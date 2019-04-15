#	true_tiebreak_multi_sam_smasher V 1.0
#	Copyright (C) 2019 Ryan C. Shean 

#		This program is free software: you can redistribute it and/or modify
#		it under the terms of the GNU General Public License as published by
#		the Free Software Foundation, either version 3 of the License, or
#		(at your option) any later version.

#		This program is distributed in the hope that it will be useful,
#		but WITHOUT ANY WARRANTY; without even the implied warranty of
#		MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#		GNU General Public License for more details.

#		You should have received a copy of the GNU General Public License
#		along with this program. If not, see <http://www.gnu.org/licenses/>.


# the final step of the CLOMP pipeline, takes all .sam output files and parses all alignments for each 
# read and outputs a final taxonomic assignment for each read, writes a pavian compatible report file
# as well as a file containing read names of all unaligned reads and a file with read_name -> taxid for 
# all assigned reads, also creates a folder for each sample and puts unique reads mapping to each taxid
# into that folder 

import subprocess
import glob
import operator
from collections import Counter
from ete3 import NCBITaxa
# loads ncbi taxonomy database into a sqlite database
ncbi = NCBITaxa()
import timeit

VERSION = '1.0'

# the location of the fully built krakenuniq-report program 
KRAKEN_LOC = '/tools/krakenuniq/krakenuniq-report'
KRAKEN_DB_LOC = '/hd2/kraktest/'


# takes a list in the form [[taxid, edit_distance], [taxid, edit_distance], ..]
# returns a single taxid as the final assignment 
def tie_break(taxid_list):
	# Figure out the best (lowest) edit distance that our read aligned with 
	score_list = []
	actual_taxid_list = []	
	for id in taxid_list:
		score_list.append(id[1])
	# we want to keep alignments that are within one of our best edit distance so add one here
	best_edit_distance = min(score_list) + 6
	
	# discard all assignments that are not within 1 of the best edit distance 
	for id in taxid_list:
		if id[1] <= best_edit_distance:
			actual_taxid_list.append(id[0])
	
	# overwrite variable to type less
	taxid_list = actual_taxid_list

	lineage_list = []
	# get full taxonomic lineages for every taxid that was assigned to this read 
	for id in taxid_list:
		# not all taxonomic assignments have lineages
		try:
			# returns a list containing the complete NCBI taxonomical lineage as integers for a given taxid 
			lineage = ncbi.get_lineage(id)
			# discard any lineages that contain 'other sequences' or 'artificial sequences' or 'environmental samples'
			if (12908 in lineage) or (28384 in lineage) or (48479 in lineage):
				lineage = [] 
		except:
			lineage = []

		# if we got a lineage add it to our lineage list
		# lists have implicit boolean(ness) in python :D 
		if lineage:
			lineage_list.append(lineage)
	
	# if we never added anything to our lineage list return '*' (unassigned)
	if not lineage_list:
		return '*'
	
	# count the number of times each node occurs in all lineages
	taxid_to_count_map = {}
	for each_lineage in lineage_list:
		for each_taxid in each_lineage:
			if each_taxid in taxid_to_count_map:
				taxid_to_count_map[each_taxid] += 1
			else:
				taxid_to_count_map[each_taxid] = 1
	
	# go through all taxid -> count map and throw out taxids that occur in less than
	# 90% of cases (with integer rounding) basically:
	# if less than 10 total assignments
	# 	thrown out all taxids appearing once or less
	# if between 10 and 20 total assignments
	# 	throw out all taxids appearing twice or less
	# and so on and so forth adding +1 every 10 
	num_assignments = len(lineage_list)
	# ---------Changing the tie breaking logic ---------
	# Swap this line for more aggressive speciation at the cost of more incorrect speciation 
	#threshold = num_assignments - ((num_assignments /10) + 1) 
	threshold = num_assignments
	surviving_taxids = []
	for taxid_key in taxid_to_count_map:
		if taxid_to_count_map[taxid_key] >= threshold:
			surviving_taxids.append(taxid_key)
	
	# if we found no taxids at above ~90% return unassigned. This should rarely happen
	# because essentially all taxids start at 1. the case of a read aligning to one phage and one
	# bacteria would result in the read getting assigned to 1 
	if len(surviving_taxids) == 0:
		return '*'
	
	# go through and find the remaining taxid with the longest lineage ( most specific taxid left)
	d = {}
	# get map of taxid lineage lists
	for the_value in surviving_taxids:
		d[the_value] = len(ncbi.get_lineage(the_value))
	
	# return the longest hit
	assigned_hit = max(d.iteritems(), key=operator.itemgetter(1))[0]
	return assigned_hit



# use the krakenuniqu report function on a full ncbi taxdmp to actually include EVERY taxid in the report 
# this lets us assign reads directly to any species that's in NCBI taxonomy 

def new_write_kraken(basename, final_counts_map, num_unassigned):
	print('Preparing output for ' + basename)
	# we write a file in the form taxid\tcount 
	l = open(basename + '_temp_kraken.tsv', 'w')
	
	# initialize with the number of unassigned, we'll need to add human host filtering in earlier
	# because some reads will get tie broken to human 
	
	l.write('0\t' + str(num_unassigned))
	# write the rest of the taxids to the file
	for key in final_counts_map.keys():
		l.write('\n' + str(key) + '\t' + str(final_counts_map[key]))
	
	# this close is required so we get an EOF before passing the file to kraken-report 
	l.close()
	
	# kraken-report creates a file that Pavian likes - we name the file base_final_report.tsv
	kraken_report_cmd = KRAKEN_LOC + ' --db ' + KRAKEN_DB_LOC + ' --taxon-counts ' + basename + '_temp_kraken.tsv > ' + basename + '_final_report.tsv'
	subprocess.call(kraken_report_cmd, shell=True)

	 
if __name__ == '__main__':
	# generate a list of all the sam files 
	file_list = glob.glob('*.sam')
	base_col = set()
	main_start_time = timeit.default_timer()
	# create a list of the unique base names 
	

	for item in file_list:
		base_col.add(item.split('_')[0])
	# iterate through each base name, keeping the same map{} files for each database, so these will get LARGE
	for file_base in base_col:
		base_start_time = timeit.default_timer()
		read_to_taxids_map = {}
		reads_seq_map = {}
		# iterate through each of the db sam output but keep the same maps
		for db_num in ['00','01','02','03','04','05','06','07','08','09','10','11','12','13']:
			file_name = file_base + '_' + db_num + '.sam'
			file_start_time = timeit.default_timer()
			print('Reading in ' + file_name)

			# Loop through every line in the sam file backwards and stop once we hit the header
			line_count = 0
			for line in open(file_name):
				line_count += 1
				if line_count > 3:

					line_list = line.split('\t')
					current_read = line_list[0]
					snap_assignment_of_current_read = line_list[2]
					sequence_of_current_read = line_list[9]
				
					# if read was unassigned append '*' to taxid list 
					if snap_assignment_of_current_read == '*':
						# higher edit distance than anything else makes sure this gets parsed out 
						current_read_taxid = [snap_assignment_of_current_read,100]
					else:
						current_read_taxid = [snap_assignment_of_current_read.split('#')[-1],int(line_list[12].split(':')[-1])]
					
					# create a map of read_name -> [[taxid, edit distance ], [taxid, edit distance], ect..]
					if current_read in read_to_taxids_map:
						read_to_taxids_map[current_read].append(current_read_taxid)
					else: 
						# if this is the first time we've seen this read add the sequence to the list
						# and initialize the read -> taxid_list map  
						read_to_taxids_map[current_read] = [current_read_taxid]
						# also store the read and the sequence, this does need to be in a map 
						reads_seq_map[current_read] = sequence_of_current_read
			file_runtime = str(timeit.default_timer() - file_start_time)
			print('Reading in file ' + file_name + ' took ' + file_runtime)
		per_base_runtime = str(timeit.default_timer() - base_start_time)
		print(file_base + ' took ' + per_base_runtime + ' in total to read')
		final_assignment_counts = {}

		print('Breaking ties for ' + file_base)
		tie_break_start = timeit.default_timer()
		# now we're done with the loop and we have a map with reads to list of taxids assigned to them
		g = open(file_base + '_assignments.txt', 'w')
		e = open(file_base + '_unassigned.txt','w')
		unass_count = 0
		taxid_to_read_set = {}
		
		# go through every read 
		for read_key in read_to_taxids_map.keys():
			# assign one taxid to read 
			tax_assignment = tie_break(read_to_taxids_map[read_key])
			
			# read is unassigned write it to the unassigned file 
			if tax_assignment == '*':
				e.write('>' + read_key + '\n' + reads_seq_map[read_key] + '\n')
				unass_count += 1
			# otherwise write it out to the read by read assignment file 
			else:
				g.write(read_key + '\t' + str(tax_assignment) + '\t' + reads_seq_map[read_key] + '\n')
				# create a mapping of taxid -> unique reads 
				if str(tax_assignment) in taxid_to_read_set:
					taxid_to_read_set[str(tax_assignment)].add(reads_seq_map[read_key],)
				else:
					taxid_to_read_set[str(tax_assignment)] = set([reads_seq_map[read_key]])
					
			
				if tax_assignment in final_assignment_counts:
					final_assignment_counts[tax_assignment] += 1
				else:
					final_assignment_counts[tax_assignment] = 1
				
		g.close()
		e.close()
		subprocess.call('mkdir ' + file_base.split('.')[0], shell=True)
		for id in taxid_to_read_set.keys():
			f = open(file_base.split('.')[0] + '/' + str(id) + '_uniques.txt', 'w')
			count = 0
			for item in taxid_to_read_set[id]:
				f.write('>' + str(count) + '\n')
				f.write(item + '\n')
				count += 1
			f.close()
		
		tie_break_time = str(timeit.default_timer() - tie_break_start)
		print('Tie breaking ' + file_base + ' took ' + tie_break_time)
		# write the final output report and then go back to the start of the loop 
		new_write_kraken(file_base, final_assignment_counts, unass_count)
