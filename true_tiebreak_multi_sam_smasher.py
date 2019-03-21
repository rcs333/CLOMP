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
ncbi = NCBITaxa()
import timeit

VERSION = '1.0'

# the location of the fully built krakenuniq-report program 
KRAKEN_LOC = '/tools/krakenuniq/krakenuniq-report'
#KRAKEN_DB_LOC = '/hd2/krakenuniq_all/krakenuniq_all2/all/'
KRAKEN_DB_LOC = '/hd2/kraktest/'

# takes a list of the form[[taxid,edit distance],[taxid,edit distance],[ect]] and breaks ties
# returning a single taxid
def tie_break(taxid_list):
	list_of_lineages = []
	list_of_assignments = []
	actual_taxid_list = []
	score_list = []
	
	# taxid list is a list of tuples so we go through and put all the edit distances into a new list
	for id in taxid_list:
		score_list.append(id[1])
	
	# then we figure out the 'best' edit distance. Edit distances are like golf, lower is better
	# except for negative numbers, which is why any edit distance that's negative gets moved to +100 
	# before passing into this function
	id_max = min(score_list)
	
	# create a new list of taxids only holding taxids with the best edit distance or one greater
	for id in taxid_list:
		if id[1] <= id_max + 1:
			actual_taxid_list.append(id[0])
	
	# then overwrite the original variable with just a list of good taxids, taxid_list is now [taxid,taxid,taxid]
	taxid_list = actual_taxid_list
	# if the read ever aligned to human just call it human and move on
	if '9606' in taxid_list:
		return 9606

	is_read_unassigned = True 
	# go through every taxid that got assigned to the current read 
	for id in taxid_list:
		# only attempt to call lineage on hits that aligned to something 
		if id != '*':
			# get list of NCBI taxonomy lineage for the taxid 
			try:
				lineage = ncbi.get_lineage(id)
				# if the taxid comes from either 'unclassified sequences' or 'other sequences' 
				# we'll not count that assignment 
				if (12908 in lineage) or (28384 in lineage):
					lineage = []
			except:
				lineage = []
			# walk backwards through the taxonomically lineage and remove leaves until 
			# we see a taxid that is not 'no rank'
			for x in range(len(lineage) -1, -1, -1):
				if ncbi.get_rank([lineage[x]])[lineage[x]] == 'no rank':
					lineage.pop(x)
				else:
					break
			
			# if we either couldn't get an NCBI linage or removed all nodes 
			# call this taxid unassigned 
			if len(lineage) == 0:
				list_of_assignments.append('*')
				
			else:
				# if we got an assignment, change boolean flag to note that we have at least one
				# good tax assignment for this read 
				is_read_unassigned = False
				# store complete lineage with 'no ranks' removed
				list_of_lineages.append(lineage)
				# store the exact assignment - we only append to the linage list if the taxid is valid
				list_of_assignments.append(lineage[-1])
	
	# if we never saw a taxid other than '*' we'll return this and handle it later
	if is_read_unassigned:
		return '*'
	# make sure that we have a list with only valid ranked NCBI taxonomies 
	for thing in list_of_assignments:
		if thing == '*':
			list_of_assignments.remove('*')

	if len(list_of_assignments) == 0:
		return '*'
	
	# count the number of times that the most common taxid was 'voted' for
	most_common_id_count = list_of_assignments.count(max(set(list_of_assignments), key=list_of_assignments.count))
	
	# if less than 90% of hits go to the best hit then we just perform an intersection
	# this covers both the case of assignments being higher up the tree, or to different species 
	# if we only have two hits just intersect them - this won't change anything if the assignments
	# are the same and if they're different it'll place them correctly
	total_counts = len(list_of_assignments)
	if most_common_id_count <= total_counts - (total_counts / 10):

		lineage_intersect = set.intersection(*map(set,list_of_lineages))
		
		# pull taxid with the longest lineage this is just because the we converted the list into a set which un-orders it
		# so we want to assign the read to the most specific taxid that survived the intersection 
		d = {}
		# get map of taxid lineage lists
		for the_value in lineage_intersect:
			d[the_value] = len(ncbi.get_lineage(the_value))
		# figure out the most specific taxid
		assigned_hit = max(d.iteritems(), key=operator.itemgetter(1))[0]
		return assigned_hit
		
	# if we didn't do an intersection just return the most common taxid 
	else: 
		return max(set(list_of_assignments), key=list_of_assignments.count)

	
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
	kraken_report_cmd = KRAKEN_LOC + ' --db ' + KRAKEN_DB_LOC + ' --taxon-counts ' + basename + '_temp_kraken.tsv > ' + basename + '_truebreak_final_report.tsv'
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
		
		for read_key in read_to_taxids_map.keys():
			tax_assignment = tie_break(read_to_taxids_map[read_key])
			if tax_assignment == '*':
				e.write(read_key + '\t' + '0' + '\t' + reads_seq_map[read_key] + '\n')
				unass_count += 1
			else:
				g.write(read_key + '\t' + str(tax_assignment) + '\t' + reads_seq_map[read_key] + '\n')
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
