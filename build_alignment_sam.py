#	build_alignment_sam V 1.0
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


# this script will build alignment .sam files of all species level tax assignments with above MIN_READ_CUTOFF reads going to that species
# This can only be run after the whole pipeline has finished, including tie-breaking. 
from Bio import Entrez 
from ete3 import NCBITaxa
ncbi = NCBITaxa()
import argparse
import glob
Entrez.email = 'uwvirongs@gmail.com'
import subprocess

# minimum number of reads that must be assigned to a species in order to assemble
MIN_READ_CUTOFF = 10
# number of chunks that you've split NT into (must be in the same format as your snap indexes
DB_LIST = ['00','01','02','03','04','05','06','07','08','09','10','11','12','12','13']
# This number is how many nodes upwards on the tree should we grab a list of reads to align to the species level reference
# Right now it's set to 2 so if we have a species assignment we'll pull all reads that were assigned to the family that the species belongs to 
ASSEMBLY_NODE_OFFSET = -2 

# traverse through every report file (for every sample) 
for file_name in glob.glob('*report.tsv'):

	base = file_name.split('_')[0]
	
	taxid_to_assemble = []
	# go through every line in the report file and pull all species level assignments with above MIN_READ_CUTOFF reads
	for line in open(file_name):
		line_list = line.split('\t')
		
		if line_list[3] == 'S' and int(line_list[2]) >= MIN_READ_CUTOFF:
			lineage = ncbi.get_lineage(line_list[4])
			# No eukaryotic assembly and no 'uncultured bacterium' assembly 
			if 2759 not in lineage and line_list[4] != '77133':
				taxid_to_assemble.append(line_list[4])

	# for every taxid we assemble a list of reads that are at the given species rank or lower 
	for taxid in taxid_to_assemble:
		taxid_search_list = [str(taxid)]
		taxid_search_list = taxid_search_list + ncbi.get_descendant_taxa(taxid, get_descendant_taxa=True)
		list_of_reads_to_pull = []
		for a_line in open(base + '_assignments.txt'):
			a_line_list = a_line.split('\t')
			if a_line_list[1] in taxid_search_list:
				list_of_reads_to_pull.append(a_line_list[0])
		
		acc_num_list = []
		# go through all the sam files and grab the accession numbers that the reads we got in the last loop aligned to 
		for db_num in db_list:
			for sam_line in open(base + '.hf.trimmed.fastq_' + db_num + '.sam'):
				sam_line_list = sam_line.split('\t')
				# only grab accession numbers that have 'complete genome' in the name to avoid assembling to partial segments
				if sam_line_list[0] in list_of_reads_to_pull and 'complete_genome' in sam_line_list[2]:
					acc_num_list.append(sam_line_list[2].split('.')[0])
		
		# if we couldn't find any accession numbers going to these reads with 'complete genome' just refuse to build an assembly for this taxid 
		if len(acc_num_list) == 0:
			print('No complete genome reference found, not assembling taxid ' + taxid)
			# this break just pops us out of this taxid and moves to the next one for the current sample 
			break
		
		# get the most common accession number for this taxid (in the event of ties this chooses completely randomly ) 
		most_common_acc_num = max(set(acc_num_list), key=acc_num_list.count)
		
		taxid_lineage = ncbi.get_lineage(taxid)
		# Change this number to align reads higher up the tree
		taxid_to_pull = taxid_lineage[-2]
		# pull all reads from the sample that are assigned at or below the taxonomic level we designated in the above line
		subprocess.call('python pull_reads.py ' + base + '_assignments.txt ' + str(taxid_to_pull) + ' --r', shell=True)
		
		print('Searching NCBI for Accession number:' + most_common_acc_num + ' for taxid ' + str(taxid_to_pull))
		record = Entrez.read(Entrez.esearch(db='nucleotide', term=most_common_acc_num))
		try:
			h2 = Entrez.efetch(db='nucleotide', id=record['IdList'][0], rettype='fasta', retmode='text')
		except:
			print(str(taxid) + ' did not return hits - not assembling')
			break
		# procedurally generate name for the reference fasta and the bowtie database we'll create 
		ref_fasta = base + '_' + str(taxid) + '_ref.fasta'
		ref_db = base + '_' + str(taxid) + '_db' 
		g = open(ref_fasta, 'w')
		g.write(h2.read())
		g.close()
		print('building bowtie2 index') 
		subprocess.call('bowtie2-build ' + ref_fasta + ' ' + ref_db + ' > /dev/null 2>&1 ', shell=True)
		print('done with index build aligning...')
		subprocess.call('bowtie2 -x ' + ref_db + ' -@ 42 -f -U ' + base + '_assignments' + str(taxid_to_pull) + '.fasta --no-unal > ' + base + '_' + str(taxid) + '.sam', shell=True)
		# delete bowtie index, I like to keep the refernce fasta just so it's easy to tell what this code pulled
		subprocess.call('rm *_db*.bt2', shell = True)
