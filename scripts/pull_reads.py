# pull reads that align to a given taxid, can pull all descendant taxa too with the --r flag 
# arguments are the _assignments.txt file with the read -> taxid mapping and then the taxid you want to pull
import argparse
from ete3 import NCBITaxa
ncbi = NCBITaxa()



if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='')
	parser.add_argument('file_to_search', help = '')
	parser.add_argument('taxid', help = '')
	parser.add_argument('--r', action='store_true', help='')
	
	args = parser.parse_args()
	search_taxid = args.taxid
	search_file = args.file_to_search
	
	taxid_search_list = [int(search_taxid)]
	
	if args.r:
		taxid_search_list =  taxid_search_list + ncbi.get_descendant_taxa(search_taxid, intermediate_nodes=True)
	
	header_list = []
	seq_list = []
	g = open(search_file.split('.')[0] + str(search_taxid) + '.fasta', 'w')
	for line in open(search_file):
		line_list = line.split('\t')
		if int(line_list[1]) in taxid_search_list:
			g.write('>' + line_list[0] + '\n')
			g.write(line_list[2])
	g.close()
	
	
