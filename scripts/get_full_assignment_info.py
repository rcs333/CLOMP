# Script that takes a 'base name' as well as a read name and parses all the sam files to get every assignment this read got
# base name by convention is everyhting before the _ in the original read files 

import argparse
import glob
if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='')
	parser.add_argument('file_base', help = 'everything before the id numbers, search will be *file_base*.sam')
	parser.add_argument('read_name', help = 'read id that you wanna pull')
	
	args = parser.parse_args()
	file_base = args.file_base
	read_name = args.read_name
	
	line_list_to_save = []
	# loop through every line in every sam file matching the wildcard expansion and save whole line
	for file_name in glob.glob('*' + file_base + '*.sam'):
		for line in open(file_name):
			if line.split('\t')[0] == read_name:
				line_list_to_save.append(line)
	
	# write all lines into a new file 
	g = open(read_name + '.txt', 'w')
	for item in line_list_to_save:
		g.write(item)