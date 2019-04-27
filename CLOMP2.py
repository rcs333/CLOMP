import ast 
import subprocess
import glob
import argparse 
import operator
from collections import Counter
from ete3 import NCBITaxa
import timeit
ncbi = NCBITaxa()


def parse_ini(ini_path):
	for line in open(ini_path):
		if line[0] and line[0] != '#' and '=' in line:
			line_list = line.split('=')
			var = line_list[0].strip()
			val = line_list[1].strip()
			if var == 'PAIRED_END':
				global PAIRED_END
				PAIRED_END = ast.literal_eval(val)
			elif var == 'SAVE_MIDDLE_FILES':
				global SAVE_MIDDLE_FILES
				SAVE_MIDDLE_FILES = ast.literal_eval(val)
			elif var == 'THREADS': 
				global THREADS
				THREADS = val
			elif var == 'BASE_DELIMITER':
				global BASE_DELIMITER
				BASE_DELIMITER = val
			elif var == 'HOST_FILTER_FIRST':
				global HOST_FILTER_FIRST
				HOST_FILTER_FIRST = ast.literal_eval(val)
			elif var == 'INPUT_SUFFIX':
				global INPUT_SUFFIX
				INPUT_SUFFIX = val
			elif var == 'TRIMMOMATIC_JAR_PATH':
				global TRIMMOMATIC_JAR_PATH
				TRIMMOMATIC_JAR_PATH = val
			elif var == 'TRIMMOMATIC_ADAPTER_PATH':
				global TRIMMOMATIC_ADAPTER_PATH
				TRIMMOMATIC_ADAPTER_PATH = val
			elif var == 'SEQUENCER':
				global SEQUENCER
				SEQUENCER = val
			elif var == 'TRIMMOMATIC_OPTIONS':
				global TRIMMOMATIC_OPTIONS
				TRIMMOMATIC_OPTIONS = val
			elif var == 'TRIM_SUFFIX':
				global TRIM_SUFFIX
				TRIM_SUFFIX = val
			elif var == 'BWT_DB_LOCATION':
				global BWT_DB_LOCATION
				BWT_DB_LOCATION = val
			elif var == 'BWT_OPTIONS':
				global BWT_OPTIONS
				BWT_OPTIONS = val
			elif var == 'BWT_SUFFIX':
				global BWT_SUFFIX
				BWT_SUFFIX = val
			elif var == 'BWT_CLEAN':
				global BWT_CLEAN
				BWT_CLEAN = ast.literal_eval(val)
			elif var == 'DB_LIST':
				global DB_LIST
				DB_LIST = ast.literal_eval(val)
			elif var == 'SNAP_ALIGNER_LOCTION':
				global SNAP_ALIGNER_LOCTION
				SNAP_ALIGNER_LOCTION = val			
			elif var == 'SNAP_OPTIONS':
				global SNAP_OPTIONS
				SNAP_OPTIONS = val
			elif var == 'KRAKEN_PATH':
				global KRAKEN_PATH
				KRAKEN_PATH = val
			elif var == 'KRAKEN_DB_PATH':
				global KRAKEN_DB_PATH
				KRAKEN_DB_PATH = val
			elif var == 'EDIT_DISTANCE_OFFSET':
				global EDIT_DISTANCE_OFFSET
				EDIT_DISTANCE_OFFSET = int(val)
			elif var == 'H_STRICT':
				global H_STRICT
				H_STRICT = ast.literal_eval(val)
			elif var == 'LOGIC':
				global LOGIC
				LOGIC = val
			elif var =='WRITE_UNIQUES':
				global WRITE_UNIQUES
				WRITE_UNIQUES = ast.literal_eval(val)
			elif var == 'BLAST_REMOVE_NONHUMAN_EUKARYOTIC_READS':
				global BLAST_REMOVE_NONHUMAN_EUKARYOTIC_READS
				BLAST_REMOVE_NONHUMAN_EUKARYOTIC_READS = ast.literal_eval(val)
			elif var == 'BUILD_SAMS':
				global BUILD_SAMS
				BUILD_SAMS = ast.literal_eval(val)
				
			

def host_filter(suffix):
	for r1 in glob.glob('*' + suffix):
		base = r1.split(BASE_DELIMITER)[0]
		if PAIRED_END:
			r2 = r1.split('R1')[0] + 'R2' + r1.split('R1')[1]
			align_cmd = 'bowtie2 ' + BWT_OPTIONS + ' --threads ' + THREADS + ' -x ' + BWT_DB_LOCATION + ' -q -1 ' + r1 + + ' -2 ' + r2 + ' -S ' + base + '_mappedSam ' + ' 2>&1 | tee -a ' + base + '.log'
		else:
			align_cmd = 'bowtie2 ' + BWT_OPTIONS + ' --threads ' + THREADS + ' -x ' + BWT_DB_LOCATION + ' -q -U ' + r1 + ' -S ' + base + '_mappedSam' + ' 2>&1 | tee -a ' + base + '.log'
		subprocess.call(align_cmd, shell=True)
		
		sort_cmd = 'samtools view -Sb -@ ' + THREADS + ' ' + base + '_mappedSam > ' + base + '_mappedBam' 
		subprocess.call(sort_cmd, shell=True)
		
		sort_cmd2 = 'samtools sort -@ ' + THREADS + ' ' + base + '_mappedBam ' + base +  '_sorted'
		subprocess.call(sort_cmd2, shell=True)
		
		index_cmd = 'samtools index ' + base + '_sorted.bam'
		subprocess.call(index_cmd, shell=True)
		
		r1_cmd = 'samtools view -@ ' + THREADS + ' -F 0x40 ' + base + '_mappedBam | awk \'{if($3 == \"*\") print \"@\" $1 \"\\n\" $10 \"\\n\" \"+\" $1 \"\\n\" $11}\' > ' + base + '_R1' + BWT_SUFFIX
		subprocess.call(r1_cmd,shell=True)
		if PAIRED_END:
			r2_cmd = 'samtools view -@ ' + THREADS + ' -f 0x40 ' + base + '_mappedBam | awk \'{if($3 == \"*\") print \"@\" $1 \"\\n\" $10 \"\\n\" \"+\" $1 \"\\n\" $11}\' > ' + base + '_R2' + BWT_SUFFIX
			subprocess.call(r2_cmd, shell=True)
		if BWT_CLEAN:
			subprocess.call('rm ' + base + '_mappedSam', shell=True)
			subprocess.call('rm ' + base + '_mappedBam', shell=True)
			subprocess.call('rm ' + base + '_sorted.bam', shell=True)
			subprocess.call('rm ' + base + '_sorted.bam.bai', shell=True)
		if not SAVE_MIDDLE_FILES:
			subprocess.call('rm ' + r1, shell=True)
			if PAIRED_END:
				subprocess.call('rm ' + r2, shell=True)
		
		
		
def trim(suffix):
	for file_name in glob.glob('*' + suffix):
		base = file_name.split(BASE_DELIMITER)[0]
		if PAIRED_END:
			r1 = glob.glob(base + '*R1*.fastq')[0]
			r2 = glob.glob(base + '*R2*.fastq')[0]
			trim_cmd = 'java -jar ' + TRIMMOMATIC_JAR_PATH + ' PE -threads ' + THREADS + ' ' + r1 + ' ' + r2 + ' ' + base + '_R1_' + TRIM_SUFFIX + ' ' + base + '_R2_' + TRIM_SUFFIX + ' ' + SEQUENCER + TRIMMOMATIC_ADAPTER_PATH + TRIMMOMATIC_OPTIONS + ' > ' + base + '.trim.log'
		else:
			trim_cmd = 'java -jar ' + TRIMMOMATIC_JAR_PATH + ' SE -threads ' + THREADS + ' ' + file_name +  ' ' + base + TRIM_SUFFIX +  ' ' + SEQUENCER + TRIMMOMATIC_ADAPTER_PATH + TRIMMOMATIC_OPTIONS + ' > ' + base + '.trim.log'
		subprocess.call(trim_cmd,shell=True)
		if not SAVE_MIDDLE_FILES:
			if PAIRED_END:
				subprocess.call('rm ' + r1, shell=True)
				subprocess.call('rm ' + r2, shell=True)
			else:
				subprocess.call('rm ' + file_name)
		

def snap(suffix):
	print('here')
	count = -1 
	for db in DB_LIST:
		print('*' + suffix)
		print(glob.glob('*' + suffix))
		count +=  1
		snap_cmd = SNAP_ALIGNER_LOCTION + ' '
		
		for file_name in glob.glob('*' + suffix):
			print(file_name)
			base = file_name.split(BASE_DELIMITER)[0]
			print('Aligning ' + base + ' to ' + db)
			
			if PAIRED_END:
				r1 = base + '_R1_' + suffix
				r2 = base + '_R2_' + suffix 
				snap_cmd += ' paired ' + db + ' ' + r1 + ' ' + r2 + ' -o ' + base + '_' + str(count) + '.sam -t ' + THREADS + ' ' + SNAP_OPTIONS + ' , '
			else:
				snap_cmd += ' single ' + db + ' ' + file_name + ' -o ' + base + '_' + str(count) + '.sam -t ' + THREADS + ' ' + SNAP_OPTIONS + ' , '
				
			subprocess.call(snap_cmd,shell=True)
	if not SAVE_MIDDLE_FILES:
		subprocess.call('rm *' + suffix)
	
			
def tie_break(taxid_list):
	score_list = [] 
	actual_taxid_list = []
	for id in taxid_list:
		score_list.append(id[1])
	best_edit_distance = min(score_list) + EDIT_DISTANCE_OFFSET
	
	for id in taxid_list:
		if id[1] <= best_edit_distance:
			actual_taxid_list.append(id[0])
			
	if H_STRICT:
		if '9606' in actual_taxid_list:
			return ['9606',False]
	taxid_list = actual_taxid_list
	
	lineage_list = []
	
	for id in taxid_list:
		try:
			lineage = ncbi.get_lineage(id)
			if (12908 in lineage) or (28384 in lineage) or (48479 in lineage):
				lineage = []
		except:
			lineage = []
		
		if lineage:
			lineage_list.append(lineage)
	
	if not lineage_list:
		return ['*',False]
		
	taxid_to_count_map = {}
	for each_lineage in lineage_list:
		for each_taxid in each_lineage:
			if each_taxid in taxid_to_count_map:
				taxid_to_count_map[each_taxid] += 1
			else:
				taxid_to_count_map[each_taxid] = 1
	
	num_assignments = len(lineage_list)
	if LOGIC == 'strict':
		threshold = num_assignments
	elif LOGIC == '90':
		threshold = num_assignments - ((num_assignments /10) + 1) 
		
	surviving_taxids = []
	for taxid_key in taxid_to_count_map:
		if taxid_to_count_map[taxid_key] >= threshold:
			surviving_taxids.append(taxid_key)
			
	if len(surviving_taxids) == 0:
		return ['*',False]

	d = {}

	for the_value in surviving_taxids:
		d[the_value] = len(ncbi.get_lineage(the_value))
	
	assigned_hit = max(d.iteritems(), key=operator.itemgetter(1))[0]
	recheck = False 
	if BLAST_REMOVE_NONHUMAN_EUKARYOTIC_READS:
		assigned_lineage = ncbi.get_lineage(assigned_hit)
		if 2759 in assigned_lineage and 9604 not in assigned_lineage:
			recheck = True
	
	return [assigned_hit, recheck]


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
	kraken_report_cmd = KRAKEN_PATH + ' --db ' + KRAKEN_DB_PATH + ' --taxon-counts ' + basename + '_temp_kraken.tsv > ' + basename + '_final_report.tsv'
	subprocess.call(kraken_report_cmd, shell=True)



	
if __name__ == '__main__':
	parse_ini('CLOMP.ini')
	
	subprocess.call('gunzip *.zip',shell=True)
	
	if HOST_FILTER_FIRST:
		host_filter(INPUT_SUFFIX)
		trim(BWT_SUFFIX)
		snap_suffix = TRIM_SUFFIX
	else:
		trim(INPUT_SUFFIX)
		host_filter(TRIM_SUFFIX)
		snap_suffix = BWT_SUFFIX
	snap(snap_suffix)
		
	file_list = glob.glob('*.sam')
	base_col = set()
	main_start_time = timeit.default_timer()
	
	
	for item in file_list:
		base_col.add(item.split(BASE_DELIMITER)[0])
	
	for file_base in base_col:
		base_start_time = timeit.default_timer()
		reads_to_taxids_map = {}
		reads_seq_map = {}
		
		for sam_file in glob.glob(file_base + '*.sam'):
			file_start_time = timeit.default_timer()
			print('Reading in ' + file_name)
			
			line_count = 0
			for line in open(file_name):
				line_count += 1
				if line_count > 3:
					line_list = line.split('\t')
					current_read = line_list[0]
					snap_assignment_of_current_read = line_list[2]
					sequence_of_current_read = line_list[9]
					if snap_assignment_of_current_read == '*':
						# higher edit distance than anything else makes sure this gets parsed out 
						current_read_taxid = [snap_assignment_of_current_read,100]
					else:
						current_read_taxid = [snap_assignment_of_current_read.split('#')[-1],int(line_list[12].split(':')[-1])]
					
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
			# assign one taxid to read 
			loaded_read = reads_seq_map[read_key]
			r_list = tie_break(read_to_taxids_map[read_key])
			tax_assignment = r_list[0]
	
			if tax_assignment == '*':
				e.write('>' + read_key + '\n' + loaded_read + '\n')
				unass_count += 1
			# otherwise write it out to the read by read assignment file 
			else:
				g.write(read_key + '\t' + str(tax_assignment) + '\t' + loaded_read + '\n')
				# create a mapping of taxid -> unique reads 
				if WRITE_UNIQUES:
					if str(tax_assignment) in taxid_to_read_set:
						taxid_to_read_set[str(tax_assignment)].add(loaded_read,)
					else:
						taxid_to_read_set[str(tax_assignment)] = set([loaded_read])
					
			
					if tax_assignment in final_assignment_counts:
						final_assignment_counts[tax_assignment] += 1
					else:
						final_assignment_counts[tax_assignment] = 1
	
		g.close()
		e.close()
		if WRITE_UNIQUES:
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
		new_write_kraken(file_base, final_assignment_counts, unass_count)
	
	
	
	
	
	
	
	
		
		
		