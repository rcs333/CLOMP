# host_filter.py 
# go through trimmed fastq files and remove all reads that align to the bowtie2 index stored at DB_LOCATION
# obviously most of the time this is just hg38 but if you wanna do giraffes or whatever just build a bowtie2
# index and point to it with the DB_LOCATION variable 

import subprocess
import glob
import argparse 

DB_LOCATION = '/db/hg38/hg38_spiro'

# go through every fastq file that has been trimmed 
for r1 in glob.glob('*_trimmed.fastq'):
  base = r1.split('_')[0]
  #r2 = r1.split('R1')[0] + 'R2' + r1.split('R1')[1]
  print('Host filtering for sample: ' + base)
  # Bowtie2 options are as follows (--very-sensitive-local
  # -D 20 -R 3 -N 0 -L 20 -i S,1,0.50
  #-D seed extensions (max number that can fail in a row) 
  #-R maximum reseed attempts (max number or reattempts)
  #-N number of mismatches per seed
  #-L seed length
  #-i interval between seeds
  # other options:
  # --threads duh 
  # -x location of database 
  # -q align fastq files
  # -U unpaired reads
  # -S output as a sam file 
  align_cmd = 'bowtie2 --very-sensitive-local --threads 42 -x ' + DB_LOCATION + ' -q -U ' + r1 + ' -S ' + base + '_mappedSam 2>&1 | tee -a ' + base + '.log'
  subprocess.call(align_cmd, shell=True)

  sort_cmd = 'samtools view -Sb -@ 42 ' + base + '_mappedSam > ' + base + '_mappedBam' 
  subprocess.call(sort_cmd, shell=True)

  sort_cmd2 = 'samtools sort -@ 42 ' + base + '_mappedBam ' + base +  '_sorted'
  subprocess.call(sort_cmd2, shell=True)

  index_cmd = 'samtools index ' + base + '_sorted.bam'
  subprocess.call(index_cmd, shell=True)

  r1_cmd = 'samtools view -@ 42 -F 0x40 ' + base + '_mappedBam | awk \'{if($3 == \"*\") print \"@\" $1 \"\\n\" $10 \"\\n\" \"+\" $1 \"\\n\" $11}\' > ' + base + '_hf.trimmed.fastq'
  #r2_cmd = 'samtools view -f 0x40 ' + base + '_mappedBam | awk \'{if($3 == \"*\") print \"@\" $1 \"\\n\" $10 \"\\n\" \"+\" $1 \"\\n\" $11}\' > ' + base + '_host_filtered_R2.fastq'
  subprocess.call(r1_cmd, shell=True)
  subprocess.call('rm ' + base + '_mappedSam', shell=True)
  subprocess.call('rm ' + base + '_mappedBam', shell=True)
  subprocess.call('rm ' + base + '_sorted.bam', shell=True)
  #subprocess.call(r2_cmd, shell=True)

# to output coverage stats
stats_cmd = 'for f in *.log;  do echo $f >> final_bowtie.txt; cat $f >> final_bowtie.txt; done'
subprocess.call(stats_cmd, shell=True)
