#	multi_snap_shell.py V 1.0
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



# host_filter.py 
# go through trimmed fastq files and host subtract the reads, output new fastq without host reads and a text file
# with statistics (% reads removed ect.) 

import subprocess
import glob
import argparse 

# this variable points to your built host genome - in our case it's hg38 + the spirometra entericia genome 
# so if you want to host filter the giraffe genome or whatever you'd build that database with bowtie2 and then
# change this variable to point to it
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
  # -q fastq files are input 
  # -U use unpaired reads
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
  subprocess.call(r1_cmd, shell=True)
  # Clean up unsued files 
  subprocess.call('rm ' + base + '_mappedSam', shell=True)
  subprocess.call('rm ' + base + '_mappedBam', shell=True)
  subprocess.call('rm ' + base + '_sorted.bam', shell=True)
  subprocess.call('rm ' + base + '_sorted.bam.bai', shell=True)

# output coverage stats
stats_cmd = 'for f in *.log;  do echo $f >> final_bowtie.txt; cat $f >> final_bowtie.txt; done'
subprocess.call(stats_cmd, shell=True)
