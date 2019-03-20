# trim_shell.py
# script to adapter and quality trim reads 

import subprocess
import glob

# location of your Trimmomatic jar
TRIM_LOC = '/tools/Trimmomatic-0.38/trimmomatic-0.38.jar'
# location of your adapters file 
ADAPTER_LOC = '/tools/Trimmomatic-0.38/adapters/adapters.fa'

for file_name in glob.glob('*.fastq'):
  base = file_name.split('_')[0]
  # MINLEN: is the minimum read length allowed, this should be bumped higher but we had to drop this to align 37bp reads, if you're doing this with 150bp runs you should up this to 50+ 
  trim_cmd = 'java -jar ' + TRIM_LOC + ' SE -threads 42 ' + file_name + ' ' + base + '_trimmed.fastq ILLUMINACLIP:' + ADAPTER_LOC + ':2:30:10 LEADING:10 TRAILING:10 SLIDINGWINDOW:4:20 MINLEN:30 > ' + base + '.trim.log'
  # SLIDINGWINDOW:size:average_quality
  subprocess.call(trim_cmd, shell=True)
subprocess.call('cat *.trim.log > final_trim_stats.txt', shell=True)

