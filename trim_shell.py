#	trim_shell.py V 1.0
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




import subprocess
import glob

# location of your Trimmomatic jar
TRIM_LOC = '/tools/Trimmomatic-0.38/trimmomatic-0.38.jar'
# location of your adapters file I recommend using the one included in the GitHub repo
# but if you know what you want to remove you can just substitute a custom file here
ADAPTER_LOC = '/tools/Trimmomatic-0.38/adapters/adapters.fa'

for file_name in glob.glob('*.fastq'):
  base = file_name.split('_')[0]
  # MINLEN: is the minimum read length allowed, this should be bumped higher but we had to drop this to align 37bp reads, if you're doing this with 150bp runs you should up this to 50+ 
  trim_cmd = 'java -jar ' + TRIM_LOC + ' SE -threads 42 ' + file_name + ' ' + base + '_trimmed.fastq ILLUMINACLIP:' + ADAPTER_LOC + ':2:30:10 LEADING:10 TRAILING:10 SLIDINGWINDOW:4:20 MINLEN:30 > ' + base + '.trim.log'
  # SLIDINGWINDOW:size:average_quality
  subprocess.call(trim_cmd, shell=True)
subprocess.call('cat *.trim.log > final_trim_stats.txt', shell=True)

