#	CLOMP - Clinically Okay Metagenomic Pipeline
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


# Exactly the same as CLOMP.sh but works on paired end input files, assumes every read has a mate
# and that both use illumina naming like our sequencer ie 'base_*R1*.fastq' and 'base_*R2*.fastq'
echo 'Unzipping fastq.gz files (if any)'
gunzip *.gz 
echo 'Done unzipping' 
python /tools/CLOMP/trim_shell.py
python /tools/CLOMP/host_filter.py
python /tools/CLOMP/multi_snap_shell.py

python /tools/CLOMP/true_tiebreak_multi_sam_smasher.py