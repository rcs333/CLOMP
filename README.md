# CLOMP - CLinically Okay Metagenomic Pipeline

A entirely open source, fast, accurate, multi-sample and highly configurable end to end metagenomics pipeline for the people.

### Table of Contents
1. [This Introduction and general description](#Introduction)

2. [Insallation and configuration guide](#Installation)

3. [Execution and run guide](#Running)

4. [Technicall details and my ramblings](#Technical)

# Introduction

This is the publically available source code and documentation for CLOMP - UW Virology's fully functional metagenomic pipeline. CLOMP takes raw sequencing reads straight off the sequencer and taxonomically assigns as many reads as possible. While CLOMP is currently relatively difficult to set up, once properly configured you can get accurate taxonomical classification of an entire MiSeq run overnight. 

Broadly the execution of the pipeline is broken down into four steps:
### 1. Read quality filtering and adapter trimming 

First, we take all reads and remove adapter and low quality sequences from the ends of the reads, and throw out reads that are generally low quality or too short. This is done using [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic) and a custom adapter file. 

### 2. Host read subtraction

Next, reads are aligned to a host geneome and thrown out. As this pipeline is designed to be 'clinically okay' and for human consuption, by default we align to the Human Genome (hg38 p.12). Depending on the sample type, this may remove almost all reads. This is great because it allows the rest of the pipeline to actually run relatively quickly. However, with very little effort you can configure CLOMP to filter out almost any host genome (as long as it's actually been sequenced). Host subtraction is performed with [bowtie2](https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.3.4.3)

### 3. Alignment of every read to NT

After we've got good quality non-host reads we take all remaining reads and align them to NCBI's NT. This is performed with the lovely tool [SNAP](http://snap.cs.berkeley.edu/). While there are a lot of important differences between SNAP and BLAST, you can basically just think of this step as running a super fast blast search of every read against NT. The fact that we align all reads to NT allows us to accuratley place many more reads than approaches that use a smaller database. 

### 4. Tiebreaking and visualization

Finally we give each read one taxonomical assignment. Again think of this as manually scrolling through some blast results and determining the ultimate origin of each read. Of course, reads that have highly conserved sequence (like ribosomal RNA) cannot be accuratley placed at the species level. A big part of the tiebreaking step is placing each read at the most specific taxonomical description that we can be confidant. For example, if a read is aligning super well to both e. coli and pseudomonas then we can't say for sure the speciees from which the read came, but we can be pretty sure it's a proteobacteria. Once this has been done for every read the final results get packed up for consumption with the amazing vizualisation tool [Pavian](https://github.com/fbreitwieser/pavian)

# Installation
## Required System Specs
Linux OS (This can theoretically work on a Mac but do yourself a favor and don't)

~3Tb of hard drive space. This can be on discontinous drives but I'm not including instructions for how that would work. This space is required just for holding your indexes and programs, so if you want to actually do processing and save the output for a sequencing run you'll need more hard drive space to hold your input data.

Between 32 and 5000 Gigabytes of RAM. The more RAM you have the quicker this will go. This was originally developed with 256Gb RAM but I've also been able to sucesfully configure this pipeline on an old iMac with 16 cores and 32Gb of RAM. (More on that later) 

An internet connection is required to do automatic alingment of results. And to download and install tools and databases.

## Installation and setup
The installation of my scripts is as simple as copying them to your computer - functionally all they are is wrappers for other programs. The first and hardest step of this will be configuring and setting up all the other programs and their associated databases. Once that's done you have to change my scripts to point to your local database and tool locations, then you need to add the main script to the PATH and you'll be good to go. This guide tries as hard as possible to make this process accessible but you will be immensely helped if you can set up and run bowtie2, snap, and Trimmomatic on your own without problems. However, once you've got all the dependencies installed correctly pretty much all you have to do is change a few lines of code to point to these and then you'll be good to go. 

### The easy dependencies
#### Approximate time: ~5 - 0 minutes (these are already installed on a lot of systems)
[Python](https://www.python.org/downloads/) (Any version is fine) Also need the ete3 module `python -m pip install ete3`

[JRE](https://www.oracle.com/technetwork/java/javase/downloads/jre8-downloads-2133155.html) (To run Trimmomatic)

libz `sudo apt install libz-dev`

samtools `sudo apt install samtools ` 

### 1. Getting Trimmomatic set up. 
#### Approximate time: ~5 minutes
You can download [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic) and unzip it to anywhere on your computer, you'll need to edit trim_shell.py to reflect the location of where Trimmomatic is living. So change line four to read `TRIM_LOC = '/your/path/to/trimmomatic-version.jar` Then you need to download the provided adapter file (adapter.fa) and change line 5 to read `ADAPTER_LOC = '/path/to/adapters.fa`


### 2. Get human host filtering set up. 
#### Approximate time: Setup ~2 hours, Index build ~1 day 
First you need to download and install [bowtie2](https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.3.4.3).

Next you need to download a copy of the human genome. I use human genome hg38 (ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.27_GRCh38.p12) You want to download the two files that end in `_genomic.fna.gz` and `_rna_from_genomic.fna.gz` This may take quite a while to download. Then you need to unzip these two files. If you're trying to add other organisms or sequences to you host filtering download them now. 
Next you need to index your host genome, the final index will take about 25Gb on your drive and requires about 4Gb of RAM to build. I like to concatenate all of my host genome files, this isn't actually necessary but I do it anyways. `cat GCA_000001405.27_GRCh38.p12_genomic.fna GCA_000001405.27_GRCh38.p12_rna_from_genomic.fna > hg38.fna`

Then build with `bowtie2-index hg38.fa hg38` This takes about 5 hours on our server so might take considerabbly longer depending on specs. Once this completes you'll need to go edit line 5 of host_filter.py to read `DB_LOCATION = '/path/to/hg38'`

### 3. SNAP Index to NT.
#### Approximate time: Setup ~ 2 hours, Index build ~3 days (might be MUCH longer)
This is absolutely the hardest part of the whole process. We use SNAP sequence aligner to quickly map reads to all of NT. However, SNAP requires that you be able to load the ENTIRE database into RAM. The final built database for NT is a bit less than 4TB in size. We have a 244Gb RAM machine and I split NT into 14 different ~12Gb fasta files and built an index for each of these (each of these indexes is about 198Gb once built). Then we sequentially run SNAP on each of these databases to produce the final output. If you have less RAM then you'll need to make more chunks, if your RAM usage starts paging during the index build it'll probobly never finish. 

##### If this moves to production I'm going to upload all taxid-appended chunks to the release page and then users can just download those and concat them as neccesary 

First you need to download NCBI's NT database (ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/) You want the one that says nt.gz

I used [pyfasta](https://pypi.org/project/pyfasta/) to split the files, but it really doesn't matter how you go about it. 
You need to split NT into however many peices you need. 14 peices load into RAM at about 200GB each but take about 220Gb to build. (Building the index takes more RAM than aligning to it). Each peice of NT takes up 14Gb of size on my computer. So if you assume that you need about (size of fasta file X 15) of RAM in order to build the database you can figure out how many peices you need to be able to fit the database into RAM. 

Next you need to process NT. I've included parse_nt.py to allow you to do this, what this does is remove sequences less than 500nt from the database as well as append the ncbi taxid to the header of each sequence. I also currate a list of 'bad' accession numbers and prune these out. This script requires that you have a locally downloaded copy of NCBI's taxonomy (ftp://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid) (You want the file called nucl_gb.accession2taxid.gz). My code works on fasta files that have already been split and are named like nt.fna.00 ect. You just need to make sure that nucl_gb.accesion2taxid is in the same folder as your nt files and parse_nt.py. Then run `python parse_nt.py` At this point if are short on disk space you can delete the nt.gz, nt.fna as well as the unprocessed chunks.

The other constraint for building the database is going to be the size on your hard drive. Currently all of NT takes around 2.5TB of hard drive space and the more chunks you make the more HD space you need.  

Next download and install [my build of SNAP](https://github.com/rcs333/snap). The original website and manual can be found [here](http://snap.cs.berkeley.edu/). The only difference my version has is that my version never writes SAM headers into the output files which when you're aligning to NT saves a huge amount of I/O time and hard drive space. 

Then use `snap index nt.00 nt.00.snapindex -tNumThreads` If snap prompts you that you need to increase the location size than do so with `snap index nt.00 nt.00.snapindex -tNumThreads -locationSize 5` and hope that the index build finishes. SNAP is nice because it'll print really helpful timing results and you can see if the index build looks like it's going to finish before the entropic heat death of the universe. Obviously you need to run this command for each of your chunks, so if you have 14 chunks you'll need to build each of these. 

Once your database is built you'll also need to modify multi_snap_shell.py. If you have more than 14 databases you need to add those numbers to line 7 - just make sure they're encased in quotes and comma seperated. Then change line 8 to represent where your databases are stored and change line 9 to point to where you installed SNAP.

### Configuring the tiebreaking
#### approximate time quick 
The tiebreaking code requires some parts of a krakenuniq tiebreaker. Either way just download and build the smallest krakenuniq database you can find and then make sure tie_break.py points to the right places. Lines 17 and 18 need to point to your krakenuniq-report script and the database respectivey.

You also need to make sure that your python installation has ete3 installed with a another copy of NCBI's taxonomy downloaded by this program as well. 
Then finally just run
`python trim_host.py; python host_filter.py; python snap_shell.py; python tie_break.py`
This assumes that you want to host trim then host filter and that all your files are in a folder with the extension .fastq. This can all be tweaked by changing the wildcard expressions in the opening for loops. Additionally, if you created more than 14 SNAP databases you'll have to modify the loops. 

# Running

# Technical 
