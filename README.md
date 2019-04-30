# CLOMP - CLinically Okay Metagenomic Pipeline

A entirely open source, fast, accurate, multi-sample and highly configurable end to end metagenomics pipeline for the people.

### Table of Contents
1. [Introduction and general description](#Introduction)

2. [Installation and configuration guide](#Installation)

3. [Execution and run guide](#Running)

4. [Technical details and other ramblings](#Technical)

# Introduction

This is the publically available source code and documentation for CLOMP - UW Virology's fully functional metagenomic pipeline. CLOMP takes raw sequencing reads straight off the sequencer and taxonomically assigns as many reads as possible. While CLOMP is currently relatively difficult to set up, once properly configured you can get accurate taxonomical classification of an entire MiSeq run overnight. 

Broadly the execution of the pipeline is broken down into four steps:
### 1. Read quality filtering and adapter trimming 

First, we take all reads and remove adapter and low quality sequences from the ends of the reads, and throw out reads that are generally low quality or too short. This is done using [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic) and a custom adapter file. 

### 2. Host read subtraction

Next, reads are aligned to a host genome and those that match are removed from further analysis. This pipeline is primarily designed for human clinical samples so by default we align to the Human Genome (hg38 p.12). Depending on the sample type, this may remove almost all reads. This is great because it allows the rest of the pipeline to run relatively quickly. However, with very little effort you can configure CLOMP to filter out almost any host genome (as long as it's actually been sequenced). Host subtraction is performed with [bowtie2](https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.3.4.3).

### 3. Alignment of every read to NT

After we've got good quality non-host reads we take all remaining reads and align them to NCBI's NT. This is performed with the lovely tool [SNAP](http://snap.cs.berkeley.edu/). While there are a lot of important differences between SNAP and BLAST, you can basically just think of this step as running a super fast blast search of every read against NT. The fact that we align all reads to NT allows us to accurately place many more reads than approaches that use a smaller database. 

### 4. Tiebreaking and visualization

Finally we give each read one taxonomical assignment. Think of this like manually scrolling through some blast results and determining the ultimate origin of each read. Of course, reads that have highly conserved sequence (like ribosomal RNA) cannot be accurately placed at the species level. A big part of the tiebreaking step is placing each read at the most specific taxonomical description that we can be sure of. For example, if a read is aligning super well to both E. coli and Pseudomonas, then we can't say for sure the species from which the read came, but we can be pretty sure it's a Proteobacteria. Once this has been done for every read the final results get packed up for consumption with the amazing visualization tool [Pavian](https://github.com/fbreitwieser/pavian).

Additional scripts are included for some optional steps, like pulling all reads aligning at or below a given taxa as well as generating coverage maps.

# Installation
## Required System Specs
Linux OS (This can theoretically work on a Mac but **do yourself a favor and don't**).

~3Tb of hard drive space. This can be on discontinuous drives but I'm not including instructions for how that would work. This space is required just for holding your indexes and programs, so if you want to actually do processing and save the output for a sequencing run you'll need more hard drive space to hold your input data.

Between 32 and 5000 Gigabytes of RAM. The more RAM you have the quicker this will go. This was originally developed with 256Gb RAM but I've also been able to successfully configure this pipeline on an old iMac with 16 cores and 32Gb of RAM. (More on that later) 

An internet connection is required to download and install tools and databases as well as for automatic alignment of results.

## Installation and setup
The installation of my scripts is as simple as copying them to your computer - functionally all they are is wrappers for other programs. The first and hardest step of this will be configuring and setting up all the other programs and their associated databases. Once that's done you have to change my scripts to point to your local database and tool locations, then you need to add the main script to the PATH and you'll be good to go. This guide tries as hard as possible to make this process accessible but you will be immensely helped if you can set up and run bowtie2, snap, and Trimmomatic on your own without problems. However, once you've got all the dependencies installed correctly pretty much all you have to do is change a few lines of code to point to these and then you'll be good to go. 

### The easy dependencies
#### Approximate time: <5 minutes (these are already installed on a lot of systems)
[Python](https://www.python.org/downloads/) (Any version is fine) Also need the ete3 module `python -m pip install ete3`

[JRE](https://www.oracle.com/technetwork/java/javase/downloads/jre8-downloads-2133155.html) (To run Trimmomatic)

libz `sudo apt install libz-dev`

samtools `sudo apt install samtools ` 

### 1. Getting Trimmomatic set up. 
#### Approximate time: ~5 minutes
You can download [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic) and unzip it to anywhere on your computer. You'll need to edit trim_shell.py to reflect the location of where Trimmomatic is living. So change line 4 to read `TRIM_LOC = '/your/path/to/trimmomatic-version.jar` Then you need to download the provided adapter file (adapter.fa) and change line 5 to read `ADAPTER_LOC = '/path/to/adapters.fa`


### 2. Get human host filtering set up. 
#### Approximate time: ~1 hour
First you need to download and install [bowtie2](https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.3.4.3).

Next you need to download a copy of the human genome. I use [human genome hg38](ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.27_GRCh38.p12). You want to download the two files that end in `_genomic.fna.gz` and `_rna_from_genomic.fna.gz` This may take quite a while to download. Then you need to unzip these two files. If you're trying to add other organisms or sequences to your host filtering download them now. 
Next you need to index your host genome. The final index will take about 25Gb on your drive and requires about 4Gb of RAM to build. I like to concatenate all of my host genome files, this isn't actually necessary but I do it anyways. `cat GCA_000001405.27_GRCh38.p12_genomic.fna GCA_000001405.27_GRCh38.p12_rna_from_genomic.fna > hg38.fna`

Then build with `bowtie2-index hg38.fa hg38` This takes about 5 hours on our server so might take considerably longer depending on specs. Once this completes you'll need to go edit line 5 of host_filter.py to read `DB_LOCATION = '/path/to/hg38'`.

### 3. SNAP Index to NT.
#### Approximate time: Setup ~ 2 hours, Index build ~3 days (might be MUCH longer)
This is absolutely the hardest part of the whole process. We use SNAP sequence aligner to quickly map reads to all of NT. However, SNAP requires that you be able to load the ENTIRE database into RAM. The final built database for NT is a bit less than 4TB in size. We have a 244Gb RAM machine and I split NT into 14 different ~12Gb fasta files and built an index for each of these (each of these indexes is about 198Gb once built). Then we sequentially run SNAP on each of these databases to produce the final output. If you have less RAM then you'll need to make more chunks, if your RAM usage starts paging during the index build it'll probably never finish. 

##### If this moves to production I'm going to upload all taxid-appended chunks to the release page and then users can just download those and concatenate them as necessary 

First you need to download NCBI's [NT database](ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/). You want the one that says `nt.gz`.

I used [pyfasta](https://pypi.org/project/pyfasta/) to split the files, but it really doesn't matter how you go about it. 
You need to split NT into however many pieces you need. 14 pieces load into RAM at about 200GB each but take about 220Gb to build. (Building the index takes more RAM than aligning to it). Each piece of NT takes up 14Gb of size on my computer. So if you assume that you need about (size of fasta file X 15) of RAM in order to build the database you can figure out how many pieces you need to be able to fit the database into RAM. 

Next you need to process NT. I've included parse_nt.py to allow you to do this, what this does is remove sequences less than 500nt from the database as well as append the NCBI taxid to the header of each sequence. I also curate a list of 'bad' accession numbers and prune these out. This script requires that you have a locally downloaded copy of [NCBI's taxonomy](ftp://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid). (You want the file called `nucl_gb.accession2taxid.gz`). My code works on fasta files that have already been split and are named like nt.fna.00 ect. You just need to make sure that nucl_gb.accesion2taxid is in the same folder as your nt files and parse_nt.py. Then run `python parse_nt.py` At this point if are short on disk space you can delete the nt.gz, nt.fna as well as the unprocessed chunks.

The other constraint for building the database is going to be the size on your hard drive. Currently all of NT takes around 2.5TB of hard drive space and the more chunks you make the more HD space you need.  

Next download and install [my build of SNAP](https://github.com/rcs333/snap). The original website and manual can be found [here](http://snap.cs.berkeley.edu/). The only difference my version has is that my version never writes SAM headers into the output files which when you're aligning to NT saves a huge amount of I/O time and hard drive space. 

Then use `snap index nt.00 nt.00.snapindex -tNumThreads` If snap prompts you to increase the location size than do so with `snap index nt.00 nt.00.snapindex -tNumThreads -locationSize 5` and hope that the index build finishes. SNAP is nice because it'll print really helpful timing results and you can see if the index build looks like it's going to finish before the entropic heat death of the universe. Obviously you need to run this command for each of your chunks, so if you have 14 chunks you'll need to build each of these. 

Once your database is built you'll also need to modify multi_snap_shell.py. If you have more than 14 databases you need to add those numbers to line 7 - just make sure they're encased in quotes and comma separated. Then change line 8 to represent where your databases are stored and change line 9 to point to where you installed SNAP.

### Configuring the tiebreaking
#### approximate time quick 
The tiebreaking code requires some parts of a krakenuniq tiebreaker. You need a seqid2taxid and a taxDB file, the best way to do this is to use the provided scripts on the github and then just point to that folder. You also need a database.idx and database.kdb file but these can be empty. Once you've got that setup up make sure tie_break.py points to the right places. Lines 17 and 18 need to point to your krakenuniq-report script and the database respectively.

You also will want to configure the main CLOMP.sh script. For this change the hard coded paths in this file to point to the place where you downloaded this repository. Then you can add CLOMP.sh to your PATH variable. To do this for all users `sudo nano /etc/environment` Then add ':/path/to/CLOMP_folder' to the end of the PATH line. Then save this file and re log in and all users will be able to clomp away by navigating to the folder which has their fastq files in the terminal and typing `CLOMP.sh`.

This assumes that you want to host trim then host filter and that all your files are in a folder with the extension .fastq. This can all be tweaked by changing the wildcard expressions in the opening for loops. Additionally, if you created more than 14 SNAP databases you'll have to modify the loops. 

# Running
Once everything is setup you can just run CLOMP.sh from inside a folder that contains all the samples that you would like to analyze, in .fastq format.

# Technical 

Adapter trimming and quality filtering is relatively self explanatory - the current build has the minimum read length set to 30. If you're doing runs with lengths any greater than that this should get increased. SNAP alignment is based off 20mers and really the read lengths should be above 50. This is also in the options for the SNAP aligner. 

Host filtering is done using bowtie2 because that's the first program I tried to use for host filtering and it seems to work pretty well. I'm also currently trying to host filter against the whole genome of spirometeria entericia because it's absolute garbage.

SNAP alignment uses a max edit distance of 11 and an om of 1. The default omax of 20 would potentially need to be tweaked depending on the number of databases you split into. 

### Basic tiebreaking logic
SNAP will report potentially multiple alignments for each read for each database. You'll get as many sam files as you have chunks of NT. true_tiebreak_multi_sam_smasher.py contains all the meat of tiebreaking. I've gone through several different iterations of the specific tiebreaking logic. The current method is optimized to reduce species level misclassifications, at the cost of underspeciating reads. Essentially we will only speciate reads which appear to be unambiguous. The basic tiebreaking logic is as follows.

1. SNAP is run with an edit distance of 11 and accepts alignments that are, at worst, one away from the best alignment. (So if the best alignment has an edit distance of 4, we would accept alignments with edit distance 5 - but not 6). 

2. Because each SNAP database has a different composition we go through all assignments and throw out any that are more than 6 edit distance away from the best hit. This essentially just throws out the case where we get one pretty tenuous hit in one database, this is much more common for less heavily sequenced organisims. 

Additionally we have used the following logic in the past, which does a good job of speciating reads, but at the cost of incorrectly speciating some reads, mainly due to disproportonate amounts of species present in the database. This code is contained in the comments of the tiebreaking file and guides for switching between the two is also in the comments. 

3. We then take the most specific taxonomical assignment that has 100% agreement. 
4. Aggregate all assignments for each read and only keep alignments that are the best edit distance or one worse than best edit distance
5. Any read that ever aligned to human within one edit of the best assignment will be assigned to human.
6. We then pick the most specific taxid that is shared by ~90% of assignments. The actual logic is you need at least number_of_assignments - (number of assignments // 10 + 1) reads. Where // is integer division. 


For the samples and sample types I've been testing this seems to strike a really good balance of accurately placing sequencing but also not over imputing. Of course weird edge cases pop up due to the crazy nature of the NCBI taxonomy tree.  

### Ask me about the licence ###

I believe in free software. I think it would be extrememly disingenuous to the bioinformatics community and the greater scientific public to allow what is essentially a parsing script for multiple sam files to be put under a restrictive license. This entire software is based off publically available databases and tools - I also think that this pipeline (as in the order of the steps and the arguments) works quite well and could be widely used, hopefully making people's lives easier. 
