# 1. Background
## What exome analysis is?
Simply saying, based on genetic knowledge, find out causative variants from not so big data.  
It also requires some computational skills.

## What are required to exome analysis?
1. knowledge of genetics
1. knowledge of target diseases, especially for genetic background
1. knowledge of associated databases (how to find the information you want)
1. knowledge of exome sequencing/high speed sequencers/molecular biology
1. computational skills
1. knowledge of ethical problems/privacy protection/IF and SF

Although not absolute, it becomes more important in order from top to bottom.  
This is completely my biased list, so I can agree with 3 or 4 or 5 in reverse/change.

# 2. Aim
This mini-course is made for a lecture which let you get a first tiny step to bioinformatic analysis.  
It assumes that this is made for **step by step hands-on style, but this could be applied to self studying.  
Do not be afraid. If you failed something, nothing will happen. No sample lost. Not wasting any tips/chips/gels/solutions/enzymes/antibodies. Just gain your experiences. Practice makes perfect**:wink:  

This course will cover:
- customise your macintosh environment for running bioinformatic programs
- experience simple target resequence data analysis
- understand what steps are there (alignment, remove duplicate, variant call..etc)
- how interpretate variants

Not cover:
- how algorithmn works
- how to handle massive amount data using cluster computers (SSH, Grid Job Scheduler, memory usage, disk I/O)
- deep/heuristic/complicated issues/knowledge/techniques/pitfalls (e.g. how can I sense abnormal data signs in heystacks)



# 3. Set up your mac :computer:
#### :point_right: Ideally, you should try and complete this section before hands-on
You have to set up your macintosh environment for informatic analyses.  
I know this is a first barricade to step in learning informatic skills, but this shold be done. I tried to make it as easy as possible.

## Install java (Java SE Development Kit 8u181 version)
Java is a kind of programming launguage.  
You need Java itself to run applications developed using Java launguage.

Download Java SE Development Kit 8u181 from here  
https://www.oracle.com/technetwork/java/javase/downloads/jdk8-downloads-2133151.html  
Then, install on your mac

After installation, confirm it by following procedures.  
Open Terminal.app (is located on /Application/Utilities), then type this shell command  
`$ java -version`

If succeeded, you will get following response

    java version "1.8.0_181"
    Java(TM) SE Runtime Environment (build 1.8.0_181-b13)
    Java HotSpot(TM) 64-Bit Server VM (build 25.181-b13, mixed mode)

## Install homebrew, a nice packge manager for macOS:beer:
A package manager maintains softwares (packages), such as install, update, and remove.  
Of course, you can manage your computer, but we usually use a package manager to make it easier.  
Go to this page, https://brew.sh/ then, follow the install instruction.  
Open Terminal.app (is located on /Application/Utilities), then type this shell command  
`$ /usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"`  
_The administrator password will be required in this process_

If succeeded, you will get following response

    ==> Next steps:
    - Run `brew help` to get started
    - Further documentation:
        https://docs.brew.sh

## Let's try to install new command (program, software, or package)
Wget command is a nice download utility for Web.  
`$ brew install wget`  
Try to download something.  
`$ wget https://www.dropbox.com/s/smuyzxmllmoctz1/test_variant_data_01.tsv`  
You can see the contents by this way  
`$ cat test_variant_data_01.tsv`  
Or  
`$ more test_variant_data_01.tsv`  
Or  
`$ less test_variant_data_01.tsv`  
_Push Q key for quit_

more/less is a viewer. Originally, there is more. Then, less was developped.

cat command means con**cat**enate. Concatenate multiple file, like this  
`$ cat fileA fileB`  
Let's test.  
`$ wget https://www.dropbox.com/s/5yfaiolgoi3cp3u/test_variant_data_02.tsv`  
`$ cat test_variant_data_01.tsv test_variant_data_02.tsv`  
To make it easier to see,  
`$ cat test_variant_data_01.tsv test_variant_data_02.tsv | less -S`  
_Push Q key for quit_  
-S is a option of less command. It change less behavior to chop-long-lines.  
Commands have their specific options. You can see like this.  
`$ less --help`  
`$ cat --help`


## Install softwares required for (minimal) sequence analysis
At first, type this, to cover scientific programs well  
`$ brew tap brewsci/bio`

### Install bwa
bwa for aligning reads to the reference genome (version 0.7.17)  
`$ brew search bwa`  
`$ brew info bwa` _see software detail information_  
`$ brew install bwa`  
Type to check the installation  
`$ bwa`  

If installation succeeded, you will get following response  

    Program: bwa (alignment via Burrows-Wheeler transformation)
    Version: 0.7.17-r1188
    Contact: Heng Li <lh3@sanger.ac.uk>
    ...

### Install samtools
SAMtools for manipulating next-generation sequencing data (version 1.9)  
`$ brew search samtools`  
`$ brew info samtools`  
`$ brew install samtools`  
Type to check the installation  
`$ samtools --version`  

If installation succeeded, you will get following response  

    samtools 1.8
    Using htslib 1.9
    Copyright (C) 2018 Genome Research Ltd.  

### Download demo data and reference genome sequence files for first analysis
I prepared small data which aquired from public sequence database. It's already modified to contain chromosome 1 reads. File size are 54M and 55M.  
`$ wget -c https://www.dropbox.com/s/eg8k4xmmw23nfnq/DRR006760_chr1_1.fastq.gz`  
`$ wget -c https://www.dropbox.com/s/b4awju0mkt8q3bn/DRR006760_chr1_2.fastq.gz`  

You also need reference genome sequence files. Totally, 8.1G will be downloaded.  
`$ wget -c https://www.dropbox.com/s/9qmtqwgq8pxyj99/human_g1k_v37_decoy.fasta`  
`$ wget -c https://www.dropbox.com/s/9dpu7ver996c8m0/human_g1k_v37_decoy.fasta.amb`  
`$ wget -c https://www.dropbox.com/s/b3rmp79xgixiyk5/human_g1k_v37_decoy.fasta.ann`  
`$ wget -c https://www.dropbox.com/s/oplswegvl68fd96/human_g1k_v37_decoy.fasta.bwt`  
`$ wget -c https://www.dropbox.com/s/4lgsboui7l01mq1/human_g1k_v37_decoy.fasta.fai`  
`$ wget -c https://www.dropbox.com/s/6dkq2f6dokddyqs/human_g1k_v37_decoy.fasta.pac`  
`$ wget -c https://www.dropbox.com/s/4braaqyewooqt4p/human_g1k_v37_decoy.fasta.sa`  

# 4. First analysis (< 10 min.)
In this section, we will analyze following public data. 
This is very first step, so we just do align and see it using viewer.
- Title: Identification of autosomal recessive spastic paraplegia with cerebellar ataxia and neuropathy
- Abstract: Objective: To identify the gene mutation responsible for a family presenting spastic paraplegia, cerebellar ataxia and neuropathy with autosomal recessive transmission. Background: Autosomal recessive hereditary spastic paraplegias (AR-HSP) constitute a heterogeneous group of neurodegenerative diseases involving pyramidal tracts dysfunction. The genes responsible for many types of AR-HSPs remain unknown.
- Methods: The present study included two patients in a Japanese consanguineous family. Their onset of symptoms was 48 and 58 years of age. Neurologic examination and DNA analysis were underwent in two patients and two normal family members. We performed a genomewide linkage analysis employing SNP arrays with two patients’ DNAs and exome sequencing using one patient’s sample.
- Results: **We identified a homozygous missense mutation in the lysosomal trafficking regulator (LYST) gene**
- DRA: http://ddbj.nig.ac.jp/DRASearch/study?acc=DRP000999
- **Causative gene: LYST c.4189T>G, p.F1397V**
- Paper: Autosomal-recessive complicated spastic paraplegia with a novel lysosomal trafficking regulator gene mutation. - PubMed - NCBI https://www.ncbi.nlm.nih.gov/pubmed/24521565

Firstly, align paired sequence reads to the 1000 genomes project-customised human reference genome build 37 (human_g1k_v37_decoy).  
```
$ bwa mem -t4 -M \
              -R "@RG\tID:FLOWCELLID\tSM:DRR006760_chr1\tPL:illumina\tLB:DRR006760_chr1_library_1" \
              human_g1k_v37_decoy.fasta \
              DRR006760_chr1_1.fastq.gz DRR006760_chr1_2.fastq.gz | \
              samtools view -@4 -1 - | samtools sort -@4 - -o - > DRR006760_chr1.aligned_reads_sorted.bam
```  

you will get following response. It will take about few min by my MacBookPro 2014 (2.2GHz).  

    [M::bwa_idx_load_from_disk] read 0 ALT contigs
    [M::process] read 396040 sequences (40000040 bp)...
    [M::process] read 396040 sequences (40000040 bp)...
    [M::mem_pestat] # candidate unique pairs for (FF, FR, RF, RR): (0, 142256, 0, 0)
    ...snip

Then, make Index file for aligned sequence file (.bam file). It takes few seconds.  
`$ samtools index -@ 4 DRR006760_chr1.aligned_reads_sorted.bam`  

To see this aligned sequence reads,  
`$ sh sh IGV_2.4.13/igv.sh DRR006760_chr1.aligned_reads_sorted.bam`  
Go to 1:235,955,287-235,955,418 to see the mutation.  
Can you see?
![](images/IGV_LYST.png "")

## tips: using memory, avoid using disk for speeding up
We did alignment using bwa program previously. As you can see, the bwa command was long, and included samtools command.  
This is done by "|". **Pipe**. "|" connects two command. In previous command, connected bwa and samtools. The result of bwa aligned data passed to samtools. This is very important. Why? If we used "|", we can use memory space instead of writing data to the very slow hard disk. Memory speed is extrem faster than hard disk.  
Let's experience.  
```
$ bwa mem -t4 -M \
              -R "@RG\tID:FLOWCELLID\tSM:DRR006760_chr1\tPL:illumina\tLB:DRR006760_chr1_library_1" \
              human_g1k_v37_decoy.fasta \
              DRR006760_chr1_1.fastq.gz DRR006760_chr1_2.fastq.gz > DRR006760_chr1.aligned_reads.sam
$ samtools view -@4 -1 DRR006760_chr1.aligned_reads.sam -o DRR006760_chr1.aligned_reads.bam
$ samtools sort -@4  DRR006760_chr1.aligned_reads.bam -o DRR006760_chr1.aligned_reads_sorted.bam
```  



あと base_dir のことが必要
