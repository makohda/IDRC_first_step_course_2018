# Aim
This mini-course is made for a lecture which let you gett a first tiny step to bioinformatic analysis.  
It is a step by step hands-on style.

This course will cover:
- customise your macintosh environment for running bioinformatic programs
- experience simple target resequence data analysis
- understand what steps are there (alignment, remove duplicate, variant call..etc)
- how interpretate variants

Not cover:
- how algorithmn works
- how to handle massive amount data using cluster computers
- deep/heuristic/complicated issues/knowledge/techniques/pitfalls (e.g. how I sense danger data signs in heystacks)

# Set up your mac

You have to set up your macintosh environment for informatic analyses.  
I know this is a first barricade to step in learning informatic skills, but this is the easiest way.

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

## Install homebrew, a nice packge manager for macOS
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

cat command is concatenate. Concatenate multiple file, like this  
`$ cat fileA fileB`  
Let's test.  
`$ wget https://www.dropbox.com/s/5yfaiolgoi3cp3u/test_variant_data_02.tsv`  
`$ cat test_variant_data_01.tsv test_variant_data_02.tsv`  
To make it easier to see,  
`$ cat test_variant_data_01.tsv test_variant_data_02.tsv | less -S`  
_Push Q key for quit_  
-S is a option of less command. It chnage less behavior to chop-long-lines.  
Commands have their specific options. You can see like this.  
`$ less --help`  
`$ cat --help`


## Install softwares required for (minimal) WES analysis
At first, type this, to cover scientific programs well  
`$ brew tap brewsci/bio`

### Install bwa
bwa for aligning reads to the reference genome (version 0.7.17)  
`$ brew search bwa`  
`$ brew info bwa` _see software detail information_  
`$ brew install bwa`

### Install samtools
SAMtools for manipulating next-generation sequencing data (version 1.9)  
`$ brew search samtools`  
`$ brew info samtools`  
`$ brew install samtools`  


