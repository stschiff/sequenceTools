# SequenceTools

[![Install with Bioconda](https://anaconda.org/bioconda/sequencetools/badges/installer/conda.svg)](https://anaconda.org/bioconda/sequencetools)

This repository contains some programs that I use for processing sequencing data.

# Installation

## Installation via precompiled executables

* [Download the latest Executable](https://github.com/stschiff/sequenceTools/releases/latest) that best matches your platform.

For example, to install `pileupCaller` in Linux, you can run the following commands to get started:

```bash
# download the current stable release binary
wget https://github.com/stschiff/sequenceTools/releases/latest/download/pileupCaller-linux
# make it executable
chmod +x pileupCaller-linux
# run it
./trident-linux -h
```

## Installation from source

1. Clone the repository: `git clone https://github.com/stschiff/sequenceTools.git`
2. Go into the repository: `cd sequenceTools`
3. Compile the executables in the repository using `stack`: `stack install`

This last step will take a while, as it not only compiles the source, but also first downloads the correct version of the Haskell compiler.

# Commands

## pileupCaller

The main tool in this repository is the program `pileupCaller` to sample alleles from low coverage sequence data. The first step is to generate a “pileup” file at all positions you wish to genotype. To do that, here is a typical command line, which restricts to mapping and base quality of 30:

    samtools mpileup -R -B -q30 -Q30 -l <list_of_positions.txt> \
        -f <reference_genome.fasta> \
        Sample1.bam Sample2.bam Sample3.bam > pileup.txt

Important Note: You should definitely use the `-B` flag, which disables base alignment quality recalibration. This mechanism is turned on by default and causes huge reference bias with low coverage ancient DNA data. This flag disables the mechanism.

In the above command line, the file "list_of_positions.txt" should either contain positions (0-based) or a bed file (see samtools manual for details). The output is a simple text file with all positions that could be genotyped in the three samples.

Next, you need to run my tool pileupCaller, which you run like this:

    pileupCaller --randomHaploid --sampleNames Sample1,Sample2,Sample3 \
        --samplePopName MyPop -f <Eigenstrat.snp> \
        -e <My_output_prefix> < pileup.txt

Here, options `--sampleNames` gives the names of the samples that is output in the Eigenstrat `*.ind` file, and and `-–samplePopName` is optional to also give the population names in that file (defaults to `Unknown`, you can also change it later in the output). Then, (required) option `-f` needs an Eigenstrat positions file. This is required for pileupCaller to know what is the reference and which the alternative allele in your reference dataset that you want to call. An Eigenstrat positions file is a line-based file format, where each line denotes a SNP position, and there are exactly six required columns, denoting in order i) SNP ID, ii) chromosome, iii) genetic position (can be set to zero), iv) physical position, v) reference allele, vi) alternate allele. Here is an example:

    rs0000  11  0.000000    0   A   C
    rs1111  11  0.001000    100000  A   G
    rs2222  11  0.002000    200000  A   T
    rs3333  11  0.003000    300000  C   A
    rs4444  11  0.004000    400000  G   A
    rs5555  11  0.005000    500000  T   A
    rs6666  11  0.006000    600000  G   T

Finally, the `-e` option specifies Eigenstrat as output format and gives the prefix for the `*.ind`, `*.pos` and `*.geno` files. Without the `-e` option, pileupCaller will output in FreqSum format,  described [here](https://rarecoal-docs.readthedocs.io/en/latest/rarecoal-tools.html#vcf2freqsum), which is useful for debugging your pipeline, since it's just a single file that is output into the terminal and can therefore easily be inspected.

You can also get some help by typing `pileupCaller -h`, which shows a lot more option, for example the sampling method, minimal coverage and other important options.

Note that you can also fuse the two steps above into one unix pipe:

    samtools mpileup -R -B -q30 -Q30 -l <list_of_positions.txt> \
        -f <reference_genome.fasta> \
        Sample1.bam Sample2.bam Sample3.bam | \
    pileupCaller --randomHaploid --sampleNames Sample1,Sample2,Sample3 \
        --samplePopName MyPop -f <Eigenstrat.snp> \
        -e <My_output_prefix>

Note that `--randomHaploid` is only one way to call genotypes. If you need stricter calling, you may want to try `--majorityCall --downSampling --minDepth 3`, which calls genotypes only on sites with at least three reads, downsamples to three if there are more, and then calls whatever of the two alleles has the majority. This will reduce errors, but also yield less data in case of lower coverage.
            
            
There is however an issue here: If you have aligned your read data to a version of the reference genome that uses `chr1`, `chr2` and so on as chromosome names, the resulting Eigenstrat file will be valid, but won't merge with other Eigenstrat datasets that use chromosome names `1`, `2` and so on. I would therefore recommend to strip the `chr` from your chromosome names if necessary. You can do that easily using a little UNIX filter using the `sed` tool. In the full pipeline, it looks like this:

    samtools mpileup -R -B -q30 -Q30 -l <list_of_positions.txt> \
        -f <reference_genome.fasta> \
        Sample1.bam Sample2.bam Sample3.bam | sed 's/chr//' | \
    pileupCaller --sampleNames Sample1,Sample2,Sample3 \
        --samplePopName MyPop -f <Eigenstrat.snp> \
        -o EigenStrat -e <My_output_prefix>

    
Note: You do not have to use a positions file in your `samtools` step. You can also just generate pileup-data for every covered position (default without passing a positions file via `-l`) and have pileupCaller filter the sites for you. This makes sense for dense genotyping, but a positions file might speed up the process for sparser genotyping.

### SingleStrandMode

pileupCaller supports a special calling mode (`--singleStrandMode`) for sequencing data generated from single-stranded libraries (Gansauge, Marie-Theres, and Matthias Meyer. 2013. “Single-Stranded DNA Library Preparation for the Sequencing of Ancient or Damaged DNA.” Nature Protocols 8 (4): 737–48.). The idea is that at C/T SNPs, forward mapping reads are discarded, and at G/A SNPs, reverse mapping reads are discarded. This will get rid of post-mortem ancient DNA damage in a conservative way, i.e. it will remove more than necessary and make sure that the remainder of the data is clean of DNA damage, improving the overall calling quality.

There is an important catch: If you have data from paired-end sequencing, and you are using _unmerged_ reads, then this approach will fail, as it will then _not_ discard potentially damaged reads.

So there are two options if you have Paired-end sequencing data:
1) Use only merged reads and `--singleStrandMode`
2) Use all reads but do _not_ use `--singleStrandMode`. Instead, in such cases I recommend to trim reads from both ends to remove ancient DNA damage. Depending on the details of the library construction, you may have UDG-treated data, in which case fewer basepairs would have to be trimmed.


## vcf2eigenstrat

Simple tool to convert a VCF file to an Eigenstrat file. Pretty self-explanatory. Please run `vcf2eigenstrat --help` to output some documentation.

## genoStats

A simple tool to get some per-individual statistics from an Eigenstrat or Freqsum-file. Run `genoStats --help` for documentation.

## Scripts
This package also contains several haskell wrapper scripts for the following [ADMIXTOOLS and EIGENSOFT](https://reich.hms.harvard.edu/software) commands: convertf, mergeit, qp3Pop, qpDstat and smartPCA. The original tools require parameter files as input, which I find tedious to use in bioinformatics pipelines. I wrote those wrapper scripts to be able to start the tools with a simple command line option interface.

If you have `stack` installed your system (see above), you should be able to run those scripts on your machine without any difficult setup. Simply clone this repository, navigate to the `scripts` subfolder and invoke any script using standard bash execution, for example

    ./convertf_wrapper.hs

If you start this the first time it may take a while, since `stack` downloads all dependencies and even the script interpreter for you, but after that it should start instantanious. If you want to use the scripts from your path, I suggest to put symbolic links into any folder that is already on your path (for example `~/.local/bin`).
