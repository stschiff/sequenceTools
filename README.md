# SequenceTools

[![Anaconda package](https://anaconda.org/bioconda/sequencetools/badges/installer/conda.svg)](https://anaconda.org/bioconda/sequencetools)

(bioconda package available thanks to [apeltzer](https://github.com/apeltzer)!)

This repository contains some programs that I use for processing sequencing data.

# Simple Installation of the main tools

1. Download stack (https://docs.haskellstack.org/en/stable/README/#how-to-install<Paste>).
2. Run `stack install sequenceTools --resolver nightly`. You should now have the executables from this package under `~/.local/bin`.
3. Add `~/.local/bin` to your PATH, for example by adding to your `~/.profile` or `~/.bash_profile` the line `PATH=$PATH:$HOME/.local/bin`. Run `source ~/.profile` or `source ~/.bash_profile`, respectively, to update your path.
4. Run `pileupCaller --version`. It should output `1.4.0`. You're all set.

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

Here, options `--sampleNames` gives the names of the samples that is output in the Eigenstrat `*.ind` file, and and `-–samplePopName` is optional to also give the population names in that file (defaults to `Unknown`, you can also change it later in the output). Then, option `-f` gives an Eigenstrat positions file. This is important because the pileup file only contains sites which could be called in at least one of your samples. In order to later merge your dataset with another Eigenstrat file, pileupCaller will check every position in the other Eigenstrat file to make sure every position is output with the correct alleles and missing genotypes if appropriate. Finally, the `-e` option specifies Eigenstrat as output format and gives the prefix for the `*.ind`, `*.pos` and `*.geno` files. Without the `-e` option, pileupCaller will output in FreqSum format,  described [here](https://rarecoal-docs.readthedocs.io/en/latest/rarecoal-tools.html#vcf2freqsum), which is useful for debugging your pipeline, since it's just a single file that is output into the terminal and can therefore easily be inspected.

You can also get some help by typing `pileupCaller -h`, which shows a lot more option, for example the sampling method, minimal coverage and other important options.

Note that you can also fuse the two steps above into one unix pipe:

    samtools mpileup -R -B -q30 -Q30 -l <list_of_positions.txt> \
        -f <reference_genome.fasta> \
        Sample1.bam Sample2.bam Sample3.bam | \
    pileupCaller --randomHaploid --sampleNames Sample1,Sample2,Sample3 \
        --samplePopName MyPop -f <Eigenstrat.snp> \
        -e <My_output_prefix>

There is however an issue here: If you have aligned your read data to a version of the reference genome that uses `chr1`, `chr2` and so on as chromosome names, the resulting Eigenstrat file will be valid, but won't merge with other Eigenstrat datasets that use chromosome names `1`, `2` and so on. I would therefore recommend to strip the `chr` from your chromosome names if necessary. You can do that easily using a little UNIX filter using the `sed` tool. In the full pipeline, it looks like this:

    samtools mpileup -R -B -q30 -Q30 -l <list_of_positions.txt> \
        -f <reference_genome.fasta> \
        Sample1.bam Sample2.bam Sample3.bam | sed 's/chr//' | \
    pileupCaller --sampleNames Sample1,Sample2,Sample3 \
        --samplePopName MyPop -f <Eigenstrat.snp> \
        -o EigenStrat -e <My_output_prefix>

## vcf2eigenstrat

Simple tool to convert a VCF file to an Eigenstrat file. Pretty self-explanatory. Please run `vcf2eigenstrat --help` to output some documentation.

## genoStats

A simple tool to get some per-individual statistics from an Eigenstrat or Freqsum-file. Run `genoStats --help` for documentation.
