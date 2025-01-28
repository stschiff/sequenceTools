# SequenceTools

[Install with Bioconda](https://anaconda.org/bioconda/sequencetools)

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

# PileupCaller

## Overview

The main tool in this repository is the program `pileupCaller` to sample alleles from low coverage sequence data. The first step is to generate a “pileup” file at all positions you wish to genotype. To do that, here is a typical command line, which restricts to mapping and base quality of 30 and uses a predefined set of positions to generate the pileup for (optional, see below):

    samtools mpileup -R -B -q30 -Q30 -l <list_of_positions.txt> \
        -f <reference_genome.fasta> \
        Sample1.bam Sample2.bam Sample3.bam > pileup.txt

Important Note: You should definitely use the `-B` flag, which disables base alignment quality recalibration. This mechanism is turned on by default and causes huge reference bias with low coverage ancient DNA data. This flag disables the mechanism.

In the above command line, if you use a positions-file, it should either contain positions (0-based) or a bed file (see samtools manual for details). The output is a simple text file with all positions that could be genotyped in the three samples.

Next, you need to run pileupCaller, which you run like this:

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

Note that you can also fuse the two steps above into one unix pipe:

    samtools mpileup -R -B -q30 -Q30 \
        -f <reference_genome.fasta> \
        Sample1.bam Sample2.bam Sample3.bam | \
    pileupCaller --randomHaploid --sampleNames Sample1,Sample2,Sample3 \
        --samplePopName MyPop -f <Eigenstrat.snp> \
        -e <My_output_prefix>

Here, I omitted the positions-file in the samtools command, because pileupCaller itself will ensure filtering for the positions listed in the Eigenstrat-Positions file given via option `-f`. Note that `--randomHaploid` is only one way to call genotypes. If you need stricter calling, you may want to try `--majorityCall --downSampling --minDepth 3`, which calls genotypes only on sites with at least three reads, downsamples to three if there are more, and then calls whatever of the two alleles has the majority. This will reduce errors, but also yield less data in case of lower coverage.
            
You will possibly encounter an issue: If you have aligned your read data to a version of the reference genome that uses `chr1`, `chr2` and so on as chromosome names, the resulting Eigenstrat file will be valid, but won't merge with other Eigenstrat datasets that use chromosome names `1`, `2` and so on. I would therefore recommend to strip the `chr` from your chromosome names if necessary. You can do that easily using a little UNIX filter using the `sed` tool. In the full pipeline, it looks like this:

    samtools mpileup -R -B -q30 -Q30 \
        -f <reference_genome.fasta> \
        Sample1.bam Sample2.bam Sample3.bam | sed 's/chr//' | \
    pileupCaller --sampleNames Sample1,Sample2,Sample3 \
        --samplePopName MyPop -f <Eigenstrat.snp> \
        -o EigenStrat -e <My_output_prefix>

## Options

You can see all options via `pileupCaller --help`, which outputs:

```
Usage: pileupCaller [--version] 
                    (--randomHaploid | --majorityCall [--downSampling] | 
                      --randomDiploid) [--keepIncongruentReads] 
                    [--seed <RANDOM_SEED>] [-d|--minDepth <DEPTH>] 
                    [--skipTransitions | --transitionsMissing | 
                      --singleStrandMode] (-f|--snpFile <FILE>) 
                    [(-e|--eigenstratOut <FILE_PREFIX>) [-z|----zip] | 
                      (-p|--plinkOut <FILE_PREFIX>) 
                      [--popNameAsPhenotype | --popNameAsBoth] [-z|----zip] |
                      --vcf] 
                    (--sampleNames NAME1,NAME2,... | --sampleNameFile <FILE>) 
                    [--samplePopName POP(s)]

  PileupCaller is a tool to create genotype calls from bam files using
  read-sampling methods. To use this tool, you need to convert bam files into
  the mpileup-format, specified at http://www.htslib.org/doc/samtools.html
  (under "mpileup"). The recommended command line to create a multi-sample
  mpileup file to be processed with pileupCaller is

      samtools mpileup -B -q30 -Q30 -l <BED_FILE> -R -f <FASTA_REFERENCE_FILE>
          Sample1.bam Sample2.bam Sample3.bam | pileupCaller ...

  You can lookup what these options do in the samtools documentation. Note that
  flag -B in samtools is very important to reduce reference bias in low coverage
  data.


  This tool is part of sequenceTools version 1.6.0.0

Available options:
  --version                Print version and exit
  -h,--help                Show this help text
  --randomHaploid          This method samples one read at random at each site,
                           and uses the allele on that read as the one for the
                           actual genotype. This results in a haploid call
  --majorityCall           Pick the allele supported by the most reads at a
                           site. If an equal numbers of alleles fulfil this,
                           pick one at random. This results in a haploid call.
                           See --downSampling for best practices for calling
                           rare variants
  --downSampling           When this switch is given, the MajorityCalling mode
                           will downsample from the total number of reads a
                           number of reads (without replacement) equal to the
                           --minDepth given. This mitigates reference bias in
                           the MajorityCalling model, which increases with
                           higher coverage. The recommendation for rare-allele
                           calling is --majorityCall --downsampling --minDepth 3
  --randomDiploid          Sample two reads at random (without replacement) at
                           each site and represent the individual by a diploid
                           genotype constructed from those two random picks.
                           This will always assign missing data to positions
                           where only one read is present, even if minDepth=1.
                           The main use case for this option is for estimating
                           mean heterozygosity across sites.
  --keepIncongruentReads   By default, pileupCaller now removes reads with
                           tri-allelic alleles that are neither of the two
                           alleles specified in the SNP file. To keep those
                           reads for sampling, set this flag. With this option
                           given, if the sampled read has a tri-allelic allele
                           that is neither of the two given alleles in the SNP
                           file, a missing genotype is generated. IMPORTANT
                           NOTE: The default behaviour has changed in
                           pileupCaller version 1.4.0. If you want to emulate
                           the previous behaviour, use this flag. I recommend
                           now to NOT set this flag and use the new behaviour.
  --seed <RANDOM_SEED>     random seed used for the random number generator. If
                           not given, use system clock to seed the random number
                           generator.
  -d,--minDepth <DEPTH>    specify the minimum depth for a call. For sites with
                           fewer reads than this number, declare Missing
                           (default: 1)
  --skipTransitions        skip transition SNPs entirely in the output,
                           resulting in a dataset with fewer sites.
  --transitionsMissing     mark transitions as missing in the output, but do
                           output the sites.
  --singleStrandMode       [THIS IS CURRENTLY AN EXPERIMENTAL FEATURE]. At C/T
                           polymorphisms, ignore reads aligning to the forward
                           strand. At G/A polymorphisms, ignore reads aligning
                           to the reverse strand. This should remove post-mortem
                           damage in ancient DNA libraries prepared with the
                           non-UDG single-stranded protocol.
  -f,--snpFile <FILE>      an Eigenstrat-formatted SNP list file for the
                           positions and alleles to call. All positions in the
                           SNP file will be output, adding missing data where
                           there is no data. Note that pileupCaller
                           automatically checks whether alleles in the SNP file
                           are flipped with respect to the human reference, and
                           in those cases flips the genotypes accordingly. But
                           it assumes that the strand-orientation of the SNPs
                           given in the SNP list is the one in the reference
                           genome used in the BAM file underlying the pileup
                           input. Note that both the SNP file and the incoming
                           pileup data have to be ordered by chromosome and
                           position, and this is checked. The chromosome order
                           in humans is 1-22,X,Y,MT. Chromosome can generally
                           begin with "chr". In case of non-human data with
                           different chromosome names, you should convert all
                           names to numbers. They will always considered to be
                           numerically ordered, even beyond 22. Finally, I note
                           that for internally, X is converted to 23, Y to 24
                           and MT to 90. This is the most widely used encoding
                           in Eigenstrat databases for human data, so using a
                           SNP file with that encoding will automatically be
                           correctly aligned to pileup data with actual
                           chromosome names X, Y and MT (or chrX, chrY and
                           chrMT, respectively).
  -e,--eigenstratOut <FILE_PREFIX>
                           Set Eigenstrat as output format. Specify the
                           filenames for the EigenStrat SNP, IND and GENO file
                           outputs: <FILE_PREFIX>.snp, <FILE_PREFIX>.ind and
                           <FILE_PREFIX>.geno. If not set, output will be
                           FreqSum (Default). Note that freqSum format,
                           described at
                           https://rarecoal-docs.readthedocs.io/en/latest/rarecoal-tools.html#vcf2freqsum,
                           is useful for testing your pipeline, since it's
                           output to standard out
  -z,----zip               GZip the output Eigenstrat or Plink genotype and SNP
                           files. Filenames will be appended with '.gz'. To zip
                           FreqSum or VCF output, just zip the standard output
                           of this program, for example `pileupCaller ... --vcf
                           | gzip -c > out.vcf.gz
  -p,--plinkOut <FILE_PREFIX>
                           Set Plink as output format. Specify the filenames for
                           the Plink BIM, FAM and BED file outputs:
                           <FILE_PREFIX>.bim, <FILE_PREFIX>.fam and
                           <FILE_PREFIX>.bed. If not set, output will be FreqSum
                           (Default). Note that freqSum format, described at
                           https://rarecoal-docs.readthedocs.io/en/latest/rarecoal-tools.html#vcf2freqsum,
                           is useful for testing your pipeline, since it's
                           output to standard out
  --popNameAsPhenotype     Only valid for Plink Output: Write the population
                           name into the last column of the fam file, as a
                           Phenotype according to the Plink Spec. By default,
                           the population name is specified as the first column
                           only (family name in the Plink spec)
  --popNameAsBoth          Only valid for Plink Output: Write the population
                           name into both the first and last column of the fam
                           file, so both as Family-ID and as a Phenotype
                           according to the Plink Spec. By default, the
                           population name is specified only as the first column
                           (family name in the Plink spec)
  -z,----zip               GZip the output Eigenstrat or Plink genotype and SNP
                           files. Filenames will be appended with '.gz'. To zip
                           FreqSum or VCF output, just zip the standard output
                           of this program, for example `pileupCaller ... --vcf
                           | gzip -c > out.vcf.gz
  --vcf                    output VCF format to stdout
  --sampleNames NAME1,NAME2,...
                           give the names of the samples as comma-separated list
                           (no spaces)
  --sampleNameFile <FILE>  give the names of the samples in a file with one name
                           per line
  --samplePopName POP(s)   specify the population name(s) of the samples, which
                           are included in the output *.ind.txt file in
                           Eigenstrat output. This will be ignored if the output
                           format is not Eigenstrat. If a single name is given,
                           it is applied to all samples, if multiple are given,
                           their number must match the the number of samples
                           (default: Left "Unknown")
```

### SingleStrandMode

pileupCaller supports a special calling mode (`--singleStrandMode`) for sequencing data generated from single-stranded libraries (Gansauge, Marie-Theres, and Matthias Meyer. 2013. “Single-Stranded DNA Library Preparation for the Sequencing of Ancient or Damaged DNA.” Nature Protocols 8 (4): 737–48.). The idea is that at C/T SNPs, forward mapping reads are discarded, and at G/A SNPs, reverse mapping reads are discarded. This will get rid of post-mortem ancient DNA damage in a conservative way, i.e. it will remove more than necessary and make sure that the remainder of the data is clean of DNA damage, improving the overall calling quality.

There is an important catch: If you have data from paired-end sequencing, and you are using _unmerged_ reads, then this approach will fail, as it will then _not_ discard potentially damaged reads.

So there are two options if you have Paired-end sequencing data:
1) Use only merged reads and `--singleStrandMode`
2) Use all reads but do _not_ use `--singleStrandMode`. Instead, in such cases I recommend to trim reads from both ends to remove ancient DNA damage. Depending on the details of the library construction, you may have UDG-treated data, in which case fewer basepairs would have to be trimmed.

### VCF output
VCF output was added in version 1.6.0.0. The VCF format is specified in detail at https://samtools.github.io/hts-specs/VCFv4.5.pdf. I just mention two specifics. First, with calling modes `--randomHaploid` and `--majorityCall`, the output genotypes will be haploid. This means that instead of genotypes like `0/0`, `0/1`, `1/1` or `./.`, you will instead just see `0`, `1` or `.`. Second, I added some possibly useful filters and statistics to the output, which are described in the header of the VCF:

```
##fileformat=VCFv4.2
##source=pileupCaller_v1.6.0.0
##command_line=pileupCaller --randomHaploid -f 1240k_eigenstrat_snp_short.snp.txt --sampleNames 1,2,3,4 --vcf
##group_names=Unknown,Unknown,Unknown,Unknown
##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of Samples With Data">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">
##FILTER=<ID=s50,Description="Less than 50% of samples have data">
##FILTER=<ID=s10,Description="Less than 10% of samples have data">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
##FORMAT=<ID=DP8,Number=8,Type=Integer,Description="Nr of Reads supporting A,C,G,T in forward strand, followed by the same quartet in reverse strand">
```

As you can see, Info fields NS, DP and AF are added and defined, as well as two filters which might come in handy. Also, beyond the required genotype `GT` tag, I added two per-sample tags `DP` and `DP8` as defined in the header. 

### Summary Statistics

PileupCaller automatically outputs a few lines of summary statistics, including the number of sites called for each sample, and the average read depth. These are output to the stderr, so do not affect stdout or file output.

# vcf2eigenstrat

Simple tool to convert a VCF file to an Eigenstrat file. Pretty self-explanatory. Please run `vcf2eigenstrat --help` to output some documentation.

# genoStats

A simple tool to get some per-individual statistics from an Eigenstrat or Freqsum-file. Run `genoStats --help` for documentation.

# Scripts
This package also contains several haskell wrapper scripts for the following [ADMIXTOOLS and EIGENSOFT](https://reich.hms.harvard.edu/software) commands: convertf, mergeit, qp3Pop, qpDstat and smartPCA. The original tools require parameter files as input, which I find tedious to use in bioinformatics pipelines. I wrote those wrapper scripts to be able to start the tools with a simple command line option interface.

If you have `stack` installed your system (see above), you should be able to run those scripts on your machine without any difficult setup. Simply clone this repository, navigate to the `scripts` subfolder and invoke any script using standard bash execution, for example

    ./convertf_wrapper.hs

If you start this the first time it may take a while, since `stack` downloads all dependencies and even the script interpreter for you, but after that it should start instantanious. If you want to use the scripts from your path, I suggest to put symbolic links into any folder that is already on your path (for example `~/.local/bin`).
