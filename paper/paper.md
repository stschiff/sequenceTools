---
title: 'PileupCaller: A command-line tool to sample genotypes from low-coverage sequencing data of ancient DNA'
tags:
  - Ancient DNA
  - Genome Sequencing
  - Bioinformatics
authors:
  - name: Stephan Schiffels
    orcid: 0000-0002-1017-9150
    affiliation: 1
affiliations:
 - name: Max Planck Institute for Evolutionary Anthropology, Leipzig, Germany
   index: 1
date: 14 November 2024
bibliography: paper.bib
---

# Summary

Next generation sequencing data is ubiquitous in medical and biological sciences. It has also become the primary tool in archaeogenetics, where ancient DNA is extracted from archaeological organic (often human skeletal) material, processed into DNA sequencing libraries and then sequenced [@Orlando2021]. As a testimony to the rapid and accellerating growth of the field, we today have close to ten thousand published ancient human genomes available in the public record [@Schmid2024; @Mallick2024], and many smaller datasets of other organisms. A key step in processing raw sequencing data is the estimation of genotypes at specific variable positions along the genome. Such positions are often pre-selected because they are informative about ancestry or of particular biological relevance [@Haak2015; @Mathieson2015; @Rohland2022]. While established tools exist for this task for high-quality modern sequencing data [@samtools; @gatk], these are often not appropriate for ancient DNA, which has often too low sequencing-coverage and a higher error rate due to post-mortem DNA damage. PileupCaller is a command-line tool written in Haskell, which randomly samples genotypes from raw alignment data at predefined bi-allelic positions. Several modes can be selected, geared towards specific input data features and research questions.

# Statement of need

Present-day DNA, for example from medical studies results in raw sequencing data with relatively low per-base error rates and sequencing-coverages of at least several multiples of 1 (for example [@1000_Genomes_Project_Consortium2015]) but in fact up to 20-30x coverage. Dedicated tools to process such data include samtools/bcftools [@samtools] and GATK [@gatk] among many other tools. Ancient DNA seuqencing data often comes with substantially lower coverage and substantially higher error rates. In terms of coverage, most ancient genomes have genome-wide coverage often below 1x and in fact very often even below 0.1x. Such low coverage means that any given genomic site is more likely not covered by a sequencing read than covered. At the same time, the low fraction of sites that is actually covered has higher error rates than modern DNA, due to ancient-DNA damage. These two factors violate the assumptions behind statistical genotype callers like `bcftools call` or `HaplotypeCaller` from GATK. 

As is widely used practice in the field, very low-coverage ancient DNA data is often "called", simply by randomly selecting reads at a given position of interest. PileupCaller is a command-line tool that does exactly that, by reading in a list of SNP positions and a stream of sequencing data, some optional filtering options, and then performs random samples at every position of interest for multiple individuals. Even before this paper, `pileupCaller` has been widely used since its creation in 2017, mostly because of its simple use and low-memory footprint thanks to streaming.

# Usage and key functionality

`pileupCaller` relies on the `pileup` format defined in [@samtools]. This format lists for each site the nucleotides from all reads covering that site, including base-qualities. An example pileup-file can be found in the [repository](). These files are rarely saved, but used as an intermediate format for streaming, specifically from `samtools mpileup`(). A typical usage command line for `pileupCaller` could be:

```bash
    samtools mpileup -R -B -q30 -Q30 -f <reference_genome.fasta> \
        Sample1.bam Sample2.bam Sample3.bam | \
    pileupCaller --randomHaploid --sampleNames Sample1,Sample2,Sample3 \
        --samplePopName MyPop -f <Eigenstrat.snp> \
        -e <My_output_prefix>
```

Among the options passed to `samtools mpileup`, we highlight the `-B` flag, which switches off base-alignment quality recalibration and is critical for avoiding reference bias with low-coverage ancient DNA. A key input ingredient is the SNP list, which needs to be given in Eigenstrat-format [@Price2006], defined in [@eigensoft_github]. Example files for Pileup and Eigenstrat-formatted data can be found under `test/testDat` in the software repository. The SNP file not only lists the positions of the variants that should be called, but also the two possible alleles at each site. This is important, as pileupCaller will ensure that only those two alleles are called. Sites at which a third allele is called will be output as missing.

In terms of output formats, pileupCaller currently supports Eigenstrat, Plink (https://www.cog-genomics.org/plink/2.0/) and VCF [@VCF], with an option to additionally compress the output in gzip-format. Command line options are documented inline via `pileupCaller --help`.

PileupCaller relies on a key Haskell dependency [`sequence-formats`](https://hackage.haskell.org/package/sequence-formats), which was developed alongside the `sequenceTools` package to which pileupCaller belongs. 

# References