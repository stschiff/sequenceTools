#!/usr/bin/env python3

import argparse

parser = argparse.ArgumentParser(description="Calculate cost-efficient additional sequencing from Preseq output. Currently for Capture data only")
parser.add_argument('-e', "--endogenous", metavar="<ENDOGENOUS_FRAC>", type=float, help="Fraction of reads (0-1 scale) that map to the reference", required=True)
parser.add_argument('-m', "--mapped_reads", metavar="<MAPPED_READS>", type=int, help="Number of mapped reads after duplicate removal", required=True)
parser.add_argument('-c', "--coverage", metavar="<COVERAGE>", help="Mean coverage on target SNPs", type=float, required=True)
parser.add_argument('-s', "--snps", metavar="<NR_SNPS>", help="Number of SNPs in the capture panel", type=int, required=True)
parser.add_argument('-g', "--goal", metavar="<GOAL>", help="Goal, given as the number of new reads you are willing to produce to get one more SNP covered. [Default=100]", default=100)
parser.add_argument("extrap_file", metavar="<FILE>", type=open, help="Preseq extrap file")

args = parser.parse_args()

