#!/usr/bin/env python3

import argparse
import math

parser = argparse.ArgumentParser(description="Calculate cost-efficient additional sequencing from Preseq output. Currently for Capture data only")
parser.add_argument('-e', "--endogenous", metavar="<ENDOGENOUS_FRAC>", type=float, help="Fraction of reads (0-1 scale) that map to the reference", required=True)
parser.add_argument('-m', "--mapped_reads", metavar="<MAPPED_READS>", type=int, help="Number of mapped reads after duplicate removal", required=True)
parser.add_argument('-c', "--coverage", metavar="<COVERAGE>", help="Mean coverage on target SNPs", type=float, required=True)
parser.add_argument('-s', "--snps", metavar="<NR_SNPS>", help="Number of SNPs in the capture panel [Default=1240000]", type=int, default=1240000)
parser.add_argument('-g', "--goal", metavar="<GOAL>", help="Goal, given as the number of new reads you are willing to produce to get one more SNP covered. [Default=100]", type=float, default=100)
parser.add_argument("extrap_file", metavar="<FILE>", type=open, help="Preseq extrap file")

args = parser.parse_args()

coverage_per_read = args.coverage / args.mapped_reads

covered_snps = None
tot_sequenced = None

print("sequenced_reads", "SNPS covered", "Cost", sep="\t")

next(args.extrap_file)
for line in args.extrap_file:
    fields = line.strip().split()
    exp_tot_mapped = float(fields[0])
    exp_unique = float(fields[1])
    exp_tot_sequenced = exp_tot_mapped / args.endogenous
    exp_coverage = exp_unique * coverage_per_read
    exp_snps_covered = (1.0 - math.exp(-exp_coverage)) * args.snps
    if covered_snps is not None:
        new_snps = exp_snps_covered - covered_snps
        new_tot_sequenced = exp_tot_sequenced - tot_sequenced
        cost = new_tot_sequenced / new_snps
        print(int(exp_tot_sequenced), int(exp_snps_covered), int(cost), sep="\t")
        if cost > args.goal:
            break
    covered_snps = exp_snps_covered
    tot_sequenced = exp_tot_sequenced
