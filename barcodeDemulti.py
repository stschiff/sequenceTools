#!/usr/bin/env python3

import argparse
import gzip
from itertools import islice, chain
import sys


################# COMMAND LINE ARGUMENT PARSING ###############
parser = argparse.ArgumentParser(description="Demultiplex Fastq files and trim barcodes")
parser.add_argument('-b', "--barcode_file", metavar="<BARCODE_FILE>", help="the file containing the barcodes. Tab separated file, where every sample gets a row. The first column contains the sample name, the second the barcode of Read 1, the third the barcode of Read 2", required=True)
parser.add_argument('-m', "--mismatches", metavar="<mismatches>", default=1, help="the number of mismatches allowed in the index- and barcode-recognition. Default=1")
parser.add_argument("-o", "--output_dir", metavar="<OUTPUT_DIRECTORY>", required=True, help="directory to write the output folders into")
parser.add_argument("-l", "--lane_spec", metavar="<lane_spec>", required=True, help="a tag to specify the lane, e.g. 001 for lane 1. Output files will be constructed as <OUTPUT_DIRECTORY/<NAME>/<NAME>_R<READ>_<lane_spec>.fastq.gz. So for example for read 1 and lane_spec 001: .../<NAME>_R1_001.fastq.gz.")
parser.add_argument("--fastqR1", required=True, metavar="<FASTQ_Read1>", help="fastq file for read1")
parser.add_argument("--fastqR2", required=True, metavar="<FASTQ_Read2>", help="fastq file for read2")
parser.add_argument("--fastqI1", required=True, metavar="<FASTQ_IndexRead1>", help="fastq file for index read1")
parser.add_argument("--fastqI2", required=True, metavar="<FASTQ_IndexRead2>", help="fastq file for index read2")


####################### DEFS #####################

class FastqEntry:
    def __init__(self, name, seq, quali):
        self.name = name
        self.seq = seq
        self.quali = quali

def iterFastQ(fn):
    f = gzip.open(fn, "rt")
    for [name, seq, dummy, quali] in islice(f, 4):
        yield FastqEntry(name, seq, quali)

def writeFastq(f, fastq):
    lines = [fastq.name, fastq.seq, "+", fastq.quali]
    f.writelines(lines)

def readBarcodeDict(barcode_file):
    ret = {}
    f = open(barcode_file, "r")
    for line in barcode_file:
        name, i1, i2, bc1, bc2 = line.strip().split()
        for letter in i1 + i2 + bc1 + bc2:
            assert letter in "ACTGN"
        ret[name] = (i1, i2, bc1, bc2)
    return ret

def generateOutputFileDict(barcodeDict, output_dir, lane_spec):
    ret = {}
    for name in chain(barcodeDict.keys(), ["undetermined"]):
        fn1 = "{o}/{n}/{n}_R1_{l}.fastq.gz".format(o=output_dir, n=name, l=lane_spec)
        fn2 = "{o}/{n}/{n}_R2_{l}.fastq.gz".format(o=output_dir, n=name, l=lane_spec)
        f1 = gzip.open(fn1, "wt")
        f2 = gzip.open(fn2, "wt")
        ret[name] = (f1, f2)
    return ret

def closeHandles(fileDict):
    for f1, f2 in fileDict.values():
        f1.close()
        f2.close()

def imatch(seq, tag, mismatches):
    nrM = sum(s != t for s, t in zip(seq, tag))
    return nrM <= mismatches

####################### MAIN ######################
args = parser.parse_args()
barcodeDict = readBarcodeDict(args.barcode_file)
output_handles = generateOutputFileDict(barcodeDict, args.output_dir, args.lane_spec)

line_nr = 1
for (r1, r2, i1, i2) in zip (iterFastQ(args.fastqR1), iterFastQ(args.fastqR2), iterFastQ(args.fastqI1), iterFastQ(args.fastqI2)):
    if line_nr % 100000 == 0:
        print("processing line {}".format(line_nr), file=sys.stderr)
    m = args.mismatches
    for name, (queryI1, queryI2, queryBC1, queryBC2) in barcodeDict.items():
        hit = imatch(i1.seq, queryI1, m) and imatch(i2.seq, queryI2, m) and imatch(r1.seq, queryBC1, m) and imatch(r2.seq, queryBC2, m)
        if hit:
            o1, o2 = output_handles[name]
            fastqEntry1 = FastqEntry(r1.name, r1.seq[len(queryBC1):], r1.quali[len(queryBC1):])
            fastqEntry2 = FastqEntry(r2.name, r2.seq[len(queryBC2):], r2.quali[len(queryBC2):])
            writeFastq(o1, fastqEntry1)
            writeFastq(o2, fastqEntry2)
            break
    else:
            o1, o2 = output_handles["undetermined"]
            fastqEntry1 = FastqEntry(r1.name, r1.seq, r1.quali)
            fastqEntry2 = FastqEntry(r2.name, r2.seq, r2.quali)
            writeFastq(o1, fastqEntry1)
            writeFastq(o2, fastqEntry2)
    line_nr += 1

closeHandles(output_handles)
