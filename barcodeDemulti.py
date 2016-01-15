#!/usr/bin/env python3

import argparse
import gzip
from itertools import islice, chain
import sys
import os

class FastqEntry:
    def __init__(self, name, seq, quali):
        self.name = name
        self.seq = seq
        self.quali = quali

def iterFastQ(fn):
    f = gzip.open(fn, "rt")
    for name, seq, dummy, quali in grouper(f, 4):
        yield FastqEntry(name, seq, quali)

def grouper(iterable, n):
    args = [iter(iterable)] * n
    return zip(*args)

def writeFastq(f, fastq):
    lines = [fastq.name, fastq.seq, "+", fastq.quali]
    f.writelines(lines)

def readSampleSheet(sample_sheet_file):
    ret = {}
    f = open(sample_sheet_file, "r")
    for line in f:
        [name, i1, i2, bc1, bc2] = line.strip().split()
        for letter in i1 + i2 + bc1 + bc2:
            assert (letter in "ACTGN-"), "illegal Index or Barcode sequence. Must contain only A, C, T, G, N or -"
        ret[name] = (i1, i2, bc1, bc2)
    return ret

def generateOutputFileDict(barcodeDict, output_dir, lane_spec):
    ret = {}
    for name in chain(barcodeDict.keys(), ["undetermined"]):
        os.makedirs("{}/{}".format(output_dir, name), exist_ok=True)
        fn1 = "{o}/{n}/{n}_R1_{l}.fastq.gz".format(o=output_dir, n=name, l=lane_spec)
        fn2 = "{o}/{n}/{n}_R2_{l}.fastq.gz".format(o=output_dir, n=name, l=lane_spec)
        print("generated output file {}".format(fn1), file=sys.stderr)
        print("generated output file {}".format(fn2), file=sys.stderr)
        f1 = gzip.open(fn1, "wt")
        f2 = gzip.open(fn2, "wt")
        ret[name] = (f1, f2)
    return ret

def closeHandles(fileDict):
    for f1, f2 in fileDict.values():
        f1.close()
        f2.close()

def imatch(seq, tag, mismatches):
    if tag == '-':
        return True
    nrM = sum(t != 'N' and s != t for s, t in zip(seq, tag))
    return nrM <= mismatches

def getIndexDict(sampleDict):
    ret = {}
    for n, indices in sampleDict.items():
        ret[indices] = n
    return ret

perfektMatch = 0
imperfectMatch = 0
def findMatch(seqs, indexDict, m):
    global perfektMatch
    global imperfectMatch
    if seqs in indexDict:
        perfektMatch += 1
        return indexDict[seq]
    elif m > 0:
        for i, n in indexDict.items():
            if imatch(seqs[0], i[0], m) and imatch(seqs[1], i[1], m) and imatch(seqs[2], i[2], m) and imatch(seqs[3], i[3], m):
                imperfectMatch += 1
                return n
        return None
    else:
        return None
            
parser = argparse.ArgumentParser(description="Demultiplex Fastq files and trim barcodes")
parser.add_argument('-s', "--sample_sheet", metavar="<SAMPLE_SHEET>", help="The file containing the sample information. Tab separated file, where every sample gets a row. The columsn are: 1) Sample_name; 2) Index 1 (P5); 3) Index2 (P7); 4) Internal Barcode Read 1; 5) Internal Barcode Read2. You can use the special symbol '-' to denote the absence of a barcode.", required=True)
parser.add_argument('-m', "--mismatches", metavar="<mismatches>", type=int, default=1, help="the number of mismatches allowed in the index- and barcode-recognition. Default=1")
parser.add_argument("-o", "--output_dir", metavar="<OUTPUT_DIRECTORY>", required=True, help="directory to write the output folders into")
parser.add_argument("-l", "--lane_spec", metavar="<lane_spec>", required=True, help="a tag to specify the lane, e.g. 001 for lane 1. Output files will be constructed as <OUTPUT_DIRECTORY/<NAME>/<NAME>_R<READ>_<lane_spec>.fastq.gz. So for example for read 1 and lane_spec 001: .../<NAME>_R1_001.fastq.gz.")
parser.add_argument("--fastqR1", required=True, metavar="<FASTQ_Read1>", help="fastq file for read1")
parser.add_argument("--fastqR2", required=True, metavar="<FASTQ_Read2>", help="fastq file for read2")
parser.add_argument("--fastqI1", required=True, metavar="<FASTQ_IndexRead1>", help="fastq file for index read1")
parser.add_argument("--fastqI2", required=True, metavar="<FASTQ_IndexRead2>", help="fastq file for index read2")

args = parser.parse_args()
sampleDict = readSampleSheet(args.sample_sheet)
indexDict = getIndexDict(sampleDict)
print("read sample sheet with {} samples".format(len(sampleDict)), file=sys.stderr)
output_handles = generateOutputFileDict(sampleDict, args.output_dir, args.lane_spec)

line_nr = 1
for (r1, r2, i1, i2) in zip (iterFastQ(args.fastqR1), iterFastQ(args.fastqR2), iterFastQ(args.fastqI1), iterFastQ(args.fastqI2)):
    if line_nr % 10000 == 0:
        print("processing line {}".format(line_nr), file=sys.stderr)
    m = args.mismatches
    name = findMatch((i1.seq, i2.seq, r1.seq, r2.seq), indexDict, m)
    if name is not None:
        o1, o2 = output_handles[name]
        queryBC1, queryBC2 = sampleDict[name][2:]
        clip1 = 0 if queryBC1 == '-' else len(queryBC1)
        clip2 = 0 if queryBC2 == '-' else len(queryBC2)
        fastqEntry1 = FastqEntry(r1.name, r1.seq[clip1:], r1.quali[clip1:])
        fastqEntry2 = FastqEntry(r2.name, r2.seq[clip2:], r2.quali[clip2:])
        writeFastq(o1, fastqEntry1)
        writeFastq(o2, fastqEntry2)
    else:
        o1, o2 = output_handles["undetermined"]
        fastqEntry1 = FastqEntry(r1.name, r1.seq, r1.quali)
        fastqEntry2 = FastqEntry(r2.name, r2.seq, r2.quali)
        writeFastq(o1, fastqEntry1)
        writeFastq(o2, fastqEntry2)
    line_nr += 1
    if line_nr > 30000:
        break

print(perfektMatch, imperfectMatch)
closeHandles(output_handles)
