#!/usr/bin/env python3

import argparse
import gzip
from itertools import islice, chain
import sys
import os

class FastqEntry:
    def __init__(self, lines):
        [name, seq, _, quali] = lines
        self.name = name
        self.seq = seq
        self.quali = quali
    
class FastQIterator:
    def __init__(self, fn):
        self.f = gzip.open(fn, "rt")
    
    def __iter__(self):
        return self
    
    def next(self):
        [name, seq, dummy, quali] = [self.f.readline() for _ in range(4)]
        return FastqEntry(name, seq, quali)
    
    def __del__(self):
        self.f.close()

class Index:
    def __init__(self, seq, mismatches):
        self.sequence = seq
        self.mismatches = mismatches
        self.hits = {}
    
    def match(self, querySeq):
        if self.sequence == '-':
            return True
        nrM = sum(t != 'N' and s != t for s, t in zip(querySeq, self.sequence))
        return nrM <= mismatches
    
    def recordHit(self, seq):
        clippedSeq = seq[:len(self.sequence)]
        if clippedSeq not in self.hits:
            self.hits[clippedSeq] = 0
        self.hits[clippedSeq] += 1

class QuadrupelIndex:
    def __init__(self, i1, i2, bc1, bc2, mismatches):
        self.i1 = Index(i1, mismatches)
        self.i2 = Index(i2, mismatches)
        self.bc1 = Index(bc1, mismatches)
        self.bc2 = Index(bc2, mismatches)
    
    def match(self, r1, r2, i1, i2):
        if self.i1.match(i1) and self.i2.match(i2) and self.bc1.match(r1) and self.bc2.match(r2):
            self.i1.recordHit(i1)
            self.i2.recordHit(i2)
            self.bc1.recordHit(r1)
            self.bc2.recordHit(r2)
            return True
        else:
            return False
        
class SampleSheet:
    def __init__(self, sample_sheet_file, outputTemplate, mismatches):
        self.samples = {}
        f = open(sample_sheet_file, "r")
        self.outputHandles = {}
        for line in f:
            [name, i1, i2, bc1, bc2] = line.strip().split()
            for letter in i1 + i2 + bc1 + bc2:
                assert (letter in "ACTGN-"), "illegal Index or Barcode sequence. Must contain only A, C, T, G, N or -"
            self.samples[name] = QuadrupelIndex(i1, i2, bc1, bc2, mismatches)
        for name in chain(self.samples.keys(), ["Undetermined"]):
            outputNameR1 = outputTemplate.replace("%name%", name).replace("%r%", "1")
            outputNameR2 = outputTemplate.replace("%name%", name).replace("%r%", "2")
            outputDir = os.path.dirname(outputNameR1)
            os.makedirs(outputDir, exist_ok=True)
            f1 = gzip.open(outputNameR1, "wt")
            f2 = gzip.open(outputNameR2, "wt")
            print("generated output file {}".format(outputNameR1), file=sys.stderr)
            print("generated output file {}".format(outputNameR2), file=sys.stderr)
            self.outputHandles[name] = (f1, f2)

    def __del__(self):
        for f1, f2 in self.outputHandles.values():
            f1.close()
            f2.close()
    
    def process(self, fastqR1, fastqR2, fastqI1, fastqI2):
        for name, qIndex in self.samples.items():
            if qIndex.match(fastqR1.seq, fastqR2.seq, fastqI1.seq, fastqI2.seq):
                queryBC1, queryBC2 = qIndex.bc1.sequence, qIndex.bc2.sequence
                clip1 = 0 if queryBC1 == '-' else len(queryBC1)
                clip2 = 0 if queryBC2 == '-' else len(queryBC2)
                fastqR1.seq = fastqR1.seq[clip1:]
                fastqR2.seq = fastqR2.seq[clip2:]
                self.writeFastq(name, (fastqR1, fastqR2))
                break
            else:
                self.writeFastq("Undetermined", (fastqR1, fastqR2))
                
    def writeFastQ(self, name, fastqPair):
        for i in [0, 1]:
            lines = [fastqPair[i].name, fastqPair[i].seq, "+", fastqPair[i].quali]
            self.outputHandles[name][i].writelines(lines) 
    
    def reportStats(self):
        for name in chain(self.samples.keys(), ["Undetermined"]):
            print(name)
            qi = self.samples[name]
            print("", "Total Reads:", sum(qi.i1.hits.values()))
            print("", "Perfect matches I1:", qi.i1.hits[qi.i1.sequence])
            print("", "Other matches I1:", [(k, v) for k, v in qi.i1.hits.items() if k != qi.i1.sequence])
            print("", "Perfect matches I2:", qi.i2.hits[qi.i2.sequence])
            print("", "Other matches I2:", [(k, v) for k, v in qi.i2.hits.items() if k != qi.i2.sequence])
            print("", "Perfect matches BC1:", qi.bc1.hits[qi.bc1.sequence])
            print("", "Other matches BC1:", [(k, v) for k, v in qi.bc1.hits.items() if k != qi.bc1.sequence])
            print("", "Perfect matches BC2:", qi.bc2.hits[qi.bc2.sequence])
            print("", "Other matches BC2:", [(k, v) for k, v in qi.bc2.hits.items() if k != qi.bc2.sequence])

def buildArgumentParser():   
    parser = argparse.ArgumentParser(description="Demultiplex Fastq files and trim barcodes")
    parser.add_argument('-s', "--sample_sheet", metavar="<SAMPLE_SHEET>", help="The file containing the sample information. Tab separated file, where every sample gets a row. The columsn are: 1) Sample_name; 2) Index 1 (P5); 3) Index2 (P7); 4) Internal Barcode Read 1; 5) Internal Barcode Read2. You can use the special symbol '-' to denote the absence of a barcode.", required=True)
    parser.add_argument('-m', "--mismatches", metavar="<mismatches>", type=int, default=1, help="the number of mismatches allowed in the index- and barcode-recognition. Default=1")
    parser.add_argument("-o", "--output_template", metavar="<OUTPUT_DIRECTORY>", required=True, help="template for files to write. Can be arbitrary paths including magic place holders %name%, which will be replaced by the sample name, and %r%, which will be replaced by the read (1 or 2)")
    parser.add_argument("--fastqR1", required=True, metavar="<FASTQ_Read1>", help="fastq file for read1")
    parser.add_argument("--fastqR2", required=True, metavar="<FASTQ_Read2>", help="fastq file for read2")
    parser.add_argument("--fastqI1", required=True, metavar="<FASTQ_IndexRead1>", help="fastq file for index read1")
    parser.add_argument("--fastqI2", required=True, metavar="<FASTQ_IndexRead2>", help="fastq file for index read2")
    return parser

if __name__ == "__main__":
    parser = buildArgumentParser()
    args = parser.parse_args()
    sampleSheet = SampleSheet(args.sample_sheet, args.output_template, args.mismatches)

    nr = 1
    for (r1, r2, i1, i2) in zip (FastQIterator(args.fastqR1), FastQIterator(args.fastqR2), FastQIterator(args.fastqI1), FastQIterator(args.fastqI2)):
        if nr % 10000 == 0:
            print("processing line {}".format(nr), file=sys.stderr)
        sampleSheet.process(r1, r2, i1, i2)
        nr += 1
        if nr > 30000:
            break
    
    