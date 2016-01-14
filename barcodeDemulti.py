#!/usr/bin/env python3

import argparse
import gzip

parser = argparse.ArgumentParser(description="Demultiplex Fastq files and trim barcodes")
parser.add_argument('-b', "--barcode_file", metavar="<BARCODE_FILE>", help="the file containing the barcodes. Tab separated file, where every sample gets a row. The first column contains the sample name, the second the barcode of Read 1, the third the barcode of Read 2", required=True)
parser.add_argument('-m', "--mismatches", metavar="<mismatches>", default=1, help="the number of mismatches allowed in the index- and barcode-recognition. Default=1")
parser.add_argument("-o", "--output_dir", metavar="<OUTPUT_DIRECTORY>", required=True, help="directory to write the output folders into")
parser.add_argument("-l", "--lane_spec", metavar="<lane_spec>", required=True, help="a tag to specify the lane, e.g. 001 for lane 1. Output files will be constructed as <OUTPUT_DIRECTORY/<NAME>/<NAME>_R<READ>_<lane_spec>.fastq.gz. So for example for read 1 and lane_spec 001: .../<NAME>_R1_001.fastq.gz.")
parser.add_argument("--fastqR1", required=True, metavar="<FASTQ_Read1>", help="fastq file for read1")
parser.add_argument("--fastqR2", required=True, metavar="<FASTQ_Read2>", help="fastq file for read2")
parser.add_argument("--fastqI1", required=True, metavar="<FASTQ_IndexRead1>", help="fastq file for index read1")
parser.add_argument("--fastqI2", required=True, metavar="<FASTQ_IndexRead2>", help="fastq file for index read2")

args = parser.parse_args()
barcodeDict = readBarcodeDict(args.barcode_file)
output_handles = generateOutputFileDict(barcodeDict, args.output_dir, args.lane_spec)



closeHandles(output_handles)


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
    for name in barcodeDict.keys():
        fn1 = "{o}/{n}/{n}_R1_{l}.fastq.gz".format(o=output_dir, n=name, l=lane_spec)
        fn2 = "{o}/{n}/{n}_R2_{l}.fastq.gz".format(o=output_dir, n=name, l=lane_spec)
        f1 = gzip.open(fn1, "w")
        f2 = gzip.open(fn2, "w")
        ret[name] = (f1, f2)

def closeHandles(fileDict):
    for f1, f2 in fileDict.values():
        f1.close()
        f2.close()
