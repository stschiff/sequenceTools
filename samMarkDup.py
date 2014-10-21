#!/usr/bin/env python

import sys
import string
import 

parser = argparse.ArgumentParser()
parser.add_argument("--remove", action="store_true", help="Remove duplicates instead of marking")
args = parser.parse_args()

lastRead = None
lastReadCount = 0
dupCounts = {}
for line in sys.stdin:
    if line[0] == '@':
        print line,
        continue
    fields = string.split(string.strip(line))
    flag = int(fields[1])
    
    if flag & 0x4: #unmapped
        print line,
        continue
    
    rname = fields[2]
    pos = int(fields[3])
    seq = fields[9]
    
    if lastRead == (rname, pos, seq):
        fields[1] = str(flag | 0x400)
        lastReadCount += 1
    else:
        if lastReadCount > 0:
            if lastReadCount not in dupCounts:
                dupCounts[lastReadCount] = 1
            else:
                dupCounts[lastReadCount] += 1
        fields[1] = str(flag & (~0x400))
        lastRead = (rname, pos, seq)
        lastReadCount = 1
    
    if not args.remove:
        print "\t".join(fields)

sys.stderr.write("Cardinality\tCounts\n")
for c in sorted(dupCounts.keys()):
    sys.stderr.write("{}\t{}\n".format(c, dupCounts[c]))

