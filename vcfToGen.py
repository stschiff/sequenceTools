#!/usr/bin/env python3

import sys

def PLtoP(PL):
    min_PL = min(PL)
    pVec = []
    for pl in PL:
        p = 10 ** (-(pl - min_PL) / 10.0)
        if p < 0.0001:
            p = 0.0
        pVec.append(p)
    norm = sum(pVec)
    return list(map(lambda p:p/norm, pVec))

line_nr = 1
for line in sys.stdin:
    if line[0] != '#':
        fields = line.strip().split()
        chr_ = fields[0]
        pos = fields[1]
        ref = fields[3]
        alt = fields[4]
        id_ = "SNP{}".format(line_nr)
        rs = "{}:{}".format(chr_, pos)


        pVecs = []
        for gen in fields[9:]:
            PL = list(map(int, gen[4:].split(",")))
            pVec = PLtoP(PL) if PL != [0, 0, 0] else PL
            pVecs.append(" ".join(map(str, pVec)))
        print(id_, rs, pos, ref, alt, " ".join(pVecs))
        line_nr += 1
