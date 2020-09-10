#!/usr/bin/python

"""usage: .py Orthogroups.csv Outfile.txt"""

import sys

sys.stdout = open(sys.argv[2], "w+")

with open(sys.argv[1], "r") as handle:
    handle.readline()
    for line in handle:
        line = line.replace(",","\t").strip().split('\t')
        for j in range(1, len(line)):
            for k in range(j + 1, len(line)):
                print(line[j].strip() + "\t" + line[k].strip())

