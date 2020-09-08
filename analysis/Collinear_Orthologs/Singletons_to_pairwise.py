#!/usr/bin/python 

import sys

with open(sys.argv[1], "r") as f:
    for line in f:
        genes = line.strip().split('\t')
        for i in range(len(genes)):
            for j in range(i+1,len(genes)):
                print(genes[i] + "\t" + genes[j])
