#!/usr/bin/python

"""This takes in as input the OrthoGroups.csv from OrthoFinder and calculates how many genes from each species are in each group. This is output in a tab-delimited table.  

usage: .py OrthoGroups.csv <outfile>"""

import sys

infile = open(sys.argv[1], "r")
outfile = open(sys.argv[2], "w")

sys.stdout = outfile

with infile as f:
    first_line = f.readline().strip()
    print(first_line)
    first_line = first_line.split('\t')
    gain_loss = {key:{"gain":[], "loss":[]} for key in first_line[1:]}
    for line in f:
        line = line.strip().split('\t')[1:]
        while len(line) < len(first_line):
            line.append(" ")
        c = [len(i.split(',')) for i in line]
        for i in range(len(first_line)):
            if line[i] == "" or line[i] == " ":
                c[i] = 0 #corrects for no genes reporting has having length of 1 (empty string)
        print(*c, sep='\t')



