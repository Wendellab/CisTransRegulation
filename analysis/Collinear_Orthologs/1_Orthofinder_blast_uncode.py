#!/usr/local/python

"""Usage: python .py SequenceIDs.txt <blastfile>

This script takes the blast file from OrthoFinder
and un-codes the sequence IDs."""

import sys

seqID = {}

with open(sys.argv[1], "r") as handle:
    for line in handle: 
        line = line.strip().split(": ")
        seqID[line[0]] = line[1]

with open(sys.argv[2], "r") as handle: 
    for line in handle: 
        line = line.strip().split("\t")
        line[0] = seqID[line[0]]
        line[1] = seqID[line[1]]
        print(*line, sep="\t")
