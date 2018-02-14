#!/usr/local/python

import sys



with open(sys.argv[1], "r") as handle: 
    first_line = handle.readline()
    first_line = first_line.strip().split('\t')
    for line in handle:
        if "," not in line:
            line = line.strip().split('\t')[1:]
            good_line = True
            for i in range(len(first_line)-1):
                if line[i] == "" or line[i] == " ":
                    good_line = False
            if len(line) == len(first_line) and good_line:
                print(*line, sep = '\t')
