#!/usr/local/python 

"""Usage: python .py <concat.gff> Orthogroups.csv 

This program takes the Orthogroups.csv results file from OrthoFinder (Emms and Kelly 2015) 
and reformats each orthogroup such that species are tab-delimited, tandemly duplicated arrays 
are surrounded by parentheses, and genes within each species and tandem array is comma-separted. 

As for the current writing, genes must be neighboring in order to be considered tandemly duplicated. 
This may be an area for further development."""


import sys

Gene_Orders = {}

with open(sys.argv[1], "r") as handle: #Bed File
    i = 0
    for line in handle:
        line = line.strip().split("\t")
        chromosome = int(line[0][2:]) * 1000000
        Gene_Orders[line[1]] = chromosome+i
        i += 1
#    print(str(float(i)) + " Genes in total. If more than 1 million, please reconsider.")

with open(sys.argv[2], "r") as handle: #Orthogroups.csv
    first_line = handle.readline().strip() #Species Identifiers
    species_num = len(first_line.split("\t"))
    print('\t' + first_line)
    for line in handle:
        line = line.strip().split("\t")
        while len(line) <= species_num:
            line.append(" ")
        lens = [] #keep track of total tandem "lengths"
        for i in range(len(line[1:]), 0, -1): #running through each species, starting at the last one
            genes = line[i].strip().split(",") 
            if len(genes) == 1: 
                if line[i] == "" or line[i] == " ": 
                    lens.append("N/A")
                else: 
                    lens.append(0)
            else: 
                #Need to find Tandem Duplicates and reinsert in list
                bed_values = []
                for gene in genes:
                    bed_values.append(Gene_Orders[gene.strip()])
                lens.append((max(bed_values) - min(bed_values)))

        lens.append(line[0])
        lens = lens[::-1]
        print(*lens, sep='\t')

