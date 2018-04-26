#!/usr/bin/python

"""This scripts takes in the all.groups.tandem output from script 2_ and the all.collinearity from MCScanX_h
and clusters each collinear block into a genome-wide set of orthogroups. This also will add any tandem duplicates 
to increase the total number of orthogroups detected. 

This program outputs one file: 
all.singletons #file in which only single copy genes are detected by MCScanX, unless that gene is part of a tandem array.

Usage: .py all.collinearity all.tandem <number of genomes in analysis>"""

import sys


def process_group():
    global block
    if block: #if string has values in it (not valid for header lines in file)
        for item in block:
            connected_components(item[0], item[1])
    block = []


def add_edge_to_graph(a,b):
    for component in components:
        if a in component:
            if b not in component:
                component.append(b)
            break # we don't have to look for other components for a
    return components


def merge_graph(a,b):
    for component in components:
        if a in component:
            for i, other_component in enumerate(components):
                if b in other_component and other_component != component: # a, and b are already in different components: merge
                    component.extend(other_component)
                    components[i:i+1] = []
                    break
            break
    return components


def connected_components(a,b):
    test_b = any(b in com for com in components)
    test_a = any(a in com for com in components)
    if test_a and test_b: #both are in componants aleady
        merge_graph(a,b)
    elif test_a: #if a in components, but b is not
        add_edge_to_graph(a,b)
    elif test_b: #if b is in components, but a is not
        add_edge_to_graph(b,a)
    else: #neither a nor b is found
        components.append([a,b])
        return(components)


def find_group(a, lists):
    for i in lists:
        if a in i:
            return i
    return []



components = []
block = []
tandem_arrays = []

with open(sys.argv[1], "r") as handle: 
    for line in handle:
        line = line.strip()
        tandem_arrays.append(line)
        line = line.split('\t')
        components.append(line)

tandem_arrays = [i.split('\t') for i in tandem_arrays]

with open(sys.argv[2], "r") as handle: #DDtAt.collinearity from MCScanX_h 
    for line in handle:
        line = line.strip()
        if line[0:2] == "##":
            process_group()
        elif line[0] != "#":
            line = line.split('\t')
            block.append([line[1], line[2]])
    process_group() #puts final collinear block through algorithm            

components.sort(key=len, reverse=True)

sys.stdout = open(sys.argv[4], "w+")

for i in components:
    i = sorted(i)
    tandems = []
    for j in range(0, len(i)):
        j_group = find_group(i[j], tandem_arrays)
        for k in range(j + 1, len(i)):
            if i[k] in j_group and i[k] not in tandems:
                tandems.append(i[k])
    #print(str(len(tandems) + int(sys.argv[3])) + "\t" + str(len(i))) 
    if len(tandems) + int(sys.argv[3]) == len(i):
        print(*i, sep = '\t')
            
        
