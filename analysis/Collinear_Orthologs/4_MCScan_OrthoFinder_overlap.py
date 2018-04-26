#!/usr/bin/python

"""Usage: .py all.singletons all.collinearity geneID.txt all.tandem

This script takes in a file with all singletons in it, 
the collinearity output from MCScanX, and a file with
species-specific geneIDs in them. Each line in the geneID.txt file
should be the begnning string for gene sequence names (e.g. Gorai, Gohir.A, Gohir.D)
"""

import sys


def process_group():
    global block
    global pos
    global neg
    if pos > 0 and block:
        #print(str(pos) + '\t' + str(neg))
        #if pos == 1: 
        for item in block:
        #        print(*item, sep='\t')
            connected_components(item[0], item[1])
    block = []
    pos = 0
    neg = 0

def check_singletons(geneA, geneB):
    if not any(geneA in sl for sl in singletons): #check if geneA is in singletons list
        return False
    elif not any(geneB in sl for sl in singletons): #check if geneA is in singletons list
        return False
    else:
        matched = [y for y in singletons if geneA in y]
        if geneB in matched:
            return True
        else: return False

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


seqID = []
singletons = []
components = []
block = []
pos = 0
neg = 0

with open(sys.argv[1], "r") as handle: #all.singletons 
    first_line = handle.readline()
    for line in handle: 
        line = line.strip().split("\t")
        singletons.append(line)
      


with open(sys.argv[2], "r") as handle: #all.collinearity
    for line in handle:
        line = line.strip()
#        print(line)
        if line[0:2] == "##":
            process_group()
        elif line[0] != "#":
            line = line.split('\t')
            block.append([line[1], line[2], line[3]])
#            print(*block, sep = '\t')
            singleton_match = find_group(line[1], singletons)
            if singleton_match:
                if line[2] in singleton_match:
                    pos += 1
                    block[-1].append("POS")
                else:
                    neg += 1
                    block[-1].append("NEG")
            else: 
               singleton_match = find_group(line[2], singletons)
               if singleton_match:
                   neg += 1
    process_group() #last group doesn't have '##', so explicitly call it here



with open(sys.argv[3], "r") as handle: #geneID.txt
    for line in handle: 
        line = line.strip()
        seqID.append(line)
    seqID = sorted(seqID)

with open(sys.argv[4], "r") as handle: #all.tandem
    for line in handle:
        line = line.strip().split(',')
        connected_components(line[0], line[1])


components.sort(key=len, reverse=True)

"""
with open("MCScanX_groups_posmore0.csv", "w+") as handle:
#with open(sys.argv[5], "w+") as handle: #outfile.txt
    for item in seqID:
        handle.write("\t%s" % item)
    handle.write("\n")
    for n,i in enumerate(components):
        index = '{0:07}'.format(n)
        final_string = "OG" + str(index)
        i = sorted(i)
        for j in seqID:
            j = [k for k in i if j in k]
            j_s = ','.join(j)
            final_string = final_string + "\t" + j_s
        handle.write(final_string + '\n')
"""



 
