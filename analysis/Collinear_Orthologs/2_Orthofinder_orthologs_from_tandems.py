#!/usr/local/python 

"""Usage: python .py <concat.tandem> Orthogroups.csv all.groups.tandem Outfile.txt 

This program takes the Orthogroups.csv results file from OrthoFinder (Emms and Kelly 2015) 
and the .tandem file from MCScanX to filter out any groups that are singleton groups, or groups
in which all paralogs are in a single tandem array with no interspacing gene(s). 

The output of this program is a pairwise list of all single-genic orthogroups, except those that are within 
a tandem array. 

As for the current writing, genes must be neighboring in order to be considered tandemly duplicated. 
This may be an area for further development."""


import sys


def add_edge_to_graph(a,b):
    for component in tandems:
        if a in component:
            if b not in component:
                component.append(b)
            break # we don't have to look for other components for a
    return tandems


def merge_graph(a,b):
    for component in tandems:
        if a in component:
            for i, other_component in enumerate(tandems):
                if b in other_component and other_component != component: # a, and b are already in different tandem array: merge
                    component.extend(other_component)
                    tandems[i:i+1] = []
                    break
            break
    return tandems


def connected_components(a,b):
    test_b = any(b in com for com in tandems)
    test_a = any(a in com for com in tandems)
    if test_a and test_b: #both are in componants aleady
        merge_graph(a,b)
    elif test_a: #if a in tandems, but b is not
        add_edge_to_graph(a,b)
    elif test_b: #if b is in tandems, but a is not
        add_edge_to_graph(b,a)
    else: #neither a nor b is found
        tandems.append([a,b])
        return(tandems)


def find_group(a, lists):
    for i in lists:
        if a in i:
            return i
    return []



tandems = []
block = []


#components.sort(key=len, reverse=True)

#for i in components:
#    i = sorted(i)
#    if len(i) == int(sys.argv[2]):
#        print(*i, sep="\t")


################################################################





with open(sys.argv[1], "r") as handle: #Concat.tandem File
    for line in handle:
        line = line.strip()
        line = line.split(',')
        connected_components(line[0], line[1])


sys.stdout = open(sys.argv[3], "w+")  #all.groups.tandem
for line in tandems:
    print(*line, sep = '\t')


sys.stdout = open(sys.argv[4], "w+") #outfile


with open(sys.argv[2], "r") as handle: #Orthogroups.csv
    first_line = handle.readline().strip() #Species Identifiers
    species_num = len(first_line.split("\t"))
    for line in handle:
        ones = [True] * species_num #If a species contains no genes in a group, the starting value is True
#        print(*ones, sep='\t')
        line = line.strip().split("\t")
        if len(line) - 1 != species_num:
            continue
        for i in range(len(line[1:]), 0, -1):
            ones[i-1] = False
#            print(*ones, sep='\t')
            genes = line[i].strip().split(",")
            if len(genes) == 1:
                ones[i-1] = True
                if line[i] == "":
                    line.pop(i)
                    ones[i-1] = False
            else: #More than one gene for a species 
                #Are all genes in tandem array list? 
                tandem_bool = True
#                print(tandem_bool)
                for gene in genes:
                    gene = gene.strip()
#                    print(gene)
                    if any(gene in com for com in tandems) != True:
                        tandem_bool = False
                    #else:
                    #    tandem_bool = False
#                print(tandem_bool)
                if tandem_bool == True:
                    ones[i-1] = True
                    for gene in genes:
                        line.append(gene.strip())
                    line.pop(i)
                else: 
                    break; #if one species has non-tandem array, move on to next Orthogroup





        if all(ones):
            line = sorted(line[1:])
            #print(*line, sep = '/t')
            for j in range(0, len(line)):
                for k in range(j + 1, len(line)):
                    print(line[j] + "\t" + line[k])



#sys.stdout = open(sys.argv[3], "w+") #outfile
#for line in tandems:
#    for j in range(0, len(line)):
#        for k in range(1, j):
#            print(line[j] + "\t" + line[k])

