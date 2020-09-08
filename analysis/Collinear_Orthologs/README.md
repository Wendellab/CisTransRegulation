# Creating a Genome-Wide Dataset of Collinear Orthologs

This pipeline is designed to take proteins sequences from several species and create a list of orthologs in which 

## Step 1: Run Orthofinder on all proteins



```
module load python/3.4.3 
orthofinder -f Collinear_Orthologs -t 9 -a 9 -S diamond
```

```
cd Collinear_Orthologs/Results_Feb14 && mkdir MCScanX
cd MCScanX && mkdir MCScanX_h
cd MCScanX_h
head -1 ../../Orthogroups.csv > Singletons.txt && awk '(NR!=1 && !/^\t|\t\t|\t.$/) {for (i=1; i<NF; i++) printf $i "\t"; print $NF}' ../../Orthogroups.csv | grep -v "," >> Singletons.txt
python /work/LAS/jfw-lab/jconover/scripts/6_Orthofinder_to_all_pairwise.py Singletons.txt DDtAt.homology
```


## Make .gff file for MCScan input    
`cp ../../../DDtAt.gff > ./DDtAt_h.gff`

```
module load mscanx/2017/-4-3
MCScanX_h MCScanX_h/DDtAt_h
#python MCScanh_OrthoFinder_singletons.py MCScanX_h/DDtAt_h.collinearity 3 > DDtAt_h.singletons
```
#### Prepare blast files for MCScanX
Since Orthofinder is also based on sequential all-vs.-all blast, we will just modify those files to fix the OrthoFinder-specific formatting gene names.    

```    
gunzip -c ../WorkingDirectory/Blast*.gz >> blast_orthoformat.txt
python OrthoFinder_blast_uncode.py ../WorkingDirectory/SequenceIDs.txt blast_orthoformat.txt > DDtAt.blast  
rm blast_orthoformat.txt  
```

#### Make gff file for MCScanX

#### Run MCScanX on full dataset with null arguments
```MCScanX -b 2 ./DDtAt
```
