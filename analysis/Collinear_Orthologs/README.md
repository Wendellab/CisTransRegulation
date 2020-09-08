# Creating a Genome-Wide Dataset of Collinear Orthologs

This pipeline is designed to take proteins sequences from several species and create a list of orthologs in which 

## Step 1: Run Orthofinder on all proteins



```
module load python/3.4.3 
orthofinder -f Collinear_Orthologs -t 9 -a 9 -S diamond
```


## Step 2: Prepare files for MCScanX_h
#### Make .homology file
```
cd Collinear_Orthologs/Results_Feb14 && mkdir MCScanX
cd MCScanX && mkdir MCScanX_h
cd MCScanX_h
head -1 ../../Orthogroups.csv > Singletons.txt && awk '(NR!=1 && !/^\t|\t\t|\t.$/) {for (i=1; i<NF; i++) printf $i "\t"; print $NF}' ../../Orthogroups.csv | grep -v "," >> Singletons.txt
python /work/LAS/jfw-lab/jconover/scripts/6_Orthofinder_to_all_pairwise.py Singletons.txt DDtAt.homology
```


#### Make gff files for MCScanX and MCScanX_h
This assumes that you already have a single gff file formatted for MCScanX named DDtAt.gff      
Please consult the [MCScanX manual](http://chibba.pgml.uga.edu/mcscan2/documentation/manual.pdf) for further instruction  
`cp ../../../DDtAt.gff > ./DDtAt_h.gff`    
`cp ../../../DDtAt.gff > ../DDtAt.gff`

## Step 3; Run MCScanX_h and extract singletons
```
module load mscanx/2017/-4-3
MCScanX_h ./DDtAt_h
python 3_Orthofinder_singletons_from_MCScanXh.py MCScanX_h/DDtAt_h.collinearity 3 > DDtAt_h.singletons
```

## Step 4: Prepare blast files for MCScanX
Since Orthofinder is also based on sequential all-vs.-all blast, we will just modify those files to fix the OrthoFinder-specific formatting gene names.    

```    
gunzip -c ../WorkingDirectory/Blast*.gz >> blast_orthoformat.txt
python 1_Orthofinder_blast_uncode.py ../WorkingDirectory/SequenceIDs.txt blast_orthoformat.txt > DDtAt.blast  
rm blast_orthoformat.txt  
```

## Step 5: Run MCScanX on full dataset with null arguments
```
MCScanX -b 2 ./DDtAt
```

## Step 6: Find all singleton genes from OrthoFinder (including those in tandem arrays)
DDtAt.tandem is the output from MCScanX       
Orthogroups.csv is the output from OrthoFinder      
all.groups.tandem is an output file (not required for input) that will contain all tandemly duplicated genes in a network format (similar to OrthoFinder Output)     
OrthoFinder_singeltons_w_tandems.txt is an output file that will have high-confidence Orthologs from OrthoFinder. This file will determine which syntenic blocks will be used in step 7    

`python 2_Orthofinder_orthologs_from_tandems.py DDtAt.tandem Orthogroups.csv all.groups.tandem OrthoFinder_singletons_w_tandems.txt `


## Step 7: Find Orthologs based on Synteny in MCScanX, bassed on singletons found in OrthoFinder 
OrthoFinder_singletons_w_tandems.txt is the output from step 6      
DDtAt.collinearity is the output from MCScanX     
Each line in the geneID.txt file should be the begnning string for gene sequence names (e.g. Gorai, Gohir.A, Gohir.D)       
DDtAt.tandem is the output from MCScanX

`python 4_MCScan_OrthoFinder_overlap.py OrthoFinder_singletons_w_tandems.txt DDtAt.collinearity geneID.txt DDtAt.tandem`


## Step 8 and 9: Check tandem Duplication Size and orthologous group sizes
This isn't necessary, but can used used as a check to see how good the pipeline did. Ideally, there should be small numbers in both files.         
Larger numbers in the group sizes file indicated large gene groups based on tandem duplications (or errors in synteny analyses)      
Larger numbers in the tandem array file indicated errors in the synteny analysis. This could happen for a number of reason, especially given genomes with a high number of polyploid events in its evolutionary past. 

```
python 5_Check_Tandems.py DDtAt.gff Orthogroups.csv #This is a naive script, and assumes the genes in the gff files are listed in sequential order as they appear in the genome
python GroupSizeTable_Orthofinder.py OrthoGroups.csv <outfile> #This simply converts any Orthogroups.csv formatted file and returns a file with orthogroup size (same format as the similar file from OrthoFinder) 
```
