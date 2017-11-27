# Analysis of cis-trans regulation underlying cotton fiber domestication

## Experimental Design
Twenty-four RNA-seq libraries were generated for 10 and 20 dpa fibers from wild and domesticated *Gossypium hirsutum* accessions and their reciprocal F1 hybrids each with three replicates: 2 fiber developmental stages x 4 accessions x 3 replicates (4 reps for Maxxa) = 26 samples.   

* fiber development: 10 dpa and 20 dpa
* cotton accessions: Maxxa, TX2094, F1_MxT, F1_TxM
* biological replicates: r1, r2, r3

## RNA-seq mapping and expression estimation
RNA-seq reads were mapped against *G. hirsutum* transcriptome to estimate gene expressions in different accessions. The latest genome sequence and annotation was used as mapping reference ([Saski et al. 2017](https://phytozome.jgi.doe.gov/pz/portal.html#!info?alias=Org_Ghirsutum_er)). In order to distinguish Maxxa(m) and TX2094(t) alleles in F1 hybrids, we used [HyLiTE](http://hylite.sourceforge.net) (Duchemin et al. 2015) to detect allelic SNPs followed by allele-specific read partitioning.



---
previous results, will be updated
* [fiberAnalysis.md](fiberAnalysis.md) - Detailed documentation of data analysis from RNA-seq to cis-trans results
* [mappingResult.md](mappingResult.md) - Summary of RNA-seq library and mapping statistics
