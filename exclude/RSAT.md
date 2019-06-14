# RSAT analysis

This protocol presents a step-by-step instruction for the task of discovering and annotating DNA motifs in clusters of ustream sequences for cotton using the [RSAT::plants](http://rsat.eead.csic.es/plants/) web server. Because the cotton genome used is not supported, input promoter sequence files need to be prepared ahead.  

## Input file preparation

In Maxxa versus TX2094 fibers, RD (regulatory divergent) genes exhibiting *trans* variants (*trans* only and *cis*&*trans*) showed higher coexpression network density than *cis* only and non-RD genes. My goal is to first discover the DNA motifs enriched in the promoter of *trans* affected gene, and then compare the discovered motifs with collections of known motifs in order to identify transcription factors that may be associated with the discovered motifs.

Input files include

1. 2kb promoter of 27816 genes surveyed for RD patterns as control: [maxxa], [tx2094], [TM1]
2. 2kb promoter of trans affected genes: [maxxa], [tx2094], [TM1]

### R srcipt

    # get promoter Grange list
    library(GenomicFeatures)
    txdb <- loadDb("/work/LAS/jfw-lab/hugj2006/cottonLeaf/refGenomes/txdb.TM1saski.sqlite")
    p2k = promoters(txdb, upstream=2000, downstream=0)
    rtracklayer::export.bed(p2k, "TM1saski.promoter2k.bed")
    
### Unix    
    module load bedtools2
    getfasta -fi /work/LAS/jfw-lab/hugj2006/cottonLeaf/refGenomes/TM1new_26.fasta -bed TM1saski.promoter2k.bed -name+ >TM1.promoter2k.fasta


### R srcipt

    library("Biostrings")
    # import TM1 ref promoter
    ref=readDNAStringSet("TM1")
    names(ref)=gsub("::.*","",names(ref))
    # get 2kb promoter of 27816 genes surveyed for RD patterns
    load("RDgenes.rdata")
    ctl = seq[names(ref)%in%rd$ID]
    writeXStringSet(ctl, "TM1_promoter_2kb_27816.fasta", format="fasta")
     # get 2kb promoter of trans affected genes
    trans = ref[names(ref)%in%rd$ID[grep("trans",rd$rd)]]
    writeXStringSet(trans, "TM1_promoter_2kb_trans_affected.fasta", format="fasta")
    
    
    ## import maxxa promoter seqs by Daojun Yuan
    A <- readDNAStringSet("/lss/research/jfw-lab/Projects/AD1_domestication_cistrans/promoterSNP_djyuan/Maxxa.At.Promoter.fa")
    D <- readDNAStringSet("/lss/research/jfw-lab/Projects/AD1_domestication_cistrans/promoterSNP_djyuan/Maxxa.Dt.Promoter.fa")
    # fix IDs
    seq=c(A,D)
    names(seq)=gsub("-.*",".1",gsub("Maxxa","Gohir",names(seq)))
    # get 2kb promoter of 27816 genes surveyed for RD patterns
    load("RDgenes.rdata")
    ctl = seq[names(seq)%in%rd$ID]
    writeXStringSet(ctl, "maxxa_promoter_2kb_27816.fasta", format="fasta")
    # get 2kb promoter of trans affected genes
    trans = seq[names(seq)%in%rd$ID[grep("trans",rd$rd)]]
    writeXStringSet(trans, "maxxa_promoter_2kb_trans_affected.fasta", format="fasta")
      

## Step-by-step instruction

### Analyzing Upstream Sequences of Co-regulated Genes
1. Open a connection to the RSAT::Plants server. It can bereached at [http://plants.rsat.eu]() and also at [http://fl oresta.eead.csic.es/rsat]().
2. In the left menu of the RSAT server, click on “NGS - ChIP-seq” and next "peak-motifs" to open the pipeline form.
3. Add a title for this job, such as "trans affected". Upload [maxxa_promoter_2kb_trans.fasta]() as **Peak sequences**, and upload [maxxa_promoter_2kb_27816.fasta]() as **Control sequences**.
4. Click on "Reduce peak sequences" and leave both fields blank ("number of top sequences to retain" and "cut peak sequences") to avoid having the sequences clipped.
5. Click on "Motif discovery parameters". Select two algorithms: "Discover over-represented words" [oligo - analysis] and “Discover over-represented spaced word pairs" [dyad - analysis]. Uncheck the program position - analysis.6. Click on “Compare discovered motifs with databases” andselect appropriate databases which will be used to annotate any found motifs: "footprintDB-plants","RSAT non-redundant plants","JASPAR core nonredundant plants".
7. Click on “Reporting Options”. Set “Origin” to “end” and“Offset” to 0.
8. Select output type (display or email) and press “GO”.
9. Once the job is complete click on the link [Download all results (peak -motifs_archive.zip)] to save the results on your computer, and click on the link [Download all matrices (transfac format)] to download the full set of discovered motifs to a local file named "cluster.motifs.tf". This file contains all motifs in the form of position- weight matrices (PWMs) in TRANSFAC format.

### Negative Control using Random Groups of Genes

1. Repeat above steps 10 times, by using random promoters as input.
2. The top motifs found by oligo-analysis anddyad-analysis are used as negative control, in comparison with motifs found from *trans* affected genes.

### Validating Motifs by Scanning Promoter Sequences

1. On the left-side menu select "Build control sets→permute-matrix".6. Upload "cluster.motifs.tf" obtained from de novo discovery, and press "GO". Save the results fi le as " cluster.motifs.perm1.tf". This list of shuffled motifs will serve negataive control7. Select “Pattern matching →matrix scan (full options)”.8. In the sequence box paste the contents of [maxxa_promoter_2kb_trans.fasta]() and optionally, [tx2094_promoter_2kb_trans.fasta]() to assess motif
conservation.9. Upload file "cluster.motifs.tf" and select “ TRANSFAC ”format.10. In the “Background model” section select Markov order 2 and choose “Organism-specific" with "Gossypium_raimondii.Graimondii2_0.40". Press “GO”.11. Save the "Scan result" file as "cluster.scan.tf" and press the "feature map" button to draw a map of the matched motif instances.12. Repeat steps 6 – 11 using the set of permuted PWMs "cluster.motifs.perm1.tf" and save the results as "cluster.perm1.scan.tf".

### Interpretation of Results

---

Maxxa vs TX2094: co-expression module member genes as input to detect over-representation of motifs as (Contreras-Moreira et al. 2016; Yu et al. 2015), compare these motifs to known plant CREs, study variants between the two.

* Collecting the full set of promoters for the genome of interest (control) and the promoter set  of co-expressed genes (test)
* Running “peak-motifs” with test and control for de novo motif discovery, and comparing discovered motifs with databases
* Negative control with random groups of genes (same gene number as test set; 50 random groups would work) by repeating Step 2
* Validating motifs by scanning promoter sequences - a. Check occurrence of discovered motifs along input test promoters in contrary to shuffled motifs (“permute-matrix”); or ortholog sequences; b. Check conservation of discovered motifs along ortholog sequences 

## Reference
* 