# Analysis of cis-trans regulation underlying cotton fiber domestication

## Experimental Design
Twenty-four RNA-seq libraries were generated for 10 and 20 dpa fibers from wild and domesticated *Gossypium hirsutum* accessions and their reciprocal F1 hybrids each with three replicates: 2 fiber developmental stages x 4 accessions x 3 replicates = 24 samples.   

* fiber development: 10 dpa and 20 dpa
* cotton accessions: Maxxa, TX2094, F1_MxT, F1_TxM
* biological replicates: r1, r2, r3

## RNA-seq mapping and expression estimation
RNA-seq reads were mapped against *G. hirsutum* transcriptome to estimate gene expressions in different accessions. In order to distinguish Maxxa(m) and TX2094(t) alleles in F1 hybrids, we used [HyLiTE](http://hylite.sourceforge.net) (Duchemin et al. 2015) to detect allelic SNPs followed by allele-specific read partitioning.

### 1.prepare reference transcriptome
Three reference genomes are now available for *G. hirstum* cultivar TM1 by @Li_genome_2015 [AD1_BGI](https://www.cottongen.org/species/Gossypium_hirsutum/bgi-AD1_genome_v1.0), @Zhang_tm1_2015 [AD1_NBI](https://www.cottongen.org/species/Gossypium_hirsutum/nbi-AD1_genome_v1.1) and Saski et al (2017) [AD1_458](https://phytozome.jgi.doe.gov/pz/portal.html#!info?alias=Org_Ghirsutum_er). The latest one (Saski et al, 2017) was used in this analysis.

    cd ~/jfw-lab/Projects/cis_trans_regulation/HyLiTE_newTM1
    ## locate reference transcriptome on server
    ln -s ~/jfw-lab/GenomicResources/archived_resources/AD1Saski/v1.1/annotation/Ghirsutum_458_v1.1.cds_primaryTranscriptOnly.fa.gz
    ## how many primary transcripts? 66,577
    zcat Ghirsutum_458_v1.1.cds_primaryTranscriptOnly.fa.gz | grep -c '>'
    zcat Ghirsutum_458_v1.1.cds_primaryTranscriptOnly.fa.gz >Ghirsutum_458_v1.1.cds_primaryTranscriptOnly.fasta

In addition to nuclear genes, we could use mitochondrial genes to check samples from reciprocal F1 hybrids. The mt genome was also sequenced ([Liu et al 2013](http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0069476#s2)), containing 35 protein genes. Go to [GenBank: JX944505.1](https://www.ncbi.nlm.nih.gov/nuccore/430728001), Send "Coding Sequences" to create "FASTA Nucleotide" file, save as `mt.transcript.fasta`. Note that only 34 coding regions were annotated.
   
    ## rename mt genes
    cp mt.transcript.fasta mt.transcript.rename.fasta 
    sed -i 's/lcl.*gene=/mt_/g' mt.transcript.rename.fasta 
    sed -i 's/\].*//g' mt.transcript.rename.fasta 
    
I prepared the TM1 reference transcriptome with both nuclear and mt genes.

    cat Ghirsutum_458_v1.1.cds_primaryTranscriptOnly.fasta mt.transcript.rename.fasta > TM1.nuclNmt.transcript.fasta
    ## load bowtie2 Version: 2.2.6 
    module load bowtie2
    ## build bowtie reference
    bowtie2-build TM1.nuclNmt.transcript.fasta TM1new

### 2.prepare RNA-seq fastq files
    ln -s ~/jfw-lab/Projects/cis_trans_regulation/HyLiTE/fastq
    
### 3.run Bowtie2 mapping
HyLiTE can take sequences (`.fastq`), mapping results (`.sam`), or mapping variants (`.pileup`) as input. It is more flexible to run mapping first. Bowtie2 mapping `-N 1` allows one mismatch, more sensitive but slower than default 0. Although not specified in HyLiTE pipeline, I used `--no-mixed --no-discordant --no-unal --dovetail` to disable mixed mode only looking for concordant alignment of PE reads, suppress SAM record for reads not aligned, and allow overlap containing one read. For example:

    bowtie2 -q -p 8 -t -N 1 --no-mixed --no-discordant --no-unal --dovetail -x TM1new -1 fastq/TX2094-3-10dpa-fiber_201_R1.fq -2 fastq/TX2094-3-10dpa-fiber_201_R2.fq -S mapping/TX2094-3-10dpa-fiber_201.sam 2>mapping/TX2094-3-10dpa-fiber_201.log
    samtools view -bS mapping/TX2094-3-10dpa-fiber_201.sam | samtools sort - -o mapping/TX2094-3-10dpa-fiber_201.bam ; samtools index mapping/TX2094-3-10dpa-fiber_201.bam

Once all mappings were finished, group alignment results from multiple BAM files. `-B -Q 0 -d 1000000` disables probabilistic realignment for the computation of base alignment quality (BAQ), set minimal base quality as 0, and maximal depth 1000000. Input `bam.files.txt` lists BAM files one file per line, and *the order must match the oder in the next protocol file `protocol_file_bam.txt`*.
    
    samtools faidx TM1.nuclNmt.transcript.fasta
    samtools mpileup -B -Q 0 -d 1000000 -f TM1.nuclNmt.transcript.fasta -b bam.files.txt > mapping/bam.pileup

Run bash command `runBowtie2.sh` [here](), keep all `.sam` files in `/mapping` directory.

    mkdir mapping
    bash runBowtie2.sh > runBowtie2.sh.out.txt 

Check mapping statistics
   
    bash checkMapping.sh


### 4.run HyLiTE
Create protocol file `sample_protocol_file_bam.txt` [here](), by listing all `.bam` files for parental and hybrid samples.

    module load python/2.7.12
    HyLiTE -v -f sample_protocol_file_bam.txt -p mapping/bam.pileup -n results052717 -r TM1.nuclNmt.transcript.fasta >results052717.log 2>&1
    # HyLiTE options:
    # -v turns on verbose runtime comments
    # -f specifies the protocol file
    # -r specifies the .fasta file containing the reference gene sequences
    # -n allows the user to provide a name for the HyLiTE analysis, and creates a directory for output files
    # -S use mapping .sam results instead of .fastq from protocol file
    # -b use pre built ref
    
    # If satart from SAM files
    # HyLiTE -v -S -f sample_protocol_file_sam.txt -r TM1.nuclNmt.transcript.fasta -b TM1new -n results053017 >results053017.log 2>&1


### 5.interpretation of HyLiTE output files

#### Result files for following analysis
* `results052717.expression.txt` : Read count table for agene expression from all samples, regardless of different alleles in F1 (as total expression of alleles used).
* `results052717.F1.F1-X-XXdpa.read.summary.txt` : Read count table for allelic gene expression in F1. On filer per 12 file one for each F1 sample12 F1 files (6 accession x 2 dpa).

#### Summary and other files
* `results052717.run.summary.txt` : Summarization of F1 read assignment and SNP identification.
  * Total number of child reads mapping on the reference:   150,089,104
  * Number of child read unambiguously assigned to a parent:        22,586,739
  * Number of child read unambiguously assigned to maxxa:   11,050,923
  * Number of child read unambiguously assigned to TX2094:  11,535,816
  * Number of child read with uninformative assignment:     48,516,887
  * Number of child read with unknown or ambiguous assignement:     78,985,478
  * Total number of SNPs identified:        597,788
  * Total number of child unique SNPs:      269,787
* `results052717.snp.txt` : Information of all SNPs detected, one SNP each line, accounting for 33 categories (presence and absence combination) of SNP
* `results052717.snp.summary.txt` : Categorization information of above SNPs tabulated by 66611 genes.
* `results052717.F1.F1-X-XXdpa.read.txt` : Specific information about every read in the child for a given biological replicate. One file per 12 F1 files (6 accession x 2 dpa). 
* `results052717.pickle` : used to save the memory content of any python object in order to be able to load it again later. 


## Cis-Trans analysis

#### 1. Preprocessing of hylite output files

    module load bioawk
    bioawk -c fastx '{ print $name, length($seq) }' <TM1.nuclNmt.transcript.fasta >TM1.nuclNmt.transcript.len.txt

    cd results052717
    module load R
    R BATCH CMD cistrans.pre.r &

    cd ..
    mkdir dataAnalysis
    mv TM1.nuclNmt.transcript.len.txt dataAnalysis/
    mv results052717/MaxxaAlleleReadCounts.txt dataAnalysis/
    mv results052717/TX2094AlleleReadCounts.txt dataAnalysis/
    mv results052717/totalReadCounts.txt dataAnalysis/
    mv results052717/totalReadCounts.txt dataAnalysis/
    
#### 2. Expression analysis by DESeq2
 
    cd dataAnalysis/
    R BATCH CMD cistrans.deseq2.r &

###
###


## Scripts and result files

* `runBowtie2.sh` [here]()
* `sample_protocol_file_bam.txt` [here]()
`cistrans.pre.r`
`cistrans.deseq2.r`


