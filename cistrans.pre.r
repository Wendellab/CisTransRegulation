# Analysis of fiber gene cis-trans effects
# Preprocessing of hylite output files
# Guanjing Hu

# setwd to the directory of hylite output files

# module load bioawk
# bioawk -c fastx '{ print $name, length($seq) }' <TM1.nuclNmt.transcript.fasta >TM1.nuclNmt.transcript.len.txt

#####################################################
############### SNP categorization ##################
#####################################################
x<-read.table(grep("snp.txt",list.files(),value=TRUE),header=TRUE,sep="\t")
dim(x) # total snps 632908
type <- as.data.frame(ftable( xtabs(~F1+Maxxa+TX2094, data=x) ))
type<-type[type$Freq>0,]
dim(type)  #33 categories
type$note <- ''
type$note[type$F1==-1 | type$Maxxa ==-1 | type$TX2094 == -1] <- "poor coverage"
type$note[ type$Maxxa == type$TX2094 & type$Maxxa==1 & type$F1==1 ] <- "common"
type$note[ type$Maxxa == type$TX2094 & type$Maxxa==1 & type$F1=="1,0" ] <- "common"
type$note[ type$Maxxa == type$TX2094 & type$Maxxa==1 & (type$F1==0 |type$F1=="0,0" ) ] <- "parents specific"
type$note[ type$Maxxa==1 & type$TX2094 ==0 & (type$F1==0 |type$F1=="0,0") ] <- "Maxxa unique1"
type$note[ type$Maxxa==0 & type$TX2094 ==1 & (type$F1==0 |type$F1=="0,0") ] <- "TX2094 unique1"
type$note[ type$Maxxa == type$TX2094 & type$Maxxa==0 ] <- "F1 unique1"
type$note[ type$Maxxa ==1 & type$TX2094==0 & (type$F1=="1,0" )] <- "allele diagnostic"
type$note[ type$Maxxa ==0 & type$TX2094==1 & (type$F1=="1,0" )] <- "allele diagnostic"
type$note[ type$Maxxa ==1 & type$TX2094==0 & (type$F1==1)] <- "allele diagnostic: xaxxa allele in F1"
type$note[ type$Maxxa ==0 & type$TX2094==1 & (type$F1==1)] <- "allele diagnostic: tx2094 allele in F1"

bg<-aggregate(type$Freq,list(type$note),sum)
total<-sum(type$Freq)#632908
bg$percentage <- bg$x/total
bg
#                                  Group.1      x   percentage
# 1                      allele diagnostic  90614 0.1431708874
# 2 allele diagnostic: tx2094 allele in F1    305 0.0004819026
# 3  allele diagnostic: xaxxa allele in F1     98 0.0001548408
# 4                                 common  50876 0.0803845109
# 5                             F1 unique1 258868 0.4090136323
# 6                          Maxxa unique1   8014 0.0126621879
# 7                       parents specific    796 0.0012576867
# 8                          poor coverage 216929 0.3427496571
# 9                         TX2094 unique1   6408 0.0101246943
#### >40% SNPs unqiue to F1
#### over 30% SNPs are uninformative due to poor read coverage
type<-type[order(type$note),]
write.table(type, "SNPtypes.txt",row.names=FALSE,sep="\t")
# what genes contain disgnostic SNPs
x$note <- ''
x$note[x$F1==-1 | x$Maxxa ==-1 | x$TX2094 == -1] <- "poor coverage"
x$note[ x$Maxxa == x$TX2094 & x$Maxxa==1 & x$F1==1 ] <- "common"
x$note[ x$Maxxa == x$TX2094 & x$Maxxa==1 & x$F1=="1,0" ] <- "common"
x$note[ x$Maxxa == x$TX2094 & x$Maxxa==1 & (x$F1==0 |x$F1=="0,0" ) ] <- "parents specific"
x$note[ x$Maxxa==1 & x$TX2094 ==0 & (x$F1==0 |x$F1=="0,0") ] <- "Maxxa unique1"
x$note[ x$Maxxa==0 & x$TX2094 ==1 & (x$F1==0 |x$F1=="0,0") ] <- "TX2094 unique1"
x$note[ x$Maxxa == x$TX2094 & x$Maxxa==0 ] <- "F1 unique1"
x$note[ x$Maxxa ==1 & x$TX2094==0 & (x$F1=="1,0" | x$F1==1)] <- "allele diagnostic: Maxxa"
x$note[ x$Maxxa ==0 & x$TX2094==1 & (x$F1=="1,0" | x$F1==1)] <- "allele diag0550nostic: tx2094"

length(unique(x$GENE))
# No. of genes contain SNPs 50500
select <- grep("diagnostic",x$note)
length(unique(x$GENE[select]))
# No. of genes contain diagnostic SNPs 27820

#####################################################
################# SNPs on genes #####################
#####################################################

y<-read.table(grep("snp.summary.txt",list.files(),value=TRUE),header=TRUE,sep="\t")
nrow(y) # total No. of genes 66611
y$snps <- rowSums(y[,-1])
table(y$snp>0) #50501 contain SNPs
table(y$F1.TX2094>0 | y$F1.Maxxa>0)  # 27820 contains SNPs diagonostic to alleles in F1
genes <- as.character( y$GENE[y$F1.TX2094>0 | y$F1.Maxxa>0] )
length(genes) #27820

#####################################################
######### Allelic specific expression ###############
#####################################################

fileL<-grep(".read.summary.txt", list.files(), value=TRUE)
# keep allele-specific counts to 4 data frames
alleleM  <- data.frame(gene=y$GENE)    # Maxxa
alleleMn <- data.frame(gene=y$GENE)    # Maxxa + Maxxa.N; N refers to F1 unique SNPs
alleleT  <- data.frame(gene=y$GENE)    # TX2094
alleleTn <- data.frame(gene=y$GENE)    # TX2094
for(file in fileL){
    ac <- read.table(file,header=TRUE,sep="\t")
    tag<-gsub("-",".",gsub(".*F1[.]|.read.summary.txt","",file))
    alleleM[,tag]   <- ac$Maxxa
    alleleMn[,tag]  <- ac$Maxxa + ac$Maxxa.N
    alleleT[,tag]   <- ac$TX2094
    alleleTn[,tag]  <- ac$TX2094 + ac$TX2094.N
}
sum(alleleM[,-1]) # No. of reads containing only Maxxa diagnostic SNPs 9552667
sum(alleleMn[,-1]) # No. of reads containing Maxxa and F1 diagnostic SNPs ** 11332177
sum(alleleTn[,-1]) # No. of reads containing TX2094 and F1 diagnostic SNPs ** 11657761
sum(alleleT[,-1]) # No. of reads containing TX2094 diagnostic SNPs 9894554

# check genes without allelic counts
check<-alleleMn$gene[rowSums(cbind(alleleMn[,-1],alleleTn[,-1]))>0]
length(check) # No. of genes containing allelic counts 29334
length(check)>length(genes) # bigger than final reported number
length(intersect(check,genes)) # 27816
length(setdiff(check,genes)) # 1518
length(setdiff(genes,check)) # 4
# check results mostly agree with genes containning disgnostic snps, the discrepancy probably came from MASKED snps, which means low coverage SNP site in some accession but probably allowing read assignment in some other accessions
# SO just focus on the common set
genes_diagnostic_with_alleles <- intersect(check,genes)
# save(alleleM, alleleMn, alleleT, alleleTn, genes, check, genes_diagnostic_with_alleles, file="Allelic_read_count.Rdata" )


#####################################################
################ Total Expression ###################
#####################################################

# import read count table:  .expression.txt file
x<-read.table(grep("expression.txt",list.files(),value=TRUE), header=TRUE, sep="\t")
dd<-which(duplicated(x$GENE))
count <- x[-dd,-1]
rownames(count) <- x$GENE[-dd]
names(count) <- gsub("^F1[.]|^TX2094[.]|^Maxxa[.]","",names(count))

#####################################################
#################### Make files  ####################
#####################################################
write.table(count, file = "totalReadCounts.txt", sep="\t")
write.table(alleleMn[alleleMn$gene %in% genes_diagnostic_with_alleles,], row.names=FALSE, file = "MaxxaAlleleReadCounts.txt", sep="\t")
write.table(alleleTn[alleleTn$gene %in% genes_diagnostic_with_alleles,], row.names=FALSE, file = "TX2094AlleleReadCounts.txt", sep="\t")

q()
n

#####################################################
################# File Discription ##################
#####################################################

# 'totalReadCounts.txt' -  counts of total gene expression, genes (row) X 26 RNA-seq libraries (column)
# 'MaxxaAlleleReadCounts.txt' - for 12 F1 libraries (column), allele-specific read count assigned to Maxxa alleles were presented for these genes contain diagnostic SNPs.
# 'TX2094AlleleReadCounts.txt' - same as above but for TX2094 allele.
# 'TM1.nuclNmt.transcript.len.txt' - gene length for all genes
# 'TM1.nuclNmt.transcript.len.diagnostic.txt' - gene length for genes with diagnostic allele-specifc SNPs