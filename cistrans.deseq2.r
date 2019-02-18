#####################################################
################# File Discription ##################
#####################################################

# 'totalReadCounts.txt' -  counts of total gene expression, genes (row) X 26 RNA-seq libraries (column)
# 'MaxxaAlleleReadCounts.txt' - for 12 F1 libraries (column), allele-specific read count assigned to Maxxa alleles were presented for these genes contain diagnostic SNPs.
# 'TX2094AlleleReadCounts.txt' - same as above but for TX2094 allele.
# 'TM1.nuclNmt.transcript.len.txt' - gene length for all genes


#####################################################
########## Finalize the sample inclusion #############
#####################################################

# load total counts
count <-read.table("totalReadCounts.txt", header=TRUE, sep="\t")
# requiring a read depth of 5 per genes, how many expressed genes in each accession?
getTotalExpressed=function(x,depth=5)
{
    length(which(rowSums(x)>=depth))
}
getTotalExpressed(count)
getTotalExpressed(count) # 56838
getTotalExpressed(count[,grep("Maxxa",names(count))]) # 51347
getTotalExpressed(count[,grep("TX2094",names(count))]) # 52953
getTotalExpressed(count[,grep("F1",names(count))]) # 52802
getTotalExpressed(count[,grep("F1.*[1-3]$",names(count))]) # 49531
getTotalExpressed(count[,grep("F1.*[4-6]$",names(count))]) # 49953
# Union of expressed fiber genes
length(which(rowSums(count[,grep("Maxxa",names(count))])>=5 | rowSums(count[,grep("TX2094",names(count))])>=5 | rowSums(count[,grep("F1.*[1-3]$",names(count))])>=5 |rowSums(count[,grep("F1.*[4-6]$",names(count))])>=5))
# 55551
length(which(rowSums(count[,grep("Maxxa",names(count))])>=5 | rowSums(count[,grep("TX2094",names(count))])>=5 | rowSums(count[,grep("F1.*[1-3]$",names(count))])>=5 |rowSums(count[,grep("F1.*[4-6]$",names(count))])>=5))/nrow(count)
# 0.8339739
length(which(rowSums(count[,grep("Maxxa",names(count))])>=5 & rowSums(count[,grep("TX2094",names(count))])>=5 & rowSums(count[,grep("F1.*[1-3]$",names(count))])>=5 &rowSums(count[,grep("F1.*[4-6]$",names(count))])>=5))
# overlap set 46411

# load total counts from bam idxstat
load("~/jfw-lab/Projects/AD1_domestication_cistrans/bowtie2_mapping/Ranalysis/counts.rdata")->l
l
# "info"       "total"      "count.rld"  "count.log2"
select<-which( (info$source=="yb"|info$source=="jj")& (info$genome=="Maxxa"|info$genome=="TX2094"|info$genome=="Yuc"|info$genome=="F1") & (info$dpa=="10dpa"|info$dpa =="20dpa"))
total<-total[,select]
info<-info[select,]
info$ID <- gsub("Yuc", "TX2094",paste(info$genome,info$dpa,info$rep, sep="."))
# need to match the sample orders
names(total)
info$ID
total<-total[,match(gsub("x$","",names(count)),info$ID)]
info<- info[match(gsub("x$","",names(count)),info$ID),]
# double check
rownames(info) == names(total)
info$ID == gsub("x$","",names(count))
rownames(info)<-names(count)
info$exclude<-factor(paste0(info$genome,gsub(".*[0-9]","_",rownames(info))))
# match gene order
total<-total[rownames(count),]
unique(rownames(total)==rownames(count))
# not exactly the same, because hylite missed the last gene, but proceed with PCA plots
colSums(total)-colSums(count)

## deseq normalization
library(DESeq2)
dds <- DESeqDataSetFromMatrix( countData = count, colData = info, design = ~ID)
# rlog transformation, note that with defalt blind=TRUE, design is actually ~1
rld <- rlog(dds)
count.rld <- as.data.frame(assay(rld))
names(count.rld) == names(total)
names(count.rld) == rownames(info)

# also possible to perform custom transformation:
dds <- estimateSizeFactors(dds)
# shifted log of normalized counts
se <- SummarizedExperiment(log2(counts(dds, normalized=TRUE) + 1), colData=colData(dds))
count.log2 <- as.data.frame(assay(se))

# RPKM transformation, Reads Per Kilobase of transcript per Million mapped reads.
# RPKM= 10^9 * number of mapped reads in gene
#             ----------------------------------------------------
#             total reads *  exon length
## rpkm normalization
len  <-read.table("~/jfw-lab/Projects/AD1_domestication_cistrans/mapping_reference/TM1.nuclNmt.transcript.len.txt",  sep="\t")
geneLen <- len$V2
names(geneLen) <- len$V1
geneLen <- geneLen[rownames(count)]
rpkm <-sweep(sweep(count*10^9, 2 ,colSums(count), FUN="/"), 1, geneLen, FUN ="/" )
rpkm.log2 <-log2(rpkm+1)

# my custom function
library(ggplot2)
library(scales)
library(ape)
plotGrouping <- function(norm_log, color, shape, text, tip, save = "plotGrouping.pdf"){
    # norm<-sweep(total,2,info$lib_size,"/")*10^6
    # norm_log <- log2(norm+1)
    pca=prcomp(t(norm_log))
    dat = as.data.frame(pca$x)
    proportion<-summary(pca)$importance[2,1:2]*100
    proportion<-paste0(names(proportion)," (", proportion, "%)")
    p<-ggplot(aes(PC1, PC2, color=color,    shape=shape),data=dat) + geom_point() +xlab(proportion[1]) + ylab(proportion[2])
    pdf(save)
    print( p + geom_text(aes_string(x = "PC1", y = "PC2", label = text), color="grey", hjust = 0, nudge_x = 0.09) )
    
    hc<-hclust( dist(t(norm_log)) )
    tre <- as.phylo(hc)
    tre$tip.label <- as.character(tip)
    # try to match ggplot color: library(scale);
    # show_col(col4<-hue_pal()(4))
    tipCol <- color
    levels(tipCol) <- hue_pal()(nlevels(tipCol))
    plot(tre,tip.col =as.character(tipCol),  type="unrooted",cex=0.6, no.margin=TRUE)
    plot(tre,tip.col =as.character(tipCol), type="fan",cex=0.6, no.margin=TRUE)
    dev.off()   }

# plots
plotGrouping(count.rld, color=info$dpa, shape=info$source, tip=rownames(info), text=info$exclude, save = "plotGrouping.rld.pdf")
plotGrouping(count.log2, color=info$dpa, shape=info$source, tip=rownames(info),text=info$exclude, save = "plotGrouping.log2.pdf")
plotGrouping(rpkm.log2, color=info$dpa, shape=info$source, tip=rownames(info),text=info$exclude, save = "plotGrouping.log2rpkm.pdf")


#####################################################
################# DGE for all genes #################
#####################################################
pairwiseDE<-function(dds, contrast,savePath)
{
    # DE analysis
    print(contrast)
    ddsPW <-dds[,dds$condition %in% contrast]
    ddsPW$condition<-droplevels(ddsPW$condition)
    res <- results(DESeq(ddsPW))
    print( summary(res,alpha=.05) ) # print results
    write.table(res, file=paste(savePath,"DE/",paste(contrast, collapse="vs"),".txt", sep=""), sep="\t")
}
getSig<-function(res,fc.threshold=0,direction=NULL){
    sig<- res[res$padj<0.05 & !is.na(res$padj) & abs(res$log2FoldChange)>=fc.threshold,]
    if(is.null(direction)){
        n<-nrow(sig)
    }else if(direction=="up"){
        n<-nrow(sig[sig$log2FoldChange>0,])
    }else if(direction=="down"){
        n<-nrow(sig[sig$log2FoldChange<0,])
    }
    return(n)
}


countT<-count[,!grepl("x$",rownames(info))]
coldataT <- info[!grepl("x$",rownames(info)),]
coldataT$condition<-factor(gsub("dpa.*","dpa",coldataT$ID) )
names(countT)==rownames(coldataT)
coldataT$lib_size - colSums(countT)  # slight difference
coldataT$lib_size <- colSums(countT) # modify

# differential expression to get expression divergence
library(DESeq2)
dds <- DESeqDataSetFromMatrix( countData = countT, colData = coldataT, design = ~ condition)
# pairwise deseq workflow
batch<- rbind(
c("Maxxa.10dpa", "TX2094.10dpa" ),
c("Maxxa.10dpa", "F1.10dpa"),
c("TX2094.10dpa", "F1.10dpa"),

c("Maxxa.20dpa", "TX2094.20dpa" ),
c("Maxxa.20dpa", "F1.20dpa"),
c("TX2094.20dpa", "F1.20dpa"),

c("Maxxa.10dpa", "Maxxa.20dpa" ),
c("TX2094.10dpa", "TX2094.20dpa" ),
c("F1.10dpa", "F1.20dpa" )  )

# make a "DE" folder
system("mkdir DE")
apply(batch,1,function(x) pairwiseDE(dds,x,savePath = ""))

# analyze MxT and TxM
coldataT2<-coldataT
mt<-which(coldataT2$genome=="F1" & coldataT2$rep %in% 1:3)
tm<-which(coldataT2$genome=="F1" & coldataT2$rep %in% 4:6)
coldataT2$condition<-as.character(coldataT2$condition)
coldataT2$condition[mt]<-gsub("F1","MxT",coldataT2$condition[mt])
coldataT2$condition[tm]<-gsub("F1","TxM",coldataT2$condition[tm])
dds2 <- DESeqDataSetFromMatrix( countData = countT, colData = coldataT2, design = ~ condition)
batch2<- rbind(
c("MxT.10dpa", "TxM.10dpa" ),
c("Maxxa.10dpa", "MxT.10dpa"),
c("Maxxa.10dpa", "TxM.10dpa"),
c("TX2094.10dpa", "MxT.10dpa"),
c("TX2094.10dpa", "TxM.10dpa"),

c("MxT.20dpa", "TxM.20dpa" ),
c("Maxxa.20dpa", "MxT.20dpa"),
c("Maxxa.20dpa", "TxM.20dpa"),
c("TX2094.20dpa", "MxT.20dpa"),
c("TX2094.20dpa", "TxM.20dpa"),

c("MxT.10dpa", "MxT.20dpa" ),
c("TxM.10dpa", "TxM.20dpa" ) )
apply(batch2,1,function(x) pairwiseDE(dds2,x,savePath = ""))


# print out results in the order of batch
#fileL<- list.files("DE")
fileL<-paste0(batch[,1],"vs",batch[,2],".txt")
fileL<-c(fileL, paste0(batch2[,1],"vs",batch2[,2],".txt"))
sigT<-c("sample 1","sample 2","DE (q<0.05)","1>2","2>1")
for(file in fileL)
{
    res<-read.table(paste0("DE/",file),sep="\t",header=TRUE)
    sigRes <- c(unlist(strsplit(gsub(".txt","",file),split="vs") ), getSig(res),getSig(res,direction="up"),getSig(res,direction="down"))
    sigT<-rbind(sigT,sigRes)
}
T<-as.data.frame(sigT[-1,], row.names=FALSE)
names(T)<-sigT[1,]
print(T)

library(gplots)
pdf("DE/checkDE.pdf")
textplot(T)
mtext("DE analysis result summary")
dev.off()

# rename this "DE" folder as "DE0"
system("mv DE/ DE-66610/")

# save total datasets
write.table(countT, file="Genes66610.raw.count.txt", sep="\t")
write.table(coldataT, file="Genes66610.sample.info.txt", sep="\t")

# normalization with rpkm
len  <-read.table("~/cistrans/TM1.nuclNmt.transcript.len.txt",  sep="\t")
geneLen <- len$V2
names(geneLen) <- len$V1
geneLen <- geneLen[rownames(countT)]
rpkm <-sweep(sweep(countT*10^9, 2 ,coldataT$lib_size, FUN="/"), 1, geneLen, FUN ="/" )
rpkm.log2 <-log2(rpkm+1)
plotGrouping(rpkm.log2, color=coldataT$dpa, shape=coldataT$genome, tip=rownames(coldataT),text=coldataT$rep, save = "plotGrouping.66610.log2rpkm.pdf")

# save total datasets
write.table(rpkm, file="Genes66610.log2rpkm.txt", sep="\t")
write.table(geneLen, file="Genes66610.gene.length.txt", sep="\t")

# read annotation
x<-read.table("Ghirsutum_458_v1.1.annotation_info.txt",sep="\t",comment.char="", quote="\"",header=TRUE)
annot<-x[grep("[.]1$",x$transcriptName),]
dim(annot)
annot<-annot[,c(3,5:13)]
desc<-annot[,c("transcriptName","arabi.defline")]
y<-read.table("mt.annotation.txt",header=TRUE,sep="\t",comment.char="", quote="\"")
y<-y[,1:2]
names(y)<-names(desc)
desc<-rbind(desc,y[,1:2])
z<-merge(desc, annot, by="transcriptName",all.x=TRUE)
rownames(z)<-z$transcriptName
annotation<-z[rownames(rpkm),]
write.table(annotation, file="Genes66610.transcrip.description.txt", sep="\t")


#####################################################
############# Finalize gene inclusion ###############
#####################################################

# import allelic read count tables
maxxa  <-read.table("MaxxaAlleleReadCounts.txt", header=TRUE, sep="\t")
tx2094 <-read.table("TX2094AlleleReadCounts.txt", header=TRUE, sep="\t")
# check column order
names(maxxa)[-1] == names(countT)[1:12]
names(maxxa) == names(tx2094)
# combine
names(maxxa)[-1] <- paste0(names(maxxa)[-1],".maxxa")
names(tx2094)[-1] <- paste0(names(tx2094)[-1],".tx2094")
allele <- merge(maxxa, tx2094, by ="gene")

# prepare for DEseq2
countA <- allele[,-1]
rownames(countA) <- allele$gene
coldataA <- rbind(coldataT[1:12,], coldataT[1:12,])
rownames(coldataA) <- names(countA)
coldataA$allele <-gsub(".*[.]", "",rownames(coldataA))
coldataA$condition <-paste0(coldataA$allele, ".", coldataA$condition)
# check percentage of allele read count in total
c(colSums(maxxa[,-1]),colSums(tx2094[,-1]))/coldataA$lib_size # ~7.6%
allele_lib_size <- c(colSums(maxxa[,-1]),colSums(tx2094[,-1]))

# add countN
names(countT) <- paste0(names(countT),".total")
rownames(coldataT) <- names(countT)
coldataT$allele <-"total"

# allele + total
countAll <- cbind(countA, countT[rownames(countA),])
coldataAll <- rbind(coldataA, coldataT)
coldataAll$condition2 <-coldataAll$condition
coldataAll$condition2[coldataAll$genome=="F1"&coldataAll$rep %in%1:3]<-gsub("F1","MxT",coldataAll$condition2[coldataAll$genome=="F1"&coldataAll$rep %in%1:3])
coldataAll$condition2[coldataAll$genome=="F1"&coldataAll$rep %in%4:6]<-gsub("F1","TxM",coldataAll$condition2[coldataAll$genome=="F1"&coldataAll$rep %in%4:6])


# normalization with rpkm
len  <-read.table("~/jfw-lab/Projects/AD1_domestication_cistrans/mapping_reference/TM1.nuclNmt.transcript.len.txt",  sep="\t")
geneLen <- len$V2
names(geneLen) <- len$V1
geneLen <- geneLen[rownames(countAll)]
rpkm <-sweep(sweep(countAll*10^9, 2 ,coldataAll$lib_size, FUN="/"), 1, geneLen, FUN ="/" )
rpkm.log2 <-log2(rpkm+1)
plotGrouping(rpkm.log2, color=coldataAll$dpa, shape=coldataAll$allele, tip=rownames(coldataAll),text=coldataAll$genome, save = "plotGrouping.27816.log2rpkm.pdf")

# read annotation
x<-read.table("Ghirsutum_458_v1.1.annotation_info.txt",sep="\t",comment.char="", quote="\"",header=TRUE)
annot<-x[grep("[.]1$",x$transcriptName),]
dim(annot)
annot<-annot[,c(3,5:13)]
desc<-annot[,c("transcriptName","arabi.defline")]
y<-read.table("mt.annotation.txt",header=TRUE,sep="\t",comment.char="", quote="\"")
y<-y[,1:2]
names(y)<-names(desc)
desc<-rbind(desc,y[,1:2])
z<-merge(desc, annot, by="transcriptName",all.x=TRUE)
rownames(z)<-z$transcriptName
annotation<-z[rownames(rpkm),]

write.table(countAll, file="Genes27816.raw.count.txt", sep="\t")
write.table(coldataAll, file="Genes27816.sample.info.txt", sep="\t")
write.table(rpkm, file="Genes27816.rpkm.txt", sep="\t")
write.table(geneLen, file="Genes27816.gene.length.txt", sep="\t")
write.table(annotation, file="Genes27816.transcrip.description.txt", sep="\t")



#####################################################
############### Differential expression #############
#####################################################
pairwiseDE<-function(dds, contrast,savePath)
{
    # DE analysis
    print(contrast)
    ddsPW <-dds[,dds$condition %in% contrast]
    ddsPW$condition<-droplevels(ddsPW$condition)
    res <- results(DESeq(ddsPW))
    print( summary(res,alpha=.05) ) # print results
    write.table(res, file=paste(savePath,"DE/",paste(contrast, collapse="vs"),".txt", sep=""), sep="\t")
}


# differential expression to get allelic expression divergence
library(DESeq2)
dds <- DESeqDataSetFromMatrix( countData = countAll, colData = coldataAll, design = ~ condition)
# if using pre-existing size factors; not necessary
# sizeFactors(dds) <- coldataAll$lib_size/mean(coldataAll$lib_size)

# pairwise deseq workflow, c(maxxa, tx2094), ratio calculated as tx2094/maxxa
batch<- rbind(
c("maxxa.F1.10dpa", "tx2094.F1.10dpa" ),
c("maxxa.F1.20dpa", "tx2094.F1.20dpa" ),
c("Maxxa.10dpa", "TX2094.10dpa" ),
c("Maxxa.20dpa", "TX2094.20dpa" ),
c("TX2094.10dpa", "F1.10dpa"),
c("TX2094.20dpa", "F1.20dpa"),
c("Maxxa.10dpa", "F1.10dpa"),
c("Maxxa.20dpa", "F1.20dpa") )
system("mkdir DE")
apply(batch,1,function(x) pairwiseDE(dds,x,savePath = ""))

### Get allelic expression divergence, test for B=0
B10<-read.table("DE/maxxa.F1.10dpavstx2094.F1.10dpa.txt", header=TRUE, sep="\t")
B20<-read.table("DE/maxxa.F1.20dpavstx2094.F1.20dpa.txt", header=TRUE, sep="\t")
dim(B10[B10$padj<0.05 & !is.na(B10$padj),]) # 3800
dim(B20[B20$padj<0.05 & !is.na(B20$padj),]) # 2887
length(unique(c(which(B10$padj<0.05 & !is.na(B10$padj)),which(B20$padj<0.05 & !is.na(B20$padj))))) #4760


#### Get parental expression divergence, test for A=0
A10<-read.table("DE/Maxxa.10dpavsTX2094.10dpa.txt", header=TRUE, sep="\t")
A20<-read.table("DE/Maxxa.20dpavsTX2094.20dpa.txt", header=TRUE, sep="\t")
dim(A10[A10$padj<0.05 & !is.na(A10$padj),]) # 5168
dim(A20[A20$padj<0.05 & !is.na(A20$padj),]) # 3528
length(unique(c(which(A10$padj<0.05 & !is.na(A10$padj)),which(A20$padj<0.05 & !is.na(A20$padj))))) #7303
length(intersect(which(A10$padj<0.05 & !is.na(A10$padj)),which(A20$padj<0.05 & !is.na(A20$padj)))) #1393

# analyze MxT and TxM
coldataAll2<-coldataAll
coldataAll2$condition<-coldataAll2$condition2
dds2 <- DESeqDataSetFromMatrix( countData = countAll, colData = coldataAll2, design = ~ condition)
# pairwise deseq workflow, c(maxxa, tx2094), ratio calculated as tx2094/maxxa
batchMT<- rbind(
c("maxxa.MxT.10dpa", "tx2094.MxT.10dpa" ),
c("maxxa.MxT.20dpa", "tx2094.MxT.20dpa" ),
c("TX2094.10dpa", "MxT.10dpa"),
c("TX2094.20dpa", "MxT.20dpa"),
c("Maxxa.10dpa", "MxT.10dpa"),
c("Maxxa.20dpa", "MxT.20dpa") )
batchTM<- rbind(
c("maxxa.TxM.10dpa", "tx2094.TxM.10dpa" ),
c("maxxa.TxM.20dpa", "tx2094.TxM.20dpa" ),
c("TX2094.10dpa", "TxM.10dpa"),
c("TX2094.20dpa", "TxM.20dpa"),
c("Maxxa.10dpa", "TxM.10dpa"),
c("Maxxa.20dpa", "TxM.20dpa") )
batchff<-rbind(
c("MxT.10dpa", "TxM.10dpa"),
c("MxT.20dpa", "TxM.20dpa"),
c("maxxa.MxT.10dpa", "maxxa.TxM.10dpa"),
c("maxxa.MxT.20dpa", "maxxa.TxM.20dpa"),
c("tx2094.MxT.10dpa", "tx2094.TxM.10dpa"),
c("tx2094.MxT.20dpa", "tx2094.TxM.20dpa")
)

apply(batchMT,1,function(x) pairwiseDE(dds2,x,savePath = ""))
apply(batchTM,1,function(x) pairwiseDE(dds2,x,savePath = ""))
apply(batchff,1,function(x) pairwiseDE(dds2,x,savePath = ""))

### Get allelic expression divergence, test for B=0
B.mt10<-read.table("DE/maxxa.MxT.10dpavstx2094.MxT.10dpa.txt", header=TRUE, sep="\t")
B.mt20<-read.table("DE/maxxa.MxT.20dpavstx2094.MxT.20dpa.txt", header=TRUE, sep="\t")
dim(B.mt10[B.mt10$padj<0.05 & !is.na(B.mt10$padj),]) # 1965    6
dim(B.mt20[B.mt20$padj<0.05 & !is.na(B.mt20$padj),]) # 1190    6

B.tm10<-read.table("DE/maxxa.TxM.10dpavstx2094.TxM.10dpa.txt", header=TRUE, sep="\t")
B.tm20<-read.table("DE/maxxa.TxM.20dpavstx2094.TxM.20dpa.txt", header=TRUE, sep="\t")
dim(B.tm10[B.tm10$padj<0.05 & !is.na(B.tm10$padj),]) # 1827    6
dim(B.tm20[B.tm20$padj<0.05 & !is.na(B.tm20$padj),]) # 1378    6

# important: A10, B10, B.mt10, B.tm10, A20, B20, B.mt20, B.tm20,
save(A10, B10, B.mt10, B.tm10, A20, B20, B.mt20, B.tm20, file="cistrans.rdata")
system("mv DE/ DE-27816/")



#####################################################
############### Visualize DE results ################
#####################################################
getSig<-function(res,fc.threshold=0,direction=NULL){
    sig<- res[res$padj<0.05 & !is.na(res$padj) & abs(res$log2FoldChange)>=fc.threshold,]
    if(is.null(direction)){
        n<-nrow(sig)
    }else if(direction=="up"){
        n<-nrow(sig[sig$log2FoldChange>0,])
    }else if(direction=="down"){
        n<-nrow(sig[sig$log2FoldChange<0,])
    }
    return(n)
}

#fileL<- list.files("DE27816")
fileL<-paste0(batch[,1],"vs",batch[,2],".txt")
fileL<-c(fileL, paste0(batchMT[,1],"vs",batchMT[,2],".txt"))
fileL<-c(fileL, paste0(batchTM[,1],"vs",batchTM[,2],".txt"))
fileL<-c(fileL, paste0(batchff[,1],"vs",batchff[,2],".txt"))

sigT<-c("sample 1","sample 2","DE (q<0.05)","1>2","2>1")
for(file in fileL)
{
    res<-read.table(paste0("DE-27816/",file),sep="\t",header=TRUE)
    sigRes <- c(unlist(strsplit(gsub(".txt","",file),split="vs") ), getSig(res),getSig(res,direction="up"),getSig(res,direction="down"))
    sigT<-rbind(sigT,sigRes)
}
T<-as.data.frame(sigT[-1,], row.names=FALSE)
names(T)<-sigT[1,]
write.table(T,file="Genes27816.DEsummary.txt",sep="\t")

library(gplots)
pdf("checkDE.pdf")
textplot(T)
mtext("DE analysis result summary")


##### Check DE results with different cutoffs
id<-rownames(B10)

fc <- seq(1,3,by=0.1)
getSig<-function(res,fc.threshold=0){no<-nrow(res[res$padj<0.05 & !is.na(res$padj) & abs(res$log2FoldChange)>=fc.threshold,]);return(no)}
AB<-data.frame(
foldChangeCutoff =fc,
B.tm10=sapply(log2(fc),function(x)getSig(B.tm10,x)),
B.tm20=sapply(log2(fc),function(x)getSig(B.tm20,x)),
B.mt10=sapply(log2(fc),function(x)getSig(B.mt10,x)),
B.mt20=sapply(log2(fc),function(x)getSig(B.mt20,x)),
B.10=sapply(log2(fc),function(x)getSig(B10,x)),
B.20=sapply(log2(fc),function(x)getSig(B20,x)),
A.10=sapply(log2(fc),function(x)getSig(A10,x)),
A.20=sapply(log2(fc),function(x)getSig(A20,x))
)
library(RColorBrewer)
n<-ncol(AB)-1
col=brewer.pal(n,"Dark2")
plot(AB[,2]~foldChangeCutoff,data=AB, type="l", col=col[1], ylim=c(1000,7000), ylab="DE genes", main="DE analysis with different cutoffs")
for(i in 1:n)
{
    lines(AB[,i+1]~foldChangeCutoff, data=AB, type="l", col=col[i])
}
#legend("topright",legend=names(AB)[2:7],col=col,pch=21)
text(1.1, AB[1,2:(n+1)]+100, names(AB)[2:(n+1)],)

par(mfrow=c(2,2))

#10dpa
# compare log2FoldChange distribution
# No apparant difference between A and B
plot(density(na.omit(A10[id,"log2FoldChange"])), main="log2 fold change, 10 dpa")
lines(density(na.omit(B10[id,"log2FoldChange"])),col="purple")
lines(density(na.omit(B.tm10[id,"log2FoldChange"])),col="blue")
lines(density(na.omit(B.mt10[id,"log2FoldChange"])),col="green")
abline(v=0,col="light grey")
legend("topright",col=c("black","purple","blue","green"),legend=c("A","B","B - TxM","B - MxT"),pch="-")
# compare standard errors
# Higher standard error in B than A
plot( A10[id,"lfcSE"], B10[id,"lfcSE"],pch=".",main="Standard error, 10 dpa",xlab="A",ylab="B")
lines(c(0,1),c(0,1),col="blue")
plot( A10[id,"lfcSE"], B.mt10[id,"lfcSE"],pch=".",main="Standard error, 10 dpa",xlab="A",ylab="B - MxT")
lines(c(0,1),c(0,1),col="blue")
plot( A10[id,"lfcSE"], B.tm10[id,"lfcSE"],pch=".",main="Standard error, 10 dpa",xlab="A",ylab="B - TxM")
lines(c(0,1),c(0,1),col="blue")

#20dpa
# compare log2FoldChange distribution
# No apparant difference between A and B
plot(density(na.omit(A20[id,"log2FoldChange"])), main="log2 fold change, 20 dpa")
lines(density(na.omit(B20[id,"log2FoldChange"])),col="purple")
lines(density(na.omit(B.tm20[id,"log2FoldChange"])),col="blue")
lines(density(na.omit(B.mt20[id,"log2FoldChange"])),col="green")
abline(v=0,col="light grey")
legend("topright",col=c("black","purple","blue","green"),legend=c("A","B","B - TxM","B - MxT"),pch="-")
# compare standard errors
# Higher standard error in B than A
plot( A20[id,"lfcSE"], B20[id,"lfcSE"],pch=".",main="Standard error, 20 dpa",xlab="A",ylab="B")
lines(c(0,1),c(0,1),col="blue")
plot( A20[id,"lfcSE"], B.mt20[id,"lfcSE"],pch=".",main="Standard error, 20 dpa",xlab="A",ylab="B - MxT")
lines(c(0,1),c(0,1),col="blue")
plot( A20[id,"lfcSE"], B.tm20[id,"lfcSE"],pch=".",main="Standard error, 20 dpa",xlab="A",ylab="B - TxM")
lines(c(0,1),c(0,1),col="blue")

### Check standard error of A and B
# volcano plot
plotVolvano<-function(res, title)
{
    plot(res[id,"log2FoldChange"], -log2(res[id,"padj"]), main=title, xlab="log2FoldChange", ylab="-log2padj",pch=".",ylim=c(0,200))
    abline(h=-log2(0.05))
}
plotVolvano(A10, "A, 10 dpa")
plotVolvano(B10, "B, F1 10 dpa")
plotVolvano(B.mt10, "B, MxT 10 dpa")
plotVolvano(B.tm10, "B, TxM 10 dpa")
plotVolvano(A20, "A, 20 dpa")
plotVolvano(B20, "B, F1 20 dpa")
plotVolvano(B.mt20, "B, MxT 20 dpa")
plotVolvano(B.tm20, "B, TxM 20 dpa")


# compare log2 Fold change
plot( A10[id,"log2FoldChange"], B10[id,"log2FoldChange"],pch=".",main="log2 Fold Change, 10 dpa",xlab="A",ylab="B - F1")
lines(c(-6,6),c(-6,6),col="blue")
plot( A10[id,"log2FoldChange"], B.mt10[id,"log2FoldChange"],pch=".",main="log2 Fold Change, 10 dpa",xlab="A",ylab="B - MxT")
lines(c(-6,6),c(-6,6),col="blue")
plot( A10[id,"log2FoldChange"], B.tm10[id,"log2FoldChange"],pch=".",main="log2 Fold Change, 10 dpa",xlab="A",ylab="B - TxM")
lines(c(-6,6),c(-6,6),col="blue")
plot( B.mt10[id,"log2FoldChange"], B.tm10[id,"log2FoldChange"],pch=".",main="log2 Fold Change, 10 dpa",xlab="B - MxT",ylab="B - TxM")
lines(c(-6,6),c(-6,6),col="blue")


plot( A20[id,"log2FoldChange"], B20[id,"log2FoldChange"],pch=".",main="log2 Fold Change, 20 dpa",xlab="A",ylab="B - F1")
lines(c(-6,6),c(-6,6),col="blue")
plot( A20[id,"log2FoldChange"], B.mt20[id,"log2FoldChange"],pch=".",main="log2 Fold Change, 20 dpa",xlab="A",ylab="B - MxT")
lines(c(-6,6),c(-6,6),col="blue")
plot( A20[id,"log2FoldChange"], B.tm20[id,"log2FoldChange"],pch=".",main="log2 Fold Change, 20 dpa",xlab="A",ylab="B - TxM")
lines(c(-6,6),c(-6,6),col="blue")
plot( B.mt20[id,"log2FoldChange"], B.tm20[id,"log2FoldChange"],pch=".",main="log2 Fold Change, 20 dpa",xlab="B - MxT",ylab="B - TxM")
lines(c(-6,6),c(-6,6),col="blue")

dev.off()


#####################################################
############### Cis/Trans Analysis ##################
#####################################################


### comparing A-B= 0 is tricky, both are log2FoldChange and its standard error lfcse
# maybe I can compare with t test
#### T test from means and standard errors ####
# m1, m2: the sample means
# s1, s2: the sample standard deviations
# n1, n2: the same sizes
# se1, se2: the sample standard errors
# se1 <- s1/sqrt(n)
# m0: the null value for the difference in means to be tested for. Default is 0.
# equal.variance: whether or not to assume equal variance. Default is FALSE.
t.test2 <- function(m1,m2,se1,se2,n1,n2,m0=0,equal.variance=FALSE)
{
    if( equal.variance==FALSE )
    {
        # se <- sqrt( (s1^2/n1) + (s2^2/n2) )
        se <- sqrt( (se1^2) + (se2^2) )
        # welch-satterthwaite df
        # df <- ( (s1^2/n1 + s2^2/n2)^2 )/( (s1^2/n1)^2/(n1-1) + (s2^2/n2)^2/(n2-1) )
        df <- ( (se1^2 + se2^2)^2 )/( (se1^2)^2/(n1-1) + (se2^2)^2/(n2-1) )
    } else
    {
        # pooled standard deviation, scaled by the sample sizes
        # se <- sqrt( (1/n1 + 1/n2) * ((n1-1)*s1^2 + (n2-1)*s2^2)/(n1+n2-2) )
        se <- sqrt( (1/n1 + 1/n2) * ((n1-1)*s1^2 + (n2-1)*s2^2)/(n1+n2-2) )
        df <- n1+n2-2
    }
    t <- (m1-m2-m0)/se
    dat <- c(m1-m2, se, t, 2*pt(-abs(t),df))
    names(dat) <- c("Difference of means", "Std Error", "t", "p-value")
    return(dat)
}

x1 = rnorm(3)
x2 = rnorm(3)
# you'll find this output agrees with that of t.test when you input x1,x2
t.test(x1,x2)
t.test2( mean(x1),  mean(x2), sd(x1)/sqrt(3), sd(x2)/sqrt(3), 3,3)

### Make a categorization table for each F1 condition: TxM10, MxT10, TxM20, MxT20,
criteria<- as.data.frame(rbind(c("A!=0;B!=0;A=B", "1.Cis only"),
c("A!=0;B=0;A!=B", "2.Trans only"),
c("A!=0;B!=0;A!=B", "Cis+Trans"),
c("A=0;B!=0;A!=B", "5.Compensatory"),
c("A=0;B=0;A=B", "6.Conserved") ))
names(criteria) <- c("class","category")

# A.res <-A10
# B.res <-B10

classCisTrans<-function(A.res, B.res, A.n, B.n, log2fc.threshold=0,plotTitle=NULL)
{
    # A = log2(TX2094/Maxxa), cis + trans
    A <- A.res[,c("log2FoldChange", "lfcSE", "padj")]
    names(A) <- c("A", "A.SE", "A.padj")
    # A = log2(F1_t/F1_m), cis
    B <- B.res[,c("log2FoldChange", "lfcSE", "padj")]
    names(B) <- c("B", "B.SE", "B.padj")
    
    A<-A[rownames(B),]
    table <- cbind(A,B)
    table$AminusB <- table$A - table$B
    table$AminusB.pvalue <- apply(table,1,function(x) t.test2(m1=x[1],m2=x[4],se1=x[2], se2=x[5], n1=A.n, n2=B.n)["p-value"])
    
    table$cisNtrans <- ifelse(abs(table$A)>=log2fc.threshold & table$A.padj<0.05 & !is.na(table$A.padj), "A!=0", "A=0")
    table$cis <- ifelse(abs(table$B)>=log2fc.threshold & table$B.padj<0.05 & !is.na(table$B.padj), "B!=0", "B=0")
    table$trans <- ifelse(abs(table$AminusB)>=log2fc.threshold & table$AminusB.pvalue<0.05, "A!=B", "A=B")
    table$class <- paste(table$cisNtrans,table$cis,table$trans,sep=";")
    table$category <- as.character(criteria$category[ match(table$class,criteria$class)])
    table$category[is.na(table$category)] <- "7.Ambiguous"
    table$category[ table$category=="Cis+Trans" & table$B*table$AminusB >0 ] <- "3.Cis+Trans: enhancing"
    table$category[ table$category=="Cis+Trans" & table$B*table$AminusB <0 ] <- "4.Cis+Trans: compensating"
    
    colors <- c("red","blue","purple","brown","green","black","grey")
    if(!is.null(plotTitle)){
        p<- ggplot( table, aes(x=A, y=B, color=category)) + geom_point(alpha=0.8) +  xlab("Cis + Trans") + ylab("Cis") + ggtitle(plotTitle) + scale_color_manual(values=colors) + geom_vline(xintercept=0) + geom_hline(yintercept=0) +  geom_abline(intercept = 0, slope = 1) + theme_bw()
        # p2<-ggplot( table, aes( factor(category), fill=category)) + geom_bar(stat="count") + scale_fill_manual(values=colors)
        print(p)
    }
    return(table)
}
res10 <- classCisTrans(A.res = A10, B.res = B10, A.n = 3, B.n=6, log2fc.threshold=0)
res.mt10 <- classCisTrans(A.res = A10, B.res = B.mt10, A.n = 3, B.n=3, log2fc.threshold=0)
res.tm10 <- classCisTrans(A.res = A10, B.res = B.tm10, A.n = 3, B.n=3, log2fc.threshold=0)
res20 <- classCisTrans(A.res = A20, B.res = B20, A.n = 3, B.n=6, log2fc.threshold=0)
res.mt20 <- classCisTrans(A.res = A20, B.res = B.mt20, A.n = 3, B.n=3, log2fc.threshold=0)
res.tm20 <- classCisTrans(A.res = A20, B.res = B.tm20, A.n = 3, B.n=3, log2fc.threshold=0)


plotCisTrans<-function(table, plotTitle="")
{
    colors <- c("red","blue","purple","brown","green","black","grey")
    p<- ggplot( table, aes(x=A, y=B, color=category)) + geom_point(alpha=0.8) +  xlab("Cis + Trans") + ylab("Cis") + ggtitle(plotTitle) + scale_color_manual(values=colors) + geom_vline(xintercept=0) + geom_hline(yintercept=0) +  geom_abline(intercept = 0, slope = 1) + theme_bw()
        # p2<-ggplot( table, aes( factor(category), fill=category)) + geom_bar(stat="count") + scale_fill_manual(values=colors)
    print(p)
}

save(A10, B10, B.mt10, B.tm10, A20, B20, B.mt20, B.tm20, res10, res.mt10, res.tm10, res20, res.mt20, res.tm20, sumT, file="cistrans.rdata")

sumT<-as.data.frame(matrix(unlist(lapply(list(res.mt10,res.tm10,res.mt20, res.tm20,res10,res20), function(x)table(x$category)) ),ncol=6))
rownames(sumT)<-names(table(res10$category))
names(sumT)<-c("10dpa MxT","10dpa TxM","20dpa MxT", "20dpa TxM","10dpa F1","20dpa F1")
print(sumT)

pdf("plotCistrans.pdf")
textplot(sumT,cex=0.6)
plotCisTrans(res.mt10, "MxT 10dpa")
plotCisTrans(res.mt20, "MxT 20dpa")
plotCisTrans(res.tm10, "TxM 10dpa")
plotCisTrans(res.tm20, "TxM 20dpa")
plotCisTrans(res10, "F1 10dpa")
plotCisTrans(res20, "F1 20dpa")
dev.off()


# plot percentage 
library(reshape);
library(plyr)
# plot 7 categories
T=sumT[,1:4]; 
T$category=rownames(T)
dat=melt(T,id.var="category")
p1<-ggplot(data=dat, aes(x=variable, y=value, fill=category)) +
  geom_bar(stat="identity")+
  scale_fill_manual(values=colors)+
  theme_minimal()
df <- ddply(dat, "variable",  transform, label_ypos=rev(cumsum(rev(value))))
df <- ddply(df, "variable",  transform, t=round(value/sum(value)*100,1))
df$value2=paste0(df$value," (",df$t, "%)")
df$value2[grep("^[1-4]",df$category)]=""
p2<-ggplot(data=df, aes(x=variable, y=value, fill=category)) +
  geom_bar(stat="identity")+
  geom_text(aes(y=label_ypos, label=value2), vjust=1.6, 
            color="white", size=3.5)+
  scale_fill_manual(values=colors)+
  theme_minimal()
# plot 4 categories
T=sumT[1:4,1:4]; 
T$category=rownames(T)
dat=melt(T,id.var="category")
df <- ddply(dat, "variable",  transform, label_ypos=rev(cumsum(rev(value))))
df <- ddply(df, "variable",  transform, t=round(value/sum(value)*100,1))
df$value2=paste0(df$value," (",df$t, "%)")
p3<-ggplot(data=df, aes(x=variable, y=value, fill=category)) +
  geom_bar(stat="identity")+
  scale_fill_manual(values=colors)+
  theme_minimal()
p4<-ggplot(data=df, aes(x=variable, y=value, fill=category)) +
  geom_bar(stat="identity")+
  geom_text(aes(y=label_ypos, label=value2), vjust=1.6, 
            color="white", size=3.5)+
  scale_fill_manual(values=colors)+
  theme_minimal()
pdf("plotCistrans.barsummary.pdf")
p1;p2;p3;p4
dev.off()



write.table(sumT, file ="cistrans.summary.txt",sep="\t")
write.table(res10, file ="cistrans.F1_10dpa.txt",sep="\t")
write.table(res.mt10, file ="cistrans.MxT_10dpa.txt",sep="\t")
write.table(res.tm10, file ="cistrans.TxM_10dpa.txt",sep="\t")
write.table(res20, file ="cistrans.F1_20dpa.txt",sep="\t")
write.table(res.mt20, file ="cistrans.MxT_20dpa.txt",sep="\t")
write.table(res.tm20, file ="cistrans.TxM_20dpa.txt",sep="\t")


###################################################################
## Change fold change cutoff and inspect cis trans summary table ##
###################################################################

library(RColorBrewer)
col=brewer.pal(7,"Dark2")
n<-length(id)
tests<-rbind(c("A10","B.mt10","MxT10"),c("A10","B.tm10","TxM10"),c("A20","B.mt20","MxT20"),c("A20","B.tm20","TxM20"), c("A10","B10","F1_10"),c("A20","B20","F1_20"))
resultT<-list()
for(i in 1:4){
    res<-sapply(log2(fc),function(x){tt<-classCisTrans(A.res = get(tests[i,1]), B.res = get(tests[i,2]), 3,3,log2fc.threshold=x ); return(table(tt$category))} )
    colnames(res)<-fc
    resultT[[tests[i,3]]]<-res
}
resultT
for(i in 5:6){
    res<-sapply(log2(fc),function(x){tt<-classCisTrans(A.res = get(tests[i,1]), B.res = get(tests[i,2]), 3,6,log2fc.threshold=x ); return(table(tt$category))} )
    colnames(res)<-fc
    resultT[[tests[i,3]]]<-res
}

colors <- c("red","blue","purple","brown","green","black","grey")
pdf("checkCistran.by.cutoff.pdf")
noCut<-cbind(resultT[[1]][,1],resultT[[2]][,1],resultT[[3]][,1],resultT[[4]][,1],resultT[[5]][,1],resultT[[6]][,1])
colnames(noCut)<-names(resultT)
barplot(noCut,col=col, xlab=names(resultT), xlim=c(0,10),legend=rownames(noCut),las=2, main="fold change cutoff: 1")
for(i in 1:6){
    res<-resultT[[i]]
    barplot(res,xlab=fc,col=col, main=tests[i,3],las=2)
}
for(i in 1:6){
    res<-resultT[[i]][1:5,]
    barplot(res,xlab=fc,col=col[1:5], main=tests[i,3],las=2)
}
dev.off()

