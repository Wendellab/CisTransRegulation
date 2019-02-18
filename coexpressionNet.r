#####################################################
################# Input Discription #################
#####################################################

# 'totalReadCounts.txt' -  counts of total gene expression, genes (row) X 26 RNA-seq libraries (column)
# 'MaxxaAlleleReadCounts.txt' - for 12 F1 libraries (column), allele-specific read count assigned to Maxxa alleles were presented for these genes contain diagnostic SNPs.
# 'TX2094AlleleReadCounts.txt' - same as above but for TX2094 allele.
# 'TM1.nuclNmt.transcript.len.txt' - gene length for all genes


#####################################################
########## Finalize the sample inclusion #############
#####################################################

setwd("/home/hugj2006/cistrans/Ranalysis/coexpressionNet/")
# load total counts from bam idxstat
load("~/jfw-lab/Projects/AD1_domestication_cistrans/bowtie2_mapping/Ranalysis/counts.rdata")->l
l
# "info"       "total"      "count.rld"  "count.log2"
select<-which(info$genome=="Maxxa"|info$genome=="TX2094"|info$genome=="Yuc")
count<-total[,select]
coldata<-info[select,]
coldata$genome <- gsub("Yuc", "TX2094",coldata$genome)
xtabs(~genome+dpa,data=coldata)

## deseq normalization
library(DESeq2)
dds <- DESeqDataSetFromMatrix( countData = count, colData = coldata, design = ~new)
# rlog transformation, note that with defalt blind=TRUE, design is actually ~1
rld <- rlog(dds)
count.rld <- as.data.frame(assay(rld))
names(count.rld) == names(count)

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
    tipCol <- factor(color)
    levels(tipCol) <- hue_pal()(nlevels(tipCol))
    plot(tre,tip.col =as.character(tipCol),  type="unrooted",cex=0.6, no.margin=TRUE)
    plot(tre,tip.col =as.character(tipCol), type="fan",cex=0.6, no.margin=TRUE)
    dev.off()   }

# plots
plotGrouping(count.rld, shape=coldata$genome, color=coldata$dpa, tip=coldata$new, text=coldata$source, save = "plotGrouping.rld.pdf")
plotGrouping(count.log2, shape=coldata$genome, color=coldata$dpa, tip=coldata$new,text=coldata$source, save = "plotGrouping.log2.pdf")
plotGrouping(rpkm.log2, shape=coldata$genome, color=coldata$dpa, tip=coldata$new,text=coldata$source, save = "plotGrouping.log2rpkm.pdf")

## Samples to exclude
# TX2094-yb-20dpa-2 appears to be 10 dpa
# TX2094-byu-20dpa-1 appear to be 10 dpa
# Yuc-jj-20dpa-2 not cluster with others
# Yuc-simon-10dpa-1 appears to be 20 dpa
coldata$type ="PE"
coldata$type[coldata$source %in% c("byu","mj","simon")]="SE"
coldata$exclude =F
coldata$exclude[coldata$new %in% c("TX2094-yb-20dpa-2","TX2094-byu-20dpa-1","Yuc-jj-20dpa-2","Yuc-simon-10dpa-1")]=T

final.rld=count.rld[,coldata$exclude==F]
coldataF = coldata[coldata$exclude==F,]
plotGrouping(final.rld, shape=coldataF$genome, color=paste(coldataF$dpa,coldataF$type), tip=coldataF$new, text=coldataF$source, save = "plotGrouping.rld.final.pdf")

save(count, coldata, count.rld, rpkm.log2,count.log2, file="input.rdata")

####################
## Clust Analysis ##
####################

load("input.rdata")
# prepare input files for clust http://clust.baselabujamous.com
# count table
names(count) = coldata$new
maxxa = count[,coldata$exclude==F & coldata$genome=="Maxxa" ]
tx2094 = count[,coldata$exclude==F & coldata$genome=="TX2094" ]
system("mkdir data")
write.table(maxxa, file="data/Maxxa.txt", sep="\t",quote=F)
write.table(tx2094, file="data/TX2094.txt", sep="\t",quote=F)
# structure
coldataF = coldata[coldata$exclude==F,]
s<-aggregate(coldataF[,"new"], by=list(paste(coldataF$genome,coldataF$dpa,coldataF$type)), print)
ss<-data.frame(s$Group.1,sapply(s$x,function(x)paste(as.character(unlist(x)),collapse=", ")))
write.table(ss, file="repStruc.txt", col.names=FALSE,row.names=FALSE, sep="\t",quote=F)
# fix order to be: 5dpaSE, 10dpaSE, 10dpaPE, 15dpaSE, 20dpaSE
# add "ID" to first rows

### run Clust
# module load python/2.7.15-ief5zfp
# module load py-pip
# pip install --user clust
# clust data/ -r repStrucM.txt -n 101 3 4 -fil-v 25

# co-expressed gene clusters
x<-read.table("coexpressionNet/Clust_Results_30_Nov_18/Clusters_Objects.tsv",header=TRUE,sep="\t")
x<-x[-1,]
library(reshape2)
y=melt(x,measure.vars=names(x))
y=y[y$value!="",]
y$cluster = gsub("[.].*","",y$variable)
clust=y

####################
## WGCNA Analysis ##
####################

library(WGCNA)
library(RColorBrewer)
library(flashClust);
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);
allowWGCNAThreads(n=10);
load("input.rdata")

###### prep multiExpr
names(count.rld) = coldata$new
maxxa = count.rld[,coldata$exclude==F & coldata$genome=="Maxxa" ]
tx2094 = count.rld[,coldata$exclude==F & coldata$genome=="TX2094" ]
shortLabels = c("TX2094","Maxxa")
nSets=2
multiExpr = vector(mode = "list", length = 2)
multiExpr[[1]] = list(data = t(tx2094))
multiExpr[[2]] = list(data = t(maxxa))
# Check that the data has the correct format for many functions operating on multiple sets:
print(checkSets(multiExpr))
# Check that all genes and samples have sufficiently low numbers of missing values.
filterMultiExpr<-function(multiExpr,nSets)
{
    gsg = goodSamplesGenesMS(multiExpr, verbose = 3);
    print(gsg$allOK)
    if (!gsg$allOK)
    {
        # Print information about the removed genes:
        if (sum(!gsg$goodGenes) > 0)
        printFlush(paste("Removing genes:", paste(names(multiExpr[[1]]$data)[!gsg$goodGenes], collapse = ", ")))
        for (set in 1:nSets)
        {
            if (sum(!gsg$goodSamples[[set]]))
            printFlush(paste("In set", setLabels[set], "removing samples",paste(rownames(multiExpr[[set]]$data)[!gsg$goodSamples[[set]]], collapse = ", ")))
            # Remove the offending genes and samples
            multiExpr[[set]]$data = multiExpr[[set]]$data[gsg$goodSamples[[set]], gsg$goodGenes];
        }
    }
    # Excluding genes from the calculation due to too many missing samples or zero variance.
    # Update exprSize
    print(checkSets(multiExpr))
    return(multiExpr)
}
multiExpr<-filterMultiExpr(multiExpr,nSets)
print(checkSets(multiExpr)$nGenes) #63665

###### choose soft power
# Choose a set of soft-thresholding powers, consider three type of adjacnecy tables, although I am going to use "signed" network for this analysis
# types<-c("unsigned", "signed", "signed hybrid")
type <- "signed"
powers = c(c(1:10), seq(from = 12, to=40, by=2))
# Initialize a list to hold the results of scale-free analysis
powerTables = vector(mode = "list", length = nSets);
# Call the network topology analysis function for each set in turn
for(set in 1:nSets){
    powerTables[[set]] = list(data = pickSoftThreshold(multiExpr[[set]]$data, powerVector=powers, verbose = 2, networkType = type, blockSize=10000)[[2]])      }
collectGarbage()
# Plot the results:
colors=brewer.pal(5,"Set1")
# Will plot these columns of the returned scale free analysis tables
plotCols = c(2,5,6,7)
colNames = c("Scale Free Topology Model Fit", "Mean connectivity", "Median connectivity", "Max connectivity");
# Get the minima and maxima of the plotted points
ylim = matrix(NA, nrow = 2, ncol = 4);
for (set in 1:nSets)
{
    for (col in 1:length(plotCols))
    {
        ylim[1, col] = min(ylim[1, col], powerTables[[set]]$data[, plotCols[col]], na.rm = TRUE);
        ylim[2, col] = max(ylim[2, col], powerTables[[set]]$data[, plotCols[col]], na.rm = TRUE);
    }
}
# Plot the quantities in the chosen columns vs. the soft thresholding power
# sizeGrWindow(8, 6)
pdf("choosePower.pdf")
par(mfcol = c(2,2));
par(mar = c(4.2, 4.2 , 2.2, 0.5))
cex1 = 0.7;
for (col in 1:length(plotCols)) for (set in 1:nSets)
{
    if (set==1)
    {
        plot(powerTables[[set]]$data[,1], -sign(powerTables[[set]]$data[,3])*powerTables[[set]]$data[,2], xlab="Soft Threshold (power)",ylab=colNames[col],type="n", ylim = ylim[, col], main = colNames[col]);
        addGrid()
    }
    if (col==1)
    {
        text(powerTables[[set]]$data[,1], -sign(powerTables[[set]]$data[,3])*powerTables[[set]]$data[,2],
        labels=powers,cex=cex1,col=colors[set]);
    } else
    text(powerTables[[set]]$data[,1], powerTables[[set]]$data[,plotCols[col]], labels=powers,cex=cex1,col=colors[set]);
    if (col==1)
    {
        legend("bottomright", legend = shortLabels, col = colors, pch = 20) ;
    } else
    legend("topright", legend = shortLabels, col = colors, pch = 20) ;
}
dev.off()
# inspect figure p=20
powerEach = 20

###### calculate individual TOMs
print("###### Calculate individual TOMs:")
iTOMs = blockwiseIndividualTOMs(
# Input data
multiExpr,
# Data checking options
checkMissingData = TRUE,
# Options for splitting data into blocks
blocks = NULL,
randomSeed = 12345,
maxBlockSize =  20000,  # as individual network # 5000 for 4G memory, 20000 for 16G, 30000 for 32 G
# Network construction arguments: correlation options, use bicor instead of default pearson
corType = "pearson",
# Adjacency and topology overlap function options
power = powerEach, networkType = "signed", TOMType = "signed",
# Save individual TOMs?
saveTOMs = TRUE,
individualTOMFileNames = "iTOM-%N-block.%b.RData"  )

###### calculate consensus modules
print("###### Construct consensus networks:")
cnet = blockwiseConsensusModules(
# Input data
multiExpr,
# Data checking options
checkMissingData = TRUE,
# Options for splitting data into blocks
blocks = NULL,
randomSeed = 12345,
maxBlockSize =  20000,  # as individual network # 5000 for 4G memory, 20000 for 16G, 30000 for 32 G
# Network construction arguments: correlation options, use bicor instead of default pearson
corType = "pearson",
# Adjacency and topology overlap function options
power = powerEach, networkType = "signed", TOMType = "signed",
# load previous TOMs
individualTOMInfo = iTOMs,
# Saving the consensus TOM
saveConsensusTOMs = TRUE,
consensusTOMFileNames = "cTOM-Block%b.RData",
# Basic tree cut options
deepSplit = 2,  #default, known to reasonable
minModuleSize = 100, #default 20, use 30 for transcriptome, or ncol(subDat)/2
pamStage = TRUE, pamRespectsDendro = TRUE, #default, known to reasonable
# Thredhold to merge modules: a height cut of 0.25 corresponding to correlation of 0.75
mergeCutHeight = 0.25,
# others
reassignThreshold = 0,
numericLabels = TRUE,
verbose = 3)

###### calculate individual modules
print("###### Construct individual networks:")
# stupid blockwiseModules only load TOM rdata file with "TOM", not "tomDS"
tomFiles<-grep("iTOM",list.files(), value=TRUE)
for(fl in tomFiles)
{
    load(fl)
    TOM<-tomDS
    save(TOM, file=fl)
}
rm(TOM)
collectGarbage()
for(i in 1:nSets)
{
    inet = blockwiseModules(
    # Input data
    multiExpr[[i]]$data,
    # Data checking options
    checkMissingData = TRUE,
    # Options for splitting data into blocks
    blocks =  iTOMs$blocks,
    #randomSeed = 12345,
    #maxBlockSize =  500,  # as individual network # 5000 for 4G memory, 20000 for 16G, 30000 for 32 G
    # Network construction arguments: correlation options, use bicor instead of default pearson
    corType = "pearson",
    # Adjacency and topology overlap function options
    power = powerEach, networkType = "signed", TOMType = "signed",
    # load previous TOMs
    loadTOM = TRUE,
    saveTOMFileBase = paste0("iTOM-",i),
    # Basic tree cut options
    deepSplit = 2,  #default, known to reasonable
    minModuleSize = 100, #default 20, use 30 for transcriptome, or ncol(subDat)/2
    pamStage = TRUE, pamRespectsDendro = TRUE, #default, known to reasonable
    # Thredhold to merge modules: a height cut of 0.25 corresponding to correlation of 0.75
    mergeCutHeight = 0.25,
    # others
    reassignThreshold = 0,
    numericLabels = TRUE,
    verbose = 3)
    assign(paste0("inet",i), inet)
}

###### Connectivity hypothesis
## TX2094
adj= adjacency(multiExpr[[1]]$data, type = "signed", power=20) # build
# whole network connectivity
# netk = rowSums(adj)-1
# kTotal, the within module connectivity kWithin, kOut=kTotal-kWithin, and kDiff=kIn-kOut=2*kIN-kTotal
k=intramodularConnectivity(adj, inet1$colors)
# module membership
kme = signedKME(multiExpr[[1]]$data, inet1$MEs)
k.tx2094 = cbind(k,kme)
##Maxxa
adj= adjacency(multiExpr[[2]]$data, type = "signed", power=20) # build
k=intramodularConnectivity(adj, inet2$colors)
kme = signedKME(multiExpr[[2]]$data, inet2$MEs)
k.maxxa = cbind(k,kme)

save(list=c( "multiExpr","shortLabels","nSets", "iTOMs", "cnet", "clust",grep("inet.",ls(), value=TRUE),"k.tx2094","k.maxxa"), file = "buildNetwork.rdata")
collectGarbage()

###### preservation test Maxxa against TX2094
# The number of permutations drives the computation time of the module preservation function. For a publication use 200 permutations.
# But for brevity, let's use a small number
nPermutations1=200
# Set it to a low number (e.g. 3) if only the medianRank statistic and other observed statistics are needed.
# Permutations are only needed for calculating Zsummary and other permutation test statistics.
# set the random seed of the permutation test analysis
set.seed(1)
system.time({
    mp = modulePreservation(multiExpr, multiColor=list(inet1$colors,inet2$colors), networkType="signed", referenceNetworks = 1, testNetworks=2, nPermutations = nPermutations1,randomSeed = 15423, quickCor = 0, verbose = 3)
})
# Save the results of the module preservation analysis
save(mp, file = "preservationZ&Medianrank.rdata")
# view MAxxa vs TX2094 median rank
mp$preservation$observed[[1]][[2]] ->mr
mr[order(mr$medianRank.pres,decreasing=TRUE),1:2]
# If needed, reload the data:
load(file = "preservationZ&Medianrank.rdata")
pdf("modulePreservation.pdf")
# specify the reference and the test networks
ref=1; test = 2
Obs.PreservationStats= mp$preservation$observed[[ref]][[test]]
Z.PreservationStats=mp$preservation$Z[[ref]][[test]]
# Look at the observed preservation statistics
Obs.PreservationStats
# Z statistics from the permutation test analysis
Z.PreservationStats
# Let us now visualize the data.
modIDs = rownames(Obs.PreservationStats)
modColors=labels2colors(order(as.numeric(modIDs) )-1 )
moduleSize = Obs.PreservationStats$moduleSize
# we will omit the grey module (background genes)
# and the gold module (random sample of genes)
selectModules = !(modColors %in% c("grey", "gold"))
# Text labels for points
point.label = modIDs[selectModules]
# Composite preservation statistics
medianRank=Obs.PreservationStats$medianRank.pres
Zsummary=Z.PreservationStats$Zsummary.pres
par(mfrow=c(1,2),mar = c(4.5,4.5,2.5,1))
# plot medianRank versus module size
plot(moduleSize[selectModules],medianRank[selectModules],col=1, bg=modColors[selectModules],
pch = 21,main=paste("medianRank -",shortLabels[ref], "vs",shortLabels[test]),
cex = 2, ylab ="medianRank",xlab="Module size", log="x")
labelPoints(moduleSize[selectModules],medianRank[selectModules],point.label,cex=1,offs=0.03)
# plot Zsummary versus module size
plot(moduleSize[selectModules],Zsummary[selectModules], col = 1, bg=modColors[selectModules],pch = 21,
main=paste("Zsummary -",shortLabels[ref], "vs",shortLabels[test]),
cex=2,ylab ="Zsummary", xlab = "Module size", log = "x")
labelPoints(moduleSize[selectModules],Zsummary[selectModules],point.label,cex=1,offs=0.03)
# Add threshold lines for Zsummary
abline(h=0); abline(h=2, col = "blue", lty = 2); abline(h=10, col = "red", lty = 2)
dev.off()

###### collect modules and check correspondence
clusters = data.frame(ID=colnames(multiExpr[[1]]$data), tx2094 =inet1$colors, maxxa=inet2$colors, concensus =cnet$colors)
rownames(clusters)=clusters$ID
# co-expressed gene clusters from Clust
x<-read.table("coexpressionNet/Clust_Results_30_Nov_18/Clusters_Objects.tsv",header=TRUE,sep="\t")
x<-x[-1,]
library(reshape2)
y=melt(x,measure.vars=names(x))
y=y[y$value!="",]
y$cluster = gsub("[.].*","",y$variable)
# combine
clusters$clust = "NC" #not clusterred
clusters[y$value,"clust"]= y$cluster
save(clusters,file ="../coexpressionNet.rdata")


#####################
## Genie3 Analysis ##
#####################

# GENIE3
# https://bioconductor.org/packages/release/bioc/html/GENIE3.html
# source("https://bioconductor.org/biocLite.R")
# biocLite("GENIE3")
if (!requireNamespace("BiocManager", quietly = TRUE))
install.packages("BiocManager")
BiocManager::install("GENIE3", version = "3.8")

## transcription factor prediction by PlantTFDB v4.0 http://planttfdb.cbi.pku.edu.cn/prediction.php
# Input 87800 sequences of "Ghirsutum_458_v1.1.protein.fa"
# Check box "link best hit in Arabidopsis"
# Output 6341 TFs identified "TFpredication.txt": 4771 are primary transcripts
# 2019-1-30
tfs <- read.table("coexpressionNet/TFs/TFpredication.txt",sep="\t")
names(tfs) = c("ID","Family","bestHitsAt","Evalue","description")
dim(tfs)
dim(tfs<-tfs[grep("[.]1[.]p",tfs$ID),])
tfs$ID = gsub("[.]p","",tfs$ID)

## TFs being RD
load("RDgenes.rdata")
rd$ID=as.character(rd$ID)
keyTFs <- tfs$ID[tfs$ID %in% rd$ID] # TFs in 27816
length(keyTFs)# 1518
table(rd$rd[rd$ID %in% keyTFs & rd$rd!="no divergence"])
# cis and trans      cis only    trans only
#         26            12            15
# a total of 53 rd TFs

## run GENIE3
library(GENIE3)
set.seed(123) # For reproducibility of results
######### TX2094
exprMatr <- t(multiExpr[[1]]$data)  # tx2094
exprMatr <- exprMatr[as.character(rd$ID),] # restrict to 27816 genes to speed up
weightMat <- GENIE3(exprMatr, regulators=keyTFs)
dim(weightMat)
#  1518 27816
linkList <- getLinkList(weightMat)
dim(linkList)
# 42217427        3
head(linkList)
#   regulatoryGene         targetGene     weight
# 1 Gohir.D08G035600.1 Gohir.A01G018300.1 0.04243262
# 2 Gohir.A11G137800.1 Gohir.A13G147500.1 0.04060390
# 3 Gohir.A11G307000.1 Gohir.D08G143800.1 0.03737583
quantile(linkList$weight, 1:10/10)
#         10%           20%           30%           40%           50%
# 0.00001555126 0.00005237248 0.00010098484 0.00016447188 0.00025094481
#          60%           70%           80%           90%          100%
# 0.00037575590 0.00056881689 0.00088935127 0.00169343318 0.04243262159
linksT = linkList
######## Maxxa
exprMatr <- t(multiExpr[[2]]$data)  # tx2094
exprMatr <- exprMatr[as.character(rd$ID),] # restrict to 27816 genes to speed up
weightMat <- GENIE3(exprMatr, regulators=keyTFs)
dim(weightMat)
#  1518 27816
linkList <- getLinkList(weightMat)
dim(linkList)
# 42217427        3
head(linkList)
quantile(linkList$weight, 1:10/10)
#           10%           20%           30%           40%           50%
# 0.00002644882 0.00006696148 0.00011577627 0.00017650869 0.00025598398
#           60%           70%           80%           90%          100%
# 0.00036693821 0.00053497221 0.00082149137 0.00157235971 0.04112548261
linksM = linkList

# filter predictions weight
weight.cutoff = 0.005
dim(tx2094 <- linksT[linksT$weight>weight.cutoff,]) # 639090
dim(maxxa  <- linksM[linksM$weight>weight.cutoff,]) # 744458
save(maxxa, tx2094, file="genie3.rdata")

# Export genie3 network edges including only RD genes to cytoscape
dim(mm <- maxxa[maxxa$targetGene%in%rd$ID[rd$rd!="no divergence"] & maxxa$regulatoryGene%in%rdTFs$ID,] ) #2055
dim(tt <- tx2094[tx2094$targetGene%in%rd$ID[rd$rd!="no divergence"] & tx2094$regulatoryGene%in%rdTFs$ID,] ) #1680
write.table(mm, row.names=FALSE, sep="\t",file="maxxa_rds_genie3_edges.txt",quote=FALSE)
write.table(tt, row.names=FALSE, sep="\t",file="tx2094_rds_genie3_edges.txt",quote=FALSE)
# Preppare node table with description
load("funcEnrich.rdata")
node=desc[rd$ID,]
node$rd=rd$rd
names(node)[1]="ID"
# is TFs or not
node$TF =node$ID %in% tfs$ID
node=merge(node, tfs[,c(1,2,3,5)],by="ID",all.x=TRUE)
# in and out degree from full network of 27816 genes
node$ID=as.character(node$ID)
node$maxxaIn = sapply(node$ID,function(x)length(which(maxxa$targetGene==x)))
node$maxxaOut = sapply(node$ID,function(x)length(which(maxxa$regulatoryGene==x)))
node$tx2094In = sapply(node$ID,function(x)length(which(tx2094$targetGene==x)))
node$tx2094Out = sapply(node$ID,function(x)length(which(tx2094$regulatoryGene==x)))
# in qtl or not
x<-read.table("/work/LAS/jfw-lab/corrinne/QTLpaper/fiber.qtl.gene.results", sep="\t")
g<- unique(paste0(gsub(".*Name=","",x$V12),".1"))
length(g)
node$qtl = node$ID %in% g
# in JJ's list of cell wall and cytoskeleton genes or not
jj=read.table("/work/LAS/jfw-lab/erdostal/JJgorai.txt")
dim(jj) #703
length(jj<-as.character(unique(jj$V1))) #402
load("homoeoPairs.rdata")
ogP$containJJ <- sapply(ogP$Gorai, function(x) length(intersect(jj, unlist(strsplit(x,split=",")))))
table(ogP$containJJ) # 0-22087; 1-307
table(ogP$Gorai %in% jj) #303 are in quadruplet
cw<-ogP[ogP$containJJ>0,]
gohir=c(cw$Gohir.A,cw$Gohir.D) # project JJ Gorai to Gohir
length(gohir <- unique(gohir[!is.na(gohir)])) #614
node$JJ = node$ID %in% gohir
# in which module or clust
load("coexpressionNet.rdata")
node <- merge(node,clusters,by="ID", all.x=TRUE)
node$maxxaC= col2hex(labels2colors(node$maxxa))
node$tx2094C= col2hex(labels2colors(node$tx2094))
node$concensusC= col2hex(labels2colors(node$concensus))
# extra annotation for rd TFs, giving rank, GO annotation of target genes
select=which(node$TF==TRUE & node$rd!="no divergence")
rdTFs = data.frame(ID= node$ID[select], rank = order(node$maxxaOut[select]+node$maxxaOut[select], decreasing=TRUE)) # rank by average out degree
node <- merge(node,rdTFs,by="ID", all.x=TRUE)
# add gorai and OG
map = data.frame(gorai =c(ogP$Gorai,ogP$Gorai), ID=c(ogP$Gohir.A,ogP$Gohir.D), OG=c(ogP$OG,ogP$OG))
dim(node<-merge(node,map,by="ID",all.x=TRUE))
# save
write.table(node, row.names=FALSE, sep="\t",file="nodes.txt",quote=FALSE)



# GO enrichment of rdTFs targets
rdTFs$maxxaTGsGO =NA
rdTFs$tx2094TGsGO =NA
refL = node$ID # 27816
library(clusterProfiler)
library(GO.db)
for(i in 1:nrow(rdTFs)){
    tf = as.character(rdTFs$ID)[i]
    # maxxa and tx2094
    print(i)
    geneL = list(maxxa=as.character(maxxa$targetGene[maxxa$regulatoryGene==tf]),tx2094=as.character(tx2094$targetGene[tx2094$regulatoryGene==tf]))
    result=GOcompare(geneL,refL,filename=paste0("rdTFs/GOcompare.",tf))
    rdTFs$maxxaTGsGO[rdTFs$ID==tf] = paste(result$Description[result$Cluster=="maxxa"][1:3],collapse="; ")
    rdTFs$tx2094TGsGO[rdTFs$ID==tf] = paste(result$Description[result$Cluster=="tx2094"][1:3],collapse="; ")
}

## what TFs have targets enriched with trans
TG.m = as.data.frame(table(maxxa$regulatoryGene))
transTG.m = as.data.frame(table(maxxa$regulatoryGene[maxxa$targetGene %in% rd$ID[grep("trans", rd$rd)]]))
unique(TG.m$Var1==transTG.m$Var1)
m=data.frame(ID=TG.m$Var1,TG=TG.m$Freq,transTG=transTG.m$Freq)
p=sum(m$transTG)/sum(m$TG)
m$pval=apply(m[,2:3],1,function(x){binom.test(x[2], n=x[1], p, alternative="greater")$p.value})
m$qval =p.adjust(m$pval,"BH")


# Maxxa vs TX2094: do in and out degree rank differ
# all nodes
with(node,wilcox.test(maxxaOut,tx2094Out,paired=TRUE))
with(node,wilcox.test(maxxaIn,tx2094In,paired=TRUE))
# TFs
with(node[node$TF==TRUE,],wilcox.test(maxxaOut,tx2094Out,paired=TRUE))
with(node[node$TF==TRUE,],wilcox.test(maxxaIn,tx2094In,paired=TRUE))
# RDs
with(node[node$rd!="no divergence",],wilcox.test(maxxaOut,tx2094Out,paired=TRUE))
with(node[node$rd!="no divergence",],wilcox.test(maxxaIn,tx2094In,paired=TRUE))
# RD TFs
with(node[node$TF==TRUE & node$rd!="no divergence",],wilcox.test(maxxaOut,tx2094Out,paired=TRUE))
with(node[node$TF==TRUE & node$rd!="no divergence",],wilcox.test(maxxaIn,tx2094In,paired=TRUE))

## TFs that exhibit regulatory divergence are more likely to cause regulatory changes of their downstream target genes?
rgn=maxxa
TG = as.data.frame(table(rgn$regulatoryGene))
rdTG = as.data.frame(table(rgn$regulatoryGene[maxxa$targetGene %in% rd$ID[rd$rd!="no divergence"]]))
unique(TG$Var1==rdTG$Var1)
df=data.frame(ID=TG$Var1,TG=TG$Freq,rdTG=rdTG$Freq, rd=rd[as.character(TG$Var1),"rd"])

 wilcox.test(df$rdTG[df$rd=="no divergence"], df$rdTG[df$rd!="no divergence"])
 wilcox.test(df$TG[df$rd=="no divergence"], df$TG[df$rd!="no divergence"])




##################
## RD x Network ##
##################
source("utilities.r")
load("coexpressionNet/buildNetwork.rdata")
load("coexpressionNet.rdata")
load("RDgenes.rdata")

## subnetwork density x RD
getNetConcept=function(expr,power=20) # density,Heterogeneity
{
    adj= adjacency(expr, type = "signed", power=20)
    # fnc = fundamentalNetworkConcepts(adj, GS = NULL)
    Size = dim(adj)[1]
    Connectivity = apply(adj, 2, sum)-1
    Density = sum(Connectivity)/(Size * (Size - 1))
    Centralization = Size * (max(Connectivity) - mean(Connectivity))/((Size - 1) * (Size - 2))
    Heterogeneity = sqrt(Size * sum(Connectivity^2)/sum(Connectivity)^2 -1)
    res=c(Density,Centralization, Heterogeneity)
    names(res)=c("Density","Centralization", "Heterogeneity")
    return(res)
}
getDensity=function(expr,power=20) # density,Heterogeneity
{
    adj= adjacency(expr, type = "signed", power=20)
    # fnc = fundamentalNetworkConcepts(adj, GS = NULL)
    Size = dim(adj)[1]
    Connectivity = apply(adj, 2, sum)-1
    Density = sum(Connectivity)/(Size * (Size - 1))
    return(Density)
}

## are RD genes more connected than a random list of genes of the same size?
expr= multiExpr[[1]]$data # tx2094
ids = as.character(rd$ID[rd$rd!="no divergence"])
no=length(ids)
den=getDensity(expr[,ids],power=20)
for(i in 1:100){
    rids=as.character(sample(rd$ID,no))
    den=c(den, getDensity(expr[,rids],power=20) )
}
t<-den
t[1];quantile(t[-1])
# [1] 0.01897275
# 0%         25%         50%         75%        100%
# 0.008143221 0.008862144 0.009186852 0.009447498 0.010197187
expr= multiExpr[[2]]$data # maxxa
ids = as.character(rd$ID[rd$rd!="no divergence"])
no=length(ids)
den=getDensity(expr[,ids],power=20)
for(i in 1:100){
    rids=as.character(sample(rd$ID,no))
    den=c(den, getDensity(expr[,rids],power=20) )
}
m[1];quantile(m[-1])
# [1] 0.009868543
# 0%         25%         50%         75%        100%
# 0.005514290 0.005968848 0.006170588 0.006376306 0.006746608

## are trans affected genes more intercent
clab = c("all","non-RD","RD","cis_only","cis&trans","trans_only")
expr1= multiExpr[[1]]$data # tx2094
expr2= multiExpr[[2]]$data # maxxa
tx2094 = c(getDensity(expr1[,as.character(rd$ID)], power=20), getDensity(expr1[,as.character(rd$ID[rd$rd=="no divergence"])],power=20), getDensity(expr1[,as.character(rd$ID[rd$rd!="no divergence"])],power=20), getDensity(expr1[,as.character(rd$ID[rd$rd=="cis only"])],power=20), getDensity(expr1[,as.character(rd$ID[rd$rd=="cis and trans"])],power=20), getDensity(expr1[,as.character(rd$ID[rd$rd=="trans only"])],power=20) )
maxxa = c(getDensity(expr2[,as.character(rd$ID)], power=20), getDensity(expr2[,as.character(rd$ID[rd$rd=="no divergence"])],power=20), getDensity(expr1[,as.character(rd$ID[rd$rd!="no divergence"])],power=20), getDensity(expr2[,as.character(rd$ID[rd$rd=="cis only"])],power=20), getDensity(expr2[,as.character(rd$ID[rd$rd=="cis and trans"])],power=20), getDensity(expr2[,as.character(rd$ID[rd$rd=="trans only"])],power=20) )
density = data.frame(subnetwork = clab,tx2094,maxxa)

# 0.008742397 0.018988961 0.018073715 0.031449437, only trans only above overall RD desnisty
expr= multiExpr[[2]]$data # maxxa
c(getDensity(expr[,as.character(rd$ID[rd$rd=="no divergence"])],power=20), getDensity(expr[,as.character(rd$ID[rd$rd=="cis only"])],power=20), getDensity(expr[,as.character(rd$ID[rd$rd=="cis and trans"])],power=20), getDensity(expr[,as.character(rd$ID[rd$rd=="trans only"])],power=20) )
# 0.005955358 0.008200147 0.011468430 0.012620720

#

expr= multiExpr[[2]]$data #maxxa
getNetConcept(expr[,as.character(rd$ID[rd$rd=="trans only"])],power=20) # 0.01262072     0.05278882     1.14557756
getNetConcept(expr[,as.character(rd$ID[rd$rd=="cis only"])],power=20) # 0.008200147    0.028780244    1.030490529
table(rd$rd)
# cis and trans      cis only no divergence    trans only
#       841           513         26161           301
den=data.frame(nrd=getNetConcept(expr[,as.character(rd$ID[rd$rd=="no divergence"])],power=20)["Density"],
cisO=getNetConcept(expr[,as.character(rd$ID[rd$rd=="cis only"])],power=20)["Density"],
cisNtrans=getNetConcept(expr[,as.character(rd$ID[rd$rd=="cis and trans"])],power=20)["Density"],
trans=getNetConcept(expr[,as.character(rd$ID[rd$rd=="trans only"])],power=20)["Density"] )
den
#         nrd        cisO  cisNtrans      trans
# 0.005955358 0.008200147 0.01146843 0.01262072
# permutation
for(i in 2:101){
    ids=as.character(sample(rd$ID[rd$rd=="no divergence"],301))
    den[i,"nrd"]=getNetConcept(expr[,ids],power=20)["Density"]
}
for(i in 2:101){
    ids=as.character(sample(rd$ID[rd$rd=="cis only"],301))
    den[i,"cisO"]=getNetConcept(expr[,ids],power=20)["Density"]
}
for(i in 2:101){
    ids=as.character(sample(rd$ID[rd$rd=="cis and trans"],301))
    den[i,"cisNtrans"]=getNetConcept(expr[,ids],power=20)["Density"]
}
boxplot(den)
points(1:4,den[1,],col="blue",pch=19)
dev.off()

## clusters x RD
node<-read.table("nodes.txt",header=TRUE,sep="\t")
pdf("network_RD.pdf")
plotFisherCorr(node$maxxa,node$rd, title="Maxxa WGCNA", fisher="greater")
plotFisherCorr(node$tx2094, node$rd, title="TX2094 WGCNA", fisher="greater")
plotFisherCorr(node$concensus,node$rd,  title="M-T consensus WGCNA", fisher="greater")
plotFisherCorr(node$clust,node$rd,  title="Clust", fisher="greater")
plotFisherCorr(node$JJ,node$rd,  title="JJ genes", fisher="greater")
plotFisherCorr(node$qtl,node$rd,  title="QTL", fisher="greater")
plotFisherCorr(node$TF,node$rd,  title="TFs", fisher="greater")
dev.off()

plotFisherCorr(clusters$module.tx2094,clusters$module.maxxa, title="", fisher="greater")
plotCorrespondence(clusters$module.tx2094,clusters$module.maxxa,file=NULL,title="", textplot=FALSE,corrplot=FALSE,fisherplot=TRUE,fisher="greater", mai=c(1.02, 0.82,0.82,0.42))
plotCorrespondence(clusters$module.tx2094,clusters$clust,file=NULL,title="", textplot=FALSE,corrplot=FALSE,fisherplot=TRUE,fisher="greater", mai=c(1.02, 0.82,0.82,0.42))
plotCorrespondence(clusters$module.maxxa,clusters$clust,file=NULL,title="", textplot=FALSE,corrplot=FALSE,fisherplot=TRUE,fisher="greater", mai=c(1.02, 0.82,0.82,0.42))
table(rd$rd,cl$clust)
plotCorrespondence(rd$rd,cl$clust,title="RD x clust", textplot=FALSE,corrplot=FALSE,fisherplot=TRUE,fisher="greater", mai=c(1.02, 0.82,0.82,0.42))
---book
