## Compare regulatory pattern between homoeologs
setwd("~/Dropbox/cistrans_Fiber/NewAnalysis")
library(gplots)
library(ggplot2)
load("cistrans.rdata")
options(scipen=999)

# all genes with allelic SNPs
ids=rownames(res10)
countADZ=function(x) {table(gsub("[0-9].*","",gsub("Gohir.1|Gohir.|_.*","",x)))}
countADZ(ids)
#     A     D    mt     Z
# 13462 14060     9   285

# load Justin's ortholog group results, v3
# smb://lss.its.iastate.edu/gluster-lss/research/jfw-lab/Projects/cis_trans_regulation/DDtAt.orthogroups.v3.csv
og<-read.table("DDtAt.orthogroups.v3.txt", sep="\t", header=TRUE)
# count how many At, Dt and D5 genes in each group
countGene <-function(x)
{
    length(unlist(strsplit(as.character(x),",")))
}
ogN<-as.data.frame(apply(og,c(1,2),countGene))
# how many categories, and how many OGs in each category
library(plyr)
unique(ogN)
z<-ddply(as.data.frame(ogN),.(Gohir.A,Gohir.D,Gorai),nrow)
# extract OGs with 1:1 At and Dt
pairs <- which(ogN$Gohir.A==1 & ogN$Gohir.D==1)
paste0("A total of ", length(pairs), " out of ",nrow(og)," OGs contain At/Dt homoeolog pairs.")
# "A total of 22394 out of 28843 OGs contain At/Dt homoeolog pairs."
ogP <-og[pairs,1:4]
ogP_snps<-ogP[ogP$Gohir.A %in% ids & ogP$Gohir.D %in% ids, ]
paste0(nrow(ogP_snps), " homoeolog pair OGs were detected with SNPs between Maxxa and TX2094.")
names(ogP)[1]<-"OG"
write.table(ogP,file="homoeolog_pairs052218.txt",sep="\t", row.names=FALSE)

# "8036 homoeolog pair OGs were detected with SNPs between Maxxa and TX2094."

# v1 - A total of 23700 out of 29636 OGs contain At/Dt homoeolog pairs.
# 8266 homoeolog pair OGs were detected with SNPs between Maxxa and TX2094.

# v2 - A total of 22048 out of 28731 OGs contain At/Dt homoeolog pairs.
# 7896 homoeolog pair OGs were detected with SNPs between Maxxa and TX2094.

# v3 - A total of 22394 out of 28843 OGs contain At/Dt homoeolog pairs.
# "8036 homoeolog pair OGs were detected with SNPs between Maxxa and TX2094."


table(ids[grep("Gohir.A",ids)] %in% ogP_snps$Gohir.A)
# FALSE  TRUE
# 5426  8036
table(ids[grep("Gohir.D",ids)] %in% ogP_snps$Gohir.D)
# FALSE  TRUE
# 6024  8036
countADZ(ids[!(ids %in% ogP_snps$Gohir.D | ids %in% ogP_snps$Gohir.A)])
#    A    D   mt    Z
# 5426 6024    9  285

# label ids by OG
og_by_id <- c(as.character(ogP$OG),as.character(ogP$OG))
names(og_by_id) <-c(as.character(ogP$Gohir.A), as.character(ogP$Gohir.D))

res=data.frame(id=ids, OG = og_by_id[ids], OG.duplicated=NA, origin = gsub("[0-9].*","",gsub("Gohir.1|Gohir.|_.*","",ids) ))
res$OG.duplicated = duplicated(res$OG) | duplicated(res$OG, fromLast=TRUE) # false means solo homoeolog
res$OG.duplicated[is.na(res$OG)] = "unknown" # not clear about OG, at least not 1:1
table(res$origin, res$OG.duplicated)
#    FALSE TRUE unknown
# A   2874 8036    2552
# D   3029 8036    2995
# mt     0    0       9
# Z      0    0     285
res$homoeo.status <- res$OG.duplicated
res$homoeo.status[res$origin=="A" & res$OG.duplicated==TRUE] <- "inPair.A"
res$homoeo.status[res$origin=="A" & res$OG.duplicated==FALSE] <- "solo.A"
res$homoeo.status[res$origin=="D" & res$OG.duplicated==TRUE] <- "inPair.D"
res$homoeo.status[res$origin=="D" & res$OG.duplicated==FALSE] <- "solo.D"
table(res$homoeo.status)
# inPair.A inPair.D   solo.A   solo.D  unknown
#     8036     8036     2874     3029     5841
write.table(res,file="homoeolog_status042618.txt",sep="\t")

# examine cis tran categorization by gene homoeolog status
library(vcd)
l<-c("res10","res20","res.mt10","res.mt20","res.tm10","res.tm20")
pdf("plotAssoc.cistrans.vs.homoeoS.pdf")
for(i in l){
    x<-get(i)
    if(unique(res$id == rownames(x))){
        cis_trans_regulatory_category = gsub("[.].*","",x$category)
        homoeolog_origin_presence = res$homoeo.status
        tbl = table(cis_trans_regulatory_category, homoeolog_origin_presence)
        textplot(tbl); title(paste0(i,": p-value ", chisq.test(tbl,simulate.p.value = TRUE)$p.value))
        assoc(tbl, shade=TRUE, main = i)
        mosaic(tbl, shade=TRUE, main = i)
    }
}
dev.off()
pdf("plotAssoc.cistrans.vs.homoeoO.pdf")
for(i in l){
    x<-get(i)
    if(unique(res$id == rownames(x))){
        cis_trans_regulatory_category = gsub("[.].*","",x$category)
        homoeolog_origin = res$origin
        tbl = table(cis_trans_regulatory_category, homoeolog_origin)
        textplot(tbl); title(paste0(i,": p-value ", chisq.test(tbl,simulate.p.value = TRUE)$p.value))
        assoc(tbl, shade=TRUE, main = i)
        mosaic(tbl, shade=TRUE, main = i)
    }
}
dev.off()
pdf("plotAssoc.cistrans.vs.homoeoD.pdf")
for(i in l){
    x<-get(i)
    if(unique(res$id == rownames(x))){
        cis_trans_regulatory_category = gsub("[.].*","",x$category)
        homoeolog_duplicated = res$OG.duplicated
        tbl = table(cis_trans_regulatory_category, homoeolog_duplicated)
        textplot(tbl); title(paste0(i,": p-value ", chisq.test(tbl,simulate.p.value = TRUE)$p.value))
        assoc(tbl, shade=TRUE, main = i)
        mosaic(tbl, shade=TRUE, main = i)
    }
}
dev.off()
library(corrplot)
library(WGCNA)
pdf("plotAssoc.cistrans.AvsD.pdf")
for(i in l){
    x<-get(i)
    if(unique(res$id == rownames(x))){
        At = x[ogP_snps$Gohir.A,"category"]
        Dt = x[ogP_snps$Gohir.D,"category"]
        OG = as.character(ogP_snps$OG)
        
        tbl = table(At,Dt)
        textplot(tbl); title(paste0(i,": p-value ", chisq.test(tbl,simulate.p.value = TRUE)$p.value))
        assoc(tbl, shade=TRUE, main = i)
        mosaic(tbl, shade=TRUE, main = i)
        
        CountTbl = as.matrix(tbl)
        residualTable  = chisq.test(tbl,simulate.p.value = TRUE)$residuals
        corrplot(residualTable, p.mat=CountTbl, insig="p-value", sig.level=0,is.corr=FALSE)
        
        pTable = CountTbl
        rowCats = levels(factor(At))
        colCats = levels(factor(Dt))
        for(rc in rowCats)
        for(cc in colCats)
        {
            rclass <- (At==rc) # row
            cclass <- (Dt==cc) # col
            pTable[rc, cc] = -log10(fisher.test(rclass, cclass, alternative = "greater")$p.value);
            CountTbl[rc, cc] = table(rclass  & cclass)["TRUE"]
        }
        CountTbl[is.na(CountTbl)]=0
        # display the p-value and counts in a color-coded table. The colors will indicate the p-value signicance
        # Truncate p values smaller than 10^{-50} to 10^{-50}
        pTable[is.infinite(pTable)] = 1.3*max(pTable[is.finite(pTable)]);
        pTable[pTable>50 ] = 50 ;
        # Marginal counts (really module sizes)
        rTotals = apply(CountTbl, 1, sum)
        cTotals = apply(CountTbl, 2, sum)
        # Use function labeledHeatmap to produce the color-coded table with all the trimmings
        labeledHeatmap( Matrix = pTable, colorLabels = TRUE,
        xLabels = paste("Dt - ", gsub("[.].*","",colCats)), yLabels = paste("At - ", gsub("[.].*","",rowCats)),
        textMatrix = CountTbl, colors = blueWhiteRed(100)[50:100],
        main = paste0(i, ": Correspondence of categories"),
        cex.text = 1, cex.lab = 1, setStdMargins = FALSE      )
    }
}
dev.off()



###################
##### Results #####
###################

x<-read.table("SNPtypes.txt",sep="\t",header=T)
x[1:4,]
# F1 Maxxa TX2094   Freq                      note
# 1    1     1      0     98  allele diagnostic: Maxxa
# 2  1,0     1      0  22065  allele diagnostic: Maxxa
# 3    1     0      1    305 allele diagnostic: tx2094
# 4  1,0     0      1  68549 allele diagnostic: tx2094
sum(x$Freq[1:4] ) # 91017
## Between Maxxa and TX2094 alleles, a total of 91,017 SNPs was detected from 27,816 gene models.

## What are genes diagostic for Maxxa and TX2094 SNPs?
# set working dir to "cistrans_Ranalysis"
load("cistrans.rdata")
# A10, B10, B.mt10, B.tm10, A20, B20, B.mt20, B.tm20, res10, res.mt10, res.tm10, res20, res.mt20, res.tm20, sumT
dim(A10) # 27816
table(gsub("..G.*|Z.*","",gsub("Gohir.","",rownames(A10))))
# 1        A        D  mt_atp6  mt_atp9 mt_ccmFN  mt_cox1  mt_matR
# 285    13462    14060        1        1        1        1        1
# mt_mttB  mt_nad3  mt_nad4 mt_rps12
# 1        1        1        1


## What is the expression divergence between wild and domesticate, among all 66,610 gene?
fn.de10<-"Maxxa.10dpavsTX2094.10dpa.txt"
fn.de20<-"Maxxa.20dpavsTX2094.20dpa.txt"
res<-read.table(paste0("DE-66610/",fn.de10),sep="\t",header=TRUE)
countADZ(rownames(res)) #   A     D    mt     Z 32295 33341    33   941
sig<- res[res$padj<0.05 & !is.na(res$padj),]
dim(sig)  # 6952
id.de10<-rownames(sig)
countADZ(id.de10)  # A    D   mt    Z, 3321 3552    4   75
res<-read.table(paste0("DE-66610/",fn.de20),sep="\t",header=TRUE)
sig<- res[res$padj<0.05 & !is.na(res$padj),]
dim(sig)  # 4878
id.de20<-rownames(sig)
countADZ(id.de20)  # A    D   mt    2368 2450   18   42
de<-unique(c(id.de10,id.de20))
countADZ(de)  # 4834 5083   18   99
n=length(de)  #10034
length(setdiff(de, rownames(A10))) # DE, no SNP 3248
x = length(intersect(de, rownames(A10))) # DE, SNP 6786
p = nrow(A10)/66610 # SNPs in all
## Is the percentage of genes containing SNPs in DE genes higher than that in all genes?
binom.test(x, n, p, alternative = c("greater"), conf.level = 0.95)
## Are these DE genes detected with MxT SNPs? Any association between DE and genetic divergence?
# make 2x2 contingency table for de x snp
M <- as.table(rbind(c(x, n-x), c(nrow(A10)-x, 66610-n-(nrow(A10)-x))))
dimnames(M) <- list(expression = c("DE", "no DE"),
sequence = c("SNPs","no SNPs"))
(Xsq <- chisq.test(M))  # Prints test summary
Xsq$observed   # observed counts (same as M)
Xsq$expected   # expected counts under the null
Xsq$residuals  # Pearson residuals
Xsq$stdres     # standardized residuals


## What is the expression divergence between wild and domesticate, among 27,816 genes?
fn.de10<-"Maxxa.10dpavsTX2094.10dpa.txt"
fn.de20<-"Maxxa.20dpavsTX2094.20dpa.txt"
res<-read.table(paste0("DE-27816/",fn.de10),sep="\t",header=TRUE)
countADZ(rownames(res))  # A     D    mt     Z 13462 14060     9   285
sig<- res[res$padj<0.05 & !is.na(res$padj),]
dim(sig)  # 5570
id.de10<-rownames(sig)
countADZ(id.de10)  # 2677 2830    3   60
res<-read.table(paste0("DE-27816/",fn.de20),sep="\t",header=TRUE)
sig<- res[res$padj<0.05 & !is.na(res$padj),]
dim(sig)  # 3578
id.de20<-rownames(sig)
countADZ(id.de20)  # 1726 1816    8   28
countADZ(de)  # 4834 5083   18   99
de<-unique(c(id.de10,id.de20))
countADZ(de)  # 3753 3951    8   77
length(de)  #7789, noting this number is bigger that 6786 from DE&SNP for 66610



