## Compare regulatory pattern between homoeologs
library(gplots)
library(ggplot2)
options(scipen=999)
load("cistrans.rdata")

###########################
## Import ortholog group ##
###########################

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
ogP_with_allele_snps<-ogP[ogP$Gohir.A %in% ids & ogP$Gohir.D %in% ids, ]
paste0(nrow(ogP_with_allele_snps), " homoeolog pair OGs were detected with SNPs between Maxxa and TX2094.")
# "8036 homoeolog pair OGs were detected with SNPs between Maxxa and TX2094."
names(ogP)[1]<-"OG"

##### in comparison with previous OG results from Justin
# v1 - A total of 23700 out of 29636 OGs contain At/Dt homoeolog pairs.
# 8266 homoeolog pair OGs were detected with SNPs between Maxxa and TX2094.
# ------
# v2 - A total of 22048 out of 28731 OGs contain At/Dt homoeolog pairs.
# 7896 homoeolog pair OGs were detected with SNPs between Maxxa and TX2094.
# ------
# v3 - A total of 22394 out of 28843 OGs contain At/Dt homoeolog pairs.
# "8036 homoeolog pair OGs were detected with SNPs between Maxxa and TX2094."

table(ids[grep("Gohir.A",ids)] %in% ogP_with_allele_snps$Gohir.A)
# FALSE  TRUE
# 5426  8036
table(ids[grep("Gohir.D",ids)] %in% ogP_with_allele_snps$Gohir.D)
# FALSE  TRUE
# 6024  8036
countADZ(ids[!(ids %in% ogP_with_allele_snps$Gohir.D | ids %in% ogP_with_allele_snps$Gohir.A)])
#    A    D   mt    Z
# 5426 6024    9  285

save(ogP, ogP_with_allele_snps, file="homoeoPairs.rdata")

###############################
## Homoeolog expression bias ##
###############################
load("homoeoPairs.rdata")
pairwiseDE<-function(dds, contrast,savePath)
{
    # DE analysis
    print(contrast)
    ddsPW <-dds[,dds$condition %in% contrast]
    ddsPW$condition<-droplevels(ddsPW$condition)
    res <- results(DESeq(ddsPW), contrast = c("condition",contrast))
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

# Differential expression of At+Dt in TX2094 or Maxxa at 10 and 20 dpa
countT = read.table( "Genes66610.raw.count.txt", sep="\t", header=TRUE)
coldataT =read.table( "Genes66610.sample.info.txt", sep="\t", header=TRUE)
tt = countT[as.character(ogP$Gohir.A),]+countT[as.character(ogP$Gohir.D),]
library(DESeq2)
dds <- DESeqDataSetFromMatrix( countData = tt, colData = coldataT, design = ~ condition)
# pairwise deseq workflow
batch<- rbind(
c("Maxxa.10dpa", "TX2094.10dpa" ),
c("Maxxa.20dpa", "TX2094.20dpa" ))
# make a "DE" folder
system("mkdir DE")
apply(batch,1,function(x) pairwiseDE(dds,x,savePath = ""))
# print out results in the order of batch
fileL<-paste0("DE/",batch[,1],"vs",batch[,2],".txt")
sigT<-c("sample 1","sample 2","DE (q<0.05)","1>2","2>1")
for(file in fileL)
{
    res<-read.table(file,sep="\t",header=TRUE)
    sigRes <- c(unlist(strsplit(gsub(".txt","",file),split="vs") ), getSig(res),getSig(res,direction="up"),getSig(res,direction="down"))
    sigT<-rbind(sigT,sigRes)
}
T<-as.data.frame(sigT[-1,], row.names=FALSE)
names(T)<-sigT[1,]
print(T)
#     sample 1     sample 2 DE (q<0.05)  1>2  2>1
# 1 Maxxa.10dpa TX2094.10dpa        3035 1610 1425 13.4% out 224394
# 2 Maxxa.20dpa TX2094.20dpa        2137 1200  937 9.5%
## out of 22394 pairs of homoeologs, 3035/2137 show expression changes in aggregated expression by domestication

# Differential expression At vs Dt in TX2094 or Maxxa at 10 and 20 dpa
a = countT[as.character(ogP$Gohir.A),]; names(a) = paste0(names(a),".At")
d = countT[as.character(ogP$Gohir.D),]; names(d) = paste0(names(d),".Dt")
countH = cbind(a,d)
coldataH = rbind(coldataT, coldataT)
rownames(coldataH)= names(countH)
coldataH$condition =gsub("dpa.*[.]","dpa.",rownames(coldataH))
# differential expression to get expression divergence
library(DESeq2)
dds <- DESeqDataSetFromMatrix( countData = countH, colData = coldataH, design = ~ condition)
# pairwise deseq workflow
batch<- rbind(
c("Maxxa.10dpa.At", "Maxxa.10dpa.Dt" ),
c("Maxxa.20dpa.At", "Maxxa.20dpa.Dt" ),
c("TX2094.10dpa.At", "TX2094.10dpa.Dt" ),
c("TX2094.20dpa.At", "TX2094.20dpa.Dt" ))
# run DE
apply(batch,1,function(x) pairwiseDE(dds,x,savePath = ""))
# print out results in the order of batch
fileL<-paste0("DE/", batch[,1],"vs",batch[,2],".txt")
sigT<-c("sample 1","sample 2","DE (q<0.05)","1>2","2>1")
for(file in fileL)
{
    res<-read.table(file,sep="\t",header=TRUE)
    sigRes <- c(unlist(strsplit(gsub(".txt","",file),split="vs") ), getSig(res),getSig(res,direction="up"),getSig(res,direction="down"))
    sigT<-rbind(sigT,sigRes)
}
T<-as.data.frame(sigT[-1,], row.names=FALSE)
names(T)<-sigT[1,]
print(T)
#         sample 1        sample 2 DE (q<0.05)  1>2  2>1
#1  Maxxa.10dpa.At  Maxxa.10dpa.Dt        7445 3670 3775
#2  Maxxa.20dpa.At  Maxxa.20dpa.Dt        4543 2246 2297
#3 TX2094.10dpa.At TX2094.10dpa.Dt        6819 3405 3414
#4 TX2094.20dpa.At TX2094.20dpa.Dt        4180 2085 2095
## out of 22394 pairs of homoeologs, homoeolog bias in either accession/dpa

#####################################################
## Domestication effects on homoeolog contribution ##
#####################################################

# FUN
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
###plot correspondecne
library(corrplot)
library(WGCNA)
options(scipen=999)
plotFisherCorr=function(c1,c2, title="", fisher="two.sided")
{
    tbl=table(data.frame(c1,c2))
    CountTbl = as.matrix(tbl)
    # residue
    chisq <- chisq.test(tbl,simulate.p.value = TRUE)
    residualTable = chisq$residuals
    ## table plot with WGCNA function
    pTable = CountTbl
    rowCats = levels(factor(c1))
    colCats = levels(factor(c2))
    for(rc in rowCats)
    for(cc in colCats)
    {
        rclass <- (c1==rc) # row
        cclass <- (c2==cc) # col
        pTable[rc, cc] = fisher.test(rclass, cclass, alternative = fisher)$p.value;
        CountTbl[rc, cc] = table(rclass  & cclass)["TRUE"]
    }
    CountTbl[is.na(CountTbl)]=0
    # Marginal counts (really module sizes)
    rP = round(rowSums(tbl)/sum(rowSums(tbl))*100,1)
    cP = round(colSums(tbl)/sum(colSums(tbl))*100,1)
    rownames(residualTable) = paste0(rownames(residualTable)," ",rP,"%")
    colnames(residualTable) = paste0(colnames(residualTable)," ",cP,"%")
    
    corrplot(residualTable, p.mat=as.matrix(pTable), insig="label_sig", sig.level=c(0.001,0.01,0.05),is.corr=FALSE, method="color",pch.col="white",tl.col="black")
}
plotCorrespondence=function(c1,c2,file=NULL,title="", textplot=TRUE,corrplot=TRUE,fisherplot=TRUE,fisher="greater", mai=c(1.02, 0.82,0.82,0.42), cex=1)
{
    if(!is.null(file)){pdf(file)}
    tbl=table(data.frame(c1,c2))
    if(textplot==TRUE){
        textplot(tbl)
        mtext(title)
    }
    ## corrplot of residule
    CountTbl = as.matrix(tbl)
    chisq <- chisq.test(tbl,simulate.p.value = TRUE)
    residualTable = chisq$residuals
    main=paste0(title, " Chi-square test P = ",chisq$p.value)
    if(corrplot==TRUE){
        corrplot(residualTable, p.mat=CountTbl, insig="p-value", sig.level=0,is.corr=FALSE, main=main, method="color")
    }
    if(fisherplot==TRUE)
    {
        ## table plot with WGCNA function
        pTable = CountTbl
        rowCats = levels(factor(c1))
        colCats = levels(factor(c2))
        for(rc in rowCats)
        for(cc in colCats)
        {
            rclass <- (c1==rc) # row
            cclass <- (c2==cc) # col
            pTable[rc, cc] = -log10(fisher.test(rclass, cclass, alternative = fisher)$p.value);
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
        mai.default=par("mai")
        par(mai=mai)
        labeledHeatmap( Matrix = pTable, colorLabels = TRUE,
        xLabels = colCats, yLabels = rowCats,
        textMatrix = CountTbl, colors = blueWhiteRed(100)[50:100],
        main = paste0(title," correspondence of categories"),
        cex.text = cex, cex.lab = cex, setStdMargins = FALSE      )
        par(mai=mai.default)
        
    }
    if(!is.null(file)){dev.off()}
}

# how does domestication changes homoeolog expression bias - ratio
### get total At+Dt changes
Tt10<-read.table("DE/Maxxa.10dpavsTX2094.10dpa.txt", header=TRUE, sep="\t")
Tt20<-read.table("DE/Maxxa.20dpavsTX2094.20dpa.txt", header=TRUE, sep="\t")
length(which( (Tt10$padj<0.05&!is.na(Tt10$padj)) | (Tt20$padj<0.05&!is.na(Tt20$padj))  ))
# 4482 show total expression change by domestication
Tt10$dom = "M=T";
Tt10$dom[Tt10$padj<0.05&!is.na(Tt10$padj)&Tt10$log2FoldChange>0]="M>T"
Tt10$dom[Tt10$padj<0.05&!is.na(Tt10$padj)&Tt10$log2FoldChange<0]="M<T"
Tt20$dom = "M=T";
Tt20$dom[Tt20$padj<0.05&!is.na(Tt20$padj)&Tt20$log2FoldChange>0]="M>T"
Tt20$dom[Tt20$padj<0.05&!is.na(Tt20$padj)&Tt20$log2FoldChange<0]="M<T"
# check numbers
table(Tt10$dom)
# M<T   M=T   M>T
# 1425 19359  1610
table(Tt20$dom)
# M<T   M=T   M>T
# 937 20257  1200
table(a=Tt10$dom,b=Tt20$dom)
#     b
# a     M<T   M=T   M>T
# M<T   255  1145    25
# M=T   667 17912   780
# M>T    15  1200   395
table(Tt10$dom!="M=T"|Tt20$dom!="M=T")
# FALSE  TRUE
# 17912  4482

### Get allelic expression divergence, test for B=0
Rm10<-read.table("DE/Maxxa.10dpa.AtvsMaxxa.10dpa.Dt.txt", header=TRUE, sep="\t")
Rm20<-read.table("DE/Maxxa.20dpa.AtvsMaxxa.20dpa.Dt.txt", header=TRUE, sep="\t")
Rt10<-read.table("DE/TX2094.10dpa.AtvsTX2094.10dpa.Dt.txt", header=TRUE, sep="\t")
Rt20<-read.table("DE/TX2094.20dpa.AtvsTX2094.20dpa.Dt.txt", header=TRUE, sep="\t")
length(which( (Rm10$padj<0.05&!is.na(Rm10$padj)) | (Rm20$padj<0.05&!is.na(Rm20$padj))  ))
# 8655 show bias in maxxa fibers
length(which( (Rt10$padj<0.05&!is.na(Rt20$padj)) | (Rt20$padj<0.05&!is.na(Rt20$padj))  ))
# 8042/22394 = 0.359114 bias in tx2094 fiber

# compare A/D between Maxxa and TX2094
compareRatios<-function(r1.res, r2.res, r1.n, r2.n,log2fc.threshold=0)
{
    r1 <- r1.res[,c("log2FoldChange", "lfcSE", "padj")]
    names(r1) <- c("R1", "R1.SE", "R1.padj")
    r2 <- r2.res[,c("log2FoldChange", "lfcSE", "padj")]
    names(r2) <- c("R2", "R2.SE", "R2.padj")
    r1<-r1[rownames(r2),]
    table <- cbind(r1,r2)
    table$R1minusR2 <- table$R1 - table$R2
    table$R1minusR2.pvalue <- apply(table,1,function(x) t.test2(m1=x[1],m2=x[4],se1=x[2], se2=x[5], n1=r1.n, n2=r2.n)["p-value"])
    # cat
    table$R1sig <- ifelse( table$R1.padj<0.05 & !is.na(table$R1.padj), ifelse(table$R1>log2fc.threshold, "R1>1","R1<1"), "R1=1")
    table$R2sig <- ifelse( table$R2.padj<0.05 & !is.na(table$R2.padj), ifelse(table$R2>log2fc.threshold, "R2>1","R2<1"), "R2=1")
    table$R1vsR2 <- ifelse( table$R1minusR2.pvalue<0.05 & !is.na(table$R1minusR2.pvalue), ifelse(table$R1minusR2>log2fc.threshold, "R1>R2","R1<R2"), "R1=R2")
    # two ways to R1 vs R1: classTri based on student t test, classDou based on cross tabulation
    table$classTri <- paste(table$R1sig,table$R2sig,table$R1vsR2,sep=";")
    table$classDou <- paste(table$R1sig,table$R2sig,sep=";")
    table$categoryDou <- "conserved"
    table$categoryDou[table$classDou == "R1>1;R2<1" | table$classDou == "R1>1;R2=1" | table$classDou == "R1=1;R2<1" ] <- "R1>R2"
    table$categoryDou[table$classDou == "R1<1;R2>1" | table$classDou == "R1<1;R2=1" | table$classDou == "R1=1;R2>1" ] <- "R1<R2"
    return(table)
}
# run
cr10=compareRatios(Rm10, Rt10, 4,3)
cr20=compareRatios(Rm20, Rt20, 3,3)
# check results from cross tabulation
xtabs(~R1sig+R2sig,data=cr10)
#   R2sig
# R1sig   R2<1  R2=1  R2>1
# R1<1  2369  1334    72
# R1=1   987 12906  1056
# R1>1    58  1335  2277
xtabs(~R1sig+R2sig,data=cr20)
#   R2sig
# R1sig   R2<1  R2=1  R2>1
# R1<1  1212  1064    21
# R1=1   857 16170   824
# R1>1    26   980  1240
# what percentage is concordant?
sum(diag(xtabs(~R1sig+R2sig,data=cr10)))/nrow(cr10)
# 0.7837814
sum(diag(xtabs(~R1sig+R2sig,data=cr20)))/nrow(cr20)
# 0.831562
# check results from direction comparisons
table(cr10$R1vsR2)/nrow(cr10)
#      R1<R2      R1=R2      R1>R2
# 0.03183889 0.93766187 0.03049924
table(cr20$R1vsR2)/nrow(cr20)
#      R1<R2      R1=R2      R1>R2
# 0.01647763 0.97021524 0.01330714
crUnion = data.frame(R1sig=(cr10$R1sig!="R1=1" |cr20$R1sig!="R1=1" ), R2sig=(cr10$R2sig!="R2=1" |cr20$R2sig!="R2=1"  ), R1vsR2=(cr10$R1vsR2!="R1=R2" |cr20$R1vsR2!="R1=R2"  ), Ttsig = Tt10$dom!="M=T"|Tt20$dom!="M=T")
rownames(crUnion)=rownames(cr10)
apply(crUnion,2,table)
#       R1sig R2sig R1vsR2 Ttsig
# FALSE 13739 14238  20694 17912
# TRUE   8655  8156   1700  4482
#       38.6% 36.4%   7.6%  20.0%
cr10_22394 =cr10
cr20_22394 =cr20
crUnion_22394=crUnion
Tt10_22394 = Tt10
Tt20_22394 = Tt20

# focus on 8036 pairs with M vs T snps
cr10 = cr10[as.character(ogP_with_allele_snps$Gohir.A),]  # restrict to 8036 pairs
cr20 = cr20[as.character(ogP_with_allele_snps$Gohir.A),]  # restrict to 8036 pairs
crUnion = crUnion[as.character(ogP_with_allele_snps$Gohir.A),]  # restrict to 8036 pairs
Tt10 = Tt10[as.character(ogP_with_allele_snps$Gohir.A),]  # restrict to 8036 pairs
Tt20 = Tt20[as.character(ogP_with_allele_snps$Gohir.A),]  # restrict to 8036 pairs
N =nrow(cr10) #8036
apply(crUnion,2,table)
#       R1sig R2sig R1vsR2 Ttsig
# FALSE  3985  4190   7085  5779
# TRUE   4051  3846    951  2257
# TRUE%  0.5041065 0.4785963 0.1183425 0.2808611
# check
xtabs(~R1sig+R2sig,data=cr10)
#      R2sig
# R1sig  R2<1 R2=1 R2>1
# R1<1 1105  565   44
# R1=1  500 3538  562
# R1>1   34  579 1109
xtabs(~R1sig+R2sig,data=cr20)
# R2sig
# R1sig  R2<1 R2=1 R2>1
#   R1<1  553  565   14
#   R1=1  327 5127  320
#   R1>1   11  514  605
table(cr10$R1vsR2)/nrow(cr10)
#  R1<R2      R1=R2      R1>R2
# 0.04790941 0.90206570 0.05002489
table(cr20$R1vsR2)/nrow(cr20)
# R1<R2      R1=R2      R1>R2
# 0.02650572 0.95246391 0.02103036

save(ogP, ogP_with_allele_snps, cr10, cr20, Tt10, Tt20, crUnion, cr10_22394, cr20_22394, crUnion_22394, Tt10_22394, Tt20_22394, file="homoeoPairs.rdata")

cr=cbind(cr10,cr20)
write.table(cr,file="homoeo.cr10n20.txt",sep="\t")



load("RDgenes.rdata")
rd$b10=ifelse(rd$mt10==rd$tm10, as.character(rd$mt10),"NA")
rd$b20=ifelse(rd$mt20==rd$tm20, as.character(rd$mt20),"NA")
rd_At = rd[as.character(ogP_with_allele_snps$Gohir.A),]
rd_Dt = rd[as.character(ogP_with_allele_snps$Gohir.D),]
rd_P = rd_At[,c("crd","trd","rd")]
rd_P$crd = rd_At$crd | rd_Dt$crd # cis variants in at least one homoeolog
rd_P$trd = rd_At$trd | rd_Dt$trd # trans variants in at least one homoeolog
rd_P$rd = rd_P$crd | rd_P$trd # cis/trans variants in at least one homoeolog

##############################################
## total expression x regulatory divergence ##
##############################################
## How does RD pattern change total or aggregated expression of homoeologs? This question was not dealt with elsewhere in the manuscript, because homoeologs are treated as individual genes.

# 
ftable(table(rd_At$b10,rd_Dt$b10,Tt10$dom))
M<T  M=T  M>T

cis and trans cis and trans     0   12    7
cis only          0    6    5
no divergence     7  218   29
trans only        0    4    4
cis only      cis and trans     2    1    5
cis only          6    4    4
no divergence    15   71   26
trans only        2    1    1
no divergence cis and trans    10  198   39
cis only         14   79   30
no divergence   265 6244  575
trans only        9   24   31
trans only    cis and trans     0    5    5
cis only          1    1    1
no divergence     8   30   23
trans only        4    8    2

ftable(table(rd_P$rd,Tt20$dom))
M<T  M=T  M>T

FALSE   265 6244  575
TRUE     78  662  212



#######################################################
## homoeolog expression bias x regulatory divergence ##
#######################################################
## Is homoeolog bias associated with RD?
### are cis/trans variants distributed randomly in homoeolog pairs regardless of bias?
pdf("homoeo.RD&bias.union.pdf")
## Union and rd
plotCorrespondence(rd_P$rd,crUnion$R1sig, title="RD, Maxxa  AtvsDt",textplot=FALSE,corrplot=FALSE,fisherplot=TRUE) # cis
plotCorrespondence(rd_P$rd,crUnion$R2sig, title="RD, TX2094 AtvsDt",textplot=FALSE,corrplot=FALSE,fisherplot=TRUE)
plotCorrespondence(rd_P$rd,crUnion$R1vsR2, title="RD, Maxxa vs TX2094 AtvsDt",textplot=FALSE,corrplot=FALSE,fisherplot=TRUE) # cis
plotCorrespondence(rd_P$rd,crUnion$Ttsig, title="RD, Maxxa vs TX2094 At+Dt",textplot=FALSE,corrplot=FALSE,fisherplot=TRUE) # cis
dev.off()
pdf("homoeo.RD&bias.pdf")
## 10dpa
plotCorrespondence(rd_P$crd,cr10$R1sig, title="10dpa, cis, Maxxa  AtvsDt",textplot=FALSE,corrplot=TRUE,fisherplot=FALSE) # cis
plotCorrespondence(rd_P$crd,cr10$R2sig, title="10dpa, cis, TX2094 AtvsDt",textplot=FALSE,corrplot=TRUE,fisherplot=FALSE)
plotCorrespondence(rd_P$trd,cr10$R1sig, title="10dpa, trans, Maxxa AtvsDt",textplot=FALSE,corrplot=TRUE,fisherplot=FALSE) # trans
plotCorrespondence(rd_P$trd,cr10$R2sig, title="10dpa, trans, TX2094 AtvsDt",textplot=FALSE,corrplot=TRUE,fisherplot=FALSE)
plotCorrespondence(rd_P$rd,cr10$R1sig, title="10dpa, cis/trans, Maxxa AtvsDt",textplot=FALSE,corrplot=TRUE,fisherplot=FALSE) # cis or trans
plotCorrespondence(rd_P$rd,cr10$R2sig, title="10dpa, cis/trans, TX2094 AtvsDt",textplot=FALSE,corrplot=TRUE,fisherplot=FALSE)
# what about altered homoeolog ratio
plotCorrespondence(rd_P$rd,cr10$R1vsR2, title="10dpa, cis/trans, TX2094 vs Maxxa AtvsDt",textplot=FALSE,corrplot=TRUE,fisherplot=TRUE)
## 20 dpa
plotCorrespondence(rd_P$crd,cr20$R1sig, title="20dpa, cis, Maxxa  AtvsDt",textplot=FALSE,corrplot=TRUE,fisherplot=FALSE) # cis
plotCorrespondence(rd_P$crd,cr20$R2sig, title="20dpa, cis, TX2094 AtvsDt",textplot=FALSE,corrplot=TRUE,fisherplot=FALSE)
plotCorrespondence(rd_P$trd,cr20$R1sig, title="20dpa, trans, Maxxa AtvsDt",textplot=FALSE,corrplot=TRUE,fisherplot=FALSE) # trans
plotCorrespondence(rd_P$trd,cr20$R2sig, title="20dpa, trans, TX2094 AtvsDt",textplot=FALSE,corrplot=TRUE,fisherplot=FALSE)
plotCorrespondence(rd_P$rd,cr20$R1sig, title="20dpa, cis/trans, Maxxa AtvsDt",textplot=FALSE,corrplot=TRUE,fisherplot=FALSE) # cis or trans
plotCorrespondence(rd_P$rd,cr20$R2sig, title="20dpa, cis/trans, TX2094 AtvsDt",textplot=FALSE,corrplot=TRUE,fisherplot=FALSE)
# what about altered homoeolog ratio
plotCorrespondence(rd_P$rd,cr20$R1vsR2, title="20dpa, cis/trans, TX2094 vs Maxxa AtvsDt",textplot=FALSE,corrplot=TRUE,fisherplot=TRUE)
dev.off()

## Given association between bias and RD, how cis-trans change the range of existing biases? 
with(cr10[cr10$R2sig!="R2=1" & rd_P$rd==TRUE,],t.test(R1,R2)) # not sig
with(cr20[cr20$R2sig!="R2=1" & rd_P$rd==TRUE,],t.test(R1,R2)) # not sig
# SO no bias in up or descrease ratio, R1 vs R2
## then how were the magnitude of existing bias change by cis-trans: bias in TX2094, contain RD, showed sig changes between ratio
with(cr10,t.test(abs(R1),abs(R2))) # not sig
with(cr10[cr10$R2sig=="R2=1" & rd_P$rd==TRUE,],t.test(abs(R1),abs(R2))) 
# sig R1>R2, increased magnitude without pre-existing bias
with(cr10[cr10$R2sig!="R2=1" & rd_P$rd==TRUE,],t.test(abs(R1),abs(R2))) 
# sig R1<R2, decreased magnitube with pre-existing bias
with(cr20,t.test(abs(R1),abs(R2))) # sig R1>R2
with(cr20[cr20$R2sig=="R2=1" & rd_P$rd==TRUE,],t.test(abs(R1),abs(R2))) # sig R1>R2 0.9349404 0.4715869
with(cr20[cr20$R2sig!="R2=1" & rd_P$rd==TRUE,],t.test(abs(R1),abs(R2))) # sig R1<R2 1.801027  2.204624
# SO reduced magnitude of ratios by RD

## is any RD pattern biased towards At or Dt ----- NO
rd_At$rd=factor(rd_At$rd, levels=c("cis only","trans only","cis and trans","no divergence"))
rd_Dt$rd=factor(rd_Dt$rd, levels=c("cis only","trans only","cis and trans","no divergence"))
pdf("homoeo.AtvsDt.pdf")
# 8036 pairs
plotCorrespondence(rd_At$rd,rd_Dt$rd, title="AtvsDt, 8036 pairs",textplot=FALSE,corrplot=TRUE,fisherplot=TRUE,cex=1.5)
# 952 pairs with RD
plotCorrespondence(rd_At$rd[rd_P$rd==TRUE],rd_Dt$rd[rd_P$rd==TRUE], title="AtvsDt, 952 pairs with RD",textplot=FALSE,corrplot=TRUE,fisherplot=TRUE, cex=1.5)
dev.off()

## For 952 pairs, how each combination of At and Dt RD change ratio and total homoeolog expression?
select = which(rd_P$rd==TRUE)
AtDt = paste0("At - ",rd_At$rd,", Dt - ", rd_Dt$rd )
pdf("homoeo.16AtDt.pdf")
plotCorrespondence(AtDt[select],cr10$R1vsR2[select], title="16 At-Dt RD groups, 10 dpa ratio",textplot=FALSE,corrplot=TRUE,fisherplot=TRUE, cex=1.5)
plotCorrespondence(AtDt[select],cr20$R1vsR2[select], title="16 At-Dt RD groups, 20 dpa ratio",textplot=FALSE,corrplot=TRUE,fisherplot=TRUE, cex=1.5)
plotCorrespondence(AtDt[select], Tt10$dom[select], title="16 At-Dt RD groups, 10 dpa total",textplot=FALSE,corrplot=TRUE,fisherplot=TRUE, cex=1.5)
plotCorrespondence(AtDt[select], Tt20$dom[select], title="16 At-Dt RD groups, 20 dpa total",textplot=FALSE,corrplot=TRUE,fisherplot=TRUE, cex=1.5)
dev.off()


ftable(table(rd_At$rd,rd_Dt$rd, cr10$R1vsR2))
ftable(table(rd_At$rd,rd_Dt$rd, cr10$R1vsR2))

## how does RD pattern change homoeolog contribution ratio?
load("cistrans.rdata")
pdf("homoe.ratioChangebyRD.pdf")
# does certain regulatory pattern of At change its ration? what about total expression?
boxplot(cr10[as.character(ogP_with_allele_snps$Gohir.A),"R1minusR2"]~res.mt10[as.character(ogP_with_allele_snps$Gohir.A),"category"], main = "At contribution by mt10 cis/trans category")
boxplot(cr10[as.character(ogP_with_allele_snps$Gohir.A),"R1minusR2"]~res.tm10[as.character(ogP_with_allele_snps$Gohir.A),"category"], main = "At contribution by tm10 cis/trans category")
boxplot(cr10[as.character(ogP_with_allele_snps$Gohir.A),"R1minusR2"]~res.mt20[as.character(ogP_with_allele_snps$Gohir.A),"category"], main = "At contribution by mt20 cis/trans category")
boxplot(cr10[as.character(ogP_with_allele_snps$Gohir.A),"R1minusR2"]~res.tm20[as.character(ogP_with_allele_snps$Gohir.A),"category"], main = "At contribution by tm20 cis/trans category")
boxplot(cr10[as.character(ogP_with_allele_snps$Gohir.A),"R1minusR2"]~rd_At$rd, main = "At contribution by RD type")
# does certain regulatory pattern of Dt change its ration? what about total expression?
boxplot(cr10[as.character(ogP_with_allele_snps$Gohir.A),"R1minusR2"]~res.mt10[as.character(ogP_with_allele_snps$Gohir.D),"category"], main = "Dt contribution by mt10 cis/trans category")
boxplot(cr10[as.character(ogP_with_allele_snps$Gohir.A),"R1minusR2"]~res.tm10[as.character(ogP_with_allele_snps$Gohir.D),"category"], main = "Dt contribution by tm10 cis/trans category")
boxplot(cr10[as.character(ogP_with_allele_snps$Gohir.A),"R1minusR2"]~res.mt20[as.character(ogP_with_allele_snps$Gohir.D),"category"], main = "Dt contribution by mt20 cis/trans category")
boxplot(cr10[as.character(ogP_with_allele_snps$Gohir.A),"R1minusR2"]~res.tm20[as.character(ogP_with_allele_snps$Gohir.D),"category"], main = "Dt contribution by tm20 cis/trans category")
boxplot(cr10[as.character(ogP_with_allele_snps$Gohir.A),"R1minusR2"]~rd_Dt$rd, main = "Dt contribution by RD type")
dev.off()
# lump together
pdf("homoe.ratioChangebyRD.union.pdf")
u=data.frame(ratio = c(cr10[as.character(ogP_with_allele_snps$Gohir.A),"R1minusR2"], cr20[as.character(ogP_with_allele_snps$Gohir.A),"R1minusR2"]), category =  c(cr10[as.character(ogP_with_allele_snps$Gohir.A),"R1vsR2"], cr20[as.character(ogP_with_allele_snps$Gohir.A),"R1vsR2"]), At =c(as.character(rd_At$rd),as.character(rd_At$rd)), Dt=c(as.character(rd_Dt$rd),as.character(rd_Dt$rd)) )
u$At=factor(u$At, levels=c("cis only","trans only","cis and trans","no divergence"))
u$Dt=factor(u$Dt, levels=c("cis only","trans only","cis and trans","no divergence"))
plotCorrespondence(u$category,u$At, title="At: RD vs homoeolog ratio change",textplot=FALSE,corrplot=TRUE,fisherplot=TRUE,cex=1)
plotCorrespondence(u$category,u$Dt, title="Dt: RD vs homoeolog ratio change",textplot=FALSE,corrplot=TRUE,fisherplot=TRUE,cex=1)
boxplot(ratio~category,data=u, main="double check")
boxplot(ratio~At,data=u, main = "At/Dt change by RD of At")
boxplot(ratio~Dt,data=u, main = "At/Dt change by RD of Dt")
# cis only of Dt preferentially decrease Dt contribution, but not change At much
# tran only of At decrease At contribution, trans only of Dt also decrease Dt contribution
# prettier boxplot
library(reshape2)
u2<-melt(u,measure.vars=c("At","Dt"))
u2$value=factor(u2$value, levels=c("cis only","trans only","cis and trans","no divergence"))
print(p<-ggplot(u2, aes(x=value,y=ratio,color=variable))+geom_boxplot()+theme_bw())
dev.off()
# whether R1vsR2 deviate from 0
t.test(u$ratio[u$At=="cis only"]) # p-value = 0.6856
t.test(u$ratio[u$Dt=="cis only"]) # p-value = 0.002579, **
t.test(u$ratio[u$At=="trans only"]) # p-value = 0.0008632 ***
t.test(u$ratio[u$Dt=="trans only"]) # p-value = 0.1442
t.test(u$ratio[u$At=="cis and trans"]) # p-value = 0.4261
t.test(u$ratio[u$Dt=="cis and trans"]) # p-value = 0.002401, **
t.test(u$ratio[u$At=="no divergence"]) # p-value = 0.1601
t.test(u$ratio[u$Dt=="no divergence"]) # p-value = 0.9187
# compare
t.test(u$ratio[u$At=="cis only"],u$ratio[u$Dt=="cis only"]) # p-value = 0.06717
t.test(u$ratio[u$At=="trans only"],u$ratio[u$Dt=="trans only"]) # p-value = 0.0006546
t.test(u$ratio[u$At=="cis and trans"],u$ratio[u$Dt=="cis and trans"]) # p-value = 0.07942
t.test(u$ratio[u$At=="no divergence"],u$ratio[u$Dt=="no divergence"]) # p-value = 0.353



## compare seq divergence between homoelogs
load("seqDivergence.rdata")
load("paml/paml.rdata") #paml
S <- seqDev[,-1]
rownames(S)<-seqDev$ID
S <- S[rownames(res.mt10),]
rownames(paml)=paml$id
R=paml[rownames(res.mt10),]
S$dNdS=R$"dN/dS"
S$dS=R$"dS"
S$dN=R$"dN"
# homoeo
Sa = S[as.character(ogP_with_allele_snps$Gohir.A),]
Sd = S[as.character(ogP_with_allele_snps$Gohir.D),]
unique(rownames(rd_At)==rownames(Sa))
Sa$rd=rd_At$rd
Sd$rd=rd_Dt$rd
Sa$genome="At"
Sd$genome="Dt"
SP=rbind(Sa,Sd)

# correlation
pdf("homoe.seqDev.pdf")
for(col in names(S)){
    print(col)
    res=cor.test(Sa[,col],Sd[,col])
    plot(Sa[,col],Sd[,col], col=alpha("black", 0.2), pch=16, main=paste0(col, ": r = ",res$estimate,", p = ",round(res$p.value,4)))
    print(wilcox.test(Sa[,col],Sd[,col])) # significant difference in s, t.indel, t.indel.size
    print(t.test(Sa[,col],Sd[,col]))
}
print(p<-ggplot(SP, aes(x=rd,y=s,color=genome))+geom_boxplot()+theme_bw())
print(p<-ggplot(SP, aes(x=rd,y=s.indel,color=genome))+geom_boxplot()+theme_bw())
print(p<-ggplot(SP, aes(x=rd,y=exon,color=genome))+geom_boxplot()+theme_bw())
print(p<-ggplot(SP, aes(x=rd,y=dS,color=genome))+geom_boxplot()+theme_bw())
print(p<-ggplot(SP, aes(x=rd,y=dN,color=genome))+geom_boxplot()+theme_bw())
print(p<-ggplot(SP, aes(x=rd,y=dN/dS,color=genome))+geom_boxplot()+theme_bw())

dev.off()
