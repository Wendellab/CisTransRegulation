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

# Differential expression At vs Dt in TX2094 or Maxxa at 10 and 20 dpa
countT = read.table( "Genes66610.raw.count.txt", sep="\t", header=TRUE)
coldataT =read.table( "Genes66610.sample.info.txt", sep="\t", header=TRUE)
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
# make a "DE" folder
system("mkdir DE")
# run DE
apply(batch,1,function(x) pairwiseDE(dds,x,savePath = ""))
# print out results in the order of batch
fileL<- list.files("DE")
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
#         sample 1        sample 2 DE (q<0.05)  1>2  2>1
#1  Maxxa.10dpa.At  Maxxa.10dpa.Dt        7445 3775 3670
#2  Maxxa.20dpa.At  Maxxa.20dpa.Dt        4543 2297 2246
#3 TX2094.10dpa.At TX2094.10dpa.Dt        6819 3414 3405
#4 TX2094.20dpa.At TX2094.20dpa.Dt        4180 2095 2085


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
plotCorrespondence=function(c1,c2,file=NULL,title="", textplot=TRUE,corrplot=TRUE,fisherplot=TRUE,fisher="greater", mai=c(1.02, 0.82,0.82,0.42))
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
        cex.text = 1, cex.lab = 1, setStdMargins = FALSE      )
        par(mai=mai.default)
        
    }
    if(!is.null(file)){dev.off()}
}

# how does domestication changes homoeolog expression bias - ratio
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
#  R1<1  2277  1335    58
#  R1=1  1056 12906   987
#  R1>1    72  1334  2369
xtabs(~R1sig+R2sig,data=cr20)
#   R2sig
# R1sig   R2<1  R2=1  R2>1
# R1<1  1240   980    26
# R1=1   824 16170   857
# R1>1    21  1064  1212
# what percentage is concordant?
sum(diag(xtabs(~R1sig+R2sig,data=cr10)))/nrow(cr10)
# 0.7837814
sum(diag(xtabs(~R1sig+R2sig,data=cr20)))/nrow(cr20)
# 0.831562
# check results from direction comparisons
table(cr10$R1vsR2)/nrow(cr10)
# R1<R2      R1=R2      R1>R2
# 0.03049924 0.93766187 0.03183889
table(cr20$R1vsR2)/nrow(cr20)
# R1<R2      R1=R2      R1>R2
# 0.01330714 0.97021524 0.01647763
crUnion = data.frame(R1sig=(cr10$R1sig!="R1=1" |cr20$R1sig!="R1=1" ), R2sig=(cr10$R2sig!="R2=1" |cr20$R2sig!="R2=1"  ), R1vsR2=(cr10$R1vsR2!="R1=R2" |cr20$R1vsR2!="R1=R2"  ))
rownames(crUnion)=rownames(cr10)
apply(crUnion,2,table)
#  R1sig R2sig R1vsR2
# FALSE 13739 14238  20694
# TRUE   8655  8156   1700
cr10_22394 =cr10
cr20_22394 =cr20
crUnion_22394=crUnion

# focus on 8036 pairs with M vs T snps
cr10 = cr10[as.character(ogP_with_allele_snps$Gohir.A),]  # restrict to 8036 pairs
cr20 = cr20[as.character(ogP_with_allele_snps$Gohir.A),]  # restrict to 8036 pairs
crUnion = crUnion[as.character(ogP_with_allele_snps$Gohir.A),]  # restrict to 8036 pairs
# check
xtabs(~R1sig+R2sig,data=cr10)
#      R2sig
# R1sig  R2<1 R2=1 R2>1
# R1<1 1109  579   34
# R1=1  562 3538  500
# R1>1   44  565 1105
xtabs(~R1sig+R2sig,data=cr20)
# R2sig
# R1sig  R2<1 R2=1 R2>1
# R1<1  605  514   11
# R1=1  320 5127  327
# R1>1   14  565  553
table(cr10$R1vsR2)/nrow(cr10)
#  R1<R2      R1=R2      R1>R2
# 0.05002489 0.90206570 0.04790941
table(cr20$R1vsR2)/nrow(cr20)
# R1<R2      R1=R2      R1>R2
# 0.02103036 0.95246391 0.02650572


save(ogP, ogP_with_allele_snps, cr10, cr20, crUnion, cr10_22394, cr20_22394, crUnion_22394, file="homoeoPairs.rdata")

#######################################################
## homoeolog expression bias x regulatory divergence ##
#######################################################

## Is homoeolog bias associated with RD?
load("RDgenes.rdata")
rd_At = rd[as.character(ogP_with_allele_snps$Gohir.A),]
rd_Dt = rd[as.character(ogP_with_allele_snps$Gohir.D),]
rd_P = rd_At[,c("crd","trd","rd")]
rd_P$crd = rd_At$crd | rd_Dt$crd # cis variants in at least one homoeolog
rd_P$trd = rd_At$trd | rd_Dt$trd # trans variants in at least one homoeolog
rd_P$rd = rd_P$crd | rd_P$trd # cis/trans variants in at least one homoeolog
### are cis/trans variants distributed randomly in homoeolog pairs regardless of bias?
pdf("homoeo.RD&bias.union.pdf")
## Union and rd
plotCorrespondence(rd_P$rd,crUnion$R1sig, title="RD, Maxxa  AtvsDt",textplot=FALSE,corrplot=FALSE,fisherplot=TRUE) # cis
plotCorrespondence(rd_P$rd,crUnion$R2sig, title="RD, TX2094 AtvsDt",textplot=FALSE,corrplot=FALSE,fisherplot=TRUE)
plotCorrespondence(rd_P$rd,crUnion$R1vsR2, title="RD, Maxxa vs TX2094 AtvsDt",textplot=FALSE,corrplot=FALSE,fisherplot=TRUE) # cis
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

## is any RD pattern biased towards At or Dt
pdf("homoeo.AtvsDt.pdf")
# 8036 pairs
plotCorrespondence(rd_At$rd,rd_Dt$rd, title="AtvsDt, 8036 pairs",textplot=FALSE,corrplot=TRUE,fisherplot=TRUE)
# 952 pairs with RD
plotCorrespondence(rd_At$rd[rd_P$rd==TRUE],rd_Dt$rd[rd_P$rd==TRUE], title="AtvsDt, 952 pairs with RD",textplot=FALSE,corrplot=TRUE,fisherplot=TRUE)
dev.off()

--book
xtabs(~crd+R1vsR2,data=cr10)
plotCorrespondence(cr10$crd,cr10$R1vsR2, title="10dpa, MxT, AtvsDt",textplot=FALSE,corrplot=TRUE,fisherplot=TRUE)
plotCorrespondence(cr10$trd,cr10$R1vsR2, title="10dpa, MxT, AtvsDt",textplot=FALSE,corrplot=TRUE,fisherplot=TRUE)
dev.off()

plotCorrespondence(cr10$rd,cr10$categoryDou, title="10dpa, MxT, AtvsDt",textplot=FALSE,corrplot=TRUE,fisherplot=TRUE)
plotCorrespondence(cr10$Adir,cr10$categoryDou, title="10dpa, MxT, AtvsDt",textplot=FALSE,corrplot=TRUE,fisherplot=TRUE)
plotCorrespondence(cr10$rd,paste(cr10$categoryDou, cr10$classDou), title="10dpa, MxT, AtvsDt",textplot=FALSE,corrplot=TRUE,fisherplot=TRUE)
plotCorrespondence(cr10$Adir,paste(cr10$categoryDou, cr10$classDou), title="10dpa, MxT, AtvsDt",textplot=FALSE,corrplot=TRUE,fisherplot=TRUE)


dev.off()

plotCorrespondence(cr10$At.mt,cr10$Dt.mt, title="10dpa, MxT, AtvsDt")
plotCorrespondence(cr10$At.tm,cr10$Dt.tm, title="10dpa, TxM, AtvsDt")
dev.off()
plotCorrespondence(cr10$At.mt[which(cr10$categoryDou=="conserved")],cr10$Dt.mt[which(cr10$categoryDou=="conserved")], title="10dpa, MxT, AtvsDt, conserved",textplot=FALSE,corrplot=TRUE,fisherplot=TRUE)
plotCorrespondence(cr10$At.tm[which(cr10$categoryDou=="conserved")],cr10$Dt.tm[which(cr10$categoryDou=="conserved")], title="10dpa, TxM, AtvsDt, conserved",textplot=FALSE,corrplot=TRUE,fisherplot=TRUE)

plotCorrespondence(cr10$At.mt[which(cr10$categoryDou=="At increase")],cr10$Dt.mt[which(cr10$categoryDou=="At increase")], title="10dpa, MxT, AtvsDt, At increase",textplot=FALSE,corrplot=TRUE,fisherplot=TRUE)
plotCorrespondence(cr10$At.tm[which(cr10$categoryDou=="At increase")],cr10$Dt.tm[which(cr10$categoryDou=="At increase")], title="10dpa, TxM, AtvsDt, At increase",textplot=FALSE,corrplot=TRUE,fisherplot=TRUE)

plotCorrespondence(cr10$At.mt[which(cr10$categoryDou=="Dt increase")],cr10$Dt.mt[which(cr10$categoryDou=="Dt increase")], title="10dpa, MxT, AtvsDt, Dt increase",textplot=FALSE,corrplot=TRUE,fisherplot=TRUE)
plotCorrespondence(cr10$At.tm[which(cr10$categoryDou=="Dt increase")],cr10$Dt.tm[which(cr10$categoryDou=="Dt increase")], title="10dpa, TxM, AtvsDt, Dt increase",textplot=FALSE,corrplot=TRUE,fisherplot=TRUE)
dev.off()


plotFisherCorr(cr10$At.mt,cr10$Dt.mt, title="10dpa, MxT, AtvsDt", fisher="two.sided")
plotFisherCorr(cr10$At.tm,cr10$Dt.tm, title="10dpa, TxM, AtvsDt", fisher="two.sided")
dev.off()




cr20 = cr20[as.character(ogP_with_allele_snps$Gohir.A),]


library(gplots)
pdf("DE/checkDE.pdf")
textplot(T)
mtext("DE analysis result summary")
dev.off()




## correlate parental expression divergence ##
corRes <-c("sample","measure","cor","pvalue")
for(i in c("res.mt10","res.tm10","res.mt20","res.tm20"))
{
    res=get(i)
    for(d in c("A","B","AminusB")){
        r = cor.test(res[ogP_with_allele_snps$Gohir.A,d], res[ogP_with_allele_snps$Gohir.D,d],use='pairwise.complete.obs')
        corRes=rbind(corRes,c(i,d,r$estimate, r$p.value))
    }
}
corRes

---book
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



