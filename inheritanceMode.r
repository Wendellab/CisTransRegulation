#options(repos="https://CRAN.R-project.org")
library(corrplot)
library(WGCNA)
options(scipen=999)

########################
## DE analysis of Mid ##
########################
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
countAll<-read.table("Genes27816.raw.count.txt", sep="\t",header=TRUE)
coldataAll <-read.table("Genes27816.sample.info.txt", sep="\t",header=TRUE)

# index of parental counts
m10 = which(coldataAll$condition2=="Maxxa.10dpa")
t10 = which(coldataAll$condition2=="TX2094.10dpa")
m20 = which(coldataAll$condition2=="Maxxa.20dpa")
t20 = which(coldataAll$condition2=="TX2094.20dpa")

# get mid
getMid=function(count, coldata, index.a,index.b){
    meanLibSize = mean(coldata$lib_size[c(index.a,index.b)])
    g = expand.grid(index.a,index.b)
    i=1
    mid=(count[,g$Var1[i]]/coldata$lib_size[g$Var1[i]] + count[,g$Var2[i]]/coldata$lib_size[g$Var2[i]])*meanLibSize/2
    for(i in 2:nrow(g)){
        mid=cbind(mid,(count[,g$Var1[i]]/coldata$lib_size[g$Var1[i]] + count[,g$Var2[i]]/coldata$lib_size[g$Var2[i]])*meanLibSize/2 )
    }
    mid=data.frame(mid)
    names(mid)=paste0("mid.",paste0("S",coldata$rep[g$Var1],"-",coldata$rep[g$Var2]))
    return(mid)
}
mid10 = getMid(countAll,coldataAll,m10,t10)
mid20 = getMid(countAll,coldataAll,m20,t20)

# F1 index
mt10 = which(coldataAll$condition2=="MxT.10dpa")
tm10 = which(coldataAll$condition2=="TxM.10dpa")
mt20 = which(coldataAll$condition2=="MxT.20dpa")
tm20 = which(coldataAll$condition2=="TxM.20dpa")

# MxT.10dpa vs mid10
count = cbind(countAll[,mt10],mid10);head(count)
info = data.frame(sample=names(count), genome = gsub("[.].*","",names(count)))
dds <- DESeqDataSetFromMatrix( countData = round(count,0), colData =info, design = ~ genome)
res = results(DESeq(dds),contrast=c("genome","F1","mid"))
print( summary(res,alpha=.05) ) # higher F1 4643, Mid 4381
write.table(res, file="DE-27816-mar2019/MxT.10dpavsMid.10dpa.txt",  sep="\t")

# TxM.10dpa vs mid10
count = cbind(countAll[,tm10],mid10);head(count)
info = data.frame(sample=names(count), genome = gsub("[.].*","",names(count)))
dds <- DESeqDataSetFromMatrix( countData = round(count,0), colData =info, design = ~ genome)
res = results(DESeq(dds),contrast=c("genome","F1","mid"))
print( summary(res,alpha=.05) ) # higher F1 2839, Mid 2497
write.table(res, file="DE-27816-mar2019/TxM.10dpavsMid.10dpa.txt",  sep="\t")

# MxT.20dpa vs mid20
count = cbind(countAll[,mt20],mid20);head(count)
info = data.frame(sample=names(count), genome = gsub("[.].*","",names(count)))
dds <- DESeqDataSetFromMatrix( countData = round(count,0), colData =info, design = ~ genome)
res = results(DESeq(dds),contrast=c("genome","F1","mid"))
print( summary(res,alpha=.05) ) # higher F1 1946, Mid 1713
write.table(res, file="DE-27816-mar2019/MxT.20dpavsMid.20dpa.txt",  sep="\t")

# TxM.20dpa vs mid20
count = cbind(countAll[,tm20],mid20);head(count)
info = data.frame(sample=names(count), genome = gsub("[.].*","",names(count)))
dds <- DESeqDataSetFromMatrix( countData = round(count,0), colData =info, design = ~ genome)
res = results(DESeq(dds),contrast=c("genome","F1","mid"))
print( summary(res,alpha=.05) ) # higher F1 5239, Mid 4955
write.table(res, file="DE-27816-mar2019/TxM.20dpavsMid.20dpa.txt",  sep="\t")


###################################
## Expression Dominance Analysis ##
###################################

# Parents - P1, P2 or A, D
# Mid parental - Mid
# Hybrid/allopolyploid - T

classDominance<-function(TvsMid, TvsP1, TvsP2, P1vsP2, log2fc.threshold=0, reverseTvsP =FALSE, Pnames=c("A","D"))
{
    # Hybrid/polyploid vs Mid parental val
    TvsMid <- data.frame(TvsMid[,c("log2FoldChange", "lfcSE", "padj")])
    names(TvsMid) <- c("TvsMid", "TvsMid.SE", "TvsMid.padj")
    # Hybrid/polyploid vs Parent 1
    TvsP1 <- data.frame(TvsP1[,c("log2FoldChange", "lfcSE", "padj")])
    names(TvsP1) <- c("TvsP1", "TvsP1.SE", "TvsP1.padj")
    # Hybrid/polyploid vs Parent 2
    TvsP2 <- data.frame(TvsP2[,c("log2FoldChange", "lfcSE", "padj")])
    names(TvsP2) <- c("TvsP2", "TvsP2.SE", "TvsP2.padj")
    # Parent 1 vs Parent 2
    P1vsP2 <- data.frame(P1vsP2[,c("log2FoldChange", "lfcSE", "padj")])
    names(P1vsP2) <- c("P1vsP2", "P1vsP2.SE", "P1vsP2.padj")

    if(reverseTvsP==TRUE){
        TvsP1$TvsP1 = -TvsP1$TvsP1
        TvsP2$TvsP2 = -TvsP2$TvsP2
    }
    
    tbl = cbind(TvsMid, TvsP1, TvsP2, P1vsP2)
    
    tbl$TvsMid[is.na(tbl$TvsMid)] =0
    tbl$TvsP1[is.na(tbl$TvsP1)] =0
    tbl$TvsP2[is.na(tbl$TvsP2)] =0
    tbl$P1vsP2[is.na(tbl$P1vsP2)] =0
    
    # judge
    tbl$additivity <- ifelse(abs(tbl$TvsMid)>log2fc.threshold & tbl$TvsMid.padj<0.05 & !is.na(tbl$TvsMid.padj), "T!=Mid", "T=Mid")
    tbl$TvsP1.reg <- ifelse(abs(tbl$TvsP1)>log2fc.threshold & tbl$TvsP1.padj<0.05 & !is.na(tbl$TvsP1.padj), "T!=P1", "T=P1")
    tbl$TvsP2.reg <- ifelse(abs(tbl$TvsP2)>log2fc.threshold & tbl$TvsP2.padj<0.05 & !is.na(tbl$TvsP2.padj), "T!=P2", "T=P2")
    tbl$P1vsP2.reg <- ifelse(abs(tbl$P1vsP2)>log2fc.threshold & tbl$P1vsP2.padj<0.05 & !is.na(tbl$P1vsP2.padj), ifelse(tbl$P1vsP2>log2fc.threshold & tbl$P1vsP2.padj<0.05 & !is.na(tbl$P1vsP2.padj), "P1>P2","P1<P2"), "P1=P2")
    
    # together
    tbl$class <- paste(tbl$P1vsP2.reg, tbl$additivity, tbl$TvsP1.reg, tbl$TvsP2.reg, sep=";")
    
    # assign category
    tbl$category = "7.Other"
    tbl$category[grep("T=Mid",tbl$class)] = "1.Additivity" # hybrid=Mid while P1!=P2
    tbl$category[grep("P1=P2;T=Mid",tbl$class)] = "6.Conserved"  # P1=P2=Mid=hybrid
    tbl$category[grep("P1>P2;T!=Mid;T=P1;T!=P2",tbl$class)] = paste0("2.",Pnames[1],"-dominant, higher")
    tbl$category[grep("P1<P2;T!=Mid;T=P1;T!=P2",tbl$class)] = paste0("2.",Pnames[1],"-dominant, lower")
    tbl$category[grep("P1>P2;T!=Mid;T!=P1;T=P2",tbl$class)] = paste0("3.",Pnames[2],"-dominant, lower")
    tbl$category[grep("P1<P2;T!=Mid;T!=P1;T=P2",tbl$class)] = paste0("3.",Pnames[2],"-dominant, higher")
    tbl$category[grepl("T!=Mid;T!=P1;T!=P2",tbl$class) & tbl$TvsP1>0 & tbl$TvsP2>0 & tbl$P1vsP2>0] = paste0("4.Transgressive Up: ",Pnames[1]," higher")
    tbl$category[grepl("T!=Mid;T!=P1;T!=P2",tbl$class) & tbl$TvsP1>0 & tbl$TvsP2>0 & tbl$P1vsP2<0] = paste0("4.Transgressive Up: ",Pnames[2]," higher")
    tbl$category[grepl("T!=Mid;T!=P1;T!=P2",tbl$class) & tbl$TvsP1<0 & tbl$TvsP2<0 & tbl$P1vsP2>0] = paste0("5.Transgressive Down: ",Pnames[1]," higher")
    tbl$category[grepl("T!=Mid;T!=P1;T!=P2",tbl$class) & tbl$TvsP1<0 & tbl$TvsP2<0 & tbl$P1vsP2<0] = paste0("5.Transgressive Down: ",Pnames[2]," higher")
    
    return(tbl)
}

### Get Parents vs F1
MvsMT10<-read.table("DE-27816-mar2019/Maxxa.10dpavsMxT.10dpa.txt", header=TRUE, sep="\t")
TvsMT10<-read.table("DE-27816-mar2019/TX2094.10dpavsMxT.10dpa.txt", header=TRUE, sep="\t")
#
MvsTM10<-read.table("DE-27816-mar2019/Maxxa.10dpavsTxM.10dpa.txt", header=TRUE, sep="\t")
TvsTM10<-read.table("DE-27816-mar2019/TX2094.10dpavsTxM.10dpa.txt", header=TRUE, sep="\t")
#
MvsMT20<-read.table("DE-27816-mar2019/Maxxa.20dpavsMxT.20dpa.txt", header=TRUE, sep="\t")
TvsMT20<-read.table("DE-27816-mar2019/TX2094.20dpavsMxT.20dpa.txt", header=TRUE, sep="\t")
#
MvsTM20<-read.table("DE-27816-mar2019/Maxxa.20dpavsTxM.20dpa.txt", header=TRUE, sep="\t")
TvsTM20<-read.table("DE-27816-mar2019/TX2094.20dpavsTxM.20dpa.txt", header=TRUE, sep="\t")
# get F1 vs Mid
MT10vsMid<-read.table("DE-27816-mar2019/MxT.10dpavsMid.10dpa.txt", header=TRUE, sep="\t")
TM10vsMid<-read.table("DE-27816-mar2019/TxM.10dpavsMid.10dpa.txt", header=TRUE, sep="\t")
MT20vsMid<-read.table("DE-27816-mar2019/MxT.20dpavsMid.20dpa.txt", header=TRUE, sep="\t")
TM20vsMid<-read.table("DE-27816-mar2019/TxM.20dpavsMid.20dpa.txt", header=TRUE, sep="\t")
# get P1 vs P2
MvsT10<-read.table("DE-27816-mar2019/Maxxa.10dpavsTX2094.10dpa.txt", header=TRUE, sep="\t")
MvsT20<-read.table("DE-27816-mar2019/Maxxa.20dpavsTX2094.20dpa.txt", header=TRUE, sep="\t")

# categorization of inheritance mode
im.mt10 = classDominance(MT10vsMid, MvsMT10, TvsMT10, MvsT10, log2fc.threshold=0, reverseTvsP=TRUE,Pnames=c("Maxxa","TX2094"))
im.tm10 = classDominance(TM10vsMid, MvsTM10, TvsTM10, MvsT10, log2fc.threshold=0, reverseTvsP=TRUE,Pnames=c("Maxxa","TX2094"))
im.mt20 = classDominance(MT20vsMid, MvsMT20, TvsMT20, MvsT20, log2fc.threshold=0, reverseTvsP=TRUE,Pnames=c("Maxxa","TX2094"))
im.tm20 = classDominance(TM20vsMid, MvsTM20, TvsTM20, MvsT20, log2fc.threshold=0, reverseTvsP=TRUE,Pnames=c("Maxxa","TX2094"))

save(im.mt10, im.tm10, im.mt20, im.tm20, file="inheritance.rdata")

# between reciprocal F1s, not as consistent as regulatory categories
table(im.mt10$category==im.tm10$category) # 8985 T-18831
table(im.mt20$category==im.tm20$category) # 8861 T-18955

# summarize category only for A!=0
class=sort(unique(matrix(apply(im[,-1],2,as.character))))
sumT=data.frame(mt10=as.numeric(table(im.mt10$category[im.mt10$P1vsP2.reg!="P1=P2"])[class]), tm10=as.numeric(table(im.tm10$category[im.tm10$P1vsP2.reg!="P1=P2"])[class]), 
mt20=as.numeric(table(im.mt20$category[im.mt20$P1vsP2.reg!="P1=P2"])[class]), 
tm20=as.numeric(table(im.tm20$category[im.tm20$P1vsP2.reg!="P1=P2"])[class]))
rownames(sumT)=class
sumT[is.na(sumT)]=0
sumT$all =rowSums(sumT)
print(sumT)
# percentage
sumP<-round(sweep(sumT,2,colSums(sumT),"/")*100,1)
names(sumP) = paste0(names(sumP),"%")
write.table(cbind(sumT,sumP), file ="inheritance.summary.txt",sep="\t")
sumP[1,]
colSums(sumP[2:5,])
colSums(sumP[6:9,])


############################################
## Relate cis trans with inheritance mode ##
############################################
load("cistrans.rdata")

### Does the composition of inheritance mode change with increased parental expression divergence?
REs<-c(res.mt10$category,res.tm10$category,res.mt20$category,res.tm20$category)
IMs<-c(im.mt10$category,im.tm10$category,im.mt20$category,im.tm20$category)
As<-c(im.mt10$P1vsP2,im.tm10$P1vsP2,im.mt20$P1vsP2,im.tm20$P1vsP2)
reg<-c(im.mt10$P1vsP2.reg,im.tm10$P1vsP2.reg,im.mt20$P1vsP2.reg,im.tm20$P1vsP2.reg)
select = which(reg!="P1=P2")
breaks <- c(0, 1, 2, 3, 4, 100)
c<-.bincode(x=abs(As[select]), b=breaks, TRUE)
cName<-c("(0-1)","(1-2)","(2-3)","(3-4)","(4+)")
res= as.matrix(table(IMs[select],c))
pdf("inheritance.with|A|.pdf")
barplot(res)
barplot(round(sweep(res,2,colSums(res),"/")*100,1) )
dev.off()
# smaller portion of dominance??? not clear

# FUN
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


pdf("inheritance.combined.pdf")
## all combined
inheritance = c(im.mt10$category, im.tm10$category,im.mt20$category, im.tm20$category)
regulation= c(res.mt10$category, res.tm10$category,res.mt20$category, res.tm20$category)
# overall plot
plotCorrespondence(inheritance,regulation,title="All",mai=c(1.72,2.32,0.82,0.42))
# restrict to parental divergence
select= grepl("^[1-5]",inheritance) & grepl("^[1-4]",regulation)
plotCorrespondence(inheritance[select],regulation[select],textplot=FALSE,title="All",mai=c(1.72,2.32,0.82,0.42))
# large category
im1= inheritance
im1[grep("dominant",im1)]="2.Dominance"
im1[grep("Transgressive",im1)]="3.Transgression"
plotFisherCorr(im1[select],regulation[select])
# finer
im2= inheritance
im2[grep("Maxxa-dominant",im2)]="2.Maxxa-dominant"
im2[grep("TX2094-dominant",im2)]="3.TX2094-dominant"
im2[grep("Transgressive Up",im2)]="4.Transgressive Up"
im2[grep("Transgressive Down",im2)]="4.Transgressive Down"
plotFisherCorr(im2[select],regulation[select])
dev.off()


pdf("inheritance.each.pdf")
for(dataset in c("mt10","tm10","mt20","tm20"))
{
    inheritance=get(paste0("im.",dataset))$category
    regulation=get(paste0("res.",dataset))$category
    plotCorrespondence(inheritance,regulation,title=dataset,mai=c(1.72,2.32,0.82,0.42))
    select= grepl("^[1-5]",inheritance) & grepl("^[1-4]",regulation)
    plotCorrespondence(inheritance[select],regulation[select],textplot=FALSE,title=dataset,mai=c(1.72,2.32,0.82,0.42))
    im2= inheritance
    im2[grep("dominant",im2)]="2.Dominance"
    im2[grep("Transgressive",im2)]="3.Transgression"
    plotFisherCorr(im2[select],regulation[select])
}
dev.off()




inheritance = c(im.mt10$category, im.tm10$category,im.mt20$category, im.tm20$category)
regulation= c(res.mt10$category, res.tm10$category,res.mt20$category, res.tm20$category)
sample=c(rep("mt10",27816),rep("tm10",27816),rep("mt20",27816),rep("tm20",27816))
im2= inheritance
im2[grep("Maxxa-dominant",im2)]="2.Maxxa-dominant"
im2[grep("TX2094-dominant",im2)]="3.TX2094-dominant"
im2[grep("Transgressive Up",im2)]="4.Transgressive Up"
im2[grep("Transgressive Down",im2)]="4.Transgressive Down"
df=as.data.frame(ftable(xtabs(~im2+regulation+sample)))
df$dpa=gsub("[mt]","",df$sample)
df$f1=gsub("[0-2]","",df$sample)
write.table(df, file ="inheritance.sub.txt",sep="\t")



# M-dominant vs T-dominant across samples, for each regulationcategory, points
df2=df[grepl("^2",df$im)&grepl("^[1-4]",df$regulation),]
df3=df[grepl("^3",df$im)&grepl("^[1-4]",df$regulation),]
df2$perc = df2$Freq/(df2$Freq+df3$Freq)*100
p1=ggplot(data=df2, aes(x=regulation, y=perc, color=dpa, shape=f1)) +
geom_point(size=4) +geom_hline(yintercept=50) + theme_minimal() +ylab("Dominant/recessive")
# Up or down transgression across samples, for each regulationcategory, points
df4=df[grepl("Transgressive Up",df$im)&grepl("^[1-4]",df$regulation),]
df5=df[grepl("Transgressive Down",df$im)&grepl("^[1-4]",df$regulation),]
df4$perc = df4$Freq/(df4$Freq+df5$Freq)*100
p2=ggplot(data=df4, aes(x=regulation, y=perc, color=dpa, shape=f1)) +
geom_point(size=4) +geom_hline(yintercept=50) + theme_minimal() +ylab("Transgressive up/down")


pdf("inheritance.sub.pdf")
df2$test="Dominance"
df4$test="Transgression"
dft=rbind(df2,df4)
ggplot(data=dft, aes(x=regulation, y=perc, color=dpa, shape=f1)) +
geom_point(size=4) +geom_hline(yintercept=50) + theme_minimal() + facet_grid(test~.)+theme(axis.text=element_text(size=12))
dev.off()


