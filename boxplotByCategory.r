## Categoryzation of cis tran
library(gplots)
library(ggplot2)
load("cistrans.rdata")

plotBox=function(x, category,ylab="", title="")
{
     regCol <- c("red","blue","purple","brown","green","black","grey")
     boxplot(x~category, main = title, col=regCol, names=NULL, ylab =ylab, xlab ="cis/trans regulatory category")
}

pdf("boxplot.|A|.pdf")
par(mfrow=c(2,3))
plotBox(abs(res.mt10$A), res.mt10$category, ylab="|A|",title="MxT 10dpa")
plotBox(abs(res.tm10$A), res.tm10$category, ylab="|A|",title="TxM 10dpa")
plotBox(abs(res.mt20$A), res.mt20$category, ylab="|A|",title="MxT 20dpa")
plotBox(abs(res.tm20$A), res.tm20$category, ylab="|A|",title="TxM 20dpa")
plotBox(abs(res10$A), res10$category, ylab="|A|",title="F1 10dpa")
plotBox(abs(res20$A), res20$category, ylab="abs(log2A)",title="F2 10dpa")
dev.off()
wilcox.test(abs(res.mt10$A[grep("^1",res.mt10$category)]),abs(res.mt10$A[grep("^2",res.mt10$category)]))

pdf("boxplot.A.pdf")
par(mfrow=c(2,3))
plotBox(res.mt10$A, res.mt10$category, ylab="A",title="MxT 10dpa")
plotBox(res.tm10$A, res.tm10$category, ylab="A",title="TxM 10dpa")
plotBox(res.mt20$A, res.mt20$category, ylab="A",title="MxT 20dpa")
plotBox(res.tm20$A, res.tm20$category, ylab="A",title="TxM 20dpa")
plotBox(res10$A, res10$category, ylab="A",title="F1 10dpa")
plotBox(res20$A, res20$category, ylab="A",title="F2 10dpa")
dev.off()
t.test(res.mt10$A[grep("^1",res.mt10$category)]) # <0 *
t.test(res.mt10$A[grep("^2",res.mt10$category)]) # >0 *
t.test(res.mt10$A[grep("^3",res.mt10$category)])
t.test(res.mt10$A[grep("^4",res.mt10$category)])
t.test(res.mt10$A[grep("^5",res.mt10$category)])
t.test(res.mt10$A[grep("^6",res.mt10$category)])
t.test(res.mt10$A[grep("^7",res.mt10$category)])

pdf("boxplot.|B|.pdf")
par(mfrow=c(2,3))
plotBox(abs(res.mt10$B), res.mt10$category, ylab="|B|",title="MxT 10dpa")
plotBox(abs(res.tm10$B), res.tm10$category, ylab="|B|",title="TxM 10dpa")
plotBox(abs(res.mt20$B), res.mt20$category, ylab="|B|",title="MxT 20dpa")
plotBox(abs(res.tm20$B), res.tm20$category, ylab="|B|",title="TxM 20dpa")
plotBox(abs(res10$B), res10$category, ylab="|B|",title="F1 10dpa")
plotBox(abs(res20$B), res20$category, ylab="|B|",title="F2 10dpa")
dev.off()
wilcox.test(abs(res.mt10$B[grep("^1",res.mt10$category)]),abs(res.mt10$B[grep("^2",res.mt10$category)]))


pdf("boxplot.B.pdf")
par(mfrow=c(2,3))
plotBox(res.mt10$B, res.mt10$category, ylab="B",title="MxT 10dpa")
plotBox(res.tm10$B, res.tm10$category, ylab="B",title="TxM 10dpa")
plotBox(res.mt20$B, res.mt20$category, ylab="B",title="MxT 20dpa")
plotBox(res.tm20$B, res.tm20$category, ylab="B",title="TxM 20dpa")
plotBox(res10$B, res10$category, ylab="B",title="F1 10dpa")
plotBox(res20$B, res20$category, ylab="B",title="F2 10dpa")
dev.off()

pdf("boxplot.Fig2C.pdf")
par(mfrow=c(2,2))
plotBox(abs(res.mt10$A), res.mt10$category, ylab="|A|",title="MxT 10dpa")
plotBox(res.mt10$A, res.mt10$category, ylab="A",title="MxT 10dpa")
plotBox(abs(res.mt10$B), res.mt10$category, ylab="|B|",title="MxT 10dpa")
plotBox(res.mt10$B, res.mt10$category, ylab="B",title="MxT 10dpa")
dev.off()

load("RDgenes.rdata")
load("seqDivergence.rdata")
## restrict to 27816 genes
S <- seqDev[,-1]
rownames(S)<-seqDev$ID
S <- S[rownames(res.mt10),]
pdf("boxplot.seqDivergence.pdf")
par(mfrow=c(2,2))
plotBox(S$transcript, res.mt10$category, ylab="snps/site",title="MxT 10dpa")
plotBox(S$p2k.s, res.mt10$category, ylab="snp/site",title="MxT 10dpa")
plotBox(S$p2k.s.indel, res.mt10$category, ylab="indel/site",title="MxT 10dpa")
plotBox(S$p2k.s.indel.size, res.mt10$category, ylab="indel_size/site",title="MxT 10dpa")
boxplot(S$transcript~rd$rd, ylab="snps/site",main="RD")
boxplot(S$p2k.s~rd$rd, ylab="snp/site",main="RD")
boxplot(S$p2k.s.indel~rd$rd, ylab="indel/site",main="RD")
boxplot(S$p2k.s.indel.size~rd$rd, ylab="indel_size/site",main="Rd")
dev.off()

load("RDgenes.rdata")
load("seqDivergence.rdata")
## restrict to 27816 genes
S <- seqDev[,-1]
rownames(S)<-seqDev$ID
S <- S[rownames(res.mt10),]
pdf("boxplot.seqDivergence.pdf")
par(mfrow=c(2,2))
plotBox(S$transcript, res.mt10$category, ylab="snps/site",title="MxT 10dpa")
plotBox(S$p2k.s, res.mt10$category, ylab="snp/site",title="MxT 10dpa")
plotBox(S$p2k.s.indel, res.mt10$category, ylab="indel/site",title="MxT 10dpa")
plotBox(S$p2k.s.indel.size, res.mt10$category, ylab="indel_size/site",title="MxT 10dpa")
boxplot(S$transcript~rd$rd, ylab="snps/site",main="RD")
boxplot(S$p2k.s~rd$rd, ylab="snp/site",main="RD")
boxplot(S$p2k.s.indel~rd$rd, ylab="indel/site",main="RD")
boxplot(S$p2k.s.indel.size~rd$rd, ylab="indel_size/site",main="Rd")
dev.off()


load("coexpressionNet/buildNetwork.rdata")
Kt <-  k.tx2094[rownames(res.mt10),]
Km <-  k.maxxa[rownames(res.mt10),]
pdf("boxplot.netConnectivity.pdf")
par(mfrow=c(2,2))
# TX2094
plotBox(Kt$kTotal, res.mt10$category, ylab="TX2094 kTotal",title="MxT 10dpa")
plotBox(Kt$kTotal, res.tm10$category, ylab="TX2094 kTotal",title="TxM 10dpa")
plotBox(Kt$kTotal, res.mt20$category, ylab="TX2094 kTotal",title="MxT 20dpa")
plotBox(Kt$kTotal, res.tm20$category, ylab="TX2094 kTotal",title="TxM 20dpa")
plotBox(Kt$kWithin, res.mt10$category, ylab="TX2094 kWithin",title="MxT 10dpa")
plotBox(Kt$kWithin, res.tm10$category, ylab="TX2094 kWithin",title="TxM 10dpa")
plotBox(Kt$kWithin, res.mt20$category, ylab="TX2094 kWithin",title="MxT 20dpa")
plotBox(Kt$kWithin, res.tm20$category, ylab="TX2094 kWithin",title="TxM 20dpa")
# MAxxa
plotBox(Km$kTotal, res.mt10$category, ylab="Maxxa kTotal",title="MxT 10dpa")
plotBox(Km$kTotal, res.tm10$category, ylab="Maxxa kTotal",title="TxM 10dpa")
plotBox(Km$kTotal, res.mt20$category, ylab="Maxxa kTotal",title="MxT 20dpa")
plotBox(Km$kTotal, res.tm20$category, ylab="Maxxa kTotal",title="TxM 20dpa")
plotBox(Km$kWithin, res.mt10$category, ylab="Maxxa kWithin",title="MxT 10dpa")
plotBox(Km$kWithin, res.tm10$category, ylab="Maxxa kWithin",title="TxM 10dpa")
plotBox(Km$kWithin, res.mt20$category, ylab="Maxxa kWithin",title="MxT 20dpa")
plotBox(Km$kWithin, res.tm20$category, ylab="Maxxa kWithin",title="TxM 20dpa")
# RD
boxplot(Kt$kTotal~rd$rd, ylab="TX2094 kTotal",title="RD")
boxplot(Kt$kWithin~rd$rd, ylab="TX2094 kWithin",title="RD")
boxplot(Km$kTotal~rd$rd, ylab="Maxxa kTotal",title="RD")
boxplot(Km$kWithin~rd$rd, ylab="Maxxa kWithin",title="RD")
dev.off()

------book

# categories with expression divergence
pdf("plotCistrans.log2A_category.pdf")
# summarize numbers in each category
tbl <- cbind(table(res10$category), table(res20$category), table(res.mt10$category), table(res.mt20$category),table(res.tm10$category), table(res.tm20$category))
colnames(tbl) <- c("F1 10dpa", "F1 20dpa", "MxT 10dpa", "MxT 20dpa","TxM 10dpa", "TxM 20dpa")
# percentages of first four
tbl4<-sweep(tbl[1:4,],2,colSums(tbl[1:4,]),"/")
textplot(tbl,cex=0.6)
mtext("Numbers of genes in each regulatory category")
textplot(round(tbl4*100,1),cex=0.6)
mtext("Percentage of genes in four regulatory category that exhibited parental divergence")
# percentages of first two
colSums(tbl4[1:2,])
par(mfrow=c(3,3))
regCol <- c("red","blue","purple","brown","green","white","grey")
boxplot(abs(res.mt10$A)~res.mt10$category, main = "MxT 10dpa", col=regCol, names=NULL,ylab ="log2 Parental Divergence", xlab ="cis/trans regulatory category")
#         legend(4.5,7.8,legend=levels(factor(res10$category)), col=regCol, pch=20)
boxplot(abs(res.mt20$A)~res.mt20$category, main = "MxT 20dpa", col=regCol,names=NULL,xlab ="cis/trans regulatory category")
boxplot(abs(res.tm10$A)~res.tm10$category, main = "TxM 10dpa", col=regCol,names=NULL, ylab ="log2 Parental Divergence", xlab ="cis/trans regulatory category")
boxplot(abs(res.tm20$A)~res.tm20$category, main = "TxM 20dpa", col=regCol, names=NULL, xlab ="cis/trans regulatory category")
boxplot(abs(res10$A)~res10$category, main = "F1 10dpa", col=regCol, names=NULL, ylab ="log2 Parental Divergence", xlab ="cis/trans regulatory category")
boxplot(abs(res20$A)~res20$category, main = "F1 20dpa", col=regCol, names=NULL, xlab ="cis/trans regulatory category")
dev.off()


pdf("plotCistrans.log2A_category.At.pdf")
# summarize numbers in each category
select=grep("Gohir.A",rownames(res10))
tbl <- cbind(table(res10$category[select]), table(res20$category[select]), table(res.mt10$category[select]), table(res.mt20$category[select]),table(res.tm10$category[select]), table(res.tm20$category[select]))
colnames(tbl) <- c("F1 10dpa", "F1 20dpa", "MxT 10dpa", "MxT 20dpa","TxM 10dpa", "TxM 20dpa")
# percentages of first four
tbl4<-sweep(tbl[1:4,],2,colSums(tbl[1:4,]),"/")
textplot(tbl,cex=0.6)
mtext("Numbers of genes in each regulatory category")
textplot(round(tbl4*100,1),cex=0.6)
mtext("Percentage of genes in four regulatory category that exhibited parental divergence")
# percentages of first two
colSums(tbl4[1:2,])
par(mfrow=c(3,3))
regCol <- c("red","blue","purple","brown","green","white","grey")
boxplot(abs(res.mt10$A)~res.mt10$category, main = "MxT 10dpa", col=regCol, names=NULL,ylab ="log2 Parental Divergence", xlab ="cis/trans regulatory category")
#         legend(4.5,7.8,legend=levels(factor(res10$category)), col=regCol, pch=20)
boxplot(abs(res.mt20$A)~res.mt20$category, main = "MxT 20dpa", col=regCol,names=NULL,xlab ="cis/trans regulatory category")
boxplot(abs(res.tm10$A)~res.tm10$category, main = "TxM 10dpa", col=regCol,names=NULL, ylab ="log2 Parental Divergence", xlab ="cis/trans regulatory category")
boxplot(abs(res.tm20$A)~res.tm20$category, main = "TxM 20dpa", col=regCol, names=NULL, xlab ="cis/trans regulatory category")
boxplot(abs(res10$A)~res10$category, main = "F1 10dpa", col=regCol, names=NULL, ylab ="log2 Parental Divergence", xlab ="cis/trans regulatory category")
boxplot(abs(res20$A)~res20$category, main = "F1 20dpa", col=regCol, names=NULL, xlab ="cis/trans regulatory category")
dev.off()


pdf("plotCistrans.log2A_category.Dt.pdf")
# summarize numbers in each category
select=grep("Gohir.D",rownames(res10))
tbl <- cbind(table(res10$category[select]), table(res20$category[select]), table(res.mt10$category[select]), table(res.mt20$category[select]),table(res.tm10$category[select]), table(res.tm20$category[select]))
colnames(tbl) <- c("F1 10dpa", "F1 20dpa", "MxT 10dpa", "MxT 20dpa","TxM 10dpa", "TxM 20dpa")
# percentages of first four
tbl4<-sweep(tbl[1:4,],2,colSums(tbl[1:4,]),"/")
textplot(tbl,cex=0.6)
mtext("Numbers of genes in each regulatory category")
textplot(round(tbl4*100,1),cex=0.6)
mtext("Percentage of genes in four regulatory category that exhibited parental divergence")
# percentages of first two
colSums(tbl4[1:2,])
par(mfrow=c(3,3))
regCol <- c("red","blue","purple","brown","green","white","grey")
boxplot(abs(res.mt10$A)~res.mt10$category, main = "MxT 10dpa", col=regCol, names=NULL,ylab ="log2 Parental Divergence", xlab ="cis/trans regulatory category")
#         legend(4.5,7.8,legend=levels(factor(res10$category)), col=regCol, pch=20)
boxplot(abs(res.mt20$A)~res.mt20$category, main = "MxT 20dpa", col=regCol,names=NULL,xlab ="cis/trans regulatory category")
boxplot(abs(res.tm10$A)~res.tm10$category, main = "TxM 10dpa", col=regCol,names=NULL, ylab ="log2 Parental Divergence", xlab ="cis/trans regulatory category")
boxplot(abs(res.tm20$A)~res.tm20$category, main = "TxM 20dpa", col=regCol, names=NULL, xlab ="cis/trans regulatory category")
boxplot(abs(res10$A)~res10$category, main = "F1 10dpa", col=regCol, names=NULL, ylab ="log2 Parental Divergence", xlab ="cis/trans regulatory category")
boxplot(abs(res20$A)~res20$category, main = "F1 20dpa", col=regCol, names=NULL, xlab ="cis/trans regulatory category")
dev.off()




# the union of each category
id= rownames(res10)
for(i in levels(factor(res10$category)))
{
    query = i
    gl = unique( id[res10$category == query | res20$category == query | res.mt10$category == query | res.mt20$category == query | res.tm10$category == query | res.tm20$category == query] )
    print(paste(i, length(gl), length(gl)/length(id)))
}
