plotBox=function(x, category,ylab="", title="")
{
     regCol <- c("red","blue","purple","brown","green","black","grey")
     boxplot(x~category, main = title, col=regCol, names=NULL, ylab =ylab, xlab ="cis/trans regulatory category")
}


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

