################################################################
## Magnitude and direction of parental and allelic divergence ##
################################################################
## Boxplot by category I-VII              
library(gplots)
library(ggplot2)
load("cistrans.rdata")

plotBox=function(x, category,ylab="", title="")
{
     regCol <- c("red","blue","purple","brown","green","black","grey")
     boxplot(x~category, main = title, col=regCol, names=NULL, ylab =ylab, xlab ="cis/trans regulatory category",medlwd=1,outpch=20)
}

pdf("boxplot.|A|.pdf")
par(mfrow=c(2,2))
plotBox(abs(res.mt10$A), res.mt10$category, ylab="|A|",title="MxT 10dpa")
plotBox(abs(res.tm10$A), res.tm10$category, ylab="|A|",title="TxM 10dpa")
plotBox(abs(res.mt20$A), res.mt20$category, ylab="|A|",title="MxT 20dpa")
plotBox(abs(res.tm20$A), res.tm20$category, ylab="|A|",title="TxM 20dpa")
#plotBox(abs(res10$A), res10$category, ylab="|A|",title="F1 10dpa")
#plotBox(abs(res20$A), res20$category, ylab="abs(log2A)",title="F2 10dpa")
dev.off()
wilcox.test(abs(res.mt10$A[grep("^1",res.mt10$category)]),abs(res.mt10$A[grep("^2",res.mt10$category)])) # p=0.12
wilcox.test(abs(res.mt10$A[grep("^1",res.mt10$category)]),abs(res.mt10$A[grep("^3",res.mt10$category)])) # p-value = 1.465e-07

pdf("boxplot.A.pdf")
par(mfrow=c(2,2))
plotBox(res.mt10$A, res.mt10$category, ylab="A",title="MxT 10dpa")
plotBox(res.tm10$A, res.tm10$category, ylab="A",title="TxM 10dpa")
plotBox(res.mt20$A, res.mt20$category, ylab="A",title="MxT 20dpa")
plotBox(res.tm20$A, res.tm20$category, ylab="A",title="TxM 20dpa")
#plotBox(res10$A, res10$category, ylab="A",title="F1 10dpa")
#plotBox(res20$A, res20$category, ylab="A",title="F2 10dpa")
dev.off()
t.test(res.mt10$A[grep("^1",res.mt10$category)]) # <0 *
t.test(res.mt10$A[grep("^2",res.mt10$category)]) # >0 *
t.test(res.mt10$A[grep("^3",res.mt10$category)])
t.test(res.mt10$A[grep("^4",res.mt10$category)])
t.test(res.mt10$A[grep("^5",res.mt10$category)])
t.test(res.mt10$A[grep("^6",res.mt10$category)])
t.test(res.mt10$A[grep("^7",res.mt10$category)])

pdf("boxplot.|B|.pdf")
par(mfrow=c(2,2))
plotBox(abs(res.mt10$B), res.mt10$category, ylab="|B|",title="MxT 10dpa")
plotBox(abs(res.tm10$B), res.tm10$category, ylab="|B|",title="TxM 10dpa")
plotBox(abs(res.mt20$B), res.mt20$category, ylab="|B|",title="MxT 20dpa")
plotBox(abs(res.tm20$B), res.tm20$category, ylab="|B|",title="TxM 20dpa")
#plotBox(abs(res10$B), res10$category, ylab="|B|",title="F1 10dpa")
#plotBox(abs(res20$B), res20$category, ylab="|B|",title="F2 10dpa")
dev.off()
wilcox.test(abs(res.mt10$B[grep("^1",res.mt10$category)]),abs(res.mt10$B[grep("^2",res.mt10$category)]))


pdf("boxplot.B.pdf")
par(mfrow=c(2,2))
plotBox(res.mt10$B, res.mt10$category, ylab="B",title="MxT 10dpa")
plotBox(res.tm10$B, res.tm10$category, ylab="B",title="TxM 10dpa")
plotBox(res.mt20$B, res.mt20$category, ylab="B",title="MxT 20dpa")
plotBox(res.tm20$B, res.tm20$category, ylab="B",title="TxM 20dpa")
#plotBox(res10$B, res10$category, ylab="B",title="F1 10dpa")
#plotBox(res20$B, res20$category, ylab="B",title="F2 10dpa")
dev.off()

pdf("boxplot.Fig2.mt10.pdf")
par(mfrow=c(2,2))
plotBox(abs(res.mt10$A), res.mt10$category, ylab="|A|",title="MxT 10dpa")
plotBox(res.mt10$A, res.mt10$category, ylab="A",title="MxT 10dpa")
plotBox(abs(res.mt10$B), res.mt10$category, ylab="|B|",title="MxT 10dpa")
plotBox(res.mt10$B, res.mt10$category, ylab="B",title="MxT 10dpa")
dev.off()


res4 <- rbind(res.mt10, res.tm10, res.mt20, res.tm20)
pdf("boxplot.Fig2.allF1s.pdf")
par(mfrow=c(2,2))
plotBox(abs(res4$A), res4$category, ylab="|A|",title="2 F1 x 2 dpa")
plotBox(res4$A, res4$category, ylab="A",title="2 F1 x 2 dpa")
plotBox(abs(res4$B), res4$category, ylab="|B|",title="2 F1 x 2 dpa")
plotBox(res4$B, res4$category, ylab="B",title="2 F1 x 2 dpa")
dev.off()
# |A|
t.test(abs(res4$A[grep("^3",res4$category)]))
wilcox.test(abs(res4$A[grep("^1",res4$category)]),abs(res4$A[grep("^2",res4$category)])) # p-value = 2.89e-05, cis only > trans only
wilcox.test(abs(res4$A[grep("^1",res4$category)]),abs(res4$A[grep("^3",res4$category)])) # p-value = 1.604e-13, enhanced highest
# A
t.test(res4$A[grep("^1",res4$category)]) # >0 *
t.test(res4$A[grep("^2",res4$category)]) # =0
t.test(res4$A[grep("^3",res4$category)]) # =0
t.test(res4$A[grep("^4",res4$category)]) # >0 *
t.test(res4$A[grep("^5",res4$category)]) # >0
t.test(res4$A[grep("^6",res4$category)]) #<0
t.test(res4$A[grep("^7",res4$category)])
# |B|

# B
t.test(res4$B[grep("^1",res4$category)]) # <0 *
t.test(res4$B[grep("^2",res4$category)]) # =0
t.test(res4$B[grep("^3",res4$category)]) # =0
t.test(res4$B[grep("^4",res4$category)]) # <0 *
t.test(res4$B[grep("^5",res4$category)]) # <0 *
t.test(res4$B[grep("^6",res4$category)]) # >0
t.test(res4$B[grep("^7",res4$category)]) # <0
