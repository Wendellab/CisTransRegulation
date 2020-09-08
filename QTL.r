options(scipen=999)
load("RDgenes.rdata")

# fiber QTL
# /work/LAS/jfw-lab/corrinne/QTLpaper/fiber.qtl.genes.results
x<-read.table("/work/LAS/jfw-lab/corrinne/QTLpaper/fiber.qtl.gene.results", sep="\t")
g<- unique(paste0(gsub(".*Name=","",x$V12),".1"))
length(g)  # 6652

# enrichment of QTL genes in RD? Kinda
rd$qtl = rd$ID %in% g
table(rd$qtl)  # 2660 qtl genes detected with SNPs, 0.09562842 out of 27816
table(rd$qtl[rd$rd!="no divergence"])  # 188 qtl genes in RD, 0.1135952 out of 27816
fisher.test(rd$qtl,rd$rd!="no divergence",alternative="greater") # p=0.006767

# then which type of RD most enriched with QTL? None
table(rd[,c("qtl","rd")])
#       rd
# qtl     cis and trans cis only no divergence trans only
#  FALSE           742      462         23689        263
#  TRUE             99       51          2472         38
chisq.test(table(rd[,c("qtl","rd")])) # p=0.03701
chisq.test(table(rd[rd$rd!="no divergence",c("qtl","rd")])) # p-value = 0.4394, no preference
fisher.test(rd$qtl,rd$crd,alternative="greater") # p=0.03081
fisher.test(rd$qtl,rd$trd,alternative="greater") # p=0.003204

source("utilities.r")
pdf("QTLvsRD.pdf")
plotFisherCorr(rd$rd!="no divergence",rd$qtl)
plotFisherCorr(rd$rd[rd$rd!="no divergence"],rd$qtl[rd$rd!="no divergence"])

dev.off()