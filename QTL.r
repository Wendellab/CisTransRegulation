options(scipen=999)

# fiber QTL
# /work/LAS/jfw-lab/corrinne/QTLpaper/fiber.qtl.genes.results
x<-read.table("/work/LAS/jfw-lab/corrinne/QTLpaper/fiber.qtl.gene.results", sep="\t")
g<- unique(paste0(gsub(".*Name=","",x$V12),".1"))
length(g)  # 6652

#
load("RDgenes.rdatadim")
rd$qtl = rd$ID %in% g
table(rd$qtl)  # 2660 qtl genes detected with SNPs, 0.09562842 out of 27816
table(rd[,c("qtl","rd")])
#       rd
# qtl     cis and trans cis only no divergence trans only
#  FALSE           742      462         23689        263
#  TRUE             99       51          2472         38
fisher.test(rd$qtl,rd$rd!="no divergence",alternative="greater") # p=0.006767
fisher.test(rd$qtl,rd$crd,alternative="greater") # p=0.03081
fisher.test(rd$qtl,rd$trd,alternative="greater") # p=0.003204
chisq.test(table(rd[,c("qtl","rd")])) # p=0.03701


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

> plotFisherCor(rd$rd,rd$qtl)
Error in plotFisherCor(rd$rd, rd$qtl) :
could not find function "plotFisherCor"
> plotFisherCorr(rd$rd,rd$qtl)
Error in corrplot(residualTable, p.mat = as.matrix(pTable), insig = "label_sig",  :
could not find function "corrplot"
> library(corrplot)
corrplot 0.84 loaded
> plotFisherCorr(rd$rd,rd$qtl)
> dev.off()
null device
1

