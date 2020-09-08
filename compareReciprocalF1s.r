## Categoryzation of cis tran
library(gplots)
library(ggplot2)
library(corrplot)
library(WGCNA)
load("cistrans.rdata")

# A
unique(res.mt10$A- res.tm10$A) #  0 NA
unique(res.mt20$A- res.tm20$A) #  0 NA

# B
cor(res.mt10$B, res.tm10$B, use="complete.obs") # 0.6495133
cor(res.mt20$B, res.tm20$B, use="complete.obs") # 0.6019386

# A-B
cor(res.mt10$AminusB, res.tm10$AminusB, use="complete.obs") # 0.672366
cor(res.mt20$AminusB, res.tm20$AminusB, use="complete.obs") # 0.6688194

# gene of same category
table(res.mt10$category==res.tm10$category)
table(res.mt20$category==res.tm20$category)

pdf("compareReciprocalF1s.pdf")

for(i in c(10,20))
{
    print(i)
    mt = get(paste0("res.mt",i))$category
    tm = get(paste0("res.tm",i))$category
    tbl=table(mt,tm)
    print(sum(diag(tbl))/sum(tbl))
    print(sum(diag(tbl[1:5,1:5]))/sum(tbl[1:5,1:5]))
    print(sum(diag(tbl[1:6,1:6]))/sum(tbl[1:6,1:6]))
    
    textplot(tbl); title(paste0(i," dpa: p-value ", chisq.test(tbl,simulate.p.value = TRUE)$p.value))
    # Pearson's Chi-squared test
    chisq <- chisq.test(tbl)
    print(chisq) # the row and the column variables are statistically significantly associated (p-value = 0).
    ## below plots not very good
    # balloonplot(tbl, show.margin=F, main = i,label=F)
    # assoc(tbl, shade=TRUE, main = i)
    # mosaic(tbl, shade=TRUE, main = i)
    ## plot matrix
    CountTbl = as.matrix(tbl)
    residualTable  = chisq.test(tbl,simulate.p.value = TRUE)$residuals
    corrplot(residualTable, p.mat=CountTbl, insig="p-value", sig.level=0,is.corr=FALSE, main=i)
    pTable = CountTbl
    rowCats = levels(factor(mt))
    colCats = levels(factor(tm))
    for(rc in rowCats)
    for(cc in colCats)
    {
        rclass <- (mt==rc) # row
        cclass <- (tm==cc) # col
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
    xLabels = paste("TxM - ", gsub("[.].*","",colCats)), yLabels = paste("MxT - ", gsub("[.].*","",rowCats)),
    textMatrix = CountTbl, colors = blueWhiteRed(100)[50:100],
    main = paste0(i, " dpa: Correspondence of categories"),
    cex.text = 1, cex.lab = 1, setStdMargins = FALSE      )
}
dev.off()

# make RD lists
rd = data.frame(ID=rownames(res.mt10), mt10 = res.mt10$category, tm10 = res.tm10$category, mt20 = res.mt20$category, tm20 = res.tm20$category)
rd$crd = (grepl("^[1,3-5]",rd$mt10)&rd$mt10==rd$tm10) | (grepl("^[1,3-5]",rd$mt20)&rd$mt20==rd$tm20) # restrict to overlap between reciprocal F1s and from category I to V at either time point
rd$trd = (grepl("^[2,3-5]",rd$mt10)&rd$mt10==rd$tm10) | (grepl("^[2,3-5]",rd$mt20)&rd$mt20==rd$tm20)# restrict to overlap between reciprocal F1s and from category I to V at either time point
with(rd,table(crd,trd))
#                 trd
# crd     FALSE  TRUE
#   FALSE 26161   301
#   TRUE    513   841
rd$rd="no divergence"
rd$rd[rd$crd==TRUE & rd$trd==FALSE]="cis only"
rd$rd[rd$trd==TRUE & rd$crd==FALSE]="trans only"
rd$rd[rd$trd==TRUE & rd$crd==TRUE]="cis and trans"


save(rd,file="RDgenes.rdata")
