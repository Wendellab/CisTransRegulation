## Compare regulatory pattern between homoeologs
library(gplots)
library(ggplot2)
options(scipen=999)

## check subgenome
countADZ=function(x) {table(gsub("[0-9].*","",gsub("Gohir.1|Gohir.|_.*","",x)))}
# countADZ(ids)
#     A     D    mt     Z
# 13462 14060     9   285



## Homoeolog expression bias ##
library(DESeq2)
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

## plot correspondecne
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
    
    corrplot(residualTable, p.mat=as.matrix(pTable), insig="label_sig", sig.level=c(0.001,0.01,0.05),is.corr=FALSE, method="color",title=title, pch.col="white",tl.col="black")
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
