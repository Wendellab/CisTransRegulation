#####################################################
################# File Discription ##################
#####################################################
library(DESeq2)
library(gtools)
# H.B. Fraser. Improving estimates of compensatory cis-trans regulatory divergence, Trends in Genetics2019
# Fraser (2019) pointed out the ASE method is intricically biased for comparing cis- and trans- regulatory divergence: because trans = parental - cis, any error in estimating cis effects will introduce an artifactual negative correlation with trans effects, hence great error in ASE estimates will lead to stronger cis-trans correlations. One solution proposed to eliminate this bias is to use different biological replicates for estimating cis and trans effects (termed as "cross-replicate comparison"), thus their errors are ramdom to each other.

## Based on the standdard approach, our previous estimates of cis an trans effects showed strong megative correlation ranging from Pearson's r= -0.76 to -0.80
load("cistrans.rdata")
with(res.mt10, cor(B, AminusB,use="pairwise.complete.obs")) # -0.7976486
with(res.tm10, cor(B, AminusB,use="pairwise.complete.obs")) # -0.7860334
with(res.mt20, cor(B, AminusB,use="pairwise.complete.obs")) # -0.7563041
with(res.tm20, cor(B, AminusB,use="pairwise.complete.obs")) # -0.7759683

## Before conducting cross-replicate comparison for each reciprocal hybrid, the cis and trans estimates between reciprocal hybrids were independently estimated and exbited weaker correlations. But these negative correlationsare still relatively strong.
cor(res.mt10$B, res.tm10$AminusB,use="pairwise.complete.obs") # -0.4474439
cor(res.tm10$B, res.mt10$AminusB,use="pairwise.complete.obs") # -0.4650338
cor(res.mt20$B, res.tm20$AminusB,use="pairwise.complete.obs") # -0.4055039
cor(res.tm20$B, res.mt20$AminusB,use="pairwise.complete.obs") # -0.4063479

## Next cross-replicate comparison
count = read.table(file="Genes27816.raw.count.txt", sep="\t",header=T)
coldata  = read.table( file="Genes27816.sample.info.txt", sep="\t",header=T)
coldata$sample= gsub("maxxa.|tx2094.","",coldata$condition2)
count=count[,coldata$allele!="total"]
coldata=coldata[coldata$allele!="total",]

# rlog transformation mirrors standard method
dds <- DESeqDataSetFromMatrix( countData = count, colData = coldata, design = ~ condition2)
rld <- rlog(dds)
r <- assay(rld)
coldata$sample= gsub("maxxa.|tx2094.","",coldata$condition2)
rcount=r[,coldata$allele!="total"]
coldata=coldata[coldata$allele!="total",]
# out of three biological replicates, pick one to estimate cis and another to estimate trans, then calculate correlation; repeat this for all six combinationslibrary(gtools)
rm(sumRes)
for(id in unique(coldata$sample)){
    cc = count[,coldata$sample==id]
    co = coldata[coldata$sample==id,]
    dds <- DESeqDataSetFromMatrix( countData = cc, colData = co, design = ~ allele)
    # standard method, all resps use
    resAll <- results(DESeq(dds), contrast = c("allele","maxxa","tx2094"))
    B = resAll$log2FoldChange
    A=get(paste0("A",gsub(".*[.]|dpa","",id)))$log2FoldChange
    AminusB = A-B
    corr=cor.test(B, AminusB)
    cisP =mean(abs(B)/(abs(B) +abs(AminusB)),na.rm=T)
    wilcox =wilcox.test(abs(B),abs(AminusB))
    corrCis = cor.test(B, A)
    corrTrans = cor.test(AminusB, A)
    res=data.frame(sample=id, cisRep=paste0(unique(co$rep),collapse=","), transRep=paste0(unique(co$rep),collapse=","), r=corr$estimate,pvalue=corr$p.value, absCisP=cisP, pvalue.wilcox=wilcox$p.value, r.BvsA=corrCis$estimate,pvalue.BvsA=corrCis$p.value,r.AminusBvsA=corrTrans$estimate,pvalue.AminusBvsA=corrTrans$p.value )
    if(exists("sumRes")){sumRes=rbind(sumRes,res)}else{sumRes=res}
    # permutation of bio rep
    combns=permutations(3,2,unique(co$rep))
    for(i in 1:nrow(combns))
    {
        # cis by combn[i,1]
        # trans by combn[i,2]
        cisC = cc[,co$rep==combns[i,1]]
        cisD = co[co$rep==combns[i,1],]
        transC = cc[,co$rep==combns[i,2]]
        transD =co[co$rep==combns[i,2],]
        # get B for cis
        B = log2(cisC[,cisD$allele=="maxxa"]/cisC[,cisD$allele=="tx2094"])
        # get B for trans
        Bt = log2(transC[,transD$allele=="maxxa"]/transC[,transD$allele=="tx2094"])
        # get A
        A=get(paste0("A",gsub(".*[.]|dpa","",id)))$log2FoldChange
        AminusB = A-Bt
        s= which(is.finite(B)&is.finite(Bt))
        corr=cor.test(B[s], AminusB[s])
        cisP =mean(abs(B)/(abs(B) +abs(AminusB)),na.rm=T)
        wilcox =wilcox.test(abs(B),abs(AminusB))
        corrCis = cor.test(B[s], A[s])
        corrTrans = cor.test(AminusB[s], A[s])
        res=data.frame(sample=id, cisRep=as.character(combns[i,1]), transRep=as.character(combns[i,2]), r=corr$estimate,pvalue=corr$p.value, absCisP=cisP, pvalue.wilcox=wilcox$p.value, r.BvsA=corrCis$estimate,pvalue.BvsA=corrCis$p.value,r.AminusBvsA=corrTrans$estimate,pvalue.AminusBvsA=corrTrans$p.value )
        if(exists("sumRes")){sumRes=rbind(sumRes,res)}else{sumRes=res}
    }
}
sumRes
write.table(sumRes,file="correlation.fraser.txt",row.names=FALSE, sep="\t")
