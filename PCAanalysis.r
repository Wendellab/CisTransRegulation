fileL <- grep("count.txt",list.files(), value=TRUE)
for(i in fileL)
{
    t<-read.table(i,sep="\t");
    if(exists("total"))
    {total<-cbind(total,t$V3)}else
    {total=t[,c(1,3)]; names(total)[1]="gene"}
    names(total)[ncol(total)]<-i
}
rownames(total)<-total$gene
total<-total[-which(rownames(total)=="*"),-1]

##### Prepare column info
info<- data.frame(sample=names(total), lib_size=colSums(total) )
# read in rename list
rinfo <- read.table("renameInfo.txt")
info$new <- rinfo[rownames(info),]
info<-cbind(info,t(as.data.frame(strsplit(as.character(info$new),"-"))))
names(info)[4:7] <-c("genome","source","dpa","rep")

## deseq normalization
library(DESeq2)
dds <- DESeqDataSetFromMatrix( countData = total, colData = info, design = ~ new)

# rlog transformation, note that with defalt blind=TRUE, design is actually ~1
rld <- rlog(dds)
count.rld <- as.data.frame(assay(rld))
names(count.rld) == names(total)

# also possible to perform custom transformation:
dds <- estimateSizeFactors(dds)
# shifted log of normalized counts
se <- SummarizedExperiment(log2(counts(dds, normalized=TRUE) + 1), colData=colData(dds))
count.log2 <- as.data.frame(assay(se))

# save
save(info, total, count.rld, count.log2, file="counts.rdata")

## noting that plotPCA code is too simple to show multiple variables
# the call to DESeqTransform() is needed to trigger our plotPCA method.
plotPCA( DESeqTransform( se ), intgroup =c("dpa") )
plotPCA(rld,intgroup =c("dpa"))

# my custom function
library(ggplot2)
library(scales)
library(ape)
plotGrouping <- function(norm_log, color, shape, text, tip, save = "plotGrouping.pdf"){
    # norm<-sweep(total,2,info$lib_size,"/")*10^6
    # norm_log <- log2(norm+1)
    pca=prcomp(t(norm_log))
    dat = as.data.frame(pca$x)
    proportion<-summary(pca)$importance[2,1:2]*100
    proportion<-paste0(names(proportion)," (", proportion, "%)")
    p<-ggplot(aes(PC1, PC2, color=color,    shape=shape),data=dat) + geom_point() +xlab(proportion[1]) + ylab(proportion[2])
    pdf(save)
    print( p + geom_text(aes_string(x = "PC1", y = "PC2", label = text), color="grey", hjust = 0, nudge_x = 0.09) )
    
    hc<-hclust( dist(t(norm_log)) )
    tre <- as.phylo(hc)
    tre$tip.label <- as.character(tip)
    # try to match ggplot color: library(scale);
    # show_col(col4<-hue_pal()(4))
    tipCol <- color
    levels(tipCol) <- hue_pal()(nlevels(tipCol))
    plot(tre,tip.col =as.character(tipCol),  type="unrooted",cex=0.6, no.margin=TRUE)
    plot(tre,tip.col =as.character(tipCol), type="fan",cex=0.6, no.margin=TRUE)
    dev.off()   }

# plots
plotGrouping(count.rld, color=info$dpa, shape=info$source, tip=info$new, text=info$genome, save = "plotGrouping.rld.pdf")
plotGrouping(count.log2, color=info$dpa, shape=info$source, text=info$genome, save = "plotGrouping.log2.pdf")

# restrict to 10 and 20 dpa
select<-which((info$genome=="Maxxa"|info$genome=="TX2094"|info$genome=="Yuc"|info$genome=="F1") & (info$dpa=="10dpa"|info$dpa =="20dpa"))
plotGrouping(count.rld[,select], color=info$dpa[select], shape=info$source[select], text=info$genome[select], tip=info$new[select], save = "plotGrouping.1020.rld.pdf")
plotGrouping(count.log2[,select], color=info$dpa[select], shape=info$source[select], text=info$genome[select], tip=info$new[select], save = "plotGrouping.1020.log2.pdf")

select<-which( (info$source=="yb"|info$source=="jj")& (info$genome=="Maxxa"|info$genome=="TX2094"|info$genome=="Yuc"|info$genome=="F1") & (info$dpa=="10dpa"|info$dpa =="20dpa"))
plotGrouping(count.rld[,select], color=info$dpa[select], shape=info$source[select], text=info$genome[select], tip=info$new[select], save = "plotGrouping.1020.rld.yj.pdf")
plotGrouping(count.log2[,select], color=info$dpa[select], shape=info$source[select], text=info$genome[select], tip=info$new[select], save = "plotGrouping.1020.log2.yj.pdf")


