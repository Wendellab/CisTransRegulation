# installation
# module load r-udunits2
# module load r/3.5.0-py2-ufvuwmm
if (!requireNamespace("BiocManager", quietly = TRUE))
install.packages("BiocManager")
BiocManager::install("clusterProfiler", version = "3.8")
BiocManager::install("AnnotationHub")
BiocManager::install("KEGG.db", version = "3.8")
BiocManager::install("KEGGREST", version = "3.8")
BiocManager::install("enrichplot", version = "3.8")

# load libaries
library("AnnotationHub")
library("GO.db")
library("clusterProfiler")
library("enrichplot")


library("biomaRt")
library("topGO")
library("Rgraphviz")
library("pathview")

########################
# check AnnotationHub ##
########################
library("AnnotationHub")
hub <- AnnotationHub::AnnotationHub()
query(hub,"Gossypium")
# AnnotationHub with 6 records
# snapshotDate(): 2018-10-24
# $dataprovider: ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/
# $species: Gossypium arboreum, Gossypium hirsutum, Gossypium hirsutum subsp...
# $rdataclass: OrgDb
# additional mcols(): taxonomyid, genome, description,
#   coordinate_1_based, maintainer, rdatadateadded, preparerclass, tags,
#   rdatapath, sourceurl, sourcetype
# retrieve records with, e.g., 'object[["AH66190"]]'
#
# title
# AH66190 | org.Gossypium_hirsutum.eg.sqlite
# AH66191 | org.Gossypium_hirsutum_subsp._mexicanum.eg.sqlite
# AH66192 | org.Gossypium_lanceolatum.eg.sqlite
# AH66193 | org.Gossypium_purpurascens.eg.sqlite
# AH66244 | org.Gossypium_raimondii.eg.sqlite
# AH66259 | org.Gossypium_arboreum.eg.sqlite
orgdb <- hub[["AH66190"]]
# list column and key types
columns(orgdb)
keytypes(orgdb)
# check content
head(keys(orgdb, "PMID"))
head(keys(orgdb, "ACCNUM"))
head(keys(orgdb, "PMID"))
id=head(keys(orgdb, "UNIGENE"))
select(orgdb, id, c("SYMBOL", "GENENAME"), "UNIGENE")
# 'select()' returned 1:1 mapping between keys and columns
# UNIGENE       SYMBOL                                      GENENAME
# 1 Ghi.10022 LOC107905551 24-methylenesterol C-methyltransferase 2-like
# 2 Ghi.10025 LOC107958018           heat shock cognate 70 kDa protein 2
# 3 Ghi.10030 LOC107934208                       protein TONNEAU 1a-like
# 4 Ghi.10068 LOC107951247         SNAP25 homologous protein SNAP33-like
# 5 Ghi.10092 LOC107925933         G2/mitotic-specific cyclin S13-7-like
# 6 Ghi.10094 LOC107898859      glutamine synthetase nodule isozyme-like
## KEGG ##
library(KEGGREST)
listDatabases()
# "pathway"  "brite"    "module"   "ko"       "genome"   "vg"
# "ag"       "compound" "glycan"   "reaction" "rclass"   "enzyme"
# "disease"  "drug"     "dgroup"   "environ"  "genes"    "ligand"
# "kegg"
org <- as.data.frame(keggList("organism"))
org[grep("Gossypium",org$species),]
# T.number organism                            species               phylogeny
# 201   T04129      gra                Gossypium raimondii       Eukaryotes;Plants;Eudicots;Mallow family
# 202   T04793      ghi Gossypium hirsutum (upland cotton)       Eukaryotes;Plants;Eudicots;Mallow family
keggList("ghi")[1:3]
# ghi:107963937
# "E3 ubiquitin-protein ligase CCNB1IP1 homolog isoform X2"
# ghi:107885929
# "THO complex subunit 4A-like"
# ghi:107885930
# "ACT domain-containing protein ACR9-like"
##### NOT the reference genome i am using


#########################
## clusterProfile prep ##
#########################
# plots: http://bioconductor.org/packages/devel/bioc/vignettes/enrichplot/inst/doc/enrichplot.html

# prep reference database
library(clusterProfiler)
library(tidyr)
desc <- read.table("Genes66610.transcript.description.txt",header=TRUE,sep="")
gomap <- desc[,c("GO","transcriptName")]
gomap<-separate_rows(gomap, GO,sep=",")
godb<-buildGOmap(gomap)
ont=go2ont(godb$GO)
godb$ont = ont$Ontology[match(godb$GO,ont$go_id)]

# prep GO term2name
library(GO.db)
goterms <- Term(GOTERM)
goterm2name = data.frame(term=names(goterms),name=goterms)

# prep GO enrichment function, plot top 10 enriched each for MF, CC and BP and save table
GOenrich = function(geneL, refL, filename="GOenrich"){
    if(!is.null(filename)){pdf(paste0(filename,".pdf"))}
    for(o in c("MF","BP","CC")){
        en <- enricher(geneL, pvalueCutoff = 0.05, pAdjustMethod = "BH", universe=refL, minGSSize = 10, maxGSSize = 500, qvalueCutoff = 0.05, TERM2GENE = godb[godb$ont==o,1:2], TERM2NAME = goterm2name)
        # print(en)
        print( dotplot(en, showCategory=10,font.size = 6, title=o) )
        if(nrow(en)>0){
            print( emapplot(en, showCategory=30,font.size = 6, title=o) )
            df=as.data.frame(en)
            df$ont=o
            df$level2 = df$ID %in% getGOLevel(o, 2)
            df$level3 = df$ID %in% getGOLevel(o, 3)
            df$level4 = df$ID %in% getGOLevel(o, 4)
            if(o=="MF"){dfL=df}
            else if(!exists("dfL")){dfL=df}
            else{dfL=rbind(dfL,df)}
        }
    }
    if(!is.null(filename)){
        dev.off()
        write.csv(dfL,file=paste0(filename,".csv"))}
    print(table(dfL$ont))
    return(dfL)
}
GOcompare = function(geneL, refL, filename="GOcompare"){
    if(!is.null(filename)){pdf(paste0(filename,".pdf"))}
    for(o in c("MF","BP","CC")){
        print(o);
        en=tryCatch(compareCluster(geneL, fun = "enricher", pvalueCutoff = 0.05, pAdjustMethod = "BH", universe=refL, minGSSize = 10, maxGSSize = 500, qvalueCutoff = 0.05, TERM2GENE = godb[godb$ont==o,1:2], TERM2NAME = goterm2name), error=function(e){cat("ERROR :",conditionMessage(e), "\n")}, finally={en=NULL})
        if(!is.null(en)){
            print( dotplot(en, showCategory=10,font.size = 6, title=o) )
            df=as.data.frame(en)
            df$ont=o
            df$level2 = df$ID %in% getGOLevel(o, 2)
            df$level3 = df$ID %in% getGOLevel(o, 3)
            df$level4 = df$ID %in% getGOLevel(o, 4)
            if(!exists("dfL")){dfL=df}else{dfL=rbind(dfL,df)}
        }
    }
    if(!exists("dfL")){dfL=NULL}
    if(!is.null(filename)){
        dev.off()
        write.csv(dfL,file=paste0(filename,".csv"))}
    print(table(dfL[,c("ont","Cluster")]))
    return(dfL)
}
GOcompare.formula = function(formula, data, filename="GOcompare"){
    if(!is.null(filename)){pdf(paste0(filename,".pdf"))}
    for(o in c("MF","BP","CC")){
        en <- compareCluster(formula, data, fun = "enricher", pvalueCutoff = 0.05, pAdjustMethod = "BH", minGSSize = 10, maxGSSize = 500, qvalueCutoff = 0.05, TERM2GENE = godb[godb$ont==o,1:2], TERM2NAME = goterm2name)
        print( dotplot(en, showCategory=10,font.size = 6, title=o) )
        df=as.data.frame(en)
        df$ont=o
        df$level2 = df$ID %in% getGOLevel(o, 2)
        df$level3 = df$ID %in% getGOLevel(o, 3)
        df$level4 = df$ID %in% getGOLevel(o, 4)
        if(o=="MF"){dfL=df}else{dfL=rbind(dfL,df)}
    }
    if(!is.null(filename)){
        dev.off()
        write.csv(dfL,file=paste0(filename,".csv"))}
    print(table(dfL[,c("ont","Cluster")]))
    return(dfL)
}

result=GOcompare.formula(ID~rd,data=rd,filename="GOcompare.RD")

# dependent functions
getGOLevel <- function(ont, level) {
    switch(ont,
    MF = {
        topNode <- "GO:0003674"
        Children <- GOMFCHILDREN
    },
    BP = {
        topNode <- "GO:0008150"
        Children <- GOBPCHILDREN
    },
    CC = {
        topNode <- "GO:0005575"
        Children <- GOCCCHILDREN
    }
    )
    
    max_level <- max(level)
    if (any(level == 1)) {
        all_nodes <- topNode
    } else {
        all_nodes <- c()
    }
    
    Node <- topNode
    for (i in seq_len(max_level-1)) {
        Node <- mget(Node, Children, ifnotfound=NA)
        Node <- unique(unlist(Node))
        Node <- as.vector(Node)
        Node <- Node[!is.na(Node)]
        if ((i+1) %in% level) {
            all_nodes <- c(all_nodes, Node)
        }
    }
    return(all_nodes)
}
dotplot.df=function(result, title="")
{
    library(DOSE) # theme_dose
    parse_ratio <- function(ratio) {
        gsize <- as.numeric(sub("/\\d+$", "", as.character(ratio)))
        gcsize <- as.numeric(sub("^\\d+/", "", as.character(ratio)))
        return(gsize/gcsize)
    }
    ratio=parse_ratio(result$GeneRatio)
    # give color to ont
    a <- ifelse(result$ont == "MF", "gold", ifelse(result$ont == "BP","steelblue","salmon"))[order(result$Description)]
    ggplot(result, aes_string(x="Cluster", y="Description", size=ratio, color="p.adjust")) +
    geom_point() +
    scale_color_continuous(low="red", high="blue", name = "p.adjust", guide=guide_colorbar(reverse=TRUE)) +
    #scale_color_gradientn(name = "p.adjust", colors=sig_palette, guide=guide_colorbar(reverse=TRUE)) +
    ylab(NULL) + ggtitle(title) + theme_dose(font.size=6) + scale_size(range=c(3, 8))+
    theme(axis.text.y = element_text(colour = a))
    #+ facet_grid(ont~.)
}


# save into rdata
save(desc, godb, goterm2name,GOenrich, GOcompare, getGOLevel, dotplot.df, file="funcEnrich.rdata")


##############
## Analysis ##
##############
library(clusterProfiler)
library(tidyr)
library(GO.db)

load("funcEnrich.rdata")

## 1. fiber expressed genes
# load total counts
count <-read.table("totalReadCounts.txt", header=TRUE, sep="\t")
# requiring a read depth of 5 per genes, how many expressed genes in each accession?
# Union of expressed fiber genes
geneL=rownames(count)[which(rowSums(count[,grep("Maxxa",names(count))])>=5 | rowSums(count[,grep("TX2094",names(count))])>=5 | rowSums(count[,grep("F1.*[1-3]$",names(count))])>=5 |rowSums(count[,grep("F1.*[4-6]$",names(count))])>=5)]
refL=rownames(count)
result=GOenrich(geneL,refL,filename="GOenrich.fiberExpressed")
#  BP  CC  MF
# 237  60  59
head(result[order(result$qvalue,decreasing=FALSE),-8],10)
result[result$level2==TRUE,-8]
result[result$level3==TRUE & result$ont=="MF",-8]

## 2. expression changes between TX2094 and Maxxa, 10 dpa vs 20 dpa
r10 <-read.table("DE-66610/Maxxa.10dpavsTX2094.10dpa.txt", header=TRUE, sep="\t")
r20 <-read.table("DE-66610/Maxxa.20dpavsTX2094.20dpa.txt", header=TRUE, sep="\t")
length(which(r10$padj<0.05&!is.na(r10$padj)))/66610  # 6883, 10.3%
length(which(r20$padj<0.05&!is.na(r20$padj)))/66610  #4829, 7.2%
length(which( (r10$padj<0.05&!is.na(r10$padj))|(r20$padj<0.05&!is.na(r20$padj)) ))/66610  # 9921, 15.0
r10id=rownames(r10)[r10$padj<0.05&!is.na(r10$padj)]
r20id=rownames(r20)[r20$padj<0.05&!is.na(r20$padj)]
geneL=list(dpa10=r10id,dpa20=r20id)
refL=rownames(count)[which(rowSums(count[,grep("Maxxa",names(count))])>=5 | rowSums(count[,grep("TX2094",names(count))])>=5 | rowSums(count[,grep("F1.*[1-3]$",names(count))])>=5 |rowSums(count[,grep("F1.*[4-6]$",names(count))])>=5)] #55551
result=GOcompare(geneL,refL,filename="GOcompare.MaxxavsTX2094")
#     Cluster
# ont  dpa10 dpa20
# BP    58    30
# CC     0    53
# MF    16    22
result[result$level3==TRUE,-9]
pdf("GOcompare.MaxxavsTX2094.all.pdf")
dotplot.df(result, title="GO enrichment of Maxxa vs TX2094 DE genes")
dev.off()

## 3. 27816 genes
load("RDgenes.rdata")
# query gene list
geneL = rownames(rd) # 27816
refL=rownames(count)[which(rowSums(count[,grep("Maxxa",names(count))])>=5 | rowSums(count[,grep("TX2094",names(count))])>=5 | rowSums(count[,grep("F1.*[1-3]$",names(count))])>=5 |rowSums(count[,grep("F1.*[4-6]$",names(count))])>=5)] #55551
result=GOenrich(geneL,refL,filename="GOenrich.27816")
#  BP  CC  MF
# 202  61 54
head(result[order(result$qvalue,decreasing=FALSE),-8],10)
result[result$level2==TRUE,-8]
result[result$level3==TRUE & result$ont=="MF",-8]

## 4. rd genes
geneL = rownames(rd)[rd$rd!="no divergence"] #1655
refL = rownames(rd) # 27816
result=GOenrich(geneL,refL,filename="GOenrich.RD")
# CC MF
# 24 17
result=GOcompare.formula(ID~rd,data=rd[rd$rd!="no divergence",],filename="GOcompare.RD")
#     Cluster
# ont  cis and trans cis only no divergence trans only
#   BP            10        0           277          1
#   CC             0       29            83         15
#   MF             1        1            97          4
pdf("GOcompare.RD.all.pdf")
dotplot.df(result, title="GO enrichment of RD genes")
dev.off()

## 5. cluster
x<-read.table("coexpressionNet/Clust_Results_30_Nov_18/Clusters_Objects.tsv",header=TRUE,sep="\t")
x<-x[-1,]
library(reshape2)
y=melt(x,measure.vars=names(x))
y=y[y$value!="",]
y$cluster = gsub("[.].*","",y$variable)
result=GOcompare.formula(value~cluster,data=y,filename="GOcompare.Clust")
# Cluster
# ont   C0  C1 C10 C11 C12 C13 C14 C15 C16 C17 C18  C2  C3  C4  C5  C6  C7  C8  C9
#   BP  17   0   1   0   1   1 100  23   0   0   4   0   0  10 114  11  14   5  21
#   CC  13   0   0   0   0   9  67  33   5   0   0   0   0   0  49   6   5   1  12
#   MF   5   0   0   4   0   1  65  13   0   0   0   6   0   5  39   6   3   9   0
## NOT making much sense because the universe is 9245 genes assigned to clusters; no enrichment when against all genes


##########
## GSEA ##
##########
## example
data(geneList, package="DOSE")
gmtfile <- system.file("extdata", "c5.cc.v5.0.entrez.gmt", package="clusterProfiler")
c5 <- read.gmt(gmtfile)
egmt2 <- GSEA(geneList, TERM2GENE=c5, verbose=FALSE)
head(egmt2)
# visualize value distributions of core enriched genes for GSEA enriched categories.
ridgeplot(egmt2)
#  visualize the distribution of the gene set and the enrichment score.
gseaplot(egmt2, geneSetID = "EXTRACELLULAR_REGION",  title = egmt2$Description[1])
gseaplot2(egmt2, geneSetID = "EXTRACELLULAR_REGION",  title = egmt2$Description[1])

## FUN
makeGSEAlist = function(value,name,decreasing=TRUE){
    #https://github.com/GuangchuangYu/DOSE/wiki/how-to-prepare-your-own-geneList
    ## feature 1: numeric vector
    geneList = value
    ## feature 2: named vector
    names(geneList) = name
    ## feature 3: decreasing order
    geneList = sort(geneList, decreasing = decreasing)
    return(geneList)
}

### rank genes by connectivity and module membership
load("coexpressionNet/buildNetwork.rdata")

# 1. given GO to see what terms are enriched in hub
pdf("GSEAnetwork.GO.pdf")
for(net in c("k.tx2094","k.maxxa")){
    k=get(net)
    for(ktype in c("kTotal","kWithin")){
        el=GSEA(geneList=makeGSEAlist(value=k[,ktype], name=rownames(k.tx2094), decreasing=TRUE), TERM2GENE=godb, TERM2NAME = goterm2name)
        print(dotplot(el, showCategory=20,font.size = 6, title=paste(net,ktype)) )
        print(emapplot(el, showCategory=30,font.size = 6, title=paste(net,ktype)) )
        print(ridgeplot(el,showCategory=30)+ggplot2::theme_bw(base_size = 6) )
    }
}
dev.off()
# no enrichment found for kWithin in Maxxa network

# 2. given cis/trans regulatory categories to see what terms are enriched in hub
load("RDgenes.rdata")
# prep gene sets
library(reshape2)
gs=melt(rd,value.name=c("ID"), measure.vars=c("rd","mt10","tm10","mt20","tm20","crd","trd"))
rd2gene = data.frame(category=paste0(gs[,2],".",gs[,3]),ID=gs[,1])
# restrict network to 27816 genes
select = as.character(rd$ID)
pdf("GSEAnetwork.RD.pdf")
for(net in c("k.tx2094","k.maxxa")){
    k=get(net)
    for(ktype in c("kTotal","kWithin")){
        print(paste(net,ktype))
        el=GSEA(geneList=makeGSEAlist(value=k[select,ktype], name=select, decreasing=TRUE), TERM2GENE=rd2gene,maxGSSize = 1500,minGSSize=5)
        print(dotplot(el, showCategory=20,font.size = 6, title=paste(net,ktype)) )
        print(emapplot(el, showCategory=30,font.size = 6, title=paste(net,ktype)) )
        print(ridgeplot(el,showCategory=30)+ggplot2::theme_bw(base_size = 6) )
        for(i in grep("^rd",el$ID))
        {
            print( gseaplot(el, geneSetID = i,  title = el$ID[i]) )
        }
    }
}
dev.off()
########
# for TX2094 network, cis only, trans only and cis+trans are all enriched in hub. looking good and suggest that domestication introduce regulatory changes to network hubs. BUT WHY no enrichment with Maxxa network. since GSEA is rank-based, the distribution difference in Maxxa vs TX2094 k doesn't matter (same result after quantile normalization).
#


# 3. given inheritance categories to see what terms are enriched in hub
load("inheritance.rdata")
# prep gene sets
rd2gene = data.frame(category=c(paste0("mt10.",im.mt10$category),paste0("tm10.",im.tm10$category),paste0("mt20.",im.mt20$category),paste0("tm20.",im.tm20$category)), ID=rep(rownames(im.mt10),4))
# restrict network to 27816 genes
select = rownames(im.mt10)
pdf("GSEAnetwork.inheritance.pdf")
for(net in c("k.tx2094","k.maxxa")){
    k=get(net)
    for(ktype in c("kTotal","kWithin")){
        print(paste(net,ktype))
        el=GSEA(geneList=makeGSEAlist(value=k[select,ktype], name=select, decreasing=TRUE), TERM2GENE=rd2gene,maxGSSize = 1500,minGSSize=5)
        print(dotplot(el, showCategory=20,font.size = 6, title=paste(net,ktype)) )
        print(emapplot(el, showCategory=30,font.size = 6, title=paste(net,ktype)) )
        print(ridgeplot(el,showCategory=30)+ggplot2::theme_bw(base_size = 6) )
        for(i in grep("^rd",el$ID))
        {
            print( gseaplot(el, geneSetID = i,  title = el$ID[i]) )
        }
    }
}
dev.off()
# more enriched categories found in TX2094 again. Not sure why. more modules 32 vs 22 in TX2094 vs Maxxa, save was found before by joeG 108 vs 47 modules.

##########################################
## JJ's cell wall and cytoskeleton genes##
##########################################
jj=read.table("/work/LAS/jfw-lab/erdostal/JJgorai.txt") #703 gorai
names(jj) ="Gorai"
load("homoeoPairs.rdata")
dim(cw<-merge(jj,ogP,by="Gorai",all.x=TRUE))
gohir=c(cw$Goria.A,cw$Gohir.D)
# project JJ Gorai to Gohir
length(gohir <- gohir[!is.na(gohir)]) #544
# how many in rd
rd$JJ = rd$ID %in% gohir

load("funcEnrich.rdata")
