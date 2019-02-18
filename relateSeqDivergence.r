options(scipen=999)
load("cistrans.rdata")

########################################
## get transcript sequence divergence ##
########################################

# import transcript SNPs
x<-read.table("/lss/research/jfw-lab/Projects/AD1_domestication_cistrans/hylite_partition/results112717/results112717.snp.txt",sep="\t",header=TRUE)
diagnostic=which((x$Maxxa ==0 & x$TX2094==1 & (x$F1=="1,0" | x$F1==1) ) | (x$Maxxa ==1 & x$TX2094==0 & (x$F1=="1,0" | x$F1==1)))
snps = x[diagnostic,]

# get ref transcript sequence
library("Biostrings")
ref <- readDNAStringSet("/lss/research/jfw-lab/Projects/AD1_domestication_cistrans/mapping_reference/TM1.nuclNmt.transcript.fasta")
names(ref) = gsub(" .*","", names(ref))
# get ref AA sequence

# get alt sequence
library(svMisc)
## maxxa
IDs = as.character(unique(snps$GENE[snps$Maxxa==1]))
nucM = ref
for(i in 1:length(IDs)){
    id=IDs[i]
    snp=snps[snps$GENE==id & snps$Maxxa==1,]
    nucM[[id]]=replaceLetterAt(nucM[[id]],snp$POS,as.character(snp$ALT))
    progress(i,length(IDs))
    if (i == length(IDs)) cat("Done!\n")
}
## TX2094
IDs = as.character(unique(snps$GENE[snps$TX2094==1]))
nucT = ref
for(i in 1:length(IDs)){
    id=IDs[i]
    snp=snps[snps$GENE==id & snps$TX2094==1,]
    nucT[[id]]=replaceLetterAt(nucT[[id]],snp$POS,as.character(snp$ALT))
    progress(i,length(IDs))
    if (i == length(IDs)) cat("Done!\n")
}

# get amino acid changes
aaRef <- Biostrings::translate(ref)
aaM <- Biostrings::translate(nucM)
aaT <- Biostrings::translate(nucT)
snps$aaPOS = ceiling(snps$POS / 3 )
snps$aaRef=NA
snps$aaM = NA
snps$aaT = NA
#library(svMisc)
for(i in 1:nrow(snps)){
    snps[i,"aaRef"] = substring( aaRef[snps[i,"GENE"] ] , snps[i,"aaPOS"] ,  snps[i,"aaPOS"])
    snps[i,"aaM"] = substring( aaM[snps[i,"GENE"] ] , snps[i,"aaPOS"] ,  snps[i,"aaPOS"])
    snps[i,"aaT"] = substring( aaT[snps[i,"GENE"] ] , snps[i,"aaPOS"] ,  snps[i,"aaPOS"])
    progress(i,nrow(snps))
    if (i == nrow(snps)) cat("Done!\n")
}
snps$MvsT= ifelse(snps$aaM==snps$aaT, "syn", "nonsyn")
table(snps$MvsT)
# nonsyn    syn
# 60612  30405
length(unique(snps$GENE))
# No. of genes contain diagnostic SNPs 27820

## KaKs calculatetion - seqinR
library(seqinr) #kaks
library(ape) #convert alignment format between biostrings to seqr
# kaks(ape::as.alignment(list(s1="ATC",s2="ATG")))
# genes with M vs T divergence
IDs = as.character(unique(snps$GENE))
# 1st gene
id=IDs[1]
aln=ape::as.alignment(list(r=as.character(nucM[[id]]),a=as.character(nucT[[id]])))
resKaKs = unlist(kaks(aln) )
kaks(aln)
for(i in 2:length(IDs)){
    id=IDs[i]
    aln=ape::as.alignment(list(r=nucM[[id]],a=nucT[[id]])) # DNAStringSet -> lists of single-character strings -> seqinr alignment
    # kaks(as.alignment(list(s1="ATC",s2="ATG")))
    res=kaks(aln)
    resKaKs = rbind(resKaKs, unlist(res))
    progress(i,length(IDs))
    if (i == length(IDs)) cat("Done!\n")
}
resKaKs = as.data.frame(resKaKs)
rownames(resKaKs) = IDs
# done and plot
pdf("plotKaKs.pdf")
plot(resKaKs[,1:2],pch=".")
boxplot(resKaKs[,1:2])
# vioplot(resKaKs[,1:2])
s=density(resKaKs$ks[resKaKs$ks>0])
ks.lg = paste0("Ks: peak at ",round(s$x[which.max(s$y)],4)) # Ks peak ar 0.0020879
a=density(resKaKs$ka[resKaKs$ka>0])
ka.lg = paste0("Ka: peak at ",round(a$x[which.max(a$y)],4)) # Ka peak ar 0.0011
a$x[which.max(a$y)] # Ka peak ar 0.001078022
plot(a,col="blue",xlim=c(0,0.02),ylim=c(0,max(max(s$y),max(a$y))+5))
lines(s,col="brown")
legend(x=0.01, y=max(max(s$y),max(a$y)), legend=c(ka.lg,ks.lg),col=c("blue","brown"),lty=1)
dev.off()

save(ref,nucM,nucT,resKaKs,snps, file="seqKaKs.rdata")


################
## PAML Ka/Ks ##
################
library(ape)
library(seqinr)
library(svMisc)
writePHYLIP=function(aln,file="",rm.stop=TRUE)
{
    seqLen= nchar(aln$seq[[1]])
    cat(paste0(" ",aln$nb," ",seqLen),file=file,sep="\n")
    for(j in 1:aln$nb){
        cat(aln$nam[j],file=file,sep="\n",append=TRUE)
        cat(substring(aln$seq[j],1,seqLen),file=file,sep="\n",append=TRUE)
    }
}
removeStop=function(aln){
    # stop codon position
    pos= sort(unique(unlist(lapply(aln$seq,function(x)which(translate(s2c(x[[1]])) =='*')))))
    # nuc pos of stop codon
    npos = as.numeric(sapply(pos, function(x)(x*3-2):(x*3)))
    # Anuc pos without stop codon
    len = nchar(aln$seq[[1]])
    ind=setdiff(1:len,npos)
    # is stop codon at the end?
    aln$seq = sapply(aln$seq, function(x)c2s(s2c(x)[ind]))
    names(aln$seq)=NULL
    print(pos)
    return(aln)
}

lastline <- function(filename) {
    ## filename is of mode character
    out <- system(sprintf("wc -l %s",filename),intern=TRUE)
    n <- as.integer(sub(sprintf("[ ]*([0-9]+)[ ]%s",filename),"\\1",out))
    # print(n)
    scan(filename,what="",skip=n-1,nlines=1,sep="\n",quiet=TRUE)
}
# module load paml
setwd("paml/")
load("seqKaKs.rdata")
aaRef <- Biostrings::translate(ref)
aaM <- Biostrings::translate(nucM)
aaT <- Biostrings::translate(nucT)
IDs = as.character(unique(snps$GENE))
df=data.frame(id=IDs,out=NA)
for(i in 1:length(IDs)){
    
    print(i)
    id=IDs[i]
    # prep alignment file
    aln=ape::as.alignment(list(maxxa=nucM[[id]],tx2094=nucT[[id]])) # DNAStringSet -> lists of single-character strings -> seqinr alignment
    aln=removeStop(aln)
    writePHYLIP(aln, file="aln.txt",rm.stop=TRUE) # remove stop codon
    #file.show("aln.txt")
    
    # prep tree file
    cat(paste0("(",paste(aln$nam,collapse=","),")"),file="treefile.txt",sep="\n")
    #file.show("treefile.txt")
    
    # prep control file
    # file.show("codeml.ctl.txt")
    # runmode=-2 to indicate that you are doing pairwise dNdS.
    
    # run codeml
    system("codeml codeml.ctl.txt")
    
    # report, silence
    o=lastline("results.txt")
    df[i,"out"]=gsub(" =","=",gsub("= ","=",gsub('\\s+'," ",o,perl=TRUE)))
    
    # save result
    # system(paste0("mv results.txt output/",id,"_report.txt"))
}
head(df)
header=
res= as.data.frame(t(data.frame(strsplit(df$out,split=" "))))
rownames(res)=NULL
names(res)= gsub("=.*","",unlist(strsplit(df$out[1],split=" ")))
res= apply(res,2,function(x){as.numeric(gsub(".*=","",x))})
paml=cbind(df,res)
save(paml,file="paml.rdata")

# pretty good cor between  geninr and paml
unique(rownames(resKaKs)==paml$id)
cor.test(resKaKs$ka,paml$dN)
cor.test(resKaKs$ks,paml$dS)
table(paml$"dN/dS">1)
# FALSE  TRUE
# 15391 12429
table(resKaKs$ka>resKaKs$ks)
# FALSE  TRUE
# 13901 13919


# plot
library(gplots)
pdf("plotKaKs_PAML.pdf")
textplot(summary(paml[,-c(1,2)]),cex=0.5)
hist(paml$"dN/dS")
plot(paml[,c("dS","dN")],pch=".")
boxplot(paml[,c("dS","dN")])
# vioplot(resKaKs[,1:2])
s=density(paml$dS[paml$dS>0])
ks.lg = paste0("Ks: peak at ",round(s$x[which.max(s$y)],4)) # Ks peak ar 0.0020879
a=density(paml$dN[paml$dN>0])
ka.lg = paste0("Ka: peak at ",round(a$x[which.max(a$y)],4)) # Ka peak ar 0.0011
a$x[which.max(a$y)] # Ka peak ar 0.001078022
plot(a,col="blue",xlim=c(0,0.02),ylim=c(0,max(max(s$y),max(a$y))+5))
lines(s,col="brown")
legend(x=0.01, y=max(max(s$y),max(a$y)), legend=c(ka.lg,ks.lg),col=c("blue","brown"),lty=1)
dev.off()




################################
## check promoter SNP pattern ##
################################

# GATK generated SNP and INDEL between Maxxa and TX2094 over 2kb promoter region
# HyLiTE generated SNPs of transcripts
# Sequence divergence was calculated as snps/site.
######## Unix command ###########
## snps
#cd /lss/research/jfw-lab/Projects/AD1_domestication_cistrans/promoterSNP_djyuan/SNP
#grep -v '^$\|^\s*#' *vcf | cut -f10,11,12 | sed 's/:\S*//g' |sort|uniq -c >pattern.txt
#sed -i 's/^ *//g' pattern.txt
#sed -i 's/ /\t/g' pattern.txt
#grep 'upstream' *vcf | cut -f8,10,11,12 | sed 's/:\S*//g' | sed 's/.*R[|]//g' >promoterSNPs.txt
## indels
#cd /lss/research/jfw-lab/Projects/AD1_domestication_cistrans/promoterSNP_djyuan/InDel
#grep -v '^$\|^\s*#' *snpEff.vcf | cut -f10,11,12 | sed 's/:\S*//g' |sort|uniq -c >pattern.txt
#sed -i 's/^ *//g' pattern.txt
#sed -i 's/ /\t/g' pattern.txt
#grep 'upstream' *snpEff.vcf | cut -f4,5,8,10,11,12 | sed 's/:\S*//g' | sed 's/\S*upstream_gene_variant[|]MODIFIER[|]//g' |sed 's/[|]\S*//g' >promoterIndels.txt

x<-read.table("/lss/research/jfw-lab/Projects/AD1_domestication_cistrans/promoterSNP_djyuan/SNP/pattern.txt",sep="\t")
t = sum(x$V1) # 566,329

# pairwise comparison
sum(x$V1[x$V2!=x$V3])/t  # maxxa != SAMN 75.1%
sum(x$V1[x$V2!=x$V4])/t  # maxxa != TX2094 72.6%
sum(x$V1[x$V3!=x$V4])/t  # SAMN != TX2094 52.1%

# consider only homozygotes
sum(x$V1[(x$V3=="0/0"&x$V2=="1/1")|(x$V2=="0/0"&x$V3=="1/1")])/t  # maxxa != SAMN 61.8%
sum(x$V1[(x$V4=="0/0"&x$V2=="1/1")|(x$V2=="0/0"&x$V4=="1/1")])/t  # maxxa != TX2094 31.0%
sum(x$V1[(x$V4=="0/0"&x$V3=="1/1")|(x$V3=="0/0"&x$V4=="1/1")])/t  # SAMN != TX2094 12.5%

# check homozygosity
sum(x$V1[(x$V2=="0/0"|x$V2=="1/1")])/t  # maxxa 84.2%
sum(x$V1[(x$V3=="0/0"|x$V3=="1/1")])/t  # SAMN 81.9%, by missisipi state
sum(x$V1[(x$V4=="0/0"|x$V4=="1/1")])/t  # TX2094 51.6%, by Udall


##################################
## check promoter Indel pattern ##
##################################

x<-read.table("/lss/research/jfw-lab/Projects/AD1_domestication_cistrans/promoterSNP_djyuan/InDel/pattern.txt",sep="\t")
t = sum(x$V1) # 143620

# pairwise comparison
sum(x$V1[x$V2!=x$V3])/t  # maxxa != SAMN 76.2%
sum(x$V1[x$V2!=x$V4])/t  # maxxa != TX2094 73.0%
sum(x$V1[x$V3!=x$V4])/t  # SAMN != TX2094 52.5%

# consider only homozygotes
sum(x$V1[(x$V3=="0/0"&x$V2=="1/1")|(x$V2=="0/0"&x$V3=="1/1")])/t  # maxxa != SAMN 61.4%
sum(x$V1[(x$V4=="0/0"&x$V2=="1/1")|(x$V2=="0/0"&x$V4=="1/1")])/t  # maxxa != TX2094 33.1%
sum(x$V1[(x$V4=="0/0"&x$V3=="1/1")|(x$V3=="0/0"&x$V4=="1/1")])/t  # SAMN != TX2094 13.4%

# check homozygosity
sum(x$V1[(x$V2=="0/0"|x$V2=="1/1")])/t  # maxxa 87.2%
sum(x$V1[(x$V3=="0/0"|x$V3=="1/1")])/t  # SAMN 83.5%
sum(x$V1[(x$V4=="0/0"|x$V4=="1/1")])/t  # TX2094 56.6%


###########################################
## get promoter variant coounts per gene ##
###########################################

#### SNP
x<-read.table("/lss/research/jfw-lab/Projects/AD1_domestication_cistrans/promoterSNP_djyuan/SNP/promoterSNPs.txt",sep="\t")
x$V1 = gsub("[|].*","",x$V1)
dim(x) # 550532 SNPs total
# SAMN
length(which((x$V3=="0/0"&x$V2=="1/1")|(x$V2=="0/0"&x$V3=="1/1"))) #340642 diagnostic snps
length(unique(x$V1[(x$V3=="0/0"&x$V2=="1/1")|(x$V2=="0/0"&x$V3=="1/1")])) # found in 53252
s = as.data.frame(table(x$V1[(x$V3=="0/0"&x$V2=="1/1")|(x$V2=="0/0"&x$V3=="1/1")]))
dim(s) #53252 genes
# TX2094
length(which((x$V4=="0/0"&x$V2=="1/1")|(x$V2=="0/0"&x$V4=="1/1"))) #171196 diagnostic snps
length(unique(x$V1[(x$V4=="0/0"&x$V2=="1/1")|(x$V2=="0/0"&x$V4=="1/1")])) # found in 33506 genes
t = as.data.frame(table(x$V1[(x$V4=="0/0"&x$V2=="1/1")|(x$V2=="0/0"&x$V4=="1/1")]))
dim(t) #33506 genes

#### INDEL
x<-read.table("/lss/research/jfw-lab/Projects/AD1_domestication_cistrans/promoterSNP_djyuan/InDel/promoterIndels.txt",sep="\t")
x$size = abs(nchar(as.character(x$V1))-nchar(as.character(x$V2)))
dim(x) # 143482 indels
# SAMN
length(which((x$V4=="0/0"&x$V5=="1/1")|(x$V5=="0/0"&x$V4=="1/1"))) # 88137 indels
length(unique(x$V3[(x$V4=="0/0"&x$V5=="1/1")|(x$V5=="0/0"&x$V4=="1/1")])) # found in 40218 genes
s.indel = as.data.frame(table(x$V3[(x$V4=="0/0"&x$V5=="1/1")|(x$V5=="0/0"&x$V4=="1/1")]))
s.indel.size = with(x[(x$V4=="0/0"&x$V5=="1/1")|(x$V5=="0/0"&x$V4=="1/1"),],aggregate(size,by=list(V3),function(x)sum(x)))
dim(s.indel) # 48284     2
dim(s.indel.size) # 40218     2# TX2094
# TX2094
length(which((x$V4=="0/0"&x$V6=="1/1")|(x$V6=="0/0"&x$V4=="1/1"))) # 47561indels
length(unique(x$V3[(x$V4=="0/0"&x$V6=="1/1")|(x$V6=="0/0"&x$V4=="1/1")])) # found in 23683 genes
t.indel = as.data.frame(table(x$V3[(x$V4=="0/0"&x$V6=="1/1")|(x$V6=="0/0"&x$V4=="1/1")]))
t.indel.size = with(x[(x$V4=="0/0"&x$V6=="1/1")|(x$V6=="0/0"&x$V4=="1/1"),],aggregate(size,by=list(V3),function(x)sum(x)))
dim(t.indel) # 48284     2
dim(t.indel.size) # 23683


#########################
## Sequence divergence ##
#########################

len = read.table("Genes66610.gene.length.txt",header=TRUE,sep="\t")
# coding seq
load("seqKaKs.rdata")
rownames(paml)=paml$id
exon = as.data.frame(table(snps$GENE))
syn = as.data.frame(table(snps$GENE[snps$MvsT=="syn"]))
nonsyn = as.data.frame(table(snps$GENE[snps$MvsT=="nonsyn"]))
seqDev =data.frame(ID=rownames(len), len=len$x, exon =0, syn=0, nonsyn=0, s = 0, t = 0, s.indel = 0, t.indel = 0, s.indel.size = 0, t.indel.size = 0)
seqDev$exon[match(exon$Var1,rownames(len))] = exon$Freq
seqDev$syn[match(syn$Var1,rownames(len))] = syn$Freq
seqDev$nonsyn[match(nonsyn$Var1,rownames(len))] = nonsyn$Freq
seqDev$s[match(paste0(s$Var1,".1"),rownames(len))] = s$Freq
seqDev$t[match(paste0(t$Var1,".1"),rownames(len))] = t$Freq
seqDev$s.indel[match(paste0(s.indel$Var1,".1"),rownames(len))] = s.indel$Freq
seqDev$t.indel[match(paste0(t.indel$Var1,".1"),rownames(len))] = t.indel$Freq
seqDev$s.indel.size[match(paste0(s.indel.size$Group.1,".1"),rownames(len))] = s.indel.size$x
seqDev$t.indel.size[match(paste0(t.indel.size$Group.1,".1"),rownames(len))] = t.indel.size$x
colSums(seqDev[,-1])
#    len         exon           syn        nonsyn            s            t
# 81346143        91017        30405        60612       340642       171196
# s.indel      t.indel s.indel.size t.indel.size
# 88137        47561       274420       145999

# number of SNPs per KB
perKB = data.frame(ID=rownames(len), promoter = seqDev$s/2, coding = seqDev$exon/seqDev$len*1000, syn=seqDev$syn/seqDev$len*1000, nonyn=seqDev$nonsyn/seqDev$len*1000)
save(s, s.indel, s.indel.size, t, t.indel, t.indel.size, seqDev,perKB, file="seqDivergence.rdata")

####################################################
## relate to expression and regulatory divergence ##
####################################################
load("cistrans.rdata")
load("RDgenes.rdata")
load("seqDivergence.rdata")
load("paml/paml.rdata")
panel.cor <- function(x, y, digits=2, prefix="", cex.cor, ...)
{
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    r <- abs(cor(x, y, use="complete.obs"))
    txt <- format(c(r, 0.123456789), digits=digits)[1]
    txt <- paste(prefix, txt, sep="")
    if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
    text(0.5, 0.5, txt, cex = cex.cor * r)
}

## restrict to 27816 genes
S <- perKB[,-1]
rownames(S)<-perKB$ID
S <- S[rownames(res.mt10),]
summary(S)
#    promoter          coding              syn              nonyn
# Min.   : 0.000   Min.   : 0.08708   Min.   : 0.0000   Min.   : 0.0000
# 1st Qu.: 1.000   1st Qu.: 1.06667   1st Qu.: 0.0000   1st Qu.: 0.6173
# Median : 2.500   Median : 1.92234   Median : 0.4845   Median : 1.2820
# Mean   : 2.919   Mean   : 2.89832   Mean   : 0.8814   Mean   : 2.0169
# 3rd Qu.: 4.000   3rd Qu.: 3.44353   3rd Qu.: 1.1779   3rd Qu.: 2.4783
# Max.   :56.000   Max.   :55.55556   Max.   :41.0256   Max.   :47.2222
S <- seqDev[,-1]
rownames(S)<-seqDev$ID
S <- S[rownames(res.mt10),]
cor(S)
#                      len        exon           Na          Ns            s
# len           1.000000000 -0.31055774 -0.313696976 -0.14138233 -0.001791317
# exon         -0.310557742  1.00000000  0.902283099  0.64612941  0.060941295
# Na           -0.313696976  0.90228310  1.000000000  0.25393047  0.028123919
# Ns           -0.141382326  0.64612941  0.253930473  1.00000000  0.086928787
# s            -0.001791317  0.06094130  0.028123919  0.08692879  1.000000000
# t             0.001266761  0.04462796  0.023176938  0.05908902  0.568157609
# s.indel       0.012226490  0.02351427  0.003220064  0.04705131  0.480819723
# t.indel       0.003421127  0.02563199  0.012659385  0.03509232  0.300853688
# s.indel.size -0.003946621  0.02257021  0.006516198  0.03909847  0.284619820
# t.indel.size -0.009138237  0.02656799  0.017084996  0.02935772  0.188698245
#                        t     s.indel     t.indel s.indel.size t.indel.size
# len          0.001266761 0.012226490 0.003421127 -0.003946621 -0.009138237
# exon         0.044627958 0.023514274 0.025631992  0.022570213  0.026567985
# Na           0.023176938 0.003220064 0.012659385  0.006516198  0.017084996
# Ns           0.059089017 0.047051306 0.035092315  0.039098470  0.029357721
# s            0.568157609 0.480819723 0.300853688  0.284619820  0.188698245
# t            1.000000000 0.297049981 0.624170559  0.175654643  0.383195528
# s.indel      0.297049981 1.000000000 0.598521888  0.528051211  0.343329059
# t.indel      0.624170559 0.598521888 1.000000000  0.316798781  0.575797031
# s.indel.size 0.175654643 0.528051211 0.316798781  1.000000000  0.607360272
# t.indel.size 0.383195528 0.343329059 0.575797031  0.607360272  1.000000000
### correlation seen between SAMN and TX2094; no correlation seen between promoter and transcript
rownames(paml)=paml$id
R=paml[rownames(res.mt10),]
S$dNdS=R$"dN/dS"
S$dS=R$"dS"
S$dN=R$"dN"

pdf("seqDivergence.pdf")
textplot(cor(S))
mtext("correlation seen between SAMN and TX2094; no correlation seen between promoter and transcript")

## sequence divergence NOT correlate with expression and regulatory divergence
#par(mfrow=c(3,2))
sumR=c("dataset","expression","seq.divergence","cor")
for(l in c("res10","res20","res.mt10","res.mt20","res.tm10","res.tm20")){
    res= get(l)
    for(d in c("A","B")){
        for(c in 2:10){
            # plot(S[,c],abs(res[,d]),xlab=names(S)[c])
            r=cor(S[,c],abs(res[,d]),use="pairwise.complete.obs")
            sumR=rbind(sumR,(c(l,d,names(S)[c],r)))
        }
    }
}
# plot results
quantile(as.numeric(sumR[-1,4]))
# 0%          25%          50%          75%         100%
# -0.044273969  0.005639656  0.015857305  0.035966835  0.102126001
textplot(sumR)
mtext("sequence divergence NOT correlate with expression and regulatory divergence")


## promoter SNPs
c=S$s
quantile(c,(1:10)/10)
# 10%  20%  30%  40%  50%  60%  70%  80%  90% 100%
# 0    1    3    4    5    6    8    9   12  112  # 90% less than 13
c[c>12]=13
table(c)
cName<-c(as.character(0:12),">12")
# |A| parental divergence
boxplot(abs(res.mt10$A)~c, axes = FALSE, xlab="Number of SNPs in promoter (2kb)", main="|A|");axis(2);axis(1,1:14, labels=cName);abline(h=median(abs(res.mt10$A),na.rm=TRUE),col="blue");text(2,10,paste0("median = ",round(median(abs(res.mt10$A),na.rm=TRUE),3)),col="blue")
# |B| cis divergence
boxplot(abs(res.mt10$B)~c, axes = FALSE, xlab="Number of SNPs in promoter (2kb)", main="|B|");axis(2);axis(1,1:14, labels=cName);abline(h=median(abs(res.mt10$B),na.rm=TRUE),col="blue");text(2,10,paste0("median = ",round(median(abs(res.mt10$B),na.rm=TRUE),3)),col="blue")
# |B|/(|AminusB|+|B|)
boxplot((abs(res.mt10$B)/(abs(res.mt10$B)+abs(res.mt10$AminusB)))~c, axes = FALSE, xlab="Number of SNPs in promoter (2kb)", main="cis%");axis(2);axis(1,1:14, labels=cName);abline(h=0.5,col="red")

##
S$rd = NA
S[rownames(rd),"rd"] = rd$rd
library(agricolae)
out=duncan.test(aov(s~rd,data=S),"rd", console=TRUE)
plot(out,variation="SE", main="Promoter SNPs of RD genes")
out=duncan.test(aov(s.indel~rd,data=S),"rd", console=TRUE)
plot(out,variation="SE",main="Promoter Indels of RD genes")
out=duncan.test(aov(exon~rd,data=S),"rd", console=TRUE)
plot(out,variation="SE",main="Coding sequence SNPs (per kb) of RD genes")
out=duncan.test(aov(syn~rd,data=S),"rd",main="Snonymous SNPs (per kb) of RD genes", console=TRUE)
plot(out,variation="SE",main="Snonymous SNPs (per kb) of RD genes")
out=duncan.test(aov(nonsyn~rd,data=S),"rd", console=TRUE)
plot(out,variation="SE",main="Nonsynonymous SNPs (per kb) of RD genes")
S$rate=S$nonsyn/S$syn
summary(S$rate[is.finite(S$rate)])
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.000   0.000   1.000   1.387   2.000  26.000
S$rate[!is.finite(S$rate)]=30
out=duncan.test(aov(rate~rd,data=S),"rd", console=TRUE) #NULL, aov not s
plot(out,variation="SE",main="Na/Ns of RD genes")

out=duncan.test(aov(dNdS~rd,data=S),"rd", console=TRUE) #NULL, aov not s
plot(out,variation="SE",main="Ka/Ks of RD genes")
out=duncan.test(aov(dN~rd,data=S),"rd", console=TRUE) #NULL, aov not s
plot(out,variation="SE",main="Ka of RD genes")
out=duncan.test(aov(dS~rd,data=S),"rd", console=TRUE) #NULL, aov not s
plot(out,variation="SE",main="Ks of RD genes")
dev.off()
# relate seq divergence with seven category, see 'boxplotByCategory.r'

ks.test(S$dN[S$rd=="no divergence"], S$dN[S$rd!="no divergence"])
ks.test(S$dS[S$rd=="no divergence"], S$dS[S$rd!="no divergence"])
ks.test(S$dNdS[S$rd=="no divergence"], S$dNdS[S$rd!="no divergence"])
