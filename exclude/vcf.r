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




##################################
## check Coding seq SNP pattern ##
##################################

# GATK generated SNP and INDEL between Maxxa and TX2094
# /work/LAS/jfw-lab/djyuan/Maxxa_TX2094/SNP/snpEff/AtDt26.PolyAll.snpEff.vcf
######## Unix command ###########
## snps
# cp /work/LAS/jfw-lab/djyuan/Maxxa_TX2094/SNP/snpEff/AtDt26.PolyAll.snpEff.vcf AtDt26.PolyAll.snpEff.vcf

# grep -v '^$\|^\s*#' AtDt26.PolyAll.snpEff.vcf | cut -f8 |cut -d "|" -f2 |sort|uniq -c
#         47260 3_prime_UTR_variant
#          3638 5_prime_UTR_premature_start_codon_gain_variant
#         23037 5_prime_UTR_variant
#        379244 downstream_gene_variant
#            28 initiator_codon_variant
#             1 initiator_codon_variant&splice_region_variant
#       8670127 intergenic_region
#        225951 intron_variant
#         91943 missense_variant
#          1438 missense_variant&splice_region_variant
#           555 splice_acceptor_variant&intron_variant
#             5 splice_acceptor_variant&splice_donor_variant&intron_variant
#             1 splice_acceptor_variant&splice_region_variant&intron_variant
#           404 splice_donor_variant&intron_variant
#             2 splice_donor_variant&splice_region_variant&intron_variant
#           886 splice_region_variant
#          6970 splice_region_variant&intron_variant
#            25 splice_region_variant&stop_retained_variant
#           886 splice_region_variant&synonymous_variant
#           168 start_lost
#             3 start_lost&splice_region_variant
#          1847 stop_gained
#            37 stop_gained&splice_region_variant
#           211 stop_lost
#            80 stop_lost&splice_region_variant
#           103 stop_retained_variant
#         53989 synonymous_variant
#        533389 upstream_gene_variant


# grep -v '^$\|^\s*#' AtDt26.PolyAll.snpEff.vcf | cut -f10,11,12 | sed 's/:\S*//g' |sort|uniq -c >pattern.txt
# sed -i 's/^ *//g' pattern.txt
# sed -i 's/ /\t/g' pattern.txt

---book

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

