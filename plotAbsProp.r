## Absolute proportion of |cis|/|cis| + |trans|
library(gplots)
library(ggplot2)
load("cistrans.rdata")

# cis% with expression change magnitude
with(res.mt10,cor.test(abs(A),abs(B)/(abs(AminusB)+abs(B)))) # -0.016
with(res.tm10,cor.test(abs(A),abs(B)/(abs(AminusB)+abs(B)))) # -0.044
with(res.mt20,cor.test(abs(A),abs(B)/(abs(AminusB)+abs(B)))) # -0.079
with(res.tm20,cor.test(abs(A),abs(B)/(abs(AminusB)+abs(B)))) # -0.076
# trans% with expression change magnitude
with(res.mt10,cor.test(abs(A),abs(AminusB)/(abs(AminusB)+abs(B)))) # 0.016
with(res.tm10,cor.test(abs(A),abs(AminusB)/(abs(AminusB)+abs(B)))) # 0.044
with(res.mt20,cor.test(abs(A),abs(AminusB)/(abs(AminusB)+abs(B)))) # 0.079
with(res.tm20,cor.test(abs(A),abs(AminusB)/(abs(AminusB)+abs(B)))) # 0.076


## Previous plot binned by A value 0,1,3,4,>4
pdf("plotCistrans.absProp.binN.I-V.pdf")
l<-c("res10","res20","res.mt10","res.mt20","res.tm10","res.tm20")
row<-c("dataset","group","prop","se","ci","geneN")
# |cis| vs |trans|
for(i in l){
    x<-get(i)
    x=x[grepl("^[1-5]",x$category),]
    print(i)
    #print(wilcox.test(abs(x$B),abs(x$AminusB)))
    #print(t.test(abs(x$B),abs(x$AminusB)))
    print(cor.test(abs(x$A),abs(x$B)/(abs(x$B)+abs(x$AminusB))))
}

    for(i in l){
    x<-get(i)
    x=x[grepl("^[1-5]",x$category),]
    p<-abs(x$B)/(abs(x$AminusB)+abs(x$B))
    # bin genes by A to 6 classes
    breaks <- c(0, 1, 2, 3, 4, 100)
    c<-.bincode(x=abs(x$A), b=breaks, TRUE)
    cName<-c("(0-1)","(1-2)","(2-3)","(3-4)","(4+)")
    # boxplot(p~c, main=i) hard to see patterns
    for(j in 1:length(cName))
    {
        pj <- p[c==j]
        mean<- mean(pj,na.rm=TRUE)
        # get error bar for 95$ CI
        error <- qt(0.975,df=length(pj)-1)*sd(pj,na.rm=TRUE)/sqrt(length(pj))
        # standard error of the mean
        se<-sd(pj,na.rm=TRUE)/sqrt(length(pj))
        row<-rbind(row,c(i,cName[j],mean,se,error,length(pj)))
    }
    mean<- mean(p,na.rm=TRUE)
    # get error bar for 95$ CI
    error <- qt(0.975,df=length(p)-1)*sd(p,na.rm=TRUE)/sqrt(length(p))
    # standard error of the mean
    se<-sd(p,na.rm=TRUE)/sqrt(length(p))
    row<-rbind(row,c(i,"all",mean,se,error,length(p)))
}
rownames(row)=NULL
colnames(row)=row[1,]
propT <-as.data.frame(row[-1,])
propT[,3:6]<-apply(propT[,3:6],2,as.numeric)
# plot results
textplot(propT)
mtext("Proportion of divergence due to cis regulatory effect grouped by overall parental divergence calculated")
# Basic line plot with points
df<-propT[propT$group!="all",]
pd <- position_dodge(0.1)
p <- ggplot(data=df, aes(x=group, y=prop, group=dataset, color=dataset)) +
geom_errorbar(aes(ymin=prop-ci, ymax=prop+ci), width=.1, position=pd) +
geom_line(position=pd) +
geom_point(position=pd) + labs(title ="|cis|/(|cis|+|trans|)", x = "log2 Parental Divergence (cis + trans)", y = "Proportion divergence due to cis")
print(p + facet_grid(dataset ~ .))
dpa<-gsub("res|[.]|m|t","", df$dataset)
print(p + facet_grid(factor(dpa) ~ .))
# seperately plot F1 and mt, tm
print(p)
print(p %+% subset(df, dataset %in% c("res10","res20")))
print(p %+% subset(df, dataset %in% c("res.mt10","res.mt20","res.tm10","res.tm20")))
dev.off()

## for figure
## Previous plot binned by A value 0,1,3,4,>4
l<-c("res.mt10","res.mt20","res.tm10","res.tm20")
row<-c("dataset","group","prop","se","ci","geneN")
# |cis| vs |trans|
for(i in l){
    x<-get(i)
    print(i)
    #print(wilcox.test(abs(x$B),abs(x$AminusB)))
    #print(t.test(abs(x$B),abs(x$AminusB)))
    print(cor.test(abs(x$A),abs(x$B)/(abs(x$B)+abs(x$AminusB))))
} # -0.016, -0.079, -0.044, -0.076
for(i in l){
    x<-get(i)
    # x=x[grepl("^[1-5]",x$category),]
    p<-abs(x$B)/(abs(x$AminusB)+abs(x$B))
    # bin genes by A to 6 classes
    breaks <- c(0, 1, 2, 3, 4, 100)
    c<-.bincode(x=abs(x$A), b=breaks, TRUE)
    cName<-c("(0-1)","(1-2)","(2-3)","(3-4)","(4+)")
    # boxplot(p~c, main=i) hard to see patterns
    for(j in 1:length(cName))
    {
        pj <- p[c==j]
        mean<- mean(pj,na.rm=TRUE)
        # get error bar for 95$ CI
        error <- qt(0.975,df=length(pj)-1)*sd(pj,na.rm=TRUE)/sqrt(length(pj))
        # standard error of the mean
        se<-sd(pj,na.rm=TRUE)/sqrt(length(pj))
        row<-rbind(row,c(i,cName[j],mean,se,error,length(pj)))
    }
    mean<- mean(p,na.rm=TRUE)
    # get error bar for 95$ CI
    error <- qt(0.975,df=length(p)-1)*sd(p,na.rm=TRUE)/sqrt(length(p))
    # standard error of the mean
    se<-sd(p,na.rm=TRUE)/sqrt(length(p))
    row<-rbind(row,c(i,"all",mean,se,error,length(p)))
}
rownames(row)=NULL
colnames(row)=row[1,]
propT <-as.data.frame(row[-1,])
propT[,3:6]<-apply(propT[,3:6],2,as.numeric)
# plot results
pdf("plotCistrans.absProp.binN.pdf")
textplot(propT)
mtext("Proportion of divergence due to cis regulatory effect grouped by overall parental divergence calculated")
# Basic line plot with points
df<-propT[propT$group!="all",]
pd <- position_dodge(0.1)
p <- ggplot(data=df, aes(x=group, y=prop, group=dataset, color=dataset)) +
geom_errorbar(aes(ymin=prop-ci, ymax=prop+ci), width=.1, position=pd) +
geom_line(position=pd) +
geom_point(position=pd) + labs(title ="|cis|/(|cis|+|trans|)", x = "log2 Parental Divergence (cis + trans)", y = "Proportion divergence due to cis")+theme_bw()
print(p + facet_grid(dataset ~ .))
dpa<-gsub("res|[.]|m|t","", df$dataset)
print(p + facet_grid(factor(dpa) ~ .))
# seperately plot F1 and mt, tm
print(p)
#print(p %+% subset(df, dataset %in% c("res10","res20")))
#print(p %+% subset(df, dataset %in% c("res.mt10","res.mt20","res.tm10","res.tm20")))
dev.off()

pdf("plotCistrans.absProp.binN.At.pdf")
## Previous plot binned by A value 0,1,3,4,>4
l<-c("res.mt10","res.mt20","res.tm10","res.tm20")
row<-c("dataset","group","prop","se","ci","geneN")
# |cis| vs |trans|
for(i in l){
    x<-get(i)
    x<-x[grep("Gohir.A",rownames(x)),]
    print(i)
    #print(wilcox.test(abs(x$B),abs(x$AminusB)))
    #print(t.test(abs(x$B),abs(x$AminusB)))
    print(cor.test(abs(x$A),abs(x$B)/(abs(x$B)+abs(x$AminusB))))
} # -0.027, -0.091, -0.056, -0.078
for(i in l){
    x<-get(i)
    x<-x[grep("Gohir.A",rownames(x)),]
    # x=x[grepl("^[1-5]",x$category),]
    p<-abs(x$B)/(abs(x$AminusB)+abs(x$B))
    # bin genes by A to 6 classes
    breaks <- c(0, 1, 2, 3, 4, 100)
    c<-.bincode(x=abs(x$A), b=breaks, TRUE)
    cName<-c("(0-1)","(1-2)","(2-3)","(3-4)","(4+)")
    # boxplot(p~c, main=i) hard to see patterns
    for(j in 1:length(cName))
    {
        pj <- p[c==j]
        mean<- mean(pj,na.rm=TRUE)
        # get error bar for 95$ CI
        error <- qt(0.975,df=length(pj)-1)*sd(pj,na.rm=TRUE)/sqrt(length(pj))
        # standard error of the mean
        se<-sd(pj,na.rm=TRUE)/sqrt(length(pj))
        row<-rbind(row,c(i,cName[j],mean,se,error,length(pj)))
    }
    mean<- mean(p,na.rm=TRUE)
    # get error bar for 95$ CI
    error <- qt(0.975,df=length(p)-1)*sd(p,na.rm=TRUE)/sqrt(length(p))
    # standard error of the mean
    se<-sd(p,na.rm=TRUE)/sqrt(length(p))
    row<-rbind(row,c(i,"all",mean,se,error,length(p)))
}
rownames(row)=NULL
colnames(row)=row[1,]
propT <-as.data.frame(row[-1,])
propT[,3:6]<-apply(propT[,3:6],2,as.numeric)
# plot results
textplot(propT)
mtext("Proportion of divergence due to cis regulatory effect grouped by overall parental divergence calculated")
# Basic line plot with points
df<-propT[propT$group!="all",]
pd <- position_dodge(0.1)
p <- ggplot(data=df, aes(x=group, y=prop, group=dataset, color=dataset)) +
geom_errorbar(aes(ymin=prop-ci, ymax=prop+ci), width=.1, position=pd) +
geom_line(position=pd) +
geom_point(position=pd) + labs(title ="|cis|/(|cis|+|trans|)", x = "log2 Parental Divergence (cis + trans)", y = "Proportion divergence due to cis")+theme_bw()
print(p + facet_grid(dataset ~ .))
dpa<-gsub("res|[.]|m|t","", df$dataset)
print(p + facet_grid(factor(dpa) ~ .))
# seperately plot F1 and mt, tm
print(p)
#print(p %+% subset(df, dataset %in% c("res10","res20")))
#print(p %+% subset(df, dataset %in% c("res.mt10","res.mt20","res.tm10","res.tm20")))
dev.off()

pdf("plotCistrans.absProp.binN.Dt.pdf")
## Previous plot binned by A value 0,1,3,4,>4
l<-c("res.mt10","res.mt20","res.tm10","res.tm20")
row<-c("dataset","group","prop","se","ci","geneN")
# |cis| vs |trans|
for(i in l){
    x<-get(i)
    x<-x[grep("Gohir.D",rownames(x)),]
    print(i)
    #print(wilcox.test(abs(x$B),abs(x$AminusB)))
    #print(t.test(abs(x$B),abs(x$AminusB)))
    print(cor.test(abs(x$A),abs(x$B)/(abs(x$B)+abs(x$AminusB))))
} # -0.0033, -0.063, -0.032, -0.074
for(i in l){
    x<-get(i)
    x<-x[grep("Gohir.D",rownames(x)),]
    # x=x[grepl("^[1-5]",x$category),]
    p<-abs(x$B)/(abs(x$AminusB)+abs(x$B))
    # bin genes by A to 6 classes
    breaks <- c(0, 1, 2, 3, 4, 100)
    c<-.bincode(x=abs(x$A), b=breaks, TRUE)
    cName<-c("(0-1)","(1-2)","(2-3)","(3-4)","(4+)")
    # boxplot(p~c, main=i) hard to see patterns
    for(j in 1:length(cName))
    {
        pj <- p[c==j]
        mean<- mean(pj,na.rm=TRUE)
        # get error bar for 95$ CI
        error <- qt(0.975,df=length(pj)-1)*sd(pj,na.rm=TRUE)/sqrt(length(pj))
        # standard error of the mean
        se<-sd(pj,na.rm=TRUE)/sqrt(length(pj))
        row<-rbind(row,c(i,cName[j],mean,se,error,length(pj)))
    }
    mean<- mean(p,na.rm=TRUE)
    # get error bar for 95$ CI
    error <- qt(0.975,df=length(p)-1)*sd(p,na.rm=TRUE)/sqrt(length(p))
    # standard error of the mean
    se<-sd(p,na.rm=TRUE)/sqrt(length(p))
    row<-rbind(row,c(i,"all",mean,se,error,length(p)))
}
rownames(row)=NULL
colnames(row)=row[1,]
propT <-as.data.frame(row[-1,])
propT[,3:6]<-apply(propT[,3:6],2,as.numeric)
# plot results
textplot(propT)
mtext("Proportion of divergence due to cis regulatory effect grouped by overall parental divergence calculated")
# Basic line plot with points
df<-propT[propT$group!="all",]
pd <- position_dodge(0.1)
p <- ggplot(data=df, aes(x=group, y=prop, group=dataset, color=dataset)) +
geom_errorbar(aes(ymin=prop-ci, ymax=prop+ci), width=.1, position=pd) +
geom_line(position=pd) +
geom_point(position=pd) + labs(title ="|cis|/(|cis|+|trans|)", x = "log2 Parental Divergence (cis + trans)", y = "Proportion divergence due to cis")+theme_bw()
print(p + facet_grid(dataset ~ .))
dpa<-gsub("res|[.]|m|t","", df$dataset)
print(p + facet_grid(factor(dpa) ~ .))
# seperately plot F1 and mt, tm
print(p)
#print(p %+% subset(df, dataset %in% c("res10","res20")))
#print(p %+% subset(df, dataset %in% c("res.mt10","res.mt20","res.tm10","res.tm20")))
dev.off()

#------------------------------------------------------------------------------------------
## binned by quantile, decreasing, not making sense
pdf("plotCistrans.absProp.binQ.pdf")
l<-c("res10","res20","res.mt10","res.mt20","res.tm10","res.tm20")
row<-c("dataset","group","prop","se","ci","geneN")
for(i in l){
    x<-get(i)
    p<-abs(x$B)/(abs(x$AminusB)+abs(x$B))
    # bin genes by A to 6 classes
    breaks <- quantile(abs(x$A),na.rm=TRUE)
    c<-.bincode(x=abs(x$A), b=breaks, TRUE)
    cName<-c("Q1","Q2","Q3","Q4")
    # boxplot(p~c, main=i) hard to see patterns
    for(j in 1:length(cName))
    {
        pj <- p[c==j]
        mean<- mean(pj,na.rm=TRUE)
        # get error bar for 95$ CI
        error <- qt(0.975,df=length(pj)-1)*sd(pj,na.rm=TRUE)/sqrt(length(pj))
        # standard error of the mean
        se<-sd(pj,na.rm=TRUE)/sqrt(length(pj))
        row<-rbind(row,c(i,cName[j],mean,se,error,length(pj)))
    }
    mean<- mean(p,na.rm=TRUE)
    # get error bar for 95$ CI
    error <- qt(0.975,df=length(p)-1)*sd(p,na.rm=TRUE)/sqrt(length(p))
    # standard error of the mean
    se<-sd(p,na.rm=TRUE)/sqrt(length(p))
    row<-rbind(row,c(i,"all",mean,se,error,length(p)))
}
rownames(row)=NULL
colnames(row)=row[1,]
propT <-as.data.frame(row[-1,])
propT[,3:6]<-apply(propT[,3:6],2,as.numeric)
# plot results
textplot(propT)
mtext("Proportion of divergence due to cis regulatory effect grouped by overall parental divergence calculated")
# Basic line plot with points
df<-propT[propT$group!="all",]
pd <- position_dodge(0.1)
p <- ggplot(data=df, aes(x=group, y=prop, group=dataset, color=dataset)) +
geom_errorbar(aes(ymin=prop-ci, ymax=prop+ci), width=.1, position=pd) +
geom_line(position=pd) +
geom_point(position=pd) + labs(title ="|cis|/(|cis|+|trans|)", x = "log2 Parental Divergence (cis + trans)", y = "Proportion divergence due to cis")
print(p + facet_grid(dataset ~ .))
dpa<-gsub("res|[.]|m|t","", df$dataset)
print(p + facet_grid(factor(dpa) ~ .))
# seperately plot F1 and mt, tm
print(p)
print(p %+% subset(df, dataset %in% c("res10","res20")))
print(p %+% subset(df, dataset %in% c("res.mt10","res.mt20","res.tm10","res.tm20")))
dev.off()
