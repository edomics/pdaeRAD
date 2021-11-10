###Code to reproduce analyses looking at the role of environment and distance in driving population structure.
###The approach follows that of Saenz-Agudelo 2015 (https://doi.org/10.1111/mec.13471) and uses the R code associated with the paper, modified for this dataset.

library(dplyr)
library(gridExtra)
setwd("/path/to/directory/")
#model data contains pairwise geographic and environmental distance (euclidian distance between pc1 & pc2 of environmental variables - min sal, max sal, min temp, max temp)
bstdata=read.table("model_data.txt",header=TRUE)
#pairwise fsts calculated in vcftools
fst=read.table("SNP.distLD.bi.fst.comb",header=FALSE)
names(fst)=c("site.comp","fst")
bstdata.fst=left_join(bstdata,fst, by = c("site"="site.comp"))

#define models
m1=lm(scale(fst)~scale(geod),data=bstdata.fst)
m2=lm(scale(fst)~scale(geod)+scale(envd),data=bstdata.fst)
m3=lm(scale(fst)~scale(geod)+scale(envd)+factor(barrier),data=bstdata.fst)
m4=lm(scale(fst)~scale(geod)+factor(barrier),data=bstdata.fst)
m5=lm(scale(fst)~scale(envd),data=bstdata.fst)
m6=lm(scale(fst)~scale(envd)+factor(barrier),data=bstdata.fst)
m7=lm(scale(fst)~factor(barrier),data=bstdata.fst)
Formula <- list(m1="geod",
                m2="geod + envd",
                m3="geod + envd + barrier",
                m4="geod + barrier",
                m5="envd",
                m6="envd + barrier",
                m7="barrier"
                )
Formula <- as.character(Formula)
RSS <- c(sum(resid(m1)^2),
         sum(resid(m2)^2),
         sum(resid(m3)^2),
         sum(resid(m4)^2),
         sum(resid(m5)^2),
         sum(resid(m6)^2),
         sum(resid(m7)^2)
         )
R2 <- c(summary(m1)$r.squared,
        summary(m2)$r.squared,
        summary(m3)$r.squared,
        summary(m4)$r.squared,
        summary(m5)$r.squared,
        summary(m6)$r.squared,
        summary(m7)$r.squared
        )
Adj.R2 <- round(c(summary(m1)$adj.r.squared,
                  summary(m2)$adj.r.squared,
                  summary(m3)$adj.r.squared,
                  summary(m4)$adj.r.squared,
                  summary(m5)$adj.r.squared,
                  summary(m6)$adj.r.squared,
                  summary(m7)$adj.r.squared
                  ),4)
AIC.table <- AIC(m1,m2,m3,m4,m5,m6,m7)
colnames(AIC.table) <- c("K","AIC")
K <- AIC.table$K
AICc <- AIC.table$AIC + 2*K*(K+1)/(nrow(bstdata.fst)-K-1)

delta.AICc <- rep(NA,nrow(AIC.table))
for(i in 1: nrow(AIC.table)){
  delta.AICc[i] <- AICc[i] - min(AICc) }
Model.likelihood  <- exp(-0.5 * delta.AICc) # of model i relative to the best K-L model
Model.probability <- Model.likelihood / sum(Model.likelihood) 
Evidence.ratio    <- rep(NA,nrow(AIC.table))
for(i in 1: nrow(AIC.table)){
  Evidence.ratio[i] <- max(Model.likelihood)/Model.likelihood[i] } # of the best K-L model relative to model i
AIC.table <- cbind(Formula, AIC.table, AICc, RSS, R2, Adj.R2, delta.AICc,Model.likelihood, Model.probability, Evidence.ratio)
AIC.table <- AIC.table[order(AIC.table$AICc),]
AIC.table[,3:ncol(AIC.table)] <- round(AIC.table[,3:ncol(AIC.table)], 4)
AIC.table
write.csv(AIC.table, file = "AIC.table_lm_bst_ALL_models.csv") # MAIN RESULTS

#Plotting figures for entire dataset
#with CIs
ggplot(bstdata.fst,aes(y=scale(fst),x=scale(envd)+scale(geod),colour=factor(barrier))) + geom_point() +scale_color_manual(values=c("orange","blue"))+stat_smooth(data=bstdata.fst,method="lm",se=T)+theme_bw()+guides(color=F) +ylab("Genetic Distance")+xlab("Geo + Env Distance")+theme(axis.text = element_text(size=14), axis.title = element_text(size=20))
#without CIs
ggplot(bstdata.fst,aes(y=scale(fst),x=scale(envd)+scale(geod),colour=factor(barrier))) + geom_point() +scale_color_manual(values=c("orange","blue"))+stat_smooth(data=bstdata.fst,method="lm",se=F)+theme_bw()+guides(color=F) +ylab("Genetic Distance")+xlab("Geo + Env Distance")+theme(axis.text = element_text(size=14), axis.title = element_text(size=20))

#final plots including for env/geo/bar only
a=ggplot(bstdata.fst,aes(y=scale(fst),x=scale(barrier))) + geom_point() +stat_smooth(data=bstdata.fst,method="lm",se=T)+theme_bw()+guides(color=F) +ylab("Genetic Distance")+xlab("Barrier")+theme(axis.text = element_text(size=14), axis.title = element_text(size=20)) + ylim(c(-2.1,2.25))
b=ggplot(bstdata.fst,aes(y=scale(fst),x=scale(geod))) + geom_point() +stat_smooth(data=bstdata.fst,method="lm",se=T)+theme_bw()+guides(color=F) +ylab("Genetic Distance")+xlab("Geo Distance")+theme(axis.text = element_text(size=14), axis.title = element_text(size=20)) + ylim(c(-2.1,2.25))
c=ggplot(bstdata.fst,aes(y=scale(fst),x=scale(envd))) + geom_point() +stat_smooth(data=bstdata.fst,method="lm",se=T)+theme_bw()+guides(color=F) +ylab("Genetic Distance")+xlab("Env Distance")+theme(axis.text = element_text(size=14), axis.title = element_text(size=20)) + ylim(c(-2.1,2.25))
d=ggplot(bstdata.fst,aes(y=scale(fst),x=scale(envd)+scale(geod),colour=factor(barrier))) + geom_point() +scale_color_manual(values=c("orange","blue"))+stat_smooth(data=bstdata.fst,method="lm",se=T)+theme_bw()+guides(color=F) +ylab("Genetic Distance")+xlab("Geo + Env Distance")+theme(axis.text = element_text(size=14), axis.title = element_text(size=20)) + ylim(c(-2.1,2.25))

grid.arrange(arrangeGrob(a,b,c,d, nrow=1, ncol=4))



