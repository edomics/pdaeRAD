#modified from https://www.analyticsvidhya.com/blog/2016/03/practical-guide-principal-component-analysis-python/
#calculation of environmental distance using PCA approach
env=read.table("env_var.txt",header=TRUE)
#use prcomp() to calc pca. By default, mean centred on - and scale=TRUE normalise to stdev = 1
prin_comp <- prcomp(env,scale.=TRUE)
#rotation give PC loadings
prin_comp$rotation
#plot with loadings
biplot(prin_comp, scale = 0)
#compute std_dev of each PC
std_dev <- prin_comp$sdev
#calc variance
pr_var <- std_dev^2
pr_var
#proportion of variance explained by each PC
prop_varex <- pr_var/sum(pr_var)
prop_varex
#table of PC values for calc of Euclidean distance
prin_comp$x

#modified from Saenz-Agudelo 2015
library(dplyr)
setwd("/Users/egs7/Documents/RAD2020/SUPPFIGS/MODEL/MODEL2020/")
bstdata=read.table("model_data.txt",header=TRUE)
fst=read.table("SNP.distLD.bi.fst.comb",header=FALSE)
names(fst)=c("site.comp","fst")
bstdata.fst=left_join(bstdata,fst, by = c("site"="site.comp"))

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

#Now we ran a signficance test based on permutations. This function was adapted from the MMRR function developed by Wang (Wang 2013). 
# Define matrices 
# best model was model 3
Y = cbind(scale(bstdata.fst$fst))
X = cbind((bstdata.fst$barrier),scale(bstdata.fst$geod),scale(bstdata.fst$envd))
colnames(X)=c("barrier","geod","envd")

#define function MMRR that compute regression coefficients and test statistics
MMRR<-function(Y,X,nperm=999){  
  nrowsY<-nrow(Y)
  y<-(Y)
  if(is.null(colnames(X)))colnames(X)<-paste("X",1:length(X),sep="")
  Xmats<-(X)
  fit<-lm(y~factor(Xmats[,1])+Xmats[,3]+Xmats[,2])
  coeffs<-fit$coefficients
  summ<-summary(fit)
  r.squared<-summ$r.squared
  tstat<-summ$coefficients[,"t value"]
  Fstat<-summ$fstatistic[1]
  tprob<-rep(1,length(tstat))
  Fprob<-1
  
  #perform permutations
  for(i in 1:nperm){
    rand<-sample(1:nrowsY)
    Yperm<-Y[rand]
    yperm<-(Yperm)
    fit<-lm(yperm~factor(Xmats[,1])+Xmats[,3]+Xmats[,2])
    summ<-summary(fit)
    Fprob<-Fprob+as.numeric(summ$fstatistic[1]>=Fstat)
    tprob<-tprob+as.numeric(abs(summ$coefficients[,"t value"])>=abs(tstat))
  }
  
  #return values
  tp<-tprob/(nperm+1)
  Fp<-Fprob/(nperm+1)
  names(r.squared)<-"r.squared"
  names(coeffs)<-c("Intercept",colnames(X))
  names(tstat)<-paste(c("Intercept",colnames(X)),"(t)",sep="")
  names(tp)<-paste(c("Intercept",colnames(X)),"(p)",sep="")
  names(Fstat)<-"F-statistic"
  names(Fp)<-"F p-value"
  return(list(r.squared=r.squared,
              coefficients=coeffs,
              tstatistic=tstat,
              tpvalue=tp,
              Fstatistic=Fstat,
              Fpvalue=Fp))
}
#run the function with 9999 permutations
MMRR (Y,X,nperm=9999)

#Plotting figures for entire dataset
#with CIs
ggplot(bstdata.fst,aes(y=scale(fst),x=scale(envd)+scale(geod),colour=factor(barrier))) + geom_point() +scale_color_manual(values=c("orange","blue"))+stat_smooth(data=bstdata,method="lm",se=T)+theme_bw()+guides(color=F) +ylab("Genetic Distance")+xlab("Geo + Env Distance")+theme(axis.text = element_text(size=14), axis.title = element_text(size=20))
#without CIs
ggplot(bstdata.fst,aes(y=scale(fst),x=scale(envd)+scale(geod),colour=factor(barrier))) + geom_point() +scale_color_manual(values=c("orange","blue"))+stat_smooth(data=bstdata,method="lm",se=F)+theme_bw()+guides(color=F) +ylab("Genetic Distance")+xlab("Geo + Env Distance")+theme(axis.text = element_text(size=14), axis.title = element_text(size=20))

