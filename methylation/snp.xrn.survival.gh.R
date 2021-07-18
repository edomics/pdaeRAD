###Pdae15888 methylation - is SNP-based differential methylation more correlated to larval survival
library(ggplot2)
#Create vector of heat survival index for individuals F1,F3,F4,A5,A6,A8
hscore=c(0.789439953, 0.998209258, 0.22654263, 0.954202839, 1.477754999, 1.300719112)

#Read S441 methlyation table (3 filters - excludes the filter for potential positions that contain SNPs)
s441=read.table("scaf441.3cfilt.tsv",header=TRUE)

#Crop table to Pdae15888 
s441.xrn=s441[which(s441$pos >= 287182 & s441$pos <= 363229),]

#Calculate means for AD and F samples
s441.xrn$amean=rowMeans(as.matrix(s441.xrn[,3:10]))
s441.xrn$fmean=rowMeans(as.matrix(s441.xrn[,11:18]))

#Identify supported/high confidence SNP positions based on the following logic
#A methylatable (for want of a better word) position that is homozygous for a CpG disruptive SNP \
#in an individual will have zero methylation for both corresponding adult and sperm samples

#Record positions that have zero methylation for both adult and sperm samples
s441.xrn$A5z=0
s441.xrn$A6z=0
s441.xrn$A7z=0
s441.xrn$A8z=0
s441.xrn$F1z=0
s441.xrn$F2z=0
s441.xrn$F3z=0
s441.xrn$F4z=0
s441.xrn$A5z[which(s441.xrn$A.A5 == 0 & s441.xrn$A.S5 == 0)]=1
s441.xrn$A6z[which(s441.xrn$A.A6 == 0 & s441.xrn$A.S6 == 0)]=1
s441.xrn$A7z[which(s441.xrn$A.A7 == 0 & s441.xrn$A.S7 == 0)]=1
s441.xrn$A8z[which(s441.xrn$A.A8 == 0 & s441.xrn$A.S8 == 0)]=1
s441.xrn$F1z[which(s441.xrn$F.A1 == 0 & s441.xrn$F.S1 == 0)]=1
s441.xrn$F2z[which(s441.xrn$F.A2 == 0 & s441.xrn$F.S2 == 0)]=1
s441.xrn$F3z[which(s441.xrn$F.A3 == 0 & s441.xrn$F.S3 == 0)]=1
s441.xrn$F4z[which(s441.xrn$F.A4 == 0 & s441.xrn$F.S4 == 0)]=1
s441.xrn$zerosum=0
s441.xrn$zerosum=s441.xrn$A5z+s441.xrn$A6z+s441.xrn$A7z+s441.xrn$A8z+s441.xrn$F1z+s441.xrn$F2z+s441.xrn$F3z+s441.xrn$F4z
#replace above with a row sum?

#Calculate correlation coefficient for each pos vs hsi

hsi=c("F.S1","F.S3","F.S4","A.S5","A.S6","A.S8")
s441.xrn.hsi=s441.xrn[hsi]
mcor=lapply(1:length(s441.xrn.hsi$F.S1), function(x) cor(as.numeric(s441.xrn.hsi[x,]),hscore, method="pearson"))
s441.xrn$hsicor=unlist(mcor)

#Calculate mean correlation coefficient for different SNP prevalence
mean(s441.xrn$hsicor[which(s441.xrn$zerosum == 0)],na.rm = TRUE) #0.1330131
mean(s441.xrn$hsicor[which(s441.xrn$zerosum >= 1)],na.rm = TRUE) #0.5720017
mean(s441.xrn$hsicor[which(s441.xrn$zerosum >= 2)],na.rm = TRUE) #0.7754107
mean(s441.xrn$hsicor[which(s441.xrn$zerosum >= 3)],na.rm = TRUE) #0.780757

#generate a box plot
ggplot(s441.xrn,aes(as.factor(zerosum),hsicor, color=as.factor(zerosum)))+geom_boxplot()+geom_jitter(alpha=0.3)+theme_bw()

#Perform permutation analysis (can perform at different SNP prevalence values).
#Test if the difference in methylation between SNP-impacted sites and non-SNP-impacted sites is greater than expected by chance. 
nperms=10000000
perm.result=numeric(nperms)
for(i in 1:nperms)
{
  s441.xrn$rs=sample(s441.xrn$hsicor,length(s441.xrn$hsicor),replace=FALSE)
  perm.result[i]=mean(s441.xrn$rs[which(s441.xrn$zerosum >= 1)],na.rm = TRUE)-mean(s441.xrn$rs[which(s441.xrn$zerosum == 0)],na.rm = TRUE)
}
hist(perm.result)


