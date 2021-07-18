#Figure 2

#XPEHH

#Create for of XPEHH for top 50 scaffolds. Highlighting scaffold 27

setwd("path/to/folder")
xpehh=read.table("xpehh.top50.2020.head",header=TRUE)
library(dplyr)

#adapted from: https://www.r-graph-gallery.com/101_Manhattan_plot.html

don <- xpehh %>% 
  
  # Compute chromosome size
  group_by(chr) %>% 
  summarise(chr_len=max(pos)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(xpehh, ., by=c("chr"="chr")) %>%
  
  # Add a cumulative position of each SNP
  arrange(chr, pos) %>%
  mutate( BPcum=pos+tot)

axisdf = don %>% group_by(chr) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

plot.xpehh=ggplot(don, aes(x=BPcum, y=normxpehh)) +
  
  # Show all points
  geom_point( aes(color=as.factor(chr)), alpha=0.8, size=1.3) +
  scale_color_manual(values = c(rep(c("grey70", "grey40"), 13 ),"orange",rep(c("grey70","grey40"),12))) +
  
  # custom X axis:
  scale_x_continuous( label = axisdf$chr, breaks= axisdf$center ) +
  scale_y_continuous(expand = c(0, 0.2) ) +     # remove space between plot area and x axis
  
  # Custom the theme:
  theme_bw() +
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  ) +
  xlab("Scaffold") + ylab("Normalised XPEHH")

s27.xpehh=xpehh[which(xpehh$chr == 27),]

ggplot(s27.xpehh, aes(pos,normxpehh))+geom_point(color="orange")+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank())+theme(axis.line = element_line(color="black"))+labs(x="Position (kb)",y="Normalised XPEHH")+ geom_hline(yintercept=2.815, linetype="dashed")+annotate(geom="text", y=3.2, x=220000, label="XPEHH Threshold")

#create histogram of xpehh values showing outlier region
xpehh.all=read.table("xpehh.all.2020",header=TRUE)
ggplot(xpehh.all,aes(normxpehh))+geom_histogram(bins=60, color="white",size=0.02,fill="black")+annotate("rect", xmin = 2.815, xmax = 6, ymin = -50, ymax = 3900,alpha = 0.1, fill="red")+theme_classic()+xlab("Normalised XPEHH")+ylab("Count")

###Nucleotide diversity

#pi plot + histogram of z-scores
setwd("path/to/folder")
#read in stacks output, separate pi values for ins and out, perform zscore, subset for s27
system("grep 'Batch' SNP.msp.recode.sumstats.tsv | sed -E 's/# //' | sed -E 's/ //g' > header.txt")
system("cat SNP.msp.recode.sumstats.tsv | grep -v '#' | grep 'INS' | cat header.txt - > INS.pi")
system("cat SNP.msp.recode.sumstats.tsv | grep -v '#' | grep 'OUT' | cat header.txt - > OUT.pi")

INS.pi=read.table("INS.pi",header=TRUE)
OUT.pi=read.table("OUT.pi",header=TRUE)

#Extract smoothed pi values for PAG (INS) and GO (OUT) populations
myvars=c("Chr","BP","SmoothedPi")
all.pi=cbind(INS.pi[myvars],OUT.pi$SmoothedPi)
colnames(all.pi)=c("Chr","BP","INS","OUT")

#Calc Z-score
all.pi$log2=log2(all.pi$INS/all.pi$OUT)
mean(all.pi$log2[is.finite(all.pi$log2)],na.rm = TRUE)
sd(all.pi$log2[is.finite(all.pi$log2)],na.rm = TRUE)
all.pi$zscore = (all.pi$log2 - (mean(all.pi$log2[is.finite(all.pi$log2)],na.rm = TRUE)))/(sd(all.pi$log2[is.finite(all.pi$log2)],na.rm = TRUE))
hist(all.pi$zscore)

#Use threshold of 5stdev from mean
NROW(all.pi$BP[which(all.pi$zscore <= -5)])
NROW(all.pi$BP[which(all.pi$zscore == -Inf)])

write.table(all.pi,"all.pi.zscores",quote = FALSE, row.names = FALSE)
outliers.pi=all.pi[which(all.pi$zscore <= -5 | all.pi$zscore == -Inf),]
write.table(outliers.pi,"outliers.pi.zscores",quote = FALSE, row.names = FALSE)

ggplot(all.pi,aes(x=zscore))+geom_histogram(bins=100,color="white",size=0.02,fill="black")+annotate("rect", xmin = -11.5, xmax = -5, ymin = -200, ymax = 10200,alpha = 0.1, fill="red")+theme_classic()+xlab(expression("Z-Score Log"[2]*"("*pi["PAG"]*"/"*pi["GO"]*")"))+ylab("Count")

#Extract Scaffold27 pi values
s27.pi=all.pi[which(all.pi$Chr == "scaffold27|size1953023"),]

#check pi outlier region covers all SNPs in window
NROW(s27.pi$BP[which(s27.pi$zscore <= -5)])
s27.pi$BP[which(s27.pi$zscore <= -5)]
NROW(s27.pi$BP[which(s27.pi$BP >= 1447130 & s27.pi$BP <= 1575802)])

#plot pi for s27 with shaded outlier region. As per Hohenlohe (2010), outlier window is extended +/- 2sigma (50kb) as SNPs in these locations can contribute to the significant average of the outliers at the limits
plot.pi=ggplot(s27.pi,aes(BP/1000))+geom_line(aes(y=INS,color="INS"))+geom_line(aes(y=OUT,color="OUT"))+geom_point(aes(y=INS,color="INS"))+geom_point(aes(y=OUT,color="OUT"))+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank())+theme(axis.line = element_line(color="black"))+labs(x="Position (kb)",y=expression("Nucleotide Diversity ("*pi*")"))+annotate("rect",xmin = 1397.130,xmax=1625.802,ymin=0,ymax=0.42,alpha=0.2)+scale_y_continuous(expand = c(0, 0.001) )+ theme(legend.position = "none")


###FST

system("cat SNP.msp.recode.fst_INS-OUT.tsv | sed -E 's/# //' | sed -E 's/ //g' > SNP.msp.recode.fst_INS-OUT.rehead.tsv")
all.fst=read.table("./SNP.msp.recode.fst_INS-OUT.rehead.tsv",header=TRUE)
s27.fst=all.fst[which(all.fst$Chr == "scaffold27|size1953023"),]
ggplot(s27.fst)+geom_line(aes(BP,SmoothedAMOVAFst))+geom_point(aes(BP,AMOVAFst))

plot.fst=ggplot(s27.fst,aes(BP/1000)) + geom_line(aes(y=SmoothedAMOVAFst, colour="Smoothed FST"))+geom_point(aes(y=AMOVAFst,colour="Point FST")) + scale_colour_manual(values=c("orange","black"))+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank())+theme(axis.line = element_line(color="black"))+labs(x="Position (kb)",y=expression("F"["ST"]))+annotate("rect",xmin = 1481.656,xmax=1625.802,ymin=-0.01,ymax=0.55,alpha=0.2)+annotate("rect",xmin = 920.751,xmax=1049.409,ymin=-0.01,ymax=0.55,alpha=0.2)+scale_y_continuous(expand = c(0, 0.001) )+ theme(legend.position = "none")

# x axis closer to zero, plan for labels, fix plot size for both plots, more ticks for fst
#check outlier window covers all SNPs - 33 SNP outliers range from 1531656-1575802 & 25 SNPs from 970751-999409, as per Hohenlohe (2010), outlier window is extended +/- 2sigma (50kb) as SNPs in these locations can contribute to the significant average of the outliers at the limits
nrow(s27.fst[which(s27.fst$BP >= 1531656 & s27.fst$BP <= 1575802),])
nrow(s27.fst[which(s27.fst$BP >= 970751 & s27.fst$BP <= 999409),])

g.pi=ggplotGrob(plot.pi)
g.fst=ggplotGrob(plot.fst)
g=rbind(g.pi,g.fst,size="first")
grid.newpage()
grid.draw(g)
pdf("s27_pi_fst_comb.pdf", width=7.5, height=7)
grid.draw(g)
dev.off()

pdf("s27.xpehh.pdf",width=7.5, height=3.5)
plot.xpehh
dev.off()


plot.xpehh=ggplot(don, aes(x=BPcum, y=normxpehh)) +
  
  # Show all points
  geom_point( aes(color=as.factor(chr)), alpha=0.8, size=1.3) +
  scale_color_manual(values = c(rep(c("grey70", "grey40"), 13 ),"orange",rep(c("grey70","grey40"),12))) +
  
  # custom X axis:
  scale_x_continuous( label = axisdf$chr, breaks= axisdf$center, expand=c(0,2000000) ) +
  scale_y_continuous(expand = c(0, 0.2) ) +     # remove space between plot area and x axis
  
  # Custom the theme:
  theme_classic() + theme(legend.position = "none")+
  xlab("Scaffold") + ylab("Normalised XPEHH")


#Venn diagram

library(eulerr)
fit=euler(c("FST"=1070,"PI"=215,"XPEHH"=244,"FST&PI"=60,"FST&XPEHH"=150,"PI&XPEHH"=9,"FST&PI&XPEHH"=33))
plot(fit)





