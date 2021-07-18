#Figure 3
setwd("path/to/folder")

#Plot Fst - can use data frame from Figure 2 (scaffold27 plots)
library(egg) #for ggarrange
#function to make sure yaxis is 2decimal places
scaleFUN <- function(x) sprintf("%.2f", x)

all.fst=read.table("./SNP.msp.recode.fst_INS-OUT.rehead.tsv",header=TRUE)

s441.fst=all.fst[which(all.fst$Chr == "scaffold441|size570437"),]

plot.s441fst=ggplot(s441.fst,aes(BP/1000)) + geom_line(aes(y=SmoothedAMOVAFst, colour="Smoothed FST"))+geom_point(aes(y=AMOVAFst,colour="Point FST")) + scale_colour_manual(values=c("orange","black"))+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank())+theme(axis.line = element_line(color="black"))+labs(x="Position (kb)",y=expression("F"["ST"]))+ theme(legend.position = "none")+scale_y_continuous(labels=scaleFUN)

#Plot EHH
ehhcomb=read.table("./EHH/s441_ehh_fig.data", header=TRUE)
ehhcomb$SITE <- factor(ehhcomb$SITE,levels=c("PAG","GO"))
plot.s441ehh=ggplot(ehhcomb,aes(x=POS/1000))+geom_line(data=ehhcomb,aes(y=GULF,colour="blue"))+geom_line(data=ehhcomb,aes(y=OUT,colour="red"))+facet_grid(~SITE)+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank())+theme(axis.line = element_line(color="black"))+labs(x="Position (kb)",y="EHH")+theme(legend.position = "none")

#Plot meth bars

meth=read.table("./BARPLOT/meth4barplot.txt",header=TRUE)
meth$SITE <- factor(meth$SITE,levels=c("PAG","GO"))
plot.meth=ggplot(meth,aes(x=SITE,y=PROP,fill=SITE))+geom_bar(stat="identity")+facet_grid(~POS)+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank())+theme(axis.line = element_line(color="black"))+labs(x="Location",y="Proportion of intact CpG sites")+scale_y_continuous(limits = c(0, 1))+theme(legend.position = "none")

pdf("s441_fst_ehh_meth_comb.pdf", width=15, height=7)
grid.draw(f)
dev.off()


