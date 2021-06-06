library(ggplot2)

setwd("/path/to/directory")
###plots for multiple K values
k2=read.table("props.k2.both.out",header=FALSE)
names(k2)=c("Index","IND","SITE","PopName","Prop")
k2$SITE=factor(k2$SITE,c("KA","KB","QB","QA","QC","WC","DH","RG","SB","RK","MB","FD","TA","TB"))
struc2=ggplot(data=k2,aes(x=Index,y=Prop,fill=PopName,width=1,ylim(0,1))) + geom_bar(stat="identity")
struc2final=struc2+facet_grid(~SITE, scales="free_x") + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour="black"),axis.line.x = element_blank(),axis.title.x = element_blank(),axis.text.x = element_blank(),axis.ticks.x = element_blank(),panel.spacing.x =unit(0.1,"lines"))+scale_y_continuous(expand=c(0,0))+scale_fill_manual(values=c("#F8766D","#00B0F6"))+ theme(legend.position = "none") + labs(title = "K = 2" ) + ylab(element_text("Proportion"))

k3=read.table("props.k3.both.out",header=FALSE)
names(k3)=c("Index","IND","SITE","PopName","Prop")
k3$SITE=factor(k3$SITE,c("KA","KB","QB","QA","QC","WC","DH","RG","SB","RK","MB","FD","TA","TB"))
struc3=ggplot(data=k3,aes(x=Index,y=Prop,fill=PopName,width=1,ylim(0,1))) + geom_bar(stat="identity")
struc3final=struc3+facet_grid(~SITE, scales="free_x") + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour="black"),axis.line.x = element_blank(),axis.title.x = element_blank(),axis.text.x = element_blank(),axis.ticks.x = element_blank(),panel.spacing.x =unit(0.1,"lines"))+scale_y_continuous(expand=c(0,0))+scale_fill_manual(values=c("#F8766D","#00B0F6","#00BF7D"))+ theme(legend.position = "none")+ labs(title = "K = 3" ) + ylab(element_text("Proportion"))

k4=read.table("props.k4.both.out",header=FALSE)
names(k4)=c("Index","IND","SITE","PopName","Prop")
k4$SITE=factor(k4$SITE,c("KA","KB","QB","QA","QC","WC","DH","RG","SB","RK","MB","FD","TA","TB"))
struc4=ggplot(data=k4,aes(x=Index,y=Prop,fill=PopName,width=1,ylim(0,1))) + geom_bar(stat="identity")
struc4final=struc4+facet_grid(~SITE, scales="free_x") + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour="black"),axis.line.x = element_blank(),axis.title.x = element_blank(),axis.text.x = element_blank(),axis.ticks.x = element_blank(),panel.spacing.x =unit(0.1,"lines"))+scale_y_continuous(expand=c(0,0))+scale_fill_manual(values=c("#F8766D","#00B0F6","#00BF7D", "#A3A500"))+ theme(legend.position = "none")+ labs(title = "K = 4" ) + ylab(element_text("Proportion"))

k5=read.table("props.k5.both.out",header=FALSE)
names(k5)=c("Index","IND","SITE","PopName","Prop")
k5$SITE=factor(k5$SITE,c("KA","KB","QB","QA","QC","WC","DH","RG","SB","RK","MB","FD","TA","TB"))
struc5=ggplot(data=k5,aes(x=Index,y=Prop,fill=PopName,width=1,ylim(0,1))) + geom_bar(stat="identity")
struc5final=struc5+facet_grid(~SITE, scales="free_x") + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour="black"),axis.line.x = element_blank(),axis.title.x = element_blank(),axis.text.x = element_blank(),axis.ticks.x = element_blank(),panel.spacing.x =unit(0.1,"lines"))+scale_y_continuous(expand=c(0,0))+scale_fill_manual(values=c("#F8766D","#00B0F6","#00BF7D", "#A3A500", "#E76BF3"))+ theme(legend.position = "none")+ labs(title = "K = 5" ) + ylab(element_text("Proportion"))

library("gridExtra")
grid.arrange(struc2final,struc3final, struc4final, struc5final, ncol=1)

##K6-K10 supp.

k6=read.table("props.k6.both.out",header=FALSE)
names(k6)=c("Index","IND","SITE","PopName","Prop")
k6$SITE=factor(k6$SITE,c("KA","KB","QB","QA","QC","WC","DH","RG","SB","RK","MB","FD","TA","TB"))
struc6=ggplot(data=k6,aes(x=Index,y=Prop,fill=PopName,width=1,ylim(0,1))) + geom_bar(stat="identity")
struc6final=struc6+facet_grid(~SITE, scales="free_x") + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour="black"),axis.line.x = element_blank(),axis.title.x = element_blank(),axis.text.x = element_blank(),axis.ticks.x = element_blank(),panel.spacing.x =unit(0.1,"lines"))+scale_y_continuous(expand=c(0,0))+scale_fill_manual(values=c("#F8766D","#00B0F6","#00BF7D", "#A3A500", "#E76BF3", "#39B600", "#00BFC4", "#9590FF", "#D89000", "#FF62BC"))+ theme(legend.position = "none")+ labs(title = "K = 6" ) + ylab(element_text("Proportion"))


k7=read.table("props.k7.both.out",header=FALSE)
names(k7)=c("Index","IND","SITE","PopName","Prop")
k7$SITE=factor(k7$SITE,c("KA","KB","QB","QA","QC","WC","DH","RG","SB","RK","MB","FD","TA","TB"))
struc7=ggplot(data=k7,aes(x=Index,y=Prop,fill=PopName,width=1,ylim(0,1))) + geom_bar(stat="identity")
struc7final=struc7+facet_grid(~SITE, scales="free_x") + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour="black"),axis.line.x = element_blank(),axis.title.x = element_blank(),axis.text.x = element_blank(),axis.ticks.x = element_blank(),panel.spacing.x =unit(0.1,"lines"))+scale_y_continuous(expand=c(0,0))+scale_fill_manual(values=c("#F8766D","#00B0F6","#00BF7D", "#A3A500", "#E76BF3", "#39B600", "#00BFC4", "#9590FF", "#D89000", "#FF62BC"))+ theme(legend.position = "none")+ labs(title = "K = 7" ) + ylab(element_text("Proportion"))


k8=read.table("props.k8.both.out",header=FALSE)
names(k8)=c("Index","IND","SITE","PopName","Prop")
k8$SITE=factor(k8$SITE,c("KA","KB","QB","QA","QC","WC","DH","RG","SB","RK","MB","FD","TA","TB"))
struc8=ggplot(data=k8,aes(x=Index,y=Prop,fill=PopName,width=1,ylim(0,1))) + geom_bar(stat="identity")
struc8final=struc8+facet_grid(~SITE, scales="free_x") + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour="black"),axis.line.x = element_blank(),axis.title.x = element_blank(),axis.text.x = element_blank(),axis.ticks.x = element_blank(),panel.spacing.x =unit(0.1,"lines"))+scale_y_continuous(expand=c(0,0))+scale_fill_manual(values=c("#F8766D","#00B0F6","#00BF7D", "#A3A500", "#E76BF3", "#39B600", "#00BFC4", "#9590FF", "#D89000", "#FF62BC"))+ theme(legend.position = "none")+ labs(title = "K = 8" ) + ylab(element_text("Proportion"))


k9=read.table("props.k9.both.out",header=FALSE)
names(k9)=c("Index","IND","SITE","PopName","Prop")
k9$SITE=factor(k9$SITE,c("KA","KB","QB","QA","QC","WC","DH","RG","SB","RK","MB","FD","TA","TB"))
struc9=ggplot(data=k9,aes(x=Index,y=Prop,fill=PopName,width=1,ylim(0,1))) + geom_bar(stat="identity")
struc9final=struc9+facet_grid(~SITE, scales="free_x") + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour="black"),axis.line.x = element_blank(),axis.title.x = element_blank(),axis.text.x = element_blank(),axis.ticks.x = element_blank(),panel.spacing.x =unit(0.1,"lines"))+scale_y_continuous(expand=c(0,0))+scale_fill_manual(values=c("#F8766D","#00B0F6","#00BF7D", "#A3A500", "#E76BF3", "#39B600", "#00BFC4", "#9590FF", "#D89000", "#FF62BC"))+ theme(legend.position = "none")+ labs(title = "K = 9" ) + ylab(element_text("Proportion"))


k10=read.table("props.k10.both.out",header=FALSE)
names(k10)=c("Index","IND","SITE","PopName","Prop")
k10$SITE=factor(k10$SITE,c("KA","KB","QB","QA","QC","WC","DH","RG","SB","RK","MB","FD","TA","TB"))
struc10=ggplot(data=k10,aes(x=Index,y=Prop,fill=PopName,width=1,ylim(0,1))) + geom_bar(stat="identity")
struc10final=struc10+facet_grid(~SITE, scales="free_x") + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour="black"),axis.line.x = element_blank(),axis.title.x = element_blank(),axis.text.x = element_blank(),axis.ticks.x = element_blank(),panel.spacing.x =unit(0.1,"lines"))+scale_y_continuous(expand=c(0,0))+scale_fill_manual(values=c("#F8766D","#00B0F6","#00BF7D", "#A3A500", "#E76BF3", "#39B600", "#00BFC4", "#9590FF", "#D89000", "#FF62BC"))+ theme(legend.position = "none")+ labs(title = "K = 10" ) + ylab(element_text("Proportion"))

grid.arrange(struc6final,struc7final, struc8final, struc9final, struc10final, ncol=1)
#saved at 10x10

#calculate admixture props for k=3
#non-PAG subpop at Offshore Qatar
1-mean(k3$Prop[which(k3$SITE == "QB" & k3$PopName == "POP1")]) 
#non-PAG subpop at Inshore Qatar
1-mean(k3$Prop[which(k3$SITE == "QA" & k3$PopName == "POP1")]) 

#calculate GO ancestry in sPAG pops
k3$loc = 0
k3$loc[which(k3$SITE == "QA" | k3$SITE == "QC" | k3$SITE == "WC" | k3$SITE == "DH" | k3$SITE == "RG")] <- "SPAG_INS"
k3$loc[which(k3$SITE == "QB" | k3$SITE == "SB" | k3$SITE == "RK")] <- "SPAG_OFF"
mean(k3$Prop[which(k3$loc == "SPAG_INS" & k3$PopName == "POP2")])
mean(k3$Prop[which(k3$loc == "SPAG_OFF" & k3$PopName == "POP2")])
max(k3$Prop[which(k3$loc == "SPAG_INS" & k3$PopName == "POP2")])
max(k3$Prop[which(k3$loc == "SPAG_OFF" & k3$PopName == "POP2")])
                                                             
#PCA plot - input from plink + pasted labels + modified from (https://speciationgenomics.github.io/pca/)
egval=read.table("plink.eigenval",header=FALSE)
egvec=read.table("plink.eigenvec.labelled",header=FALSE)
egvec$V2=factor(egvec$V2,c("KA","KB","QB","QA","QC","WC","DH","RG","SB","RK","MB","FD","TA","TB"))
pve <- data.frame(PC = 1:20, pve = egval/sum(egval)*100)
#plot for all sites
ggplot(egvec,aes(x=V6,y=V7,color=V2))+geom_point()+theme_bw()+coord_fixed()+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+xlab(paste0("PC1 (", signif(pve$V1[1], 4), "%)")) + ylab(paste0("PC2 (", signif(pve$V1[2], 4), "%)")) +labs(color = "Site")
#plot for all regions
egvec$V3=factor(egvec$V3,c("AG","KW","GO","QO"))
ggplot(egvec,aes(x=V6,y=V7,color=V3))+geom_point()+theme_bw()+coord_fixed()+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+xlab(paste0("PC1 (", signif(pve$V1[1], 4), "%)")) + ylab(paste0("PC2 (", signif(pve$V1[2], 4), "%)")) +labs(color = "Location")
#save as 4x4
ggplot(egvec,aes(x=V6,y=V7,color=V3))+geom_point()+theme_bw()+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+xlab(paste0("PC1 (", signif(pve$V1[1], 4), "%)")) + ylab(paste0("PC2 (", signif(pve$V1[2], 4), "%)")) +labs(color = "Location")+theme(aspect.ratio = 1)+theme(legend.position = "none")



