#stairway plot based on Misha Matz code @ https://github.com/z0on/Adaptive-pathways-of-coral-populations-on-the-Great-Barrier-Reef/blob/master/GBR_stairway.R)
library(ggplot2)
#read in output from stairway plot
stair=read.table("path/to/summaryfile",header = TRUE)
#plot
ggplot(stair,aes(x=year))+geom_line(aes(y=Ne_median), color="red", lwd=2)+geom_ribbon(aes(ymin=Ne_12.5.,ymax=Ne_87.5.), alpha=0.2) + geom_ribbon(aes(ymin=Ne_2.5.,ymax=Ne_97.5.), alpha=0.2) + scale_x_continuous(trans='log10',labels=c("100","1k","10k","100k","1M"),breaks=c(1e+2,1e+3,1e+4,1e+5,1e+6))+scale_y_continuous(trans='log10', limits=c(1800, 210000),breaks=c(2e+3,2e+4,2e+5),labels=c("2k","20k","200k"))+theme_bw()+xlab("Years Before Present")+ylab("Effective Population Size")






