#calc of bifurcation plots
library(rehh)
library(vcfR)

hho <- data2haplohh(hap_file = "scaf27.out.head.vcf",polarize_vcf = FALSE)
#central SNP for bifurcation analyses is scaffold27_1539277
#Because RAD SNPs are patchy, present on SNP by SNP basis rather than distance
hho@positions=c(1:356)
furcationo <- calc_furcation(hho,mrk = 271)

plot(furcationo, xlim=c(241,295),lwd=0.4)

hhi<- data2haplohh(hap_file = "scaf27.ins.head.vcf",polarize_vcf = FALSE)
hhi@positions=c(1:356)
furcationi <- calc_furcation(hhi,mrk = 271)

plot(furcationi, xlim=c(241,295),lwd=0.4)
