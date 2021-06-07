###Three pop calcs
#threepop ran in treemix package "threepop -i input.treemix.gz > threepop.raw"
#tidied up in command line to tab separate and label pops to regions (i.e. NPAG, SPAG...)
#Purpose is to calculate prop of significant to total test for each region_region comparison
setwd("/path/to/directory")
pop3=read.table("three.raw.poplab",header=FALSE)

pop3$BC="xxx"
pop3$BC[which(pop3$V8 == "NPAG" & pop3$V9 == "NPAG")] = "NPAG_NPAG"
pop3$BC[which(pop3$V8 == "SPAG" & pop3$V9 == "SPAG")] = "SPAG_SPAG"
pop3$BC[which(pop3$V8 == "GO" & pop3$V9 == "GO")] = "GO_GO"
pop3$BC[which(pop3$V8 == "QO" & pop3$V9 == "QO")] = "QO_QO"
pop3$BC[which(pop3$V8 == "NPAG" & pop3$V9 == "SPAG")] = "NPAG_SPAG"
pop3$BC[which(pop3$V9 == "NPAG" & pop3$V8 == "SPAG")] = "NPAG_SPAG"
pop3$BC[which(pop3$V8 == "NPAG" & pop3$V9 == "GO")] = "NPAG_GO"
pop3$BC[which(pop3$V9 == "NPAG" & pop3$V8 == "GO")] = "NPAG_GO"
pop3$BC[which(pop3$V8 == "NPAG" & pop3$V9 == "QO")] = "NPAG_QO"
pop3$BC[which(pop3$V9 == "NPAG" & pop3$V8 == "QO")] = "NPAG_QO"
pop3$BC[which(pop3$V8 == "SPAG" & pop3$V9 == "GO")] = "SPAG_GO"
pop3$BC[which(pop3$V9 == "SPAG" & pop3$V8 == "GO")] = "SPAG_GO"
pop3$BC[which(pop3$V8 == "SPAG" & pop3$V9 == "QO")] = "SPAG_QO"
pop3$BC[which(pop3$V9 == "SPAG" & pop3$V8 == "QO")] = "SPAG_QO"
pop3$BC[which(pop3$V8 == "GO" & pop3$V9 == "QO")] = "GO_QO"
pop3$BC[which(pop3$V9 == "GO" & pop3$V8 == "QO")] = "GO_QO"

pop3$pval=pnorm(pop3$V7)
pop3$pval.adj=p.adjust(pop3$pval,method="BH")
pop3$BC=as.factor(pop3$BC)
pop05=pop3[which(pop3$pval.adj < 0.05),]

site.names=unique(pop3$V1)
for (x in 1:length(site.names)){
  name=site.names[x]
  print(summary(pop3$BC[which(pop3$V1 == name)]))
}

for (x in 1:length(site.names)){
  name=site.names[x]
  print(summary(pop05$BC[which(pop05$V1 == name)]))
}
















