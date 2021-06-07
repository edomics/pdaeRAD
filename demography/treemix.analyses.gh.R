
###Code to analyse and plot multiple treemix reps
###Requires plotting functions R code provided with Treemix
###https://github.com/joepickrell/pophistory-tutorial/blob/master/example2/plotting_funcs.R
###Modifies some code from here (https://speciationgenomics.github.io/Treemix/)

###Includes plots for dealing with failed Treemix reps.
###These seem to occur when runs get stuck in suboptimal solutions.
###Alternatively, you could repeat failed runs with new seed. I kept them for transparency.

#Plot of max likelihood
setwd("/path/to/directory")
source("plotting_funcs.R")
tmxpref="prefix used for treemix"
#completed.list is a list of completed runs
com=read.table("completed.list",header=FALSE)
#failed.list is a list of failed runs
fail=read.table("failed.list",header=FALSE)
likemat=matrix(0,dim(com)[1],11)
for(line in 1:dim(com)[1]){
  r=com[line,1]
  for(m in 0:10){
    filename=(paste0(tmxpref,m,".r.",r,".llik"))
    ll=readLines(filename)
    x=paste0("Exiting ln\\(likelihood\\) with ",m," migration events\\: ")
    lik=as.numeric(gsub(x,"",ll[2]))
    likemat[line,m+1]=lik
  }
  
}
row.names(likemat)=rownames(likemat, do.NULL = FALSE, prefix = "rep")
colnames(likemat) <- c("m0","m1","m2","m3","m4","m5","m6","m7","m8","m9","m10")

likelong=melt(likemat)
ggplot(likelong,aes(x=Var2,group=Var2,y=value))+geom_boxplot()

#make failed matrix
likematf=matrix(0,dim(fail)[1],6)
for(line in 1:dim(fail)[1]){
  r=fail[line,1]
  for(m in 0:5){
    filename=(paste0(tmxpref,m,".r.",r,".llik"))
    ll=readLines(filename)
    x=paste0("Exiting ln\\(likelihood\\) with ",m," migration events\\: ")
    lik=as.numeric(gsub(x,"",ll[2]))
    likematf[line,m+1]=lik
  }
  
}
row.names(likematf)=rownames(likematf, do.NULL = FALSE, prefix = "rep")
colnames(likematf) <- c("m0","m1","m2","m3","m4","m5")

likelongf=melt(likematf)

#make failed matrix - but dummy up to m10
likematfd=matrix(NA,dim(fail)[1],11)
for(line in 1:dim(fail)[1]){
  r=fail[line,1]
  for(m in 0:5){
    filename=(paste0(tmxpref,m,".r.",r,".llik"))
    ll=readLines(filename)
    x=paste0("Exiting ln\\(likelihood\\) with ",m," migration events\\: ")
    lik=as.numeric(gsub(x,"",ll[2]))
    likematfd[line,m+1]=lik
  }
  
}
row.names(likematfd)=rownames(likematfd, do.NULL = FALSE, prefix = "rep")
colnames(likematfd) <- c("m0","m1","m2","m3","m4","m5","m6","m7","m8","m9","m10")

likelongfd=melt(likematfd)





par(mfrow=c(3,10))
for(rep in 1:30){
  plot_tree(cex=0.8,paste0(tmxpref,"1.r.",rep))
  title(paste(round((c[rep,2]*100),digits=2),"% exp"))
}

#% explained

c=matrix(0,dim(com)[1],11)
row.names(c)=rownames(c, do.NULL = FALSE, prefix = "rep")
colnames(c) <- c("m0","m1","m2","m3","m4","m5","m6","m7","m8","m9","m10")

for(line in 1:dim(com)[1]){
  r=com[line,1]
  for(m in 0:10){
    pexplain=get_f(paste0(tmxpref, m,".r.", r))
    c[line,m+1]=pexplain
  }
}
long=melt(c)
ggplot(long,aes(x=Var2,group=Var2,y=value))+geom_boxplot()

#max likelihood trees for m0 to m10
par(mfrow=c(3,4))
for(m in 0:10){
  plot_tree(cex=0.8,paste0(tmxpref,m,".r.28"))
  title(paste("m = ",m,":",round((c[17,m+1]*100),digits=1),"% exp"))
}

#or... with residuals
par(mfrow=c(3,4))
for(m in 0:5){
  plot_tree(cex=1,paste0(tmxpref,m,".r.28"))
  title(paste("m = ",m,":",round((c[17,m+1]*100),digits=1),"% exp"))
  plot_resid(paste0(tmxpref, m,".r.28"),pop_order = "pop.order")
}

#all rep plots for a single m with % exp

c.all=matrix(0,30,6)
row.names(c.all)=rownames(c.all, do.NULL = FALSE, prefix = "rep")
colnames(c.all) <- c("m0","m1","m2","m3","m4","m5","m6","m7","m8","m9","m10")

for(j in 1:30){
  for(i in 0:5){
    pexplain=get_f(paste0(tmxpref, i,".r.", j))
    c.all[j,i+1]=pexplain
  }
}

#set desired m for tree plots
mview=0
par(mfrow=c(3,10))
for(rep in 1:30){
  plot_tree(cex=0.8,paste0(tmxpref, mview,".r.",rep))
  title(paste(round((c.all[rep,mview+1]*100),digits=2),"% exp"))
}

#set desired m for residual plots
mview=1
par(mfrow=c(3,10))
for(rep in 1:30){
  plot_resid(paste0(tmxpref, mview,".r.",rep),pop_order = "pop.order")
  title(paste(round((c.all[rep,mview+1]*100),digits=2),"% exp"))
}
