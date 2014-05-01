require(phytools)
require(geiger)
ntips <- 100
node <- ntips+1
L <- 2*ntips-2
edge=matrix(ncol=2,nrow=2*ntips-2)
edge[,1]=sort(rep((ntips+1):(2*ntips-1),2))
edge[seq(1,2*ntips-2,2),2]=1:(ntips-1)
edge[seq(2,2*ntips-2,2),2]=(ntips+2):(2*ntips)
edge[2*ntips-2,2]=ntips
edge.length=rep(1/ntips,L)
edge.length[c(seq(1,L,2),L)]=1-1/ntips*0:(ntips-1)
tree=list('edge'=edge,'tip.label'=1:ntips,'edge.length'=edge.length,'Nnode'=ntips-1)
class(tree)='phylo'
is.ultrametric(tree)
try(plot(tree))
plot(tree)
X=fastBM(tree,nsim=100)
pc <- princomp(X)
plot(pc$scores[,1],pc$scores[,2])
nsims <- 10
OUT <- NULL
OU.OUT <- NULL

for(i in 1:nsims){
  X <- fastBM(tree,nsim=100)
  pc <- princomp(X)
  apply(X,1,function(x) (abs(x)-min(x))/(max(x)-min(x)))
  fit1 <- fitContinuous(tree,pc$scores[,1:7],model="BM")
  fit2 <- fitContinuous(tree,pc$scores[,1:7],model="OU")
  #spars <- sapply(1:10,function(x) fit1[[x]]$opt$sigsq)
  bm.lik <- sapply(1:7,function(x) fit1[[x]]$opt$aic)
  ou.lik <- sapply(1:7,function(x) fit2[[x]]$opt$aic)
  OUT <- rbind(OUT,bm.lik)
  OU.OUT <- rbind(OU.OUT,ou.lik)
}


plot(OUT[,1],ylim=c(0,5))
plot(1:ncol(OUT),apply(OUT,2,mean))
