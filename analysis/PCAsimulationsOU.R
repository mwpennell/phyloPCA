## # OU simulations

## Load packages
require(phytools)
require(geiger)
require(phylolm)
require(TreeSim)
require(foreach)
require(doMC)
require(MASS)
require(nlme)

## Read in functions
source("PCAfunctions.R")

## Set the number of cores for parallel analysis
registerDoMC(cores=4)

## Specify the number of tips, traits and a vector specifying which traits to fit models to.
ntips <- 50
ntraits <- 30
trait.seq= c(1, 2, 5, 10, 20, 30)

## This function simulates a phylogenetic tree using a birth-death model with extinction = 0. The tree is rescaled to unit height. It then simulates a trait variance-covariance matrix by simulating eigenvalues from an exponential distribution with lambda = 1/100. These are then used to simulate under multivariate Brownian motion on the phylogeny. Note that simulating the data for a large number of tips and traits is computationally intensive (requires a ntraits x ntips covariance matrix).

#Skipped running and just load saved simulations
#res <- foreach(i=1:100) %dopar% simandfitPCs.ind.ou(ntips, ntraits, alpha=2, sig2=1, trait.seq=trait.seq)
#save(res, file=paste("simsou", ntips, ntraits,".rds", sep="_"))
load(paste("simsou", ntips, ntraits,".rds", sep="_"))

## Rearrange the output to be usable for plotting
res.all <- lapply(1:length(res[[1]]), function(x) sapply(1:length(res), function(y) res[[y]][[x]]))
res.all <- lapply(res.all, function(x){x <- t(x); colnames(x) <- colnames(res[[1]][[1]]); x})

## Create boxplots of Akaike weights for each transformation (raw, pca, and phylo  pca) for each of the traits in trait.seq, for each of the models (BM, OU, EB).
par(mfrow=c(3,3), mar=c(4, 4.1, 2.1, 0.1))
boxplot.prep <- function(x){
  tmp <- sapply(res.all, function(y) y[,x])
  colnames(tmp) <- paste("",trait.seq, sep="")
  tmp
}
bx <- lapply(1:ncol(res.all[[1]]), function(x) boxplot(boxplot.prep(x), ylab="AICw",xlab="Trait", main=colnames(res.all[[1]])[x]))

## Create plot of contrasts through time
nsims = 50
alpha=10
sig2=1

dat <- lapply(1:nsims, function(x) sim.tree.pcs.ind.ou(ntips, ntraits, alpha, sig2))
btimes <- lapply(1:nsims, function(x) 1 - branching.times(dat[[x]]$tree))
picrawdat <- lapply(1:nsims, function(x) sapply(1:ntraits, function(y) pic(dat[[x]]$raw[,y], dat[[x]]$tree)))
picpcdat <- lapply(1:nsims, function(x) sapply(1:ntraits, function(y) pic(dat[[x]]$pc[,y], dat[[x]]$tree)))
picppcdat <- lapply(1:nsims, function(x) sapply(1:ntraits, function(y) pic(dat[[x]]$ppc[,y], dat[[x]]$tree)))

raw.bytrait <- lapply(1:ntraits, function(x) lapply(1:nsims, function(y) dat[[y]]$raw[,x]))
dttraw <- lapply(1:ntraits, function(y) lapply(1:nsims, function(x) dtt(dat[[x]]$tree, raw.bytrait[[y]][[x]], plot=FALSE)))
dispraw <- lapply(1:ntraits, function(x) unlist(lapply(1:nsims, function(y) dttraw[[x]][[y]]$dtt)))
times <- lapply(1:ntraits, function(x) unlist(lapply(1:nsims, function(y) dttraw[[x]][[y]]$times)))

pc.bytrait <- lapply(1:ntraits, function(x) lapply(1:nsims, function(y) dat[[y]]$pc[,x]))
dttpc <- lapply(1:ntraits, function(y) lapply(1:nsims, function(x) dtt(dat[[x]]$tree, pc.bytrait[[y]][[x]], plot=FALSE)))
disppc <- lapply(1:ntraits, function(x) unlist(lapply(1:nsims, function(y) dttpc[[x]][[y]]$dtt)))

ppc.bytrait <- lapply(1:ntraits, function(x) lapply(1:nsims, function(y) dat[[y]]$ppc[,x]))
dttppc <- lapply(1:ntraits, function(y) lapply(1:nsims, function(x) dtt(dat[[x]]$tree, ppc.bytrait[[y]][[x]], plot=FALSE)))
dispppc <- lapply(1:ntraits, function(x) unlist(lapply(1:nsims, function(y) dttppc[[x]][[y]]$dtt)))




par(mfrow=c(2,3))
pal = rainbow
alph=10
cex=0.5

plot(c(0,1), c(0,3), type='n', main="raw", xlab="time", ylab="contrasts")
garbage <- lapply(1:ntraits, function(x) lapply(1:nsims, function(y) points(btimes[[y]], abs(picrawdat[[y]][,x]), col=makeTransparent(pal(ntraits)[x], alpha=alph),bg=makeTransparent(pal(ntraits)[x], alpha=alph), cex=cex, pch=21)))
picraw.bytrait <- lapply(1:ntraits, function(x) sapply(1:nsims, function(y) abs(picrawdat[[y]][,x])))
garbage <- lapply(1:ntraits,function(x) abline(lm(unlist(as.data.frame(picraw.bytrait[[x]]))~unlist(btimes)), col=pal(ntraits)[x], lwd=1.5))

plot(c(0,1), c(0,3), type='n', main="pc", xlab="time", ylab="contrasts")
garbage <- lapply(1:ntraits, function(x) lapply(1:nsims, function(y) points(btimes[[y]], abs(picpcdat[[y]][,x]), col=makeTransparent(pal(ntraits)[x], alpha=alph),bg=makeTransparent(pal(ntraits)[x], alpha=alph), cex=cex, pch=21)))
picpc.bytrait <- lapply(1:ntraits, function(x) sapply(1:nsims, function(y) abs(picpcdat[[y]][,x])))
garbage <- lapply(1:ntraits,function(x) abline(lm(unlist(as.data.frame(picpc.bytrait[[x]]))~unlist(btimes)), col=pal(ntraits)[x], lwd=1.5))


plot(c(0,1), c(0,3), type='n', main="ppc", xlab="time", ylab="contrasts")
garbage <- lapply(1:ntraits, function(x) lapply(1:nsims, function(y) points(btimes[[y]], abs(picppcdat[[y]][,x]), col=makeTransparent(pal(ntraits)[x], alpha=alph),bg=makeTransparent(pal(ntraits)[x], alpha=alph), cex=cex, pch=21)))
picppc.bytrait <- lapply(1:ntraits, function(x) sapply(1:nsims, function(y) abs(picppcdat[[y]][,x])))
garbage <- lapply(1:ntraits,function(x) abline(lm(unlist(as.data.frame(picppc.bytrait[[x]]))~unlist(btimes)), col=pal(ntraits)[x], lwd=1.5))

plot(c(0,1), c(0,1.5), type="n", main="raw", xlab="time", ylab="disparity")
garbage <- lapply(1:ntraits, function(x) lines.loess(times[[x]], dispraw[[x]], col=pal(ntraits)[x]))

plot(c(0,1), c(0,1.5), type="n", main="pc", xlab="time", ylab="disparity")
garbage <- lapply(1:ntraits, function(x) lines.loess(times[[x]], disppc[[x]], col=pal(ntraits)[x]))

plot(c(0,1), c(0,1.5), type="n", main="ppc", xlab="time", ylab="disparity")
garbage <- lapply(1:ntraits, function(x) lines.loess(times[[x]], dispppc[[x]], col=pal(ntraits)[x]))
