## # Simulations exploring the effects of different types of PCA on model-based phylogenetic inference.

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
#res <- foreach(i=1:100) %dopar% simandfitPCs(ntips, ntraits, trait.seq=trait.seq)
#save(res, file=paste("sims", ntips, ntraits,".rds", sep="_"))
load(paste("sims", ntips, ntraits,".rds", sep="_"))

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


##
registerDoMC(cores=8)
trait.seq <- c(2, 3, 5, 10, 15, 20, 30)
pglsres <- foreach(i=1:50) %dopar% simandfitPGLS(100, 30, trait.seq, foc2trans="BM", rescalepar=1, me=0)
rawpgls <- sapply(pglsres, function(x) x$raw[,2])
pcpgls  <- sapply(pglsres, function(x) x$pc[,2])
ppcpgls <- sapply(pglsres, function(x) x$ppc[,2])
simdat <- lapply(pglsres, function(x) x$simdat)
pairs((simdat[[1]]$pc[1:50, 1:10]))
plot(simdat[[10]]$pc[,2],simdat[[1]]$pc[,3])

par(mfrow=c(1,3))
plot(ecdf(rawpgls[1,]), ylim=c(0,1), xlim=c(0,1))
lapply(1:length(trait.seq), function(x) lines(ecdf(rawpgls[x,]), col=heat.colors(13)[x]))
curve(1*x, add=TRUE)
plot(ecdf(pcpgls[1,]), ylim=c(0,1), xlim=c(0,1))
lapply(1:length(trait.seq), function(x) lines(ecdf(pcpgls[x,]), col=heat.colors(13)[x]))
curve(1*x, add=TRUE)
plot(ecdf(rawpgls[1,]), ylim=c(0,1), xlim=c(0,1))
lapply(1:length(trait.seq), function(x) lines(ecdf(ppcpgls[x,]), col=heat.colors(13)[x]))
curve(1*x, add=TRUE)

par(mfrow=c(1,3))
yl <- c(-0, 1)
boxplot(t(rawpgls), ylim=yl)
boxplot(t(pcpgls ), ylim=yl)
boxplot(t(ppcpgls), ylim=yl)

plot(pcpgls[,1], ylim=c(0,1), type="n")
lapply(1:20, function(x) lines(pcpgls[,x], col=x))
plot(ppcpgls[,1], ylim=c(0,1), type="n")
lapply(1:20, function(x) lines(ppcpgls[,x], col=x))

