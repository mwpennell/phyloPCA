## # PGLS simulations

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
