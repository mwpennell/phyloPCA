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
require(reshape)

## Read in functions
source("PCAfunctions.R")

## Set the number of cores for parallel analysis
registerDoMC(cores=3)

## Specify the number of tips, traits and a vector specifying which traits to fit models to.
ntips <- 50
ntraits <- 20

## Select which of the traits to study
trait.seq= c(1, 2, 3, 5, 7, 10, 13, 16, 20)

## This function simulates a phylogenetic tree using a birth-death model with extinction = 0. The tree is rescaled to unit height. It then simulates a trait variance-covariance matrix by simulating eigenvalues from an exponential distribution with lambda = 1/100. These are then used to simulate under multivariate Brownian motion on the phylogeny. Note that simulating the data for a large number of tips and traits is computationally intensive (requires a ntraits x ntips covariance matrix).
## Skipped running and just load saved simulations
#simdat <- foreach(i=1:100) %dopar% sim.tree.pcs.mv(ntips, ntraits, sig2dist=rexp, lambda=1/100)
#res <- foreach(i=1:100) %dopar% fitPCs(simdat[[i]], trait.seq, models=c("BM", "OUfixedRoot", "EB"))
#res <- do.call(rbind, res)
#save(res, file=paste("sims","mv", ntips, ntraits, ".rds", sep="_"))
load(paste("sims","mv", ntips, ntraits, ".rds", sep="_"))
par(mfrow=c(3,3), mar=c(4, 4.1, 2.1, 0.1))
for(i in 2:ncol(res)){
  boxplot(res[,i]~trait, data=res, ylim=c(0,1))
}

## Create plot of contrasts through time, for these simulations, we are going to use independent traits rather than mvBM.
nsims = 50
bmdat <- lapply(1:nsims, function(x) sim.tree.pcs.ind(ntips, ntraits, sig2dist=function(x){0.25}))
oudat <- lapply(1:nsims, function(x) sim.tree.pcs.ind.ou(ntips, ntraits, alpha=2, sig2=1))
ebdat <- lapply(1:nsims, function(x) sim.tree.pcs.ind.eb(ntips, ntraits, a=log(0.02), sig2=1))

dat <- list(bm=bmdat, ou=oudat, eb=ebdat)
models=c("BM", "OU", "EB")
dispdat <- lapply(dat, get.dtt)
dispdat[[1]]$simmodel = "BM"
dispdat[[2]]$simmodel = "OU"
dispdat[[3]]$simmodel = "EB"
dispdat <- do.call(rbind, dispdat)

bmcontrasts <- get.contrasts(bmdat)
bmcontdat <- do.call(rbind, bmcontrasts)
bmcontmelt <- melt(bmcontdat, id=c("type", "times"))
bmcontmelt$value <- abs(bmcontmelt$value)
bmcontmelt$simmodel <- "BM"

oucontrasts <- get.contrasts(oudat)
oucontdat <- do.call(rbind, oucontrasts)
oucontmelt <- melt(oucontdat, id=c("type", "times"))
oucontmelt$value <- abs(oucontmelt$value)
oucontmelt$simmodel <- "OU"

ebcontrasts <- get.contrasts(ebdat)
ebcontdat <- do.call(rbind, ebcontrasts)
ebcontmelt <- melt(ebcontdat, id=c("type", "times"))
ebcontmelt$value <- abs(ebcontmelt$value)
ebcontmelt$simmodel <- "EB"

contmelt <- rbind(bmcontmelt, oucontmelt, ebcontmelt)

#save(dispdat, file="disparitydata.rds")
#save(contmelt, file="contrastsdata.rds")


## Simulate under OU, BM and EB uncorrelated, fit models and calculate AICweights and compare parameter estimates.
trait.seq <- 1:20
#registerDoMC(cores=3)
#simdat <- foreach(i=1:100) %dopar% sim.tree.pcs.ind.ou(ntips, ntraits, alpha=2, sig2=1)
#res <- foreach(i=1:100) %dopar% fitPCOU(simdat[[i]], trait.seq, models=c("BM", "OUfixedRoot", "EB"))

#aicws <- lapply(res, function(x) x$aicw.table)
#oures <- do.call(rbind, aicws)

#parsdf <- get.parsOU(res)
#OUmeanAICw <- get.meanAICws(res)

#save(OUmeanAICw, file="OUmeanAICw.rds")
#save(parsdf, file="OUsimParameterEstimates.rds")
#save(oures, file=paste("simsOU","ind", ntips, ntraits, ".rds", sep="_"))

load("OUmeanAICw.rds")
load("OUsimParameterEstimates.rds")
load("simsOU_ind_50_20_.rds")

#simdat <- foreach(i=1:100) %dopar% sim.tree.pcs.ind.eb(ntips, ntraits, a=log(0.02), sig2=1)
#ebres <- foreach(i=1:100) %dopar% fitPCs(simdat[[i]], trait.seq, models=c("BM", "OUfixedRoot", "EB"))
#ebres <- do.call(rbind, ebres)
#save(ebres, file=paste("simsEB","ind", ntips, ntraits, ".rds", sep="_"))
load(paste("simsEB","ind", ntips, ntraits, ".rds", sep="_"))

#simdat <- foreach(i=1:100) %dopar% sim.tree.pcs.ind(ntips, ntraits, sig2dist=function(q, x){x}, x=1)
#bmres <- foreach(i=1:100) %dopar% fitPCs(simdat[[i]], trait.seq, models=c("BM", "OUfixedRoot", "EB"))
#bmres <- do.call(rbind, bmres)
#save(bmres, file=paste("simsBM","ind", ntips, ntraits, ".rds", sep="_"))
load(paste("simsBM","ind", ntips, ntraits, ".rds", sep="_"))
