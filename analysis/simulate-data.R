## Script for generating and analyzing simulated datasets

## source the helper file
## contains functions for performing the analysis
source("R/simulate-helper.R")


## Set the number of cores for parallel analysis
registerDoMC(cores=3)

## Specify the number of tips, traits and a vector specifying which traits to fit models to.
ntips <- 50
ntraits <- 20

## Simulate mv Brownian motion (correlated characters)

## This function simulates a phylogenetic tree using a birth-death model with extinction = 0. The tree is rescaled to unit height. It then simulates a trait variance-covariance matrix by simulating eigenvalues from an exponential distribution with lambda = 1/100. These are then used to simulate under multivariate Brownian motion on the phylogeny. Note that simulating the data for a large number of tips and traits is computationally intensive (requires a ntraits x ntips covariance matrix).

simdat <- foreach(i=1:100) %dopar% sim.tree.pcs.mv(ntips, ntraits,
                      sig2dist=rexp, lambda=1/100)
res <- foreach(i=1:100) %dopar% fitPCs(simdat[[i]], c(1:20),
                   models=c("BM", "OUfixedRoot", "EB"))
saveRDS(res, "output/sim-res/bm-cor.rds")



## Simulate uncorrelated data under BM, OU and EB and analyze using model selection, node-height test, and disparity through time.

## Simulate the data
nsims <- 100
bmdat <- lapply(1:nsims, function(x) sim.tree.pcs.ind(ntips, ntraits,sig2dist=function(x){0.25}))
oudat <- lapply(1:nsims, function(x) sim.tree.pcs.ind.ou(ntips, ntraits, alpha=2, sig2=1))
ebdat <- lapply(1:nsims, function(x) sim.tree.pcs.ind.eb(ntips, ntraits, a=log(0.02), sig2=1))

## Calculate the contrasts for each of the first 50 simulated datasets.
bmcont <- get.contrasts(bmdat[1:50], "BM")
oucont <- get.contrasts(oudat[1:50], "OU")
ebcont <- get.contrasts(ebdat[1:50], "EB")

## Combine all and save
contrasts <- rbind(bmcont, oucont, ebcont)
saveRDS(contrasts, file="output/sim-res/cont-height.rds")

## Calculate the disparity through time for each simulated dataset
bmdtt <- get.dtt(bmdat[1:50], "BM") 
oudtt <- get.dtt(oudat[1:50], "OU")
ebdtt <- get.dtt(ebdat[1:50], "EB")

dispdat <- do.call(rbind, list(bmdtt, oudtt, ebdtt))
saveRDS(dispdat, file="output/sim-res/disp-time.rds")



## Simulate under OU, BM and EB uncorrelated, fit models and calculate AICweights and compare parameter estimates.
registerDoMC(cores=3)
oufits <- foreach(i=1:100) %dopar% fitPCOU(oudat[[i]], 1:20, models=c("BM", "OUfixedRoot", "EB"))
ebfits <- foreach(i=1:100) %dopar% fitPCs(ebdat[[i]], 1:20, models=c("BM", "OUfixedRoot", "EB"))
bmfits <- foreach(i=1:100) %dopar% fitPCs(bmdat[[i]], 1:20, models=c("BM", "OUfixedRoot", "EB"))

## Get just the aicws for the OU results, filtering out the parameter estimates
ouaicws <- lapply(oufits, function(x) x$aicw.table)
oures <- do.call(rbind, ouaicws)

## Collect other AIC weights for EB and BM
ebres <- do.call(rbind, ebfits)
bmres <- do.call(rbind, bmfits)

## Save output
saveRDS(bmres, "output/sim-res/bm-uncor.rds")
saveRDS(oures, "output/sim-res/ou-uncor.rds")
saveRDS(ebres, "output/sim-res/eb-uncor.rds")



## For OU simulations, extract the parameter estimates
parsdf <- get.parsOU(oufits)
saveRDS(parsdf, "output/sim-res/OU-param.rds")




## Rank simulations
nsims <- 25
ntraits <- 20

## This sequence gives a reasonably smooth visualization across sets of eigenvalues where the proportion of variance explained by PC1 varies evenly between 1/20 and 1.
seqa <- c(-5, -1.5, -1, -0.75, -.5, -0.2, 0, 0.15, seq(0.3,2.1, 0.1), 2.25, 2.5, 5)
rank.seq <- exp(seqa)

## Look at the simulated values
#plot(seqa, sapply(rank.seq, function(x) 1/sum(ev^x)), xlim=c(-5, 5), ylim=c(0, 1))
#abline(h=sapply(rank.seq, function(x) 1/sum(ev^x)))

## Function for generating matrices of different rank
#base vector of eigenvalues that represent the mean relative values for a random set of 20 values drawn from an exponential distribution with rate 1/100
ev <-  rev(c(0.01, 0.03, 0.04, 0.06, 0.08, 0.10, 0.12, 0.14, 0.16, 0.19, 0.21, 0.24, 0.28, 0.32, 0.37, 0.42, 0.49, 0.58, 0.72, 1.00))
ev.rank <- function(n, p){
  ev^p
}


varEV1 <- sapply(1:length(rank.seq), function(x) 1/(sum(ev^rank.seq[x])))
rankdat <- lapply(rank.seq, function(y) lapply(1:nsims, function(x) sim.tree.pcs.mv(ntips, ntraits, sig2dist=ev.rank, p=y)))
rankcont <- lapply(1:length(rank.seq), function(x) get.contrasts(rankdat[[x]], x))
rankcont <- do.call(rbind, rankcont)
rankslopes <- ddply(rankcont, .(type, rep, simmodel, variable), summarize, slope=lmslope(value, times))
## store the values of varEV1 so that these can be used for plotting
ranks <- list(rankslopes=rankslopes, exp.val=varEV1)
saveRDS(ranks, file="output/sim-res/rankslopes.rds")








