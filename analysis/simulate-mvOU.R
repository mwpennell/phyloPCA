## Here simulations for the mvOU case using mvSLOUCH:

source("R/simulate-helper.R")

## Parameters:
registerDoMC(cores=2)
ntips <- 50
ntraits <- 20
alpha <- 2
sig2 <- 1
nsims <- 50

## ## Sims using a diag R matrix and diag A matrix.
## Results should be equivalent to the ones simulated under independent traits (using same parameters here).
R <- diag(sig2, nrow = ntraits)
A <- diag(alpha, nrow = ntraits)
simou.uncor <- foreach(i=1:nsims) %dopar% sim.tree.pcs.mv.ou(ntips, ntraits, alpha=A, sig2=R)
res.ou.uncor <- foreach(i=1:nsims) %dopar% fitPCs(simou.uncor[[i]], c(1:20), models=c("BM", "OUfixedRoot", "EB"))
res.ou.uncor <- do.call(rbind, res.ou.uncor)
## Using 'diag' flag to remember is a A and R mat diag sims.
## Results from this sims should be equal to the ones with independent OU.
saveRDS(res.ou.uncor, "output/sim-res/mv.ou-uncor-diag.rds")

## ## Using a Posdef() R matrix:
R1 <- Posdef(20, ev=rexp(20, 1/sig2) )
A1 <- diag(alpha, nrow = ntraits)
simdat <- foreach(i=1:nsims) %dopar% sim.tree.pcs.mv.ou(ntips, ntraits, alpha=A1, sig2=R1)
res <- foreach(i=1:nsims) %dopar% fitPCs(simdat[[i]], c(1:20), models=c("BM", "OUfixedRoot", "EB"))
res <- do.call(rbind, res)
saveRDS(res, "output/sim-res/mv.ou-cor-sig.rds")

## ## Using a Posdef() A and R matrix:
## A has diag = rexp(1/2), then mean = 2.
R2 <- Posdef(20, ev=rexp(20, 1/sig2) )
A2 <- Posdef(20, ev=rexp(20, 1/alpha) )
simdat <- foreach(i=1:nsims) %dopar% sim.tree.pcs.mv.ou(ntips, ntraits, alpha=A2, sig2=R2)
res <- foreach(i=1:nsims) %dopar% fitPCs(simdat[[i]], c(1:20), models=c("BM", "OUfixedRoot", "EB"))
res <- do.call(rbind, res)
saveRDS(res, "output/sim-res/mv.ou-cor-sig-alpha.rds")
