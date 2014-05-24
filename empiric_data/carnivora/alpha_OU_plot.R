load("./results_fitted.RData")

library(phylolm)

## Checking parameter estimates for the Ou model:

alpha.raw <- vector()
alpha.pc <- vector()
alpha.ppc <- vector()

for(i in 1:15){
    alpha.raw[i] <- rawfits[[2]][[i]]$optpar
    alpha.pc[i] <- pcfits[[2]][[i]]$optpar
    alpha.ppc[i] <- ppcfits[[2]][[i]]$optpar
}

## Even if OU is the best selected model for all pcs (independent of transformations) we expect the values of the
## alpha to be lower in the first set of pcs and decrease towards the last pcs.
## This trend should not be observed in the ppcs.

par(mfrow = c(2,2))
plot(1:15, alpha.raw, ylab = "alpha", xlab = "Climate", main = "raw")
plot(1:15, alpha.pc, ylab = "alpha", xlab = "PCs", main = "pc")
plot(1:15, alpha.ppc, ylab = "alpha", xlab = "PCs", main = "ppc")

dev.copy2pdf()
