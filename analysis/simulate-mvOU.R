## Simulations with a diag alpha matrix of value 2 and a positive definite matrix for sigma like in the original simulations:

library(mvSLOUCH)
source("R/simulate-helper.R")

## Simulate correlated data:
registerDoMC(cores=2)
ntips <- 50
ntraits <- 20
alpha <- 2

## Using a Posdef() R matrix:
simdat <- foreach(i=1:50) %dopar% sim.tree.pcs.mv.ou(ntips, ntraits, alpha=alpha, cor=TRUE,
                      sig2dist=rexp, lambda=1/100)
res <- foreach(i=1:50) %dopar% fitPCs(simdat[[i]], c(1:20),
                   models=c("BM", "OUfixedRoot", "EB"))
res <- do.call(rbind, res)
saveRDS(res, "output/sim-res/mv.ou-cor.rds")

## Using a diag() R matrix:
simdat.uncor <- foreach(i=1:50) %dopar% sim.tree.pcs.mv.ou(ntips, ntraits, alpha=alpha, cor=FALSE,
                      sig2dist=rexp, lambda=1/100)
res.uncor <- foreach(i=1:50) %dopar% fitPCs(simdat.uncor[[i]], c(1:20),
                   models=c("BM", "OUfixedRoot", "EB"))
res.uncor <- do.call(rbind, res.uncor)
saveRDS(res.uncor, "output/sim-res/mv.ou-uncor.rds")

## Check for correlation among traits. Pairs plot:
panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...)
{
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    r <- abs(cor(x, y))
    txt <- format(c(r, 0.123456789), digits = digits)[1]
    txt <- paste0(prefix, txt)
    if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
    text(0.5, 0.5, txt, cex = cex.cor * r)
}
pdf(file = "pairs.uncor.mv.OU.pdf")
pairs(simdat.uncor[[1]]$raw[,1:10], lower.panel = panel.smooth, upper.panel = panel.cor)
dev.off()

###########################################################
## Rank simulations -- From 'simulate-data.R'
## Create data using a continuum of proportion of variance explained by the first PC axis.
## Lower values of prop explained by first axis means less correlation, higher values
##       means more correlation.

nsims <- 5
ntips <- 50
ntraits <- 20

seqa <- c(-5, -1.5, -1, -0.75, -.5, -0.2, 0, 0.15, seq(0.3,2.1, 0.1), 2.25, 2.5, 5)
rank.seq <- exp(seqa)
ev <-  rev(c(0.01, 0.03, 0.04, 0.06, 0.08, 0.10, 0.12, 0.14, 0.16, 0.19, 0.21, 0.24, 0.28, 0.32, 0.37, 0.42, 0.49, 0.58, 0.72, 1.00))
ev.rank <- function(n, p){
  ev^p
}

varEV1 <- sapply(1:length(rank.seq), function(x) 1/(sum(ev^rank.seq[x])))
rankdat <- lapply(rank.seq, function(y)
    lapply(1:nsims, function(x) sim.tree.pcs.mv.ou(ntips, ntraits, alpha=2, cor = "FALSE"
                                                 , sig2dist=ev.rank, p=y))
                  )
rankcont <- lapply(1:length(rank.seq), function(x) get.contrasts(rankdat[[x]], x))
rankcont <- do.call(rbind, rankcont)
rankslopes <- ddply(rankcont, .(type, rep, simmodel, variable), summarize, slope=lmslope(value, times))
ranks.uncor <- list(rankslopes=rankslopes, exp.val=varEV1)
saveRDS(ranks, file="output/sim-res/rankslopes.uncor_mvOU.rds")
