# Simulations exploring the effects of different types of PCA on model-based phylogenetic inference.

Load packages

```r
require(phytools)
require(geiger)
require(phylolm)
require(TreeSim)
require(foreach)
require(doMC)
require(MASS)
```


Read in functions

```r
source("PCAfunctions.R")
```


Set the number of cores for parallel analysis

```r
registerDoMC(cores = 4)
```


Specify the number of tips, traits and a vector specifying which traits to fit models to.

```r
ntips <- 50
ntraits <- 30
trait.seq = c(1, 2, 5, 10, 20, 30)
```


This function simulates a phylogenetic tree using a birth-death model with extinction = 0. The tree is rescaled to unit height. It then simulates a trait variance-covariance matrix by simulating eigenvalues from an exponential distribution with lambda = 1/100. These are then used to simulate under multivariate Brownian motion on the phylogeny. Note that simulating the data for a large number of tips and traits is computationally intensive (requires a ntraits x ntips covariance matrix).


```r
# Skipped running and just load saved simulations res <- foreach(i=1:100)
# %dopar% simandfitPCs(ntips, ntraits, trait.seq=trait.seq) save(res,
# file=paste('sims', ntips, ntraits,'.rds', sep='_'))
load(paste("sims", ntips, ntraits, ".rds", sep = "_"))
```


Rearrange the output to be usable for plotting

```r
res.all <- lapply(1:length(res[[1]]), function(x) sapply(1:length(res), function(y) res[[y]][[x]]))
res.all <- lapply(res.all, function(x) {
    x <- t(x)
    colnames(x) <- colnames(res[[1]][[1]])
    x
})
```


Create boxplots of Akaike weights for each transformation (raw, pca, and phylo  pca) for each of the traits in trait.seq, for each of the models (BM, OU, EB).

```r
par(mfrow = c(3, 3), mar = c(4, 4.1, 2.1, 0.1))
boxplot.prep <- function(x) {
    tmp <- sapply(res.all, function(y) y[, x])
    colnames(tmp) <- paste("", trait.seq, sep = "")
    tmp
}
bx <- lapply(1:ncol(res.all[[1]]), function(x) boxplot(boxplot.prep(x), ylab = "AICw", 
    xlab = "Trait", main = colnames(res.all[[1]])[x]))
```

![plot of chunk unnamed-chunk-7](figure/unnamed-chunk-7.png) 
