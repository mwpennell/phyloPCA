# Simulations exploring the effects of different types of PCA on model-based phylogenetic inference.

Load packages

```r
rm(list = ls(all = TRUE))
require(phytools)
```

```
## Loading required package: phytools
## Loading required package: ape
## Loading required package: maps
## Loading required package: rgl
```

```r
require(geiger)
```

```
## Loading required package: geiger
```

```r
require(phylolm)
```

```
## Loading required package: phylolm
```

```r
require(foreach)
```

```
## Loading required package: foreach
```

```r
require(doMC)
```

```
## Loading required package: doMC
## Loading required package: iterators
## Loading required package: parallel
```

```r
require(MASS)
```

```
## Loading required package: MASS
```

```r
require(nlme)
```

```
## Loading required package: nlme
```

```r
require(reshape2)
```

```
## Loading required package: reshape2
```

```r
require(ggplot2)
```

```
## Loading required package: ggplot2
```

```r
require(plyr)
```

```
## Loading required package: plyr
```


Read in functions

```r
source("./R/analysis-helper.R")
```


Set the number of cores for parallel analysis

```r
registerDoMC(cores = 3)
```


Specify the number of tips, traits and a vector specifying which traits to fit models to.

```r
ntips <- 50
ntraits <- 20
```


Select which of the traits to study

```r
trait.seq = c(1, 2, 5, 10, 15, 20)
```


This function simulates a phylogenetic tree using a birth-death model with extinction = 0. The tree is rescaled to unit height. It then simulates a trait variance-covariance matrix by simulating eigenvalues from an exponential distribution with lambda = 1/100. These are then used to simulate under multivariate Brownian motion on the phylogeny. Note that simulating the data for a large number of tips and traits is computationally intensive (requires a ntraits x ntips covariance matrix).
Skipped running and just load saved simulations

```r
# simdat <- foreach(i=1:100) %dopar% sim.tree.pcs.mv(ntips, ntraits,
# sig2dist=rexp, lambda=1/100) res <- foreach(i=1:100) %dopar%
# fitPCs(simdat[[i]], trait.seq, models=c('BM', 'OUfixedRoot', 'EB')) res <-
# do.call(rbind, res) save(res, file=paste('./output/sims','mv', ntips,
# ntraits, '.rds', sep='_'))

load(paste("./output/sims", "mv", ntips, ntraits, ".rds", sep = "_"))

bm.cor <- build.sim.data.table(res, trait.set = trait.seq)
```

```
## Using trait as id variables
```


Figure 1- Distribution of Akaike weights for the Brownian motion (top row), Ornstein–Uhlenbeck (middle row) or Early–Burst model (bottom row).
Data simulated under correlated multivariate Brownian Motion

```r
fig.box.aicw(bm.cor)
```

![plot of chunk unnamed-chunk-7](figure/unnamed-chunk-7.png) 




 # Simulate uncorrelated data under BM, OU and EB and analyze using model selection, node-height test, and disparity through time.

Simulate the data

```r
nsims = 100
# bmdat <- lapply(1:nsims, function(x) sim.tree.pcs.ind(ntips, ntraits,
# sig2dist=function(x){0.25})) oudat <- lapply(1:nsims, function(x)
# sim.tree.pcs.ind.ou(ntips, ntraits, alpha=2, sig2=1)) ebdat <-
# lapply(1:nsims, function(x) sim.tree.pcs.ind.eb(ntips, ntraits,
# a=log(0.02), sig2=1))
```


Calculate the contrasts for each of the first 50 simulated datasets.

```r
# bmcont <- get.contrasts(bmdat[1:50], 'BM') oucont <-
# get.contrasts(oudat[1:50], 'OU') ebcont <- get.contrasts(ebdat[1:50],
# 'EB')
```


Combine all and save

```r
# contrasts <- rbind(bmcont, oucont, ebcont) save(contrasts,
# file='./output/contrastsdata.rds')
```


Calculate the disparity through time for each simulated dataset

```r
# bmdtt <- get.dtt(bmdat[1:50], 'BM') oudtt <- get.dtt(oudat[1:50], 'OU')
# ebdtt <- get.dtt(ebdat[1:50], 'EB')

# dispdat <- do.call(rbind, list(bmdtt, oudtt, ebdtt)) save(dispdat,
# file='./output/disparitydata.rds')
```


Load back in saved output for contrasts and disparity

```r
load("./output/contrastsdata.rds")
load("./output/disparitydata.rds")
```


Figure 3- Relationship between the average phylogenetic independent contrasts and the height of the node across 100 datasets simulated under either

```r
## a BM (top row), OU (middle row) or EB (bottom row) model of evolution.
fig.nh.3panel(contrasts)
```

![plot of chunk unnamed-chunk-13](figure/unnamed-chunk-13.png) 


Figure 4- Relative disparity through time for the same simulated datasets as in Figure 3.

```r
fig.dtt.3panel(dispdat)
```

```
## Error: cannot allocate vector of size 128.0 Mb
```


Simulate under OU, BM and EB uncorrelated, fit models and calculate AICweights and compare parameter estimates.

```r
registerDoMC(cores = 3)
# oufits <- foreach(i=1:100) %dopar% fitPCOU(oudat[[i]], 1:20,
# models=c('BM', 'OUfixedRoot', 'EB')) ebfits <- foreach(i=1:100) %dopar%
# fitPCs(ebdat[[i]], trait.seq, models=c('BM', 'OUfixedRoot', 'EB')) bmfits
# <- foreach(i=1:100) %dopar% fitPCs(bmdat[[i]], trait.seq, models=c('BM',
# 'OUfixedRoot', 'EB'))
```


Get just the aicws for the OU results, filtering out the parameter estimates

```r
# ouaicws <- lapply(oufits, function(x) x$aicw.table) oures <-
# do.call(rbind, ouaicws)
```

Collect other AIC weights for EB and BM

```r
# ebres <- do.call(rbind, ebfits) bmres <- do.call(rbind, bmfits)
```


Save output

```r
# save(oures, file=paste('./output/simsOU','ind', ntips, ntraits, '.rds',
# sep='_')) save(ebres, file=paste('./output/simsEB','ind', ntips, ntraits,
# '.rds', sep='_')) save(bmres, file=paste('./output/simsBM','ind', ntips,
# ntraits, '.rds', sep='_'))
```


Load fits back into the workspace

```r
load(paste("./output/simsOU", "ind", ntips, ntraits, ".rds", sep = "_"))
load(paste("./output/simsBM", "ind", ntips, ntraits, ".rds", sep = "_"))
load(paste("./output/simsEB", "ind", ntips, ntraits, ".rds", sep = "_"))
```


Get the parameter estimates from the OU simulations, as well as the mean AIC weights.

```r
# parsdf <- get.parsOU(oufits) OUmeanAICw <- get.meanAICws(oufits)

# save(OUmeanAICw, file='./output/OUmeanAICw.rds') save(parsdf,
# file='./output/OUsimParameterEstimates.rds')

load("./output/OUmeanAICw.rds")
load("./output/OUsimParameterEstimates.rds")
```


Multivariate BM

```r
bm.ind <- build.sim.data.table(bmres, trait.seq)
```

```
## Using trait as id variables
```


Boxplots for AIC weights for models for uncorrelated mvBM

```r
fig.box.aicw(bm.ind)
```

![plot of chunk unnamed-chunk-22](figure/unnamed-chunk-22.png) 


Boxplots for AIC weights for models for uncorrelated mvOU

```r
ou.ind <- build.sim.data.table(oures, trait.seq)
```

```
## Using trait as id variables
```

```r
fig.box.aicw(ou.ind)
```

![plot of chunk unnamed-chunk-23](figure/unnamed-chunk-23.png) 


Boxplots for AIC weights for models for uncorrelated mvEB

```r
eb.ind <- build.sim.data.table(ebres, trait.seq)
```

```
## Using trait as id variables
```

```r
fig.box.aicw(eb.ind)
```

![plot of chunk unnamed-chunk-24](figure/unnamed-chunk-24.png) 


do a bit of processing

```r
parsdf$trait <- as.character(parsdf$trait)
parsdf$trait <- sapply(parsdf$trait, function(x) sub(".", "", x, fixed = TRUE))
parsdf$trait <- factor(parsdf$trait)
pars.tmp <- melt(parsdf)
```

```
## Using trait, variable as id variables
```

```r
colnames(pars.tmp) <- colnames(OUmeanAICw)
```


combine the dataframes

```r
est <- rep("AICw", nrow(OUmeanAICw))
OUmeanAICw <- cbind(OUmeanAICw, est)
est <- rep("Par", nrow(pars.tmp))
pars.tmp <- cbind(pars.tmp, est)
oudf <- rbind(OUmeanAICw, pars.tmp)
```


prune down dataset to raw and phylo pca only

```r
oudf <- subset(oudf, type %in% c("raw", "ppc"))
oudf$type <- factor(oudf$type)
oudf$type <- factor(oudf$type, levels = c("raw", "ppc"), labels = c("Original data", 
    "Phylogenetic PCA"))
```


prune down dataset to exclude only BM, OU, EB, alpha

```r
oudf <- subset(oudf, simmodel %in% c("BM", "EB", "OUfixedRoot", "alpha"))
oudf$simmodel <- factor(oudf$simmodel)
oudf$simmodel <- factor(oudf$simmodel, levels = c("BM", "OUfixedRoot", "EB", 
    "alpha"), labels = c("BM", "OU", "EB", "alpha"))

oudf$trait <- factor(oudf$trait, levels = unique(oudf$trait), labels = c(1:20))
```



Figure showing detailed model fit results and parameter estimates when simulated under uncorrelated mvOU and fit to each of the different models.

```r
fig.model.support.alpha(oudf)
```

```
## Loading required package: gridExtra
## Loading required package: grid
```

```
## Warning: Removed 1 rows containing non-finite values (stat_boxplot).
```

![plot of chunk unnamed-chunk-29](figure/unnamed-chunk-29.png) 


## Empirical analyses of morphological dataset for Felids

Reference for the Felids phylogeny, trees available in on-line supplementary file:
Nyakatura, K., and O. R. Bininda-Emonds. 2012. Updating the evolutionary history of Carnivora (Mammalia): A new species-level supertree complete with divergence time estimates. BMC biology 10:12.

Get only the first dated phylogeny in the file. Tree has best estimate for branch lengths.

```r
phy <- read.nexus("./datasets/1741-7007-10-12-s5.nex")[[1]]
```


Morphological data manually compiled from:
Slater, G. J., and B. Van Valkenburgh. 2009. Allometry and performance: The evolution of skull form and function in felids. Journal of Evolutionary Biology 22:2278–2287.
Sakamoto, M., G. T. Lloyd, and M. J. Benton. 2010. Phylogenetically structured variance in felid bite force: The role of phylogeny in the evolution of biting performance. Journal of Evolutionary Biology 23:463–478.

Get dataset:

```r
data <- read.csv("./datasets/data.csv", header = TRUE, sep = "\t")
rownames(data) <- data[, 1]
```


Drop from the tree the species we do not have data for:

```r
phy.fel <- drop.tip(phy, tip = phy$tip.label[which(!phy$tip.label %in% data[, 
    1])])
phy.fel$edge.length <- phy.fel$edge.length/max(branching.times(phy.fel))
```


Model selection analysis:

```r
mm <- match(phy.fel$tip.label, data[, 1])
dt <- log(data[mm, 2:8])
felidPC <- princomp(dt, cor = TRUE)
felidPPC <- phyl.pca(phy.fel, dt, mode = "corr")

felidae.dat <- list(tree = phy.fel, raw = dt, pc = felidPC$scores, ppc = felidPPC$S, 
    pcall = felidPC, ppcall = felidPPC)
felidaeFit <- fitPCs(felidae.dat, 1:ncol(dt), models = c("BM", "OUfixedRoot", 
    "EB"))

felidae.ind <- build.sim.data.table(felidaeFit, 1:7)
```

```
## Using trait as id variables
```


Felid dataset is highly correlated. PC1 explains 96.9% of the variance under PCA 93.7 under pPCA.

```r
round((felidPC$sdev^2)[1]/sum(felidPC$sdev^2), digits = 3)
```

```
## Comp.1 
##  0.969
```

```r
round(diag(felidPPC$Eval)[1]/sum(diag(felidPPC$Eval)), digits = 3)
```

```
##   PC1 
## 0.937
```


Model support across datasets showing AIC weights:

```r
fig.emp.aicw(felidae.ind)
```

![plot of chunk unnamed-chunk-35](figure/unnamed-chunk-35.png) 


Analysis of contrasts and disparity:

```r
nsims = 1
ntraits = ncol(dt)
felidCont <- get.contrasts(list(felidae.dat), "contrasts")
felidDisp <- get.dtt(list(felidae.dat), "disparity")
```


Node height and disparity through time plots for the felid dataset.

```r
fig.felidae.contrasts(felidCont)
```

![plot of chunk unnamed-chunk-37](figure/unnamed-chunk-371.png) 

```r
fig.felidae.dtt(felidDisp)
```

```
## geom_smooth: method="auto" and size of largest group is <1000, so using loess. Use 'method = x' to change the smoothing method.
## geom_smooth: method="auto" and size of largest group is <1000, so using loess. Use 'method = x' to change the smoothing method.
## geom_smooth: method="auto" and size of largest group is <1000, so using loess. Use 'method = x' to change the smoothing method.
```

![plot of chunk unnamed-chunk-37](figure/unnamed-chunk-372.png) 


## Empirical analyses of morphological dataset for Cyprinodon.

Read the dataset and phylogeny.

```r
tr <- read.nexus("./datasets/fish_finaltree_underscore.nex")
dt <- read.csv("./datasets/fish_dataset.csv")

dt$species <- gsub(" ", "_", dt$species)
tr$tip.label <- gsub("'", "", tr$tip.label)
```


Keep only the data for Cyprinodon species.

```r
sp.split <- sapply(dt$species, FUN = function(x) strsplit(x, split = "_"))
cyp <- sapply(1:length(sp.split), FUN = function(x) sp.split[[x]][1])
id <- which(cyp == "Cyprinodon")
dt <- dt[id, ]
```


Some species names need correction.

```r
dt$species[which(dt$species == "Cyprinodon_sp._'durophage'_San_Salvador_Island")] <- "Cyprinodon_bozo"
dt$species[which(dt$species == "Cyprinodon_sp._'scale-eater'_San_Salvador_Island")] <- "Cyprinodon_bulldog"
dt$species[which(dt$species == "Cyprinodon_sp._'detritivore'_San_Salvador_Island")] <- "Cyprinodon_normal_Crescent_Pond"
```


Check if all species in the dataset are in the phylogeny.

```r
index <- which(tr$tip.label %in% dt$species)
sum(dt$species %in% tr$tip.label[index]) == dim(dt)[1]
```

```
## [1] TRUE
```


Keep only the species in the dataset in the tree (Drop the non-Cyprinodon).

```r
tr <- drop.tip(phy = tr, tip = tr$tip.label[-index])
```


Match dataset and tree. Keep only the morphological measurements.
Data is already size corrected.

```r
mm <- match(tr$tip.label, dt$species)
dt <- dt[mm, -c(17:25)]
rownames(dt) <- dt$species
dt <- dt[, -1]
```


Repeat analysis done with the Felid:

Model selection:

```r
cypriPC <- princomp(dt, cor = TRUE)
cypriPPC <- phyl.pca(tr, dt, mode = "corr")

cypri.dat <- list(tree = tr, raw = dt, pc = cypriPC$scores, ppc = cypriPPC$S, 
    pcall = cypriPC, ppcall = cypriPPC)
cypriFit <- fitPCs(cypri.dat, 1:ncol(dt), models = c("BM", "OUfixedRoot", "EB"))

cypri.ind <- build.sim.data.table(cypriFit, 1:15)
```

```
## Using trait as id variables
```


Show correlation. PC1 explains 37.8% and pPC1 32.3%.

```r
round((cypriPC$sdev^2)[1]/sum(cypriPC$sdev^2), digits = 3)
```

```
## Comp.1 
##  0.378
```

```r
round(diag(cypriPPC$Eval)[1]/sum(diag(cypriPPC$Eval)), digits = 3)
```

```
##   PC1 
## 0.323
```


Model support across datasets showing AIC weights:

```r
fig.emp.aicw(cypri.ind)
```

![plot of chunk unnamed-chunk-46](figure/unnamed-chunk-46.png) 


Analysis of contrasts and disparity:

```r
nsims = 1
ntraits = ncol(dt)
cypriCont <- get.contrasts(list(cypri.dat), "contrasts")
cypriDisp <- get.dtt(list(cypri.dat), "disparity")
```


Node height and disparity through time plots.

```r
fig.felidae.contrasts(cypriCont)
```

![plot of chunk unnamed-chunk-48](figure/unnamed-chunk-481.png) 

```r
fig.felidae.dtt(cypriDisp)
```

```
## geom_smooth: method="auto" and size of largest group is <1000, so using loess. Use 'method = x' to change the smoothing method.
## geom_smooth: method="auto" and size of largest group is <1000, so using loess. Use 'method = x' to change the smoothing method.
## geom_smooth: method="auto" and size of largest group is <1000, so using loess. Use 'method = x' to change the smoothing method.
```

![plot of chunk unnamed-chunk-48](figure/unnamed-chunk-482.png) 


## Need to put this part in the same order of the manuscript.
## Maybe use the same headings to organize this.
# Matrix Rank Simulations

```r
nsims = 25
ntraits = 20

# This sequence gives a reasonably smooth visualization across sets of
# eigenvalues where the proportion of variance explained by PC1 varies
# evenly between 1/20 and 1.
seqa <- c(-5, -1.5, -1, -0.75, -0.5, -0.2, 0, 0.15, seq(0.3, 2.1, 0.1), 2.25, 
    2.5, 5)
rank.seq <- exp(seqa)
# plot(seqa, sapply(rank.seq, function(x) 1/sum(ev^x)), xlim=c(-5, 5),
# ylim=c(0, 1)) abline(h=sapply(rank.seq, function(x) 1/sum(ev^x)))

varEV1 <- sapply(1:length(rank.seq), function(x) 1/(sum(ev^rank.seq[x])))
# rankdat <- lapply(rank.seq, function(y) lapply(1:nsims, function(x)
# sim.tree.pcs.mv(ntips, ntraits, sig2dist=ev.rank, p=y))) rankcont <-
# lapply(1:length(rank.seq), function(x) get.contrasts(rankdat[[x]], x))
# rankcont <- do.call(rbind, rankcont) rankslopes <- ddply(rankcont, .(type,
# rep, simmodel, variable), summarize, slope=lmslope(value, times))
# save(rankslopes, file='./output/rankslopes.rds')

load("./output/rankslopes.rds")

fig.rankslopes(rankslopes)
```

```
## geom_smooth: method="auto" and size of largest group is <1000, so using loess. Use 'method = x' to change the smoothing method.
## geom_smooth: method="auto" and size of largest group is <1000, so using loess. Use 'method = x' to change the smoothing method.
## geom_smooth: method="auto" and size of largest group is <1000, so using loess. Use 'method = x' to change the smoothing method.
```

![plot of chunk unnamed-chunk-49](figure/unnamed-chunk-49.png) 
