library(geiger)
library(phytools)
library(phylolm)

## Bininda-Edmunds 2012 BMC Felidae supertree.
## Getting only the first dated phylo in the file that has the best estimate for branch lengths.
phy <- read.nexus("../../datasets/1741-7007-10-12-s5.nex")[[1]]

## Get data:
data <- read.csv("../../datasets/data.csv", header = TRUE, sep = "\t")
rownames(data) <- data[,1]

## Drop from the tree the species we do not have data for:
phy.fel <- drop.tip(phy, tip = phy$tip.label[which(!phy$tip.label %in% data[,1])])

## Analysis:
mm <- match(phy.fel$tip.label, data[,1])
dt <- log(data[mm,2:8])

pc <- princomp(dt)
ppc <- phyl.pca(phy.fel,dt)

## Fitting the models:
## models = c("BM", "OUrandomRoot", "EB")
models <- c("BM", "OUfixedRoot", "EB")
rawfits <- list()
pcfits <- list()
ppcfits <- list()
for(i in 1:length(models)){
    rawfits[[i]] <- lapply(1:dim(dt)[2], function(x) phylolm(dt[,x]~1, phy=phy.fel, model=models[i]))
    pcfits[[i]] <-  lapply(1:dim(dt)[2], function(x) phylolm(pc$scores[,x]~1, phy=phy.fel, model=models[i]))
    ppcfits[[i]] <-  lapply(1:dim(dt)[2], function(x) phylolm(ppc$S[,x]~1, phy=phy.fel, model=models[i]))
}

## Get the AICw:
bm.aic <- lapply(1:dim(dt)[2], function(x) (c(rawfits[[1]][[x]]$aic, pcfits[[1]][[x]]$aic, ppcfits[[1]][[x]]$aic)))
ou.aic <- lapply(1:dim(dt)[2], function(x) (c(rawfits[[2]][[x]]$aic, pcfits[[2]][[x]]$aic, ppcfits[[2]][[x]]$aic)))
eb.aic <- lapply(1:dim(dt)[2], function(x) (c(rawfits[[3]][[x]]$aic, pcfits[[3]][[x]]$aic, ppcfits[[3]][[x]]$aic)))
all.aicw <- lapply(1:dim(dt)[2], function(x) lapply(1:3, function(y) aicw(c(bm.aic[[x]][y], ou.aic[[x]][y], eb.aic[[x]][y]))$w))


## Make result tables:
bm.table <- sapply(all.aicw, function(x) sapply(x, function(y) y[1]))
ou.table <- sapply(all.aicw, function(x) sapply(x, function(y) y[2]))
eb.table <- sapply(all.aicw, function(x) sapply(x, function(y) y[3]))
rownames(bm.table) <- rownames(ou.table) <- rownames(eb.table) <- c("raw","pc","ppc")

## Save results:
##save.image("felidae_results.RData")
