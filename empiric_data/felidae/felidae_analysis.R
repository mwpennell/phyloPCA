library(geiger)
library(phytools)
library(phylolm)

## Bininda-Edmunds 2012 BMC Felidae supertree.
## Getting only the first dated phylo in the file that has the best estimate for branch lengths.
phy <- read.nexus("./1741-7007-10-12-s5.nex")[[1]]

## Get data:
data <- read.csv("./data.csv", header = TRUE, sep = "\t")
rownames(data) <- data[,1]

## Drop from the tree the species we do not have data for:
phy.fel <- drop.tip(phy, tip = phy$tip.label[which(!phy$tip.label %in% data[,1])])

## Analysis:
mm <- match(phy.fel$tip.label, data[,1])
dt <- log(data[mm,2:8])

pc <- princomp(dt)
ppc <- phyl.pca(phy.fel,dt)

## Fitting the models:
models = c("BM", "OUrandomRoot", "EB")
rawfits <- list()
pcfits <- list()
ppcfits <- list()
for(i in 1:length(models)){
    rawfits[[i]] <- lapply(1:dim(dt)[2], function(x) phylolm(dt[,x]~1, phy=phy.fel, model=models[i]))
    pcfits[[i]] <-  lapply(1:dim(dt)[2], function(x) phylolm(pc$scores[,x]~1, phy=phy.fel, model=models[i]))
    ppcfits[[i]] <-  lapply(1:dim(dt)[2], function(x) phylolm(ppc$S[,x]~1, phy=phy.fel, model=models[i]))
}

## Get the AICw:
bm.aicw <- lapply(1:dim(dt)[2], function(x) aicw(c(rawfits[[1]][[x]]$aic, pcfits[[1]][[x]]$aic, ppcfits[[1]][[x]]$aic))$w)
ou.aicw <- lapply(1:dim(dt)[2], function(x) aicw(c(rawfits[[2]][[x]]$aic, pcfits[[2]][[x]]$aic, ppcfits[[2]][[x]]$aic))$w)
eb.aicw <- lapply(1:dim(dt)[2], function(x) aicw(c(rawfits[[3]][[x]]$aic, pcfits[[3]][[x]]$aic, ppcfits[[3]][[x]]$aic))$w)

## Make result tables:
bm.table <- do.call(rbind, bm.aicw)
colnames(bm.table) <- c("raw","pc","ppc")
ou.table <- do.call(rbind, ou.aicw)
colnames(ou.table) <- c("raw","pc","ppc")
eb.table <- do.call(rbind, eb.aicw)
colnames(eb.table) <- c("raw","pc","ppc")

## Plots:
pdf("Felidae_pcs.pdf")
par(mfrow = c(3,3))
for(i in 1:3){
    plot(1:dim(dt)[2], bm.table[,i], main = paste(models[1],colnames(bm.table)[i],sep="_")
         , ylim = c(0.0,1.0), ylab = "AICw", xlab = "PCs")
    plot(1:dim(dt)[2], ou.table[,i], main = paste(models[2],colnames(ou.table)[i],sep="_")
         , ylim = c(0.0,1.0), ylab = "AICw", xlab = "PCs")
    plot(1:dim(dt)[2], eb.table[,i], main = paste(models[3],colnames(eb.table)[i],sep="_")
         , ylim = c(0.0,1.0), ylab = "AICw", xlab = "PCs")
}
dev.off()
