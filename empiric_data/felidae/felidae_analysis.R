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

## Taking a look at the principal component analysis:
pc$loadings
pc$scores
ppc$L
ppc$S

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
raw.aicw <- list()
pc.aicw <- list()
ppc.aicw <- list()
for(i in 1:length(models)){
    raw.aicw[[i]] <- aicw(sapply(1:dim(dt)[2], function(x) rawfits[[i]][[x]]$aic))$w
    pc.aicw[[i]] <- aicw(sapply(1:dim(dt)[2], function(x) pcfits[[i]][[x]]$aic))$w
    ppc.aicw[[i]] <- aicw(sapply(1:dim(dt)[2], function(x) ppcfits[[i]][[x]]$aic))$w
}

## Plots:
pdf("Felidae_pcs.pdf")
par(mfrow = c(3,3))
for(i in 1:length(models)){
    plot(1:dim(dt)[2], raw.aicw[[i]], main = paste(models[i],"raw",sep="_")
         , ylim = c(0.0,1.0), ylab = "AICw", xlab = "PCs")
    plot(1:dim(dt)[2], pc.aicw[[i]], main = paste(models[i],"pc",sep="_")
         , ylim = c(0.0,1.0), ylab = "AICw", xlab = "PCs")
    plot(1:dim(dt)[2], ppc.aicw[[i]], main = paste(models[i],"ppc",sep="_")
         , ylim = c(0.0,1.0), ylab = "AICw", xlab = "PCs")
}
dev.off()
