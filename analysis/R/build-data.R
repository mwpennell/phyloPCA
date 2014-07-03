## Little script for building data
## This is only intended to be run once.
## Purpose of the script is to have a record of how datasets were modified

library(geiger)
## Felidae
fel.phy <- read.nexus("datasets/1741-7007-10-12-s5.nex")[[1]]
fel.dat <- read.csv("datasets/felidae_data.csv", header=TRUE, sep="\t", row.names=1)
fel.dat <- fel.dat[,c(1:7)]
fel.dat <- log(fel.dat)

## match tips
fel <- treedata(fel.phy, fel.dat)

## write to file
saveRDS(fel, "datasets/felidae.rds")


## Cyprinodon
cyp.phy <- read.nexus("datasets/cypri_underscore.nex")
cyp.phy$tip.label <- gsub("'", "", cyp.phy$tip.label)
cyp.dat <- read.csv("datasets/cypri_dataset.csv", row.names=1)
rownames(cyp.dat) <- gsub(" ", "_", rownames(cyp.dat))

## Remove outgroup
cyp.gen <- sapply(cyp.phy$tip.label, function(x) {strsplit(x, split="_")[[1]][1]})
cyp.phy <- drop.tip(cyp.phy, tip=cyp.phy$tip.label[which(cyp.gen != "Cyprinodon")])

## fix a couple of names in dataset
fix.cyp.names <- function(x){
    if (x == "Cyprinodon_sp._'durophage'_San_Salvador_Island")
        x <- "Cyprinodon_bozo"

    if (x == "Cyprinodon_sp._'scale-eater'_San_Salvador_Island")
        x <- "Cyprinodon_bulldog"

    if (x == "Cyprinodon_sp._'detritivore'_San_Salvador_Island")
        x <- "Cyprinodon_normal_Crescent_Pond"

    x
}

rownames(cyp.dat) <- sapply(rownames(cyp.dat), function(x) fix.cyp.names(x))

cyp <- treedata(cyp.phy, cyp.dat)

## write to file
saveRDS(cyp, "datasets/cyprinodon.rds")

