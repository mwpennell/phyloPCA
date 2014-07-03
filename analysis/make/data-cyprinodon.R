## Cyprinodon
library(geiger)
cyp.phy <- read.nexus("datasets/cypri_underscore.nex")
cyp.phy$tip.label <- gsub("'", "", cyp.phy$tip.label)
cyp.dat <- read.csv("datasets/cypri_dataset.csv", row.names=1)
rownames(cyp.dat) <- gsub(" ", "_", rownames(cyp.dat))

## Remove outgroup
cyp.gen <- sapply(cyp.phy$tip.label, function(x) {strsplit(x, split="_")[[1]][1]})
cyp.phy <- drop.tip(cyp.phy, tip=cyp.phy$tip.label[which(cyp.gen != "Cyprinodon")])


## Rescale phylogeny to unit height
cyp.phy$edge.length <- cyp.phy$edge.length / max(branching.times(cyp.phy))

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

## remove PC socres
cyp.dat <- cyp.dat[,c(1:16)]

cyp <- treedata(cyp.phy, cyp.dat)

## write to file
## create directory if necessary
tmp <- dir("output")
if (!"data" %in% tmp)
    dir.create("output/data")

saveRDS(cyp, "output/data/cyprinodon.rds")
