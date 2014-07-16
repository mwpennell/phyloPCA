## Anolis
library(geiger)
anoles.phy <- read.tree("datasets/GA_Anolis.tre")
anoles.dat <- read.csv("datasets/Anolis_Mahler2010.csv", row.names=1)

## Rescale phylogeny to unit height
anoles.phy$edge.length <- anoles.phy$edge.length / max(branching.times(anoles.phy))

## Match tree and data
anoles <- treedata(anoles.phy, anoles.dat)

## write to file
## create directory if necessary
tmp <- dir("output")
if (!"data" %in% tmp)
    dir.create("output/data")

saveRDS(anoles, "output/data/anoles.rds")
