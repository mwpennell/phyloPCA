library(geiger)
## Felidae
fel.phy <- read.nexus("datasets/1741-7007-10-12-s5.nex")[[1]]

## rescale tree to unit height
fel.phy$edge.length <- fel.phy$edge.length/(max(branching.times(fel.phy)))

fel.dat <- read.csv("datasets/felidae_data.csv", header=TRUE, sep="\t", row.names=1)
fel.dat <- fel.dat[,c(1:7)]
fel.dat <- log(fel.dat)

## match tips
fel <- treedata(fel.phy, fel.dat)

## write to file
saveRDS(fel, "output/data/felidae.rds")
