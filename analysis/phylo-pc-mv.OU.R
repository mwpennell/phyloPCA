## This is for plotting results and results analysis:

library(mvSLOUCH)
source("R/simulate-helper.R")
source("R/plotting-helper.R")

## Build colour palettes
pal.bl <- colorRampPalette(brewer.pal(9, "Blues"))
pal.gn <- colorRampPalette(brewer.pal(9, "Greens"))
cols.bl <- rev(pal.bl(20))
cols.gn <- rev(pal.gn(20))
cols.line <- "#d73027"
cols <- c(cols.bl, cols.gn, cols.line)
## subset cols for model support figures
cols.aic <- c(cols.bl[2], cols.gn[2], cols.line)

## ## Model support for simulated data
## Note that this results have 10 instead of 100 sims.
ou.cor <- readRDS("output/sim-res/mv.ou-cor.rds")

## Build matrix for plotting

## As an overall summary of model support, we calculate the difference in AIC weight between an OU model and an EB model.
ou.cor.df <- build.sim.data.step(ou.cor)

## Figure 1 -- Model support across PC and pPC axes
jpeg(file = "ou.cor.df.jpeg", quality = 90)
fig.aicw(ou.cor.df, cols.aic)
dev.off()

## ## Effect of dimensionality
## We conducted an additional set of simulations to investigate the effect of matrix dimensionality on the patterns we found in the other simulations

rank.sim <- readRDS("output/sim-res/rankslopes_mvOU.rds")

## What the look of the AIC.w for each model depending on the domensionality
##     of the OU dataset?
class(rank.sim[[1]])
dim(rank.sim[[1]])
names(rank.sim[[1]])
summary(rank.sim[[1]])
ss <- rank.sim[[1]]
head(ss)
length(rank.sim[[1]][[1]])
jpeg(file = "ou.cor.df.jpeg", quality = 90)
fig.aicw(ou.cor.df, cols.aic)
dev.off()

## Figure 2 -- "Onion" plots showing slope of node height test as a function of the proportion of variance explained by leading eigenvector
jpeg(file = "rankslope.onion.mvOU.jpeg", quality = 90, width = 960)
fig.rankslopes(rank.sim$rankslopes, rank.sim$exp.val, cols)
dev.off()

## ## Uncorrelated data:
## Uncorrelated data is the R matrix with 0 as diagonals.
ou.uncor <- readRDS("output/sim-res/mv.ou-uncor.rds")
ou.uncor.df <- build.sim.data.step(ou.uncor)
jpeg(file = "ou.uncor.df.jpeg", quality = 90)
fig.aicw(ou.uncor.df, cols.aic)
dev.off()

## ## Effect of dimensionality
rank.uncor <- readRDS("output/sim-res/rankslopes.uncor_mvOU.rds")
jpeg(file = "rankslope.onion.uncor.mvOU.jpeg", quality = 90, width = 960)
fig.rankslopes(rank.uncor$rankslopes, rank.uncor$exp.val, cols)
dev.off()
