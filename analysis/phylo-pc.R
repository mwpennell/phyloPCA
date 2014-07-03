## # Statistical and conceptual challenges in the comparative analysis of principal components

## Read in functions for analyzing and plotting results
source("R/simulate-helper.R")
source("R/plotting-helper.R")

## Build colour palettes
pal.bl <- colorRampPalette(brewer.pal(9, "Blues"))
pal.gn <- colorRampPalette(brewer.pal(9, "Greens"))
cols.bl <- rev(pal.bl(20))
cols.gn <- rev(pal.gn(20))
cols.line <- "#d73027"
cols <- c(cols.bl, cols.gn, cols.line)


## ## Model support for simulated data

## Note: All simulations can be run with the script `simulate-data.R`. As this takes a relatively long time (a couple of hours on a laptop), the results are stored as output.

## ### Multivariate BM (correlated traits)

## Multivariate BM was simulated on a pure-birth tree and rescaled to unit height. To obtain a variance-covariance matrix for the traits, we drew eigenvalues from an exponential distribution with rate 1/100.

## Read in data
bm.cor <- readRDS("output/bm-cor.rds")

## Build matrix for plotting

## As an overall summary of model support, we calculate the difference in AIC weight between an OU model and an EB model.
bm.cor.df <- build.sim.data.step(bm.cor)

## Figure 1 -- Model support across PC and pPC axes
fig.aicw(bm.cor.df)


## ### Multivariate BM (uncorrelated traits)

## Simulated in the same way as the correlated traits scenario (above) but setting all eigenvalues to be equal.

## Read in data
bm.uncor <- readRDS("output/bm-uncor.rds")

## Build matrix for plotting, using the difference in AICw between OU and EB as a summary statistic.
bm.uncor.df <- build.sim.data.step(bm.uncor)

## Supplementary Figure 1 -- Model support across PC and pPC axes
fig.aicw(bm.uncor.df)


## ### Multivariate OU (uncorrelated traits)

## All traits are evolved independently

## Read in data
ou.uncor <- readRDS("output/ou-uncor.rds")

## Build matrix for plotting, using the difference in AICw between OU and EB as a summary statistic.
ou.uncor.df <- build.sim.data.step(ou.uncor)

## Figure 3 -- Model support across PC and pPC axes
fig.aicw(ou.uncor.df)


## ### Multivariate EB (uncorrelated traits)

## All traits are evolved independently

## Read in data
eb.uncor <- readRDS("output/eb-uncor.rds")

## Build matrix for plotting, using the difference in AICw between OU and EB as a summary statistic
eb.uncor.df <- build.sim.data.step(eb.uncor)

## Supplementary Figure 2 -- Model support across PC and pPC axes
fig.aicw(eb.uncor.df)



## ## Node Height test

## Read in data containing size of contrasts and the height above the root that each contrast was calculated
cont <- readRDS("output/cont-height.rds")

## Figure 4 -- plot average slope of node-height state for all PC and pPC axes
fig.nh.2panel(cont, cols)



## ## Disparity through time

## Read in data containing disparity (sensu Harmon et al. 2003) and time at which it was calculated
disp <- readRDS("output/disp-time.rds")

## Figure 5 -- plot average loess function for all PC and pPC axes
fig.dtt.2panel(disp, cols)



## ## Parameter esitmation (OU model only)
## For the case of the uncorrelated OU model, we have obtained the estimated alpha parameter for each simulation
par.alpha <- readRDS("output/OU-param.rds")

## Only looking at the results from the phylogenetic PCA
ppca.alpha <- subset(par.alpha, variable == "ppc")

## Plot the estimated alpha values for each of the PCs
fig.alpha.est(ppca.alpha)


## ## Effect of dimensionality
## We conducted an additional set of simulations to investigate the effect of matrix dimensionality on the patterns we found in the other simulations

rank.sim <- readRDS("output/rankslopes.rds")

## Figure 2 -- "Onion" plots showing slope of node height test as a function of the proportion of variance explained by leading eigenvector
fig.rankslopes(rank.sim)




## ## Empirical studies

## ## Felidae
## Read in data
fel <- readRDS("datasets/felidae.rds")

## Compute principal components using the correlation approach as units of measurement differ between traits
fel.pca <- princomp(fel$data, cor=TRUE)

## Compute phylogenetic principal components again, using the corrleation approach
fel.ppca <- phyl.pca(fel$phy, fel$data, mode="corr")

## Prepare dataset for model fitting
fel.dat <- list(tree=fel$phy, raw=fel$data, pc=fel.pca$scores, ppc=fel.ppca$S,
                pcall=fel.pca, ppcall=fel.ppca)

## ### Fit models and compute AICw across all traits
fel.fit <- fitPCs(fel.dat, seq_len(ncol(fel$data)))

## Supplementary Figure 4 -- Model support across all PC and pPC axes
fel.df <- build.emp.data.step(fel.fit)

fig.aicw.empirical(fel.df)


## ### Node height test and Disparity through time
fel.cont <- get.contrasts(list(fel.dat), "contrats")
fel.disp <- get.dtt(list(fel.dat), "disparity")

## Supplementary Figure 5 -- Felidae contrasts and dtt
fig.nh.dtt(fel.cont, fel.disp)



## ## Cyprinodon
## Read in data
cyp <- readRDS("datasets/cyprinodon.rds")

## Compute principal components using the correlation approach as units of measurement differ between traits
cyp.pca <- princomp(cyp$data, cor=TRUE)

## Compute phylogenetic principal components again, using the corrleation approach
cyp.ppca <- phyl.pca(cyp$phy, cyp$data, mode="corr")

## Prepare dataset for model fitting
cyp.dat <- list(tree=cyp$phy, raw=cyp$data, pc=cyp.pca$scores, ppc=cyp.ppca$S,
                pcall=cyp.pca, ppcall=cyp.ppca)

## ### Fit models and compute AICw across all traits
cyp.fit <- fitPCs(cyp.dat, seq_len(ncol(cyp$data)))

## Supplementary Figure 6 -- Model support across all PC and pPC axes
cyp.df <- build.emp.data.step(cyp.fit)

fig.aicw.empirical(cyp.df)


## ### Node height test and Disparity through time
cyp.cont <- get.contrasts(list(cyp.dat), "contrats")
cyp.disp <- get.dtt(list(cyp.dat), "disparity")

## Supplementary Figure 7 -- Cyprinodon contrasts and dtt
fig.nh.dtt(cyp.cont, cyp.disp)



















