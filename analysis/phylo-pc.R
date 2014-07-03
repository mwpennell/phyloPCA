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
## subset cols for model support figures
cols.aic <- c(cols.bl[2], cols.gn[2], cols.line)


## ## Model support for simulated data

## Note: All simulations can be run with the script `simulate-data.R`. As this takes a relatively long time (a couple of hours on a laptop), the results are stored as output.

## ### Multivariate BM (correlated traits)

## Multivariate BM was simulated on a pure-birth tree and rescaled to unit height. To obtain a variance-covariance matrix for the traits, we drew eigenvalues from an exponential distribution with rate 1/100.

## Read in data
bm.cor <- readRDS("output/sim-res/bm-cor.rds")

## Build matrix for plotting

## As an overall summary of model support, we calculate the difference in AIC weight between an OU model and an EB model.
bm.cor.df <- build.sim.data.step(bm.cor)

## Figure 1 -- Model support across PC and pPC axes
fig.aicw(bm.cor.df, cols.aic)


## ### Multivariate BM (uncorrelated traits)

## Simulated in the same way as the correlated traits scenario (above) but setting all eigenvalues to be equal.

## Read in data
bm.uncor <- readRDS("output/sim-res/bm-uncor.rds")

## Build matrix for plotting, using the difference in AICw between OU and EB as a summary statistic.
bm.uncor.df <- build.sim.data.step(bm.uncor)

## Supplementary Figure 1 -- Model support across PC and pPC axes
fig.aicw(bm.uncor.df, cols.aic)


## ### Multivariate OU (uncorrelated traits)

## All traits are evolved independently

## Read in data
ou.uncor <- readRDS("output/sim-res/ou-uncor.rds")

## Build matrix for plotting, using the difference in AICw between OU and EB as a summary statistic.
ou.uncor.df <- build.sim.data.step(ou.uncor)

## Figure 3 -- Model support across PC and pPC axes
fig.aicw(ou.uncor.df, cols.aic)


## ### Multivariate EB (uncorrelated traits)

## All traits are evolved independently

## Read in data
eb.uncor <- readRDS("output/sim-res/eb-uncor.rds")

## Build matrix for plotting, using the difference in AICw between OU and EB as a summary statistic
eb.uncor.df <- build.sim.data.step(eb.uncor)

## Supplementary Figure 2 -- Model support across PC and pPC axes
fig.aicw(eb.uncor.df, cols.aic)



## ## Node Height test

## Read in data containing size of contrasts and the height above the root that each contrast was calculated
cont <- readRDS("output/sim-res/cont-height.rds")

## Figure 4 -- plot average slope of node-height state for all PC and pPC axes
fig.nh.2panel(cont, cols)



## ## Disparity through time

## Read in data containing disparity (sensu Harmon et al. 2003) and time at which it was calculated
disp <- readRDS("output/sim-res/disp-time.rds")

## Figure 5 -- plot average loess function for all PC and pPC axes
fig.dtt.2panel(disp, cols)



## ## Parameter esitmation (OU model only)
## For the case of the uncorrelated OU model, we have obtained the estimated alpha parameter for each simulation
par.ou <- readRDS("output/sim-res/OU-param.rds")

## Only looking at the results from the phylogenetic PCA
ppca.ou <- subset(par.ou, variable == "ppc")

## Plot the estimated alpha values for each of the PCs
fig.alpha.est(ppca.ou, col.pt=cols.gn[2], col.line=cols.line)


## ## Effect of dimensionality
## We conducted an additional set of simulations to investigate the effect of matrix dimensionality on the patterns we found in the other simulations

rank.sim <- readRDS("output/sim-res/rankslopes.rds")

## Figure 2 -- "Onion" plots showing slope of node height test as a function of the proportion of variance explained by leading eigenvector
fig.rankslopes(rank.sim$rankslopes, rank.sim$exp.val, cols)




## ## Empirical studies

## ## Felidae
## Read in data
fel <- readRDS("output/data/felidae.rds")

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

fig.aicw.empirical(fel.df, cols=c(cols.bl[2], cols.gn[2], cols.line))


## Felidae dataset is highly correlated: PC1 explains 96.9% of the variance with PCA
round((fel.pca$sdev^2)[1] / sum(fel.pca$sdev^2), digits=3)
## and 93.7% of the variance with pPCA
round(diag(fel.ppca$Eval)[1] / sum(diag(fel.ppca$Eval)), digits=3)



## ### Node height test and Disparity through time
## add dummy variables to use functions built for simulation
nsims <- 1
ntraits <- ncol(fel$data)
fel.cont <- get.contrasts(list(fel.dat), "contrasts")
fel.disp <- get.dtt(list(fel.dat), "disparity")

## Prepare data for plotting
fel.cont <- fel.cont[,which(colnames(fel.cont) != "rep")]
fel.all <- rbind(fel.cont, fel.disp)
cols.fel <- c(cols.bl[1:ntraits], cols.gn[1:ntraits], cols.line)

## Supplementary Figure 5 -- Felidae contrasts and dtt
fig.nh.dtt.emp(fel.all, cols.fel)



## ## Cyprinodon
## Read in data
cyp <- readRDS("output/data/cyprinodon.rds")

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

fig.aicw.empirical(cyp.df, cols=c(cols.bl[2], cols.gn[2], cols.line))


## ### Node height test and Disparity through time
ntraits <- ncol(cyp$data)
cyp.cont <- get.contrasts(list(cyp.dat), "contrasts")
cyp.disp <- get.dtt(list(cyp.dat), "disparity")

## Prepare data for plotting
cyp.cont <- cyp.cont[,which(colnames(cyp.cont) != "rep")]
cyp.all <- rbind(cyp.cont, cyp.disp)
cols.cyp <- c(cols.bl[1:ntraits], cols.gn[1:ntraits], cols.line)

## Supplementary Figure 5 -- Cypridon contrasts and dtt
fig.nh.dtt.emp(cyp.all, cols.cyp)

## Cyprinodon dataset is much less correlated: PC1 explains 39.8% of the variance with PCA
round((cyp.pca$sdev^2)[1] / sum(cyp.pca$sdev^2), digits=3)
## and 32.0% of the variance with pPCA
round(diag(cyp.ppca$Eval)[1] / sum(diag(cyp.ppca$Eval)), digits=3)





















