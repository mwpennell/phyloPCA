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


## ### Multivariate OU
## 1) All traits are evolved independently:

## Read in data
ou.ind <- readRDS("output/sim-res/ou-res-ind.rds")

## Build matrix for plotting, using the difference in AICw between OU and EB as a summary statistic.
ou.ind.df <- build.sim.data.step(ou.ind)

## Figure 3 -- Model support across PC and pPC axes
#pdf("Figure3_uncor.pdf", width = 10.5, height = 7)
fig.aicw(ou.ind.df, cols.aic)
#dev.off()

## 2) All traits are evolved correlated:
## Same steps of previous block (1):
ou.cor <- readRDS("output/sim-res/ou-res-cor.rds")
ou.cor.df <- build.sim.data.step(ou.cor)
#pdf("Figure3_corr.pdf", width = 10.5, height = 7)
fig.aicw(ou.cor.df, cols.aic)
#dev.off()


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
## 1) For the uncorrelated OU traits:
disp <- readRDS("output/sim-res/disp-time.rds")
## Figure 5 -- plot average loess function for all PC and pPC axes

fig.dtt.2panel(disp, cols)

## ## Parameter esitmation (OU model only)
## 1) For the uncorrelated OU traits:
## For the case of the uncorrelated OU model, we have obtained the estimated alpha parameter for each simulation
par.ou.ind <- readRDS("output/sim-res/OU-param-ind.rds")

## Only looking at the results from the phylogenetic PCA
ppca.ou.ind <- subset(par.ou.ind, variable == "ppc")

## Plot the estimated alpha values for each of the PCs
fig.alpha.est(ppca.ou.ind, col.pt=cols.gn[2], col.line=cols.line)

## 1) For the correlated OU traits:
## Same steps as above:
par.ou.cor <- readRDS("output/sim-res/OU-param-cor.rds")
ppca.ou.cor <- subset(par.ou.cor, variable == "ppc")
fig.alpha.est(ppca.ou.cor, col.pt=cols.gn[2], col.line=cols.line)


## Effect of mixing different models - ACDC models with parameter drawn from normal distribution
## Replace with better plotting function
acdcres <- readRDS("output/sim-res/acdc.rds")

fig.acdc(acdcres, cols.aic)

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
fel.df <- build.emp.data.stack(fel.fit)

fig.aicw.empirical(fel.df, cols.aic[c(1,3,2)])

## Felidae dataset is highly correlated: PC1 explains 96.9% of the variance with PCA
round((fel.pca$sdev^2)[1] / sum(fel.pca$sdev^2), digits=3)
## and 93.7% of the variance with pPCA
round(diag(fel.ppca$Eval)[1] / sum(diag(fel.ppca$Eval)), digits=3)

## ### Node height test and Disparity through time
## add dummy variables to use functions built for simulation
nsims <- 1
ntraits <- ncol(fel$dat)
fel.cont <- get.contrasts(list(fel.dat), "contrasts")
fel.disp <- get.dtt(list(fel.dat), "disparity")

## Prepare data for plotting
fel.cont <- fel.cont[,which(colnames(fel.cont) != "rep")]
fel.all <- rbind(fel.cont, fel.disp)
cols.fel <- c(cols.bl[1:ntraits], cols.gn[1:ntraits], cols.line)

## Supplementary Figure 5 -- felidae contrasts and dtt
fig.nh.dtt.emp(fel.all, cols.fel)


## ## Anolis (Mahler et al. 2010)
anoles <- readRDS("output/data/anoles.rds")

## Compute principal components using the correlation approach as units of measurement differ between traits
anoles.pca <- princomp(anoles$dat, cor=TRUE)
anoles.ppca <- phyl.pca(anoles$phy, anoles$dat, mode="corr")

## Prepare dataset for model fitting
anoles.dat <- list(tree=anoles$phy, raw=anoles$dat, pc=anoles.pca$scores, ppc=anoles.ppca$S,
                   pcall=anoles.pca, ppcall=anoles.ppca)

## ### Fit models and compute AICw across all traits
anoles.fit <- fitPCs(anoles.dat, seq_len(ncol(anoles$dat)))

## Supplementary Figure 4 -- Model support across all PC and pPC axes
anoles.df <- build.emp.data.stack(anoles.fit)

fig.aicw.empirical(anoles.df, cols.aic[c(1,3,2)])

## Anolis dataset is also highly correlated: PC1 explains 92.6% of the variance with PCA
round((anoles.pca$sdev^2)[1] / sum(anoles.pca$sdev^2), digits=3)
## and 90.0% of the variance with pPCA
round(diag(anoles.ppca$Eval)[1] / sum(diag(anoles.ppca$Eval)), digits=3)

## ### Node height test and Disparity through time
## add dummy variables to use functions built for simulation
nsims <- 1
ntraits <- ncol(anoles$dat)
anoles.cont <- get.contrasts(list(anoles.dat), "contrasts")
anoles.disp <- get.dtt(list(anoles.dat), "disparity")

## Prepare data for plotting
anoles.cont <- anoles.cont[,which(colnames(anoles.cont) != "rep")]
anoles.all <- rbind(anoles.cont, anoles.disp)
cols.anoles <- c(cols.bl[1:ntraits], cols.gn[1:ntraits], cols.line)

## Supplementary Figure 5 -- anoles contrasts and dtt
fig.nh.dtt.emp(anoles.all, cols.anoles)


## ## Cyprinodon
## Read in data
#cyp <- readRDS("output/data/cyprinodon.rds")

## Compute principal components using the correlation approach as units of measurement differ between traits
#cyp.pca <- princomp(cyp$data, cor=TRUE)

## Compute phylogenetic principal components again, using the corrleation approach
#cyp.ppca <- phyl.pca(cyp$phy, cyp$data, mode="corr")

## Prepare dataset for model fitting
#cyp.dat <- list(tree=cyp$phy, raw=cyp$data, pc=cyp.pca$scores, ppc=cyp.ppca$S,
#                pcall=cyp.pca, ppcall=cyp.ppca)

## ### Fit models and compute AICw across all traits
#cyp.fit <- fitPCs(cyp.dat, seq_len(ncol(cyp$data)))

## Supplementary Figure 6 -- Model support across all PC and pPC axes
#cyp.df <- build.emp.data.step(cyp.fit)

#fig.aicw.empirical(cyp.df, cols.aic[c(1,3,2)])


## ### Node height test and Disparity through time
#ntraits <- ncol(cyp$data)
#cyp.cont <- get.contrasts(list(cyp.dat), "contrasts")
#cyp.disp <- get.dtt(list(cyp.dat), "disparity")

## Prepare data for plotting
#cyp.cont <- cyp.cont[,which(colnames(cyp.cont) != "rep")]
#cyp.all <- rbind(cyp.cont, cyp.disp)
#cols.cyp <- c(cols.bl[1:ntraits], cols.gn[1:ntraits], cols.line)

## Supplementary Figure 7 -- Cypridon contrasts and dtt
#fig.nh.dtt.emp(cyp.all, cols.cyp)

## Cyprinodon dataset is much less correlated: PC1 explains 39.8% of the variance with PCA
#round((cyp.pca$sdev^2)[1] / sum(cyp.pca$sdev^2), digits=3)
## and 32.0% of the variance with pPCA
#round(diag(cyp.ppca$Eval)[1] / sum(diag(cyp.ppca$Eval)), digits=3)


## ## Produce figure output

if (!interactive()){

    dev.off()
    ## create figure directory
    tmp <- dir("output")
    if (!"figs" %in% tmp)
        dir.create("output/figs")
    
    pdf("output/figs/mv-bm-aic.pdf", height=7, width=9)
    fig.aicw(bm.cor.df, cols.aic)
    dev.off()

    pdf("output/figs/uncor-bm-aic.pdf", height=7, width=9)
    fig.aicw(bm.uncor.df, cols.aic)
    dev.off()

    pdf("output/figs/uncor-ou-aic.pdf", height=7, width=9)
    fig.aicw(ou.ind.df, cols.aic)
    dev.off()
    
    pdf("output/figs/corr-ou-aic.pdf", height=7, width=9)
    fig.aicw(ou.cor.df, cols.aic)
    dev.off()

    pdf("output/figs/uncor-eb-aic.pdf", height=7, width=9)
    fig.aicw(eb.uncor.df, cols.aic)
    dev.off()

    pdf("output/figs/nh-2panel.pdf", height=7, width=9)
    fig.nh.2panel(cont, cols)
    dev.off()

    pdf("output/figs/dtt-2panel.pdf", height=7, width=9)
    fig.dtt.2panel(disp, cols)
    dev.off()

    pdf("output/figs/alpha-est-uncor.pdf", height=7, width=9)
    fig.alpha.est(ppca.ou.ind, col.pt=cols.gn[2], col.line=cols.line)
    dev.off()

    pdf("output/figs/alpha-est-corr.pdf", height=7, width=9)
    fig.alpha.est(ppca.ou.cor, col.pt=cols.gn[2], col.line=cols.line)
    dev.off()

    pdf("output/figs/onion.pdf", height=7, width=9)
    fig.rankslopes(rank.sim$rankslopes, rank.sim$exp.val, cols)
    dev.off()

    pdf("output/figs/felidae-aicw.pdf", height=7, width=9)
    fig.aicw.empirical(fel.df, cols.aic[c(1,3,2)])
    dev.off()

    pdf("output/figs/felidae_nh-dtt.pdf", height=7, width=9)
    fig.nh.dtt.emp(fel.all, cols.fel)
    dev.off()

 #   pdf("output/figs/cypri_aicw.pdf", height=7, width=9)
 #   fig.aicw.empirical(cyp.df, cols.aic[c(1,3,2)])
 #   dev.off()

 #   pdf("output/figs/cypri_nh-dtt.pdf", height=7, width=9)
 #   fig.nh.dtt.emp(cyp.all, cols.cyp)
 #   dev.off()
    
    pdf("output/figs/acdc_slopes.pdf", height=5, width=12)
    fig.acdc(acdcres, cols.aic)
    dev.off()
      
    pdf("output/figs/anoles_aicw.pdf", height=7, width=9)
    fig.aicw.empirical(anoles.df, cols.aic[c(1,3,2)])
    dev.off()
    
    pdf("output/figs/anoles_nh-dtt.pdf", height=7, width=9)
    fig.nh.dtt.emp(anoles.all, cols.anoles)
    dev.off()
}
