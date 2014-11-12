## This is for plotting results and results analysis:

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
## Uncorrelated diag results.
ou.diag <- readRDS("output/sim-res/mv.ou-uncor-diag.rds")

## Figure 1 -- Model support across PC and pPC axes
ou.diag.df <- build.sim.data.step(ou.diag)
jpeg(file = "ou.diag.jpeg", quality = 90)
fig.aicw(ou.diag.df, cols.aic)
dev.off()

## Correlated mvBM (alpha matrix still diag). R matrix Postdef.
ou.cor.sig <- readRDS("output/sim-res/mv.ou-cor-sig.rds")
ou.cor.sig.df <- build.sim.data.step(ou.cor.sig)

## Figure 1 -- Model support across PC and pPC axes
jpeg(file = "ou.cor.sig.df.jpeg", quality = 90)
fig.aicw(ou.cor.sig.df, cols.aic)
dev.off()

## Correlated mvBM and alpha. A and R are Postdef.
ou.cor.sig.alp <- readRDS("output/sim-res/mv.ou-cor-sig-alpha.rds")
ou.cor.sig.alp.df <- build.sim.data.step(ou.cor.sig.alp)

## Figure 1 -- Model support across PC and pPC axes
jpeg(file = "ou.cor.sig.alp.df.jpeg", quality = 90)
fig.aicw(ou.cor.sig.alp.df, cols.aic)
dev.off()

## Now the figures for the DTT plots:
