## I am going to try this very nice package.
## Check this website: http://ropensci.org/blog/2014/04/22/rwbclimate-sp/
## library(devtools)
## devtools::install_github("spocc", "ropensci")

library(geiger)
library(spocc)
library(maptools)
library(raster)
library(plyr)
library(phytools)
library(phylolm)

load("climate.data.RData")

## phy <- read.nexus("1741-7007-10-12-s5.nex")[[1]]
## sp.total <- phy$tip.label
## ## Take out the outgroup species.
## sp <- sp.total[which(!sp.total %in% c("Equus_caballus","Bos_taurus","Artibeus_jamaicensis"
##                                      ,"Mystacina_tuberculata","Tadarida_brasiliensis","Homo_sapiens"
##                                      ,"Rattus_norvegicus","Mus_musculus"))]

## Get occurrence data for the species:
## dt <- occ(query = sp, from = c("gbif"), limit = 100)
## dtfix <- fixnames(dt, how = "query")
## dt.df <- occ2df(dtfix)
## head(dt.df) ## Nice!
## dim(dt.df)

## Plot to see where the species are:
data(wrld_simpl)
## plot(wrld_simpl, col = "light yellow", axes = T)
## points(dt.df$longitude, dt.df$latitude, col = "red", cex = 0.25)

## Get climate data:

## Data downloaded from: http://biogeo.ucdavis.edu/
## Using the 2-5m (accuracy) bio_clim data.

## BIO1 = Annual Mean Temperature (amT)
## BIO2 = Mean Diurnal Range (Mean of monthly (max temp - min temp)) (mdrT)
## BIO3 = Isothermality (BIO2/BIO7) (* 100) (isoT)
## BIO4 = Temperature Seasonality (standard deviation *100) (Tseas)
## BIO5 = Max Temperature of Warmest Month (TmaxWM)
## BIO6 = Min Temperature of Coldest Month (TminCM)
## BIO7 = Temperature Annual Range (BIO5-BIO6) (Tar)
## BIO8 = Mean Temperature of Wettest Quarter (TmeanWeQ)
## BIO9 = Mean Temperature of Driest Quarter (TmeanDQ)
## BIO10 = Mean Temperature of Warmest Quarter (TmeanWaQ)
## BIO11 = Mean Temperature of Coldest Quarter (TmeanCQ)
## BIO12 = Annual Precipitation (aP)
## BIO13 = Precipitation of Wettest Month (PWeM)
## BIO14 = Precipitation of Driest Month (PDM)
## BIO15 = Precipitation Seasonality (Coefficient of Variation) (Pseas)
## BIO16 = Precipitation of Wettest Quarter (PWeQ)
## BIO17 = Precipitation of Driest Quarter (PDQ)
## BIO18 = Precipitation of Warmest Quarter (PWaQ)
## BIO19 = Precipitation of Coldest Quarter (PCQ)

## nms <- c("amT","mdrT","isoT","Tseas","TmaxWM","TminCM","Tar","TmeanWeQ","TmeanDQ","TmeanWaQ","TmeanCQ","aP"
##          ,"PWeM","PDM","Pseas","PWeQ","PDQ","PWaQ","PCQ")
## index <- 1:19

## climate <- stack()
## for(i in index){
##     r <- raster(paste("~/Documents/Academicos/Harmon_Lab/Projects/PhyloPCA_paper/bio_2-5m_bil/bio",index[i]
##                       ,".bil",sep = ""))
##     climate <- stack(climate, r)
## }
## nlayers(climate)

## ## plot(climate@layers[[12]], main = nms[12]) ## Annual Precipitation plot.
## ## points(dt.df$longitude, dt.df$latitude, col = "red", cex = 0.25)

## ## Extract the bioclimate data for each occurrence:
## bioclim <- extract(climate, dt.df[,2:3], method='simple', df=TRUE)
## bioclim <- bioclim[,-1]
## names(bioclim) <- c(nms)

## biobuff <- extract(climate, dt.df[,2:3], method='simple', buffer=500, fun=mean, df=TRUE)
## This will get the climate data for a 5km radius, then give the mean for each point.

## Make data.frame with everything:
bio.dt <- data.frame(species=dt.df$name, bioclim)
names(bio.dt)

## Get mean bioclimate data by species:
mean.c <- colwise(mean, na.rm = TRUE)
bio.mean <- ddply(bio.dt, "species", function(x) mean.c(x) )

## Check the data (it has NA's):
summary(bio.mean) ## Yep.

## Take out Na:
ss <- apply(bio.mean[,-1], 2, function(x) is.na(x) )
ss[,1]
bio.mean <- bio.mean[!ss[,1],]
summary(bio.mean)

## Drop tips that are not in the dataset:
bio.mean$species <- gsub(" ", "_", bio.mean$species)
phy.car <- drop.tip(phy, tip = phy$tip.label[which(!phy$tip.label %in% bio.mean[,1])])
rownames(bio.mean) <- bio.mean$species ## phyl.pca requirement.
bio.mean <- bio.mean[,-1]

## Find PCs and Phylogenetic PCs.
phy.car <- multi2di(phy.car)

pc <- princomp(bio.mean)
ppc <- phyl.pca(phy.car, bio.mean)

#################################################################################

## Fitting the models with log-transformation:
models <- c("BM", "OUrandomRoot", "EB")
rawfits <- list()
pcfits <- list()
ppcfits <- list()

index <- 1:dim(bio.mean)[2]

for(i in 1:length(models)){
    rawfits[[i]] <- lapply(index, function(x) phylolm(bio.mean[,x]~1, phy=phy.car, model=models[i]))
    pcfits[[i]] <-  lapply(index, function(x) phylolm(pc$scores[,x]~1, phy=phy.car, model=models[i]))
    ppcfits[[i]] <-  lapply(index, function(x) phylolm(ppc$S[,x]~1, phy=phy.car, model=models[i]))
}

## Get the AICw:
bm.aic <- lapply(index, function(x) (c(rawfits[[1]][[x]]$aic, pcfits[[1]][[x]]$aic, ppcfits[[1]][[x]]$aic)))
ou.aic <- lapply(index, function(x) (c(rawfits[[2]][[x]]$aic, pcfits[[2]][[x]]$aic, ppcfits[[2]][[x]]$aic)))
eb.aic <- lapply(index, function(x) (c(rawfits[[3]][[x]]$aic, pcfits[[3]][[x]]$aic, ppcfits[[3]][[x]]$aic)))
all.aicw <- lapply(index, function(x) lapply(1:3, function(y) aicw(c(bm.aic[[x]][y]
                                                                     , ou.aic[[x]][y], eb.aic[[x]][y]))$w))

## Make result tables:
bm.table <- sapply(all.aicw, function(x) sapply(x, function(y) y[1]))
ou.table <- sapply(all.aicw, function(x) sapply(x, function(y) y[2]))
eb.table <- sapply(all.aicw, function(x) sapply(x, function(y) y[3]))
rownames(bm.table) <- rownames(ou.table) <- rownames(eb.table) <- c("raw","pc","ppc")

## Plots:
pdf("Carnivora_2-5_bioclim.pdf")
par(mfrow = c(3,3))
for(i in 1:3){
    plot(index, bm.table[i,], main = paste(models[1],rownames(bm.table)[i],sep="_")
         , ylim = c(0.0,1.0), ylab = "AICw", xlab = "PCs")
    plot(index, ou.table[i,], main = paste(models[2],rownames(ou.table)[i],sep="_")
         , ylim = c(0.0,1.0), ylab = "AICw", xlab = "PCs")
    plot(index, eb.table[i,], main = paste(models[3],rownames(eb.table)[i],sep="_")
         , ylim = c(0.0,1.0), ylab = "AICw", xlab = "PCs")
}
dev.off()

### --- ### --- ### --- ### --- ### --- ### --- ### --- ### --- ###
## Some graphs to visualize the data:

par(mfrow = c(3,4))
for(i in 1:12){
    r <- bio.mean[,i]
    names(r) <- rownames(bio.mean)
    phenogram(phy.car, r)
}
dev.copy2pdf()

par(mfrow = c(3,4))
for(i in 1:12){
    phenogram(phy.car, pc$scores[,i])
}
dev.copy2pdf()

par(mfrow = c(3,4))
for(i in 1:12){
    phenogram(phy.car, ppc$S[,i])
}
dev.copy2pdf()
