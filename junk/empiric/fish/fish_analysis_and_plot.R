library(geiger)
library(phytools)
library(phylolm)

library(ggplot2)
library(reshape2)
library(gridExtra)

tr <- read.nexus("../../datasets/fish_finaltree_underscore.nex")
dt <- read.csv("../../datasets/fish_dataset.csv")

dt$species <- gsub(" ", "_", dt$species)
tr$tip.label <- gsub("'", "", tr$tip.label)

## Keep only the data for Cyprinodon species:
sp.split <- sapply(dt$species, FUN = function(x) strsplit(x, split = "_"))
cyp <- sapply(1:length(sp.split), FUN = function(x) sp.split[[x]][1])
id <- which(cyp == "Cyprinodon")
dt <- dt[id,]

## Fix some species names:
index <- which(tr$tip.label %in% dt$species)
which(!dt$species %in% tr$tip.label[index])
dt$species[9] <- "Cyprinodon_bozo"
dt$species[10] <- "Cyprinodon_bulldog"
dt$species[25] <- "Cyprinodon_normal_Crescent_Pond"

## Drop some tips from the tree:
tr <- drop.tip(phy = tr, tip = tr$tip.label[-index])

mm <- match(tr$tip.label, dt$species)

dt <- dt[,-c(17:25)]
rownames(dt) <- dt$species
dt <- dt[,-1]

## Plot correlation among traits (heatmap):
cr <- cor(dt)
heat <- melt(cr)

h <- ggplot(heat, aes(x = Var2, y=Var1))
h <- h + geom_tile(aes(fill = value))
h <- h + scale_fill_gradient2(low = "blue", mid = "white", high = "red")
h <- h + theme(strip.background=element_rect(fill="white"),
               plot.background=element_blank(),
               panel.grid.major=element_blank(),
               panel.grid.minor=element_blank(),
               axis.text.x=element_text(angle = -45, hjust = 0, vjust = 1))
h <- h + xlab("")
h <- h + ylab("")

pdf("fish-cor.pdf", width = 10, height = 10)
h
dev.off()

## Analysis:
pc <- princomp(dt)
ppc <- phyl.pca(tr, dt[mm,])

## Variance explained by PC1:
plot(pc)

## Fitting the models:
models <- c("BM", "OUfixedRoot", "EB")
rawfits <- list()
pcfits <- list()
ppcfits <- list()
for(i in 1:length(models)){
    rawfits[[i]] <- lapply(1:dim(dt)[2], function(x) phylolm(dt[mm,x]~1, phy=tr, model=models[i]))
    pcfits[[i]] <-  lapply(1:dim(dt)[2], function(x) phylolm(pc$scores[,x]~1, phy=tr, model=models[i]))
    ppcfits[[i]] <-  lapply(1:dim(dt)[2], function(x) phylolm(ppc$S[,x]~1, phy=tr, model=models[i]))
}

## Get the AICw:
bm.aic <- lapply(1:dim(dt)[2], function(x) (c(rawfits[[1]][[x]]$aic, pcfits[[1]][[x]]$aic, ppcfits[[1]][[x]]$aic)))
ou.aic <- lapply(1:dim(dt)[2], function(x) (c(rawfits[[2]][[x]]$aic, pcfits[[2]][[x]]$aic, ppcfits[[2]][[x]]$aic)))
eb.aic <- lapply(1:dim(dt)[2], function(x) (c(rawfits[[3]][[x]]$aic, pcfits[[3]][[x]]$aic, ppcfits[[3]][[x]]$aic)))
all.aicw <- lapply(1:dim(dt)[2], function(x) lapply(1:3, function(y) aicw(c(bm.aic[[x]][y], ou.aic[[x]][y], eb.aic[[x]][y]))$w))


## Make result tables:
bm.table <- sapply(all.aicw, function(x) sapply(x, function(y) y[1]))
ou.table <- sapply(all.aicw, function(x) sapply(x, function(y) y[2]))
eb.table <- sapply(all.aicw, function(x) sapply(x, function(y) y[3]))
rownames(bm.table) <- rownames(ou.table) <- rownames(eb.table) <- c("raw","pc","ppc")

## Now making all the plots:

## set colors
col <- c("#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD", "#117777", "#44AAAA", "#77CCCC", "#117744", "#44AA77", "#88CCAA", "#777711", "#AAAA44", "#DDDD77", "#774411", "#AA7744", "#DDAA77", "#771122", "#AA4455", "#DD7788")

## Figure 1 NH for 3 models (uncorrelated)
fig.nh.3panel <- function(df){
       .e <- environment()
       
    p <- ggplot(df, aes(times, value, colour=factor(variable)), environment=.e)
    p <- p + stat_smooth(method="lm",se=FALSE)
    p <- p + facet_grid(.~type, scales="free_y")
    p <- p + scale_color_manual(values=col, name="Trait/PC axis")
    p <- p + theme_bw()
    p <- p + theme(strip.background=element_rect(fill="white"),
                   plot.background=element_blank(),
                   panel.grid.major=element_blank(),
                   panel.grid.minor=element_blank(),
                   axis.text.y=element_blank(),
                   axis.ticks.y=element_blank(),
                   axis.text.x=element_blank(),
                   axis.ticks.x=element_blank())
    p <- p + xlab("Time")
    p <- p + ylab("Contrasts")   
    p
}

## Figure 2 DTT for 3 models (uncorrelated)
fig.dtt.3panel <- function(df){
       .e <- environment()

    p <- ggplot(df, aes(times, value, colour=factor(variable)), environment=.e)
    p <- p + stat_smooth(se=FALSE)
    p <- p + facet_grid(.~type, scales="free_y")
    p <- p + scale_color_manual(values=col, name="Trait/PC axis")
    p <- p + theme_bw()
    p <- p + theme(strip.background=element_rect(fill="white"),
                   plot.background=element_blank(),
                   panel.grid.major=element_blank(),
                   panel.grid.minor=element_blank(),
                   axis.text.y=element_blank(),
                   axis.ticks.y=element_blank(),
                   axis.text.x=element_blank(),
                   axis.ticks.x=element_blank())
    p <- p + xlab("Time")
    p <- p + ylab("Disparity")   
    p
  
}

## Get contrast data:

get.contrasts <- function(data, tree, type){
  btimes <- 1 - branching.times(tree)
  ntraits <- dim(data)[2]
  pic <- sapply(1:ntraits, function(y) pic(data[,y], tree))
  tr <- as.character(1:dim(data)[2])
  ll <- dim(pic)[1]
  variable <- unlist(data.frame(sapply(1:length(tr), function(x) rep(tr[x], times = ll))))
  times <- rep(btimes, times = length(tr))
  ty <- rep(type, times = length(times))
  value <- unlist(data.frame(pic))
  df <- data.frame(type=ty, times=times, variable=variable, value=value)
  return(df)
}

emp.rawcont <- get.contrasts(dt[mm,], tr, "Original data")
emp.pccont <- get.contrasts(pc$scores, tr, "PCA")
emp.ppccont <- get.contrasts(ppc$S, tr, "Phylogenetic PCA")

emp.df <- rbind(emp.rawcont, emp.pccont, emp.ppccont)

## Get disparity data:
get.dtt <- function(data, tree){
    ## Need the names of the species in row.names of the data.
    ntraits <- dim(data)[2]
    bytrait <- lapply(1:ntraits, function(x) data[,x])
    dtt <- lapply(1:ntraits, function(y) {tr <- bytrait[[y]]; names(tr) <- row.names(data); dtt(tree, tr, plot=FALSE)})
    ## dtt <- lapply(1:ntraits, function(y) dtt(tree, bytrait[[y]], plot=FALSE))
    disp <- lapply(1:ntraits, function(x) unlist(dtt[[x]]$dtt))
    return(list(dtt=dtt, disp=disp))
}

emp.rawdtt <- get.dtt(dt[mm,], tr)
emp.pcdtt <- get.dtt(pc$scores, tr)
emp.ppcdtt <- get.dtt(ppc$S, tr)

times <- unlist(lapply(1:dim(dt)[2], function(x) unlist(emp.rawdtt$dtt[[1]]$times)))

dttraw <- lapply(1:length(emp.rawdtt$dtt), function(x) emp.rawdtt$dtt[[x]]$dtt)
dttpc <- lapply(1:length(emp.pcdtt$dtt), function(x) emp.pcdtt$dtt[[x]]$dtt)
dttppc <- lapply(1:length(emp.ppcdtt$dtt), function(x) emp.ppcdtt$dtt[[x]]$dtt)

disp <- list(do.call(cbind, dttraw), do.call(cbind, dttpc), do.call(cbind, dttppc))
types <- c("Original data", "PCA", "Phylogenetic PCA")
dispdf <- lapply(1:3, function(x) data.frame(type=types[x], times=times, trait=disp[[x]]))
dispdf <- do.call(rbind, dispdf)
dispmelt <- melt(dispdf, id=c('type', 'times'))

## Making the plot with both data:

## Node heigth figure:
p <- fig.nh.3panel(emp.df)
## Disparity figure:
q <- fig.dtt.3panel(dispmelt)

g_legend <- function(a.gplot){
    tmp <- ggplot_gtable(ggplot_build(a.gplot))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend <- tmp$grobs[[leg]]
    return(legend)
}

legend <- g_legend(p)
lwidth <- sum(legend$width)

pdf("fish-nh-dtt.pdf", height = 10, width = 10)
grid.arrange(arrangeGrob(p + theme(legend.position="none"), q + theme(legend.position="none")),
                         legend, widths=unit.c(unit(1, "npc") - lwidth, lwidth), nrow=1)
dev.off()

## Make the box with the aicw:
## Prepare the dataset:

bm <- data.frame(bm.table)
bm$type <- row.names(bm)
bm$model <- "BM"

ou <- data.frame(ou.table)
ou$type <- row.names(ou)
ou$model <- "OU"

eb <- data.frame(eb.table)
eb$type <- row.names(eb)
eb$model <- "EB"

emp.fit <- rbind(bm, ou, eb)
names(emp.fit) <- c(1:dim(dt)[2],"type","model")
emp.fit$type <- factor(emp.fit$type, levels = c("raw","pc","ppc"),
                       labels = c("Original data","PCA","Phylogenetic PCA"))
melt.fit <- melt(emp.fit, id=c('type','model'))

p <- ggplot(melt.fit, aes(x = variable,y = value))
p <- p + geom_point()
p <- p + facet_grid(model~type)
p <- p + scale_fill_manual(values="black")
p <- p + theme_bw()
p <- p + theme(strip.background=element_rect(fill="white"),
               plot.background=element_blank(),
               panel.grid.major=element_blank(),
               panel.grid.minor=element_blank(),
               legend.position="none")
p <- p + scale_x_discrete(breaks=c(seq(1, dim(dt)[2], by = 1)))
p <- p + xlab("Trait/PC axis")
p <- p + ylab("AICw")

pdf("fish-models.pdf", height = 10, width = 10)
p
dev.off()

## The phenograms:
for(i in 1:5){
    pdf(paste("pheno.pc.",i,".pdf", sep = ""))
    phenogram(tr, setNames(pc$scores[,i], tr$tip.label), ftype="off")
    dev.off()
    
    pdf(paste("pheno.dt.",i,".pdf", sep = ""))       
    phenogram(tr, setNames(dt[,i], tr$tip.label), ftype="off")
    dev.off()
    
    pdf(paste("pheno.ppc.",i,".pdf", sep = ""))
    phenogram(tr, ppc$S[,i], ftype="off")
    dev.off()
}

## Save results:
save.image("fish_analysis_plot.RData")
