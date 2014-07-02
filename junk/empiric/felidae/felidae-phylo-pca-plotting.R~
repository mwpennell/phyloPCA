library(ggplot2)
library(reshape2)
library(gridExtra)
library(geiger)

## set colors
col <- c("#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD", "#117777", "#44AAAA", "#77CCCC", "#117744", "#44AA77", "#88CCAA", "#777711", "#AAAA44", "#DDDD77", "#774411", "#AA7744", "#DDAA77", "#771122", "#AA4455", "#DD7788")

## read in data
load("results_fitted.RData")

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

emp.rawcont <- get.contrasts(bio.mean, phy.car, "Original data")
emp.pccont <- get.contrasts(pc$scores, phy.car, "PCA")
emp.ppccont <- get.contrasts(ppc$S, phy.car, "Phylogenetic PCA")

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

emp.rawdtt <- get.dtt(bio.mean, phy.car)
emp.pcdtt <- get.dtt(pc$scores, phy.car)
emp.ppcdtt <- get.dtt(ppc$S, phy.car)

times <- unlist(lapply(1:ntraits, function(x) unlist(emp.rawdtt$dtt[[1]]$times)))

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

pdf("carnivora-nh-dtt.pdf", height = 10, width = 10)
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
names(emp.fit) <- c(1:19,"type","model")
emp.fit$type <- factor(emp.fit$type, levels = c("raw","pc","ppc"),
                       labels = c("Original data","PCA","Phylogenetic PCA"))
melt.fit <- melt(emp.fit, id=c('type','model'))
head(melt.fit)

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
p <- p + scale_x_discrete(breaks=c(seq(1, 19, by = 3)))
p <- p + xlab("Trait/PC axis")
p <- p + ylab("AICw")

pdf("carnivora-2.5-bio.clim.pdf", height = 10, width = 10)
p
dev.off()

## Heat map of the correlation among variables:
cr <- cor(bio.mean)
heat <- melt(cr)
head(heat)

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

pdf("carnivora-cor.pdf", width = 10, height = 10)
h
dev.off()

## Plot with the value of the alpha parameter:
alpha.raw <- vector()
alpha.pc <- vector()
alpha.ppc <- vector()

for(i in 1:19){
    alpha.raw[i] <- rawfits[[2]][[i]]$optpar
    alpha.pc[i] <- pcfits[[2]][[i]]$optpar
    alpha.ppc[i] <- ppcfits[[2]][[i]]$optpar
}

alpha <- data.frame(alpha.raw, alpha.pc, alpha.ppc, trait=1:19)
names(alpha) <- c("Original data", "PCA", "Phylogenetic PCA", 'trait')
alpha.m <- melt(alpha, id = 'trait')
alpha.m

p <- ggplot(alpha.m, aes(x = trait,y = value))
p <- p + geom_point()
p <- p + facet_grid(.~variable)
p <- p + scale_fill_manual(values="black")
p <- p + theme_bw()
p <- p + theme(strip.background=element_rect(fill="white"),
               plot.background=element_blank(),
               panel.grid.major=element_blank(),
               panel.grid.minor=element_blank(),
               legend.position="none")
p <- p + scale_x_discrete(breaks=c(seq(1, 19, by = 3)))
p <- p + xlab("Trait/PC axis")
p <- p + ylab("alpha paramter")

pdf("carnivora_alpha.pdf", width = 10, height = 5)
p
dev.off()
