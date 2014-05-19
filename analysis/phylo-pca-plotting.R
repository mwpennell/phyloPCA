library(ggplot2)
library(reshape2)
library(gridExtra)

## set colors
col <- c("#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD", "#117777", "#44AAAA", "#77CCCC", "#117744", "#44AA77", "#88CCAA", "#777711", "#AAAA44", "#DDDD77", "#774411", "#AA7744", "#DDAA77", "#771122", "#AA4455", "#DD7788")

## read in data
load("contrastsdata.rds")
load("disparitydata.rds")


## little function for changing trait names
tr.nm <- function(x){
    sapply(x, function(y) {paste("trait", y, sep=".")})
}


## Figure 1 NH for 3 models (uncorrelated)
fig.nh.3panel <- function(df){
       .e <- environment()

    df$type <- factor(df$type, levels=c("raw", "pc", "ppc"),
                      labels=c("Original data", "PCA", "Phylogenetic PCA"))

    df$simmodel <- factor(df$simmodel, levels=c("BM", "OU", "EB"))

    ll <- length(levels(df$variable))
    df$variable <- factor(df$variable, levels=tr.nm(seq_len(ll)), labels=seq_len(ll))

    p <- ggplot(df, aes(times, value, colour=factor(variable)), environment=.e)
    p <- p + stat_smooth(method="lm",se=FALSE)
    p <- p + facet_grid(simmodel~type, scales="free_y")
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

    df$type <- factor(df$type, levels=c("raw", "pc", "ppc"),
                      labels=c("Original data", "PCA", "Phylogenetic PCA"))

    df$simmodel <- factor(df$simmodel, levels=c("BM", "OU", "EB"))

    ll <- length(levels(df$variable))
    df$variable <- factor(df$variable, levels=tr.nm(seq_len(ll)), labels=seq_len(ll))

    p <- ggplot(df, aes(times, value, colour=factor(variable)), environment=.e)
    p <- p + stat_smooth(se=FALSE)
    p <- p + facet_grid(simmodel~type, scales="free_y")
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

pdf("nh-3models.pdf", height=9, width=10)
fig.nh.3panel(contmelt)
dev.off()

pdf("dtt-3models.pdf", height=9, width=10)
fig.dtt.3panel(dispdat)
dev.off()


## Figure 3

## function to restructure sim data
build.sim.data.table <- function(x){
    df <- melt(x)
    type <- rep(NA, nrow(df))
    df <- cbind(df, type)
    df$variable <- as.character(df$variable)

    for (i in 1:nrow(df)){
        tmp <- strsplit(df[i,"variable"], split="_")
        df[i,"variable"] <- tmp[[1]][1]
        df[i,"type"] <- tmp[[1]][2]
    }
    df
}

    

fig.box.aicw <- function(df){

    df$type <- factor(df$type, levels=c("BM", "OUfixedRoot", "EB"),
                      labels=c("BM", "OU", "EB"))
    df$variable <- factor(df$variable, levels=c("raw", "pc", "ppc"),
                      labels=c("Original data", "PCA", "Phylogenetic PCA"))

    df$trait <- factor(df$trait, levels=levels(df$trait), labels=c(1,2,3,5,7,10,13,16, 20))
    
    p <- ggplot(df, aes(factor(trait), value))
    p <- p + geom_boxplot(aes(fill=factor(trait)))
    p <- p + facet_grid(type~variable)
    p <- p + scale_fill_manual(values=col[c(1,2,3,5,7,10,13,16, 20)])
    p <- p + theme_bw()
    p <- p + theme(strip.background=element_rect(fill="white"),
                   plot.background=element_blank(),
                   panel.grid.major=element_blank(),
                   panel.grid.minor=element_blank(),
                   legend.position="none")
    p <- p + xlab("Trait/PC axis")
    p <- p + ylab("AICw")
    p <- p 
    p
}

## Multivariate BM
load("simsBM_ind_50_20_.rds")
bm.ind <- build.sim.data.table(bmres)
pdf("box-aicw-mvbm.pdf", height=8, width=9)
fig.box.aicw(bm.ind)
dev.off()


## Multivariate OU
load("simsOU_ind_50_20_.rds")
ou.ind <- build.sim.data.table(oures)
pdf("box-aicw-mvou.pdf", height=8, width=9)
fig.box.aicw(ou.ind)
dev.off()

## Multivariate EB
load("simsEB_ind_50_20_.rds")
eb.ind <- build.sim.data.table(ebres)
pdf("box-aicw-mveb.pdf", height=8, width=9)
fig.box.aicw(eb.ind)
dev.off()

## Note trait numbers differ so I have used a different form of the figure fxn for now

## Multivariate and Correlated BM
load("sims_mv_50_20_.rds")
bm.cor <- build.sim.data.table(res)






## Figure 4 -- Model support and alpha estimation
## Will need to modified once we figure out structure of data

## Supporting function to remove strips
strip.remover <- function(ggp, what="x") {
  require(gridExtra)

  zeroGrob <- function() {
    g0 <- grob(name="NULL")
    class(g0) <- c("zeroGrob",class(g0))
    g0
  }

  g <- ggplotGrob(ggp)

  g$grobs <- lapply(g$grob, function(gr) {
    if (any(grepl(paste0("strip.text.", what),names(gr$children)))) {
      gr$children[[grep("strip.background",names(gr$children))]] <- zeroGrob()
      gr$children[[grep("strip.text",names(gr$children))]] <- zeroGrob()
    }
    return(gr)
  }
  )

  class(g) = c("arrange", "ggplot",class(g)) 
  g
}

## Supporting function to place legend properly
g.legend<-function(a.gplot){
tmp <- ggplot_gtable(ggplot_build(a.gplot))
leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
legend <- tmp$grobs[[leg]]
return(legend)}




fig.model.support.alpha <- function(df){
    .e <- environment()

    df$simmodel <- factor(df$simmodel,
                          levels=c("BM", "OUfixedRoot", "EB", "alpha", "sig2", "halflife"),
                          labels=c("BM", "OU", "EB", "alpha", "sig2", "halflife"))
    df$trait <- factor(df$trait, levels=c(paste("trait", 1:20, sep="")), labels=1:20)
    df$type <- factor(df$type, levels=c("raw", "ppc"),
                      labels=c("Original data", "Phylogenetic PCA"))

    p <- q <- ggplot(df, aes(trait, value, fill=simmodel), environment = .e)
    p <- p +  geom_bar(data=subset(df, est == "AICw"), stat="identity", position="stack")
    p <- p + scale_y_continuous(name="AICw")
    p <- p + scale_fill_manual(values=col[c(2,11,14)], name="Model")
    p <- p + facet_grid(.~type)
    p <- p + theme_bw()
    p <- p + theme(strip.background=element_rect(fill="white"),
                   plot.background=element_blank(),
                   panel.grid.major=element_blank(),
                   panel.grid.minor=element_blank(),
                   axis.ticks.y=element_blank(),
                   axis.text.y=element_blank(),
                   axis.title.x=element_blank())
    legend <- g.legend(p)
    lwidth <- sum(legend$width)
    

    q <- q + geom_hline(aes(yintercept=2), color=col[2])
    q <- q + geom_point(data=subset(df, est=="Alpha" & simmodel=="alpha"),
                        alpha=0.75, color=col[11])
    q <- q + ylim(c(0, 20))
    q <- q + facet_grid(.~type)
    q <- q + theme_bw()
    q <- q + theme(strip.background=element_rect(fill="white"),
                   plot.background=element_blank(),
                   panel.grid.major=element_blank(),
                   panel.grid.minor=element_blank(),
                   legend.position="none")
    q <- q + ylab(expression(alpha))
    q <- q + xlab("Trait/PC axis")
    q <- strip.remover(q, "x")

    grid.arrange(arrangeGrob(p + theme(legend.position="none"),legend,q,
                             widths=unit.c(unit(1, "npc") - lwidth,lwidth), nrow=2))
    
}


load("OUmeanAICw.rds")
load("OUsimParameterEstimates.rds")
## do a bit of processing
pars.tmp <- parsdf[,-ncol(parsdf)]
colnames(pars.tmp) <- colnames(OUmeanAICw)
pars.tmp$value <- pars.tmp$simmodel
pars.tmp$simmodel <- "alpha"
pars.tmp$trait <- as.character(pars.tmp$trait)
pars.tmp$trait <- sapply(pars.tmp$trait, function(x) sub(".", "", x, fixed=TRUE))

## combine the dataframes
est <- rep("AICw", nrow(OUmeanAICw))
OUmeanAICw <- cbind(OUmeanAICw, est)
est <- rep("Alpha", nrow(pars.tmp))
pars.tmp <- cbind(pars.tmp, est)
oudf <- rbind(OUmeanAICw, pars.tmp)

## remove standard PC
oudf$type <- as.character(oudf$type)
oudf <- oudf[-which(oudf$type == "pc"),]
oudf$type <- as.factor(oudf$type)
oudf$trait <- as.factor(oudf$trait)

pdf("model-support-alpha.pdf", height=7, width=9, onefile = FALSE)
fig.model.support.alpha(oudf)
dev.off()











