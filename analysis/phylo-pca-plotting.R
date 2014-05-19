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
build.sim.data.table <- function(x, trait.set){

    ## subset by desired trait set
    trait.names <- levels(x$trait)[trait.set]
    subx <- subset(x, trait %in% trait.names)
    subx$trait <- factor(subx$trait)

    ## creat data frame
    df <- melt(subx)

    ## add extra row to sort datasets
    type <- rep(NA, nrow(df))
    df <- cbind(df, type)

    ## break column names up
    df$variable <- as.character(df$variable)

    for (i in 1:nrow(df)){
        tmp <- strsplit(df[i,"variable"], split="_")
        df[i,"variable"] <- tmp[[1]][1]
        df[i,"type"] <- tmp[[1]][2]
    }

    ## recreate as factors
    df$trait <- factor(df$trait, levels=levels(df$trait), labels=trait.set)
    df$variable <- factor(df$variable, levels=c("raw", "pc", "ppc"),
                      labels=c("Original data", "PCA", "Phylogenetic PCA"))
    df$type <- factor(df$type, levels=c("BM", "OUfixedRoot", "EB"),
                      labels=c("BM", "OU", "EB"))
    df
}

    

fig.box.aicw <- function(df){

    p <- ggplot(df, aes(factor(trait), value))
    p <- p + geom_boxplot(aes(fill=factor(trait)))
    p <- p + facet_grid(type~variable)
    p <- p + scale_fill_manual(values=col[as.numeric(levels(df$trait))])
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



## set subset
trait.set <- c(1,2,5,10,15,20)


## Multivariate BM
load("simsBM_ind_50_20_.rds")
bm.ind <- build.sim.data.table(bmres, trait.set)
pdf("box-aicw-mvbm.pdf", height=8, width=9)
fig.box.aicw(bm.ind)
dev.off()


## Multivariate OU
load("simsOU_ind_50_20_.rds")
ou.ind <- build.sim.data.table(oures, trait.set)
pdf("box-aicw-mvou.pdf", height=8, width=9)
fig.box.aicw(ou.ind)
dev.off()

## Multivariate EB
load("simsEB_ind_50_20_.rds")
eb.ind <- build.sim.data.table(ebres, trait.set)
pdf("box-aicw-mveb.pdf", height=8, width=9)
fig.box.aicw(eb.ind)
dev.off()

## Multivariate and Correlated BM
## Note that different trait numbers were used so I am subsetting differently for now
load("sims_mv_50_20_.rds")
bm.cor <- build.sim.data.table(res, trait.set=seq_len(length(levels(res$trait))))

## traits 1,2,3,4,5,10,15,20
cor.tr <- c(1,2,5,6,7,8)
bm.cor <- subset(bm.cor, trait %in% cor.tr)
bm.cor$trait <- factor(bm.cor$trait)
bm.cor$trait <- factor(bm.cor$trait, levels=levels(bm.cor$trait),
                       labels=c(1,2,5,10,15,20))
pdf("box-aicw-corbm.pdf", height=8, width=9)
fig.box.aicw(bm.cor)
dev.off()






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
                   axis.text.y=element_text(color="white", size=3.75),
                   axis.title.x=element_blank())
    legend <- g.legend(p)
    lwidth <- sum(legend$width)
    


    q <- q + geom_hline(aes(yintercept=2), color=col[5], size=2)
    q <- q + geom_boxplot(data=subset(df, est=="Par" & simmodel=="alpha"),fill=col[11], color="black", outlier.size = 1)
    q <- q + facet_grid(.~type)
    q <- q + theme_bw()
    q <- q + theme(strip.background=element_rect(fill="white"),
                   plot.background=element_blank(),
                   panel.grid.major=element_blank(),
                   panel.grid.minor=element_blank(),
                   legend.position="none")
    q <- q + ylim(c(0,20))
    q <- q + ylab(expression(alpha))
    q <- q + xlab("Trait/PC axis")
    q <- strip.remover(q, "x")

    grid.arrange(arrangeGrob(p + theme(legend.position="none"),legend,q,
                             widths=unit.c(unit(1, "npc") - lwidth,lwidth), nrow=2))
    
}


load("OUmeanAICw.rds")
load("OUsimParameterEstimates.rds")
## do a bit of processing
parsdf$trait <- as.character(parsdf$trait)
parsdf$trait <- sapply(parsdf$trait, function(x) sub(".", "", x, fixed=TRUE))
parsdf$trait <- factor(parsdf$trait)
pars.tmp <- melt(parsdf)
colnames(pars.tmp) <- colnames(OUmeanAICw)

## combine the dataframes
est <- rep("AICw", nrow(OUmeanAICw))
OUmeanAICw <- cbind(OUmeanAICw, est)
est <- rep("Par", nrow(pars.tmp))
pars.tmp <- cbind(pars.tmp, est)
oudf <- rbind(OUmeanAICw, pars.tmp)

## prune down dataset to raw and phylo pca only
oudf <- subset(oudf, type %in% c("raw", "ppc"))
oudf$type <- factor(oudf$type)
oudf$type <- factor(oudf$type, levels=c("raw", "ppc"), labels=c("Original data", "Phylogenetic PCA"))

## prune down dataset to exclude only BM, OU, EB, alpha
oudf <- subset(oudf, simmodel %in% c("BM", "EB", "OUfixedRoot", "alpha"))
oudf$simmodel <- factor(oudf$simmodel)
oudf$simmodel <- factor(oudf$simmodel, levels=c("BM", "OUfixedRoot", "EB", "alpha"),
                        labels=c("BM", "OU", "EB", "alpha"))

oudf$trait <- factor(oudf$trait, levels=unique(oudf$trait), labels=c(1:20))


pdf("model-support-alpha.pdf", height=7, width=10, onefile = FALSE)
fig.model.support.alpha(oudf)
dev.off()











