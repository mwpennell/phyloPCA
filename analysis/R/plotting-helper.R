## Functions for plotting the data

## libraries for plotting
require(ggplot2)
require(reshape2)
require(RColorBrewer)
require(scales)
require(lattice)


## Function for building table for AICw figures
build.sim.data.step <- function(x){

    x$trait <- as.character(x$trait)
    raw <- data.frame()
    pca <- data.frame()
    ppc <- data.frame()
    for (i in seq_len(nrow(x))){
        raw.tmp <- cbind.data.frame(x[i,"trait"], (x[i,"raw_OUfixedRoot"] - x[i, "raw_EB"]))
        raw <- rbind(raw, raw.tmp)
        pca.tmp <- cbind.data.frame(x[i,"trait"], (x[i,"pc_OUfixedRoot"] - x[i, "pc_EB"]))
        pca <- rbind(pca, pca.tmp)
        ppc.tmp <- cbind.data.frame(x[i,"trait"], (x[i,"ppc_OUfixedRoot"] - x[i, "ppc_EB"]))
        ppc <- rbind(ppc, ppc.tmp)
    }

    colnames(raw) <- colnames(pca) <- colnames(ppc) <- c("trait", "dAIC")

    tr <- unique(raw$trait)
    m.raw <- sapply(tr, function(x) mean(raw$dAIC[raw$trait == x]))
    m.pca <- sapply(tr, function(x) mean(pca$dAIC[pca$trait == x]))
    m.ppc <- sapply(tr, function(x) mean(ppc$dAIC[ppc$trait == x]))

    q1.pca <- sapply(tr, function(x) {quantile(pca$dAIC[pca$trait == x], 0.75)})
    q2.pca <- sapply(tr, function(x) {quantile(pca$dAIC[pca$trait == x], 0.25)})
    q1.ppc <- sapply(tr, function(x) {quantile(ppc$dAIC[ppc$trait == x], 0.75)})
    q2.ppc <- sapply(tr, function(x) {quantile(ppc$dAIC[ppc$trait == x], 0.25)})

    list(m.raw=m.raw, m.pca=m.pca, m.ppc=m.ppc, q1.pca=q1.pca,
         q2.pca=q2.pca, q1.ppc=q1.ppc, q2.ppc=q2.ppc)

}





## Function for making plots of model suppport
fig.aicw <- function(x,cols){

    par(mar=c(5,4,2,2))
    plot(x$m.pca, type="s", lwd=1.5, las=1, xlab="PC axis", ylab="",
     col=cols[1], ylim=c(-1,1), axes=FALSE)
    ## add additional line edge to last point
    segments(20, x$m.pca[20], 21, x$m.pca[20], col=cols[1], lwd=1.5)
    lines(x$q1.pca, lwd=1, type="s", col=alpha(cols[1], 0.5))
    segments(20, x$q1.pca[20], 21, x$q1.pca[20], col=alpha(cols[1], 0.5), lwd=1)
    lines(x$q2.pca, lwd=1, type="s", col=alpha(cols[1], 0.5))
    segments(20, x$q2.pca[20], 21, x$q2.pca[20], col=alpha(cols[1], 0.5), lwd=1)
    tmp <- sapply(c(1:20), function(i) polygon.step(i, x$m.pca[i], x$q1.pca[i],
                                                col=cols[1]))
    tmp <- sapply(c(1:20), function(i) polygon.step(i, x$q2.pca[i], x$m.pca[i],
                                                col=cols[1]))


    lines(x$m.ppc, type="s", col=cols[2], lwd=1.5)
    ## add additional line edge to last point
    segments(20, x$m.ppc[20], 21, x$m.ppc[20], col=cols[2], lwd=1.5)
    lines(x$q1.ppc, lwd=1, type="s", col=alpha(cols[2], 0.5))
    segments(20, x$q1.ppc[20], 21, x$q1.ppc[20], col=alpha(cols[2], 0.5), lwd=1)
    lines(x$q2.ppc, lwd=1, type="s", col=alpha(cols[2], 0.5))
    segments(20, x$q2.ppc[20], 21, x$q2.ppc[20], col=alpha(cols[2], 0.5), lwd=1)
    tmp <- sapply(c(1:20), function(i) polygon.step(i, x$m.ppc[i], x$q1.ppc[i],
                                                col=cols[2]))
    tmp <- sapply(c(1:20), function(i) polygon.step(i, x$q2.ppc[i], x$m.ppc[i],
                                                col=cols[2]))

    abline(h=mean(x$m.raw), lwd=1.5, col=cols[3])
    axis(side=1, at=c(1:20))
    axis(side=2, at=c(-1, 1), labels=c("",""), las=1)
    legend("bottomright", c("PCA", "Phylogenetic PCA"),
       lwd=1.5, col=cols[1:2], bty="n", cex=1)
    parsave <- par(new=TRUE, mar=c(0,0,0,0))
    plot(NA, xlim=c(0,1), ylim=c(0,1), type="n", xlab="", ylab="", col="white")
    arrows(0.01, 0.65, 0.01, 0.9, length=0.1)
    arrows(0.01, 0.45, 0.01, 0.2, length=0.1)
    text(-0.01, 0.55, "BM", srt=90)
    text(-0.01, 0.775, "Support for OU", srt=90)
    text(-0.01, 0.325, "Support for EB", srt=90)
}


## Little function to help generate correct polygons
polygon.step <- function(i, y1, y2, col){
    px <- c(i, i+1, i+1, i)
    py <- c(y1, y1, y2, y2)
    polygon(px, py, col=alpha(col, 0.25), border=alpha(col, 0))
}



## Function for building table for AICw figures for empirical data
build.emp.data.step <- function(x){

    x$trait <- as.character(x$trait)
    raw <- data.frame()
    pca <- data.frame()
    ppc <- data.frame()
    for (i in seq_len(nrow(x))){
        raw.tmp <- cbind.data.frame(x[i,"trait"], (x[i,"raw_OUfixedRoot"] - x[i, "raw_EB"]))
        raw <- rbind(raw, raw.tmp)
        pca.tmp <- cbind.data.frame(x[i,"trait"], (x[i,"pc_OUfixedRoot"] - x[i, "pc_EB"]))
        pca <- rbind(pca, pca.tmp)
        ppc.tmp <- cbind.data.frame(x[i,"trait"], (x[i,"ppc_OUfixedRoot"] - x[i, "ppc_EB"]))
        ppc <- rbind(ppc, ppc.tmp)
    }

    colnames(raw) <- colnames(pca) <- colnames(ppc) <- c("trait", "dAIC")

    tr <- unique(raw$trait)
    m.raw <- sapply(tr, function(x) mean(raw$dAIC[raw$trait == x]))
    m.pca <- sapply(tr, function(x) mean(pca$dAIC[pca$trait == x]))
    m.ppc <- sapply(tr, function(x) mean(ppc$dAIC[ppc$trait == x]))

    list(m.raw=m.raw, m.pca=m.pca, m.ppc=m.ppc)

}

## Function for building table for AICw figures for empirical data
build.emp.data.stack <- function(x){
  x$trait <- as.character(x$trait)
  df <- melt(x)
  ids <- apply(do.call(rbind, strsplit(as.character(df[,2]), "_")), 2, factor)
  ids[ids[,2]=="OUfixedRoot",2] <- "OU"
  df <- data.frame(df[,1], ids, df[,3])
  colnames(df) <- c("trait", "variable", "type", "value")
  df[,1] <- factor(df[,1], levels = paste("trait", 1:length(unique(df[,1])), sep=""))
  return(df)
}

## Function for making plots of model suppport for the empirical datasets
fig.aicw.empirical <- function(df){
  .e <- environment()
  df$variable <- as.character(df$variable)
  df$variable[df$variable=="pc"] <- "PCA"; df$variable[df$variable=="ppc"] <- "Phylogenetic PCA"; df$variable[df$variable=="raw"] <- "Original Data"
  df$variable <- factor(df$variable, levels=c("Original Data", "PCA", "phylogenetic PCA"))
  p <- q <- ggplot(df, aes(trait, value, fill=type), environment = .e)
  p <- p +  geom_bar(data=df, stat="identity", position="stack")
  p <- p + scale_y_continuous(name="AICw")
  #p <- p + scale_fill_manual(values=col[c(6,21,12)], name="Model")
  p <- p + facet_grid(.~variable)
  p <- p + theme_bw()
  p <- p + theme(strip.background=element_rect(fill="white"),
                 plot.background=element_blank(),
                 panel.grid.major=element_blank(),
                 panel.grid.minor=element_blank(),
                 axis.text.x=element_blank(),
                 axis.ticks.x=element_blank())
  p <- p + xlab("Trait/PC axis")
  p  
}



## Function for plotting 2 panel NH test
fig.nh.2panel <- function(df, cols){
    .e <- environment()

    df1 <- subset(df, type != "raw")
    df2 <- subset(df, type == "raw")
    ## combine all the traits
    df2[,"variable"] <- factor("raw")
    dfpca <- dfppca <- df2
    dfpca[,"type"] <- factor("pc")
    dfppca[,"type"] <- factor("ppc")

    df.plot <- rbind(df1, dfpca, dfppca)

    df.plot$color.use <- apply(df.plot[,c("type", "variable")], 1, paste,
                               collapse="-")

    ## reorder
    df.plot <- df.plot[order(df.plot$type),]

    df.plot$color.use <- factor(df.plot$color.use,
                                levels=unique(df.plot$color.use))

    df.plot$type <- factor(df.plot$type, levels=c("pc", "ppc"),
                      labels=c("PCA", "Phylogenetic PCA"))

    df.plot$simmodel <- factor(df.plot$simmodel, levels=c("BM", "OU", "EB"))

    p <- ggplot(df.plot, aes(times, value, colour=color.use), environment=.e)
    p <- p + stat_smooth(method="lm",se=FALSE)
    p <- p + facet_grid(type~simmodel, scales="free_y")
    p <- p + scale_color_manual(values=c(cols[1:20], cols[41], cols[21:40], cols[41]))
    p <- p + theme_bw()
    p <- p + theme(strip.background=element_rect(fill="white"),
                   plot.background=element_blank(),
                   panel.grid.major=element_blank(),
                   panel.grid.minor=element_blank(),
                   axis.text.y=element_blank(),
                   axis.ticks.y=element_blank(),
                   axis.text.x=element_blank(),
                   axis.ticks.x=element_blank(),
                   legend.position="none")
    p <- p + xlab("Time")
    p <- p + ylab("Contrasts")   
    p
}


## Function for plotting 2-panel DTT plots
fig.dtt.2panel <- function(df, cols){
    .e <- environment()
    
    df1 <- subset(df, type != "raw")
    df2 <- subset(df, type == "raw")
    rm(df) ## save memory
    df2 <- subset(df2, variable == "trait.1")
    
    ## only use first trait
    df2[,"variable"] <- factor("raw")
    dfpca <- dfppca <- df2
    rm(df2) ## save memory
  
    
    dfpca[,"type"] <- factor("pc")
    dfppca[,"type"] <- factor("ppc")

    df.plot <- rbind(df1, dfpca, dfppca)

    df.plot$color.use <- apply(df.plot[,c("type", "variable")], 1, paste,
                               collapse="-")
    ## reorder
    df.plot <- df.plot[order(df.plot$type),]

    df.plot$color.use <- factor(df.plot$color.use,
                                levels=unique(df.plot$color.use))

    df.plot$type <- factor(df.plot$type, levels=c("pc", "ppc"),
                      labels=c("PCA", "Phylogenetic PCA"))

    df.plot$simmodel <- factor(df.plot$simmodel, levels=c("BM", "OU", "EB"))


                                   

    p <- ggplot(df.plot, aes(times, value, colour=color.use), environment=.e)
    p <- p + stat_smooth(se=FALSE)
    p <- p + facet_grid(type~simmodel, scales="free_y")
    p <- p + scale_color_manual(values=c(cols[1:20], cols[41], cols[21:40], cols[41]))
    p <- p + theme_bw()
    p <- p + theme(strip.background=element_rect(fill="white"),
                   plot.background=element_blank(),
                   panel.grid.major=element_blank(),
                   panel.grid.minor=element_blank(),
                   axis.text.y=element_blank(),
                   axis.ticks.y=element_blank(),
                   axis.text.x=element_blank(),
                   axis.ticks.x=element_blank(),
                   legend.position="none")
    p <- p + xlab("Time")
    p <- p + ylab("Contrasts")   
    p
}




## Figure for plotting the empirical NH and DTT together
fig.nh.dtt.emp <- function(df, cols){
    .e <- environment()
    
    df1 <- subset(df, type != "raw")
    df2 <- subset(df, type == "raw")
    rm(df) ## save memory
    df2 <- subset(df2, variable == "trait.1")
    
    ## only use first trait
    df2[,"variable"] <- factor("raw")
    dfpca <- dfppca <- df2
    
    dfpca[,"type"] <- factor("pc")
    dfppca[,"type"] <- factor("ppc")

    df.plot <- rbind(df1, dfpca, dfppca)

    df.plot$color.use <- apply(df.plot[,c("type", "variable")], 1, paste,
                               collapse="-")
    ## reorder
    df.plot <- df.plot[order(df.plot$type),]

    df.plot$color.use <- factor(df.plot$color.use,
                                levels=unique(df.plot$color.use))

    df.plot$type <- factor(df.plot$type, levels=c("pc", "ppc"),
                      labels=c("PCA", "Phylogenetic PCA"))

    df.plot$simmodel <- factor(df.plot$simmodel, levels=c("contrasts", "disparity"),
                               labels=c("Node height test", "Disparity through time"))


    ll <- length(unique(df.plot$variable)) - 1

    p <- ggplot(df.plot, aes(times, value, colour=color.use), environment=.e)
    p <- p + stat_smooth(data=subset(df.plot, simmodel=="Node height test"), method="lm", se=FALSE)
    p <- p + stat_smooth(data=subset(df.plot, simmodel=="Disparity through time"), se=FALSE)
    p <- p + facet_grid(type~simmodel, scales="free")
    p <- p + scale_color_manual(values=c(cols[1:ll], cols[(2*ll+1)], cols[(ll+1):(2*ll)], cols[(2*ll+1)]))
    p <- p + theme_bw()
    p <- p + theme(strip.background=element_rect(fill="white"),
                   plot.background=element_blank(),
                   panel.grid.major=element_blank(),
                   panel.grid.minor=element_blank(),
                   axis.text.y=element_blank(),
                   axis.ticks.y=element_blank(),
                   axis.text.x=element_blank(),
                   axis.ticks.x=element_blank(),
                   axis.title.y=element_blank(),
                   legend.position="none")
    p <- p + xlab("Time")  
    p
}




fig.rankslopes <- function(df, exp.val, cols){
    .e <- environment()
    
    df1 <- subset(df, type != "raw")
    df2 <- subset(df, type == "raw")
    rm(df) ## save memory
    
    ## only use first trait
    df2[,"variable"] <- factor("raw")
    dfpca <- dfppca <- df2
    rm(df2) ## save memory
  
    
    dfpca[,"type"] <- factor("pc")
    dfppca[,"type"] <- factor("ppc")

    df.plot <- rbind(df1, dfpca, dfppca)

    df.plot$color.use <- apply(df.plot[,c("type", "variable")], 1, paste,
                               collapse="-")
    ## reorder
    df.plot <- df.plot[order(df.plot$type),]

    df.plot$color.use <- factor(df.plot$color.use,
                                levels=unique(df.plot$color.use))

    df.plot$type <- factor(df.plot$type, levels=c("pc", "ppc"),
                      labels=c("PCA", "Phylogenetic PCA"))

   df.plot$simmodel <- as.numeric(as.character(factor(df.plot$simmodel, levels=1:length(exp.val), labels=exp.val)))
    p <- ggplot(df.plot, aes(simmodel, slope, colour=color.use), environment=.e)
    p <- p + stat_smooth(se=FALSE)
    p <- p + facet_grid(.~type, scales="free_y")
    p <- p + scale_color_manual(values=c(cols[1:20], cols[41], cols[21:40], cols[41]))
    p <- p + theme_bw()
    p <- p + theme(strip.background=element_rect(fill="white"),
                   plot.background=element_blank(),
                   panel.grid.major=element_blank(),
                   panel.grid.minor=element_blank(),
                   legend.position="none")
  p <- p + xlab("Proportion of variance explained by leading eigenvector")
  p <- p + ylab("Node-height test slope")   
  p
}



## Functions for plotting alpha parameter estimates
## Tufte style boxplots
panel.tuftebxp <- 
function (x, y, box.ratio = 1, box.width = box.ratio/(1 + box.ratio), horizontal=FALSE,
    pch = box.dot$pch, col = box.dot$col, 
    alpha = box.dot$alpha, cex = box.dot$cex, font = box.dot$font, 
    fontfamily = box.dot$fontfamily, fontface = box.dot$fontface, 
    fill = box.rectangle$fill, varwidth = FALSE, notch = FALSE, 
    notch.frac = 0.5, ..., levels.fos = if (horizontal) sort(unique(y)) else sort(unique(x)), 
    stats = boxplot.stats, coef = 1.5, do.out = TRUE, identifier = "bwplot") 
{
    if (all(is.na(x) | is.na(y))) 
        return()
    x <- as.numeric(x)
    y <- as.numeric(y)
    box.dot <- trellis.par.get("box.dot")
    box.rectangle <- trellis.par.get("box.rectangle")
    box.umbrella <- trellis.par.get("box.umbrella")
    plot.symbol <- trellis.par.get("plot.symbol")
    fontsize.points <- trellis.par.get("fontsize")$points
    cur.limits <- current.panel.limits()
    xscale <- cur.limits$xlim
    yscale <- cur.limits$ylim
    if (!notch) 
        notch.frac <- 0
    #removed horizontal code
     blist <- tapply(y, factor(x, levels = levels.fos), stats, 
            coef = coef, do.out = do.out)
        blist.stats <- t(sapply(blist, "[[", "stats"))
        blist.out <- lapply(blist, "[[", "out")
        blist.height <- box.width
        if (varwidth) {
            maxn <- max(table(x))
            blist.n <- sapply(blist, "[[", "n")
            blist.height <- sqrt(blist.n/maxn) * blist.height
        }
        blist.conf <- if (notch) 
            sapply(blist, "[[", "conf")
        else t(blist.stats[, c(2, 4), drop = FALSE])
        ybnd <- cbind(blist.stats[, 3], blist.conf[2, ], blist.stats[, 
            4], blist.stats[, 4], blist.conf[2, ], blist.stats[, 
            3], blist.conf[1, ], blist.stats[, 2], blist.stats[, 
            2], blist.conf[1, ], blist.stats[, 3])
        xleft <- levels.fos - blist.height/2
        xright <- levels.fos + blist.height/2
        xbnd <- cbind(xleft + notch.frac * blist.height/2, xleft, 
            xleft, xright, xright, xright - notch.frac * blist.height/2, 
            xright, xright, xleft, xleft, xleft + notch.frac * 
                blist.height/2)
        xs <- cbind(xbnd, NA_real_)
        ys <- cbind(ybnd, NA_real_)
        panel.segments(rep(levels.fos, 2), c(blist.stats[, 2], 
            blist.stats[, 4]), rep(levels.fos, 2), c(blist.stats[, 
            1], blist.stats[, 5]), col = box.umbrella$col, alpha = box.umbrella$alpha, 
            lwd = box.umbrella$lwd, lty = box.umbrella$lty, identifier = paste(identifier, 
                "whisker", sep = "."))

        if (all(pch == "|")) {
            mult <- if (notch) 
                1 - notch.frac
            else 1
            panel.segments(levels.fos - mult * blist.height/2, 
                blist.stats[, 3], levels.fos + mult * blist.height/2, 
                blist.stats[, 3], lwd = box.rectangle$lwd, lty = box.rectangle$lty, 
                col = box.rectangle$col, alpha = alpha, identifier = paste(identifier, 
                  "dot", sep = "."))
        }
        else {
            panel.points(x = levels.fos, y = blist.stats[, 3], 
                pch = pch, col = col, alpha = alpha, cex = cex, 
                 identifier = paste(identifier, 
                  "dot", sep = "."))
        }
        panel.points(x = rep(levels.fos, sapply(blist.out, length)), 
            y = unlist(blist.out), pch = plot.symbol$pch, col = plot.symbol$col, 
            alpha = plot.symbol$alpha, cex = plot.symbol$cex*0.5, 
            identifier = paste(identifier, "outlier", sep = "."))

}



fig.alpha.est <- function(df, col.pt, col.line){
    df$trait <- factor(df$trait, levels=as.character(unique(df$trait)),
                       labels=c(1:20))
    my.theme <- list(
        box.umbrella = list(col=col.pt),
        box.rectangle = list(fill=rep(c(col.pt,col.pt),2)),
        box.dot = list(col=col.pt, pch=19, cex=1),
        plot.symbol = list(cex=1, col=col.pt, pch=1))

    bwplot(alpha~trait, data=df, xlab="pPC axis", ylab=expression(alpha),
           par.settings=my.theme, panel=function(x,y,...){
               panel.abline(h=2, lwd=1.5, col=col.line)
               panel.tuftebxp(x=x,y=y,...)
           })
}
    








