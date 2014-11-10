## Functions to help with simulating and analyzing data

## libraries
require(phytools)
require(geiger)
require(phylolm)
require(foreach)
require(doMC)
require(MASS)
require(nlme)
require(reshape2)
require(plyr)
require(mvSLOUCH)

## Function for simulating positive definite covariance matrices
Posdef <- function (n, ev = rexp(n, 1/100)) {
  Z <- matrix(ncol=n, rnorm(n^2))
  decomp <- qr(Z)
  Q <- qr.Q(decomp) 
  R <- qr.R(decomp)
  d <- diag(R)
  ph <- d / abs(d)
  O <- Q %*% diag(ph)
  Z <- t(O) %*% diag(ev) %*% O
  return(Z)
}

## Function for simulating correlated multivariate Brownian motion
sim.tree.pcs.mv <- function(ntips, traits, sig2dist=rexp, lambda=0.1, mu=0, ...){
  tree <- pbtree(b=lambda, d=mu, n=ntips)
  tree$edge.length <- tree$edge.length/max(branching.times(tree))
  Sigma <- vcv.phylo(tree)
  R <- Posdef(traits, ev=sig2dist(traits, ...))
  simdat <- t(matrix(mvrnorm(1, rep(0, traits*ntips), R %x% Sigma), nrow=traits, ncol=ntips, byrow=TRUE))
  pc <- princomp(simdat)
  ppc <- phyl.pca(tree,simdat)
  row.names(simdat) <- tree$tip.label
  rownames(pc$scores) <- tree$tip.label
  rownames(ppc$S) <- tree$tip.label
  return(list(tree=tree, raw=simdat, pc=pc$scores, ppc=ppc$S, pcall=pc, ppcall=ppc, R=R))
}


## Function for simulating uncorrelated Brownian motion data
sim.tree.pcs.ind <- function(ntips, traits, sig2dist=rexp, lambda=0.1, mu=0, ...){
  tree <- pbtree(b=lambda, d=mu, n=ntips)
  tree$edge.length <- tree$edge.length/max(branching.times(tree))
  #outree <- rescale(tree,alpha=log(2),model="OU")
  X <- cbind(sapply(1:traits,function(x) fastBM(tree, sig2=sig2dist(1, ...), nsim=1)))
  pc <- princomp(X)
  ppc <- phyl.pca(tree,X)
  return(list(tree=tree, raw=X, pc=pc$scores, ppc=ppc$S, pcall=pc, ppcall=ppc))
}

## Function for simulating uncorrelated (independent) Ornstein-Uhlenbeck data
## This give the same results from 'sim.tree.pcs.mv.ou.uncor' using same values for
##   alpha, sig2 and root = 0.
sim.tree.pcs.ind.ou <- function(ntips, traits, alpha, sig2, lambda=0.1, mu=0, ...){
  tree <- pbtree(b=lambda, d=mu, n=ntips)
  tree$edge.length <- tree$edge.length/max(branching.times(tree))
  tr <- rescale(tree, model="OU", alpha)
  #outree <- rescale(tree,alpha=log(2),model="OU")
  ## Default for fastBM root state is 0.
  X <- sapply(1:traits, function(x) fastBM(tr, sig2=sig2, nsim=1))
  pc <- princomp(X)
  ppc <- phyl.pca(tree,X)
  return(list(tree=tree, raw=X, pc=pc$scores, ppc=ppc$S, pcall=pc, ppcall=ppc))
}

## Function for simulating multivariate uncorrelated Ornstein-Uhlenbeck data under mvSLOUCH.
## Results from this function should be the same as using 'sim.tree.pcs.ind.ou'.
sim.tree.pcs.mv.ou.uncor <- function(ntips, traits, sig2, lambda=0.1, mu=0, root=0, alpha, ...){
  ## A single mv optima set to be equal to root (same optima value for every trait).
  ## Alpha ans sig2 matrices are diag() of alpha. Off-diag equal to 0.
  tree <- pbtree(b=lambda, d=mu, n=ntips)
  tree$edge.length <- tree$edge.length/max(branching.times(tree))
  ouchtree <- ape2ouch(tree, scale=1)
  ouchtree@nodelabels[1:(ouchtree@nnodes-ouchtree@nterm)] <- as.character(1:(ouchtree@nnodes-ouchtree@nterm))
  ## Note here both R and A matrices are diag.
  R <- diag(sig2, nrow=traits)
  A <- diag(alpha, nrow=traits)
  Y0 <- rep(root, traits)
  mPsi <- matrix(rep(root, traits), ncol=1)
  simdat <- simulOUCHProcPhylTree(ouchtree, list(vY0=Y0, Syy=R, A=A, mPsi=mPsi))
  simdat <- simdat[-c(1:(ouchtree@nnodes-ouchtree@nterm)),]
  pc <- princomp(simdat)
  ppc <- phyl.pca(tree, simdat)
  row.names(simdat) <- tree$tip.label
  rownames(pc$scores) <- tree$tip.label
  rownames(ppc$S) <- tree$tip.label
  return(list(tree=tree, raw=simdat, pc=pc$scores, ppc=ppc$S, pcall=pc, ppcall=ppc, R=R))
}

## Function for simulating multivariate correlated Ornstein-Uhlenbeck data.
sim.tree.pcs.mv.ou <- function(ntips, traits, sig2dist=rexp, lambda=0.1, mu=0, root=0, alpha, diag=FALSE, ...){
  ## A single mv optima set to be equal to root (same optima value for every trait).
  ## Alpha matrix is diag() of alpha. Off-diag equal to 0.
  ## diag = FALSE to use a Posdef() R matrix; diag = TRUE to use a diag() R matrix.
  tree <- pbtree(b=lambda, d=mu, n=ntips)
  tree$edge.length <- tree$edge.length/max(branching.times(tree))
  ouchtree <- ape2ouch(tree, scale=1)
  ouchtree@nodelabels[1:(ouchtree@nnodes-ouchtree@nterm)] <- as.character(1:(ouchtree@nnodes-ouchtree@nterm))
  R <- Posdef(traits, ev=sig2dist(traits, ...))
  R[lower.tri(R)] <- 0
  if(diag == TRUE){
      R[upper.tri(R)] <- 0
  }
  Y0 <- rep(root, traits)
  mPsi <- matrix(rep(root, traits), ncol=1)
  A <- diag(alpha, nrow=traits)
  simdat <- simulOUCHProcPhylTree(ouchtree, list(vY0=Y0, Syy=R, A=A, mPsi=mPsi))
  simdat <- simdat[-c(1:(ouchtree@nnodes-ouchtree@nterm)),]
  pc <- princomp(simdat)
  ppc <- phyl.pca(tree, simdat)
  row.names(simdat) <- tree$tip.label
  rownames(pc$scores) <- tree$tip.label
  rownames(ppc$S) <- tree$tip.label
  return(list(tree=tree, raw=simdat, pc=pc$scores, ppc=ppc$S, pcall=pc, ppcall=ppc, R=R))
}

## Function for simulating uncorrelated Early-Burst data
sim.tree.pcs.ind.eb <- function(ntips, traits, a, sig2, lambda=0.1, mu=0, ...){
  tree <- pbtree(b=lambda, d=mu, n=ntips)
  tree$edge.length <- tree$edge.length/max(branching.times(tree))
  tr <- rescale(tree, model="EB", a)
  #outree <- rescale(tree,alpha=log(2),model="OU")
  X <- sapply(1:traits, function(x) fastBM(tr,sig2=sig2, nsim=1))
  pc <- princomp(X)
  ppc <- phyl.pca(tree,X)
  return(list(tree=tree, raw=X, pc=pc$scores, ppc=ppc$S, pcall=pc, ppcall=ppc))
}

## # Analysis functions
## Function for fitting models to raw, pc and ppca datasets
fitPCs <- function(simdat, trait.seq=1, models = c("BM", "OUfixedRoot", "EB")){
  transforms <- c("raw", "pc", "ppc")
  lower.bounds <- c(NULL, 10^-12, NULL)
  rawfits <- foreach(i=trait.seq) %dopar% suppressWarnings(lapply(1:length(models), function(x) phylolm(simdat$raw[,i]~1, phy=simdat$tree, model=models[x])))
  pcfits <-  foreach(i=trait.seq) %dopar% suppressWarnings(lapply(1:length(models), function(x) phylolm(simdat$pc[,i]~1, phy=simdat$tree, model=models[x])))
  ppcfits <-  foreach(i=trait.seq) %dopar% suppressWarnings(lapply(1:length(models), function(x) phylolm(simdat$ppc[,i]~1, phy=simdat$tree, model=models[x])))
  raw.aicw <- lapply(1:length(trait.seq), function(y) aicw(sapply(rawfits[[y]], function(x) x$aic))$w)
  pc.aicw <- lapply(1:length(trait.seq), function(y) aicw(sapply(pcfits[[y]], function(x) x$aic))$w)
  ppc.aicw <- lapply(1:length(trait.seq), function(y) aicw(sapply(ppcfits[[y]], function(x) x$aic))$w)
  aicw.table <- lapply(1:length(trait.seq), function(y) matrix(c(raw.aicw[[y]], pc.aicw[[y]], ppc.aicw[[y]]), nrow=1))
  aicw.table <- lapply(1:length(trait.seq), function(y){ colnames(aicw.table[[y]]) <- as.vector(t(outer(transforms, models, paste, sep="_"))); aicw.table[[y]]})
  aicw.table <- do.call(rbind, aicw.table)
  aicw.table <- as.data.frame(aicw.table)
  aicw.table <- data.frame("trait"=factor(paste("trait", trait.seq, sep=""),levels=paste("trait", trait.seq, sep="")), aicw.table)
  return(aicw.table)
}

## Special version of fitPCs that saves more output, like parameter estimates, for the OU model
fitPCOU <- function(simdat, trait.seq=1, models = c("BM", "OUfixedRoot", "EB")){
  transforms <- c("raw", "pc", "ppc")
  lower.bounds <- c(NULL, 10^-12, NULL)
  rawfits <- foreach(i=trait.seq) %dopar% suppressWarnings(lapply(1:length(models), function(x) phylolm(simdat$raw[,i]~1, phy=simdat$tree, model=models[x])))
  pcfits <-  foreach(i=trait.seq) %dopar% suppressWarnings(lapply(1:length(models), function(x) phylolm(simdat$pc[,i]~1, phy=simdat$tree, model=models[x])))
  ppcfits <-  foreach(i=trait.seq) %dopar% suppressWarnings(lapply(1:length(models), function(x) phylolm(simdat$ppc[,i]~1, phy=simdat$tree, model=models[x])))
  raw.aicw <- lapply(1:length(trait.seq), function(y) aicw(sapply(rawfits[[y]], function(x) x$aic))$w)
  pc.aicw <- lapply(1:length(trait.seq), function(y) aicw(sapply(pcfits[[y]], function(x) x$aic))$w)
  ppc.aicw <- lapply(1:length(trait.seq), function(y) aicw(sapply(ppcfits[[y]], function(x) x$aic))$w)
  aicw.table <- lapply(1:length(trait.seq), function(y) matrix(c(raw.aicw[[y]], pc.aicw[[y]], ppc.aicw[[y]]), nrow=1))
  aicw.table <- lapply(1:length(trait.seq), function(y){ colnames(aicw.table[[y]]) <- as.vector(t(outer(transforms, models, paste, sep="_"))); aicw.table[[y]]})
  aicw.table <- do.call(rbind, aicw.table)
  aicw.table <- as.data.frame(aicw.table)
  aicw.table <- data.frame("trait"=factor(paste("trait", trait.seq, sep=""),levels=paste("trait", trait.seq, sep="")), aicw.table)
  alpha <-  sapply(list(rawfits, pcfits, ppcfits), function(x) sapply(1:length(trait.seq), function (y) sapply(x[[y]][2], function(z) z$optpar)))
  sig2 <- sapply(list(rawfits, pcfits, ppcfits), function(x) sapply(1:length(trait.seq), function (y) sapply(x[[y]][2], function(z) z$sigma2)))
  colnames(alpha)<- colnames(sig2) <- c("raw", "pc", "ppc")
  alpha <- data.frame(trait=paste("trait", trait.seq, sep="."), as.data.frame(alpha))
  sig2 <- data.frame(trait=paste("trait", trait.seq, sep="."), as.data.frame(sig2))
  return(list(alpha=alpha, sig2=sig2, aicw.table=aicw.table))
}


## # Utility & processing functions
## Get contrasts from a simulated dataset
get.contrasts <- function(dat, simmodel, melt=TRUE){
  btimes <- lapply(1:length(dat), function(x) 1 - branching.times(dat[[x]]$tree))
  picraw <- lapply(1:length(dat), function(x) sapply(1:ncol(dat[[x]]$raw), function(y) pic(dat[[x]]$raw[,y], dat[[x]]$tree)))
  picpc <- lapply(1:length(dat), function(x) sapply(1:ncol(dat[[x]]$pc), function(y) pic(dat[[x]]$pc[,y], dat[[x]]$tree)))
  picppc <- lapply(1:length(dat), function(x) sapply(1:ncol(dat[[x]]$ppc), function(y) pic(dat[[x]]$ppc[,y], dat[[x]]$tree)))
  rawdf <- lapply(1:length(dat), function(x) data.frame(type="raw", times=btimes[[x]], trait=picraw[[x]]))
  pcdf <- lapply(1:length(dat), function(x) data.frame(type="pc", times=btimes[[x]], trait=picpc[[x]]))
  ppcdf <- lapply(1:length(dat), function(x) data.frame(type="ppc", times=btimes[[x]], trait=picppc[[x]]))
  df <- lapply(1:length(dat), function(x) rbind(rawdf[[x]], pcdf[[x]], ppcdf[[x]]))
  df <- lapply(1:length(dat), function(x) data.frame(df[[x]], rep=x))
  contdat <- do.call(rbind, df)
  contmelt <- melt(contdat, id=c("type", "times", "rep"))
  contmelt$value <- abs(contmelt$value)
  contmelt$simmodel <- simmodel
  return(contmelt)
}

## Get slopes of the node-height test
lmslope <- function(y, x){
    ff <- lm(y~x)
    ff$coef[2]
  }

## Get the disparity through time from a simulated dataset
get.dtt <- function(dat, simmodel){
  raw.bytrait <- lapply(1:ntraits, function(x) lapply(1:length(dat), function(y) dat[[y]]$raw[,x]))
  pc.bytrait <- lapply(1:ntraits, function(x) lapply(1:length(dat), function(y) dat[[y]]$pc[,x]))
  ppc.bytrait <- lapply(1:ntraits, function(x) lapply(1:length(dat), function(y) dat[[y]]$ppc[,x]))
  dttraw <- lapply(1:ntraits, function(y) lapply(1:length(dat), function(x) dtt(dat[[x]]$tree, setNames(raw.bytrait[[y]][[x]], dat[[x]]$tree$tip.label), plot=FALSE)))
  dispraw <- lapply(1:ntraits, function(x) unlist(lapply(1:length(dat), function(y) dttraw[[x]][[y]]$dtt)))
  times <- unlist(lapply(1:ntraits, function(x) unlist(lapply(1:length(dat), function(y) dttraw[[x]][[y]]$times))))
  
  dttpc <- lapply(1:ntraits, function(y) lapply(1:length(dat), function(x) dtt(dat[[x]]$tree, setNames(pc.bytrait[[y]][[x]], dat[[x]]$tree$tip.label), , plot=FALSE)))
  disppc <- lapply(1:ntraits, function(x) unlist(lapply(1:length(dat), function(y) dttpc[[x]][[y]]$dtt)))
  dttppc <- lapply(1:ntraits, function(y) lapply(1:length(dat), function(x) dtt(dat[[x]]$tree, setNames(ppc.bytrait[[y]][[x]], dat[[x]]$tree$tip.label), plot=FALSE)))
  dispppc <- lapply(1:ntraits, function(x) unlist(lapply(1:length(dat), function(y) dttppc[[x]][[y]]$dtt)))
  disp <- list(do.call(cbind, dispraw), do.call(cbind, disppc), do.call(cbind, dispppc))
  types <- c("raw", "pc", "ppc")
  dispdf <- lapply(1:3, function(x) data.frame(type=types[x], times=times, trait=disp[[x]]))
  dispdf <- do.call(rbind, dispdf)
  dispmelt <- melt(dispdf, id=c('type', 'times'))
  dispmelt$simmodel <- simmodel
  return(dispmelt)
}

## Get the parameter estimates for the OU simulated datasets
get.parsOU <- function(res){
  alphares <- lapply(res, function(x) x[['alpha']])
  alphares <- do.call(rbind, alphares)
  sig2res <- lapply(res, function(x) x[['sig2']])
  sig2res <- do.call(rbind, sig2res)
  sig2melt <- melt(sig2res, id.vars="trait")
  alphamelt <- melt(alphares, id.vars="trait")
  parsdf <- alphamelt
  colnames(parsdf)[3] <- "alpha"
  parsdf$sig2 <- sig2melt$value
  parsdf$halflife <- log(2)/parsdf$alpha
  #head(parsdf)
  #parsmelt <- melt(parsdf, id.vars=c("trait", "variable"))
  #colnames(parsmelt)[3] = "parameter"
  #parsmelt$traitno <- trait.seq
  return(parsdf)
}

## Get the average AIC weights for the OU fit results
get.meanAICws <- function(res){
  aicws <- lapply(res, function(x) x$aicw.table)
  aicws <- do.call(rbind, aicws)
  meanaicws <- apply(aicws[2:ncol(aicws)], 2, function(x) tapply(x, aicws[,1], mean))
  meanaicws <- as.data.frame(meanaicws)
  meanaicws$trait <- rownames(meanaicws)
  meanaicws <- melt(meanaicws, id.vars="trait")
  meanaicws$type <- factor(sapply(strsplit(as.character(meanaicws$variable), split="_"), function(x) x[[1]]))
  meanaicws$simmodel <- factor(sapply(strsplit(as.character(meanaicws$variable), split="_"), function(x) x[[2]]))
  meanaicws <- meanaicws[,-2]
  meanaicws <- meanaicws[,c(1,3,4,2)]
}


## Function for simulating uncorrelated ACDC data
sim.tree.pcs.ind.acdc <- function(ntips, traits, asd, sig2, lambda=0.1, mu=0, cor=TRUE, ...){
  tree <- pbtree(b=lambda, d=mu, n=ntips)
  tree$edge.length <- tree$edge.length/max(branching.times(tree))
  a <- rnorm(traits, 0, asd)
  trees <- lapply(a, function(x) geiger:::rescale.phylo(tree, model="EB", x))
  #outree <- rescale(tree,alpha=log(2),model="OU")
  X <- sapply(1:traits, function(x) fastBM(trees[[x]],sig2=sig2, nsim=1))
  if(cor){
    pc <- princomp(X, cor=TRUE)
    ppc <- phyl.pca(tree,X, mode="corr")
  } else {
    pc <- princomp(X)
    ppc <- phyl.pca(tree,X)  
  }
  return(list(tree=tree, raw=X, pc=pc$scores, ppc=ppc$S, pcall=pc, ppcall=ppc, a=a))
}

## Prepare ACDC dataset for plotting
prepare.acdc <- function(dat){
  aLL <- lapply(dat, function(x) cbind(x$a, x$pcall$loadings, x$ppcall$L))
  aLL.slopes <- lapply(aLL, function(x) sapply(2:ncol(x), function(y) lm(abs(x[,y]) ~ x[,1])$coeff[2]))
  aLL.slopes <- do.call(rbind, aLL.slopes)
  colnames(aLL.slopes) <- c(paste("PC", 1:20, sep=""), paste("PPC", 1:20, sep=""))
  aLL <- do.call(rbind, aLL)
  colnames(aLL) <- c("a", paste("PC", 1:20, sep=""), paste("PPC", 1:20, sep=""))
  aLL <- as.data.frame(aLL)
  aLL.melt <- melt(aLL, id.vars="a")
  return(list(melt=aLL.melt, slopes=aLL.slopes))
}

