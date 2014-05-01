sim.tree.pcs.ind <- function(ntips, traits, sig2dist=rexp, lambda=0.1, mu=0, ...){
  tree <- sim.bd.taxa(ntips,1,lambda=0.1,mu=0)[[1]][[1]]
  tree$edge.length <- tree$edge.length/max(branching.times(tree))
  #outree <- rescale(tree,alpha=log(2),model="OU")
  X <- sapply(1:traits,function(x) fastBM(tree,sig2=sig2dist(1, ...), nsim=1))
  pc <- princomp(X)
  ppc <- phyl.pca(tree,X)
  return(list(tree=tree, raw=X, pc=pc$scores, ppc=ppc$scores, pcall=pc, ppcall=ppc))
}

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


sim.tree.pcs.mv <- function(ntips, traits, sig2dist=rexp, lambda=0.1, mu=0, ...){
  tree <- sim.bd.taxa(ntips,1,lambda=0.1,mu=0)[[1]][[1]]
  tree$edge.length <- tree$edge.length/max(branching.times(tree))
  Sigma <- vcv.phylo(tree)
  #rates <- c(1,10,100,1000)#sig2dist(traits)#, ...)
  #R <- matrix(0, ncol=traits, nrow=traits)
  #corr.n <- (traits*traits-traits)/2
  #diag(R) <- 1
  #corr.R <- c(0.8, 0.8, -0.5, 0.8, -0.5, -0.5)#rep(0.5, corr.n)#runif(length(R[lower.tri(R)]), 0, 1)
  #R[lower.tri(R)] <- corr.R
  #R[upper.tri(R)] <- (t(R))[upper.tri(R)]
  #R <- (sqrt(rates) %*% t(sqrt(rates))) * R
  R <- Posdef(traits, ev=sig2dist(traits, ...))
  simdat <- t(matrix(mvrnorm(1, rep(0, traits*ntips), R %x% Sigma), nrow=traits, ncol=ntips, byrow=TRUE))
  pc <- princomp(simdat)
  ppc <- phyl.pca(tree,simdat)
  row.names(simdat) <- tree$tip.label
  rownames(pc$scores) <- tree$tip.label
  rownames(ppc$S) <- tree$tip.label
  return(list(tree=tree, raw=simdat, pc=pc$scores, ppc=ppc$S, pcall=pc, ppcall=ppc))
}

simandfitPCs <- function(ntips, traits, trait.seq=1, models = c("BM", "OUrandomRoot", "EB")){
  simdat <- sim.tree.pcs.mv(ntips, traits, sig2dist=rexp, lambda=1/100)
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
  names(aicw.table) <- paste("trait", trait.seq, sep="")
  return(aicw.table)
}