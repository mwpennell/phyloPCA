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
  return(list(tree=tree, raw=simdat, pc=pc$scores, ppc=ppc$S, pcall=pc, ppcall=ppc, R=R))
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

#Simulate tree and data for pgls using mvBM. Removes first trait from the pca.
sim.tree.pcs.pgls <- function(ntips, traits, foc2trans="BM", rescalepar=1, sig2dist=rexp, me=0, lambda=0.1, mu=0, ...){
  #tree <- sim.bd.taxa(ntips,1,lambda=0.1,mu=0)[[1]][[1]]
  #tree$edge.length <- tree$edge.length/max(branching.times(tree))
  tree <- laddertree(ntips)
  Sigma <- vcv.phylo(tree)
  R <- Posdef(traits, ev=sig2dist(traits, 1/100))
  R[1, 2:ncol(R)] <- 0
  R[2:nrow(R), 1] <- 0
  focdat <- fastBM(rescale(tree, model=foc2trans, rescalepar), a=0, sig2=R[1,1])+rnorm(ntips, 0, me)
  simdat <- t(matrix(mvrnorm(1, rep(0, traits*ntips), R %x% Sigma), nrow=traits, ncol=ntips, byrow=TRUE))
  diag(simdat) <- diag(simdat)+rnorm(length(diag(simdat)), 0, me)
  pc <- princomp(simdat[,-1])
  ppc <- phyl.pca(tree,simdat[,-1])
  row.names(simdat) <- tree$tip.label
  rownames(pc$scores) <- tree$tip.label
  rownames(ppc$S) <- tree$tip.label
  return(list(tree=tree, raw=as.data.frame(cbind(focdat, simdat[,-1])), pc=as.data.frame(cbind(focdat, pc$scores)), ppc=as.data.frame(cbind(focdat, ppc$S)), pcall=pc, ppcall=ppc, R=R))
}

simandfitPGLS <- function(ntips, traits, trait.seq=1, models = c("BM", "OUrandomRoot", "EB"), foc2trans="BM", rescalepar=1, me=0){
  simdat <- sim.tree.pcs.pgls(ntips, traits, sig2dist=rexp, foc2trans=foc2trans, rescalepar=rescalepar)
  transforms <- c("raw", "pc", "ppc")
  lower.bounds <- c(NULL, 10^-12, NULL)
  corBM <- corBrownian(1, simdat$tree)
  rawfits <- foreach(i=trait.seq) %dopar% {fmla = as.formula(paste(colnames(simdat$raw)[1], "~", colnames(simdat$raw)[i])); 
                                           gls(fmla, correlation=corBM, data=simdat$raw)}
  pcfits <- foreach(i=trait.seq) %dopar% {fmla = as.formula(paste(colnames(simdat$pc)[1], "~", colnames(simdat$pc)[i])); 
                                           gls(fmla, correlation=corBM, data=simdat$pc)}
  ppcfits <- foreach(i=trait.seq) %dopar% {fmla = as.formula(paste(colnames(simdat$ppc)[1], "~", colnames(simdat$ppc)[i])); 
                                           gls(fmla, correlation=corBM, data=simdat$ppc)}
  rawcoef <- t(sapply(1:length(trait.seq), function(y) c(rawfits[[y]]$coeff[2], summary(rawfits[[y]])$tTable[2,4])))
  pccoef <- t(sapply(1:length(trait.seq), function(y) c(pcfits[[y]]$coeff[2], summary(pcfits[[y]])$tTable[2,4])))
  ppccoef <- t(sapply(1:length(trait.seq), function(y) c(ppcfits[[y]]$coeff[2], summary(ppcfits[[y]])$tTable[2,4])))
  #rawsigma <- sapply(1:length(trait.seq), function(y) rawfits[[y]]$sigma)
  #pcsigma <- sapply(1:length(trait.seq), function(y) pcfits[[y]]$sigma)
  #ppcsigma <- sapply(1:length(trait.seq), function(y) ppcfits[[y]]$sigma)
  return(list(raw=rawcoef, pc=pccoef, ppc=ppccoef, simdat=simdat))
}

phenogramn <- function(tree, dat){
  names(dat) <- tree$tip.label
  phenogram(tree, dat)
}

laddertree <- function(ntips){
  node <- ntips+1
  L <- 2*ntips-2
  edge=matrix(ncol=2,nrow=2*ntips-2)
  edge[,1]=sort(rep((ntips+1):(2*ntips-1),2))
  edge[seq(1,2*ntips-2,2),2]=1:(ntips-1)
  edge[seq(2,2*ntips-2,2),2]=(ntips+2):(2*ntips)
  edge[2*ntips-2,2]=ntips
  edge.length=rep(1/ntips,L)
  edge.length[c(seq(1,L,2),L)]=1-1/ntips*0:(ntips-1)
  edge.length[2*ntips-2] <- 0.02
  tree=list('edge'=edge,'tip.label'=1:ntips,'edge.length'=edge.length,'Nnode'=ntips-1)
  class(tree)='phylo'
  tree$edge.length <- tree$edge.length/(max(branching.times(tree)))
  return(tree)
}

#cherrytree <- function()
#cherrymaker <- function(n1){
#  paste("(", n1, ":1,", n1+1 ,":1)",sep="")
#}
#cherries <- sapply(seq(1, 64, 2), cherrymaker)
##cherrymaker2 <- function(n1, cherries=cherries){
#  paste("(", cherries[n1], ":1,", cherries[n1+1] ,":1)",sep="")
#}
#cherry2 <- sapply(seq(1, 32, 2), cherrymaker2, cherries=cherries)
#cherry3 <- sapply(seq(1, 16, 2), cherrymaker2, cherries=cherry2)
#cherry4 <- sapply(seq(1, 8, 2), cherrymaker2, cherries=cherry3)
#cherry5 <- sapply(seq(1, 4, 2), cherrymaker2, cherries=cherry4)
#cherry6 <- sapply(1, cherrymaker2, cherries=cherry5)
#cherry6 <- paste(cherry6, ";", sep="")
#write(cherry6, file="balancedtree.tre")
#setwd("~/BayesOUFit/shiftestimation")
#tree <- read.tree("balancedtree.tre")
##tree <- reorder(tree, "postorder")
#tree$tip.label <- sapply(tree$tip.label, function(x) paste("t", x, sep=""))
#tree$edge.length <- tree$edge.length/max(nodeHeights(tree))
