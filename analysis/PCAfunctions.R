sim.tree.pcs.ind <- function(ntips, traits, sig2dist=rexp, lambda=0.1, mu=0, ...){
  tree <- sim.bd.taxa(ntips,1,lambda=0.1,mu=0)[[1]][[1]]
  tree$edge.length <- tree$edge.length/max(branching.times(tree))
  #outree <- rescale(tree,alpha=log(2),model="OU")
  X <- sapply(1:traits,function(x) fastBM(tree,sig2=sig2dist(1, ...), nsim=1))
  pc <- princomp(X)
  ppc <- phyl.pca(tree,X)
  return(list(tree=tree, raw=X, pc=pc$scores, ppc=ppc$S, pcall=pc, ppcall=ppc))
}

sim.tree.pcs.ind.ou <- function(ntips, traits, alpha, sig2, lambda=0.1, mu=0, ...){
  tree <- sim.bd.taxa(ntips,1,lambda=0.1,mu=0)[[1]][[1]]
  tree$edge.length <- tree$edge.length/max(branching.times(tree))
  tr <- rescale(tree, model="OU", alpha)
  #outree <- rescale(tree,alpha=log(2),model="OU")
  X <- sapply(1:traits, function(x) fastBM(tr, sig2=sig2, nsim=1))
  pc <- princomp(X)
  ppc <- phyl.pca(tree,X)
  return(list(tree=tree, raw=X, pc=pc$scores, ppc=ppc$S, pcall=pc, ppcall=ppc))
}

sim.tree.pcs.ind.eb <- function(ntips, traits, a, sig2, lambda=0.1, mu=0, ...){
  tree <- sim.bd.taxa(ntips,1,lambda=0.1,mu=0)[[1]][[1]]
  tree$edge.length <- tree$edge.length/max(branching.times(tree))
  tr <- rescale(tree, model="EB", a)
  #outree <- rescale(tree,alpha=log(2),model="OU")
  X <- sapply(1:traits, function(x) fastBM(tr,sig2=sig2, nsim=1))
  pc <- princomp(X)
  ppc <- phyl.pca(tree,X)
  return(list(tree=tree, raw=X, pc=pc$scores, ppc=ppc$S, pcall=pc, ppcall=ppc))
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

get.contrasts <- function(dat){
  btimes <- lapply(1:nsims, function(x) 1 - branching.times(dat[[x]]$tree))
  picraw <- lapply(1:nsims, function(x) sapply(1:ntraits, function(y) pic(dat[[x]]$raw[,y], dat[[x]]$tree)))
  picpc <- lapply(1:nsims, function(x) sapply(1:ntraits, function(y) pic(dat[[x]]$pc[,y], dat[[x]]$tree)))
  picppc <- lapply(1:nsims, function(x) sapply(1:ntraits, function(y) pic(dat[[x]]$ppc[,y], dat[[x]]$tree)))
  rawdf <- lapply(1:nsims, function(x) data.frame(type="raw", times=btimes[[x]], trait=picraw[[x]]))
  pcdf <- lapply(1:nsims, function(x) data.frame(type="pc", times=btimes[[x]], trait=picpc[[x]]))
  ppcdf <- lapply(1:nsims, function(x) data.frame(type="ppc", times=btimes[[x]], trait=picppc[[x]]))
  df <- lapply(1:nsims, function(x) rbind(rawdf[[x]], pcdf[[x]], ppcdf[[x]]))
}

get.dtt <- function(dat){
  raw.bytrait <- lapply(1:ntraits, function(x) lapply(1:nsims, function(y) dat[[y]]$raw[,x]))
  pc.bytrait <- lapply(1:ntraits, function(x) lapply(1:nsims, function(y) dat[[y]]$pc[,x]))
  ppc.bytrait <- lapply(1:ntraits, function(x) lapply(1:nsims, function(y) dat[[y]]$ppc[,x]))
  dttraw <- lapply(1:ntraits, function(y) lapply(1:nsims, function(x) dtt(dat[[x]]$tree, raw.bytrait[[y]][[x]], plot=FALSE)))
  dispraw <- lapply(1:ntraits, function(x) unlist(lapply(1:nsims, function(y) dttraw[[x]][[y]]$dtt)))
  times <- unlist(lapply(1:ntraits, function(x) unlist(lapply(1:nsims, function(y) dttraw[[x]][[y]]$times))))
  
  dttpc <- lapply(1:ntraits, function(y) lapply(1:nsims, function(x) dtt(dat[[x]]$tree, pc.bytrait[[y]][[x]], plot=FALSE)))
  disppc <- lapply(1:ntraits, function(x) unlist(lapply(1:nsims, function(y) dttpc[[x]][[y]]$dtt)))
  dttppc <- lapply(1:ntraits, function(y) lapply(1:nsims, function(x) dtt(dat[[x]]$tree, ppc.bytrait[[y]][[x]], plot=FALSE)))
  dispppc <- lapply(1:ntraits, function(x) unlist(lapply(1:nsims, function(y) dttppc[[x]][[y]]$dtt)))
  disp <- list(do.call(cbind, dispraw), do.call(cbind, disppc), do.call(cbind, dispppc))
  types <- c("raw", "pc", "ppc")
  dispdf <- lapply(1:3, function(x) data.frame(type=types[x], times=times, trait=disp[[x]]))
  dispdf <- do.call(rbind, dispdf)
  dispmelt <- melt(dispdf, id=c('type', 'times'))
  return(dispmelt)
}

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

simandfitPCs.ind.ou <- function(ntips, traits, trait.seq=1, alpha, sig2, models = c("BM", "OUrandomRoot", "EB")){
  simdat <- sim.tree.pcs.ind.ou(ntips, traits, alpha, sig2, lambda=1/100)
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

simandfitPCs.ind.eb <- function(ntips, traits, trait.seq=1, a, sig2, models = c("BM", "OUrandomRoot", "EB")){
  simdat <- sim.tree.pcs.ind.eb(ntips, traits, a, sig2, lambda=1/100)
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

makeTransparent <- function (someColor, alpha = 100){
  newColor <- col2rgb(someColor)
  apply(newColor, 2, function(curcoldata) {
    rgb(red = curcoldata[1], green = curcoldata[2], blue = curcoldata[3], 
        alpha = alpha, maxColorValue = 255)
  })
}

lines.loess <- function(x, y, ...){
  df <- as.data.frame(x=x, y=y)
  tmp <- loess(y~x, data=df, )
  pred <- predict(tmp, data.frame(x=seq(0,1,0.01)), se=FALSE)
  lines(seq(0,1,0.01), pred, ...)
}