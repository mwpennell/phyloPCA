require(bayou)
require(MASS)
require(corrplot)
require(phytools)
require(foreach)
require(doMC)
data(chelonia)
phy <- chelonia$phy
##Example dataset
tree <- chelonia$phy
dat <- chelonia$dat
tree <- reorder(tree, "postorder")
dat <- dat[tree$tip.label]
##Define prior function. This will be the trickiest part. I've chosen log-normal priors for alpha and sigma^2, 
##normal prior for theta, a conditional Poisson for the number of shifts, and "dsb" controls how many shifts can be per branch (either 0, 1 or Inf) and the probability of a shift being on that branch
prior <- make.prior(tree, dists=list(dalpha="dlnorm", dsig2="dlnorm",dsb="dsb", dk="cdpois", dtheta="dnorm"), 
                    param=list(dalpha=list(meanlog=-1, sdlog=1), dsig2=list(meanlog=-1, sdlog=2), 
                               dk=list(lambda=15, kmax=200), dsb=list(bmax=1,prob=1), dtheta=list(mean=mean(dat), sd=2)))

simpars <- priorSim(prior, tree, plot=FALSE, nsim=10)
simdats <- lapply(simpars$pars, function(x) sapply(1:30, function(y) dataSim(x, model="OU", tree, phenogram=FALSE)$dat))

ppcs <- lapply(simdats, function(x) phyl.pca(tree, x, method="BM"))

ouwiepars.raw <- lapply(1:30, function(y) bayou2OUwie(simpars$pars[[3]], tree, simdats[[3]][,y]))
ouwiepars.ppc <- lapply(1:30, function(y) bayou2OUwie(simpars$pars[[3]], tree, ppcs[[3]]$S[,y]))

registerDoMC(cores=8)
require(OUwie)
res.ppc <- foreach(i=c(1, 2, 5, 10, 15, 20, 25, 30)) %dopar% {
  OUwie(ouwiepars.ppc[[i]]$tree, ouwiepars.ppc[[i]]$dat, model="OUM")
}

res.raw <- foreach(i=c(1, 2, 5, 10, 15, 20, 25, 30)) %dopar% {
  OUwie(ouwiepars.raw[[i]]$tree, ouwiepars.raw[[i]]$dat, model="OUM")
}

thalf <- sapply(res.ppc, function(x) log(2)/x$solution[1,1])
sig2 <- sapply(res.ppc, function(x) x$solution[2,1])


plot(c(1, 2, 5, 10, 15, 20, 25, 30), thalf)
abline(h=log(2)/simpars$pars[[3]]$alpha)
plot(c(1, 2, 5, 10, 15, 20, 25, 30), sig2)

par(mfrow=c(4,2))
for(i in c(1,2,10,30)){
tmp <- ppcs[[3]]$S[,i]
names(tmp) <- tree$tip.label
phenogram(tree, tmp, ftype="off")
#tmp <- ppcs[[3]]$S[,1]
#names(tmp) <- tree$tip.label
phenogram(tree, simdats[[3]][,i], ftype="off")
}


str(res.ppc[[1]])



