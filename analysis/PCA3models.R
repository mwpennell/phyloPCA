require(phytools)
require(geiger)
require(phylolm)
require(TreeSim)
ntips <- 100
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
is.ultrametric(tree)
try(plot(tree))
plot(tree)
X=fastBM(tree,nsim=100)
pc <- princomp(X)
plot(pc$scores[,1],pc$scores[,2])


akaike.weights <- function(AIC){
  D <- AIC-min(AIC)
  w <- exp(-D/2)/sum(exp(-D/2))
  return(w)
}


ntips <- 100
traits <- 100
nsims <- 10
OUT <- list("BM"=list(PC=NULL,RAW=NULL),"OU"=list("PC"=NULL,"RAW"=NULL),"EB"=list("PC"=NULL,"RAW"=NULL))
OUT.AICw <- matrix(NA,ncol=12,nrow=nsims)
colnames(OUT.AICw) <- c("BMpc","OUpc","EBpc","BMr","OUr","EBr","BMppc","OUppc","EBppc","BMsc","OUsc","EBsc")
pc.i <- 1
for(i in 1:nsims){
  #tree <- sim.bd.taxa(ntips,1,lambda=0.1,mu=0)[[1]][[1]]
  #tree$edge.length <- tree$edge.length/max(nodeHeights(tree))
  #outree <- rescale(tree,alpha=log(2),model="OU")
  X <- sapply(1:traits,function(x) fastBM(tree,sig2=rexp(1,1/100),nsim=1))
  pc <- princomp(X)
  ppc <- phyl.pca(tree,X)
  sc <- scale(X,center=FALSE)
  scm <- apply(sc,1,mean)
                                        #apply(X,1,function(x) (abs(x)-min(x))/(max(x)-min(x)))
  pcdat <- data.frame(pc$scores)
  ppcdat <- data.frame(ppc$S)
  fitBM <- phylolm(pcdat[,pc.i]~1,data=pcdat,phy=tree,model="BM")
  fitOU <- phylolm(pcdat[,pc.i]~1,data=pcdat,phy=tree,model="OUfixedRoot",lower.bound=("alpha"=10^-12))
  fitEB <- phylolm(pcdat[,pc.i]~1,data=pcdat,phy=tree,model="EB")
  fitBMr <- phylolm(X[,1]~1,phy=tree,model="BM")
  fitOUr <- phylolm(X[,1]~1,phy=tree,model="OUfixedRoot",lower.bound=("alpha"=10^-12))
  fitEBr <- phylolm(X[,1]~1,phy=tree,model="EB")
  fitBMppc <- phylolm(ppcdat[,pc.i]~1,dat=ppcdat,phy=tree,model="BM")
  fitOUppc <- phylolm(ppcdat[,pc.i]~1,dat=ppcdat,phy=tree,model="OUfixedRoot",lower.bound=("alpha"=10^-12))
  fitEBppc <- phylolm(ppcdat[,pc.i]~1,dat=ppcdat,phy=tree,model="EB")
  fitBMsc <- phylolm(scm~1,phy=tree,model="BM")
  fitOUsc <- phylolm(scm~1,phy=tree,model="OUfixedRoot",lower.bound=("alpha"=10^-12))
  fitEBsc <- phylolm(scm~1,phy=tree,model="EB")
  OUT$BM$PC[[i]] <- fitBM
  OUT$OU$PC[[i]] <- fitOU
  OUT$EB$PC[[i]] <- fitEB
  OUT$BM$RAW[[i]] <- fitBMr
  OUT$OU$RAW[[i]] <- fitOUr
  OUT$EB$RAW[[i]] <- fitEBr
  AIC <- sapply(list(fitBM,fitOU,fitEB),function(x) x$aic)
  AICr <- sapply(list(fitBMr,fitOUr,fitEBr),function(x) x$aic)
  AICppc <- sapply(list(fitBMppc,fitOUppc,fitEBppc),function(x) x$aic)
  AICsc <- sapply(list(fitBMsc,fitOUsc,fitEBsc),function(x) x$aic)
  OUT.AICw[i,1:3] <- akaike.weights(AIC)
  OUT.AICw[i,4:6] <- akaike.weights(AICr)
  OUT.AICw[i,7:9] <- akaike.weights(AICppc)
  OUT.AICw[i,10:12] <- akaike.weights(AICsc)
  print(OUT.AICw[i,])
}
boxplot(OUT.AICw, range=50)

plot(OUT[,1],ylim=c(0,5))
plot(1:ncol(OUT),apply(OUT,2,mean))







require(geiger)
require(TreeSim)
ntaxa <- 64
tree <- sim.bd.taxa(ntaxa,1,lambda=0.1,mu=0)[[1]][[1]]
X <- sim.char(tree,1,nsim=100,model="BM")
X <- data.frame(X)
pc <- princomp(X)
pc1 <- pc$scores[,1]
X.bar <- apply(X,1,mean)

bm.fit <- fitContinuous(tree,pc1,model="BM")
ou.fit <- fitContinuous(tree,pc1,model="OU")
eb.fit <- fitContinuous(tree,pc1,model="EB")

bm.fit2 <- fitContinuous(tree,X.bar,model="BM")
ou.fit2 <- fitContinuous(tree,X.bar,model="OU")
eb.fit2 <- fitContinuous(tree,X.bar,model="EB")


aic.scores <- c(bm.fit$opt$aic,ou.fit$opt$aic,eb.fit$opt$aic)
aic.scores2 <- c(bm.fit2$opt$aic,ou.fit2$opt$aic,eb.fit2$opt$aic)
