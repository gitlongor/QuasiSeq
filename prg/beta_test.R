base.dir=file.path('D:/work/dev/Rstudio/QuasiSeq')
out.dir=file.path(base.dir, 'out')
prg.dir=file.path(base.dir, 'prg')
dat.dir=file.path(base.dir, 'dat')
img.dir=file.path(base.dir, 'img')

if(FALSE){
	library(parallel)
	RNGkind("L'Ecuyer-CMRG")
	cl=makeCluster(8L)
	clusterEvalQ(cl, {
		RNGkind("L'Ecuyer-CMRG"); 
		options(contrasts=c('contr.sum', 'contr.sum'))
	})
	clusterSetRNGStream(cl , 9992722L)
}
set.seed(9992722L)

prg.name=paste('beta_test', Sys.Date(), sep='_')
img.name=file.path(img.dir, paste(prg.name, 'RData', sep='.'))
out.file=function(suffix, subdir) if(missing(subdir)) file.path(out.dir, paste(prg.name, suffix, sep='_')) else file.path(out.dir, subdir, paste(prg.name, suffix, sep='_')) 

setwd(out.dir)

options(contrasts=c('contr.sum', 'contr.sum'))

.sessionInfo=list(sessionInfo=sessionInfo(), Sys.info=Sys.info(),
  Sys.getenv=Sys.getenv(names=TRUE), capabilities=capabilities(),
  options=options(), RNGkind=RNGkind(), .libPaths=.libPaths(), 
  date=date(), wd=getwd(),.Random.seed=.Random.seed)


####################################### END OF HEADER ################################################

load(file.path(dat.dir, 'testdata.RData'))

### Load Bioconductor Packages
require(edgeR)

### Load CRAN Packages
require(QuasiSeq)

source(file.path(prg.dir, "bnbReg.R"))

    ### Create Design Matrices
    design.list <- list(X)
    
    ### Compute Normalization Factors
    nf.nb <- calcNormFactors(counts.nb) * apply(counts.nb, 2, sum)
    
    ### fit for SimSeq data
    fit <- QL.fit(counts.nb, design.list = design.list, log.offset = log(nf.nb), Model = "NegBin", print.progress = FALSE)
				  

phis=pmax(fit$phi.hat.dev,1+1e-3)
taus=fit$NB.disp
offs=log(nf.nb)

G=nrow(counts.nb)
p=ncol(X)

betahats=matrix(NA_real_, G, p)

for(g in seq(G)) {
	tmp=Mlebeta(counts.nb[g,], x=X,o=offs, phi=phis[g], tau=taus[g], start=fit$coefficients[g,])
	betahats[g,]=tmp$bet
}

plot(betahats[,1], fit$coef[,1]); abline(0,1)
plot(betahats[,2], fit$coef[,2]); abline(0,1)


trt.sums=cbind(rowSums(counts.nb[,1:5]), rowSums(counts.nb[,1:5+5L]))
non0.idx=which(apply(trt.sums,1L, min)>0L)


plot(betahats[non0.idx,1], fit$coef[non0.idx,1]); abline(0,1)
plot(betahats[non0.idx,2], fit$coef[non0.idx,2]); abline(0,1)

plot(betahats[-non0.idx,1], fit$coef[-non0.idx,1]); abline(0,1)
plot(betahats[-non0.idx,2], fit$coef[-non0.idx,2]); abline(0,1)
