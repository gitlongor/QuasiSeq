base.dir=file.path('C:/Users/Tharshan-PC/Desktop/Research/QuasiSeq')
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

prg.name=paste('mle_test', Sys.Date(), sep='_')
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

source(file.path(prg.dir, "bnbRegMle.R"))

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
tauPhihats=matrix(NA_real_, G, 2L)
niters=numeric(G)

for(g in seq(G)) {
	tmp=bnbRegMle(counts.nb[g,], x=X,o=offs, start=c(fit$coefficients[g,], phis[g], taus[g]) , iter.max=5e3L)
	betahats[g,]=tmp$bet
	tauPhihats[g,]=c(tmp$phi, tmp$tau)
	niters[g]=tmp$iter
}

table(niters)/G

plot(betahats[,1], fit$coef[,1]); abline(0,1)
plot(betahats[,2], fit$coef[,2]); abline(0,1)

plot(log10(tauPhihats[,1]), log10(phis)); abline(0,1)
cor(log10(tauPhihats[,1]), log10(phis))
plot(log10(tauPhihats[,2L]+1), log10(taus+1)); abline(0,1)
cor(log10(tauPhihats[,2L]+1), log10(taus+1))

quasi.means=exp(sweep(tcrossprod(fit$coefficients, X),2L, offs, "+"))
quasi.vars=phis*quasi.means*(1+taus*quasi.means)

bnb.means=exp(sweep(tcrossprod(betahats, X),2L, offs, "+"))
bnb.vars=tauPhihats[,1L]*bnb.means*(1+tauPhihats[,2L]*bnb.means)

