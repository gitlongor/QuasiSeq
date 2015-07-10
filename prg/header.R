base.dir=file.path('C:/Users/Long/XXXXXXXXXXXXXXXXXXXXXXXXX')
out.dir=file.path(base.dir, 'out')
prg.dir=file.path(base.dir, 'prg')
dat.dir=file.path(base.dir, 'dat')
img.dir=file.path(base.dir, 'img')

if(TRUE){
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

prg.name=paste('XXXXXXXXXXXXXXXXXXXXXXXXX', Sys.Date(), sep='_')
img.name=file.path(img.dir, paste(prg.name, 'RData', sep='.'))
out.file=function(suffix, subdir) if(missing(subdir)) file.path(out.dir, paste(prg.name, suffix, sep='_')) else file.path(out.dir, subdir, paste(prg.name, suffix, sep='_')) 

setwd(out.dir)

options(contrasts=c('contr.sum', 'contr.sum'))

.sessionInfo=list(sessionInfo=sessionInfo(), Sys.info=Sys.info(),
  Sys.getenv=Sys.getenv(names=TRUE), capabilities=capabilities(),
  options=options(), RNGkind=RNGkind(), .libPaths=.libPaths(), 
  date=date(), wd=getwd(),.Random.seed=.Random.seed)


####################################### END OF HEADER ################################################
	