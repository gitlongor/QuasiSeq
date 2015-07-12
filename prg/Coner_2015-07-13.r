base.dir=file.path('C:/Users/Tharshan-PC/Desktop/Research/QuasiSeq')
out.dir=file.path(base.dir, 'out')
prg.dir=file.path(base.dir, 'prg')
dat.dir=file.path(base.dir, 'dat')
img.dir=file.path(base.dir, 'img')

if(F){
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

prg.name=paste('Conver', Sys.Date(), sep='_')
img.name=file.path(img.dir, paste(prg.name, 'RData', sep='.'))
out.file=function(suffix, subdir) if(missing(subdir)) file.path(out.dir, paste(prg.name, suffix, sep='_')) else file.path(out.dir, subdir, paste(prg.name, suffix, sep='_')) 

setwd(out.dir)

options(contrasts=c('contr.sum', 'contr.sum'))

.sessionInfo=list(sessionInfo=sessionInfo(), Sys.info=Sys.info(),
  Sys.getenv=Sys.getenv(names=TRUE), capabilities=capabilities(),
  options=options(), RNGkind=RNGkind(), .libPaths=.libPaths(), 
  date=date(), wd=getwd(),.Random.seed=.Random.seed)


####################################### END OF HEADER ################################################


#####Graphical comparison QLNB and LBNB

my_QLNB<-function(y, mu,phi,tau) {
	exp(1/phi)*((1/(1+tau*mu))^(1/tau)*(tau*mu/(tau*mu+1))^y)
}

my_LBNB<-function(y, mu,phi,tau){
	term1 = mu*(tau*phi+1)/(phi-1)
	term2 = (tau*(2*phi-1)+1)/(tau*(phi-1))
	numerator = lbeta(term1+term2,y+1/tau)
	term3 = lbeta(term2,1/tau)
	ans = exp(		
		numerator-lbeta(y, term1)-term3
	)/y
	ans[y==0]=exp(numerator-term3)
	ans
}

		
##################END####################################

#This is the function to find exp(h(y,tau,pi))=H,where 
#exp(h(y,tau=t,pi))=integration(f(QLNB)-f(LBNB)) over the mean values
#f(QLNB)=exp(1/pi)*((1/(1+t*x))^(1/t)*(t*x/(t*x+1))^y)
#f(LBNB)=(exp(lgamma(y+x*(t*pi+1)/(pi-1))+lbeta((x*(t*pi+1))/(pi-1)
#             +(t*(2*pi-1)+1)/(t*(pi-1)),y+1/t)-log(factorial(y))-lgamma(x*(t*pi+1)/(pi-1))
#             -lbeta((t*(2*pi-1)+1)/t*(pi-1),1/t)))

get.h<-function(phi,y,tau){
	f = function(mu) my_QLNB(y, mu, phi, tau)
	fm = optimize(f,c(0,1e6), maximum=TRUE)$objective
	f = function(mu) my_LBNB(y, mu, phi, tau)
	fm2 = optimize(f,c(0,1e6), maximum=TRUE)$objective
	fm2/fm
}

get.h(2,10,0.2)

logLikQNB=function(y, mu, phi, tau)
{
	h=get.h(phi, y, tau)
	my_QLNB(y, mu, phi, tau) * h
}
logLikQNB = Vectorize(logLikQNB, c('phi','y','tau'))


curve(log(logLikQNB(y=10, mu=x, phi=1.1, tau=.5)), 8, 12)
curve(log(my_LBNB(y=10, mu=x, phi=1.1, tau=.5)), 8, 12, add=T, col=4, lty=3, lwd=3)
curve(log(logLikQNB(y=10, mu=x, phi=1.1, tau=.5)), 200, 2000)
curve(log(my_LBNB(y=10, mu=x, phi=1.1, tau=.5)), 200, 2000, add=T, col=4, lty=3, lwd=3)
abline(v=10)


curve(log(logLikQNB(y=0, mu=x, phi=1.1, tau=.5)), 0, 20)
curve(log(my_LBNB(y=0, mu=x, phi=1.1, tau=.5)), 0, 20, add=T, col=4, lty=3, lwd=3)
curve(log(logLikQNB(y=0, mu=x, phi=1.1, tau=.5)), 200, 2000)
curve(log(my_LBNB(y=0, mu=x, phi=1.1, tau=.5)), 200, 2000, add=T, col=4, lty=3, lwd=3)

curve(log(logLikQNB(y=1000, mu=x, phi=1.1, tau=.5)), 900, 1100)
curve(log(my_LBNB(y=1000, mu=x, phi=1.1, tau=.5)), 900, 1100, add=T, col=4, lty=3, lwd=3)
curve(log(logLikQNB(y=1000, mu=x, phi=1.1, tau=.5)), 0, 500)
curve(log(my_LBNB(y=1000, mu=x, phi=1.1, tau=.5)), 0, 500, add=T, col=4, lty=3, lwd=3)

save.image(img.name)
