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

library(QuasiSeq)

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
	y0=which(y==0)
	ans[y0]=exp(numerator-term3)[y0]
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

likQNB=function(y, mu, phi, tau)
{
	h=get.h(phi, y, tau)
	my_QLNB(y, mu, phi, tau) * h
}
likQNB = Vectorize(logLikQNB, c('phi','y','tau'))


#Simulation
my.simu<-function(n, mu,phi,tau) {
	alpha<-(tau*(2*phi-1)+1)/(tau*(phi-1))
	beta<-1/tau
	r=mu*(alpha-1)/beta
	p<-sort(rbeta(n,alpha,beta))
	y=rnbinom(n,r,prob=p)
	attr(y, 'phi')=phi
	attr(y, 'tau')=tau
	attr(y, 'mu')=mu
	y
}

my.simu(8, 10,2,.5)


curve(log(likQNB(y=10, mu=x, phi=1.1, tau=.5)), 8, 12)
curve(log(my_LBNB(y=10, mu=x, phi=1.1, tau=.5)), 8, 12, add=T, col=4, lty=3, lwd=3)
curve(log(likQNB(y=10, mu=x, phi=1.1, tau=.5)), 200, 2000)
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

#Score functions comparison with mean
my_BNB_Score1<-function(y,mu,phi,tau){
	term1 = mu*(tau*phi+1)/(phi-1)
	term2 = (tau*(2*phi-1)+1)/(tau*(phi-1))
    term3= term1+term2+y+1/tau
	term4 = (tau*phi+1)/(phi-1)
	ans = term4*(digamma(y+term1)+digamma(term1+term2)-
                   digamma(term3)-digamma(term1))
ans
}
my_QNB_Score1<-function(y, mu,phi,tau){
   (y-mu)/(phi*(mu+tau*mu^2))
}

curve(my_BNB_Score1(y=10, mu=x, phi=1.1, tau=.5), 9, 12)
curve(my_QNB_Score1(y=10, mu=x, phi=1.1, tau=.5), 9, 12, add=T, col=4, lty=3, lwd=3)
abline(h=0)

curve(my_BNB_Score1(y=0, mu=x, phi=1.1, tau=.5), 0, 1)
curve(my_QNB_Score1(y=0, mu=x, phi=1.1, tau=.5), 0, 1, add=T, col=4, lty=3, lwd=3)
abline(h=0)

bnbScore=function(mu, y, phi, tau)sum(my_BNB_Score1(y,mu,phi,tau))
qnbScore=function(mu, y, phi, tau)sum(my_QNB_Score1(y,mu,phi,tau))
bnbNegLogLik=function(mu, y, phi, tau)-sum(log(my_LBNB(y,mu,phi,tau)))


bnbMle=function(y, phi, tau)
	tryCatch(optim(max(1e-9,mean(y)), bnbNegLogLik, method='L-BFGS-B', lower=0, y=y, phi=phi, tau=tau)$par, error=function(...)NA_real_)
	

	
qnbMle=function(y, phi, tau) mean(y)	


phis=seq(1+1e-3, 3, length=20)
taus=seq(1e-3, 3, length=20)
mus =seq(.1, 1e4, length=20)
ns = c(2, 4, 8, 16, 32)
cases=expand.grid(n=ns, mu=mus, tau=taus, phi=phis)
reps=1e3L

getBnbQnb=function(n, mu, tau, phi)
{replicate(reps,
	c(BNB=bnbMle(my.simu(n, mu, phi, tau), phi, tau),
	  QNB=qnbMle(my.simu(n, mu, phi, tau), phi, tau)
	)
)
}
getBnbQnb2=function(x) getBnbQnb(x[1],x[2],x[3],x[4])
est = apply(cases, 1, getBnbQnb2)


save.image(img.name)


 