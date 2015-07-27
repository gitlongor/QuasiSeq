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

prg.name=paste('Coner', Sys.Date(), sep='_')
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

likQNB = Vectorize(likQNB, c('phi','y','tau'))



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

jpeg(file = "plot10.jpg")
curve(log(likQNB(y=10000, mu=x, phi=1.1, tau=.9)), 8000, 12000,main="y=10000,phi=1.1,tau=0.9",xlab="mean",ylab="Log likelihood")
curve(log(my_LBNB(y=10000, mu=x, phi=1.1, tau=.9)), 8000, 12000, add=T, col=4, lty=3, lwd=3,main="y=10000,phi=1.1,tau=0.9",xlab="mean",ylab="Log likelihood")
legend("topright", c("LogLQNB","LogLBNB"), cex=0.8,col=c("black","blue"),pch=21:22,lwd=2:3)
dev.off()

jpeg(file = "plot11.jpg")
curve(log(likQNB(y=10, mu=x, phi=1.1, tau=.5)), 200, 2000)
curve(log(my_LBNB(y=10, mu=x, phi=1.1, tau=.5)), 200, 2000, add=T, col=4, lty=3, lwd=3)
legend("topright", c("LogLQNB","LogLBNB"), cex=0.8,col=c("black","blue"),pch=21:22,lwd=2:3)
abline(v=10)
dev.off()

curve(log(likQNB(y=0, mu=x, phi=1.1, tau=.5)), 0, 20)
curve(log(my_LBNB(y=0, mu=x, phi=1.1, tau=.1)), 0, 20, add=T, col=4, lty=3, lwd=3)


curve(log(likQNB(y=0, mu=x, phi=1.1, tau=.5)), 200, 2000)
curve(log(my_LBNB(y=0, mu=x, phi=1.1, tau=.5)), 200, 2000, add=T, col=4, lty=3, lwd=3)

curve(log(likQNB(y=1000, mu=x, phi=1.1, tau=.5)), 900, 1100)
curve(log(my_LBNB(y=1000, mu=x, phi=1.1, tau=.5)), 900, 1100, add=T, col=4, lty=3, lwd=3)
curve(log(likQNB(y=1000, mu=x, phi=1.1, tau=.5)), 0, 500)
curve(log(my_LBNB(y=1000, mu=x, phi=1.1, tau=.5)), 0, 500, add=T, col=4, lty=3, lwd=3)



#Score functions comparison with mean
my_BNB_Score1<-function(y,mu,phi,tau){
	term1 = mu*(tau*phi+1)/(phi-1)
	term2 = (tau*(2*phi-1)+1)/(tau*(phi-1))
    term3= term1+term2+1/tau
	term4 = (tau*phi+1)/(phi-1)
	k=0:(y-1)
	terms1=NULL
	terms2=NULL
	for(i in 1:y){
	terms1[i]=(1/(k[i]+term1))
	terms2[i]=(1/(term3+k[i]))
}
    term5=sum(terms1)
    term6=sum(terms2)
	ans = term4*(term5+digamma(term1+term2)-
                   digamma(term3)-term6)
	y0=which(y==0)
	ans[y0]=term4*(digamma(term1+term2)-digamma(term3))[y0]
	ans	   
}
my_QNB_Score1<-function(y, mu,phi,tau){
   (y-mu)/(phi*(mu+tau*mu^2))
}
jpeg(file = "plot28.jpg")
curve(my_BNB_Score1(y=10, mu=x, phi=1.1, tau=0.5), 8, 20,col=4, lty=3, lwd=3,main="y=10,phi=1.1,tau=0.5",xlab="mean",ylab="Score")
curve(my_QNB_Score1(y=10, mu=x, phi=1.1, tau=0.5), 8, 20, add=T,main="y=10,phi=1.1,tau=0.5",xlab="mean",ylab="Score")
legend("topright", c("ScoreBNB","ScoreQNB"), cex=0.8,col=c("blue","black"),pch=21:22,lwd=3:4)
abline(h=0)
dev.off()

curve(my_BNB_Score1(y=0, mu=x, phi=1.1, tau=.5), 0, 1)
curve(my_QNB_Score1(y=0, mu=x, phi=1.1, tau=.5), 0, 1, add=T, col=4, lty=3, lwd=3)
abline(h=0)

#Hessian value
my_BNB_Hess<-function(y,mu,phi,tau){
	term1 = mu*(tau*phi+1)/(phi-1)
	term2 = (tau*(2*phi-1)+1)/(tau*(phi-1))
    term3= term1+term2+y+1/tau
	term4 = (tau*phi+1)/(phi-1)
	k=0:(y-1)
	terms1=NULL
	terms2=NULL
	for(i in 1:y){
	terms1[i]=(-1/(k[i]+term1)^2)
	terms2[i]=(-1/(term3+k[i])^2)
}
    term5=sum(terms1)
    term6=sum(terms2)
	ans=term4^2*(term5+trigamma(term1+term2)-
                   trigamma(term3)-term6)
	y0=which(y==0)
	ans[y0]=term4^2*(trigamma(term1+term2)-trigamma(term3))[y0]
	ans	
}

bnbScore=function(mu, y, phi, tau)-sum(my_BNB_Score1(y,mu,phi,tau))
qnbScore=function(mu, y, phi, tau)sum(my_QNB_Score1(y,mu,phi,tau))
#bnbHess=function(mu, y, phi, tau)structure(-sum(my_BNB_Hess(y,mu,phi,tau)), dim=c(1L,1L))
bnbHess=function(mu, y, phi, tau)-sum(my_BNB_Hess(y,mu,phi,tau))
#bnbNegLogLik=function(mu, y, phi, tau)min(.Machine$double.xmax, -sum(log(my_LBNB(y,mu,phi,tau))))
bnbNegLogLik=function(mu, y, phi, tau)-sum(log(my_LBNB(y,mu,phi,tau)))

#optim
bnbMle=function(y, phi, tau){
	#bnbNegLogLik.logmu=function(logmu)min(.Machine$double.xmax, -sum(log(my_LBNB(y,exp(logmu),phi,tau))))
	#bnbGrad.logmu=function(logmu)bnbScore(exp(logmu),y,phi,tau)/exp(logmu)
	
	par.hess=tryCatch(optim((max(1e-9,mean(y))), bnbNegLogLik,bnbScore,method='L-BFGS-B',lower=0, y=y, phi=phi, tau=tau,hessian=TRUE), error=function(...)NA_real_)
    par=par.hess$par
	hessian=par.hess$hessian
	variance=(hessian)
	var2=bnbHess(par, y=y, phi=phi, tau=tau)
	return(c(mu=par,var.numerical=variance,
			 var.theoretical=var2))
	}

qnbMle=function(y, phi, tau) mean(y)
	
getBnbQnb=function(reps,n, mu, tau, phi)#reps=number of replication
{replicate(reps,
	c(BNB=bnbMle(my.simu(n, mu, phi, tau), phi, tau),
      QNB=qnbMle(my.simu(n, mu, phi, tau), phi, tau)
	)
)
}	


#nlminb
bnbMle1=function(y, phi, tau)
nlminb(max(1e-9,mean(y)),bnbNegLogLik,bnbScore,scale=1,lower=0,y=y,phi=phi, tau=tau)

getBnbQnb1=function(reps,n, mu, tau, phi)
{replicate(reps,
	c(BNB=bnbMle1(my.simu(n, mu, phi, tau), phi, tau)$par,
	  QNB=qnbMle(my.simu(n, mu, phi, tau), phi, tau)
	)
)
}

	
#simulation
simulation=function(simu,n,mu,tau,phi)#simu=desired number of simulations  
{replicate(simu,my.simu(n, mu, phi, tau))}

	

##############################
phis=seq(1+1e-3, 3, length=20)
taus=seq(1e-3, 3, length=20)
mus =seq(.1, 1e4, length=20)
ns = c(2, 4, 8, 16, 32)
cases=expand.grid(n=ns, mu=mus, tau=taus, phi=phis)
reps=1e3L


getBnbQnb2=function(x) getBnbQnb(x[1],x[2],x[3],x[4])
est = apply(cases, 1, getBnbQnb2)
dim(est)=c(2L, reps, nrow(cases))



save.image(img.name)


 
