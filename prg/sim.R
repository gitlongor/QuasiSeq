
my.simu<-function(n, beta,phi,tau) {
mu <- exp(x%*%beta+o)
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
getBnbQnb=function(reps,n, beta, tau, phi)#reps=number of replication
{replicate(reps,
	Beta=Mlebeta(my.simu(n, beta, phi, tau), phi, tau)
      
)
}	
