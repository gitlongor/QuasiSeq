
my_LBNB<-function(beta,y,x,o,phi,tau){
mu <- exp(x%*%beta+o)
	term1 = mu*(tau*phi+1)/(phi-1)
	term2 = (tau*(2*phi-1)+1)/(tau*(phi-1))
	numerator = lbeta(term1+term2,y+1/tau)
	term3 = lbeta(term2,1/tau)
	ans = exp(		
		numerator-lbeta(y, term1)-term3
	)/y
	y0=which(y==0)
	ans[y0]=exp(numerator-term3)[y0]
	-sum(log(ans))
}


My_score_Beta<-function(beta, y, x, o, tau, phi)
{
	n=NROW(x)

	mu <- exp(x%*%beta+o)
	mu.x=as.vector(mu)*x
	

	mu.tau.phip1.dphin1 = mu*(tau*phi+1)/(phi-1)
	t.2phin1p1.dtphin1=(tau*(2*phi-1)+1)/(tau*(phi-1))
	t.2phin1p1.dtphin1.p1tau=t.2phin1p1.dtphin1 +1/tau
	tau.phip1.dphin1=(tau*phi+1)/(phi-1)

	pos.psi0.arg=mu.tau.phip1.dphin1+t.2phin1p1.dtphin1
	pos.psi0=digamma(pos.psi0.arg)
	neg.psi0=digamma(pos.psi0.arg+1/tau)
	
	sums=numeric(n)
	for(i in seq_len(n)){
		k=seq_len(y[i])-1L
		pos.denom=k+mu.tau.phip1.dphin1[i]
		neg.denom=pos.denom+ t.2phin1p1.dtphin1.p1tau
		sums[i] = sum(1/pos.denom-1/neg.denom)
	}
	dmu=(sums+pos.psi0-neg.psi0)*tau.phip1.dphin1
	
	score=colSums(as.vector(dmu) * mu.x)
	
	-score
}


 
My_hess_beta<-function(beta, y, x, o, tau, phi)
{
	n=NROW(x)
	mu <- exp(x%*%beta+o)
	mu.x<-as.vector(mu)*x
	
	
	mu.tau.phip1.dphin1 = mu*(tau*phi+1)/(phi-1)
	t.2phin1p1.dtphin1=(tau*(2*phi-1)+1)/(tau*(phi-1))
	t.2phin1p1.dtphin1.p1tau=t.2phin1p1.dtphin1 +1/tau
	tau.phip1.dphin1=(tau*phi+1)/(phi-1)
	
	pos.psi0.arg=mu.tau.phip1.dphin1+t.2phin1p1.dtphin1
	pos.psi0=digamma(pos.psi0.arg)
	neg.psi0=digamma(pos.psi0.arg+1/tau)
	pos.psi2=trigamma(pos.psi0.arg)
	neg.psi2=trigamma(pos.psi0.arg+1/tau)
	
	sums1=numeric(n)
	sums2=numeric(n)
	
	for(i in seq_len(n)){
		k=seq_len(y[i])-1L
		pos.denom1=k+mu.tau.phip1.dphin1[i]
		neg.denom1=pos.denom1+ t.2phin1p1.dtphin1.p1tau
		sums1[i] = sum(1/pos.denom1-1/neg.denom1)
		neg.denom2=(k+mu.tau.phip1.dphin1[i])
		squ.neg.denom2=(k+mu.tau.phip1.dphin1[i])^2
		pos.denom2=(neg.denom2+ t.2phin1p1.dtphin1.p1tau)^2
		sums2[i] = sum(1/pos.denom2-1/squ.neg.denom2)
	}
	dmu=(sums1+pos.psi0-neg.psi0)*tau.phip1.dphin1
	dbetaj.dmuilogl.ai<-as.vector((sums2+pos.psi2-neg.psi2)*(tau.phip1.dphin1)^2*mu)
	dmu.xi.mu<-t(as.vector(dmu)*x)%*%(as.vector(mu)*x)
	Hess<-dmu.xi.mu+t(dbetaj.dmuilogl.ai*x)%*%(as.vector(mu)*x)
	-Hess
}



Mlebeta=function(y, x,o,phi,tau,start=rep(1,NCOL(x))){
	
	par.hess=#tryCatch(
			#optim(start, my_LBNB,My_score_Beta,method='BFGS',x=x,o=o, y=y, phi=phi, tau=tau,hessian=TRUE)
			nlminb(start, objective=my_LBNB, gradient = My_score_Beta, hessian = My_hess_beta,x=x,o=o, y=y, phi=phi, tau=tau)	
		#, error=function(...)NA_real_)

	par=par.hess$par
	#hessian=par.hess$hessian
	#variance=solve(hessian)
	var2=solve(My_hess_beta(par, y=y, x=x,o=o, phi=phi, tau=tau))
	return(list(bet=par, #var.numerical=variance,
			 var.theoretical=var2, optim=par.hess))
}


if(FALSE){
	library(numDeriv)
	set.seed(234234)
	x=cbind(1, matrix(runif(10*3, 0, 1), nc=3))
	bet=runif(4, 0,1)
	xb=x%*%bet
	y=rnbinom(10, 1/(tau=2), mu=drop(exp(xb)))

	grad(my_LBNB, x=bet,,,, y=y, x,o=rep(0,10), phi=1+1e-3, tau=2)
	My_score_Beta(bet, y=y, x=x, o=rep(0,10), tau=2, phi=1+1e-3)

	hessian(function(bet)my_LBNB(beta=bet,x=x,o=rep(0,10), y=y, phi=1+1e-3, tau=2), bet) 
	My_hess_beta(beta=bet,y=y,x=x,o=rep(0,10), phi=1+1e-3, tau=2)

	
	(tmp=Mlebeta(y, x=x,o=rep(0,length(y)), phi=1+1e-3, tau=2))

	My_score_Beta(tmp$bet, y=y, x=x, o=rep(0,10), tau=2, phi=1+1e-3)
}
