
my_LBNB<-function(y, beta,x,o,phi,tau){
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
	ans
}


My_score_Beta<-function(beta, y, x, o, tau, phi)
{
	n=NROW(x)

	mu <- exp(x%*%beta+o)
	mu.x=as.vector(mu)*x
	

	mu.tau.phip1.dphin1 = mu*(tau*phi+1)/(phi-1)
	t.2phin1p1.dtphin1=(tau*(2*phi-1)+1)/tau/(phi-1)
	t.2phin1p1.dtphin1.p1tau=t.2phin1p1.dtphin1 +1/tau
	tau.phip1.dphin1=(tau*phi+1)/(phi-1)

	pos.psi0.arg=mu.tau.phip1.dphin1+t.2phin1p1.dtphin1
	pos.psi0=digamma(pos.psi0.arg)
	neg.psi0=digamma(pos.psi0.arg+1/tau)
	
	sums=numeric(n)
	for(i in seq_len(n)){
		k=seq_len(y)-1L
		pos.denom=k+mu.tau.phip1.dphin1[i]
		neg.denom=pos.denom+ t.2phin1p1.dtphin1.p1tau
		sums[i] = sum(1/pos.denom-1/neg.denom)
	}
	dmu=(sums+pos.psi0-neg.psi0)*tau.phip1.dphin1
	
	score=colSums(as.vector(dmu) * mu.x)
	
	score
}
library(numDeriv)
grad(function(beta)log(my_LBNB(beta=beta, y=10,x=matrix(1:3,nrow=1),o=1, phi=2, tau=.5)), matrix(1:3,nrow=3))
My_score_Beta(matrix(1:3,nrow=3), y=10, x=matrix(1:3,nrow=1), o=1, tau=.5, phi=2)

grad(function(beta)log(my_LBNB(beta=beta, y=0,x=matrix(1:3,nrow=1),o=1, phi=2, tau=.5)), matrix(1:3,nrow=3))
My_score_Beta(matrix(1:3,nrow=3), y=0, x=matrix(1:3,nrow=1), o=1, tau=.5, phi=2)

 
My_diag_hess<-function(beta, y, x, o, tau, phi)
{
	n=NROW(x)
	mu <- exp(x%*%beta+o)
	mu.x<-as.vector(mu)*x
	dmu2<-as.vector(mu)*(x*x)
	#dmu2.dbi.dbj<-as.vector(mu)*combn(x)
	
	mu.tau.phip1.dphin1 = mu*(tau*phi+1)/(phi-1)
	t.2phin1p1.dtphin1=(tau*(2*phi-1)+1)/tau/(phi-1)
	t.2phin1p1.dtphin1.p1tau=t.2phin1p1.dtphin1 +1/tau
	tau.phip1.dphin1=(tau*phi+1)/(phi-1)
	
	pos.psi0.arg=mu.tau.phip1.dphin1+t.2phin1p1.dtphin1
	pos.psi0=trigamma(pos.psi0.arg)
	neg.psi0=trigamma(pos.psi0.arg+1/tau)
	
	sums1=numeric(n)
	sums2=numeric(n)
	
	for(i in seq_len(n)){
		k=seq_len(y)-1L
		pos.denom1=k+mu.tau.phip1.dphin1[i]
		neg.denom1=pos.denom1+ t.2phin1p1.dtphin1.p1tau
		sums1[i] = sum(1/pos.denom1-1/neg.denom1)
		neg.denom2=(k+mu.tau.phip1.dphin1[i])^2
		pos.denom2=(neg.denom2+ t.2phin1p1.dtphin1.p1tau)^2
		sums2[i] = sum(1/pos.denom2-1/neg.denom2)
	}
	dmu=(sums1+pos.psi0-neg.psi0)*tau.phip1.dphin1
	dbetaj.dmuilogl<-as.vector((sums2+pos.psi0-neg.psi0)*(tau.phip1.dphin1)^2)*mu.x
	
	diago<-colSums(as.vector(dmu)*dmu2+mu.x*dbetaj.dmuilogl)
	#offdiag<-
	diago
}

hessian(function(beta)log(my_LBNB(beta=beta, y=10,x=matrix(1:3,nrow=1),o=1, phi=2, tau=.5)), matrix(1:3,nrow=3))
grad(function(beta)My_score_Beta(beta=beta, y=10, x=matrix(1:3,nrow=1), o=1, tau=.5, phi=2),matrix(1:3,nrow=3))
My_diag_hess(matrix(1:3,nrow=3), y=10, x=matrix(1:3,nrow=1), o=1, tau=.5, phi=2)

