
my_LBNB_beta<-function(beta,y,x,o,phi,tau){
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
	sum(log(ans))
}

my_LBNB_mu<-function(y,mu,phi,tau)
{
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

dtau=function(y, mu, phi, tau)
{
-sum(
	(1/(tau^2 * (-1 + phi))) * (
		(-1 + phi) * digamma(1/tau) - 
		phi * digamma((tau - phi - 2 * tau * phi)/(tau - tau * phi)) - 
		mu * tau^2 * phi * digamma((mu * (1 + tau * phi))/(-1 + phi)) + 
		digamma(( 1 - tau + 2 * tau * phi)/(tau * (-1 + phi))) - 
		digamma(( 1 + mu * tau^2 * phi + tau * (-1 + mu + 2 * phi))/(tau * (-1 + phi))) + 
		mu * tau^2 * phi * digamma(( 1 + mu * tau^2 * phi + tau * (-1 + mu + 2 * phi))/(tau * (-1 + phi)))

		-(-1 + phi) * digamma(y + 1/tau) + 
		mu * tau^2 * phi * digamma(y + (mu + mu * tau * phi)/(-1 + phi)) + 
		phi * digamma(y + (phi + mu * tau^2 * phi + tau * (-1 + mu + 2 * phi))/(tau * (-1 + phi))) - 
		mu * tau^2 * phi * digamma(y + (phi + mu * tau^2 * phi + tau * (-1 + mu + 2 * phi))/(tau * (-1 + phi)))
	)
)
} 



My_phi_score<-function(y, mu,phi,tau)
{
    term1=mu*(tau*phi+1)/(phi-1)
	term2=(tau*(2*phi-1)+1)/(tau*(phi-1))
	term3= -mu*(tau+1)/(phi-1)^2
	term4=-(tau+1)/(tau*(phi-1)^2)
	
	digam.term1.term2=digamma(term1+term2)
	digam.term2.tau=digamma(term2+1/tau)
	digam.term2=digamma(term2)
	digam.term1.term2.tau=digamma(term1+term2+1/tau)
	
	sums1=sums2=numeric(length(y))
	for(i in seq_along(y)){
		k=seq_len(y[i])-1
		pos.denom=k+term1[i]
		neg.denom=k+term1[i]+term2+1/tau
		sums1[i] = sum(1/pos.denom)
		sums2[i]=sum(1/neg.denom)
	}

	dphi=sums1*term3+digam.term1.term2*(term3+term4)+digam.term2.tau*term4-
	      (digam.term1.term2.tau+sums2)*(term3+term4)-digam.term2*term4
	-sum(dphi)
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
	
	score
}

library(numDeriv)
grad(function(tau)my_LBNB_mu(tau=tau, y=10, phi=2, mu=10.5), .5)
dtau(y=10, mu=10.5, phi=2, tau=.5)

grad(function(phi)my_LBNB_mu(phi=phi, y=10, tau=.5, mu=10.5), 2)
My_phi_score(y=10, mu=10.5, phi=2, tau=.5)


grad(function(tau)my_LBNB_mu(tau=tau, y=10:20, phi=2, mu=10.5:20.5), .5)
dtau(y=10:20, mu=10.5:20.5, phi=2, tau=.5)

grad(function(phi)my_LBNB_mu(phi=phi, y=10:20, tau=.5, mu=10.5:20.5), 2)
My_phi_score(y=10:20, mu=10.5:20.5, phi=2, tau=.5)

## need: 
mle=function(y, x, o)
{
beta=function(y, x, o,phi=2,tau=.5){
par.beta=tryCatch(optim(max(0,10), my_LBNB_beta,My_score_Beta,method='L-BFGS-B',lower=-inf, y=y, phi=phi, tau=tau,hessian=TRUE), error=function(...)NA_real_)
beta=par.beta$par
beta
}
phi=function(y,x,o,tau=.5){
beta=runif(4, 0,1)
mu = exp(x%*%beta+o)
par.phi=tryCatch(optim(max(0,10), my_LBNB_beta,My_phi_score,method='L-BFGS-B',lower=-inf, y=y, mu=mu, tau=tau,hessian=TRUE), error=function(...)NA_real_)
phi=par.phi$par
phi
}
tau=function(y,x,o,phi=2){
beta=runif(4, 0,1)
mu = exp(x%*%beta+o)
par.tau=tryCatch(optim(max(0,10), my_LBNB_beta,dtau,method='L-BFGS-B',lower=-inf, y=y, phi=phi, mu=mu,hessian=TRUE), error=function(...)NA_real_)
tau=par.tau$par
tau
}
beta=beta(y, x, o,phi=2,tau=.5)
phi=phi(y,x,o,tau=.5)
tau=tau(y,x,o,phi=2)
return(list(beta=beta, phi=phi, tau=tau))
}