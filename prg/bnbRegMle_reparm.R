bnbRegMle_reparm=function(y, x, o, start=NULL, iter.max=1000L, eps.bet=1e-6, eps.tp=1e-3)
{
	n=NROW(x)
	p=NCOL(X)

	NLLmu = function(y,mu,phi,tau)
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
	betObj=function(bet, theta){
		NLLmu(y, exp(x%*%bet+o), theta[1L], theta[2L])
	}
	betGrad=function(bet, theta){
		phi=theta[1L]; tau=theta[2L]
		mu <- exp(x%*%bet+o)
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
		
		-colSums(as.vector(dmu) * mu.x)
	}
	betHess=function(bet, theta){
		phi=theta[1L]; tau=theta[2L]
		mu <- exp(x%*%bet+o)
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
	dtau.phi=function(theta, mu)
	{
		phi=theta[1L]; tau=theta[2L]
		const1=digamma((1-tau+2*tau*phi)/(tau*(phi-1)))
		const2=digamma((tau-phi-2*tau*phi)/(tau-tau*phi))
		const3=digamma(1/tau)

		dg1=digamma((1+mu*tau^2*phi+tau*(mu-1+2*phi))/(tau*(phi-1)))
		dg2=digamma((phi+mu*tau^2*phi+tau*(mu-1+y*(phi-1)+2*phi))/(tau*(phi-1)))
		dg3=digamma((mu*(1+tau*phi))/(phi-1))
		dg4=digamma(y+(mu+mu*tau*phi)/(phi-1))
		dg5=digamma(y+1/tau)

		dtau=sum(
			(1/(tau^2*(phi-1)))*(
				const1 -dg1
				+(phi-1)*(const3-dg5)
				+phi*(dg2-const2)
				+mu*tau^2*phi*(dg1-dg2-dg3+dg4)
			))

		dphi=sum(
			(1/(tau * (phi-1)^2))*(1 + tau) * (
				dg2 -const2 + const1 -dg1
				+ mu * tau * (dg3 -dg4 +dg2 -dg1)
			))
	 
		-c(dphi, dtau)
	}
	rhoxi2theta=function(rhoxi)rhoxi/(1-rhoxi)
	theta2rhoxi=function(theta)theta/(1+theta)
	drho.xi=function(rhoxi, mu)
	{
		theta=rhoxi2theta(rhoxi)
		dphiTau=dtau.phi(theta, mu)
		dphiTau/(1-rhoxi)^2
	}
	tpObj=function(theta, mu){
		NLLmu(y, mu, theta[1L], theta[2L])
	}
	rxObj=function(rhoxi, mu){
		theta=rhoxi2theta(rhoxi)
		NLLmu(y, mu, theta[1L], theta[2L])
	}
	
	optBet=function(start, theta){
		#tryCatch(
			nlminb(start, objective=betObj, gradient = betGrad, hessian = betHess,theta=theta)$par
		#, error=function(...)rep(NA_real_, length(start)))
	}
	optTheta=function(start, mu){
		#tryCatch(
			nlminb(start, objective=tpObj, gradient = dtau.phi, lower=1:0+1e-6,mu=mu)$par
		#, error=function(...)rep(NA_real_, length(start)))
	}
	optRhoxi=function(start, mu){
		#tryCatch(
			nlminb(start, objective=rxObj, gradient = drho.xi, lower=c(.5,0)+1e-6, upper=c(1, 1)-1e-6, mu=mu)$par
		#, error=function(...)rep(NA_real_, length(start)))
	}
	
	if(is.null(start)){
		.NotYetImplemented()
	}
	
	old.bet=start[1:p]
	old.tp =start[1:2+p]
	
	for(iter in seq_len(iter.max)){
		bet=optBet(old.bet, old.tp)
		mu=exp(x%*%bet+o)
		theta=rhoxi2theta(optRhoxi(theta2rhoxi(old.tp), mu))
		if(max(abs(bet-old.bet))<= eps.bet && 
		   max(abs(theta-old.tp))<=eps.tp) break
		old.bet=bet
		old.tp=theta
	}
	grad=c(betGrad(bet, theta), dtau.phi(theta, mu))
	
	if(iter==iter.max ) {
		warning('Maximum iterations reached')
	}
	list(beta=bet, tau=theta[2L], phi=theta[1L], 
		grad=grad, hess.bet=betHess(bet, theta),
		iter=iter)
}
