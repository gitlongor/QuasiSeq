bnbRegMle=function(y, x, o, start=NULL, iter.max=1000L, eps.tp=1e-3)
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
	tpObj=function(theta, mu){
		NLLmu(y, mu, theta[1L], theta[2L])
	}
	optTheta=function(start, mu){
		#tryCatch(
			nlminb(start, objective=tpObj, gradient = dtau.phi, lower=1:0+1e-6,mu=mu)$par
		#, error=function(...)rep(NA_real_, length(start)))
	}
	if(is.null(start)){
		.NotYetImplemented()
	}
	old.tp =start[1:2+p]
	for(iter in seq_len(iter.max)){
		mu=exp(x%*%bet+o)
		theta=optTheta(old.tp, mu)
		if(max(abs(theta-old.tp))<=eps.tp) break
		old.tp=theta
	}
	grad= dtau.phi(theta, mu)
	if(iter==iter.max ) {
		warning('Maximum iterations reached')
	}
	list(tau=theta[2L], phi=theta[1L], 
		grad=grad,iter=iter)
}

