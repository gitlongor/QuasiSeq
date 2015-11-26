bnbRegMle_CommonTau=function(y, x, o, beta.start=NULL, phitau.start=NULL, iter.max=1000L, eps.mu=1e-6, eps.tp=1e-3, verbose=FALSE)
{
	G=NROW(y)
	n=NROW(x)
	p=NCOL(X)

	NLLmu = function(y,mu,phi,tau)
	{
		term1 = mu*(tau*phi+1)/(phi-1)
		term2 = (tau*(2*phi-1)+1)/(tau*(phi-1))
		numerator = lbeta(term1+term2,y+1/tau)
		term3 = lbeta(term2,1/tau)
		ans = 		
			numerator-lbeta(y, term1)-term3 - log(y)
		
		y0=which(y==0)
		ans[y0]=(numerator-term3)[y0]
		-sum(ans)
	}
	betObj=function(bet, theta,y){
		NLLmu(y, exp(x%*%bet+o), theta[1L], theta[2L])
	}
	betGrad=function(bet, theta,y){
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
	betHess=function(bet, theta,y){
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
	
	dphi=function(phi, tau, mu, y)
	{
		const1=digamma((1-tau+2*tau*phi)/(tau*(phi-1)))
		const2=digamma((tau-phi-2*tau*phi)/(tau-tau*phi))
		const3=digamma(1/tau)

		dg1=digamma((1+mu*tau^2*phi+tau*(mu-1+2*phi))/(tau*(phi-1)))
		dg2=digamma((phi+mu*tau^2*phi+tau*(mu-1+y*(phi-1)+2*phi))/(tau*(phi-1)))
		dg3=digamma((mu*(1+tau*phi))/(phi-1))
		dg4=digamma(y+(mu+mu*tau*phi)/(phi-1))
		dg5=digamma(y+1/tau)

		dphi=sum(
			(1/(tau * (phi-1)^2))*(1 + tau) * (
				dg2 -const2 + const1 -dg1
				+ mu * tau * (dg3 -dg4 +dg2 -dg1)
			))
	 
		-dphi
	}
	dtau=function(tau, phi, mu, y)
	{
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
 
		-dtau
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

		dphi=rowSums(
			(1/(tau * (phi-1)^2))*(1 + tau) * (
				dg2 -const2 + const1 -dg1
				+ mu * tau * (dg3 -dg4 +dg2 -dg1)
			))
	 
		-c(dphi, dtau)
	}
	tpObj=function(phi, tau, y, mu){
		NLLmu(y, mu, phi, tau)
	}
	
	optBet=function(start, phi, tau){
		ans=start
		for(g in seq_len(G)) {
			ans[g,]=
			nlminb(start[g,], objective=betObj, gradient = betGrad, hessian = betHess, theta=c(phi[g],tau), y=y[g,])$par
		}#, error=function(...)rep(NA_real_, length(start)))
		ans
	}
	optTheta=function(start, mu, iter.max=15L){
		obj=function(theta, mu)tpObj(theta[-(G+1L)], theta[G+1L], y, mu)
		nlminb(start, objective=obj, gradient = dtau.phi, lower=rep(1:0+1e-6,c(G,1L)), upper=rep(1e6,G+1L), mu=mu, control=list(iter.max=iter.max))$par
	}	
	optTheta=function(start, mu, iter.max=15L){
		old.ans=ans=start
		iter.inner=1L
		repeat{
			for(g in seq_len(G)){
				ans[g]=nlminb(old.ans[g], objective=tpObj, gradient = dphi, lower=1+1e-6, upper=1e6, tau=old.ans[G+1], y=y[g,], mu=mu[g,])$par
			}
			ans[G+1L]=nlminb(old.ans[G+1L], objective=tpObj, gradient = dtau, lower=1e-6, upper=1e6, phi=ans[-(G+1L)], y=y, mu=mu)$par
			
			if(	iter.inner>=iter.max+iter || 
				max(abs(old.ans-ans))<=eps.tp
			){
				print(iter.inner)
				return(ans)
			} 
			ans=old.ans
			iter.inner=iter.inner+1L
		}
	}
	
	if(is.null(start)){
		.NotYetImplemented()
	}
	
	old.bet=beta.start
	old.tp =phitau.start
	
	oo=o[rep(seq_len(n), each=G)]
	old.mu=exp(tcrossprod(old.bet,x)+oo)
	old.ll=-NLLmu(y, old.mu, old.tp[-(G+1L)], old.tp[G+1L])
	
	start.taus=10^seq(-6,6, length=50)
	start.phiTaus=as.matrix(expand.grid( start.taus))
	start.lls=numeric(nrow(start.phiTaus))

	if(verbose )
		cat("iter=", 0, "\tphi=", quantile(old.tp[-(G+1L)]), "\ttau=", old.tp[G+1L], "\n")
	for(iter in seq_len(iter.max)){

		bet=optBet(old.bet, old.tp[-(G+1L)], old.tp[G+1L])
		mu=exp(tcrossprod(bet,x)+oo)
		
		for(i in seq_along(start.lls)){
			start.lls[i]=-NLLmu(y, mu, old.tp[-(G+1L)], start.phiTaus[i,1L])
		}
		theta=optTheta(c(old.tp[-(G+1L)], start.phiTaus[which.max(start.lls),1]), mu)
		ll=-NLLmu(y, mu, theta[-(G+1L)], theta[G+1L])

		if((#max(abs(bet-old.bet))<= eps.bet ||
			max(abs(mu-old.mu)/pmax(1,abs(mu)))<= eps.mu     ) &&
			max(abs(ll-old.ll))<= eps.mu      &&
		   max(abs(theta-old.tp))<=eps.tp) break
		   
		if(verbose && iter%%verbose==0L){
			cat("iter=", iter, "\tphi=", quantile(abs(old.tp[-(G+1L)]-theta[-(G+1L)])), "\ttau=", old.tp[G+1L],  "\tdiffMu=", quantile(apply(abs(mu-old.mu),1L,max)), "\timaxDiffMu", which.max(apply(abs(mu-old.mu)/pmax(1,abs(mu)),1L,max)), "\tLL=", ll, "\n")
		}		
		old.bet=bet
		old.tp=theta
		old.mu=mu
		old.ll=ll
	}
	grads=matrix(NA_real_, G, p)
	for(g in seq_len(G)) grads[g,]=betGrad(bet[g,], c(theta[g],theta[G+1L]), y[g,])
	grad.phi.tau=dtau.phi(theta, mu)
	
	if(iter==iter.max ) {
		warning('Maximum iterations reached')
	}
	list(beta=bet, tau=theta[G+1L], phi=theta[-(G+1L)], 
		beta.grad=grads, phitau.grad=grad.phi.tau,
		# hess.bet=betHess(bet, theta),
		mu=mu, logLik=ll, 
		iter=iter)
}


if(FALSE){

	debug(bnbRegMle_CommonTau)

	set.seed(234345L)
	tmpidx=seq_len(NROW(counts.nb))
	#tmpidx=sort(sample(nrow(counts.nb),  500L) )
	
	tmp1=bnbRegMle_CommonTau(y=counts.nb[tmpidx,], x=X, o=offs, beta.start=fit$coeff[tmpidx,], phitau.start=c(pmin(1e6,tauPhihats[tmpidx,1]),median(tauPhihats[tmpidx,2])), iter.max=1000L, eps.mu=1e-6, eps.tp=1e-3, verbose=1L)

	tmp5=bnbRegMle_CommonTau(y=counts.nb[tmpidx,], x=X, o=offs, beta.start=fit$coeff[tmpidx,], phitau.start=c(1+1e-5, 1e-5), iter.max=1000L, eps.mu=1e-6, eps.tp=1e-3, verbose=1L)

	NLLmu = function(y,mu,phi,tau)
	{
		term1 = mu*(tau*phi+1)/(phi-1)
		term2 = (tau*(2*phi-1)+1)/(tau*(phi-1))
		numerator = lbeta(term1+term2,y+1/tau)
		term3 = lbeta(term2,1/tau)
		ans = 		
			numerator-lbeta(y, term1)-term3 - log(y)
		
		y0=which(y==0)
		ans[y0]=(numerator-term3)[y0]
		-sum(ans)
	}
	start.phis=1+10^(-5:5)
	start.taus=10^(-5:5)
	start.phiTaus=expand.grid(start.phis, start.taus)
	start.mus=fit$fitted.values+1e-3
	start.lls=numeric(nrow(start.phiTaus))
	for(i in seq_along(start.lls)){
		start.lls[i]=-NLLmu(counts.nb[tmpidx,], start.mus[tmpidx,], start.phiTaus[i,1L], start.phiTaus[i,2L])
	}
	
}


	tmp1=bnbRegMle_CommonTau(y=counts.nb[tmpidx,], x=X, o=offs, beta.start=fit$coeff[tmpidx,], phitau.start=c(pmin(1e6,tauPhihats[tmpidx,1]),median(tauPhihats[tmpidx,2])), iter.max=1000L, eps.mu=1e-3, eps.tp=1e-3, verbose=1L)
