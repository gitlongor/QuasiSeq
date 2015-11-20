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

dtau.phi=function(theta, y, mu)#####
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
 
	-c(dtau, dphi)
}


if(FALSE){
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
#recrrence relationship for tau?
	library(numDeriv)
	grad(function(phi)my_LBNB_mu(phi=phi, y=10:20, tau=.5, mu=10.5:20.5), 2)
	My_phi_score(y=10:20, mu=10.5:20.5, phi=2, tau=.5)

	grad(function(tau)my_LBNB_mu(tau=tau, y=10:20, phi=2, mu=10.5:20.5), .5)
	grad(function(phi)my_LBNB_mu(phi=phi, y=10:20, tau=.5, mu=10.5:20.5), 2)
	dtau.phi(y=10:20, mu=10.5:20.5, c(phi=2, tau=.5))
}
