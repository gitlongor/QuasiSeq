#bartlettNBglm=function(x, y, weights = rep(1, length(y)), offset = rep(0, length(y)), family, link='log')
#{
bartlettFactor=function(glmFit)
{
  good.weights=glmFit$weights>0  ## ??
  yy=glmFit$y[good.weights]; xx=glmFit$x[good.weights,,drop=FALSE]; oo=glmFit$offset[good.weights]

  family=glmFit$family

  mu.eta=family$mu.eta
  variance=family$variance
  linkinv=family$linkinv
  d2linkfun=family$d2linkfun

  odisp=1/family$getTheta()
  this.beta=glmFit$coef
  this.eta=as.vector(xx%*%this.beta+oo)
  this.mu=linkinv(this.eta)
  this.mu.eta=mu.eta(this.eta)
  this.var=variance(this.mu)
  this.d2g=d2linkfun(this.mu)
  this.mu2.eta=-this.d2g*(this.mu.eta)^3

  V=diag(this.var); VI=solve(V)
  D1=diag(this.mu.eta); d1vec=this.mu.eta
  D2=diag(this.mu2.eta - (this.mu.eta)^2*(1+2*odisp*this.mu)/(this.mu+odisp*(this.mu)^2)); 
	d2vec = this.mu2.eta - (this.mu.eta)^2*(1+2*odisp*this.mu)/(this.mu+odisp*(this.mu)^2)
  W=D1%*%VI%*%D1; wvec=d1vec^2/this.var

  Qmat=xx%*%solve(t(xx)%*%W%*%xx)%*%t(xx) ## switch to QR decomp
  P=D1%*%Qmat%*%D1%*%VI

  rho4vec= - 6/this.mu + 1/this.var+6*this.var/this.mu^2
  a = sum(diag(P)^2 * rho4vec)
  if(FALSE){
	a=0
	for(i in 1:nrow(P))
	{
	a = a + P[i,i]^2*(-6/this.mu[i]+1/this.var[i]+6*this.var[i]/this.mu[i]^2)
	}
  }

  k3=(-this.var+2*this.var^2/this.mu)
  if(FALSE){
	  b=0
	  for(i in 1:nrow(P))
	  {
		for(j in 1:ncol(P))
		{
		  b = b + (P[i,i]*k3[i]/this.var[i])/(this.var[i])*P[i,j]*(P[j,j]*k3[j]/this.var[j])
		}
	  }
  }else b = drop((diag(P)*k3/this.var^2)%*%P%*%(diag(P)*k3/this.var))

  if(FALSE){
	  cval=0
	  for(i in 1:nrow(P))
	  {
		for(j in 1:ncol(P))
		{
		  cval = cval + (solve(this.var[i])*P[i,j])^3*k3[i]*k3[j]
		}
	  }
  }else  cval=sum((k3%o%k3)*P^3 / this.var^3)


  qvec=diag(Qmat)
  if(FALSE){
		H=D2%*%VI%*%(diag(nrow(P))-P)
  }else H=(d2vec/this.var) * (diag(nrow(P))-P)

  if(FALSE){
		d = t(qvec)%*%H%*%D2%*%qvec
  }else d = drop(qvec%*%H%*%(d2vec*qvec))
  
  if(FALSE){
		e = sum(Qmat*Qmat*(H%*%D2))
  }else e = sum(Qmat*Qmat*(H*D2[col(H)]))

  qstar = qvec * wvec * k3 / this.var
  if(FALSE){
	  
  }else f = drop(qvec%*%H%*%qstar)
  
  ans.MN = -1/4*a+1/4*b+1/6*cval-1/4*d+1/2*e-1/2*f  ## epsilon_p
  
  
  ans.C = 0
  
  list(nelder=ans.MN, cordeiro = ans.C)
}
