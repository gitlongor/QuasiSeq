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
  this.dvar=family$dvar(this.mu)
  this.d2var=family$d2var(this.mu)
  this.d2g=d2linkfun(this.mu)
  this.mu2.eta=-this.d2g*(this.mu.eta)^3
  this.cum3 = family$cumulant3(mu=this.mu, var=this.var)
  this.cum4 = family$cumulant4(mu=this.mu, var=this.var)

  if(FALSE){
	V=diag(this.var); VI=solve(V)
	D1=diag(this.mu.eta);
	D2=diag(this.mu2.eta - (this.mu.eta)^2*(1+2*odisp*this.mu)/(this.mu+odisp*(this.mu)^2));
	W=D1%*%VI%*%D1;
	Qmat=xx%*%solve(t(xx)%*%W%*%xx)%*%t(xx)
	P=D1%*%Qmat%*%D1%*%VI
  }else {
	d1vec=this.mu.eta
	d2vec = this.mu2.eta - (this.mu.eta)^2 / this.var * this.dvar ##
	wvec=d1vec^2/this.var; wvec.half=sqrt(wvec)
	QQ = (tcrossprod(qr.Q(glmFit$qr)))
	Qmat = 1/wvec.half * QQ / wvec.half[col(QQ)]
	P = d1vec * Qmat * (d1vec/this.var)[col(Qmat)]
  }

  rho4vec= this.cum4 / this.var^2    ## rho4vec = - 6/this.mu + 1/this.var+6*this.var/this.mu^2
  if(FALSE){
	a=0
	for(i in 1:nrow(P))
	{
	a = a + P[i,i]^2*(-6/this.mu[i]+1/this.var[i]+6*this.var[i]/this.mu[i]^2)
	}
  }else a = sum(diag(P)^2 * rho4vec)

  k3= this.cum3		## k3 = (-this.var+2*this.var^2/this.mu)
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
  }else e = sum(Qmat*Qmat*(H*d2vec[col(H)]))

  qstar = qvec * wvec * k3 / this.var
  if(FALSE){

  }else f = drop(qvec%*%H%*%qstar)

  ans.MN = -1/4*a+1/4*b+1/6*cval-1/4*d+1/2*e-1/2*f  ## epsilon_p

  wvec.cord=(this.mu.eta)^2/this.var; wvec.half.cord=sqrt(wvec.cord); wmat.cord=diag(wvec.cord)
  ZZ.cord=(tcrossprod(qr.Q(glmFit$qr)))
  phi.cord=diag(nrow(wmat.cord))
  Z.cord=1/wvec.half.cord*ZZ.cord*phi.cord/wvec.half.cord[col(ZZ.cord)]
  Z3.cord=Z.cord^3

  H.cord=1/this.var*this.mu2.eta*(this.mu2.eta-4*w.cord*this.dvar)+w.cord^2*(2/this.var*this.dvar^2-this.d2var)
  F.cord=1/this.var*this.mu.eta*this.mu2.eta
  G.cord=1/this.var*this.mu.eta*this.mu2.eta-1/this.var^2*this.dvar*this.mu.eta^3
  Zd.cord=diag(Z.cord)

  a.cord=sum(diag(phi.cord*H.cord*Zd.cord^2))
  b.cord=sum(phi.cord*G.cord*Z3.cord*(F.cord+G.cord)*phi.cord)
  c.cord=drop(sum(phi.cord*F.cord*(2*Z3.cord+3*Zd.cord %*% Z.cord %*% Zd.cord)*F.cord*phi.cord))
  ans.C = 1/4*a.cord-1/3*b.cord+1/12*c.cord ## epsilon_k

  list(nelder=ans.MN, cordeiro = ans.C)
}
