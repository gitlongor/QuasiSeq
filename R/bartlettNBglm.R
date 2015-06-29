#bartlettNBglm=function(x, y, weights = rep(1, length(y)), offset = rep(0, length(y)), family, link='log')
#{
barlettFactor=function(glmFit)
{
  good.weights=glmFit$weights>0
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

  D1=diag(this.mu.eta)
  D2=diag(this.mu2.eta - (this.mu.eta)^2*(1+2*odisp*this.mu)/(this.mu+odisp*(this.mu)^2))
  W=D1%*%solve(this.var)%*%D1

  Q=xx%*%solve(t(xx)%*%W%*%xx)%*%t(xx)
  P=D1%*%Q%*%D1%*%solve(this.var)

  a=0
  for(i in 1:nrow(P))
  {
    a = a + P[i,i]^2*(-6/this.mu[i]+1/this.var[i]+6*this.var[i]/this.mu[i]^2)
  }

  b=0
  for(i in 1:nrow(P))
  {
    for(j in 1:ncol(P))
    {
      b = b + (P[i,i]*(-this.var[i]+2*this.var[i]^2/this.mu[i])/(this.mu[i]))*solve(this.var[i])*P[i,j]*(P[j,j]*(-this.var[j]+2*this.var[j]^2/this.mu[j])/(this.mu[j]))
    }
  }

  c=0
  for(i in 1:nrow(P))
  {
    for(j in 1:ncol(P))
    {
      c = c + (solve(this.var[i])*P[i,j])^3*(-this.var[i]+2*this.var[i]^2/this.mu[i])*(-this.var[j]+2*this.var[j]^2/this.mu[j])
    }
  }

  q=diag(Q)
  H=D2%*%solve(this.var)%*%(diag(nrow(P))-P)

  d = t(q)%*%H%*%D2%*%q

  e = sum(Q%*%Q%*%H%*%D2)

  f=0
  for(j in 1:ncol(W))
  {
    f = f + t(q)%*%H*q[j]*W[j,j]*(-this.var[j]+2*this.var[j]^2/this.mu[j])/(this.mu[j])
  }

  ep=-1/4*a+1/4*b+1/6*c-1/4*d+1/2*e-1/2*f
}
