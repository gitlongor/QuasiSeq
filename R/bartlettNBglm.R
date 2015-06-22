bartlettNBglm=function(x, y, weights = rep(1, length(y)), offset = rep(0, length(y)), family, link='log')
{
  mu.eta=family$mu.eta
  variance=family$variance

  good.weights=weights>0
  xx=x[good.weights,,drop=FALSE]

  D1=diag(mu.eta)
  D2=diag()
  W=D1*inv(variance)*D1

  Q=xx*inv(t(xx)*W*xx)*t(xx)
  P=D1*Q*D1*inv(variance)

  a=0
  for(i in 1:nrow(P))
  {
    a = a+P[i,i]^2*mean(y[i]-mean(y[i]))^4/variance[i]-3
  }

  b=0
  for(i in 1:nrow(P))
  {
    for(j in 1:ncol(P))
    {
      b = b+(P[i,i]*mean(y[i]-mean(y[i]))*inv(variance)[i]*P[i,j]*(P[j,j]*mean(y[j]-mean(y[j]))))
    }
  }

  c=0
  for(i in 1:nrow(P))
  {
    for(j in 1:ncol(P))
    {
      c = c+(inv(variance)[i]*P[i,j])^3*mean(y[i]-mean(y[j]))^5
    }
  }

  q=as.matrix(1:nrow(Q))
  for(i in 1:nrow(Q))
  {
    q[1,i]=Q[i,i]
  }

  C=D2*inv(variance)*(diag(nrow(P))-P)*D2
  d=t(q)*C*q

  e=0
  for(i in 1:nrow(Q))
  {
    for(j in 1:ncol(Q))
    {
      e = e+Q[i,j]*Q[i,j]*C[i,j]
    }
  }

  f=0
  for(j in 1:ncol(W))
  {
    f = f + t(q)*D2*inv(variance)*(diag(nrow(P))-P)*q[j]*W[j]*mean(y[j]-mean(y[j]))
  }

  ep=-1/4*a+1/4*b+1/6*c-1/4*d+1/2*e-1/2*f
}
