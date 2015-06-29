coef.glm=function(object, type=c('raw','bias','corrected'), ...)
{
	type=match.arg(type)
	if(type=='raw') return(object$coefficients)
	
	good.wt=object$weights>0
	this.qr=object$qr
	this.hatd=.rowSums(qr.Q(this.qr)[,seq_len(this.qr$rank), drop=FALSE]^2, sum(good.wt), this.qr$rank)## not affected by pivoting 
	## this.hatd / object$weights[good.wt] will be diag(Q matrix) in GLM book section 15.2

	this.mu=object$fitted.values[good.wt]
	this.eta=object$family$linkfun(this.mu)
	this.w2ksi= 0.5* this.hatd / sqrt(object$weights[good.wt]) *object$family$d2linkfun(this.mu)*object$family$mu.eta(this.eta)^2  ## this requires existence of d2linkfun components on the family
	## mu'' is the same as -d2linkfun(mu) * (du/deta)^3
	## thus mu''/mu' is -d2linkfun(mu) * (du/deta)^2
	
	this.w2ksi = this.w2ksi * object$prior.weights ## Cordeiro & McCullagh (JRSSB,1991)
	
	this.bias=qr.coef(this.qr, this.w2ksi)
	if(type=='bias') return(this.bias)
	
	ans=object$coefficient-drop(this.bias)
	attr(ans, 'bias')=this.bias
	ans
}

