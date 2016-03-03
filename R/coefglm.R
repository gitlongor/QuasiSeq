#' Coefficient and Bias From Generalized Linear Model Fit
#' 
#' \code{coef.glm} is a \code{S3} method computing the raw coefficients, 
#' leading order bias, and bias corrected coefficients from a 
#' generalized linear model (glm) fit. 
#' 
#' @param object An object with similar structure as returned from 
#'  \code{\link[stats:glm]{glm}}.
#' @param type A character scalar, being one of \code{'raw'}, 
#'  \code{'bias'}, and \code{'corrected'}. \code{'raw'} returns the 
#'  usual coefficients; \code{'bias'} returns the estimated leading order bias;
#'  and \code{'corrected'} returns the coefficients after subtracting
#'  the estimated bias.
#' @param ... Not used.
#'  
#' @return A numeric vector of requested components. 
#'  When \code{type='corrected'}, the \code{'bias'} attribute
#'  will be set to a numeric vector of estimated biases 
#'  being subtracted from the raw coefficients.
#' @author Long Qu <long.qu@wright.edu>
#' @references McCullagh and Nelder (1989) "Generalized Linear Models", 2nd edition. London: Chapman and Hall.
#' 
#' Cordeiro and McCullah (1991) "Bias Correction in Generalized Linear Models", \emph{Journal of the Royal Statistical Society: Series B}, \bold{53}, 629-643.
#'  
#' @examples 
#'  x=1:30
#'  y=rpois(30L, x/10)
#'  glmfit=glm(y~x, poisson('log'))
#'  coef(glmfit)
#'  coef(glmfit, type='bias')
#'  coef(glmfit, type='corrected')
#'  
#' @keywords methods
#' @concept bias correction
#' @export
coef.glm=function(object, type=c('raw','bias','corrected'), ...)
{
	type=match.arg(type)
	if(type=='raw') return(object$coefficients)
	
	good.wt=object$weights>0
	this.qr=object$qr
	this.hatd=.rowSums(qr.Q(this.qr)[,seq_len(this.qr$rank), drop=FALSE]^2, sum(good.wt), this.qr$rank)## not affected by pivoting

	this.mu=object$fitted.values[good.wt]
	this.eta=object$family$linkfun(this.mu)
	d2link = if(!is.null(object$family$`d2linkfun`)){
		object$family$d2linkfun(this.mu) ## this requires existence of d2linkfun components on the family
	}else switch(object$family$link, ## named links in stats::make.link
		'log'= -1/this.mu^2,
		'sqrt'=	-.25/this.mu^1.5,
		'logit'= (2*this.mu - 1)/(this.mu*(1-this.mu))^2,
		'probit'=local({.a=qnorm(this.mu); .a/dnorm(.a)^2}), 
		'cauchit'=-2*(base::pi)^2/tanpi(this.mu)/sinpi(this.mu)^2,
		'cloglog'=local({.a=1-this.mu; .b=log(.a); -(1+log(.a))/.a^2/.b^2}),
		'identity'=this.mu,
		'1/mu^2'=6/this.mu^4,
		'inverse'=2/this.mu^3,
		numDeriv::hessian(object$family$linkfun, this.mu)
	)
	this.w2ksi= 0.5* this.hatd / sqrt(object$weights[good.wt]) *d2link*object$family$mu.eta(this.eta)^2  
	this.bias=qr.coef(this.qr, this.w2ksi)
	if(type=='bias') return(this.bias)
	
	ans=object$coefficient-drop(this.bias)
	attr(ans, 'bias')=this.bias
	ans
}

