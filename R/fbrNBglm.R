#' Firth-Type Bias-Reduced Negative Binomial Log-Linear Model
#' 
#' This is the estimation algorithm for generalized linear model with negative binomial responses
#' and log link function. Sore functions are modified properly such that 
#' bias in the coefficients are reduced. 
#' 
#' @param x,y,weights,offset Defined the same as in \code{\link{glm.fit}}.
#' @param family The same as in \code{\link{glm.fit}}, but the default 
#'      (and the only currently supported) choice is \code{negbin{'log', odisp}}, 
#'      where \code{odisp} is defined below.
#' @param odisp A numeric scalar of negative binomial over-dispersion parameter. 
#'      This is the same as \code{1/size}, where \code{size} is the parameter 
#'      as in \code{\link{dnbinom}}.
#' @param control A list returned from \code{fbrNBglm.control}.
#' @return It depends.
#' @export
#' @name fbrNBglm
#' @rdname fbrNBglm
#' @author Long Qu <long.qu@wright.edu>
#' @keywords models regression iteration
#' @concept bias reduction
#' @concept Firth
#' @examples 
#'  ## prepare example data
#'  data(mockRNASeqData)
#'  x=mockRNASeqData$design.matrix
#'  y=mockRNASeqData$counts[3462,]
#'  offset = log(mockRNASeqData$estimated.normalization)
#'  overDisp = mockRNASeqData$estimated.nbdisp[3462]
#'  nbfam = negbin('log', overDisp)
#'  
#'  ## usual maximum likelihood estimate
#'  coef(glm.fit(x, y, offset=offset, family=nbfam))
#'  
#'  ## 2nd-order biased reduced fit with observed information
#'  ctrl= fbrNBglm.control(infoParms=list(j=1,k=1,m=1), order=2L, coefOnly=TRUE)
#'  fbrNBglm.fit(x, y, offset=offset, family=nbfam, control=ctrl)
#'
#'  ## 2nd-order biased reduced fit with expected information
#'  ctrl= fbrNBglm.control(infoParms=list(j=0,k=1,m=1), order=2L, coefOnly=TRUE)
#'  fbrNBglm.fit(x, y, offset=offset, family=nbfam, control=ctrl)
#'  
#'  ## 3rd-order biased reduced fit with observed information
#'  ## Not available yet if offsets are non-constants with a treatment
#'  offset.avg = ave(offset, mockRNASeqData$treatment)
#'  ctrl= fbrNBglm.control(infoParms=list(j=1,k=1,m=1), order=3L, coefOnly=TRUE)
#'  fbrNBglm.fit(x, y, offset=offset.avg, family=nbfam, control=ctrl)
#'
fbrNBglm.fit=function(x, y, weights = rep(1, length(y)), 
                      offset = rep(0, length(y)), family, 
                      odisp, control = fbrNBglm.control()) 
{
	if(missing(family)) family=negbin('log', odisp)
	if(family$link!='log') stop('Currently only log link has been implemented.')
	
	nobs=NROW(y)
	good.weights=weights>0
	ngood=length(good.weights)
	if(ngood==0L) stop('None of the weights is positive')
	ww=weights[good.weights]
	if(any(ww!=ww[1L])) {
		warning('Bias reduction in the presence of non-equal positive prior weights has not been thoroughly tested. Use this option with caution.')
		# ww[]=1
	}

	if(!is.matrix(x)) x=as.matrix(x)
	ncolx=NCOL(x)
	yy=y[good.weights]; xx=x[good.weights,,drop=FALSE]; oo=offset[good.weights]
	start=control$start
	if(isTRUE(control$standardizeX)){
		x.norm=.colSums(xx*xx, ngood, ncolx) # sqrt will be taken later
		#x.stdCols=apply(xx,2L,sd)>0   ## apply is slow
		# x.stdCols=sqrt(diag(var(xx)))>0 ## still slow
		x.stdCols=x.norm > .colSums(xx, ngood, ncolx)^2 / ngood
		x.norm=sqrt(x.norm)
		if(any(x.stdCols)){
			xx[,x.stdCols]=sweep(xx[,x.stdCols,drop=FALSE],2L,x.norm[x.stdCols],'/')
			if(!is.null(start)) {
				start[x.stdCols] = start[x.stdCols] * x.norm[x.stdCols]
				control$start = start
			}
		}
	}

	odisp=1/family$getTheta()	
    variance <- family$variance
    linkfun <- family$linkfun
    linkinv <- family$linkinv
    mu.eta <- family$mu.eta	
	d2link = family$d2link
	dvar = family$dvar
	

	xxqr=qr(xx); rk=xxqr$rank
	if(rk<ncolx){
		#xx=qr.Q(xxqr)[,seq_len(rk)]%*%qr.R(xxqr)[seq_len(rk),seq_len(rk)]  ## working x mat always full rank
		xx=xx[,xxqr$pivot[seq_len(rk)],drop=FALSE]  ## this should avoid round-off errors
		if(!is.null(control$start)) {
			start=start[xxqr$pivot[seq_len(rk)]]
			control$start = start
		}
	}
	xxuniq=unique(xx)
	
	infoParmsj=control$infoParms$j
	infoParmsk=control$infoParms$k
	infoParmsm=control$infoParms$m
	fbrOrder =control$order
	
	.C(C_initQRdecomp, ngood, rk)
	on.exit(.C(C_finalQRdecomp))  ## should not call adjScoreFunc anymore
	adjScoreFunc=function(bet, approxJacob=FALSE) ## xx, oo, yy in evalFrame; j,k,m
	{## G =O.I.;  R=d{X'W(Y-mu)}
		this.eta=as.vector(xx%*%bet+oo)
		this.mu=linkinv(this.eta)
		this.mu.eta=mu.eta(this.eta)

		this.var=variance(this.mu)
		#this.d2g=d2link(this.mu)
		#this.dvar=dvar(this.mu)
		#this.dg=1/this.mu.eta
		
		this.weight=this.mu.eta^2/this.var
		this.resid=(yy-this.mu)
		this.wresid=this.resid*this.weight
		score=crossprod(xx, this.wresid/this.mu)
		
		this.w2x=sqrt(this.weight)*xx
		if(approxJacob) return(-crossprod(this.w2x))
		
		this.bias=try(.Call(C_getGlmBias, rtwx = this.w2x, wrt = sqrt(this.weight), ngood, rk))
		if(inherits(this.bias, 'try-error')) { ## original R version that is known working. 
			## the next three rows consume 60% of the time
			this.qr=qr(this.w2x,  tol=control$qr.tol)
			this.hatd=.rowSums(qr.Q(this.qr)[,seq_len(this.qr$rank), drop=FALSE]^2, ngood, this.qr$rank)## not affected by pivoting
			this.bias=qr.coef(this.qr, -0.5*this.hatd/sqrt(this.weight))
			this.bias[is.na(this.bias)]=0
		}
		
		this.adjWt=this.weight*(
			this.resid*infoParmsk*(this.var*d2link(this.mu)+dvar(this.mu)/this.mu.eta)^infoParmsm/(this.var/this.mu.eta)^infoParmsj +1
			## k=0    : this is EI weight
			## k=j=m=1:  this is OI weight
		)
		this.adjInfo=crossprod(xx, this.adjWt*xx) 
		
		adjScore=-this.adjInfo%*%this.bias
		as.vector(score + adjScore)
	}
	approxJacob=NULL
	attr(adjScoreFunc, 'getApproxJacob')=function(...)approxJacob
	
	test.1stepFF=function()
	{
		if(rk!=NROW(xxuniq)) stop('this function should only be used for one-way designs')
		#group=grpDuplicated(xx, factor=TRUE)
		#constOffset = sapply(split(oo, group), function(ooo)max(ooo)==min(ooo)) # faster than by, tapply, ave
		
		group=grpDuplicated(xx)
		constOffset = all(ave(oo, group) == oo)
		oneWayX=1 * (group==rep(seq_len(rk), each=ngood)); dim(oneWayX)=c(ngood, rk)
	
		#exact = (infoParmsk==0 || infoParmsj==1) && all(constOffset)
		exact =  constOffset && (
			(fbrOrder==2L && (infoParmsk==0 || infoParmsj==1) ) || 
			(fbrOrder==3L && infoParmsj==1) 
			)
		#attr(exact, 'group')=unclass(group)
		attr(exact, 'group')=group
		attr(exact, 'oneWayX')=oneWayX
		attr(exact, 'constOffset')=constOffset
		exact
	}
	
	fullFactorial1Step=function( ## yy, odisp, ngood, rk, oo
		groupX, 
		ns=.colSums(groupX, ngood, rk), 
		ybars=crossprod(groupX, yy)/ns,
		off=crossprod(groupX, oo)/ns, 
		fitted.mean    ## only for shuting up CRAN checker
		) 
	{	#ns=.colSums(groupX, ngood, rk)
		#fitted.mean=ybars=crossprod(groupX, yy)/ns  
		#off=crossprod(groupX, oo)/ns  
		eval(expr.1step)   # defined below
		fitted.coef=xxuniqInv %*% ( linkfun(fitted.mean)-off)
	}
	expr.1step=
	if(infoParmsk==0){  ## expected information
		expression(
			fitted.mean <- (ns*ybars+0.5)/(ns-odisp*.5)
		)
	}else if(infoParmsj==1 && fbrOrder ==2L){ # linear solution (including observed information)
		expression({
			.tmp=2*ns + infoParmsk * odisp^infoParmsm
			fitted.mean=(.tmp*ybars+1)/(.tmp-odisp)
		})
	}else if (infoParmsj==1 && fbrOrder ==3L){ # 3rd order: quadratic solution 
		expression({
			.tmpktaum=infoParmsk*odisp^infoParmsm
			.tmpc=1+.tmpktaum*ybars;
			.tmpb=12*ns+24*ns^2*ybars+5*.tmpktaum+12*ns*ybars*.tmpktaum+6*ybars*.tmpktaum^2-ybars*.tmpktaum*odisp;
			.tmpa=-24*ns^2+12*ns*odisp-odisp^2-12*ns*.tmpktaum-6*.tmpktaum^2+7*.tmpktaum*odisp;
			fitted.mean = (-.tmpb - sqrt(.tmpb^2 - 4 * .tmpa * .tmpc))/2/.tmpa
		})
	}else if(infoParmsj==0){ # quadratic solution
		expression({
			.tmpktaum=infoParmsk*odisp^infoParmsm
			.tmp=.5*(ybars+(odisp-2*ns)/.tmpktaum/odisp - 1/odisp)
			fitted.mean = .tmp + sqrt(.tmp^2 +(1+(2*ns+.tmpktaum)*ybars)/.tmpktaum/odisp)
		})
	}else { ##CAUTION: this is not exact; simply a default (Jeffreys) to allow later iterations
		expression({
			fitted.mean <- ybars + .5/ns
		})
	}
	getMuStart=expression(
	{
		eval(family$initialize)
		etastart=linkfun(mustart[good.weights])
		start=lm.fit(xx, etastart, offset=oo)$coef
	})
	if(rk<NROW(xxuniq)){
		doIteration=TRUE
		if(is.null(start)){
			eval(getMuStart, envir=sys.frame(sys.nframe()))
			startAdjscore=adjScoreFunc(start,approxJacob=FALSE)
			if(FALSE){
				## one-step starting values : not very effective
				xxkm=kmeans(qr.Q(xxqr), rk)		## slow
				approx1wayx=model.matrix(~0+as.factor(xxkm$cluster))  ## slow
				xxuniqInv=diag(1, rk, rk)
					yy.bak=yy; oo.bak=oo
					yy=yy/exp(oo-mean(oo)); oo[]=mean(oo)
				fff=fullFactorial1Step(approx1wayx)
					yy=yy.bak; oo=oo.bak
				fff=qr.coef(xxqr, fff[xxkm$cluster])
				fffAdjscore=adjScoreFunc(fff, approxJacob=FALSE)
				if(sum(fffAdjscore^2)<sum(startAdjscore^2)){
					start=fff
				}
			}
		}
		if(!all(is.finite(start)) # || !all(is.finite(startAdjscore)) 
		){
			eval(getMuStart, envir=sys.frame(sys.nframe()))
			# startAdjscore=adjScoreFunc(start,approxJacob=FALSE)
		}
	}else 
	if(rk==NROW(xxuniq)){
		FFtestRslt=test.1stepFF()
		oneWayGroup=attr(FFtestRslt, 'group')
		oneWayX=attr(FFtestRslt, 'oneWayX')
		oneWayN=.colSums(oneWayX, ngood, rk)
		xxuniqInv=if(rk==1L) matrix(1/xxuniq) else if(rk==2L) solve22(xxuniq) else solve(xxuniq)
		
		if(FFtestRslt){
			doIteration=FALSE
			start=fullFactorial1Step(oneWayX, ns=oneWayN)
			attr(start, 'method')='exact'
			attr(start, 'success')=TRUE
			attr(start, 'iter')=0L
		}else{
			doIteration=TRUE
			if(is.null(start)) eval(getMuStart)
						
			ss=exp(oo)
			ssOdisp=ss*odisp
			if(infoParmsj==1 && infoParmsk==1 && infoParmsm==1 && fbrOrder==2L){ ## OI full factorial iteratation equation
				workMat=cbind(ss, yy, ss*(1+yy*odisp))
				rhs=function(this.mu){
					onePlusOdispMuScale=1+ssOdisp*this.mu[oneWayGroup]
					tmpMat=workMat/onePlusOdispMuScale
					tmpMat[,3L]=tmpMat[,3L]/onePlusOdispMuScale
					ssSyssY=crossprod(oneWayX, tmpMat)
					ans=ssSyssY[, 2L]/ssSyssY[, 1L]+.5*ssSyssY[, 3L]/ssSyssY[,1L]^2
					if(any(is.na(ans))) ans=exp(log(ssSyssY[, 2L])-log(ssSyssY[, 1L]))+.5*exp(log(ssSyssY[, 3L])-2*log(ssSyssY[,1L]))
					ans
				}
			}else if(infoParmsj==1 && fbrOrder==3L){ ## 3rd order correction
				.NotYetImplemented()
			}else{	## full factorial iteration (general) equation
				tmp=ss*infoParmsk*odisp^infoParmsm
				workMat=cbind(ss, yy, ss2=ss*tmp, sy=yy*tmp)
				rhs=function(this.mu){
					onePlusOdispMuScale=1+ssOdisp*this.mu[oneWayGroup]
					tmpMat=workMat/onePlusOdispMuScale
					tmpMat[,3L:4L]=tmpMat[,3L:4L]/onePlusOdispMuScale^infoParmsj
					ssSyssY=crossprod(oneWayX, tmpMat)
					( ssSyssY[, 1L]*ssSyssY[, 2L] + .5*(ssSyssY[, 4L] +ssSyssY[, 1L]  )  )/
						( ssSyssY[, 1L]^2 + .5*ssSyssY[, 3L]  )
				}
			}
			
			tryCatch({ ## fixed-point algorithm with bisection bracketing and Steffensen acceleration
				it=1L
				iterMax=control$maxit
				startExpXb0=startExpXb1=as.vector(crossprod(oneWayX, exp(xx%*%start))/oneWayN)

				## try better starting values using the cheaper rhs-lhs as surrogate adjusted score
				startAdjscore = rhs(startExpXb0) - startExpXb0
					yy.bak=yy; oo.bak=oo
					yy=yy/exp(oo-mean(oo)); oo[]=mean(oo)
				fff=fullFactorial1Step(oneWayX, ns=oneWayN)
					yy=yy.bak; oo=oo.bak
				fffmu=as.vector(crossprod(oneWayX, exp(xx%*%fff))/oneWayN)
				fffAdjscore=rhs(fffmu) - fffmu
				if(!all(is.finite(startAdjscore)) || sum(startAdjscore^2)>sum(fffAdjscore^2))	{
					start=fff
					startExpXb0=startExpXb1=fffmu
					startAdjscore = fffAdjscore
				}
				
				bracketFactor = rep(0.618, rk)
				expXbLower =expXbUpper = startExpXb0
				idx0 = startAdjscore < 0
				if(any( idx0 ) ){
					bf=bracketFactor; bf[!idx0]=1 
					repeat{
						expXbLower = expXbLower * bf
						tmpIdx = rhs(expXbLower)[idx0] > expXbLower[idx0]
						if(all(tmpIdx)) { break
						}else bf[idx0 & tmpIdx] = 1
					}
				}
				idx0=!idx0
				if ( any(idx0) ){
					bf=bracketFactor; bf[!idx0]=1
					repeat{
						expXbUpper= expXbUpper / bf
						tmpIdx = rhs(expXbUpper)[idx0] < expXbUpper[idx0]
						if(all(tmpIdx)) { break
						}else bf[idx0 & tmpIdx] = 1
					}
				}
				
				rhsNextExpXb=rhs(startExpXb1)
				while(it<=iterMax){ ## bisection-protected Steffensen-accelerated fixed-point iterations
					nextExpXb=rhsNextExpXb
					if(it%%2L==0L && all((steffDenom<-nextExpXb-2*startExpXb1+startExpXb0)!=0)) {
						bak=nextExpXb
						nextExpXb=startExpXb0-(startExpXb1-startExpXb0)^2/steffDenom  #  Steffensen's method (Aitken acceleration)
						if(any(!is.finite(nextExpXb))) {
							nextExpXb=bak
						}
					}
					
					outsideIdx =(nextExpXb < expXbLower | nextExpXb > expXbUpper) 
					if(any(outsideIdx)) { ## fixed point steps outside bracketing interval
						nextExpXb[outsideIdx]=(expXbLower + .5*(expXbUpper-expXbLower))[outsideIdx]
						rhsNextExpXb = rhs(nextExpXb)
							tmpIdx = rhsNextExpXb[outsideIdx] > nextExpXb[outsideIdx]
							if(any(tmpIdx)) 	expXbLower[outsideIdx][tmpIdx] =  nextExpXb[outsideIdx][tmpIdx]
							tmpIdx = !tmpIdx
							if(any(tmpIdx)) 	expXbUpper[outsideIdx][tmpIdx] =  nextExpXb[outsideIdx][tmpIdx]
					}else rhsNextExpXb = rhs(nextExpXb)
					
					insideIdx = !outsideIdx
					if(any(insideIdx)) { ## fixed point stays inside bracketing interval
						tmpIdx = rhsNextExpXb[insideIdx] > nextExpXb[insideIdx]
						if(any(tmpIdx)) expXbLower[insideIdx][tmpIdx] = nextExpXb[insideIdx][tmpIdx]
						tmpIdx = !tmpIdx
						if(any(tmpIdx)) expXbUpper[insideIdx][tmpIdx] = nextExpXb[insideIdx][tmpIdx]
					}
					
					if( sqrt(sum((nextExpXb-startExpXb1)^2))   <= control$tol 
					 || sqrt(sum((expXbLower - expXbUpper)^2)) <= control$tol
					) {
						diffbet=abs((nextExpXb-startExpXb1)/nextExpXb)
						if(all(abs(xxuniqInv%*%(diffbet-.5*diffbet^2+diffbet^3/3)) <=control$tol) )
							break  ## 3-term Taylor approx log(startExpXb)-log(nextExpXb)
					}

					it=it+1L
					startExpXb0=startExpXb1
					startExpXb1=nextExpXb
				}
				start= as.vector(xxuniqInv %*% linkfun(nextExpXb))
				startAdjscore=adjScoreFunc(start,approxJacob=FALSE)
				if(sqrt(sum(startAdjscore^2))<=control$tol) {
					doIteration=FALSE
					attr(start, 'method')='fullFactorialIter'
					attr(start, 'success')=TRUE
					attr(start, 'iter')=it
				}
			}, error=function(e){doIteration <<- TRUE; print(e); NULL})  ## tryCatch
		}
	}else stop("rank of x is larger than number of unique rows of x")
	

	if(doIteration){
		approxJacob=adjScoreFunc(start,approxJacob=TRUE)
		ans=suppressWarnings(nlsolve(start, adjScoreFunc, control$solvers, control))
		if(!attr(ans, 'nlsolve.success')){
			if(!is.null(control$start) && any(start!=control$start)) {
				start=control$start  ## restart using user specified starting value
				ans=nlsolve(start, adjScoreFunc, control$solvers, control)
			}else warning('None of the non-linear equation solvers succeeded.')
		}
		attrs=attributes(ans)
		names(attrs)=gsub('^nlsolve\\.', '', names(attrs))
		attributes(ans)=attrs
	}else ans=start

	if(!isTRUE(control$coefOnly)) finalAdjScore = adjScoreFunc(ans)
	
	#.C(C_finalQRdecomp)  ## should not call adjScoreFunc anymore
	
	if(rk<ncolx){
		ans0=ans
		ans=rep(NA_real_, ncolx)
		ans[xxqr$pivot[seq_len(rk)]]=ans0
		attributes(ans)=attributes(ans0)
	}
	
	if(isTRUE(control$standardizeX)){
		if(any(x.stdCols))
			ans[x.stdCols]=ans[x.stdCols]/x.norm[x.stdCols]
		xx=x[good.weights,,drop=FALSE]
	}

	if(control$coefOnly) {
		ans
	}else{
		ans0=ans
		ans0[is.na(ans0)]=0
		this.ctrl=control
		this.ctrl$maxit=0L
		ans=withCallingHandlers(glm.fit3(x=x, y=y, weights = weights, start = ans0, offset = offset, family = family,  control = this.ctrl, intercept = TRUE), simpleWarning=ignorableWarnings) 
		ans$converged = attr(ans0, 'success')
		ans$method = attr(ans0, 'method')
		ans$iter = attr(ans0, 'iter')
		ans$adjusted.score = finalAdjScore
		
		ans
	}
}

#' @rdname fbrNBglm
#' @param standardizeX A logical scalar. If \code{TRUE}, columns of 
#'      design matrix will be standardized with norm 1 during the 
#'      fitting process, except for columns with identical values. 
#' @param coefOnly A logical scalar. If \code{TRUE}, only the regression
#'      coefficients will be returned. This is useful when being called
#'      from other functions, e.g. \code{\link{NBDev}}. If \code{coef=FALSE},
#'      a \code{glm} object will be returned. 
#' @param solvers The non-linear equation solvers to be used if iterative 
#'      fitting is necessary.
#' @param verbose A logical scalar, indicating whether intermediate messages 
#'      should be printed.
#' @param maxit A positive integer, the maximum number of iterations allowed
#'      if iterative fitting is necessary. 
#' @param start A numeric vector of starting values, with length being the
#'      same as the number of columns of \code{x}.
#' @param infoParms A list of three components, named \code{j}, \code{k}, 
#'      and \code{m}, respectively. These parameters control the type of 
#'      adjustment made to the score function. When all j, k, and m are 1, 
#'      the observed information is used in the adjustment. When \code{k=0},
#'      the expected information is used. See reference for details. 
#' @param order A positive integer. Usually this should be set to 2, 
#'      indicating the second order, i.e., \eqn{O(n^{-1})}{O(n^{-1})} bias 
#'      being reduced by the adjustment. For one-way design, if \code{infoParms$j=1} and offsets are constants within each treatment, the third
#'      order bias reduction with \code{order=3} is also support, resulting 
#'      in an estimate with both \eqn{O(n^{-1})}{O(n^{-1})} order and 
#'      \eqn{O(n^{-2})}{O(n^{-2})} order biases being reduced. 
#' @param tol Small positive integer, indicating the desired accuracy
#'      of parameter estimates.
#' @param qr.tol The same as the \code{tol} argument of \link[base:qr]{qr}.
#' @export
fbrNBglm.control=
function (standardizeX = TRUE, coefOnly=TRUE, solvers=nlSolvers, 
          verbose = FALSE, maxit=25L, start = NULL, 
          infoParms=list(j=1,k=1,m=1), order=2L, 
          tol=sqrt(.Machine$double.eps), qr.tol=tol) 
{
	stopifnot(all(sort(names(infoParms))==c('j','k','m')))
	if (order == 3L){
		if(infoParms$j!=1) stop(message("Third order bias reduction with 'infoParms$j != 1' is not implemented yet"))
	}else if (order != 2L){
		stop(message('Currently "order" needs to be 2 or 3.'))
	}
    structure(list(standardizeX = standardizeX, coefOnly = coefOnly, infoParms=infoParms, order=order, solvers = solvers, verbose = verbose, maxit=maxit, start = start, tol = tol, qr.tol=qr.tol), class = "fbrNBglm.control")
}

