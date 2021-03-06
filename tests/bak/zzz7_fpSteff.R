	fbrNBglm.fit=function(x, y, weights = rep(1, length(y)), offset = rep(0, length(y)), family, link='log', odisp, control = fbrNBglm.control()) 
{
	if(missing(family)) family=negbin(link, odisp)
	if(family$link!='log') stop('Currently only log link has been implemented.')
	
	nobs=NROW(y)
	good.weights=weights>0
	ngood=length(good.weights)
	if(ngood==0L) stop('None of the weights is positive')
	ww=weights[good.weights]
	if(any(ww!=ww[1L])) stop('Bias reduction in the presence of non-equal positive prior weights is currently not available')
	ww[]=1

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
			this.qr=qr(this.w2x,  tol=qr.tol)
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
		exact = (infoParmsk==0 || infoParmsj==1) && constOffset
		#attr(exact, 'group')=unclass(group)
		attr(exact, 'group')=group
		attr(exact, 'oneWayX')=oneWayX
		attr(exact, 'constOffset')=constOffset
		exact
	}
	
	fullFactorial1Step=function(groupX) ## yy, odisp
	{	ns=.colSums(groupX, ngood, rk)
		fitted.mean=ybars=crossprod(groupX, yy)/ns  
		off=crossprod(groupX, oo)/ns  
		eval(expr.1step)   # defined below
		fitted.coef=xxuniqInv %*% ( linkfun(fitted.mean)-off)
	}
	expr.1step=
	if(infoParmsk==0){  ## expected information
		expression(
			fitted.mean <- (ns*ybars+0.5)/(ns-odisp*.5)
		)
	}else if(infoParmsj==1){ # linear solution (including observed information)
		expression({
			.tmp=2*ns + infoParmsk * odisp^infoParmsm
			fitted.mean=(.tmp*ybars+1)/(.tmp-odisp)
		})
	}else if(infoParsmj==0){ # quadratic solution
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
			start=fullFactorial1Step(oneWayX)
			attr(start, 'method')='exact'
			attr(start, 'success')=TRUE
			attr(start, 'iter')=1L
		}else{
			doIteration=TRUE
			if(is.null(start)) eval(getMuStart)
						
			ss=exp(oo)
			ssOdisp=ss*odisp
			if(infoParmsj==1 && infoParmsk==1 && infoParmsm==1){ ## OI full factorial iteratation equation
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
				rhs1=function(this.mu, grpId){
					onePlusOdispMuScale=1+ssOdisp*this.mu #[oneWayGroup]
					tmpMat=workMat/onePlusOdispMuScale
					tmpMat[,3L]=tmpMat[,3L]/onePlusOdispMuScale
					ssSyssY=crossprod(oneWayX[,grpId], tmpMat)
					ans=ssSyssY[, 2L]/ssSyssY[, 1L]+.5*ssSyssY[, 3L]/ssSyssY[,1L]^2
					if(any(is.na(ans))) ans=exp(log(ssSyssY[, 2L])-log(ssSyssY[, 1L]))+.5*exp(log(ssSyssY[, 3L])-2*log(ssSyssY[,1L]))
					ans
				}
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
				rhs1=function(this.mu, grpId){
					onePlusOdispMuScale=1+ssOdisp*this.mu #[oneWayGroup]
					tmpMat=workMat/onePlusOdispMuScale
					tmpMat[,3L:4L]=tmpMat[,3L:4L]/onePlusOdispMuScale^infoParmsj
					ssSyssY=crossprod(oneWayX[,grpId], tmpMat)
					( ssSyssY[, 1L]*ssSyssY[, 2L] + .5*(ssSyssY[, 4L] +ssSyssY[, 1L]  )  )/
						( ssSyssY[, 1L]^2 + .5*ssSyssY[, 3L]  )
				}
			}
			
			tryCatch({ ## fixed-point with Steffensen algorithm
				it=1L
				iterMax=control$maxit
				startExpXb0=startExpXb1=as.vector(crossprod(oneWayX, exp(xx%*%start))/oneWayN)

				## try better starting values using the cheaper rhs-lhs as surrogate adjusted score
				startAdjscore = rhs(startExpXb0) - startExpXb0
					yy.bak=yy; oo.bak=oo
					yy=yy/exp(oo-mean(oo)); oo[]=mean(oo)
				fff=fullFactorial1Step(oneWayX)
					yy=yy.bak; oo=oo.bak
				fffmu=as.vector(crossprod(oneWayX, exp(xx%*%fff))/oneWayN)
				fffAdjscore=rhs(fffmu) - fffmu
				if(!all(is.finite(startAdjscore)) || sum(startAdjscore^2)>sum(fffAdjscore^2))	{
					start=fff
					startExpXb0=startExpXb1=fffmu
					startAdjscore = fffAdjscore
				}
				
				if(FALSE){
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
				}
				## Fixed Point With Steffensen Iterations
				nextExpXb = startExpXb1
				for(grpId in seq_len(rk)){ 
					it=1L
					mu0=mu=startExpXb1[grpId]
					steff=FALSE
					while(it <= control$maxit){
						nextmu=rhs1(mu, grpId)
						if (abs(mu - nextmu) <= control$tol ) break
						if(steff && (steffDenom<-nextmu-2*mu+mu0) !=0) {
							bak=nextmu
							nextmu=mu0-(mu-mu0)^2/steffDenom  #  Steffensen's method (Aitken acceleration)
							if(any(!is.finite(nextmu))) {
								nextmu=bak
							}
						}
						steff = !steff
						it=it+1L
						mu0=mu
						mu=nextmu
					}
					nextExpXb[grpId]=nextmu
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
	

	if(doIteration){ print("doIteration")
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


fbrNBglm.control=
function (standardizeX = TRUE, coefOnly=TRUE, solvers=nlSolvers, verbose = FALSE, maxit=25L, start = NULL, infoParms=list(j=1,k=1,m=1), tol=sqrt(.Machine$double.eps)) 
{
	stopifnot(all(sort(names(infoParms))==c('j','k','m')))
    structure(list(standardizeX = standardizeX, coefOnly = coefOnly, infoParms=infoParms, solvers = solvers, verbose = verbose, maxit=maxit, start = start, tol = tol), class = "fbrNBglm.control")
}

