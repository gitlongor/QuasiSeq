ignorableWarnings=function(w)
{
	harmless=c(
		"step size truncated due to divergence",
		"step size truncated: out of bounds", 
		"step size truncated due to increasing deviance", 
		"inner loop 3: cannot correct step size",
		"inner loop 3; cannot correct step size",
		"glm.fit: algorithm did not converge",
		"glm.fit: algorithm stopped at boundary value",
		"glm.fit3: algorithm did not converge. Try increasing the maximum iterations",
		"glm.fit3: algorithm stopped at boundary value", 
		"glm.fit3: fitted rates numerically 0 occurred"
	)
	if (!(w$message %in% harmless)) return(w)
	
	invokeRestart("muffleWarning")
}

#' @rdname glmSolve
#' @export
#' @format \code{glmSolvers} is currently a list of length 4. See 
#' \code{solvers} argument for details and the example. 
#' @examples
#' str(glmSolvers) ## available solvers stored in glmSolvers object
glmSolvers=list(
	glm=list(solve=expression(withCallingHandlers(glm(formula=formula, family=family, data=data, control=control, ...),simpleWarning=ignorableWarnings	)),
			 success=expression(ans$converged), 
			 result=expression(ans)),
			 
	glm.fit3=list(solve=expression(withCallingHandlers(glm(formula=formula, family=family, data=data, control=control, method=glm.fit3, ...),simpleWarning=ignorableWarnings)),
			 success=expression(ans$converged), 
			 result=expression(ans)), 
	nlminb=list(solve=expression({
		this.ctrl=control; 
		this.ctrl$maxit=1L; this.initGlm=withCallingHandlers(glm(formula=formula, family=family, data=data, control=this.ctrl, ...), simpleWarning=ignorableWarnings);
		this.wt=this.initGlm$weights; 
		this.idxn=which(this.wt>0); this.idxp=which(!is.na(this.initGlm$coefficients)); 
		this.off=if(is.null(this.initGlm$offset)) 0 else this.initGlm$offset; 
		this.x=model.matrix(this.initGlm);  		
		dev.func=function(bet){
			sum(family$dev.resids(this.initGlm$y[this.idxn], family$linkinv(this.x[this.idxn,this.idxp,drop=FALSE]%*%bet[this.idxp]+ this.off[this.idxn]), this.wt[this.idxn]))
		}
		nlminb(this.initGlm$coefficients, dev.func, control=list(iter.max=control$maxit, eval.max=round(control$maxit/3L*4L), abs.tol=control$epsilon^2, rel.tol=control$epsilon, x.tol=control$epsilon, trace=FALSE))
		}), 
		success=expression(ans$convergence==0), 
		result=expression({
			this.rslt=ans$par; this.rslt[is.na(this.rslt)]=0; 
			this.ans=withCallingHandlers(tryCatch(glm(formula=formula, family=family, data=data, control=this.ctrl, start=this.rslt, method=glm.fit3, ...), simpleError=function(message, call = NULL)message), simpleWarning=ignorableWarnings)
			if(!inherits(this.ans, 'glm')) next
			this.ans
		})), 
	BFGS=list(solve=expression({
		this.ctrl=control; 
		this.ctrl$maxit=1L; this.initGlm=withCallingHandlers(glm(formula=formula, family=family, data=data, control=this.ctrl, ...), simpleWarning=ignorableWarnings);
		this.wt=this.initGlm$weights; 
		this.idxn=which(this.wt>0); this.idxp=which(!is.na(this.initGlm$coefficients)); 
		this.off=if(is.null(this.initGlm$offset)) 0 else this.initGlm$offset; 
		this.x=model.matrix(this.initGlm);  		
		dev.func=function(bet){
			sum(family$dev.resids(this.initGlm$y[this.idxn], family$linkinv(this.x[this.idxn,this.idxp,drop=FALSE]%*%bet[this.idxp]+ this.off[this.idxn]), this.wt[this.idxn]))
		}
		this.start=this.initGlm$coefficients; if(length(this.idxp)<length(this.start)) this.start[-this.idxp]=0
		optim(this.start, dev.func, method='BFGS', control=list(maxit=control$maxit, reltol=control$epsilon, trace=FALSE))
		}), 
		success=expression(ans$convergence==0), 
		result=expression({
			this.rslt=ans$par; this.rslt[is.na(this.rslt)]=0 ;
			this.ans=withCallingHandlers(tryCatch(glm(formula=formula, family=family, data=data, control=this.ctrl, start=this.rslt, method=glm.fit3, ...), simpleError=function(message, call = NULL)message), simpleWarning=ignorableWarnings)
			if(!inherits(this.ans, 'glm'))	next
			this.ans
		}))
)

#' Generalized linear model (glm) solvers
#' 
#' \code{glmsolve} sequencially tries multiple glm solvers and returns the 
#' result from the success. This tries to avoid numerical instabilities
#' that occassionally affect some glm solving routines. If one fails, the 
#' next solver will be tried. 
#' 
#' Currently, the supported are solvers are \code{c('glm', 'glm.fit3', 
#' 'nlminb', 'BFGS')}. \code{'glm'} is the one that comes with the default
#' \code{stats} package. \code{'glm.fit3'} is a modification with better 
#' stability (but slightly slower). \code{'nlminb'} uses the general 
#' optimization routine \code{\link[stats]{nlminb}}. \code{'BFGS'} uses 
#' the \code{method='BFGS'} option provided by \code{\link{optim}}.
#' 
#' @param formula,family,data Identical to those of \code{\link{glm}}.
#' @param control A list of named control options that specifies the 
#' details for each solver. Named elements not used by any solver will 
#' be ignored. 
#' @param solvers A named list of glm solvers, default to \code{glmSolvers}. 
#' Each element is itself 
#' a list with components named \code{solve}, \code{success}, and 
#' \code{result}, each of which is an expression. \code{solve} is an 
#' expression used to fit the glm model, the result of which is assigned
#' to an object \code{ans}. \code{success} is an expression that returns a 
#' logical scalar, indicating whether the solver has succeeded. It may 
#' refer to the \code{ans} returned by \code{solve}. Often, it will check 
#' the \code{converged} element of \code{ans}, or something similar. 
#' \code{result} is an expression that returns a \code{glm} object from 
#' the success solver. If the \code{solver} does not return a \code{glm} 
#' object, it is the duty of \code{result} expression to convert it to a 
#' \code{glm} object. 
#' @param ... Additional arguments passed to solvers. 
#' 
#' @return If at least one of the solvers succeed, an \code{glm} object
#' will be returned from the success. If all solvers fail but at least 
#' one solver did not throw an error, the result from the last solver 
#' that did not throw an error will be returned. In this case, 
#' \code{glmsolve} will throw a warning, stating that none of the solvers
#' succeeded. If all solvers ended up throwing errors, \code{glmsolve} 
#' will also throw an error, again stating that non of the solvers 
#' succeeded.
#' 
#' @author Long Qu
#' @keywords iteration
#' @seealso \code{\link[stats]{glm}}, \code{\link[stats]{nlminb}}, 
#' \code{\link[stats]{optim}}, \code{\link{nlsolve}}
#' 
#' @examples
#' ## Taken from stats::glm:
#' counts <- c(18,17,15,20,10,20,25,13,12)
#' outcome <- gl(3,1,9)
#' treatment <- gl(3,3)
#' print(d.AD <- data.frame(treatment, outcome, counts))
#' glm.D93 <- glmsolve(counts ~ outcome + treatment, family = poisson())
#' 
#' @export 
#' @name glmSolve
glmsolve=function(formula, family=gaussian, data, control=list(...), solvers=glmSolvers, ...)
{
	nsolvers=length(solvers)
	if(missing(data)) data=parent.frame()
	if(nsolvers==0L || !is.list(solvers))stop('solvers must be supplied as a list')
	solver.names=names(solvers)
	if(is.null(solver.names)) solver.names=paste('anonymous',seq_len(nsolvers),sep='.')
	tmp=solver.names==''
	if(any(tmp)) solver.names[tmp]=paste('anonymous',seq_len(sum(tmp)),sep='.')
	
	thisEnv=sys.frame(sys.nframe())
	lastNonError=NULL
	for(i in seq_len(nsolvers)){
		solver.name=solver.names[i]
		ans=try( eval(solvers[[i]]$solve, envir=thisEnv), silent=TRUE)
		if(inherits(ans,'try-error')) next 
		success=eval(solvers[[i]]$success, envir=thisEnv)
		ans=eval(solvers[[i]]$result, envir=thisEnv)
		attr(ans, 'glmsolve.method')=solver.name
		if(isTRUE(success)){
			attr(ans, 'glmsolve.success')=TRUE
			return(ans)
		}else {
			lastNonError=ans
			attr(lastNonError, 'glmsolve.success')=FALSE
		}
	}
	
	if(inherits(family, 'fbrNBfamily')) {
		family=tryCatch(reinitialize.fbrNBfamily(family), simpleWarning=function(w)w$message)
		return(Recall(formula, family, data, control, solvers, ...))
	}
	
	do.call(if(is.null(lastNonError)) 'stop' else 'warning', list('None of the glm solvers succeeded.'))
	lastNonError
}
