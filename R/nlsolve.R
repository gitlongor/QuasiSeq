#' @rdname nlSolve
#' @export
#' @format \code{nlSolvers} is currently a named list of length 8. 
#' @examples
#' str(nlSolvers) ## list of existing solvers
nlSolvers=list(
	broyden=list(solve=expression({J0=attr(func, 'getApproxJacob')(); pracma::broyden(func, start, J0=J0,maxiter=control$maxit, tol=control$tol, ...)}),
				 success=expression(ans$fnorm<1e-6),
				 result=expression(structure(ans$zero, iter=ans$niter))),
	nleqslv.Broyden=list(solve=expression({jac=attr(func, 'getApproxJacob'); nleqslv::nleqslv(start, func, jac=jac, method='Broyden', control=list(trace=0, xtol=control$tol, ftol=control$tol, maxit=control$maxit),...)}),
				 success=expression(ans$termcd<=2), 
				 result=expression(structure(ans$x, iter=ans$iter))),
	nleqslv.Newton=list(solve=expression(nleqslv::nleqslv(start, func, method='Newton',control=list(trace=0, xtol=control$tol, ftol=control$tol, maxit=control$maxit),...)),
				 success=expression(ans$termcd<=2), 
				 result=expression(structure(ans$x, iter=ans$iter))),
	fsolve=list(solve=expression(pracma::fsolve(func, start,maxiter=control$maxit, tol=control$tol,...)), 
				 success=expression(max(abs(ans$fval))<1e-6),
				 result=expression(structure(ans$x, iter=NA_integer_))),
	newtonsys=list(solve=expression(pracma::newtonsys(func, start,maxiter=control$maxit, tol=control$tol,...)), 
				 success=expression(ans$fnorm<1e-6),
				 result=expression(structure(ans$zero, iter=ans$iter))),
	sane=list(solve=expression(BB::sane(start, func, quiet=TRUE, control=list(maxit=control$maxit, tol=control$tol, trace=FALSE), ...)),
				 success=expression(ans$convergence==0), 
				 result=expression(structure(ans$par, iter=ans$iter))),
	dfsane=list(solve=expression(BB::dfsane(start, func, quiet=TRUE, control=list(maxit=control$maxit, tol=control$tol, trace=FALSE), ...)),
				 success=expression(ans$convergence==0), 
				 result=expression(structure(ans$par, iter=ans$iter))),
	BBsolve=list(solve=expression(BB::BBsolve(start, func, quiet=TRUE, control=list(maxit=control$maxit, tol=control$tol, trace=FALSE),  ...)),
				 success=expression(ans$convergence==0), 
				 result=expression(structure(ans$par, iter=ans$iter)))
	## dfsane & BBsolve are not quiet
)

#' Nonlinear equation solvers
#' 
#' \code{nlsolve} sequentially tries multiple nonliear equation solvers 
#' until a solution is found or all solvers have been tried. 
#' 
#' Currently supported solvers are \code{c('broyden', 'nleqslv.Broyden',
#' 'nleqslv.Newton','fsolve','newtonsys','sane','dfsane','BBsolve')}. 
#' They are respectively \code{\link[pracma]{broyden}} in the \code{pracma}
#' package, \code{\link[nleqslv]{nleqslv}} in the \code{nleqslv} with 
#' \code{method='Broyden'} or with \code{method='Newton'}, functions
#' \code{\link[pracma]{fsolve}} and \code{\link[pracma]{newtonsys}} 
#' in the \code{pracma} package, and functions \code{\link[BB]{sane}}, 
#' \code{\link[BB]{dfsane}}, and \code{\link[BB]{BBsolve}} in the 
#' \code{BB} package.
#' 
#' @param start A numeric vector of starting values. 
#' @param  func An objective function of which the zeros are sought.
#' @param solvers A (named) list, with each element being itself 
#' a list with components named \code{solve}, \code{success} and 
#' \code{result}. See the \code{solvers} argument of \code{\link{glmsolve}}. 
#' @param control See the \code{control} argument of \code{\link{glmsolve}}. 
#' @param ... Additional arguments.
#'  
#' @return If at least one solver succeeds, a numeric vector of solutions
#' will be returned, with attribute \code{attr(*, 'nlsolve.metho')} 
#' set to a string indicating the name of the success solver being used, 
#' and with attribute \code{attr(*, 'nlsolve.success')} being a logical
#' scalar indicating whether any solver succeeded. When all solvers 
#' fail, the last one is returned, but with \code{nlsolve.success} being 
#' set to \code{FALSE}. 
#'  
#' @author Long Qu
#' @keywords iteration
#' @seealso\code{\link[pracma]{broyden}}, \code{\link[nleqslv]{nleqslv}},
#' \code{\link[pracma]{fsolve}},  \code{\link[pracma]{newtonsys}},
#' \code{\link[BB]{sane}}, \code{\link[BB]{dfsane}}, \code{\link[BB]{BBsolve}},
#' \code{\link{glmsolve}}.
#' 
#' @export 
#' @name nlSolve
nlsolve=function(start, func, solvers=nlSolvers, control, ...)
{
	force(start)
	force(func)
	force(control)
	nsolvers=length(solvers)
	if(nsolvers==0L || !is.list(solvers))stop('solvers must be supplied as a list')
	solver.names=names(solvers)
	if(is.null(solver.names)) solver.names=paste('anonymous',seq(nsolvers),sep='.')
	tmp=solver.names==''
	if(any(tmp)) solver.names[tmp]=paste('anonymous',seq(sum(tmp)),sep='.')
	
	thisEnv=sys.frame(sys.nframe())
	cur.obj=sqrt(mean(func(start)^2)); cur.method='start'
	for(i in seq(nsolvers)){
		solver.name=solver.names[i]
		ans=try( suppressWarnings(eval(solvers[[i]]$solve, envir=thisEnv)), silent=TRUE)
		if(inherits(ans,'try-error')) next
		success=eval(solvers[[i]]$success, envir=thisEnv)
		ans=eval(solvers[[i]]$result, envir=thisEnv)
		if(isTRUE(success)){
			attr(ans, 'nlsolve.method')=solver.name
			attr(ans, 'nlsolve.success')=TRUE
			return(ans)
		}else{
			new.obj=sqrt(mean(func(ans)^2))
			if(new.obj<cur.obj){
				start=ans
				cur.obj=new.obj
				cur.method=solver.name
			}
		}
	}
	warning('None of the non-linear equation solvers succeeded.')
	attr(start, 'nlsolve.method')=cur.method
	attr(start, 'nlsolve.success')=FALSE
	start
}
