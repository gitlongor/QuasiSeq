duplicated.matrix = function (x, incomparables = FALSE, MARGIN = 1L, fromLast = FALSE, ...)
{
	if (!is.matrix(x) || !(is.numeric(x) || is.logical(x)) || !identical(incomparables, FALSE) || (MARGIN!=1L && MARGIN!=2L) )
		return(base::duplicated.matrix(x, incomparables, MARGIN, fromLast, ...))
	.Call(C_dupNumMat, x, as.integer(MARGIN), as.logical(fromLast))
}

unique.matrix=function (x, incomparables = FALSE, MARGIN = 1, fromLast = FALSE, ...)
{
	if (!is.matrix(x) || !(is.numeric(x) || is.logical(x)) || !identical(incomparables, FALSE) || (MARGIN!=1L && MARGIN!=2L) )
		return(base::unique.matrix(x, incomparables, MARGIN, fromLast, ...))
	if(MARGIN==1){
		x[!.Call(C_dupNumMat, x, as.integer(MARGIN), as.logical(fromLast)),,drop=FALSE]
	}else if(MARGIN==2){
		x[,!.Call(C_dupNumMat, x, as.integer(MARGIN), as.logical(fromLast)),drop=FALSE]
	}else stop('MARGIN should be 1L or 2L')
}
