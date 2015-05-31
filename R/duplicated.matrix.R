duplicated.matrix = function (x, incomparables = FALSE, MARGIN = 1L, fromLast = FALSE, ...)
{
<<<<<<< HEAD
	if (!is.matrix(x) || !is.atomic(x) || !identical(incomparables, FALSE) || (MARGIN[1L]!=1L && MARGIN[1L]!=2L) || length(MARGIN)!=1L || anyNA(x) )
=======
	if (!is.matrix(x) || !is.atomic(x) || !identical(incomparables, FALSE) || (MARGIN[1L]!=1L && MARGIN[1L]!=2L) || length(MARGIN)!=1L || any(is.na(x)) )
>>>>>>> origin/master
		return(base::duplicated.matrix(x, incomparables, MARGIN, fromLast, ...))
	.Call(C_dupAtomMat, x, as.integer(MARGIN), as.logical(fromLast))
}

unique.matrix=function (x, incomparables = FALSE, MARGIN = 1, fromLast = FALSE, ...)
{
<<<<<<< HEAD
	if (!is.matrix(x) || !is.atomic(x) || !identical(incomparables, FALSE) || (MARGIN[1L]!=1L && MARGIN[1L]!=2L) || length(MARGIN)!=1L || anyNA(x) )
=======
	if (!is.matrix(x) || !is.atomic(x) || !identical(incomparables, FALSE) || (MARGIN[1L]!=1L && MARGIN[1L]!=2L) || length(MARGIN)!=1L || any(is.na(x)) )
>>>>>>> origin/master
		return(base::unique.matrix(x, incomparables, MARGIN, fromLast, ...))
	dups=.Call(C_dupAtomMat, x, as.integer(MARGIN), as.logical(fromLast))
	if(MARGIN==1L) x[!dups,,drop=FALSE] else x[,!dups,drop=FALSE]
}
