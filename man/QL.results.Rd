% File src/library/QuasiSeq/man/QL.results.Rd
\name{QL.results}
\alias{QL.results}
\title{Obtain p-values and q-values using results from \code{QL.fit}}
\description{
   Obtain significance results for quasi-likelihood model fits to RNA-seq expression count data using the methods detailed in Lund, Nettleton, McCarthy, and Smyth (2012).
}

\usage{	
QL.results(fit,Dispersion="Deviance",spline.df=NULL,Plot=TRUE)
}

\arguments{
    \item{fit}{The list returned by the function QL.fit}

	\item{Dispersion}{Must be one of "Deviance" or "Pearson", specifying which type of estimator should be used for estimating quasi-likelihood dispersion parameter.}

	\item{spline.df}{Optional.  User may specify the degrees of freedom to use when fitting a cubic spline to log-scale(estimated dispersion) versus the log(average count).  Default uses cross-validation in \code{sreg} function to pick appropriate degrees of freedom.} 	

	\item{Plot}{logical.  If TRUE, the estimated dispersion versus the average count are plotted on a log-scale with the corresponding cubic spline fit overlaid. }
}

\value{list containing:
	\item{"P.values"}{list of matrices providing p-values for the QL, QLShrink and QLSpline methods, respectively.  The i^th column of each element of \code{pvals} corresponds to the hypothesis test assigned in the i^th row of \code{test.mat}.}

	\item{"Q.values"}{list of matrices providing q-values for the QL, QLShrink and QLSpline methods, respectively.  The i^th column of each element of \code{qvals} corresponds to the hypothesis test assigned in the i^th row of \code{test.mat}. Q-values are computed using the methods of Nettleton et al. (2006) JABES 11, 337-356.}
	
	\item{"F.stat"}{list of matrices providing F-statistics for the QL, QLShrink and QLSpline methods, respectively.  The i^th column of each element of \code{F.stat} corresponds to the hypothesis test assigned in the i^th row of \code{test.mat}.}
	
	\item{"m0"}{matrix providing estimated number of true null hypotheses for each test(arranged by row) under each of the three methods(arranged by column). \code{m0} values are computed using the methods of Nettleton et al. (2006) JABES 11, 337-356.}
	
	\item{"d0"}{ vector containing estimated additional denominator degrees of freedom  gained from shrinking dispersion estimates in the QLShrink and QLSpline procedures, respectively.}
}

\author{Steve Lund \email{lundsp@gmail.com}}

\references{
Lund, Nettleton, McCarthy and Smyth (2012) "Detecting differential expression in RNA-sequence data using quasi-likelihood with shrunken dispersion estimates" \emph{SAGMB}, \bold{11}(5).
}

\seealso{\code{\link{QL.fit}}, \code{\link{NBDev}}, \code{\link{mockRNASeqData}}}

\examples{  
## see examples for QL.fit()
}

\keyword{htest}
\keyword{regression}
\concept{RNA-seq}
\concept{quasi-likelihood}
\concept{differential expression}
\concept{negative binomial}
\concept{poisson}



