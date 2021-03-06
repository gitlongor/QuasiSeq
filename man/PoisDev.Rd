% File src/library/QuasiSeq/man/PoisDev.Rd
\name{PoisDev}
\alias{PoisDev}
\title{Compute Poisson deviances (up to a constant) for given design matrix}

\description{A function called within \code{QL.fit} to compute Poisson deviances of each gene for given design matrix}

\usage{
PoisDev(counts,design,log.offset,print.progress=TRUE)
}

\arguments{
    \item{counts}{RNA-seq data matrix of integer expression counts.  Each row contains observations from a single gene. Each column contains observations from a single experimental unit.}

	\item{design}{A single element from the \code{design.list} argument given to \code{QL.fit}.}

	\item{log.offset}{A vector of log-scale, additive factors used to adjust estimated log-scale means for differences in library sizes across samples.  Commonly used offsets include,\code{log.offset=log(colSums(counts))} or \code{log.offset=log(apply(counts,2,quantile,.75))}.  The default setting in \code{QL.fit} makes no adjustment for library sizes (i.e. log.offset=0).}

	\item{print.progress}{logical. If \code{TRUE}, a text progress bar will be displayed during the running of this function.}
}
	
\value{list containing:
	\item{dev}{vector containing the deviance for each gene under a Poisson model fit to design matrix (or vector, for one-factor experiments) specified by \code{design}. This vector is passed along within the \code{QL.fit} function.}
  
	\item{means}{matrix of fitted mean values for each gene}

	\item{parms}{matrix of estimated coefficients for each gene}	
	
	\item{bartlett.epsilon}{A numeric vector containing the epsilon values that can be later used for Bartllet correction of likelihood ratio test.}
}

\author{Originally authored by Steve Lund \email{lundsp@gmail.com}; later modified by Long Qu \email{rtistician@gmail.com}}

\examples{
\dontrun{

pois.fit =
with(mockRNASeqData,
	PoisDev(counts, design.matrix, log(estimated.normalization))
)

}
}
\keyword{distribution}
\keyword{iteration}
\concept{RNA-seq}
\concept{quasi-likelihood}
\concept{differential expression}
\concept{Poisson}





