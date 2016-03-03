% File src/library/QuasiSeq/man/NBDev.Rd
\name{NBDev}
\alias{NBDev}
\title{Fit a negative binomial GLM for given design matrix}
\description{A function called within \code{QL.fit} to fit a negative binomial GLM of each gene for given design matrix}

\usage{
NBDev(counts,design,log.offset,nb.disp,print.progress=TRUE, bias.fold.tolerance=1.10)
}

\arguments{
      \item{counts}{RNA-seq data matrix of integer expression counts.  Each row contains observations from a single gene. Each column contains observations from a single experimental unit.}
      \item{design}{A single element from the \code{design.list} argument given to \code{QL.fit}.}
	\item{log.offset}{A vector of log-scale, additive factors used to adjust estimated log-scale means for differences in library sizes across samples.  Commonly used offsets include,\code{log.offset=log(colSums(counts))} or \code{log.offset=log(apply(counts,2,quantile,.75))}.  The default setting  in \code{QLfit} makes no adjustment for library sizes (i.e. log.offset=0).}
	\item{nb.disp}{estimated negative binomial dispersion parameters obtained from either \code{estimateGLMTrendedDisp} or \code{estimateGLMCommonDisp} in package \code{edgeR}.  These estimates are treated as known and are used to compute deviances.}
	\item{print.progress}{logical. If TRUE, the function will provide an update on what gene (row number) is being analyzed.  Updates occur frequently to start then eventually occur every 5000 genes.}
	\item{bias.fold.tolerance}{A numerical value no smaller than 1. If the bias reduction of maximum likelihood estimates of (log) fold change is likely to result in a ratio of fold changes greater than this value, then bias reduction will be performed on such genes. Setting \code{bias.fold.tolerance=Inf} will completely disable bias reduction;  setting \code{bias.fold.tolerance=1} will always perform bias reduction. See details. }
}
	
\value{list containing:
\item{dev}{vector containing the deviance for each gene under a negative binomial model fit to design matrix specified by \code{design}. This vector is passed along within the \code{QL.fit} function.}
  \item{means}{matrix of fitted mean values for each gene}
\item{parms}{matrix of estimated coefficients for each gene. Note that these are given on the log scale. (i.e. intercept coefficient reports log(average count) and non-intercept coefficients report estimated (and bias reduced, when appropriate) log fold-changes.)  
%Genes with at least one zero count and initial absolute coefficient estimates greater than three undergo the bias correction of Kosmidis & Firth (2009) to moderate extreme coefficient estimates.
}	
}

\details{
This functions fits, for each row of \code{counts}, a negative binomial log-linear model through GLM framework with the over-dispersion parameter fixed. Asymptotic biases of regression coefficients (i.e., log fold changes) are then estimated by a plug-in estimate [eqn. (15.4) of McCullagh and Nelder, 1989] from the last iteration of iteratively re-weighted least squares (IWLS) procedure. The fitted response values are then compared with or without such a bias term. If the ratio of fitted response values are larger than \code{bias.fold.tolerance} for any observation, then the bias-reduction (not bias-correction) procedure according to Firth (1993) and Kosmidis & Firth (2009) is applied to such rows of \code{counts}, by adjusting the score equation with a term based on the observed information. Such bias-reduced estimates are more advantageous than directly subtracting the estimated bias from the maximum likelihood estimates as the latter may not exist (e.g., when all counts in one treatment group are zeros). 
}

\references{
Firth (1993) "Bias reduction of maximum likelihood estimates" \emph{Biometrika}, \bold{80}, 27--38. 

Kosmidis and Firth (2009) "Bias reduction in exponential family nonlinear models" \emph{Biometrika}, \bold{96}, 793--804.

Lund, Nettleton, McCarthy and Smyth (2012) "Detecting differential expression in RNA-sequence data using quasi-likelihood with shrunken dispersion estimates" emph{SAGMB}, \bold{11}(5).

McCullagh and Nelder (1989) "Generalized Linear Models", 2nd edition. London: Chapman and Hall. 
}

\examples{
#\dontrun{

## no bias reduction
noReduction = 
with(mockRNASeqData,
	NBDev(counts, design.matrix, log(estimated.normalization),
		  estimated.nbdisp, bias.fold.tolerance=Inf)
)

## bias reduction for all genes
allReduction = 
with(mockRNASeqData,
	NBDev(counts, design.matrix, log(estimated.normalization),
		  estimated.nbdisp, bias.fold.tolerance=1)
)

## default: bias reduction for genes showing large bias
someReduction = 
with(mockRNASeqData,
	NBDev(counts, design.matrix, log(estimated.normalization),
		  estimated.nbdisp, bias.fold.tolerance=1.1)
)

#}
}
\author{Steve Lund (\email{lundsp@gmail.com}), Long Qu (\email{long.qu@wright.edu})}

\keyword{distribution, iteration}




