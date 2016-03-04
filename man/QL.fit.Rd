% File src/library/QuasiSeq/man/QL.fit.Rd
     \name{QL.fit}
     \alias{QL.fit}
     \title{Fit quasi-likelihood models to matrix of RNA-seq expression count data}
     \description{
       Fit a quasi-likelihood model to RNA-seq expression count data using the methods detailed in Lund, Nettleton, McCarthy, and Smyth (2012).}
     \usage{
    QL.fit(counts,design.list, test.mat=NULL, log.offset=NULL, 
           Model="NegBin", print.progress=TRUE, NBdisp="trend", 
           bias.fold.tolerance=1.10, ...)
     }
     \arguments{
      \item{counts}{RNA-seq data matrix of integer expression counts.  Each row contains observations from a single gene. Each column contains observations from a single experimental unit.}
      \item{design.list}{List of design matrices for the full model and reduced model(s). The first element of \code{design.list} must describe the overall full model, as this design is used to compute deviance residuals for estimating dispersion.  One-factor designs may be specified as vectors. The number of rows in each design matrix (or the length of each design vector) must be \code{ncol(counts)}.  Means are modeled with a log link function.}
	\item{test.mat}{T by 2 matrix dictating which designs are to be compared, where T is the total number of desired hypothesis tests for each gene.  Each row contains two integers, which provide the indices within \code{design.list} of the two designs to be compared.  If \code{test.mat} is not specified, the default is compare the first design (the full model) to each of the other designs provided in \code{design.list} in order (i.e. first design compared to second design, first design compared to third design, first design compared to fourth design, etc.).}
	\item{log.offset}{A vector of log-scale, additive factors used to adjust estimated log-scale means for differences in library sizes across samples.  Commonly used offsets include \code{log.offset=log(colSums(counts))} and \code{log.offset=log(apply(counts[rowSums(counts)!=0,],2,quantile,.75))}.  The default setting makes no adjustment for library sizes (i.e. log.offset=0).}
	\item{Model}{Must be one of "Poisson" or "NegBin", specifying use of a quasi-Poisson or quasi-negative binomial model, respectively.}
	\item{print.progress}{logical. If \code{TRUE}, updates are provided regard what gene (row number) is being analyzed.  Updates occur frequently to start then eventually occur every 5000 genes.}
      \item{NBdisp}{Used only when \code{Model="NegBin"}. Must be one of "trend", "common" or a vector of non-negative real numbers with length equal to \code{nrow(counts)}. Specifying \code{NBdisp="trend"} or \code{NBdisp="common"} will use {\link{estimateGLMTrendedDisp}} or \code{\link{estimateGLMCommonDisp}}, respectively, from the package \code{edgeR} to estimate negative binomial dispersion parameters for each gene.  Estimates obtained from other sources can be used by entering \code{NBdisp} as a vector containing the negative binomial dispersion value to use for each gene when fitting quasi-likelihood model.}
	\item{bias.fold.tolerance}{A numerical value no smaller than 1. If the bias reduction of maximum likelihood estimates of (log) fold change is likely to result in a ratio of fold changes greater than this value, then bias reduction will be performed on such genes. Setting \code{bias.fold.tolerance=Inf} will completely disable bias reduction;  setting \code{bias.fold.tolerance=1} will always perform bias reduction.  See \code{\link{NBDev}} for details. } 
	\item{...}{arguments to be passed to the function \code{\link{estimateGLMTrendedDisp}} or \code{\link{estimateGLMCommonDisp}} from the package \code{edgeR}.}
	}

\value{list containing:
	\item{"LRT"}{matrix providing unadjusted likelihood ratio test statistics.  Each column contains statistics from a single hypothesis test, applied separately to each gene.}
	\item{"LRT.Bart"}{The same as \code{LRT}, except that the LRT test statistics are adjusted by the Bartlett factor to better match their asymptotic moments.}
	\item{"phi.hat.dev"}{vector providing unshrunken, deviance-based estimates of quasi-dispersion (phi) for each gene.}
	\item{"phi.hat.pearson"}{vector providing unshrunken, Pearson estimates of quasi-dispersion (phi) for each gene.}
	\item{"mn.cnt"}{vector providing average count for each gene.}
	\item{"den.df"}{denominator degrees of freedom. Equal to the number of samples minus the number of fitted parameters in the full model, which is specified by the first element of \code{design.list}}
	\item{"num.df"}{vector of numerator degrees of freedom for each test, computed as the difference in the number of fitted parameters between the full and reduced models for each test.}
	\item{"Model"}{Either "Poisson" or "NegBin", specifying which model (quasi-Poisson or quasi-negative binomial, respectively) was used.}
	\item{"nb.disp"}{Only appears when \code{Model="NegBin"}.  Vector providing negative binomial dispersion parameter estimate used during fitting of quasi-negative binomial model for each gene.} 
	\item{fitted.values}{matrix of fitted mean values}
	\item{coefficients}{matrix of estimated coefficients for each gene.  Note that these are given on the log scale. (i.e. intercept coefficient reports log(average count) and non-intercept coefficients report estimated (and bias reduced, when appropriate) log fold-changes.)  %  Genes with at least one zero count and initial absolute coefficient estimates greater than three undergo the bias correction of Kosmidis & Firth (2009) to moderate extreme coefficient estimates.
	}	
}

\references{
Kosmidis and Firth (2009) "Bias reduction in exponential family nonlinear models" \emph{Biometrika}, \bold{96}, 793--804.

Lund, Nettleton, McCarthy and Smyth (2012) "Detecting differential expression in RNA-sequence data using quasi-likelihood with shrunken dispersion estimates" \emph{SAGMB}, \bold{11}(5).

McCarthy, Chen and Smyth (2012) "Differential expression analysis of multifactor RNA-Seq experiments with respect to biological variation" \emph{Nucleic Acids Res.} \bold{40}(10), 4288--97.

Cordeiro (1983) "Improved Likelihood Ratio Statistics for Generalized Linear Models" \emph{Journal of the Royal Statistical Society. Series B (Methodological)}, \bold{45}, 404--413.
}

\seealso{\code{\link{QL.results}}, \code{\link{NBDev}}, \code{\link{PoisDev}}, \code{\link{mockRNASeqData}}}

\author{Originally authored by Steve Lund \email{lundsp@gmail.com}; later modified by Long Qu \email{long.qu@wright.edu}}

\examples{ 
data(mockRNASeqData)

### Create list of designs describing model under null and alternative hypotheses
design.list=list(
	mockRNASeqData$design.matrix, 
	matrix(1, nrow=rep(mockRNASeqData$nsamples))
)

#\dontrun{
### Analyze using QL, QLShrink and QLSpline methods applied to quasi-Poisson model
fit = 
with(mockRNASeqData,
	QL.fit(counts, design.list,
	log.offset=log(estimated.normalization),
	 Model="Poisson")
)
results = QL.results(fit)

### How many significant genes at FDR=.05 from QLSpline method?
apply(results$Q.values[[3]]<.05,2,sum)

### Indexes for Top 10 most significant genes from QLSpline method
head(order(results$P.values[[3]]), 10)
#}


### Analyze using QL, QLShrink and QLSpline methods 
### applied to quasi-negative binomial model
fit2 = 
with(mockRNASeqData,
	QL.fit(counts, design.list,
		log.offset=log(estimated.normalization),
		Model="NegBin")
)
##########################################################
## Note: 'NBdisp' typically will not be specified when  ##
## calling QL.fit while analyzing real data. Providing  ##
## numeric values for 'NBdisp' prevents neg binomial    ##
## dispersions from being estimated from the data.      ##
##########################################################
results2<-QL.results(fit2)

### How many significant genes at FDR=.05 for QLSpline method?
apply(results2$Q.values[[3]]<.05,2,sum)

### Indexes for Top 10 most significant genes from QLShrink method
head(order(results2$P.values[[2]]), 10)

}


\keyword{models}
\keyword{regression}
\keyword{iteration} 
\keyword{multivariate}
\concept{RNA-seq}
\concept{quasi-likelihood}
\concept{differential expression}
\concept{negative binomial}
\concept{poisson}

