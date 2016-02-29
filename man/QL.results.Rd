% File src/library/QuasiSeq/man/QL.results.Rd
     \name{QL.results}
     \alias{QL.results}
     \title{Obtain p-values and q-values using results from \code{QL.fit}}
     \description{
       Obtain significance results for quasi-likelihood model fits to RNA-seq expression count data using the methods detailed in Lund, Nettleton, McCarthy, and Smyth (2012).}
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
	\item{"Q.values"}{list of matrices providing q-values for the QL, QLShrink and QLSpline methods, respectively.  The i^th column of each element of \code{qvals} corresponds to the hypothesis test assigned in the i^th row of \code{test.mat}.
 Q-values are computed using the methods of Nettleton et al. (2006) JABES 11, 337-356.}
	\item{"F.stat"}{list of matrices providing F-statistics for the QL, QLShrink and QLSpline methods, respectively.  The i^th column of each element of \code{F.stat} corresponds to the hypothesis test assigned in the i^th row of \code{test.mat}.}
	\item{"m0"}{matrix providing estimated number of true null hypotheses for each test(arranged by row) under each of the three methods(arranged by column). \code{m0} values are computed using the methods of Nettleton et al. (2006) JABES 11, 337-356.}
	\item{"d0"}{ vector containing estimated additional denominator degrees of freedom  gained from shrinking dispersion estimates in the QLShrink and QLSpline procedures, respectively.}
}

\author{Steve Lund \email{lundsp@gmail.com}}

\references{
Lund, Nettleton, McCarthy and Smyth (2012) "Detecting differential expression in RNA-sequence data using quasi-likelihood with shrunken dispersion estimates" \emph{SAGMB}, \bold{11}(5).
}


\examples{  
set.seed(234092)

n.genes<-100
n.de<-round(.5*n.genes)
trt<-rep(1:2,each=4)
n.samp<-length(trt)
mu<-rgamma(n.genes,1.5,.01)

## specify gene specific negative binomial dispersions
size<-(log(mu+exp(1))-1)/mu  ### Var(Y)=E(Y)log(E(Y)+exp(1))

## add noise to gene specific negative binomial dispersions
size<-size*4.5/rchisq(n.genes,4.5)

sim.mn<-matrix(mu,n.genes,2)

### Simulate fold changes
B<-exp((2*rbinom(n.de,1,.5)-1)*(.25+rbeta(n.de,1,2)))

sim.mn[1:n.de,1]<-sim.mn[1:n.de,1]*B^(.5)+5
sim.mn[1:n.de,2]<-sim.mn[1:n.de,2]*B^(-.5)

### Simulate library size factors
sim.offset<-2^(rnorm(n.samp,0,.15))

### Compute final means
sim.mn2<-t(t(sim.mn[,trt])*sim.offset)

### Simulate data
simdat<-matrix(rnbinom(n.samp*n.genes,mu=sim.mn2,size=1/size),n.genes,n.samp)

### Simulate estimated dispersions to save time
################################################################
## THIS STEP SHOULD NOT BE PERFORMED WHEN ANALYZING REAL DATA ##
################################################################
est.nb.disp<-size*rchisq(n.genes,n.samp-2)/(n.samp-2)
est.nb.disp<-est.nb.disp

### Keep genes with at least 10 total counts
est.nb.disp<-est.nb.disp[rowSums(simdat)>9]
simdat<-simdat[rowSums(simdat)>9,]


### Create list of designs describing model under null and alternative hypotheses
design.list<-vector("list",2)
design.list[[1]]<-model.matrix(~as.factor(trt))  #This also could have just been ``trt''.
design.list[[2]]<-rep(1,length(trt))

log.offset<-log(apply(simdat,2,quantile,.75))
\dontrun{
### Analyze using QL, QLShrink and QLSpline methods applied to quasi-Poisson model
fit<-QL.fit(simdat, design.list,log.offset=log.offset, Model="Poisson")

results<-QL.results(fit)

### How many significant genes at FDR=.05 from QLSpline method?
apply(results$Q.values[[3]]<.05,2,sum)

### Indexes for Top 10 most significant genes from QLSpline method
head(order(results$P.values[[3]]), 10)
}
### Analyze using QL, QLShrink and QLSpline methods 
### applied to quasi-negative binomial model
fit2<-QL.fit(simdat, design.list,log.offset=log.offset, nb.disp=est.nb.disp, Model="NegBin")

##########################################################
## Note: 'nb.disp' typically will not be specified when ##
## calling QL.fit while analyzing real data. Providing  ##
## numeric values for 'nb.disp' prevents neg binomial   ##
## dispersions from being estimated from the data.      ##
##########################################################

results2<-QL.results(fit2)

### How many significant genes at FDR=.05 for QLSpline method?
apply(results2$Q.values[[3]]<.05,2,sum)

### Indexes for Top 10 most significant genes from QLShrink method
head(order(results2$P.values[[2]]), 10)
}
     \keyword{RNA-seq, quasi-likelihood, differential expression}



