\name{mockRNASeqData}
\docType{data}
\alias{mockRNASeqData}
\title{A Simulated RNA-Seq Data Set}
\description{
  This is a simulated RNA-Seq data set using a negative binomial model with 10000 genes and 8 experimental unit, under a balanced two-treatment comparison design.
}
\usage{mockRNASeqData}
\format{This is a list with the following components: 
\describe{
\item{counts}{This is a numeric data matrix with 10000 rows and 8 columns, containing counts for each gene (row) and each experimental unit (column). 
}
\item{treatment}{This is a factor with 2 levels, indicating the treatment group of each column of \code{counts}.
}
\item{design.matrix}{This is an example of design matrix corresponding to \code{treatment}.
}
\item{true.offset}{This is a numeric vector of normalizing factors actually used to simulate the data matrix. 
}
\item{estimated.offset}{This is a numeric vector of normalizing factors estimated from the data matrix, using the so-called "TMM" method.
}
\item{true.nbdisp}{This is a numeric vector of negative binomial over-dispersion parameters actually used to simulate the data. This is using the parameterization such that \code{true.nbdisp = 1/size}, where \code{size} is the parameter used in \code{\link{rnbinom}}. 
}
\item{estimated.nbdisp}{This is a numeric vector of estimated negative binomial over-dispersion parameters, using the "TrendedDisp" method from the \link[edgeR:estimateGLMTrendedDisp]{edgeR} package.
}
\item{ngenes}{Integer scalar 10000, the number of rows of \code{counts}.
}
\item{nsamples}{Integer scalar 8, the number of columns of \code{counts}.
}
\item{true.DEgenes}{An integer vector of length 3500, indicating the correct row indices of differentially expressed genes, i.e., rows whose means differ across the two treatments.
}
\item{true.foldChanges}{A numeric vector of length 3500, indicating the true ratio of means for each differentially expressed genes. 
}
\item{simulation.expression}{This is a \code{R} expression that was used to simulate the \code{mockRNASeqData} data set itself. \code{eval(mockRNASeqData$simulation.expression)} should generate an identical data set, except for the \code{simulation.expression} component itself.
}

}
}


\keyword{datasets}