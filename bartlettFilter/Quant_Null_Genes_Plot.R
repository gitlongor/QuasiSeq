print(getwd())
setwd("C:/Users/Klirk/Desktop/Simulation_Code/Simulation_Code/KIRC_Simulations/bartlettFilter")
### Load CRAN packages
require(ggplot2)

###Non-Cooks analysis

QuantThresh <- function(pvals, null.genes, specificity){
  quantile(pvals[null.genes], specificity)
}

k.ind <- 5
source("benchmarks_quant.R")
#num.deseq2.simseq <- apply(pvals.deseq2.simseq, 1, QuantThresh, null.genes, 0.05)
#num.edger.simseq <- apply(pvals.edger.simseq, 1, QuantThresh, null.genes, 0.05)
num.quasiseq.simseq.ql <- apply(pvals.quasiseq.simseq.ql, 1, QuantThresh, null.genes, 0.05)
num.quasiseq.simseq.spline <- apply(pvals.quasiseq.simseq.spline, 1, QuantThresh, null.genes, 0.05)
num.quasiseq.simseq.bart <- apply(pvals.quasiseq.simseq.bart, 1, QuantThresh, null.genes, 0.05)
#num.samseq.simseq <- apply(pvals.samseq.simseq, 1, QuantThresh, null.genes, 0.05)
num.voom.simseq <- apply(pvals.voom.simseq, 1, QuantThresh, null.genes, 0.05)
#num.deseq2.nb <- apply(pvals.deseq2.nb, 1, QuantThresh, null.genes, 0.05)
#num.edger.nb <- apply(pvals.edger.nb, 1, QuantThresh, null.genes, 0.05)
num.quasiseq.nb.ql <- apply(pvals.quasiseq.nb.ql, 1, QuantThresh, null.genes, 0.05)
num.quasiseq.nb.spline <- apply(pvals.quasiseq.nb.spline, 1, QuantThresh, null.genes, 0.05)
num.quasiseq.nb.bart <- apply(pvals.quasiseq.nb.bart, 1, QuantThresh, null.genes, 0.05)
#num.samseq.nb <- apply(pvals.samseq.nb, 1, QuantThresh, null.genes, 0.05)
num.voom.nb <- apply(pvals.voom.nb, 1, QuantThresh, null.genes, 0.05)

l <- length(num.quasiseq.simseq.ql)
num1 <- data.frame(num = c(num.quasiseq.simseq.ql, num.quasiseq.simseq.spline, num.quasiseq.simseq.bart, num.voom.simseq), 
                   Method = as.factor(c(rep("QuasiSeqQL", l), rep("QuasiSeqQLSpline", l), rep("QuasiSeqBart", l), rep("Voom", l))),
                   Simulation = rep("SimSeq", l*4),
                   SS = rep("Sample Size: 5",l*4))

l <- length(num.quasiseq.nb.ql)
num2 <- data.frame(num = c(num.quasiseq.nb.ql, num.quasiseq.nb.spline, num.quasiseq.nb.bart, num.voom.nb), 
                   Method = as.factor(c(rep("QuasiSeqQL", l), rep("QuasiSeqQLSpline", l), rep("QuasiSeqBart", l), rep("Voom", l))),
                   Simulation = rep("NB", l*4),
                   SS = rep("Sample Size: 5",l*4))

k.ind <- 10
source("benchmarks_quant.R")
#num.deseq2.simseq <- apply(pvals.deseq2.simseq, 1, QuantThresh, null.genes, 0.05)
#num.edger.simseq <- apply(pvals.edger.simseq, 1, QuantThresh, null.genes, 0.05)
num.quasiseq.simseq.ql <- apply(pvals.quasiseq.simseq.ql, 1, QuantThresh, null.genes, 0.05)
num.quasiseq.simseq.spline <- apply(pvals.quasiseq.simseq.spline, 1, QuantThresh, null.genes, 0.05)
num.quasiseq.simseq.bart <- apply(pvals.quasiseq.simseq.bart, 1, QuantThresh, null.genes, 0.05)
#num.samseq.simseq <- apply(pvals.samseq.simseq, 1, QuantThresh, null.genes, 0.05)
num.voom.simseq <- apply(pvals.voom.simseq, 1, QuantThresh, null.genes, 0.05)
#num.deseq2.nb <- apply(pvals.deseq2.nb, 1, QuantThresh, null.genes, 0.05)
#num.edger.nb <- apply(pvals.edger.nb, 1, QuantThresh, null.genes, 0.05)
num.quasiseq.nb.ql <- apply(pvals.quasiseq.nb.ql, 1, QuantThresh, null.genes, 0.05)
num.quasiseq.nb.spline <- apply(pvals.quasiseq.nb.spline, 1, QuantThresh, null.genes, 0.05)
num.quasiseq.nb.bart <- apply(pvals.quasiseq.nb.bart, 1, QuantThresh, null.genes, 0.05)
#num.samseq.nb <- apply(pvals.samseq.nb, 1, QuantThresh, null.genes, 0.05)
num.voom.nb <- apply(pvals.voom.nb, 1, QuantThresh, null.genes, 0.05)

l <- length(num.quasiseq.simseq.ql)
num3 <- data.frame(num = c(num.quasiseq.simseq.ql, num.quasiseq.simseq.spline, num.quasiseq.simseq.bart, num.voom.simseq), 
                   Method = as.factor(c(rep("QuasiSeqQL", l), rep("QuasiSeqQLSpline", l), rep("QuasiSeqBart", l), rep("Voom", l))),
                   Simulation = rep("SimSeq", l*4),
                   SS = rep("Sample Size: 10",l*4))

l <- length(num.quasiseq.nb.ql)
num4 <- data.frame(num = c(num.quasiseq.nb.ql, num.quasiseq.nb.spline, num.quasiseq.nb.bart, num.voom.nb), 
                   Method = as.factor(c(rep("QuasiSeqQL", l), rep("QuasiSeqQLSpline", l), rep("QuasiSeqBart", l), rep("Voom", l))),
                   Simulation = rep("NB", l*4),
                   SS = rep("Sample Size: 10",l*4))

k.ind <- 20
source("benchmarks_quant.R")
#num.deseq2.simseq <- apply(pvals.deseq2.simseq, 1, QuantThresh, null.genes, 0.05)
#num.edger.simseq <- apply(pvals.edger.simseq, 1, QuantThresh, null.genes, 0.05)
num.quasiseq.simseq.ql <- apply(pvals.quasiseq.simseq.ql, 1, QuantThresh, null.genes, 0.05)
num.quasiseq.simseq.spline <- apply(pvals.quasiseq.simseq.spline, 1, QuantThresh, null.genes, 0.05)
num.quasiseq.simseq.bart <- apply(pvals.quasiseq.simseq.bart, 1, QuantThresh, null.genes, 0.05)
#num.samseq.simseq <- apply(pvals.samseq.simseq, 1, QuantThresh, null.genes, 0.05)
num.voom.simseq <- apply(pvals.voom.simseq, 1, QuantThresh, null.genes, 0.05)
#num.deseq2.nb <- apply(pvals.deseq2.nb, 1, QuantThresh, null.genes, 0.05)
#num.edger.nb <- apply(pvals.edger.nb, 1, QuantThresh, null.genes, 0.05)
num.quasiseq.nb.ql <- apply(pvals.quasiseq.nb.ql, 1, QuantThresh, null.genes, 0.05)
num.quasiseq.nb.spline <- apply(pvals.quasiseq.nb.spline, 1, QuantThresh, null.genes, 0.05)
num.quasiseq.nb.bart <- apply(pvals.quasiseq.nb.bart, 1, QuantThresh, null.genes, 0.05)
#num.samseq.nb <- apply(pvals.samseq.nb, 1, QuantThresh, null.genes, 0.05)
num.voom.nb <- apply(pvals.voom.nb, 1, QuantThresh, null.genes, 0.05)

l <- length(num.quasiseq.simseq.ql)
num5 <- data.frame(num = c(num.quasiseq.simseq.ql, num.quasiseq.simseq.spline, num.quasiseq.simseq.bart, num.voom.simseq), 
                   Method = as.factor(c(rep("QuasiSeqQL", l), rep("QuasiSeqQLSpline", l), rep("QuasiSeqBart", l), rep("Voom", l))),
                   Simulation = rep("SimSeq", l*4),
                   SS = rep("Sample Size: 20",l*4))

l <- length(num.quasiseq.nb.ql)
num6 <- data.frame(num = c(num.quasiseq.nb.ql, num.quasiseq.nb.spline, num.quasiseq.nb.bart, num.voom.nb), 
                   Method = as.factor(c(rep("QuasiSeqQL", l), rep("QuasiSeqQLSpline", l), rep("QuasiSeqBart", l), rep("Voom", l))),
                   Simulation = rep("NB", l*4),
                   SS = rep("Sample Size: 20",l*4))

num <- rbind(num1, num2, num3, num4, num5, num6)
qplot(Simulation, num, data = num, geom = "boxplot", facets = SS ~ Method, 
      ylab = "0.05 Quantile of p-values from null genes") + theme_bw()+ 
  theme(axis.text.x = element_text(angle=90))

if(FALSE){
### Cook's Filtering
QuantThreshCooks <- function(pvals, filt, null.genes, specificity, pvals.list = FALSE)
{
  if(pvals.list){
    n.iter <- length(pvals)
  } else {
    n.iter <- nrow(pvals)
  }
  res <- rep(NA, n.iter)
  for(i in 1:n.iter)
  {
    if(pvals.list){
      pvals.temp <- pvals[[i]]
    } else {
      pvals.temp <- pvals[i, ][filt[i, ]] 
    }
    null.genes.temp <- null.genes[filt[i, ]]
    res[i] <- QuantThresh(pvals.temp, null.genes.temp, specificity)
  }
  return(res)
}

### Cooks Analysis
k.ind <- 5
source("benchmarks_quant.R")
num.deseq2.simseq <- QuantThreshCooks(pvals.deseq2.cooks.simseq, filt.cooks.simseq, null.genes, specificity = 0.05, pvals.list = TRUE)
num.edger.simseq <- QuantThreshCooks(pvals.edger.cooks.simseq, filt.cooks.simseq, null.genes, specificity = 0.05, pvals.list = TRUE)
num.quasiseq.simseq <- QuantThreshCooks(pvals.quasiseq.cooks.simseq, filt.cooks.simseq, null.genes, specificity = 0.05, pvals.list = TRUE)
num.quasiseq.simseq.bart <- QuantThreshCooks(pvals.quasiseq.cooks.simseq.bart, filt.cooks.simseq, null.genes, specificity = 0.05, pvals.list = TRUE)
num.samseq.simseq <- QuantThreshCooks(pvals.samseq.cooks.simseq, filt.cooks.simseq, null.genes, specificity = 0.05, pvals.list = TRUE)
num.voom.simseq <- QuantThreshCooks(pvals.voom.cooks.simseq, filt.cooks.simseq, null.genes, specificity = 0.05, pvals.list = TRUE)
num.deseq2.nb <- QuantThreshCooks(pvals.deseq2.cooks.nb, filt.cooks.nb, null.genes, specificity = 0.05, pvals.list = TRUE)
num.edger.nb <- QuantThreshCooks(pvals.edger.cooks.nb, filt.cooks.nb, null.genes, specificity = 0.05, pvals.list = TRUE)
num.quasiseq.nb <- QuantThreshCooks(pvals.quasiseq.cooks.nb, filt.cooks.nb, null.genes, specificity = 0.05, pvals.list = TRUE)
num.quasiseq.nb.bart <- QuantThreshCooks(pvals.quasiseq.cooks.nb.bart, filt.cooks.nb, null.genes, specificity = 0.05, pvals.list = TRUE)
num.samseq.nb <- QuantThreshCooks(pvals.samseq.cooks.nb, filt.cooks.nb, null.genes, specificity = 0.05, pvals.list = TRUE)
num.voom.nb <- QuantThreshCooks(pvals.voom.cooks.nb, filt.cooks.nb, null.genes, specificity = 0.05, pvals.list = TRUE)

l <- length(num.deseq2.simseq)
num1 <- data.frame(num = c(num.deseq2.simseq, num.edger.simseq, num.quasiseq.simseq,num.quasiseq.simseq.bart, num.samseq.simseq, num.voom.simseq), 
                   Method = as.factor(c(rep("DESeq2", l), rep("EdgeR", l),
                                        rep("QuasiSeq", l),rep("QuasiSeqBart", l), rep("SAMseq", l),rep("Voom", l))),
                   Simulation = rep("SimSeq", l*6),
                   SS = rep("Sample Size: 5",l*6))

l <- length(num.deseq2.nb)
num2 <- data.frame(num = c(num.deseq2.nb, num.edger.nb, num.quasiseq.nb,num.quasiseq.nb.bart, num.samseq.nb, num.voom.nb), 
                   Method = as.factor(c(rep("DESeq2", l), rep("EdgeR", l),
                                        rep("QuasiSeq", l),rep("QuasiSeqBart", l), rep("SAMseq", l),rep("Voom", l))),
                   Simulation = rep("NB", l*6),
                   SS = rep("Sample Size: 5",l*6))

k.ind <- 10
source("benchmarks_quant.R")
num.deseq2.simseq <- QuantThreshCooks(pvals.deseq2.simseq, filt.cooks.simseq, null.genes, specificity = 0.05, pvals.list = FALSE)
num.edger.simseq <- QuantThreshCooks(pvals.edger.simseq, filt.cooks.simseq, null.genes, specificity = 0.05, pvals.list = FALSE)
num.quasiseq.simseq <- QuantThreshCooks(pvals.quasiseq.simseq, filt.cooks.simseq, null.genes, specificity = 0.05, pvals.list = FALSE)
num.quasiseq.simseq.bart <- QuantThreshCooks(pvals.quasiseq.simseq.bart, filt.cooks.simseq, null.genes, specificity = 0.05, pvals.list = FALSE)
num.samseq.simseq <- QuantThreshCooks(pvals.samseq.simseq, filt.cooks.simseq, null.genes, specificity = 0.05, pvals.list = FALSE)
num.voom.simseq <- QuantThreshCooks(pvals.voom.simseq, filt.cooks.simseq, null.genes, specificity = 0.05, pvals.list = FALSE)
num.deseq2.nb <- QuantThreshCooks(pvals.deseq2.nb, filt.cooks.nb, null.genes, specificity = 0.05, pvals.list = FALSE)
num.edger.nb <- QuantThreshCooks(pvals.edger.nb, filt.cooks.nb, null.genes, specificity = 0.05, pvals.list = FALSE)
num.quasiseq.nb <- QuantThreshCooks(pvals.quasiseq.nb, filt.cooks.nb, null.genes, specificity = 0.05, pvals.list = FALSE)
num.quasiseq.nb.bart <- QuantThreshCooks(pvals.quasiseq.nb.bart, filt.cooks.nb, null.genes, specificity = 0.05, pvals.list = FALSE)
num.samseq.nb <- QuantThreshCooks(pvals.samseq.nb, filt.cooks.nb, null.genes, specificity = 0.05, pvals.list = FALSE)
num.voom.nb <- QuantThreshCooks(pvals.voom.nb, filt.cooks.nb, null.genes, specificity = 0.05, pvals.list = FALSE)

l <- length(num.deseq2.simseq)
num3 <- data.frame(num = c(num.deseq2.simseq, num.edger.simseq, num.quasiseq.simseq,num.quasiseq.simseq.bart, num.samseq.simseq, num.voom.simseq), 
                   Method = as.factor(c(rep("DESeq2", l), rep("EdgeR", l),
                                        rep("QuasiSeq", l),rep("QuasiSeqBart", l), rep("SAMseq", l),rep("Voom", l))),
                   Simulation = rep("SimSeq", l*6),
                   SS = rep("Sample Size: 10",l*6))

l <- length(num.deseq2.nb)
num4 <- data.frame(num = c(num.deseq2.nb, num.edger.nb, num.quasiseq.nb, num.quasiseq.nb.bart,num.samseq.nb, num.voom.nb), 
                   Method = as.factor(c(rep("DESeq2", l), rep("EdgeR", l),
                                        rep("QuasiSeq", l),rep("QuasiSeqBart", l), rep("SAMseq", l),rep("Voom", l))),
                   Simulation = rep("NB", l*6),
                   SS = rep("Sample Size: 10",l*6))

k.ind <- 20
source("benchmarks_quant.R")
num.deseq2.simseq <- QuantThreshCooks(pvals.deseq2.simseq, filt.cooks.simseq, null.genes, specificity = 0.05, pvals.list = FALSE)
num.edger.simseq <- QuantThreshCooks(pvals.edger.simseq, filt.cooks.simseq, null.genes, specificity = 0.05, pvals.list = FALSE)
num.quasiseq.simseq <- QuantThreshCooks(pvals.quasiseq.simseq, filt.cooks.simseq, null.genes, specificity = 0.05, pvals.list = FALSE)
num.quasiseq.simseq.bart <- QuantThreshCooks(pvals.quasiseq.simseq.bart, filt.cooks.simseq, null.genes, specificity = 0.05, pvals.list = FALSE)
num.samseq.simseq <- QuantThreshCooks(pvals.samseq.simseq, filt.cooks.simseq, null.genes, specificity = 0.05, pvals.list = FALSE)
num.voom.simseq <- QuantThreshCooks(pvals.voom.simseq, filt.cooks.simseq, null.genes, specificity = 0.05, pvals.list = FALSE)
num.deseq2.nb <- QuantThreshCooks(pvals.deseq2.nb, filt.cooks.nb, null.genes, specificity = 0.05, pvals.list = FALSE)
num.edger.nb <- QuantThreshCooks(pvals.edger.nb, filt.cooks.nb, null.genes, specificity = 0.05, pvals.list = FALSE)
num.quasiseq.nb <- QuantThreshCooks(pvals.quasiseq.nb, filt.cooks.nb, null.genes, specificity = 0.05, pvals.list = FALSE)
num.quasiseq.nb.bart <- QuantThreshCooks(pvals.quasiseq.nb.bart, filt.cooks.nb, null.genes, specificity = 0.05, pvals.list = FALSE)
num.samseq.nb <- QuantThreshCooks(pvals.samseq.nb, filt.cooks.nb, null.genes, specificity = 0.05, pvals.list = FALSE)
num.voom.nb <- QuantThreshCooks(pvals.voom.nb, filt.cooks.nb, null.genes, specificity = 0.05, pvals.list = FALSE)

l <- length(num.deseq2.simseq)
num5 <- data.frame(num = c(num.deseq2.simseq, num.edger.simseq, num.quasiseq.simseq, num.quasiseq.simseq.bart,num.samseq.simseq, num.voom.simseq), 
                   Method = as.factor(c(rep("DESeq2", l), rep("EdgeR", l),
                                        rep("QuasiSeq", l),rep("QuasiSeqBart", l), rep("SAMseq", l),rep("Voom", l))),
                   Simulation = rep("SimSeq", l*6),
                   SS = rep("Sample Size: 20",l*6))

l <- length(num.deseq2.nb)
num6 <- data.frame(num = c(num.deseq2.nb, num.edger.nb, num.quasiseq.nb,num.quasiseq.nb.bart, num.samseq.nb, num.voom.nb), 
                   Method = as.factor(c(rep("DESeq2", l), rep("EdgeR", l),
                                        rep("QuasiSeq", l),rep("QuasiSeqBart", l), rep("SAMseq", l),rep("Voom", l))),
                   Simulation = rep("NB", l*6),
                   SS = rep("Sample Size: 20",l*6))

num <- rbind(num1, num2, num3, num4, num5, num6)
qplot(Simulation, num, data = num, geom = "boxplot", facets = SS ~ Method, 
      ylab = "0.05 Quantile of p-values from null genes") + theme_bw()+ 
  theme(axis.text.x = element_text(angle=90))
}
save.image('Quant_Null_Genes_Plot.RData')