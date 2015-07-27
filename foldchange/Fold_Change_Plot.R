print(getwd())
setwd("C:/Users/Klirk/Desktop/Simulation_Code/Simulation_Code/KIRC_Simulations/foldchange")
require(ggplot2)
require(SimSeq)

### Set simulation variables
n.genes <- 8500  # No of genes in each simulated matrix
n.diff <- 2000   # No of DE genes in each simulated matrix
n.genes.trim <- 5000 # No of genes to trim down to
n.diff.trim <- 1000 # No of DE genes to trim down to
filter.mean <- 0 # lower bound of average read count for simulated genes
filter.nonzero <- 1 # lower bound for nonzero read counts for simulated genes
qval.cuts <- seq(0, 0.15, by = 0.001)
null.genes <- c( rep(TRUE, n.genes.trim - n.diff.trim), rep(FALSE, n.diff.trim))

### Load Data
bak.seed=.Random.seed
RNGkind("Mersenne-Twister", "Inversion")
set.seed(156165L)
data(kidney); 
kidney$counts=kidney$counts[sample(nrow(kidney$counts)),]
.Random.seed = bak.seed
counts <- kidney$counts
tumor <- kidney$treatment
replic <- kidney$replic

### Remove low count genes
keep.counts <- ( rowMeans(counts) >= filter.mean ) & ( rowSums(counts > 0) >= filter.nonzero )
counts <- counts[keep.counts, ]

lib.sizes <- apply(counts, 2, sum)
nf <- calcNormFactors(counts) * lib.sizes

lambdas <- matrix(NA, nrow = nrow(counts), ncol = 2)

sum.nf.nontumor <- sum(nf[tumor == "Non-Tumor"])
sum.nf.tumor <- sum(nf[tumor == "Tumor"])
lambdas[, 1] <- rowSums(counts[, tumor == "Non-Tumor"])/sum.nf.nontumor
lambdas[, 2] <- rowSums(counts[, tumor == "Tumor"])/sum.nf.tumor

logFoldChange = lambdas[, 2] - lambdas[, 1]

QuantThresh <- function(pvals, null.genes, specificity){
  quantile(pvals[null.genes], specificity)
}

k.ind <- 5

mainDir <- getwd()
subDir <- "pval_output"

fold.quasiseq.simseq.fit1 <- readRDS(file.path(mainDir, subDir, paste0("ss", k.ind), "fold_quasiseq_simseq_fit1.RDS"))
fold.quasiseq.simseq.fit2 <- readRDS(file.path(mainDir, subDir, paste0("ss", k.ind), "fold_quasiseq_simseq_fit2.RDS"))
fold.quasiseq.simseq.fit3 <- readRDS(file.path(mainDir, subDir, paste0("ss", k.ind), "fold_quasiseq_simseq_fit3.RDS"))
fold.voom.simseq <- readRDS(file.path(mainDir, subDir, paste0("ss", k.ind), "fold_voom_simseq.RDS"))

fold.quasiseq.nb.fit1 <- readRDS(file.path(mainDir, subDir, paste0("ss", k.ind), "fold_quasiseq_nb_fit1.RDS"))
fold.quasiseq.nb.fit2 <- readRDS(file.path(mainDir, subDir, paste0("ss", k.ind), "fold_quasiseq_nb_fit2.RDS"))
fold.quasiseq.nb.fit3 <- readRDS(file.path(mainDir, subDir, paste0("ss", k.ind), "fold_quasiseq_nb_fit3.RDS"))
fold.voom.nb <- readRDS(file.path(mainDir, subDir, paste0("ss", k.ind), "fold_voom_nb.RDS"))

fold.trulyDEgeneList <- readRDS(file.path(mainDir, subDir, paste0("ss", k.ind), "fold_trulyDEgeneList.RDS"))

#convert to log2
fold.quasiseq.simseq.fit1 <- log2(exp(fold.quasiseq.simseq.fit1*2))
fold.quasiseq.simseq.fit2 <- log2(exp(fold.quasiseq.simseq.fit2*2))
fold.quasiseq.simseq.fit3 <- log2(exp(fold.quasiseq.simseq.fit3*2))
fold.voom.simseq <- log2(exp(fold.voom.simseq*2))

fold.quasiseq.nb.fit1 <- log2(exp(fold.quasiseq.nb.fit1*2))
fold.quasiseq.nb.fit2 <- log2(exp(fold.quasiseq.nb.fit2*2))
fold.quasiseq.nb.fit3 <- log2(exp(fold.quasiseq.nb.fit3*2))
fold.voom.nb <- log2(exp(fold.voom.nb*2))


##end of work
num.quasiseq.simseq.fit1 <- apply(fold.quasiseq.simseq.fit1, 1, QuantThresh, null.genes, 0.05)
num.quasiseq.simseq.fit2 <- apply(fold.quasiseq.simseq.fit2, 1, QuantThresh, null.genes, 0.05)
num.quasiseq.simseq.fit3 <- apply(fold.quasiseq.simseq.fit3, 1, QuantThresh, null.genes, 0.05)
num.voom.simseq <- apply(fold.voom.simseq, 1, QuantThresh, null.genes, 0.05)

num.quasiseq.nb.fit1 <- apply(fold.quasiseq.nb.fit1, 1, QuantThresh, null.genes, 0.05)
num.quasiseq.nb.fit2 <- apply(fold.quasiseq.nb.fit2, 1, QuantThresh, null.genes, 0.05)
num.quasiseq.nb.fit3 <- apply(fold.quasiseq.nb.fit3, 1, QuantThresh, null.genes, 0.05)
num.voom.nb <- apply(fold.voom.nb, 1, QuantThresh, null.genes, 0.05)

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
  
save.image('Fold_Change_Plot.RData')