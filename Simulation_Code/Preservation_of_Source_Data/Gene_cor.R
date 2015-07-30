### Load CRAN packages
require(SimSeq)
require(fdrtool)
require(matrixStats)
require(ggplot2)
require(gridExtra)

### Load Bioconductor packages
require(edgeR)

### Set seed  
set.seed(392313412)

### Set simulation variables
n.iter <- 200  # Number of iterations
k.ind <- 10    # Sample size in each simulated treatment group 
n.genes <- 8500  # No of genes in each simulated matrix
n.diff <- 2000   # No of DE genes in each simulated matrix
n.genes.trim <- 5000 # No of genes to trim down to
n.diff.trim <- 1000 # No of DE genes to trim down to
filter.mean <- 10 # lower bound of average read count for simulated genes
filter.nonzero <- 2 # lower bound for nonzero read counts for simulated genes

### Load Data
data(kidney)
counts <- kidney$counts
tumor <- kidney$treatment
replic <- kidney$replic

### Remove low count genes
keep.counts <- ( rowMeans(counts) >= filter.mean ) & ( rowSums(counts > 0) >= filter.nonzero )
counts <- counts[keep.counts, ]

### Preprocessing Steps to speed up SimData function from SimSeq package

### Compute normalization factors to use in SimData function
### Effective library size is product of library size and size factors
### from calcNormFactors
lib.sizes <- apply(counts, 2, sum)
nf <- calcNormFactors(counts) * lib.sizes
nf <- nf/exp(mean(log(nf))) #normalize so that geometric mean is 1

### Compute weights to sample DE genes in SimData function
probs <- CalcPvalWilcox(counts, treatment = tumor, replic = replic,
                        sort.method = "paired", sorted = TRUE, nf,
                        exact = FALSE)
wghts <- 1 - fdrtool(probs, statistic = "pvalue", plot = FALSE, verbose = FALSE)$lfdr

### Compute test statistics on source dataset

counts.norm <- t(t(counts[, tumor == "Tumor"]) * 1/nf[tumor == "Tumor"])
cor.source <- cor(t(counts.norm), method = "spearman")

### Preprocessing steps for simulating NB data
lambdas <- matrix(NA, nrow = nrow(counts), ncol = 2)

sum.nf.nontumor <- sum(nf[tumor == "Non-Tumor"])
sum.nf.tumor <- sum(nf[tumor == "Tumor"])
lambdas[, 1] <- rowSums(counts[, tumor == "Non-Tumor"])/sum.nf.nontumor
lambdas[, 2] <- rowSums(counts[, tumor == "Tumor"])/sum.nf.tumor

### rnbinom simulates counts according to Var(Y) = mu + mu^2/size
### dispersions from edgeR given from model with Var(Y) = mu + mu^2 * size
nb.disp.tumor <- 1 / readRDS("nbdisp_tumor_edgeR.RDS")
nb.disp.nontumor <- 1 / readRDS("nbdisp_nontumor_edgeR.RDS")


### This code was used to find a pair of genes involved in epithelial cell differentiation  with
### high Spearman Rank correlation in tumor group
### does not need to be run

# GO.Hs <- readRDS("Gene_Ontologies_Hs.RDS")
# kp.samp <- sapply(GO.Hs, function(x) "GO:0030855" %in% x)[keep.counts] ### GO:0030855 epithelial cell differentiation
# cor.source[kp.samp, kp.samp] > 0.8
# 
# cor.source[kp.krt, kp.krt][16, 55] #0.8167406 correlation - find large sample correlation
# which(kp.krt)[c(16, 55)] #row numbers for pair of genes are 3393, 14299

### storage for spearman rank correlation between genes in row 3393 and 14299
cor.simseq <- rep(NA, n.iter)
cor.nb <- rep(NA, n.iter)

for(i in 1:n.iter){
  print(i)
  
  ### repeat simulations until gene pair 3393 and 14299 are included
  repeat{
  ### Simulate matrix of read counts from SimSeq
  counts.simseq.list <- SimData(counts = counts, treatment = tumor, replic = replic, 
                                sort.method = "paired", k.ind = k.ind, n.genes = n.genes,
                                n.diff = n.diff, norm.factors = nf, weights = wghts,
                                switch.trt = TRUE)
  
  counts.simseq <- counts.simseq.list$counts # Simulated Count matrix from SimSeq
  genes.samp <- counts.simseq.list$genes.subset # Genes sampled from source matrix
  de.genes <- counts.simseq.list$DE.genes # DE genes sampled from source matrix
  ee.genes <- genes.samp[ ! (genes.samp %in% de.genes) ] # EE genes sampled from source matrix
  samp.col <- counts.simseq.list$col # Columns sampled in SimSeq algorithm
  de.genes.sim <- counts.simseq.list$genes.subset %in% de.genes # logical vector giving which genes are DE in simulted matrix
  
  ### Compute matrix of means for simulating from NB model
  mu.samp <- matrix(NA, nrow = n.genes, ncol = 2 * k.ind)
  nf.samp <- nf[samp.col]
  
  ### Use normalization factors from Tumor group to match SimSeq algorithm
  mu.samp[de.genes.sim, 1:k.ind] <- lambdas[de.genes, 1, drop = FALSE] %*% nf.samp[ (k.ind + 1):(2 * k.ind) ]
  mu.samp[de.genes.sim, (k.ind + 1):(2 * k.ind)] <- lambdas[de.genes, 2, drop = FALSE] %*% nf.samp[ (2 * k.ind + 1):(3 * k.ind) ]
  mu.samp[ !de.genes.sim, ] <-  lambdas[ee.genes, 2, drop = FALSE] %*% nf.samp[ (k.ind + 1):(3 * k.ind) ]
  
  ### Set dispersion estimates
  disp.samp <- matrix(NA, nrow = n.genes, ncol = 2)
  disp.samp[!de.genes.sim, 1] <- disp.samp[!de.genes.sim, 2] <- nb.disp.tumor[ee.genes]
  disp.samp[de.genes.sim, 1] <- nb.disp.nontumor[de.genes]
  disp.samp[de.genes.sim, 2] <- nb.disp.tumor[de.genes]
  
  ### Simulate matrix of read counts from NB model
  counts.nb <- matrix(NA, nrow = n.genes, ncol = 2 * k.ind)
  for(jj in 1:n.genes){
    counts.nb[jj, 1:k.ind] <- rnbinom(k.ind, size = disp.samp[jj, 1], mu = mu.samp[jj, 1:k.ind])
    counts.nb[jj, -(1:k.ind)] <- rnbinom(k.ind, size = disp.samp[jj, 2], mu = mu.samp[jj, -(1:k.ind)])
  }
  
  ### Apply filtering rules to both simulated datasets and only keep genes who pass both filters
  keep.genes.simseq <- ( rowMeans(counts.simseq) >= filter.mean ) & ( rowSums(counts.simseq > 0) >= filter.nonzero )
  keep.genes.nb <- ( rowMeans(counts.nb) >= filter.mean ) & ( rowSums(counts.nb > 0) >= filter.nonzero )
  keep <- keep.genes.simseq & keep.genes.nb
  
  ee.genes <- sample(which(!de.genes.sim & keep), n.genes.trim - n.diff.trim)
  de.genes <- sample(which(de.genes.sim & keep), n.diff.trim)
  counts.simseq <- counts.simseq[c(ee.genes, de.genes), ]
  counts.nb <- counts.nb[c(ee.genes, de.genes), ]
  
  
    if( all(c(3393, 14299) %in% genes.samp[c(ee.genes, de.genes)]) ){
      break
    }
  }
  
  ### Compute norm factors, scale to geometric mean of 1
  nf.simseq <- calcNormFactors(counts.simseq) * apply(counts.simseq, 2, sum)
  nf.nb <- calcNormFactors(counts.nb) * apply(counts.nb, 2, sum)
  nf.simseq <- nf.simseq/exp(mean(log(nf.simseq)))
  nf.nb <- nf.nb/exp(mean(log(nf.nb)))
  
  counts.simseq.norm <- t(t(counts.simseq[, -(1:k.ind)]) * 1/nf.simseq[-(1:k.ind)])
  counts.nb.norm <- t(t(counts.nb[, -(1:k.ind)]) * 1/nf.nb[-(1:k.ind)])
  
  kp <- genes.samp[c(ee.genes, de.genes)] %in% c(3393, 14299)
  cor.simseq[i] <- cor(t(counts.simseq.norm[kp,]), method = "spearman")[2,1]
  cor.nb[i] <- cor(t(counts.nb.norm[kp,]), method = "spearman")[2,1]
  
}

### Make plots
p <- ggplot() + geom_histogram(aes(x = cor.simseq), colour = "black", fill = "grey", binwidth = 0.1) + 
  xlab("Spearman Rank Correlation") + theme_bw() + ggtitle("A") + xlim(-1,1)
q <- ggplot() + geom_histogram(aes(x = cor.nb), colour = "black", fill = "grey", binwidth = 0.1) + 
  xlab("Spearman Rank Correlation") + theme_bw() + ggtitle("B") + xlim(-1,1)
grid.arrange(p, q, ncol = 2)
