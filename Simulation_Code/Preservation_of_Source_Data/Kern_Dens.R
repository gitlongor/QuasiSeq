### Load in packages required for simulation/analysis

### CRAN packages
require(SimSeq)
require(fdrtool)
require(matrixStats)
require(ggplot2)
require(gridExtra)

### Bioconductor packages
require(edgeR)

### Set seed  
set.seed(392313412)

### Set simulation variables
n.iter <- 200   # Number of iterations
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

n.k <- 2^12
st <- 0
fn <- 50000
mean.dens.samp <- matrix(nrow = n.k, ncol = n.iter)
var.dens.samp <- matrix(nrow = n.k, ncol = n.iter)
lfc.dens.samp <- matrix(nrow = n.k, ncol = n.iter)
mean.dens.grid <- density(log2(rowMeans(counts[, tumor == "Tumor"])),
                     n = n.k, from = st, to = log2(fn))$x
mean.dens.source <- density(log2(rowMeans(counts[, tumor == "Tumor"])),
                            n = n.k, from = st, to = log2(fn))$y
var.dens.grid <- density(log2(apply(counts[, tumor == "Tumor"], 1, var)),
                          n = n.k, from = st, to = 30)$x
var.dens.source <- density(log2(apply(counts[, tumor == "Tumor"], 1, var)),
                            n = n.k, from = st, to = 30)$y

lfc.source <- rowMeans(log2(counts[, tumor == "Tumor"] + 1)) 
- rowMeans(log2(counts[, tumor == "Non-Tumor"] + 1)) 

lfc.dens.grid <- density(lfc.source, n = n.k, from = -3, to = 3)$x

### for log fold change dist of the source, we average the estimates
### using the sampled DE genes from each simulation
lfc.dens.source <- matrix(nrow = n.k, ncol = n.iter)

for(i in 1:n.iter){
  print(i)
  ### Simulate matrix of read counts from SimSeq
  counts.simseq.list <- SimData(counts = counts, treatment = tumor, replic = replic, 
                                sort.method = "paired", k.ind = k.ind, n.genes = n.genes,
                                n.diff = n.diff, norm.factors = nf, weights = wghts, switch.trt = TRUE)
  
  counts.simseq <- counts.simseq.list$counts # Simulated Count matrix from SimSeq
  genes.samp <- counts.simseq.list$genes.subset # Genes sampled from source matrix
  de.genes <- counts.simseq.list$DE.genes # DE genes sampled from source matrix
  ee.genes <- genes.samp[ ! (genes.samp %in% de.genes) ] # EE genes sampled from source matrix
  samp.col <- counts.simseq.list$col # Columns sampled in SimSeq algorithm
  de.genes.sim <- counts.simseq.list$genes.subset %in% de.genes # logical vector giving which genes are DE in simulted matrix
  
  ### Apply filtering rules to both simulated datasets and only keep genes who pass both filters
  keep <- ( rowMeans(counts.simseq) >= filter.mean ) & ( rowSums(counts.simseq > 0) >= filter.nonzero )
  
  ee.genes <- sample(which(!de.genes.sim & keep), n.genes.trim - n.diff.trim)
  de.genes <- sample(which(de.genes.sim & keep), n.diff.trim)
  counts.simseq <- counts.simseq[c(ee.genes, de.genes), ]
  
  row.means.samp <- log2(rowMeans(counts.simseq[, -(1:k.ind)]))
  row.var.samp <- log2(apply(counts.simseq[, -(1:k.ind)], 1, var))
  row.lfc.samp <- rowMeans(log2((counts.simseq[-(1:(n.genes.trim - n.diff.trim)), 1:k.ind] + 1))) - 
    rowMeans(log2((counts.simseq[-(1:(n.genes.trim - n.diff.trim)), -(1:k.ind)] + 1)))
  mean.dens.samp[, i] <- density(row.means.samp, n = n.k, from = st, to = log2(fn))$y
  var.dens.samp[, i] <- density(row.var.samp, n = n.k, from = st, to = 30)$y
  lfc.dens.samp[, i] <- density(row.lfc.samp, n = n.k, from = -3, to = 3)$y
  
  
  lfc.source <- rowMeans(log2(counts[genes.samp[de.genes], tumor == "Tumor"] + 1)) - 
    rowMeans(log2(counts[genes.samp[de.genes], tumor == "Non-Tumor"] + 1)) 

  lfc.dens.source[, i] <- density(lfc.source, n = n.k, from = -3, to = 3)$y
 }

p1 <- ggplot() + 
geom_line(aes(mean.dens.grid, mean.dens.source), colour = "red") +
geom_line(aes(mean.dens.grid, rowMeans(mean.dens.samp))) +
geom_line(aes(mean.dens.grid, apply(mean.dens.samp, 1, quantile, 0.025)), linetype = "dashed") +
geom_line(aes(mean.dens.grid, apply(mean.dens.samp, 1, quantile, 0.975)), linetype = "dashed") +
theme_bw() + xlab("Log Base 2 of Mean Expression") + ylab("Density") + ggtitle("A")

p2 <- ggplot() + 
  geom_line(aes(var.dens.grid, var.dens.source), colour = "red") +
  geom_line(aes(var.dens.grid, rowMeans(var.dens.samp))) +
  geom_line(aes(var.dens.grid, apply(var.dens.samp, 1, quantile, 0.025)), linetype = "dashed") +
  geom_line(aes(var.dens.grid, apply(var.dens.samp, 1, quantile, 0.975)), linetype = "dashed") +
  theme_bw() + xlab("Log Base 2 of Sample Variance") + ylab("Density") + ggtitle("B")

p3 <- ggplot() + 
  geom_line(aes(lfc.dens.grid, rowMeans(lfc.dens.source)), colour = "red") +
  geom_line(aes(lfc.dens.grid, rowMeans(lfc.dens.samp))) +
  geom_line(aes(lfc.dens.grid, apply(lfc.dens.samp, 1, quantile, 0.025)), linetype = "dashed") +
  geom_line(aes(lfc.dens.grid, apply(lfc.dens.samp, 1, quantile, 0.975)), linetype = "dashed") +
  theme_bw() + xlab("Log Base 2 Fold Change") + ylab("Density") + ggtitle("C")

grid.arrange(p1, p2, p3, ncol = 3)
