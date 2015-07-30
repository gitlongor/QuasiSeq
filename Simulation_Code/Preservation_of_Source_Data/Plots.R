### Load Bioconductor packages
require(edgeR)
require(org.Hs.eg.db)
require(annotate)

### Load CRAN packages
require(SimSeq)
require(fdrtool)
require(stringr)
require(ggplot2)
require(gridExtra)


### Simulate one SimSeq and negative binomial dataset

### Set seed  
set.seed(783623211)

### Set simulation variables
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

counts.simseq.list <- SimData(counts = counts, treatment = tumor, replic = replic, 
                              sort.method = "paired", k.ind = k.ind, n.genes = n.genes,
                              n.diff = n.diff, norm.factors = nf, weights = wghts, switch.trt = TRUE)

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

nf.simseq <- calcNormFactors(counts.simseq) * apply(counts.simseq, 2, sum)
nf.simseq <- nf.simseq/exp(mean(log(nf.simseq))) #normalize so that geometric mean is 1
nf.nb <- calcNormFactors(counts.nb) * apply(counts.nb, 2, sum)
nf.nb <- nf.nb/exp(mean(log(nf.nb))) #normalize so that geometric mean is 1

y.source <- DGEList(counts = counts[genes.samp[c(ee.genes, de.genes)], ], group = tumor)
y.source <- calcNormFactors(y.source)
design <- model.matrix(~tumor)
y.source <- estimateDisp(y.source, design)

y.nb <- DGEList(counts = counts.nb, group = c(rep("Trt 1",k.ind), rep("Trt 2", k.ind)))
y.nb <- calcNormFactors(y.nb)
design <- model.matrix(~c(rep("Trt 1",k.ind), rep("Trt 2", k.ind)))
y.nb <- estimateDisp(y.nb)

y.simseq <- DGEList(counts = counts.simseq, group = c(rep("Trt 1",k.ind), rep("Trt 2", k.ind)))
y.simseq <- calcNormFactors(y.simseq)
design <- model.matrix(~c(rep("Trt 1",k.ind), rep("Trt 2", k.ind)))
y.simseq <- estimateDisp(y.simseq)

### Make Mean Variance Plot
### Old
#par(mfrow = c(1,3))
#plotMeanVar(y.source, main = "A", show.raw.vars = TRUE, show.ave.raw.vars = FALSE, pch = 19)
#plotMeanVar(y.simseq, main = "B", show.raw.vars = TRUE, show.ave.raw.vars = FALSE, pch = 19)
#plotMeanVar(y.nb, main = "C", show.raw.vars = TRUE, show.ave.raw.vars = FALSE, pch = 19)
#par(mfrow = c(1,1))

### Make Mean Variance Plot
par(mfrow = c(1,3))
counts.source <- counts[genes.samp[c(ee.genes, de.genes)], ]
x <- log2(c(rowMeans(counts.source[, tumor == "Non-Tumor"]), 
            rowMeans(counts.source[, tumor == "Tumor"])))
y <- log2(c(rowVars(counts.source[, tumor == "Non-Tumor"]), 
            rowVars(counts.source[, tumor == "Tumor"])))
smoothScatter(x, y, xlab = "Log Base 2 of Mean Expression", ylab = "Log Base 2 of Sample Variance",
              xlim = c(-1,18), ylim = c(0, 35), main = "A")

x <- log2(c(rowMeans(counts.simseq[, 1:k.ind]), 
            rowMeans(counts.simseq[, -(1:k.ind)])))
y <- log2(c(rowVars(counts.simseq[, 1:k.ind]), 
            rowVars(counts.simseq[, -(1:k.ind)])))
smoothScatter(x, y, xlab = "Log Base 2 of Mean Expression", ylab = "Log Base 2 of Sample Variance",
              xlim = c(-1,18), ylim = c(0, 35), main = "B")

x <- log2(c(rowMeans(counts.nb[, 1:k.ind]), 
            rowMeans(counts.nb[, -(1:k.ind)])))
y <- log2(c(rowVars(counts.nb[, 1:k.ind]), 
            rowVars(counts.nb[, -(1:k.ind)])))
smoothScatter(x, y, xlab = "Log Base 2 of Mean Expression", ylab = "Log Base 2 of Sample Variance",
              xlim = c(-1,18), ylim = c(0, 35), main = "C")
par(mfrow = c(1,1))


### Subset down to just DE genes for MA plot
y.source$counts <- y.source$counts[-(1:(n.genes.trim - n.diff.trim)), ]
y.nb$counts <- y.nb$counts[-(1:(n.genes.trim - n.diff.trim)), ]
y.simseq$counts <- y.simseq$counts[-(1:(n.genes.trim - n.diff.trim)), ]

### Make MA plot
par(mfrow = c(1,3))
plotSmear(y.source, ylim = c(-10, 8), xlim = c(-3,15), ylab = "Log Fold Change", main = "A")
plotSmear(y.nb, ylim = c(-10, 8), xlim = c(-3,15), ylab = "Log Fold Change", main = "B")
plotSmear(y.simseq, ylim = c(-10, 8), xlim = c(-3,15), ylab = "Log Fold Change", main = "C")
par(mfrow = c(1,1))

### Load GO's for each gene annotation
GO.Hs <- readRDS("Gene_Ontologies_Hs.RDS")
kp.genes.source <- sapply(GO.Hs, function(x) "GO:0030855" %in% x)[keep.counts]
#"GO:0030855" genes involved in epithelial cell differentiation
#gene pair 3393, 14299 used in Figure 3

kp.samp <- genes.samp[c(ee.genes, de.genes)]

### Normalize counts and compute Spearman's rank correlation in Tumor group
cor.source <- cor(t(counts[, tumor == "Tumor"]) * 1/nf[tumor == "Tumor"] , method = "spearman")
cor.simseq <- cor(t(counts.simseq[,-(1:k.ind)]) * 1/nf.simseq[-(1:k.ind)] , method = "spearman")
cor.nb <- cor(t(counts.nb[, -(1:k.ind)]) * 1/nf.nb[-(1:k.ind)] , method = "spearman")

GO.samp <- kp.samp %in% which(kp.genes.source)

### Shape Spearman's rank correlation from matrix to vector form
kp.fin <- as.vector(lower.tri((cor.source[kp.samp, kp.samp])[GO.samp, GO.samp])) #remove dup's on upper triangle
c.org <- as.vector((cor.source[kp.samp, kp.samp])[GO.samp, GO.samp])[kp.fin]
c.sim <- as.vector(cor.simseq[GO.samp, GO.samp])[kp.fin]
c.nb <- as.vector(cor.nb[GO.samp, GO.samp])[kp.fin]

# gene pair 3393, 14299
hglight <- c.org %in% cor.source[3393, 14299]
hglight <- ifelse(hglight, "red", "black")

### Make plot of Spearman Correlation for simulated data vs source data
p <- ggplot() + geom_point(aes(x = c.org, y = c.sim, colour = hglight)) + 
  scale_color_identity() + 
  xlab("Spearman Rank Correlation, KIRC") + ylab("Spearman Rank Correlation, SimSeq") +
  theme_bw() + ggtitle("A")
q <- ggplot() + geom_point(aes(x = c.org, y = c.nb, colour = hglight)) + 
  scale_color_identity() + 
  xlab("Spearman Rank Correlation, KIRC") + ylab("Spearman Rank Correlation, Negative Binomial") +
  theme_bw() + ggtitle("B")
grid.arrange(p, q, ncol = 2)