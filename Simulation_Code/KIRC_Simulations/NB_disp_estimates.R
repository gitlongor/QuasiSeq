### Load Bioconductor Packages
require(edgeR)

### Load CRAN Packages
require(SimSeq)

### Load Data
data(kidney)
counts <- kidney$counts
treatment <- kidney$treatment

### Set simulation variables
filter.mean <- 10 # lower bound of average read count for simulated genes
filter.nonzero <- 2 # lower bound for nonzero read counts for simulated genes

### Remove low count genes
keep.counts <- ( rowMeans(counts) >= filter.mean ) & ( rowSums(counts > 0) >= filter.nonzero )
counts <- counts[keep.counts, ]

### Subset down to just tumor counts for dispersion estimation in tumor group
counts.tumor <- counts[ , treatment == "Tumor"]

### Use edgeR package's tagwise dispersion estimate
### Calculate tagwise estimates from edgeR. Not that edgeR uses the
### parameterizion: Var(y) = \mu + \omega * \mu^2 where E(Y) = \mu and
### \omega is a dispersion parameter for the NB model.
y <- DGEList(counts = counts.tumor)
y <- calcNormFactors(y)
y <- estimateDisp(y)
nbdisp.tumor <- y$tagwise.dispersion

### Subset down to just non-tumor counts for dispersion estimation in non-tumor group
counts.nontumor <- counts[ , treatment == "Non-Tumor"]

### Calculate tagwise estimates from edgeR. Not that edgeR uses the
### parameterizion: Var(y) = \mu + \omega * \mu^2 where E(Y) = \mu and
### \omega is a dispersion parameter for the NB model.
y <- DGEList(counts = counts)
y <- calcNormFactors(y)
y <- estimateDisp(y)
nbdisp.nontumor <- y$tagwise.dispersion

### Save Results
saveRDS(nbdisp.tumor, file = "nbdisp_tumor_edgeR.RDS")
saveRDS(nbdisp.nontumor, file = "nbdisp_nontumor_edgeR.RDS")
