### Load Bioconductor Packages
require(edgeR)

### Load Data
counts <- read.table("bottomly_count_table.txt", header = TRUE, as.is = TRUE)
counts <- as.matrix(counts[, -1])
metadata <- read.table("bottomly_phenodata.txt", header = TRUE, as.is = TRUE)
metadata$strain <- as.factor(metadata$strain)

### Set simulation variables
filter.mean <- 10 # lower bound of average read count for simulated genes
filter.nonzero <- 2 # lower bound for nonzero read counts for simulated genes

### Remove low count genes
keep.counts <- ( rowMeans(counts) >= filter.mean ) & ( rowSums(counts > 0) >= filter.nonzero )
counts <- counts[keep.counts, ]

### Subset down to counts from C57BL/6J group for dispersion estimation in C57BL/6J group
counts.C57BL <- counts[ , metadata$strain == "C57BL/6J"]

### Use edgeR package's tagwise dispersion estimate
### Calculate tagwise estimates from edgeR. Not that edgeR uses the
### parameterizion: Var(y) = \mu + \omega * \mu^2 where E(Y) = \mu and
### \omega is a dispersion parameter for the NB model.
y <- DGEList(counts = counts.C57BL)
y <- calcNormFactors(y)
y <- estimateDisp(y)
nbdisp.C57BL <- y$tagwise.dispersion

### Repeat analysis for DBA/2J estimates
### Subset down to just DBA/2J counts for dispersion estimation in DBA/2J group
counts.DBA <- counts[ , metadata$strain == "DBA/2J"]

### Calculate tagwise estimates from edgeR. Not that edgeR uses the
### parameterizion: Var(y) = \mu + \omega * \mu^2 where E(Y) = \mu and
### \omega is a dispersion parameter for the NB model.
y <- DGEList(counts = counts.DBA)
y <- calcNormFactors(y)
y <- estimateDisp(y)
nbdisp.DBA <- y$tagwise.dispersion

### Save Results
saveRDS(nbdisp.C57BL, file = "nbdisp_C57BL_edgeR.RDS")
saveRDS(nbdisp.DBA, file = "nbdisp_DBA_edgeR.RDS")
