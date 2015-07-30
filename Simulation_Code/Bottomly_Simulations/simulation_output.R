### Load Bioconductor Packages
require(DESeq2)
require(edgeR)
require(limma)

### Load CRAN Packages
require(QuasiSeq)
require(samr)
require(fdrtool)
require(SimSeq)

### Set seed  
set.seed(3258986)

### Set simulation variables
n.iter <- 200     # Number of iterations
k.ind <- 3    # Sample size in each simulated treatment group 
n.genes <- 8500  # No of genes in each simulated matrix
n.diff <- 2000   # No of DE genes in each simulated matrix
n.genes.trim <- 5000 # No of genes to trim down to
n.diff.trim <- 1000 # No of DE genes to trim down to
filter.mean <- 10 # lower bound of average read count for simulated genes
filter.nonzero <- 2 # lower bound for nonzero read counts for simulated genes

### Load Data
counts <- read.table("bottomly_count_table.txt", header = TRUE, as.is = TRUE)
counts <- as.matrix(counts[, -1])
metadata <- read.table("bottomly_phenodata.txt", header = TRUE, as.is = TRUE)
treatment <- as.factor(metadata$strain)

### Remove low count genes
keep.counts <- ( rowMeans(counts) >= filter.mean ) & ( rowSums(counts > 0) >= filter.nonzero )
counts <- counts[keep.counts, ]

### Preprocessing Steps to speed up SimData function from SimSeq package

### Compute normalization factors to use in SimData function
### Effective library size is product of library size and size factors
### from calcNormFactors
lib.sizes <- apply(counts, 2, sum)
nf <- calcNormFactors(counts) * lib.sizes

### Compute weights to sample DE genes in SimData function
probs <- CalcPvalWilcox(counts, treatment = treatment,
                        sort.method = "unpaired", sorted = TRUE, norm.factors = nf,
                        exact = FALSE)
wghts <- 1 - fdrtool(probs, statistic = "pvalue", plot = FALSE, verbose = FALSE)$lfdr

### Initialize matrix of p-value output for each statistical method
pvals.samseq.simseq <- matrix(NA, nrow = n.iter, ncol = n.genes.trim)
pvals.quasiseq.simseq <- matrix(NA, nrow = n.iter, ncol = n.genes.trim)
pvals.edger.simseq <- matrix(NA, nrow = n.iter, ncol = n.genes.trim)
pvals.deseq2.simseq <- matrix(NA, nrow = n.iter, ncol = n.genes.trim)
pvals.voom.simseq <- matrix(NA, nrow = n.iter, ncol = n.genes.trim)

pvals.samseq.nb <- matrix(NA, nrow = n.iter, ncol = n.genes.trim)
pvals.quasiseq.nb <- matrix(NA, nrow = n.iter, ncol = n.genes.trim)
pvals.edger.nb <- matrix(NA, nrow = n.iter, ncol = n.genes.trim)
pvals.deseq2.nb <- matrix(NA, nrow = n.iter, ncol = n.genes.trim)
pvals.voom.nb <- matrix(NA, nrow = n.iter, ncol = n.genes.trim)

### if sample size in simulated matrices is greater than or equal to 
### 7, Cook's distance filtering replaces outlier counts which creates a new
### dataset so that another set of p-values is needed.
if(k.ind >= 7)
{
  pvals.samseq.cooks.simseq <- vector("list", n.iter)
  pvals.quasiseq.cooks.simseq <- vector("list", n.iter)
  pvals.edger.cooks.simseq <- vector("list", n.iter)
  pvals.deseq2.cooks.simseq <- vector("list", n.iter)
  pvals.voom.cooks.simseq <- vector("list", n.iter)
  
  pvals.samseq.cooks.nb <- vector("list", n.iter)
  pvals.quasiseq.cooks.nb <- vector("list", n.iter)
  pvals.edger.cooks.nb <- vector("list", n.iter)
  pvals.deseq2.cooks.nb <- vector("list", n.iter)
  pvals.voom.cooks.nb <- vector("list", n.iter)
}


### Preprocessing steps for simulating NB data
lambdas <- matrix(NA, nrow = nrow(counts), ncol = 2)

sum.nf.C57BL <- sum(nf[treatment == "C57BL/6J"])
sum.nf.DBA <- sum(nf[treatment == "DBA/2J"])
lambdas[, 1] <- rowSums(counts[, treatment == "C57BL/6J"])/sum.nf.C57BL
lambdas[, 2] <- rowSums(counts[, treatment == "DBA/2J"])/sum.nf.DBA

### rnbinom simulates counts according to Var(Y) = mu + mu^2/size
### dispersions from edgeR given from model with Var(Y) = mu + mu^2 * size
nb.disp.DBA <- 1 / readRDS("nbdisp_DBA_edgeR.RDS")
nb.disp.C57BL <- 1 / readRDS("nbdisp_C57BL_edgeR.RDS")

### matrix where each row gives a logical vector for the Cook's distance filter
filt.cooks.simseq <- matrix(NA, nrow = n.iter, ncol = n.genes.trim)
filt.cooks.nb <- matrix(NA, nrow = n.iter, ncol = n.genes.trim)

for(i in 1:n.iter){
  
  repeat{
    ### Simulate matrix of read counts from SimSeq
    ### On very rare occaisions, QuasiSeq errors out due to computational issues
    ### Repeat simulation when this occurs
  
    ### Simulate matrix of read counts from SimSeq
    counts.simseq.list <- SimData(counts = counts, treatment = treatment, 
                         sort.method = "unpaired", k.ind = k.ind, n.genes = n.genes,
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
    
    ### Use normalization factors from DBA/2J group to match SimSeq algorithm
    mu.samp[de.genes.sim, 1:k.ind] <- lambdas[de.genes, 1, drop = FALSE] %*% nf.samp[ (k.ind + 1):(2 * k.ind) ]
    mu.samp[de.genes.sim, (k.ind + 1):(2 * k.ind)] <- lambdas[de.genes, 2, drop = FALSE] %*% nf.samp[ (2 * k.ind + 1):(3 * k.ind) ]
    mu.samp[ !de.genes.sim, ] <-  lambdas[ee.genes, 2, drop = FALSE] %*% nf.samp[ (k.ind + 1):(3 * k.ind) ]
  
    ### Set dispersion estimates
    disp.samp <- matrix(NA, nrow = n.genes, ncol = 2)
    disp.samp[!de.genes.sim, 1] <- disp.samp[!de.genes.sim, 2] <- nb.disp.DBA[ee.genes]
    disp.samp[de.genes.sim, 1] <- nb.disp.C57BL[de.genes]
    disp.samp[de.genes.sim, 2] <- nb.disp.DBA[de.genes]
    
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
        
    ### DESeq Analysis
    colData <- data.frame(trt = c(rep(0, k.ind), rep(1, k.ind)))
    
    ### fit for SimSeq data
    des <- DESeqDataSetFromMatrix(countData = counts.simseq, colData = colData, design = formula(~trt))
    des <- DESeq(des)
    des.res <- results(des, cooksCutoff = FALSE, independentFiltering = FALSE)
    pvals.deseq2.simseq[i, ] <- des.res$pvalue
    if(k.ind < 7){
    filt.cooks.simseq[i, ] <- !(is.na(results(des, cooksCutoff = TRUE, independentFiltering = FALSE)$pvalue))
    }
    
    ### fit for NB data
    des <- DESeqDataSetFromMatrix(countData = counts.nb, colData = colData, design = formula(~trt))
    des <- DESeq(des)
    des.res <- results(des, cooksCutoff = FALSE, independentFiltering = FALSE)
    pvals.deseq2.nb[i, ] <- des.res$pvalue
    if(k.ind < 7){
    filt.cooks.nb[i, ] <- !(is.na(results(des, cooksCutoff = TRUE, independentFiltering = FALSE)$pvalue))
    }
    
    if(k.ind >= 7) #reccomended if ss >= 7
    {
      ### Replace outliers with trimmed mean via Cook's filtering rule for SimSeq
      des.cooks <- DESeqDataSetFromMatrix(countData = counts.simseq, colData = colData, design = formula(~trt))
      des.cooks <- DESeq(des.cooks) 
      des.cooks <- replaceOutliersWithTrimmedMean(des.cooks)
      counts.cooks.simseq <- counts(des.cooks)
      filt.cooks.simseq[i, ] <- ( rowMeans(counts.cooks.simseq) >= filter.mean ) & ( rowSums(counts.cooks.simseq > 0) >= filter.nonzero )
      counts.cooks.simseq <- counts.cooks.simseq[filt.cooks.simseq[i, ], ]
      
      ### Replace outliers with trimmed mean via Cook's filtering rule for NB
      des.cooks <- DESeqDataSetFromMatrix(countData = counts.nb, colData = colData, design = formula(~trt))
      des.cooks <- DESeq(des.cooks) 
      des.cooks <- replaceOutliersWithTrimmedMean(des.cooks)
      counts.cooks.nb <- counts(des.cooks)
      filt.cooks.nb[i, ] <- ( rowMeans(counts.cooks.nb) >= filter.mean ) & ( rowSums(counts.cooks.nb > 0) >= filter.nonzero )
      counts.cooks.nb <- counts.cooks.nb[filt.cooks.nb[i, ], ]
      
      ### fit for SimSeq data
      des <- DESeqDataSetFromMatrix(countData = counts.cooks.simseq, colData = colData, design = formula(~trt))
      des <- DESeq(des)
      des.res <- results(des, cooksCutoff = FALSE, independentFiltering = FALSE)
      pvals.deseq2.cooks.simseq[[i]] <- des.res$pvalue
  
      ### fit for NB data
      des <- DESeqDataSetFromMatrix(countData = counts.cooks.nb, colData = colData, design = formula(~trt))
      des <- DESeq(des)
      des.res <- results(des, cooksCutoff = FALSE, independentFiltering = FALSE)
      pvals.deseq2.cooks.nb[[i]] <- des.res$pvalue
   }
  
    ### QuasiSeq Analysis
    ### Create Design Matrices
    design.list <- vector("list", 2)
    trt <- c(rep(0, k.ind), rep(1, k.ind))
    design.list[[1]] <- model.matrix(~trt)  
    design.list[[2]] <- rep(1, length(trt))
    
    ### Compute Normalization Factors
    nf.simseq <- calcNormFactors(counts.simseq) * apply(counts.simseq, 2, sum)
    nf.nb <- calcNormFactors(counts.nb) * apply(counts.nb, 2, sum)
    
    ### fit for SimSeq data
    fit <- tryCatch(QL.fit(counts.simseq, design.list = design.list, log.offset = log(nf.simseq),
                          Model = "NegBin", method = "optim", print.progress = FALSE), error = function(w) NA)
    if( !is.list(fit) ) next
   
    res.fit <- QL.results(fit, Plot = FALSE)
    pvals.quasiseq.simseq[i, ] <- res.fit$P.values[[3]]
    
    ### fit for NB data
    fit <- tryCatch(QL.fit(counts.nb, design.list = design.list, log.offset = log(nf.nb),
                          Model = "NegBin", method = "optim", print.progress = FALSE), error = function(w) NA)
    if( !is.list(fit) ) next
    
    res.fit <- QL.results(fit, Plot = FALSE)
    pvals.quasiseq.nb[i, ] <- res.fit$P.values[[3]]
    
    
    if(k.ind >= 7){
      ### Compute Normalization Factors
      nf.cooks.simseq <- calcNormFactors(counts.cooks.simseq) * apply(counts.cooks.simseq, 2, sum)
      nf.cooks.nb <- calcNormFactors(counts.cooks.nb) * apply(counts.cooks.nb, 2, sum)
      
      ### fit for SimSeq data
      fit <- tryCatch(QL.fit(counts.cooks.simseq, design.list = design.list, log.offset = log(nf.simseq),
                             Model = "NegBin", method = "optim", print.progress = FALSE), error = function(w) NA)
      if( !is.list(fit) ) next
      
      res.fit <- QL.results(fit, Plot = FALSE)
      pvals.quasiseq.cooks.simseq[[i]] <- res.fit$P.values[[3]]
      
      ### fit for NB data
      fit <- tryCatch(QL.fit(counts.cooks.nb, design.list = design.list, log.offset = log(nf.nb),
                             Model = "NegBin", method = "optim", print.progress = FALSE), error = function(w) NA)
      if( !is.list(fit) ) next
      
      res.fit <- QL.results(fit, Plot = FALSE)
      pvals.quasiseq.cooks.nb[[i]] <- res.fit$P.values[[3]]
    }
 
    break
  }
 
 
  ### edgeR Analysis
  
  ### fit for SimSeq data
  y <- DGEList(counts = counts.simseq, group = trt)
  y <- calcNormFactors(y)
  design <- model.matrix(~trt)
  y <- estimateDisp(y, design)
  fit.edgeR <- glmFit(y, design)
  lrt <- glmLRT(fit.edgeR, coef = 2)
  pvals.edger.simseq[i, ] <- lrt$table$PValue
  
  ### fit for NB data
  y <- DGEList(counts = counts.nb, group = trt)
  y <- calcNormFactors(y)
  design <- model.matrix(~trt)
  y <- estimateDisp(y, design)
  fit.edgeR <- glmFit(y, design)
  lrt <- glmLRT(fit.edgeR, coef = 2)
  pvals.edger.nb[i, ] <- lrt$table$PValue

  if(k.ind >= 7){
    ### fit for SimSeq data
    y <- DGEList(counts = counts.cooks.simseq, group = trt)
    y <- calcNormFactors(y)
    design <- model.matrix(~trt)
    y <- estimateDisp(y, design)
    fit.edgeR <- glmFit(y, design)
    lrt <- glmLRT(fit.edgeR, coef = 2)
    pvals.edger.cooks.simseq[[i]] <- lrt$table$PValue
    
    ### fit for NB data
    y <- DGEList(counts = counts.cooks.nb, group = trt)
    y <- calcNormFactors(y)
    design <- model.matrix(~trt)
    y <- estimateDisp(y, design)
    fit.edgeR <- glmFit(y, design)
    lrt <- glmLRT(fit.edgeR, coef = 2)
    pvals.edger.cooks.nb[[i]] <- lrt$table$PValue
  }
  
 ### Limma Voom Analysis
 
 ### fit for SimSeq data
 y <- DGEList(counts = counts.simseq, group = trt)
 y <- calcNormFactors(y)
 design <- model.matrix(~trt)
 v <- voom(y, design)
 fit <- lmFit(v,design)
 fit <- eBayes(fit, proportion = 0.8)
 pvals.voom.simseq[i, ] <- fit$p.value[, 2]
 
 ### fit for NB data
 y <- DGEList(counts = counts.nb, group = trt)
 y <- calcNormFactors(y)
 design <- model.matrix(~trt)
 v <- voom(y, design)
 fit <- lmFit(v,design)
 fit <- eBayes(fit, proportion = 0.8)
 pvals.voom.nb[i, ] <- fit$p.value[, 2]
 
 if(k.ind >= 7){
   ### fit for SimSeq data
   y <- DGEList(counts = counts.cooks.simseq, group = trt)
   y <- calcNormFactors(y)
   design <- model.matrix(~trt)
   v <- voom(y, design)
   fit <- lmFit(v,design)
   fit <- eBayes(fit, proportion = 0.8)
   pvals.voom.cooks.simseq[[i]] <- fit$p.value[, 2]
   
   ### fit for NB data
   y <- DGEList(counts = counts.cooks.nb, group = trt)
   y <- calcNormFactors(y)
   design <- model.matrix(~trt)
   v <- voom(y, design)
   fit <- lmFit(v,design)
   fit <- eBayes(fit, proportion = 0.8)
   pvals.voom.cooks.nb[[i]] <- fit$p.value[, 2]
 } 
 
  #SAMSeq Analysis

  ### fit for SimSeq data
  trt <- c(rep(1, k.ind), rep(2, k.ind))
  samfit <- SAMseq(counts.simseq, trt, resp.type = "Two class unpaired")
  pvals.samseq.simseq[i, ] <- samr.pvalues.from.perms(samfit$samr.obj$tt, samfit$samr.obj$ttstar)

  ### fit for NB data
  trt <- c(rep(1, k.ind), rep(2, k.ind))
  samfit <- SAMseq(counts.nb, trt, resp.type = "Two class unpaired")
  pvals.samseq.nb[i, ] <- samr.pvalues.from.perms(samfit$samr.obj$tt, samfit$samr.obj$ttstar)
  
  if(k.ind >= 7){
    ### fit for SimSeq data
    trt <- c(rep(1, k.ind), rep(2, k.ind))
    samfit <- SAMseq(counts.cooks.simseq, trt, resp.type = "Two class unpaired")
    pvals.samseq.cooks.simseq[[i]] <- samr.pvalues.from.perms(samfit$samr.obj$tt, samfit$samr.obj$ttstar)
    
    ### fit for NB data
    trt <- c(rep(1, k.ind), rep(2, k.ind))
    samfit <- SAMseq(counts.cooks.nb, trt, resp.type = "Two class unpaired")
    pvals.samseq.cooks.nb[[i]] <- samr.pvalues.from.perms(samfit$samr.obj$tt, samfit$samr.obj$ttstar)
  }

}

mainDir <- getwd()
subDir <- "pval_output"
if( !file.exists(file.path(mainDir, subDir)) )
{
  dir.create(file.path(mainDir, subDir)) 
}
if( !file.exists(file.path(mainDir, subDir, paste0("ss", k.ind))) )
{
  dir.create(file.path(mainDir, subDir, paste0("ss", k.ind))) 
}
saveRDS(pvals.deseq2.simseq, file.path(mainDir, subDir, paste0("ss", k.ind), "p_deseq2_simseq.RDS"))
saveRDS(pvals.quasiseq.simseq, file.path(mainDir, subDir, paste0("ss", k.ind), "p_quasiseq_simseq.RDS"))
saveRDS(pvals.edger.simseq, file.path(mainDir, subDir, paste0("ss", k.ind), "p_edger_simseq.RDS"))
saveRDS(pvals.samseq.simseq, file.path(mainDir, subDir, paste0("ss", k.ind), "p_samseq_simseq.RDS"))
saveRDS(pvals.voom.simseq, file.path(mainDir, subDir, paste0("ss", k.ind), "p_voom_simseq.RDS"))

saveRDS(pvals.deseq2.nb, file.path(mainDir, subDir, paste0("ss", k.ind), "p_deseq2_nb.RDS"))
saveRDS(pvals.quasiseq.nb, file.path(mainDir, subDir, paste0("ss", k.ind), "p_quasiseq_nb.RDS"))
saveRDS(pvals.edger.nb, file.path(mainDir, subDir, paste0("ss", k.ind), "p_edger_nb.RDS"))
saveRDS(pvals.samseq.nb, file.path(mainDir, subDir, paste0("ss", k.ind), "p_samseq_nb.RDS"))
saveRDS(pvals.voom.nb, file.path(mainDir, subDir, paste0("ss", k.ind), "p_voom_nb.RDS"))

if(k.ind >= 7)
{
  saveRDS(pvals.deseq2.cooks.simseq, file.path(mainDir, subDir, paste0("ss", k.ind), "p_deseq2_cooks_simseq.RDS"))
  saveRDS(pvals.quasiseq.cooks.simseq, file.path(mainDir, subDir, paste0("ss", k.ind), "p_quasiseq_cooks_simseq.RDS"))
  saveRDS(pvals.edger.cooks.simseq, file.path(mainDir, subDir, paste0("ss", k.ind), "p_edger_cooks_simseq.RDS"))
  saveRDS(pvals.samseq.cooks.simseq, file.path(mainDir, subDir, paste0("ss", k.ind), "p_samseq_cooks_simseq.RDS"))
  saveRDS(pvals.voom.cooks.simseq, file.path(mainDir, subDir, paste0("ss", k.ind), "p_voom_cooks_simseq.RDS"))
  
  saveRDS(pvals.deseq2.cooks.nb, file.path(mainDir, subDir, paste0("ss", k.ind), "p_deseq2_cooks_nb.RDS"))
  saveRDS(pvals.quasiseq.cooks.nb, file.path(mainDir, subDir, paste0("ss", k.ind), "p_quasiseq_cooks_nb.RDS"))
  saveRDS(pvals.edger.cooks.nb, file.path(mainDir, subDir, paste0("ss", k.ind), "p_edger_cooks_nb.RDS"))
  saveRDS(pvals.samseq.cooks.nb, file.path(mainDir, subDir, paste0("ss", k.ind), "p_samseq_cooks_nb.RDS"))
  saveRDS(pvals.voom.cooks.nb, file.path(mainDir, subDir, paste0("ss", k.ind), "p_voom_cooks_nb.RDS"))
}

saveRDS(filt.cooks.simseq, file.path(mainDir, subDir, paste0("ss", k.ind), "filt_cooks_simseq.RDS"))
saveRDS(filt.cooks.nb, file.path(mainDir, subDir, paste0("ss", k.ind), "filt_cooks_nb.RDS"))
