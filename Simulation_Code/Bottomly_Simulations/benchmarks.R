bh.q=function(p)
{
  ### The following function computes q-values given a set of p-values 
  ### using the approach of Benjamini-Hochberg (1995)
  
  m = length(p)
  k = 1:m
  ord = order(p)
  p[ord] = (p[ord] * m)/(1:m)
  qval = p
  qval[ord]=rev(cummin(rev(qval[ord])))
  return(qval)
}

pau <- function(pvals, null.genes, specificity)
{
  ### The following function computes partial area under the ROC curve
  ### over the interval from zero to the level of specificitiy given.
  ### pvals is a vector of p-values for testing the null hypothesis of 
  ### no differential expression for each gene, null.genes is a logical vector
  ### which genes are null genes, and specificity is the upper endpoint of the
  ### interval over which to calculate the area
  specificity <- quantile(pvals[null.genes], specificity)
  disc <- which(pvals <= specificity)
  disc <- disc[order(pvals[disc])]
  disc <- !null.genes[disc]
  tot.null <- sum(null.genes)
  tot.de <- sum(!null.genes)
  height <- cumsum(disc)[disc == 0]/tot.de
  pau <- sum(height/tot.null)
  return(pau)
}

pau.cooks <- function(pvals, filt, null.genes, specificity, pvals.list = FALSE)
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
    res[i] <- pau(pvals.temp, null.genes.temp, specificity)
  }
  return(res)
}


fdp <- function(pvals, null.genes, qval.cuts){
  qvals <- bh.q(pvals)  
  res <- sapply(qval.cuts, function(x) mean(null.genes[which(qvals < x)]) )
  res[is.na(res)] <- 0
  res <- res - qval.cuts
  return(res)
}

fdp.cooks <- function(pvals, filt, null.genes, qval.cuts, pvals.list = FALSE)
{
  if(pvals.list){
    n.iter <- length(pvals)
  } else {
    n.iter <- nrow(pvals)
  }
  res <- matrix(NA, nrow = length(qval.cuts), ncol = n.iter)
  for(i in 1:n.iter)
  {
    if(pvals.list){
      pvals.temp <- pvals[[i]]
    } else {
      pvals.temp <- pvals[i, ][filt[i, ]] 
    }
    null.genes.temp <- null.genes[filt[i, ]]
    res[, i] <- fdp(pvals.temp, null.genes.temp, qval.cuts)
  }
  return(res)
}

k.ind <- 3
n.genes.trim <- 5000 # No of genes to trim down to
n.diff.trim <- 1000 # No of DE genes to trim down to
qval.cuts <- seq(0, 0.15, by = 0.001)
  
mainDir <- getwd()
subDir <- "pval_output"

pvals.deseq2.simseq <- readRDS(file.path(mainDir, subDir, paste0("ss", k.ind), "p_deseq2_simseq.RDS"))
pvals.quasiseq.simseq <- readRDS(file.path(mainDir, subDir, paste0("ss", k.ind), "p_quasiseq_simseq.RDS"))
pvals.edger.simseq <- readRDS(file.path(mainDir, subDir, paste0("ss", k.ind), "p_edger_simseq.RDS"))
pvals.samseq.simseq <- readRDS(file.path(mainDir, subDir, paste0("ss", k.ind), "p_samseq_simseq.RDS"))
pvals.voom.simseq <- readRDS(file.path(mainDir, subDir, paste0("ss", k.ind), "p_voom_simseq.RDS"))

pvals.deseq2.nb <- readRDS(file.path(mainDir, subDir, paste0("ss", k.ind), "p_deseq2_nb.RDS"))
pvals.quasiseq.nb <- readRDS(file.path(mainDir, subDir, paste0("ss", k.ind), "p_quasiseq_nb.RDS"))
pvals.edger.nb <- readRDS(file.path(mainDir, subDir, paste0("ss", k.ind), "p_edger_nb.RDS"))
pvals.samseq.nb <- readRDS(file.path(mainDir, subDir, paste0("ss", k.ind), "p_samseq_nb.RDS"))
pvals.voom.nb <- readRDS(file.path(mainDir, subDir, paste0("ss", k.ind), "p_voom_nb.RDS"))

if( k.ind >= 7)
{
  pvals.deseq2.cooks.simseq <- readRDS(file.path(mainDir, subDir, paste0("ss", k.ind), "p_deseq2_cooks_simseq.RDS"))
  pvals.quasiseq.cooks.simseq <- readRDS(file.path(mainDir, subDir, paste0("ss", k.ind), "p_quasiseq_cooks_simseq.RDS"))
  pvals.edger.cooks.simseq <- readRDS(file.path(mainDir, subDir, paste0("ss", k.ind), "p_edger_cooks_simseq.RDS"))
  pvals.samseq.cooks.simseq <- readRDS(file.path(mainDir, subDir, paste0("ss", k.ind), "p_samseq_cooks_simseq.RDS"))
  pvals.voom.cooks.simseq <- readRDS(file.path(mainDir, subDir, paste0("ss", k.ind), "p_voom_cooks_simseq.RDS"))
  
  pvals.deseq2.cooks.nb <- readRDS(file.path(mainDir, subDir, paste0("ss", k.ind), "p_deseq2_cooks_nb.RDS"))
  pvals.quasiseq.cooks.nb <- readRDS(file.path(mainDir, subDir, paste0("ss", k.ind), "p_quasiseq_cooks_nb.RDS"))
  pvals.edger.cooks.nb <- readRDS(file.path(mainDir, subDir, paste0("ss", k.ind), "p_edger_cooks_nb.RDS"))
  pvals.samseq.cooks.nb <- readRDS(file.path(mainDir, subDir, paste0("ss", k.ind), "p_samseq_cooks_nb.RDS"))
  pvals.voom.cooks.nb <- readRDS(file.path(mainDir, subDir, paste0("ss", k.ind), "p_voom_cooks_nb.RDS"))  
}

filt.cooks.simseq <- readRDS(file.path(mainDir, subDir, paste0("ss", k.ind), "filt_cooks_simseq.RDS"))
filt.cooks.nb <- readRDS(file.path(mainDir, subDir, paste0("ss", k.ind), "filt_cooks_nb.RDS"))


### Compute PAU
null.genes <- c( rep(TRUE, n.genes.trim - n.diff.trim), rep(FALSE, n.diff.trim))
pau.deseq2.simseq <- apply(pvals.deseq2.simseq, 1, pau, null.genes, specificity = 0.05)
pau.edger.simseq <- apply(pvals.edger.simseq, 1, pau, null.genes, specificity = 0.05)
pau.quasiseq.simseq <- apply(pvals.quasiseq.simseq, 1, pau, null.genes, specificity = 0.05)
pau.samseq.simseq <- apply(pvals.samseq.simseq, 1, pau, null.genes, specificity = 0.05)
pau.voom.simseq <- apply(pvals.voom.simseq, 1, pau, null.genes, specificity = 0.05)

pau.deseq2.nb <- apply(pvals.deseq2.nb, 1, pau, null.genes, specificity = 0.05)
pau.edger.nb <- apply(pvals.edger.nb, 1, pau, null.genes, specificity = 0.05)
pau.quasiseq.nb <- apply(pvals.quasiseq.nb, 1, pau, null.genes, specificity = 0.05)
pau.samseq.nb <- apply(pvals.samseq.nb, 1, pau, null.genes, specificity = 0.05)
pau.voom.nb <- apply(pvals.voom.nb, 1, pau, null.genes, specificity = 0.05)

if(k.ind < 7)
{
  pau.deseq2.cooks.simseq <- pau.cooks(pvals.deseq2.simseq, filt.cooks.simseq, null.genes, specificity = 0.05)
  pau.edger.cooks.simseq <- pau.cooks(pvals.edger.simseq, filt.cooks.simseq, null.genes, specificity = 0.05)
  pau.quasiseq.cooks.simseq <- pau.cooks(pvals.quasiseq.simseq, filt.cooks.simseq, null.genes, specificity = 0.05)
  pau.samseq.cooks.simseq <- pau.cooks(pvals.samseq.simseq, filt.cooks.simseq, null.genes, specificity = 0.05)
  pau.voom.cooks.simseq <- pau.cooks(pvals.voom.simseq, filt.cooks.simseq, null.genes, specificity = 0.05)
  
  pau.deseq2.cooks.nb <- pau.cooks(pvals.deseq2.nb, filt.cooks.nb, null.genes, specificity = 0.05)
  pau.edger.cooks.nb <- pau.cooks(pvals.edger.nb, filt.cooks.nb, null.genes, specificity = 0.05)
  pau.quasiseq.cooks.nb <- pau.cooks(pvals.quasiseq.nb, filt.cooks.nb, null.genes, specificity = 0.05)
  pau.samseq.cooks.nb <- pau.cooks(pvals.samseq.nb, filt.cooks.nb, null.genes, specificity = 0.05)
  pau.voom.cooks.nb <- pau.cooks(pvals.voom.nb, filt.cooks.nb, null.genes, specificity = 0.05)
} else {
  pau.deseq2.cooks.simseq <- pau.cooks(pvals.deseq2.cooks.simseq, filt.cooks.simseq, null.genes, specificity = 0.05, pvals.list = TRUE)
  pau.edger.cooks.simseq <- pau.cooks(pvals.edger.cooks.simseq, filt.cooks.simseq, null.genes, specificity = 0.05, pvals.list = TRUE)
  pau.quasiseq.cooks.simseq <- pau.cooks(pvals.quasiseq.cooks.simseq, filt.cooks.simseq, null.genes, specificity = 0.05, pvals.list = TRUE)
  pau.samseq.cooks.simseq <- pau.cooks(pvals.samseq.cooks.simseq, filt.cooks.simseq, null.genes, specificity = 0.05, pvals.list = TRUE)
  pau.voom.cooks.simseq <- pau.cooks(pvals.voom.cooks.simseq, filt.cooks.simseq, null.genes, specificity = 0.05, pvals.list = TRUE)
  
  pau.deseq2.cooks.nb <- pau.cooks(pvals.deseq2.cooks.nb, filt.cooks.nb, null.genes, specificity = 0.05, pvals.list = TRUE)
  pau.edger.cooks.nb <- pau.cooks(pvals.edger.cooks.nb, filt.cooks.nb, null.genes, specificity = 0.05, pvals.list = TRUE)
  pau.quasiseq.cooks.nb <- pau.cooks(pvals.quasiseq.cooks.nb, filt.cooks.nb, null.genes, specificity = 0.05, pvals.list = TRUE)
  pau.samseq.cooks.nb <- pau.cooks(pvals.samseq.cooks.nb, filt.cooks.nb, null.genes, specificity = 0.05, pvals.list = TRUE)
  pau.voom.cooks.nb <- pau.cooks(pvals.voom.cooks.nb, filt.cooks.nb, null.genes, specificity = 0.05, pvals.list = TRUE)
}

fdp.deseq2.simseq <- apply(pvals.deseq2.simseq, 1, fdp, null.genes, qval.cuts)
fdp.edger.simseq <- apply(pvals.edger.simseq, 1, fdp, null.genes, qval.cuts)
fdp.quasiseq.simseq <- apply(pvals.quasiseq.simseq, 1, fdp, null.genes, qval.cuts)
fdp.samseq.simseq <- apply(pvals.samseq.simseq, 1, fdp, null.genes, qval.cuts)
fdp.voom.simseq <- apply(pvals.voom.simseq, 1, fdp, null.genes, qval.cuts)

fdp.deseq2.nb <- apply(pvals.deseq2.nb, 1, fdp, null.genes, qval.cuts)
fdp.edger.nb <- apply(pvals.edger.nb, 1, fdp, null.genes, qval.cuts)
fdp.quasiseq.nb <- apply(pvals.quasiseq.nb, 1, fdp, null.genes, qval.cuts)
fdp.samseq.nb <- apply(pvals.samseq.nb, 1, fdp, null.genes, qval.cuts)
fdp.voom.nb <- apply(pvals.voom.nb, 1, fdp, null.genes, qval.cuts)

if(k.ind < 7)
{
  fdp.deseq2.cooks.simseq <- fdp.cooks(pvals.deseq2.simseq, filt.cooks.simseq, null.genes, qval.cuts)
  fdp.edger.cooks.simseq <- fdp.cooks(pvals.edger.simseq, filt.cooks.simseq, null.genes, qval.cuts)
  fdp.quasiseq.cooks.simseq <- fdp.cooks(pvals.quasiseq.simseq, filt.cooks.simseq, null.genes, qval.cuts)
  fdp.samseq.cooks.simseq <- fdp.cooks(pvals.samseq.simseq, filt.cooks.simseq, null.genes, qval.cuts)
  fdp.voom.cooks.simseq <- fdp.cooks(pvals.voom.simseq, filt.cooks.simseq, null.genes, qval.cuts)
  
  fdp.deseq2.cooks.nb <- fdp.cooks(pvals.deseq2.nb, filt.cooks.nb, null.genes, qval.cuts)
  fdp.edger.cooks.nb <- fdp.cooks(pvals.edger.nb, filt.cooks.nb, null.genes, qval.cuts)
  fdp.quasiseq.cooks.nb <- fdp.cooks(pvals.quasiseq.nb, filt.cooks.nb, null.genes, qval.cuts)
  fdp.samseq.cooks.nb <- fdp.cooks(pvals.samseq.nb, filt.cooks.nb, null.genes, qval.cuts)
  fdp.voom.cooks.nb <- fdp.cooks(pvals.voom.nb, filt.cooks.nb, null.genes, qval.cuts)
} else {
  fdp.deseq2.cooks.simseq <- fdp.cooks(pvals.deseq2.cooks.simseq, filt.cooks.simseq, null.genes, qval.cuts, pvals.list = TRUE)
  fdp.edger.cooks.simseq <- fdp.cooks(pvals.edger.cooks.simseq, filt.cooks.simseq, null.genes, qval.cuts, pvals.list = TRUE)
  fdp.quasiseq.cooks.simseq <- fdp.cooks(pvals.quasiseq.cooks.simseq, filt.cooks.simseq, null.genes, qval.cuts, pvals.list = TRUE)
  fdp.samseq.cooks.simseq <- fdp.cooks(pvals.samseq.cooks.simseq, filt.cooks.simseq, null.genes, qval.cuts, pvals.list = TRUE)
  fdp.voom.cooks.simseq <- fdp.cooks(pvals.voom.cooks.simseq, filt.cooks.simseq, null.genes, qval.cuts, pvals.list = TRUE)
  
  fdp.deseq2.cooks.nb <- fdp.cooks(pvals.deseq2.cooks.nb, filt.cooks.nb, null.genes, qval.cuts, pvals.list = TRUE)
  fdp.edger.cooks.nb <- fdp.cooks(pvals.edger.cooks.nb, filt.cooks.nb, null.genes, qval.cuts, pvals.list = TRUE)
  fdp.quasiseq.cooks.nb <- fdp.cooks(pvals.quasiseq.cooks.nb, filt.cooks.nb, null.genes, qval.cuts, pvals.list = TRUE)
  fdp.samseq.cooks.nb <- fdp.cooks(pvals.samseq.cooks.nb, filt.cooks.nb, null.genes, qval.cuts, pvals.list = TRUE)
  fdp.voom.cooks.nb <- fdp.cooks(pvals.voom.cooks.nb, filt.cooks.nb, null.genes, qval.cuts, pvals.list = TRUE)
}

mainDir <- getwd()
subDir <- "fdp_output"
if( !file.exists(file.path(mainDir, subDir)) )
{
  dir.create(file.path(mainDir, subDir)) 
}
if( !file.exists(file.path(mainDir, subDir, paste0("ss", k.ind))) )
{
  dir.create(file.path(mainDir, subDir, paste0("ss", k.ind))) 
}
saveRDS(fdp.deseq2.simseq, file.path(mainDir, subDir, paste0("ss", k.ind), "fdp_deseq2_simseq.RDS"))
saveRDS(fdp.quasiseq.simseq, file.path(mainDir, subDir, paste0("ss", k.ind), "fdp_quasiseq_simseq.RDS"))
saveRDS(fdp.edger.simseq, file.path(mainDir, subDir, paste0("ss", k.ind), "fdp_edger_simseq.RDS"))
saveRDS(fdp.samseq.simseq, file.path(mainDir, subDir, paste0("ss", k.ind), "fdp_samseq_simseq.RDS"))
saveRDS(fdp.voom.simseq, file.path(mainDir, subDir, paste0("ss", k.ind), "fdp_voom_simseq.RDS"))

saveRDS(fdp.deseq2.nb, file.path(mainDir, subDir, paste0("ss", k.ind), "fdp_deseq2_nb.RDS"))
saveRDS(fdp.quasiseq.nb, file.path(mainDir, subDir, paste0("ss", k.ind), "fdp_quasiseq_nb.RDS"))
saveRDS(fdp.edger.nb, file.path(mainDir, subDir, paste0("ss", k.ind), "fdp_edger_nb.RDS"))
saveRDS(fdp.samseq.nb, file.path(mainDir, subDir, paste0("ss", k.ind), "fdp_samseq_nb.RDS"))
saveRDS(fdp.voom.nb, file.path(mainDir, subDir, paste0("ss", k.ind), "fdp_voom_nb.RDS"))

saveRDS(fdp.deseq2.cooks.simseq, file.path(mainDir, subDir, paste0("ss", k.ind), "fdp_deseq2_cooks_simseq.RDS"))
saveRDS(fdp.quasiseq.cooks.simseq, file.path(mainDir, subDir, paste0("ss", k.ind), "fdp_quasiseq_cooks_simseq.RDS"))
saveRDS(fdp.edger.cooks.simseq, file.path(mainDir, subDir, paste0("ss", k.ind), "fdp_edger_cooks_simseq.RDS"))
saveRDS(fdp.samseq.cooks.simseq, file.path(mainDir, subDir, paste0("ss", k.ind), "fdp_samseq_cooks_simseq.RDS"))
saveRDS(fdp.voom.cooks.simseq, file.path(mainDir, subDir, paste0("ss", k.ind), "fdp_voom_cooks_simseq.RDS"))

saveRDS(fdp.deseq2.cooks.nb, file.path(mainDir, subDir, paste0("ss", k.ind), "fdp_deseq2_cooks_nb.RDS"))
saveRDS(fdp.quasiseq.cooks.nb, file.path(mainDir, subDir, paste0("ss", k.ind), "fdp_quasiseq_cooks_nb.RDS"))
saveRDS(fdp.edger.cooks.nb, file.path(mainDir, subDir, paste0("ss", k.ind), "fdp_edger_cooks_nb.RDS"))
saveRDS(fdp.samseq.cooks.nb, file.path(mainDir, subDir, paste0("ss", k.ind), "fdp_samseq_cooks_nb.RDS"))
saveRDS(fdp.voom.cooks.nb, file.path(mainDir, subDir, paste0("ss", k.ind), "fdp_voom_cooks_nb.RDS"))

subDir <- "pau_output"
if( !file.exists(file.path(mainDir, subDir)) )
{
  dir.create(file.path(mainDir, subDir)) 
}
if( !file.exists(file.path(mainDir, subDir, paste0("ss", k.ind))) )
{
  dir.create(file.path(mainDir, subDir, paste0("ss", k.ind))) 
}
saveRDS(pau.deseq2.simseq, file.path(mainDir, subDir, paste0("ss", k.ind), "pau_deseq2_simseq.RDS"))
saveRDS(pau.quasiseq.simseq, file.path(mainDir, subDir, paste0("ss", k.ind), "pau_quasiseq_simseq.RDS"))
saveRDS(pau.edger.simseq, file.path(mainDir, subDir, paste0("ss", k.ind), "pau_edger_simseq.RDS"))
saveRDS(pau.samseq.simseq, file.path(mainDir, subDir, paste0("ss", k.ind), "pau_samseq_simseq.RDS"))
saveRDS(pau.voom.simseq, file.path(mainDir, subDir, paste0("ss", k.ind), "pau_voom_simseq.RDS"))

saveRDS(pau.deseq2.nb, file.path(mainDir, subDir, paste0("ss", k.ind), "pau_deseq2_nb.RDS"))
saveRDS(pau.quasiseq.nb, file.path(mainDir, subDir, paste0("ss", k.ind), "pau_quasiseq_nb.RDS"))
saveRDS(pau.edger.nb, file.path(mainDir, subDir, paste0("ss", k.ind), "pau_edger_nb.RDS"))
saveRDS(pau.samseq.nb, file.path(mainDir, subDir, paste0("ss", k.ind), "pau_samseq_nb.RDS"))
saveRDS(pau.voom.nb, file.path(mainDir, subDir, paste0("ss", k.ind), "pau_voom_nb.RDS"))

saveRDS(pau.deseq2.cooks.simseq, file.path(mainDir, subDir, paste0("ss", k.ind), "pau_deseq2_cooks_simseq.RDS"))
saveRDS(pau.quasiseq.cooks.simseq, file.path(mainDir, subDir, paste0("ss", k.ind), "pau_quasiseq_cooks_simseq.RDS"))
saveRDS(pau.edger.cooks.simseq, file.path(mainDir, subDir, paste0("ss", k.ind), "pau_edger_cooks_simseq.RDS"))
saveRDS(pau.samseq.cooks.simseq, file.path(mainDir, subDir, paste0("ss", k.ind), "pau_samseq_cooks_simseq.RDS"))
saveRDS(pau.voom.cooks.simseq, file.path(mainDir, subDir, paste0("ss", k.ind), "pau_voom_cooks_simseq.RDS"))

saveRDS(pau.deseq2.cooks.nb, file.path(mainDir, subDir, paste0("ss", k.ind), "pau_deseq2_cooks_nb.RDS"))
saveRDS(pau.quasiseq.cooks.nb, file.path(mainDir, subDir, paste0("ss", k.ind), "pau_quasiseq_cooks_nb.RDS"))
saveRDS(pau.edger.cooks.nb, file.path(mainDir, subDir, paste0("ss", k.ind), "pau_edger_cooks_nb.RDS"))
saveRDS(pau.samseq.cooks.nb, file.path(mainDir, subDir, paste0("ss", k.ind), "pau_samseq_cooks_nb.RDS"))
saveRDS(pau.voom.cooks.nb, file.path(mainDir, subDir, paste0("ss", k.ind), "pau_voom_cooks_nb.RDS"))