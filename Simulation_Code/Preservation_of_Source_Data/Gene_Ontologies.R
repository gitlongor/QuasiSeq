### Load Bioconductor packages
require(org.Hs.eg.db)
require(annotate)

### Load CRAN packages
require(stringr)

### Load Gene annotations
gene_annot <- readRDS("gene_annotations.RDS")
gene_fam <- strsplit(gene_annot, "|", fixed = TRUE)
gene_first <- sapply(gene_fam, function(x) x[1])
gene_sec <- sapply(gene_fam, function(x) x[2])

### Create Gene Ontologies for each gene annotation
l <- length(gene_sec)
GO.Hs <- vector("list", l)
for(w in 1:l){
  nam <- names(getGO(gene_sec[w], 'org.Hs.eg')[[1]])
  if(!is.null(nam)) GO.Hs[[w]] <- nam
}

### Save Results
saveRDS(GO.Hs, "Gene_Ontologies_Hs.RDS")
