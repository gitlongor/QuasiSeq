print(getwd())
setwd("C:/Users/Klirk/Desktop/Simulation_Code/Simulation_Code/KIRC_Simulations/foldchange")
require(ggplot2)
require(SimSeq)
require(edgeR)

### Set seed  
set.seed(1039245)

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

logLambdaFoldChange = lambdas[, 2] - lambdas[, 1]
logLambdaFoldChangeSim <- matrix(0, 200, 5000)

mainDir <- getwd()
subDir <- "pval_output"

k.ind <- 5
fold.quasiseq.simseq.fit1 <- readRDS(file.path(mainDir, subDir, paste0("ss", k.ind), "fold_quasiseq_simseq_fit1.RDS"))
fold.quasiseq.simseq.fit2 <- readRDS(file.path(mainDir, subDir, paste0("ss", k.ind), "fold_quasiseq_simseq_fit2.RDS"))
fold.quasiseq.simseq.fit3 <- readRDS(file.path(mainDir, subDir, paste0("ss", k.ind), "fold_quasiseq_simseq_fit3.RDS"))
fold.voom.simseq <- readRDS(file.path(mainDir, subDir, paste0("ss", k.ind), "fold_voom_simseq.RDS"))
fold.quasiseq.nb.fit1 <- readRDS(file.path(mainDir, subDir, paste0("ss", k.ind), "fold_quasiseq_nb_fit1.RDS"))
fold.quasiseq.nb.fit2 <- readRDS(file.path(mainDir, subDir, paste0("ss", k.ind), "fold_quasiseq_nb_fit2.RDS"))
fold.quasiseq.nb.fit3 <- readRDS(file.path(mainDir, subDir, paste0("ss", k.ind), "fold_quasiseq_nb_fit3.RDS"))
fold.voom.nb <- readRDS(file.path(mainDir, subDir, paste0("ss", k.ind), "fold_voom_nb.RDS"))
fold.trulyDEgeneList <- readRDS(file.path(mainDir, subDir, paste0("ss", k.ind), "fold_trulyDEgeneList.RDS"))

#finding de genes for each simulation in logLambdaFoldChange
for(ii in 1:200){
	for(jj in 1:1000){
		logLambdaFoldChangeSim[ii,4000+jj] = logLambdaFoldChange[fold.trulyDEgeneList[ii,jj]]
	}
}

#convert to log2
fold.quasiseq.simseq.fit1 <- log2(exp(fold.quasiseq.simseq.fit1*2))
fold.quasiseq.simseq.fit2 <- log2(exp(fold.quasiseq.simseq.fit2*2))
fold.quasiseq.simseq.fit3 <- log2(exp(fold.quasiseq.simseq.fit3*2))
fold.voom.simseq <- fold.voom.simseq*2
fold.quasiseq.nb.fit1 <- log2(exp(fold.quasiseq.nb.fit1*2))
fold.quasiseq.nb.fit2 <- log2(exp(fold.quasiseq.nb.fit2*2))
fold.quasiseq.nb.fit3 <- log2(exp(fold.quasiseq.nb.fit3*2))
fold.voom.nb <- fold.voom.nb*2

#calculate estimate errors
fold.quasiseq.simseq.fit1.estError = fold.quasiseq.simseq.fit1 - logLambdaFoldChangeSim
fold.quasiseq.simseq.fit2.estError = fold.quasiseq.simseq.fit2 - logLambdaFoldChangeSim
fold.quasiseq.simseq.fit3.estError = fold.quasiseq.simseq.fit3 - logLambdaFoldChangeSim
fold.voom.simseq.estError = fold.voom.simseq - logLambdaFoldChangeSim
fold.quasiseq.nb.fit1.estError = fold.quasiseq.nb.fit1 - logLambdaFoldChangeSim
fold.quasiseq.nb.fit2.estError = fold.quasiseq.nb.fit2 - logLambdaFoldChangeSim
fold.quasiseq.nb.fit3.estError = fold.quasiseq.nb.fit3 - logLambdaFoldChangeSim
fold.voom.nb.estError = fold.voom.nb - logLambdaFoldChangeSim

#separate ee genes
fold.quasiseq.simseq.fit1.estErrorEE = fold.quasiseq.simseq.fit1[,1:4000] - logLambdaFoldChangeSim[,1:4000]
fold.quasiseq.simseq.fit2.estErrorEE = fold.quasiseq.simseq.fit2[,1:4000] - logLambdaFoldChangeSim[,1:4000]
fold.quasiseq.simseq.fit3.estErrorEE = fold.quasiseq.simseq.fit3[,1:4000] - logLambdaFoldChangeSim[,1:4000]
fold.voom.simseq.estErrorEE = fold.voom.simseq[,1:4000] - logLambdaFoldChangeSim[,1:4000]
fold.quasiseq.nb.fit1.estErrorEE = fold.quasiseq.nb.fit1[,1:4000] - logLambdaFoldChangeSim[,1:4000]
fold.quasiseq.nb.fit2.estErrorEE = fold.quasiseq.nb.fit2[,1:4000] - logLambdaFoldChangeSim[,1:4000]
fold.quasiseq.nb.fit3.estErrorEE = fold.quasiseq.nb.fit3[,1:4000] - logLambdaFoldChangeSim[,1:4000]
fold.voom.nb.estErrorEE = fold.voom.nb[,1:4000] - logLambdaFoldChangeSim[,1:4000]

#separate de genes
fold.quasiseq.simseq.fit1.estErrorDE = fold.quasiseq.simseq.fit1[,4001:5000] - logLambdaFoldChangeSim[,4001:5000]
fold.quasiseq.simseq.fit2.estErrorDE = fold.quasiseq.simseq.fit2[,4001:5000] - logLambdaFoldChangeSim[,4001:5000]
fold.quasiseq.simseq.fit3.estErrorDE = fold.quasiseq.simseq.fit3[,4001:5000] - logLambdaFoldChangeSim[,4001:5000]
fold.voom.simseq.estErrorDE = fold.voom.simseq[,4001:5000] - logLambdaFoldChangeSim[,4001:5000]
fold.quasiseq.nb.fit1.estErrorDE = fold.quasiseq.nb.fit1[,4001:5000] - logLambdaFoldChangeSim[,4001:5000]
fold.quasiseq.nb.fit2.estErrorDE = fold.quasiseq.nb.fit2[,4001:5000] - logLambdaFoldChangeSim[,4001:5000]
fold.quasiseq.nb.fit3.estErrorDE = fold.quasiseq.nb.fit3[,4001:5000] - logLambdaFoldChangeSim[,4001:5000]
fold.voom.nb.estErrorDE = fold.voom.nb[,4001:5000] - logLambdaFoldChangeSim[,4001:5000]

#calculate row means all
fold.quasiseq.simseq.fit1.mean = rowMeans(fold.quasiseq.simseq.fit1.estError)
fold.quasiseq.simseq.fit2.mean = rowMeans(fold.quasiseq.simseq.fit2.estError)
fold.quasiseq.simseq.fit3.mean = rowMeans(fold.quasiseq.simseq.fit3.estError)
fold.voom.simseq.mean = rowMeans(fold.voom.simseq.estError)
fold.quasiseq.nb.fit1.mean = rowMeans(fold.quasiseq.nb.fit1.estError)
fold.quasiseq.nb.fit2.mean = rowMeans(fold.quasiseq.nb.fit2.estError)
fold.quasiseq.nb.fit3.mean = rowMeans(fold.quasiseq.nb.fit3.estError)
fold.voom.nb.mean = rowMeans(fold.voom.nb.estError)

#calculate row means ee genes
fold.quasiseq.simseq.fit1.mean.ee = rowMeans(fold.quasiseq.simseq.fit1.estErrorEE)
fold.quasiseq.simseq.fit2.mean.ee = rowMeans(fold.quasiseq.simseq.fit2.estErrorEE)
fold.quasiseq.simseq.fit3.mean.ee = rowMeans(fold.quasiseq.simseq.fit3.estErrorEE)
fold.voom.simseq.mean.ee = rowMeans(fold.voom.simseq.estErrorEE)
fold.quasiseq.nb.fit1.mean.ee = rowMeans(fold.quasiseq.nb.fit1.estErrorEE)
fold.quasiseq.nb.fit2.mean.ee = rowMeans(fold.quasiseq.nb.fit2.estErrorEE)
fold.quasiseq.nb.fit3.mean.ee = rowMeans(fold.quasiseq.nb.fit3.estErrorEE)
fold.voom.nb.mean.ee = rowMeans(fold.voom.nb.estErrorEE)

#calculate row means de genes
fold.quasiseq.simseq.fit1.mean.de = rowMeans(fold.quasiseq.simseq.fit1.estErrorDE)
fold.quasiseq.simseq.fit2.mean.de = rowMeans(fold.quasiseq.simseq.fit2.estErrorDE)
fold.quasiseq.simseq.fit3.mean.de = rowMeans(fold.quasiseq.simseq.fit3.estErrorDE)
fold.voom.simseq.mean.de = rowMeans(fold.voom.simseq.estErrorDE)
fold.quasiseq.nb.fit1.mean.de = rowMeans(fold.quasiseq.nb.fit1.estErrorDE)
fold.quasiseq.nb.fit2.mean.de = rowMeans(fold.quasiseq.nb.fit2.estErrorDE)
fold.quasiseq.nb.fit3.mean.de = rowMeans(fold.quasiseq.nb.fit3.estErrorDE)
fold.voom.nb.mean.de = rowMeans(fold.voom.nb.estErrorDE)

l <- length(fold.quasiseq.simseq.fit1.mean)
num1 <- data.frame(num = c(fold.quasiseq.simseq.fit1.mean, fold.quasiseq.simseq.fit2.mean, fold.quasiseq.simseq.fit3.mean, fold.voom.simseq.mean), 
                   Method = as.factor(c(rep("Bias Fold Tolerance = 1.10", l), rep("Bias Fold Tolerance = 1", l), rep("Bias Fold Tolerance = Inf", l)
							, rep("Voom", l))),
                   Simulation = rep("SimSeq", l*4),
                   SS = rep("Sample Size: 5",l*4))

l <- length(fold.quasiseq.nb.fit1.mean)
num2 <- data.frame(num = c(fold.quasiseq.nb.fit1.mean, fold.quasiseq.nb.fit2.mean, fold.quasiseq.nb.fit3.mean, fold.voom.nb.mean), 
                   Method = as.factor(c(rep("Bias Fold Tolerance = 1.10", l), rep("Bias Fold Tolerance = 1", l), rep("Bias Fold Tolerance = Inf", l)
							, rep("Voom", l))),
                   Simulation = rep("NB", l*4),
                   SS = rep("Sample Size: 5",l*4)) 
				   
l <- length(fold.quasiseq.simseq.fit1.mean)
num1EE <- data.frame(num = c(fold.quasiseq.simseq.fit1.mean.ee, fold.quasiseq.simseq.fit2.mean.ee, fold.quasiseq.simseq.fit3.mean.ee, fold.voom.simseq.mean.ee), 
                   Method = as.factor(c(rep("Bias Fold Tolerance = 1.10", l), rep("Bias Fold Tolerance = 1", l), rep("Bias Fold Tolerance = Inf", l)
							, rep("Voom", l))),
                   Simulation = rep("SimSeq", l*4),
                   SS = rep("Sample Size: 5",l*4))

l <- length(fold.quasiseq.nb.fit1.mean)
num2EE <- data.frame(num = c(fold.quasiseq.nb.fit1.mean.ee, fold.quasiseq.nb.fit2.mean.ee, fold.quasiseq.nb.fit3.mean.ee, fold.voom.nb.mean.ee), 
                   Method = as.factor(c(rep("Bias Fold Tolerance = 1.10", l), rep("Bias Fold Tolerance = 1", l), rep("Bias Fold Tolerance = Inf", l)
							, rep("Voom", l))),
                   Simulation = rep("NB", l*4),
                   SS = rep("Sample Size: 5",l*4)) 
				   
l <- length(fold.quasiseq.simseq.fit1.mean)
num1DE <- data.frame(num = c(fold.quasiseq.simseq.fit1.mean.de, fold.quasiseq.simseq.fit2.mean.de, fold.quasiseq.simseq.fit3.mean.de, fold.voom.simseq.mean.de), 
                   Method = as.factor(c(rep("Bias Fold Tolerance = 1.10", l), rep("Bias Fold Tolerance = 1", l), rep("Bias Fold Tolerance = Inf", l)
							, rep("Voom", l))),
                   Simulation = rep("SimSeq", l*4),
                   SS = rep("Sample Size: 5",l*4))

l <- length(fold.quasiseq.nb.fit1.mean)
num2DE <- data.frame(num = c(fold.quasiseq.nb.fit1.mean.de, fold.quasiseq.nb.fit2.mean.de, fold.quasiseq.nb.fit3.mean.de, fold.voom.nb.mean.de), 
                   Method = as.factor(c(rep("Bias Fold Tolerance = 1.10", l), rep("Bias Fold Tolerance = 1", l), rep("Bias Fold Tolerance = Inf", l)
							, rep("Voom", l))),
                   Simulation = rep("NB", l*4),
                   SS = rep("Sample Size: 5",l*4)) 

k.ind <- 10
fold.quasiseq.simseq.fit1 <- readRDS(file.path(mainDir, subDir, paste0("ss", k.ind), "fold_quasiseq_simseq_fit1.RDS"))
fold.quasiseq.simseq.fit2 <- readRDS(file.path(mainDir, subDir, paste0("ss", k.ind), "fold_quasiseq_simseq_fit2.RDS"))
fold.quasiseq.simseq.fit3 <- readRDS(file.path(mainDir, subDir, paste0("ss", k.ind), "fold_quasiseq_simseq_fit3.RDS"))
fold.voom.simseq <- readRDS(file.path(mainDir, subDir, paste0("ss", k.ind), "fold_voom_simseq.RDS"))
fold.quasiseq.nb.fit1 <- readRDS(file.path(mainDir, subDir, paste0("ss", k.ind), "fold_quasiseq_nb_fit1.RDS"))
fold.quasiseq.nb.fit2 <- readRDS(file.path(mainDir, subDir, paste0("ss", k.ind), "fold_quasiseq_nb_fit2.RDS"))
fold.quasiseq.nb.fit3 <- readRDS(file.path(mainDir, subDir, paste0("ss", k.ind), "fold_quasiseq_nb_fit3.RDS"))
fold.voom.nb <- readRDS(file.path(mainDir, subDir, paste0("ss", k.ind), "fold_voom_nb.RDS"))
fold.trulyDEgeneList <- readRDS(file.path(mainDir, subDir, paste0("ss", k.ind), "fold_trulyDEgeneList.RDS"))

#finding de genes for each simulation in logLambdaFoldChange
for(ii in 1:200){
	for(jj in 1:1000){
		logLambdaFoldChangeSim[ii,4000+jj] = logLambdaFoldChange[fold.trulyDEgeneList[ii,jj]]
	}
}

#convert to log2
fold.quasiseq.simseq.fit1 <- log2(exp(fold.quasiseq.simseq.fit1*2))
fold.quasiseq.simseq.fit2 <- log2(exp(fold.quasiseq.simseq.fit2*2))
fold.quasiseq.simseq.fit3 <- log2(exp(fold.quasiseq.simseq.fit3*2))
fold.voom.simseq <- fold.voom.simseq*2
fold.quasiseq.nb.fit1 <- log2(exp(fold.quasiseq.nb.fit1*2))
fold.quasiseq.nb.fit2 <- log2(exp(fold.quasiseq.nb.fit2*2))
fold.quasiseq.nb.fit3 <- log2(exp(fold.quasiseq.nb.fit3*2))
fold.voom.nb <- fold.voom.nb*2

#calculate estimate errors
fold.quasiseq.simseq.fit1.estError = fold.quasiseq.simseq.fit1 - logLambdaFoldChangeSim
fold.quasiseq.simseq.fit2.estError = fold.quasiseq.simseq.fit2 - logLambdaFoldChangeSim
fold.quasiseq.simseq.fit3.estError = fold.quasiseq.simseq.fit3 - logLambdaFoldChangeSim
fold.voom.simseq.estError = fold.voom.simseq - logLambdaFoldChangeSim
fold.quasiseq.nb.fit1.estError = fold.quasiseq.nb.fit1 - logLambdaFoldChangeSim
fold.quasiseq.nb.fit2.estError = fold.quasiseq.nb.fit2 - logLambdaFoldChangeSim
fold.quasiseq.nb.fit3.estError = fold.quasiseq.nb.fit3 - logLambdaFoldChangeSim
fold.voom.nb.estError = fold.voom.nb - logLambdaFoldChangeSim

#separate ee genes
fold.quasiseq.simseq.fit1.estErrorEE = fold.quasiseq.simseq.fit1[,1:4000] - logLambdaFoldChangeSim[,1:4000]
fold.quasiseq.simseq.fit2.estErrorEE = fold.quasiseq.simseq.fit2[,1:4000] - logLambdaFoldChangeSim[,1:4000]
fold.quasiseq.simseq.fit3.estErrorEE = fold.quasiseq.simseq.fit3[,1:4000] - logLambdaFoldChangeSim[,1:4000]
fold.voom.simseq.estErrorEE = fold.voom.simseq[,1:4000] - logLambdaFoldChangeSim[,1:4000]
fold.quasiseq.nb.fit1.estErrorEE = fold.quasiseq.nb.fit1[,1:4000] - logLambdaFoldChangeSim[,1:4000]
fold.quasiseq.nb.fit2.estErrorEE = fold.quasiseq.nb.fit2[,1:4000] - logLambdaFoldChangeSim[,1:4000]
fold.quasiseq.nb.fit3.estErrorEE = fold.quasiseq.nb.fit3[,1:4000] - logLambdaFoldChangeSim[,1:4000]
fold.voom.nb.estErrorEE = fold.voom.nb[,1:4000] - logLambdaFoldChangeSim[,1:4000]

#separate de genes
fold.quasiseq.simseq.fit1.estErrorDE = fold.quasiseq.simseq.fit1[,4001:5000] - logLambdaFoldChangeSim[,4001:5000]
fold.quasiseq.simseq.fit2.estErrorDE = fold.quasiseq.simseq.fit2[,4001:5000] - logLambdaFoldChangeSim[,4001:5000]
fold.quasiseq.simseq.fit3.estErrorDE = fold.quasiseq.simseq.fit3[,4001:5000] - logLambdaFoldChangeSim[,4001:5000]
fold.voom.simseq.estErrorDE = fold.voom.simseq[,4001:5000] - logLambdaFoldChangeSim[,4001:5000]
fold.quasiseq.nb.fit1.estErrorDE = fold.quasiseq.nb.fit1[,4001:5000] - logLambdaFoldChangeSim[,4001:5000]
fold.quasiseq.nb.fit2.estErrorDE = fold.quasiseq.nb.fit2[,4001:5000] - logLambdaFoldChangeSim[,4001:5000]
fold.quasiseq.nb.fit3.estErrorDE = fold.quasiseq.nb.fit3[,4001:5000] - logLambdaFoldChangeSim[,4001:5000]
fold.voom.nb.estErrorDE = fold.voom.nb[,4001:5000] - logLambdaFoldChangeSim[,4001:5000]

#calculate row means all
fold.quasiseq.simseq.fit1.mean = rowMeans(fold.quasiseq.simseq.fit1.estError)
fold.quasiseq.simseq.fit2.mean = rowMeans(fold.quasiseq.simseq.fit2.estError)
fold.quasiseq.simseq.fit3.mean = rowMeans(fold.quasiseq.simseq.fit3.estError)
fold.voom.simseq.mean = rowMeans(fold.voom.simseq.estError)
fold.quasiseq.nb.fit1.mean = rowMeans(fold.quasiseq.nb.fit1.estError)
fold.quasiseq.nb.fit2.mean = rowMeans(fold.quasiseq.nb.fit2.estError)
fold.quasiseq.nb.fit3.mean = rowMeans(fold.quasiseq.nb.fit3.estError)
fold.voom.nb.mean = rowMeans(fold.voom.nb.estError)

#calculate row means ee genes
fold.quasiseq.simseq.fit1.mean.ee = rowMeans(fold.quasiseq.simseq.fit1.estErrorEE)
fold.quasiseq.simseq.fit2.mean.ee = rowMeans(fold.quasiseq.simseq.fit2.estErrorEE)
fold.quasiseq.simseq.fit3.mean.ee = rowMeans(fold.quasiseq.simseq.fit3.estErrorEE)
fold.voom.simseq.mean.ee = rowMeans(fold.voom.simseq.estErrorEE)
fold.quasiseq.nb.fit1.mean.ee = rowMeans(fold.quasiseq.nb.fit1.estErrorEE)
fold.quasiseq.nb.fit2.mean.ee = rowMeans(fold.quasiseq.nb.fit2.estErrorEE)
fold.quasiseq.nb.fit3.mean.ee = rowMeans(fold.quasiseq.nb.fit3.estErrorEE)
fold.voom.nb.mean.ee = rowMeans(fold.voom.nb.estErrorEE)

#calculate row means de genes
fold.quasiseq.simseq.fit1.mean.de = rowMeans(fold.quasiseq.simseq.fit1.estErrorDE)
fold.quasiseq.simseq.fit2.mean.de = rowMeans(fold.quasiseq.simseq.fit2.estErrorDE)
fold.quasiseq.simseq.fit3.mean.de = rowMeans(fold.quasiseq.simseq.fit3.estErrorDE)
fold.voom.simseq.mean.de = rowMeans(fold.voom.simseq.estErrorDE)
fold.quasiseq.nb.fit1.mean.de = rowMeans(fold.quasiseq.nb.fit1.estErrorDE)
fold.quasiseq.nb.fit2.mean.de = rowMeans(fold.quasiseq.nb.fit2.estErrorDE)
fold.quasiseq.nb.fit3.mean.de = rowMeans(fold.quasiseq.nb.fit3.estErrorDE)
fold.voom.nb.mean.de = rowMeans(fold.voom.nb.estErrorDE)

l <- length(fold.quasiseq.simseq.fit1.mean)
num3 <- data.frame(num = c(fold.quasiseq.simseq.fit1.mean, fold.quasiseq.simseq.fit2.mean, fold.quasiseq.simseq.fit3.mean, fold.voom.simseq.mean), 
                   Method = as.factor(c(rep("Bias Fold Tolerance = 1.10", l), rep("Bias Fold Tolerance = 1", l), rep("Bias Fold Tolerance = Inf", l)
							, rep("Voom", l))),
                   Simulation = rep("SimSeq", l*4),
                   SS = rep("Sample Size: 10",l*4))

l <- length(fold.quasiseq.nb.fit1.mean)
num4 <- data.frame(num = c(fold.quasiseq.nb.fit1.mean, fold.quasiseq.nb.fit2.mean, fold.quasiseq.nb.fit3.mean, fold.voom.nb.mean), 
                   Method = as.factor(c(rep("Bias Fold Tolerance = 1.10", l), rep("Bias Fold Tolerance = 1", l), rep("Bias Fold Tolerance = Inf", l)
							, rep("Voom", l))),
                   Simulation = rep("NB", l*4),
                   SS = rep("Sample Size: 10",l*4)) 
				   
l <- length(fold.quasiseq.simseq.fit1.mean.ee)
num3EE <- data.frame(num = c(fold.quasiseq.simseq.fit1.mean.ee, fold.quasiseq.simseq.fit2.mean.ee, fold.quasiseq.simseq.fit3.mean.ee, fold.voom.simseq.mean.ee), 
                   Method = as.factor(c(rep("Bias Fold Tolerance = 1.10", l), rep("Bias Fold Tolerance = 1", l), rep("Bias Fold Tolerance = Inf", l)
							, rep("Voom", l))),
                   Simulation = rep("SimSeq", l*4),
                   SS = rep("Sample Size: 10",l*4))

l <- length(fold.quasiseq.nb.fit1.mean.ee)
num4EE <- data.frame(num = c(fold.quasiseq.nb.fit1.mean.ee, fold.quasiseq.nb.fit2.mean.ee, fold.quasiseq.nb.fit3.mean.ee, fold.voom.nb.mean.ee), 
                   Method = as.factor(c(rep("Bias Fold Tolerance = 1.10", l), rep("Bias Fold Tolerance = 1", l), rep("Bias Fold Tolerance = Inf", l)
							, rep("Voom", l))),
                   Simulation = rep("NB", l*4),
                   SS = rep("Sample Size: 10",l*4)) 
				   
l <- length(fold.quasiseq.simseq.fit1.mean.de)
num3DE <- data.frame(num = c(fold.quasiseq.simseq.fit1.mean.de, fold.quasiseq.simseq.fit2.mean.de, fold.quasiseq.simseq.fit3.mean.de, fold.voom.simseq.mean.de), 
                   Method = as.factor(c(rep("Bias Fold Tolerance = 1.10", l), rep("Bias Fold Tolerance = 1", l), rep("Bias Fold Tolerance = Inf", l)
							, rep("Voom", l))),
                   Simulation = rep("SimSeq", l*4),
                   SS = rep("Sample Size: 10",l*4))

l <- length(fold.quasiseq.nb.fit1.mean.de)
num4DE <- data.frame(num = c(fold.quasiseq.nb.fit1.mean.de, fold.quasiseq.nb.fit2.mean.de, fold.quasiseq.nb.fit3.mean.de, fold.voom.nb.mean.de), 
                   Method = as.factor(c(rep("Bias Fold Tolerance = 1.10", l), rep("Bias Fold Tolerance = 1", l), rep("Bias Fold Tolerance = Inf", l)
							, rep("Voom", l))),
                   Simulation = rep("NB", l*4),
                   SS = rep("Sample Size: 10",l*4)) 

				   
k.ind <- 20
fold.quasiseq.simseq.fit1 <- readRDS(file.path(mainDir, subDir, paste0("ss", k.ind), "fold_quasiseq_simseq_fit1.RDS"))
fold.quasiseq.simseq.fit2 <- readRDS(file.path(mainDir, subDir, paste0("ss", k.ind), "fold_quasiseq_simseq_fit2.RDS"))
fold.quasiseq.simseq.fit3 <- readRDS(file.path(mainDir, subDir, paste0("ss", k.ind), "fold_quasiseq_simseq_fit3.RDS"))
fold.voom.simseq <- readRDS(file.path(mainDir, subDir, paste0("ss", k.ind), "fold_voom_simseq.RDS"))
fold.quasiseq.nb.fit1 <- readRDS(file.path(mainDir, subDir, paste0("ss", k.ind), "fold_quasiseq_nb_fit1.RDS"))
fold.quasiseq.nb.fit2 <- readRDS(file.path(mainDir, subDir, paste0("ss", k.ind), "fold_quasiseq_nb_fit2.RDS"))
fold.quasiseq.nb.fit3 <- readRDS(file.path(mainDir, subDir, paste0("ss", k.ind), "fold_quasiseq_nb_fit3.RDS"))
fold.voom.nb <- readRDS(file.path(mainDir, subDir, paste0("ss", k.ind), "fold_voom_nb.RDS"))
fold.trulyDEgeneList <- readRDS(file.path(mainDir, subDir, paste0("ss", k.ind), "fold_trulyDEgeneList.RDS"))

#finding de genes for each simulation in logLambdaFoldChange
for(ii in 1:200){
	for(jj in 1:1000){
		logLambdaFoldChangeSim[ii,4000+jj] = logLambdaFoldChange[fold.trulyDEgeneList[ii,jj]]
	}
}

#convert to log2
fold.quasiseq.simseq.fit1 <- log2(exp(fold.quasiseq.simseq.fit1*2))
fold.quasiseq.simseq.fit2 <- log2(exp(fold.quasiseq.simseq.fit2*2))
fold.quasiseq.simseq.fit3 <- log2(exp(fold.quasiseq.simseq.fit3*2))
fold.voom.simseq <- fold.voom.simseq*2
fold.quasiseq.nb.fit1 <- log2(exp(fold.quasiseq.nb.fit1*2))
fold.quasiseq.nb.fit2 <- log2(exp(fold.quasiseq.nb.fit2*2))
fold.quasiseq.nb.fit3 <- log2(exp(fold.quasiseq.nb.fit3*2))
fold.voom.nb <- fold.voom.nb*2

#calculate estimate errors
fold.quasiseq.simseq.fit1.estError = fold.quasiseq.simseq.fit1 - logLambdaFoldChangeSim
fold.quasiseq.simseq.fit2.estError = fold.quasiseq.simseq.fit2 - logLambdaFoldChangeSim
fold.quasiseq.simseq.fit3.estError = fold.quasiseq.simseq.fit3 - logLambdaFoldChangeSim
fold.voom.simseq.estError = fold.voom.simseq - logLambdaFoldChangeSim
fold.quasiseq.nb.fit1.estError = fold.quasiseq.nb.fit1 - logLambdaFoldChangeSim
fold.quasiseq.nb.fit2.estError = fold.quasiseq.nb.fit2 - logLambdaFoldChangeSim
fold.quasiseq.nb.fit3.estError = fold.quasiseq.nb.fit3 - logLambdaFoldChangeSim
fold.voom.nb.estError = fold.voom.nb - logLambdaFoldChangeSim

#separate ee genes
fold.quasiseq.simseq.fit1.estErrorEE = fold.quasiseq.simseq.fit1[,1:4000] - logLambdaFoldChangeSim[,1:4000]
fold.quasiseq.simseq.fit2.estErrorEE = fold.quasiseq.simseq.fit2[,1:4000] - logLambdaFoldChangeSim[,1:4000]
fold.quasiseq.simseq.fit3.estErrorEE = fold.quasiseq.simseq.fit3[,1:4000] - logLambdaFoldChangeSim[,1:4000]
fold.voom.simseq.estErrorEE = fold.voom.simseq[,1:4000] - logLambdaFoldChangeSim[,1:4000]
fold.quasiseq.nb.fit1.estErrorEE = fold.quasiseq.nb.fit1[,1:4000] - logLambdaFoldChangeSim[,1:4000]
fold.quasiseq.nb.fit2.estErrorEE = fold.quasiseq.nb.fit2[,1:4000] - logLambdaFoldChangeSim[,1:4000]
fold.quasiseq.nb.fit3.estErrorEE = fold.quasiseq.nb.fit3[,1:4000] - logLambdaFoldChangeSim[,1:4000]
fold.voom.nb.estErrorEE = fold.voom.nb[,1:4000] - logLambdaFoldChangeSim[,1:4000]

#separate de genes
fold.quasiseq.simseq.fit1.estErrorDE = fold.quasiseq.simseq.fit1[,4001:5000] - logLambdaFoldChangeSim[,4001:5000]
fold.quasiseq.simseq.fit2.estErrorDE = fold.quasiseq.simseq.fit2[,4001:5000] - logLambdaFoldChangeSim[,4001:5000]
fold.quasiseq.simseq.fit3.estErrorDE = fold.quasiseq.simseq.fit3[,4001:5000] - logLambdaFoldChangeSim[,4001:5000]
fold.voom.simseq.estErrorDE = fold.voom.simseq[,4001:5000] - logLambdaFoldChangeSim[,4001:5000]
fold.quasiseq.nb.fit1.estErrorDE = fold.quasiseq.nb.fit1[,4001:5000] - logLambdaFoldChangeSim[,4001:5000]
fold.quasiseq.nb.fit2.estErrorDE = fold.quasiseq.nb.fit2[,4001:5000] - logLambdaFoldChangeSim[,4001:5000]
fold.quasiseq.nb.fit3.estErrorDE = fold.quasiseq.nb.fit3[,4001:5000] - logLambdaFoldChangeSim[,4001:5000]
fold.voom.nb.estErrorDE = fold.voom.nb[,4001:5000] - logLambdaFoldChangeSim[,4001:5000]

#calculate row means all
fold.quasiseq.simseq.fit1.mean = rowMeans(fold.quasiseq.simseq.fit1.estError)
fold.quasiseq.simseq.fit2.mean = rowMeans(fold.quasiseq.simseq.fit2.estError)
fold.quasiseq.simseq.fit3.mean = rowMeans(fold.quasiseq.simseq.fit3.estError)
fold.voom.simseq.mean = rowMeans(fold.voom.simseq.estError)
fold.quasiseq.nb.fit1.mean = rowMeans(fold.quasiseq.nb.fit1.estError)
fold.quasiseq.nb.fit2.mean = rowMeans(fold.quasiseq.nb.fit2.estError)
fold.quasiseq.nb.fit3.mean = rowMeans(fold.quasiseq.nb.fit3.estError)
fold.voom.nb.mean = rowMeans(fold.voom.nb.estError)

#calculate row means ee genes
fold.quasiseq.simseq.fit1.mean.ee = rowMeans(fold.quasiseq.simseq.fit1.estErrorEE)
fold.quasiseq.simseq.fit2.mean.ee = rowMeans(fold.quasiseq.simseq.fit2.estErrorEE)
fold.quasiseq.simseq.fit3.mean.ee = rowMeans(fold.quasiseq.simseq.fit3.estErrorEE)
fold.voom.simseq.mean.ee = rowMeans(fold.voom.simseq.estErrorEE)
fold.quasiseq.nb.fit1.mean.ee = rowMeans(fold.quasiseq.nb.fit1.estErrorEE)
fold.quasiseq.nb.fit2.mean.ee = rowMeans(fold.quasiseq.nb.fit2.estErrorEE)
fold.quasiseq.nb.fit3.mean.ee = rowMeans(fold.quasiseq.nb.fit3.estErrorEE)
fold.voom.nb.mean.ee = rowMeans(fold.voom.nb.estErrorEE)

#calculate row means de genes
fold.quasiseq.simseq.fit1.mean.de = rowMeans(fold.quasiseq.simseq.fit1.estErrorDE)
fold.quasiseq.simseq.fit2.mean.de = rowMeans(fold.quasiseq.simseq.fit2.estErrorDE)
fold.quasiseq.simseq.fit3.mean.de = rowMeans(fold.quasiseq.simseq.fit3.estErrorDE)
fold.voom.simseq.mean.de = rowMeans(fold.voom.simseq.estErrorDE)
fold.quasiseq.nb.fit1.mean.de = rowMeans(fold.quasiseq.nb.fit1.estErrorDE)
fold.quasiseq.nb.fit2.mean.de = rowMeans(fold.quasiseq.nb.fit2.estErrorDE)
fold.quasiseq.nb.fit3.mean.de = rowMeans(fold.quasiseq.nb.fit3.estErrorDE)
fold.voom.nb.mean.de = rowMeans(fold.voom.nb.estErrorDE)

l <- length(fold.quasiseq.simseq.fit1.mean)
num5 <- data.frame(num = c(fold.quasiseq.simseq.fit1.mean, fold.quasiseq.simseq.fit2.mean, fold.quasiseq.simseq.fit3.mean, fold.voom.simseq.mean), 
                   Method = as.factor(c(rep("Bias Fold Tolerance = 1.10", l), rep("Bias Fold Tolerance = 1", l), rep("Bias Fold Tolerance = Inf", l)
							, rep("Voom", l))),
                   Simulation = rep("SimSeq", l*4),
                   SS = rep("Sample Size: 20",l*4))

l <- length(fold.quasiseq.nb.fit1.mean)
num6 <- data.frame(num = c(fold.quasiseq.nb.fit1.mean, fold.quasiseq.nb.fit2.mean, fold.quasiseq.nb.fit3.mean, fold.voom.nb.mean), 
                   Method = as.factor(c(rep("Bias Fold Tolerance = 1.10", l), rep("Bias Fold Tolerance = 1", l), rep("Bias Fold Tolerance = Inf", l)
							, rep("Voom", l))),
                   Simulation = rep("NB", l*4),
                   SS = rep("Sample Size: 20",l*4)) 
				   
l <- length(fold.quasiseq.simseq.fit1.mean.ee)
num5EE <- data.frame(num = c(fold.quasiseq.simseq.fit1.mean.ee, fold.quasiseq.simseq.fit2.mean.ee, fold.quasiseq.simseq.fit3.mean.ee, fold.voom.simseq.mean.ee), 
                   Method = as.factor(c(rep("Bias Fold Tolerance = 1.10", l), rep("Bias Fold Tolerance = 1", l), rep("Bias Fold Tolerance = Inf", l)
							, rep("Voom", l))),
                   Simulation = rep("SimSeq", l*4),
                   SS = rep("Sample Size: 20",l*4))

l <- length(fold.quasiseq.nb.fit1.mean.ee)
num6EE <- data.frame(num = c(fold.quasiseq.nb.fit1.mean.ee, fold.quasiseq.nb.fit2.mean.ee, fold.quasiseq.nb.fit3.mean.ee, fold.voom.nb.mean.ee), 
                   Method = as.factor(c(rep("Bias Fold Tolerance = 1.10", l), rep("Bias Fold Tolerance = 1", l), rep("Bias Fold Tolerance = Inf", l)
							, rep("Voom", l))),
                   Simulation = rep("NB", l*4),
                   SS = rep("Sample Size: 20",l*4)) 
				   
l <- length(fold.quasiseq.simseq.fit1.mean.de)
num5DE <- data.frame(num = c(fold.quasiseq.simseq.fit1.mean.de, fold.quasiseq.simseq.fit2.mean.de, fold.quasiseq.simseq.fit3.mean.de, fold.voom.simseq.mean.de), 
                   Method = as.factor(c(rep("Bias Fold Tolerance = 1.10", l), rep("Bias Fold Tolerance = 1", l), rep("Bias Fold Tolerance = Inf", l)
							, rep("Voom", l))),
                   Simulation = rep("SimSeq", l*4),
                   SS = rep("Sample Size: 20",l*4))

l <- length(fold.quasiseq.nb.fit1.mean.de)
num6DE <- data.frame(num = c(fold.quasiseq.nb.fit1.mean.de, fold.quasiseq.nb.fit2.mean.de, fold.quasiseq.nb.fit3.mean.de, fold.voom.nb.mean.de), 
                   Method = as.factor(c(rep("Bias Fold Tolerance = 1.10", l), rep("Bias Fold Tolerance = 1", l), rep("Bias Fold Tolerance = Inf", l)
							, rep("Voom", l))),
                   Simulation = rep("NB", l*4),
                   SS = rep("Sample Size: 20",l*4))  

num <- rbind(num1, num2, num3, num4, num5, num6)
qplot(Simulation, num, data = num, geom = "boxplot", facets = SS ~ Method, 
      ylab = "Across All Genes") + theme_bw()+ 
  theme(axis.text.x = element_text(angle=90))

num <- rbind(num1EE, num2EE, num3EE, num4EE, num5EE, num6EE)
qplot(Simulation, num, data = num, geom = "boxplot", facets = SS ~ Method, 
      ylab = "EE Genes Only") + theme_bw()+ 
  theme(axis.text.x = element_text(angle=90))

num <- rbind(num1DE, num2DE, num3DE, num4DE, num5DE, num6DE)
qplot(Simulation, num, data = num, geom = "boxplot", facets = SS ~ Method, 
      ylab = "DE Genes Only") + theme_bw()+ 
  theme(axis.text.x = element_text(angle=90))  
save.image('Fold_Change_Plot.RData')