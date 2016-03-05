PoisDev<-function(counts,design,log.offset,print.progress=TRUE)
{
	n<-NCOL(counts)
	if(missing(log.offset)) log.offset=rep(0, n)
	if(!is.matrix(counts)) counts=as.matrix(counts)

	### For one-way equivalent designs, MLE's for Poisson GLM have simple closed form, so we can avoid one-gene-at-a-time GLM fitting.
	### This part was rewritten by Long Qu since version 1.0-9. 
	### Even if "design" is not a vector 
	if(is.vector(design) || (is.factor(design) && is.null(dim(design))) )
		design = model.matrix(~as.factor(design))
	stopifnot(is.matrix(design) && is.numeric(design))
	rk=qr(design)$rank
	uniqx=unique(design)
	uniqxTInv = if(rk==2) solve22(t(uniqx)) else solve(t(uniqx))
	if(NROW(uniqx) == rk) {  
		group = grpDuplicated(design)
		oneWayX=1 * (group==rep(seq_len(rk), each=n)); dim(oneWayX)=c(n, rk)
		## oneWayX is an equivalent design matrix in 0-1 coding, w/o intercept
		oneWayN = colSums(oneWayX)
		
		##offset should NOT be given on log scale here
		offset <- exp(log.offset)
		
		sum.counts = counts %*% oneWayX
		sum.scales = drop(crossprod(oneWayX, offset))
		fitted.expXb = sweep(sum.counts, 2L, sum.scales, "/")
		
		means <- fitted.expXb[ , match(group, seq_len(rk)),drop=FALSE]
		means = sweep(means, 2, offset, "*")
		linear.predictors = log(means) ## for use in bartlett; might be -Inf; but -Inf won't be used when computing bartlett factor, since they won't have good.weights

		parms = pmax(-.5*.Machine$double.xmax,log(fitted.expXb)) ## eta under new parameterization
		dim(parms)=dim(fitted.expXb)
		parms = parms %*% uniqxTInv ## beta under original parameterization
		colnames(parms) = colnames(design)
		
		### 0's require special attention since 0^0=1, but 0*log(0)=NaN
		deviance<-means-counts
		deviance[counts!=0]<-deviance[counts!=0]+(counts*log(counts/means))[counts!=0]
		deviance<-2*rowSums(deviance)
		
		BartEpsilons=rep(NA_real_, NROW(counts))
		glmFit = list(
			prior.weights = rep(1, n), 
			family=fix.family.link(fix.family.var(poisson('log')))
		)
		if(print.progress) progressBar = txtProgressBar(max = NROW(counts),style=3L)
		for(i in seq_len(nrow(counts))){
			#if(i%in%c(2,10,100,500,1000,2500,4000,5000*(1:200))&print.progress) print(paste("Analyzing Gene #",i))
			if(print.progress) setTxtProgressBar(progressBar, value=i)
			glmFit$linear.predictors=linear.predictors[i,]
			glmFit$weights=means[i,]
			glmFit$qr=qr(sqrt(means[i,])*design)
			BartEpsilons[i]=bartlettFactor(glmFit)
		}
		if(print.progress) close(progressBar)
	}else{
		### For general designs, the first column of each element (matrix) of design.list should be a column of 1's, pertaining to the intercept. 
		## Note by LQ: There seems no reason to require an intercept... Changed model formula below in glm.
		
		deviance<-rep(NA_real_,nrow(counts)); 
		means<-matrix(NA_real_,nrow(counts),ncol(counts)); 
		parms<-matrix(NA_real_,nrow(counts),ncol(design))
		BartEpsilons = deviance ## NA's
		
		### For each gene and given design matrix, fit GLM to find model parameters (for mean structure) that optimize quasi-likelihood
		if(print.progress) progressBar = txtProgressBar(max = NROW(counts),style=3L)
		for(i in 1:nrow(counts)){
			### If wanted, provide running progress update (eventually once every 5000 genes) 
			#if(i%in%c(2,10,100,500,1000,2500,4000,5000*(1:200))&print.progress) print(paste("Analyzing Gene #",i))
			if(print.progress) setTxtProgressBar(progressBar, value=i)
			
			### Fit GLM
			res<-withCallingHandlers(
				glm(counts[i,]~design-1,family="poisson",offset=log.offset,method=glm.fit3)	##offset should be given on log scale here
				, simpleWarning=ignorableWarnings
			)
			
			### Save optimized means (used in Pearson's dispersion estimator)
			means[i,]<-res$fitted.values
			parms[i,]<-res$coefficients
			BartEpsilons[i]=bartlettFactor(res)

			### Save deviance (used to compute LRT for comparing models and also deviance dispersion estimator)
			deviance[i]<-res$deviance
		}
		if(print.progress) close(progressBar)
	}

return(list(dev=deviance,means=means,parms=parms, bartlett.epsilon = BartEpsilons))
}