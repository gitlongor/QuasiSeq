if(FALSE){ ## only for reference
	mockRNASeqData.expression=expression({
		set.seed(9992722L, kind="Mersenne-Twister")
		
		n.genes<-10000L
		n.de<-round(.35 * n.genes)
		trt<-gl(2L,4L)
		design = model.matrix(~trt)
		n.samp<-length(trt)
		mu<-rgamma(n.genes, 1.5, .01)

		minTotalCount = 0
		## specify gene specific negative binomial dispersions
		size<-(log(mu+exp(1))-1)/mu  ### Var(Y)=E(Y)log(E(Y)+exp(1))

		## add noise to gene specific negative binomial dispersions
		size<-size*4.5/rchisq(n.genes,4.5)

		sim.mn<-matrix(mu,n.genes,2)

		### Simulate fold changes
		B<-exp((2*rbinom(n.de,1,.5)-1)*(.25+rbeta(n.de,1,2)))

		sim.mn[1:n.de,1]<-sim.mn[1:n.de,1]*B^(.5) ## there was an extra  +5 for this line
		sim.mn[1:n.de,2]<-sim.mn[1:n.de,2]*B^(-.5)

		### Simulate library size factors
		sim.offset <- 2^(rnorm(n.samp,0,.15))
		sim.offset = exp(scale(log(sim.offset), center=TRUE, scale=FALSE))
		attributes(sim.offset)=NULL

		### Compute final means
		sim.mn2<-(sim.mn[,trt])*sim.offset[rep(seq(n.samp), each=n.genes)]

		### Simulate data
		simdat<-matrix(rnbinom(n.samp*n.genes,mu=sim.mn2,size=1/size),n.genes,n.samp)
		rowsToKeep = which(rowSums(simdat)>minTotalCount)
		ngenes = length(rowsToKeep)

		## estimated log.offsets
		est.offset <- edgeR::calcNormFactors(simdat, method='TMM')

		## estimated nb.disp
		d <- edgeR::DGEList(counts = simdat, group = trt, norm.factors = est.offset)
		nb.disp <- edgeR::estimateGLMTrendedDisp(d, design)$trended.dispersion

		list(
			counts = simdat[rowsToKeep,],
			treatment = trt, 
			design.matrix = design, 
			true.normalization = sim.offset, 
			estimated.normalization = est.offset, 
			true.nbdisp = size, 
			estimated.nbdisp = nb.disp, 
			ngenes = ngenes, 
			nsamples = n.samp,
			true.DEgenes = which(rowsToKeep <= n.de), 
			true.foldChanges = B[rowsToKeep <= n.de],
			simulation.expression = NULL
		)
	})
	mockRNASeqData = eval(mockRNASeqData.expression)
	mockRNASeqData$simulation.expression = mockRNASeqData.expression
	save(mockRNASeqData, file='data/mockRNASeqData.RData',compress='xz', compression_level=9L)
}
