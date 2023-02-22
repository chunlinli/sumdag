#estimate key quantities from summary statistics
#give an option to enter a genotype txt file to calculate LD: an example is "chr22.26SNPs.ref.3000_3.txt"
#N is a matrix correspond to the sample size for each SNP-trait pair, it can also be a number
summary_stats_preprocess <- function( fullbetas, fullses, fullP, N, xtxdiag=0, ref=T, geno.ref.txt=NULL){
	
	if (length(N)==1) N <- matrix(rep(N, p*q), nrow=q, ncol=p)
		
	#read in a txt file and calculate its LD matrix
	if (ref) {geno = read.table(geno.ref.txt)	
	n.ref = nrow(geno)

	#calculate LD
	genoMat.std <- apply(geno,2,function(x) scale(x,scale=T, center=T))
	tempMat = t(genoMat.std) %*% (genoMat.std)/n.ref
	#calculate var of SNP j
	sds = diag(tempMat)
    rm(geno)
    #if no genotype is provided, then we assume all SNPs are pruned to be independent and used a diagonal matrix for LD
    } else {
    	#for standardized genotype
    	sds = rep(1,nrow(fullbetas))
    	tempMat = diag(sds)
    }

	#yty/N
	ytyVec <- NULL
	for (j in 1:ncol(fullbetas)){
		tempEst <- vector()
	for (i in 1:nrow(fullbetas)) {
		ses = fullses[,j]
		betas = fullbetas[,j]
        tempEst[i] = N[i,j] * sds[i] * ses[i]^2 + sds[i] * betas[i]^2
    }
    ytyEst = median(tempEst)
    ytyVec = c(ytyVec, ytyEst)
    }
    
    #calculate xty/N
    xtyMat <- NULL
	for (j in 1:ncol(fullbetas)){
	betas = fullbetas[,j]
    xtyEst = sds[j] * betas
    xtyMat = cbind(xtyMat, xtyEst)
    }
    
    #xtx/n
    xtx = tempMat + diag(rep(xtxdiag, nrow(tempMat)))
    		
	return(list(xtx=xtx, ytyVec=ytyVec, xtyMat=xtyMat, sds=sds))
	}  
	
#estimating V matrix based on summary statistic; y is standardized in the TLPsum
estimateV_sum <- function(ss, N, ss_back){
	ss$fullP[ss$fullP==0] = 1e-100
	
	p <- ncol(ss$fullP)
	q <- nrow(ss$fullP)
	if (length(N)==1) N <- matrix(rep(N, p*q), nrow=q, ncol=p)

  	V <- matrix(0, nrow=q, ncol = p)

	for (i in 1:p){
	#standardize y
	xty <- ss_back$xtyMat[,i]/sqrt(ss_back$ytyVec[i])
	
    x1 <- glmtlp::sumtlp(XX=ss_back$xtx, Xy=xty, method='tlp-c',kappa=1:ncol(ss_back$xtx), nobs=median(N[,i]))
    
    # Extracting result for pseudo AIC/BIC function
    penBetas <- x1$beta    
    uniBetas <- ss$fullBetas[,i]
    uniSD <- ss$fullBetaSD[,i]
    
    # Obtaining psuedo AIC/BIC values
    x2 <- pseudoBIC.internal(penalizedBetas = penBetas , betas = uniBetas, ses = uniSD, Ni =  N[,i], sds=ss_back$sds, xtx=ss_back$xtx, ytyEst=1, xtyEst=xty)
                
    # Saving V from minimum BIC
    V[,i] <- penBetas[,which.min(x2$bic)] # * sqrt(ss_back$ytyVec[i])    
  }
  return(Vest = V)
}

#pesudoBIC
#N here is a vector of sample size 
pseudoBIC.internal <- function(penalizedBetas, betas, ses, Ni, sds, xtx, ytyEst, xtyEst){
    SSEvec = NULL
    qVec = NULL
    aicVec = NULL
    bicVec = NULL
    bxxb = NULL
    bxy = NULL
    SSEvec = NULL
    
    #because it is low dimension, so all the SNPs are used for estimating the sigma^2
    xtxinv = solve(xtx)
    
    #obtain a consistent estimator for sigma^2 tilde
    qTemp = nrow(xtx)
    sigSqTilde = (ytyEst - t(xtyEst) %*% xtxinv %*% xtyEst) * median(Ni)/(median(Ni) - qTemp)
         
    #calculate SSE 
    for (k in 1:ncol(penalizedBetas)) {
      	penalizedBetasTemp = penalizedBetas[, k]        
        
        #q here is the number of parameters not equal to 0
        qVec = sum(penalizedBetasTemp != 0)
                
        #bxxb for SSE (correct a mistake here on Oct 23)
        bxxb = t(penalizedBetasTemp) %*% xtx %*% penalizedBetasTemp
        #bxy for SSE
        bxyTemp = t(penalizedBetasTemp) %*% xtyEst
        
        #sse
        SSEest = (ytyEst - (2 * bxyTemp) + bxxb)* median(Ni)
        SSEvec = c(SSEvec, SSEest[1, 1])
        logLik = -SSEest /(2 * sigSqTilde)
        
        aicTemp = (2 * qVec) - (2 * logLik)
        bicTemp = (log(median(Ni)) * qVec) - (2 * logLik) 
        aicVec = c(aicVec, aicTemp[1, 1])
        bicVec = c(bicVec, bicTemp[1, 1])
    }
    
    toReturn = structure(list(aic = aicVec, bic = bicVec, SSE = SSEvec,
        q = qVec, logLik = logLik, bxxb = bxxb,  sigSqTilde = sigSqTilde))
    return(toReturn)
}

peeling <- function(V, thresh){
	source('./topological_order.R')
	peel.res <- topological_order(v=V, thresh = thresh)
  	PiHat <- as.matrix(peel.res$an_mat) #ancestry
  	col.names(PiHat) <- colnames(V)
  	row.names(PiHat) <- colnames(V)
  	PhiHat <- as.matrix(peel.res$iv_mat) #IV
  	col.names(PhiHat) <- colnames(V)
  	row.names(PhiHat) <- rownames(V) 
  	
  	#calculate the effective number of tests
  	 ntest = 0
     for(j in 1:nrow(PiHat)){
        ancestor <- which(PiHat[, j] != 0)
        intervention <- which(PhiHat[, j] != 0)
		ntest <- ntest + sum(length(ancestor)!=0 )
        }
    p <- ncol(PiHat)
    effect.test = p^2 - p - ntest    
        
 	return(list(PiHat=PiHat, PhiHat=PhiHat, effect.test=effect.test))
}

#TLP after peeling
TLP_afterpeel_U <- function(PiHat, PhiHat, ss_back, cormat, N){
	
	PiHat.new <- matrix(0, nrow=nrow(PiHat) , ncol=ncol(PiHat)) #15x15
	PhiHat.new <- matrix(0, nrow=nrow(PhiHat) , ncol=ncol(PhiHat)) 
	#covmat = diag(sqrt(output_prepare$ytyVec)) %*% cormat %*% diag(sqrt(output_prepare$ytyVec))
	covmat=cormat

	for (j in 1:ncol(PiHat)){
	PiHatj = PiHat[,j]
	PhiHatj = PhiHat[,j]
	Jp <- which(PhiHatj != 0) #IV set
	Ap <- which(PiHatj != 0) #ancestry set
	#standardize y
	xtyMat <- t(t(ss_back$xtyMat)/sqrt(ss_back$ytyVec))
	if (length(Jp) > 0) {
		XjpXjp = ss_back$xtx[Jp,Jp, drop=F]
		XjpYp = xtyMat[Jp, j]
		Nvec = N[Jp,j]
		#when there is ancestor
		if (length(Ap) > 0) {
	XjpYap = xtyMat[Jp,Ap, drop=F]
	YapYap = covmat[Ap, Ap, drop=F]
	YapYp = covmat[Ap, j]
	XTX = rbind(cbind(XjpXjp, XjpYap), cbind(t(XjpYap),YapYap))
	XTY = c(XjpYp, YapYp)
	x1 <- glmtlp::sumtlp(XX=XTX, Xy=XTY, method='tlp-c',kappa=1:ncol(XTX), nobs=median(Nvec))
	
	penBetas <- x1$beta        
    # Obtaining psuedo AIC/BIC values
    x2 <- pseudoBIC.internal1(penalizedBetas = penBetas , N = median(Nvec),  xtx=XTX, ytyEst=1, xtyEst=XTY)
       
    PhiHat.new[Jp, j] = penBetas[1:nrow(XjpXjp),which.min(x2$bic)]
    PiHat.new[Ap, j] = penBetas[(nrow(XjpXjp)+1):ncol(XTX), which.min(x2$bic)]
  	} 
	} else if (length(Jp)==0 & length(Ap) > 0){
		YapYap = covmat[Ap, Ap, drop=F]
		YapYp = covmat[Ap, j]
		BETAS = solve(YapYap) %*% YapYp
		x1 <- glmtlp::sumtlp(XX=YapYap, Xy=YapYp, method='tlp-c',kappa=1:ncol(YapYap), nobs=median(Nvec))
		penBetas <- x1$beta        
    # Obtaining psuedo AIC/BIC values
    	x2 <- pseudoBIC.internal1(penalizedBetas = penBetas , N = median(Nvec),  xtx=YapYap, ytyEst=1, xtyEst=XTY)
       
    PhiHat.new[Jp, j] = penBetas[,which.min(x2$bic)]
    PiHat.new[Ap, j] = 0		
	}
}
	return(list(W.est=PhiHat.new, U.est=PiHat.new))
	}

