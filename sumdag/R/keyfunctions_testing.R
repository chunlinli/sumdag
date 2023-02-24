#testing

#estimating Sigma using summary statistics in the regression of Y ~ X(IV) + Y(ancestor) on j's y


#' Title
#'
#' @param output_prepare
#' @param PhiHatj
#' @param PiHatj
#' @param covmat
#' @param N
#' @param j
#'
#' @return
#' @export
#'
#' @examples
estimate_Sigma <- function(output_prepare, PhiHatj, PiHatj, covmat, N, j){

		Jp <- which(PhiHatj != 0) #IV set
		Ap <- which(PiHatj != 0) #ancestry set
		#standardize y
	    xtyMat <- t(t(ss_back$xtyMat)/sqrt(ss_back$ytyVec))
		if (length(Jp) > 0 & length(Ap) > 0){
		XjpXjp = output_prepare$xtx[Jp,Jp, drop=F]
		XjpYp = xtyMat[Jp, j]
		XjpYap = xtyMat[Jp,Ap, drop=F]
		YapYap = covmat[Ap, Ap, drop=F]
		YapYp = covmat[Ap, j]
		INV = solve(rbind(cbind(XjpXjp, XjpYap), cbind(t(XjpYap),YapYap)))
		BETAS = INV %*% c(XjpYp, YapYp)
		Nvec = N[Jp,j]
		RSS = (1 - t(c(XjpYp, YapYp))%*%BETAS)*median(Nvec)
		sigma2p = RSS/(median(Nvec)-length(Ap)-length(Jp) + 1)
    	} else if (length(Jp) == 0 & length(Ap) > 0){
		YapYap = covmat[Ap, Ap, drop=F]
		YapYp = covmat[Ap, j]
		BETAS = solve(YapYap) %*% YapYp
		Nvec = N[,j]
		RSS = (1 - t(YapYp)%*%BETAS)*median(Nvec)
		sigma2p = RSS/(median(Nvec) - length(Ap) + 1)
    	} else if (length(Jp) > 0 & length(Ap)==0){
    		Nvec = N[Jp,j]
    		XjpXjp = output_prepare$xtx[Jp,Jp, drop=F]
		XjpYp = xtyMat[Jp, j]
		BETAS = solve(XjpXjp) %*% XjpYp
		RSS = (1 - t(XjpYp)%*%BETAS)*median(Nvec)
		sigma2p = RSS/(median(Nvec) - length(Jp) + 1)
    	} else if (length(Jp)==0 & length(Ap)==0){
    		Nvec = N[,j]
    		RSS =  median(Nvec)
    		sigma2p = RSS/(median(Nvec) - 1)
    	}
    	return(list(sigma2hat=sigma2p, RSS=RSS))
}

#specify H and check acyclicity constraint
#f_mat is a U matrix specifying what to detect
#' Title
#'
#' @param PiHat
#' @param testj
#' @param j
#'
#' @return
#' @export
#'
#' @examples
check_regularity <- function(PiHat, testj, j){
	  	kk <- setdiff(which(PiHatj != 0), testj)
	  	if (length(kk) > 0){
	  		for ( i in 1:length(kk)) {
	  			if (PiHat[j,kk[i]]!=0) {paste0('The path from ', kk[i], ' to ', j, 'is not testable \n')
	  				testj[kk[i]] = PiHatj[kk[i]]
	  				}
	  		}
	  		}
	  return(testj)
	  }

#' Title
#'
#' @param ss_back
#' @param PiHat
#' @param PhiHat
#' @param pairs
#' @param covmat
#' @param N
#' @param method
#' @param test_type
#'
#' @return
#' @export
#'
#' @examples
causal_inference <- function(ss_back, PiHat, PhiHat, pairs, covmat, N,
                             method = "asymptotic", test_type = c("edge", "path")) {
    # based on an_mat, in_mat, f_mat, generate causal inference p-value
    if (test_type == "edge" && method == "asymptotic") {
        likelihood_ratios <- asymptotic_inference_internal(ss_back, PiHat, PhiHat, pairs, covmat, N)
        likelihood_ratio <- sum(likelihood_ratios$likelihood_ratios)
        df <- sum(likelihood_ratios$df)
        p_value <- pchisq(likelihood_ratio, df = df, lower.tail = FALSE)
        list(likelihood_ratio = likelihood_ratio, df = df, p_value = p_value)
    }  else if (test_type == "path" && method == "asymptotic") {
        likelihood_ratios <- asymptotic_inference_internal(ss_back, PiHat, PhiHat, pairs, covmat, N)
        likelihood_ratio <- likelihood_ratios$likelihood_ratios
		#if (!all(likelihood_ratio==0))  likelihood_ratio <- likelihood_ratio[which(likelihood_ratio!=0)]
        p_value <- max(pchisq(likelihood_ratio, df = 1, lower.tail = FALSE))
        list(likelihood_ratios = likelihood_ratios, p_value = p_value)

    } else {
        stop("Invalid test_type or method.")
    }
}


#' Title
#'
#' @param ss_back
#' @param PiHat
#' @param PhiHat
#' @param pairs
#' @param covmat
#' @param N
#'
#' @return
#' @export
#'
#' @examples
asymptotic_inference_internal <- function(ss_back, PiHat, PhiHat, pairs, covmat, N) {

    # based on an_mat, in_mat, f_mat, generate likelihood ratios
    d_mat <- PiHat
    alt_mat <- PiHat
    jj <- NULL
    for (j in 1:length(pairs)){
    	    loc <- pairs[[j]]
    	    if (PiHat[loc[2], loc[1]] == 0) {
    		d_mat[loc[1],loc[2]] <- 0
      	alt_mat[loc[1],loc[2]] <- 1
    		jj <- c(jj, loc[2])
    		} else {print(paste('Edge from', loc[1],'to', loc[2], 'is degenerate'))}
    }
    jj <- unique(jj)

    likelihood_ratios <- NULL
    df <- NULL
    p <- ncol(PiHat)

	if (length(jj) > 0){
    for (jjj in 1:length(jj)) {
    	#print(j)
    		j <- jj[jjj]
   	 	#PiHatj <- PiHat[,j]
        PhiHatj <- PhiHat[,j]
        testj <- d_mat[,j]
        testaltj <- alt_mat[,j]
        ancestor <- setdiff(which(testaltj!= 0), which(testj!=0))
        #intervention <- which(PhiHat[, j] != 0)
        dfj <- 0

        if ( median(N[,j]) >  sum(testaltj!=0) + sum(PhiHatj!=0) + 1 ) {
                if (length(ancestor) > 0)  {
                an_mat <- testaltj
                iv_mat <- PhiHatj
                sigma_h1 <- estimate_Sigma(ss_back, iv_mat, an_mat, covmat, N, j)

                an_mat <- testj
                sigma_h0 <- estimate_Sigma(ss_back, iv_mat, an_mat, covmat, N, j)

				likelihood_ratio <- (sigma_h0$RSS - sigma_h1$RSS) / sigma_h1$sigma2hat
				likelihood_ratios <- c(likelihood_ratios, likelihood_ratio)
				dfj <- dfj + 1
				df <- c(df, dfj)
				#print(likelihood_ratio)
				#print(j)
        } else likelihood_ratios <- c(likelihood_ratios, 0)
        }
    }
    }

    if (is.null(likelihood_ratios)) {likelihood_ratios = 0   ; df=0}
    return(list(likelihood_ratios=likelihood_ratios, df=df))
}


