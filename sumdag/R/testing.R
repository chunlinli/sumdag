

#' Internal function: estimating Sigma using summary statistics in the regression of Y ~ X(IV) + Y(ancestor) on j's y
#'
#' @param output_prepare output_prepare
#' @param PhiHatj PhiHatj
#' @param PiHatj PiHatj
#' @param cormat cormat
#' @param N N
#' @param j j
#' @param ss_back ss_back
#'
#' @import stats
#' @examples
#' \dontrun{
#'
#' }
estimate_Sigma <- function(output_prepare, PhiHatj, PiHatj, cormat, N, j, ss_back) {
    Jp <- which(PhiHatj != 0) # IV set
    Ap <- which(PiHatj != 0) # ancestry set
    # standardize y
    xtyMat <- t(t(ss_back$xtyMat) / sqrt(ss_back$ytyVec))
    if (length(Jp) > 0 & length(Ap) > 0) {
        XjpXjp <- output_prepare$xtx[Jp, Jp, drop = F]
        XjpYp <- xtyMat[Jp, j]
        XjpYap <- xtyMat[Jp, Ap, drop = F]
        YapYap <- cormat[Ap, Ap, drop = F]
        YapYp <- cormat[Ap, j]
        INV <- solve(rbind(cbind(XjpXjp, XjpYap), cbind(t(XjpYap), YapYap)))
        BETAS <- INV %*% c(XjpYp, YapYp)
        Nvec <- N[Jp, j]
        RSS <- (1 - t(c(XjpYp, YapYp)) %*% BETAS) * median(Nvec)
        sigma2p <- RSS / (median(Nvec) - length(Ap) - length(Jp) + 1)
    } else if (length(Jp) == 0 & length(Ap) > 0) {
        YapYap <- cormat[Ap, Ap, drop = F]
        YapYp <- cormat[Ap, j]
        BETAS <- solve(YapYap) %*% YapYp
        Nvec <- N[, j]
        RSS <- (1 - t(YapYp) %*% BETAS) * median(Nvec)
        sigma2p <- RSS / (median(Nvec) - length(Ap) + 1)
    } else if (length(Jp) > 0 & length(Ap) == 0) {
        Nvec <- N[Jp, j]
        XjpXjp <- output_prepare$xtx[Jp, Jp, drop = F]
        XjpYp <- xtyMat[Jp, j]
        BETAS <- solve(XjpXjp) %*% XjpYp
        RSS <- (1 - t(XjpYp) %*% BETAS) * median(Nvec)
        sigma2p <- RSS / (median(Nvec) - length(Jp) + 1)
    } else if (length(Jp) == 0 & length(Ap) == 0) {
        Nvec <- N[, j]
        RSS <- median(Nvec)
        sigma2p <- RSS / (median(Nvec) - 1)
    }
    return(list(sigma2hat = sigma2p, RSS = RSS))
}

#' Likelihood ratio tests
#'
#' @param PiHat  A pxp superset of the directed relations among phenotypes; output from superDAG function
#' @param PhiHat A qxp superset of the interventional relation; output from superDAG function
#' @param ss_back output from the summary_stats_preprocess function
#' @param pairs A list of vectors of length two to be tested, with each vector represents a directed node from the first to the second item of the vector
#' @param cormat pxp correlation matrix of phenotypes; this matrix can be estimated by empricial correlation of null SNP.
#' @param N A qxp matrix of corresponding sample size for each SNP-phenotype pair. A scalar is also acceptable when all SNP-phentype pairs have the same sample size
#' @param test_type "edge" refers to the test of multiple directed relations; "path" refers to the test of a directed pathway
#'
#' @import stats
#' @return A list of likelihood ratio, and p-value
#' @export
#'
#' @examples
#' \dontrun{
#'
#' }
causal_inference <- function(PiHat, PhiHat, ss_back, pairs, cormat, N,
                             test_type = c("edge", "path")) {
    # based on an_mat, in_mat, f_mat, generate causal inference p-value
    if (test_type == "edge") {
        likelihood_ratios <- asymptotic_inference_internal(ss_back, PiHat, PhiHat, pairs, cormat, N)
        likelihood_ratio <- sum(likelihood_ratios$likelihood_ratios)
        df <- sum(likelihood_ratios$df)
        p_value <- pchisq(likelihood_ratio, df = df, lower.tail = FALSE)
        list(likelihood_ratio = likelihood_ratio, df = df, p_value = p_value)
    } else if (test_type == "path") {
        likelihood_ratios <- asymptotic_inference_internal(ss_back, PiHat, PhiHat, pairs, cormat, N)
        likelihood_ratio <- likelihood_ratios$likelihood_ratios
        p_value <- max(pchisq(likelihood_ratio, df = 1, lower.tail = FALSE))
        list(likelihood_ratios = likelihood_ratios, p_value = p_value)
    } else {
        stop("Invalid test_type or method.")
    }
}

#' Internal function: asymptotic inference
#'
#' @param ss_back ss_back
#' @param PiHat PiHat
#' @param PhiHat PhiHat
#' @param pairs pairs
#' @param cormat cormat
#' @param N N
#'
#' @import stats
#' @examples
#' \dontrun{
#'
#' }
asymptotic_inference_internal <- function(ss_back, PiHat, PhiHat, pairs, cormat, N) {

    # based on an_mat, in_mat, f_mat, generate likelihood ratios
    d_mat <- PiHat
    alt_mat <- PiHat
    jj <- NULL
    for (j in 1:length(pairs)) {
        loc <- pairs[[j]]
        if (PiHat[loc[2], loc[1]] == 0) {
            d_mat[loc[1], loc[2]] <- 0
            alt_mat[loc[1], loc[2]] <- 1
            jj <- c(jj, loc[2])
        } else {
            print(paste("Edge from", loc[1], "to", loc[2], "is degenerate"))
        }
    }
    jj <- unique(jj)

    likelihood_ratios <- NULL
    df <- NULL
    p <- ncol(PiHat)

    if (length(jj) > 0) {
        for (jjj in 1:length(jj)) {
            s <- jj[jjj]
            PhiHatj <- PhiHat[, s]
            testj <- d_mat[, s]
            testaltj <- alt_mat[, s]
            ancestor <- setdiff(which(testaltj != 0), which(testj != 0))
            dfj <- 0

            if (median(N[, s]) > sum(testaltj != 0) + sum(PhiHatj != 0) + 1) {
                if (length(ancestor) > 0) {
                    an_mat <- testaltj
                    iv_mat <- PhiHatj
                    sigma_h1 <- estimate_Sigma(ss_back, iv_mat, an_mat, cormat, N, j)

                    an_mat <- testj
                    sigma_h0 <- estimate_Sigma(ss_back, iv_mat, an_mat, cormat, N, j)

                    likelihood_ratio <- (sigma_h0$RSS - sigma_h1$RSS) / sigma_h1$sigma2hat
                    likelihood_ratios <- c(likelihood_ratios, likelihood_ratio)
                    dfj <- dfj + 1
                    df <- c(df, dfj)
                    # print(likelihood_ratio)
                    # print(j)
                } else {
                    likelihood_ratios <- c(likelihood_ratios, 0)
                }
            }
        }
    }

    if (is.null(likelihood_ratios)) {
        likelihood_ratios <- 0
        df <- 0
    }
    return(list(likelihood_ratios = likelihood_ratios, df = df))
}
