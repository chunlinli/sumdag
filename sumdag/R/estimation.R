
#' Title Estimation of key quantities from GWAS summary statistics
#'
#' @param fullbetas A qxp matrix of estimated regression coefficient (beta) of outcome on each genotype, where q is the number of SNPs and p is the number of phenotype.
#' @param fullses A qxp matrix of the standard error of matrix 'fullbetas'
#' @param N A qxp matrix of corresponding sample size for each SNP-phenotype pair. A scalar is also acceptable when all SNP-phentype pairs have the same sample size
#' @param xtxdiag A constant number used to regularize the genetic correlation matrix; by default, no penalization is used.
#' @param ref ref = T if users provide a reference panel; ref = F if all SNPs are pruned and standardized and no reference panel is needed
#' @param geno.ref.txt The location of the user-provided reference panel with the genotype of the q SNPs; only txt files are acceptable; The order SNPs needs to be matched with fullbetas and fullses.
#'
#' @import utils stats
#' @return A list of quantities: qxq matrix xtx, p-vector ytyVec, qxp matrix xtyMat, q vector sds
#' @export
#'
#' @examples
#' \dontrun{
#'
#' }
summary_stats_preprocess <- function(fullbetas, fullses, N, xtxdiag = 0, ref = FALSE, geno.ref.txt = NULL) {
    if (length(N) == 1) paste("All SNP-phenotype pairs are assumed to have the same sample size", N)
    p <- ncol(fullbetas)
    q <- ncol(fullses)
    if (length(N) == 1) N <- matrix(rep(N, p * q), nrow = q, ncol = p)
    if (ref == T & is.null(geno.ref.txt)) stop("Please provide a reference txt file.")
    if (ref == F) stop("Genotypes are assumed to be independent and standardized.")
    if (nrow(fullbetas) != nrow(fullses)) stop("Dimension of beta and se(beta) do not match.")
    if (ncol(fullbetas) != ncol(fullses)) stop("Dimension of beta and se(beta) do not match.")

    # read in a txt file and calculate its LD matrix
    if (ref) {
        geno <- read.table(file = geno.ref.txt)
        n.ref <- nrow(geno)

        # calculate LD
        genoMat.std <- apply(geno, 2, function(x) scale(x, scale = T, center = T))
        tempMat <- t(genoMat.std) %*% (genoMat.std) / n.ref
        # calculate var of SNPs
        sds <- diag(tempMat)
        rm(geno)
        # if no genotype is provided, then we assume all SNPs are pruned to be independent and used a diagonal matrix for LD
    } else {
        # for standardized genotype
        sds <- rep(1, nrow(fullbetas))
        tempMat <- diag(sds)
    }

    # yty/N
    ytyVec <- NULL
    for (j in 1:ncol(fullbetas)) {
        tempEst <- vector()
        for (i in 1:nrow(fullbetas)) {
            ses <- fullses[, j]
            betas <- fullbetas[, j]
            tempEst[i] <- N[i, j] * sds[i] * ses[i]^2 + sds[i] * betas[i]^2
        }
        ytyEst <- median(tempEst)
        ytyVec <- c(ytyVec, ytyEst)
    }

    # calculate xty/N
    xtyMat <- NULL
    for (j in 1:ncol(fullbetas)) {
        betas <- fullbetas[, j]
        xtyEst <- sds[j] * betas
        xtyMat <- cbind(xtyMat, xtyEst)
    }

    # xtx/n
    xtx <- tempMat + diag(rep(xtxdiag, nrow(tempMat)))

    return(list(xtx = xtx, ytyVec = ytyVec, xtyMat = xtyMat, sds = sds))
}


#' Internal function: estimating V matrix based on summary statistic
#'
#' @param fullbetas A qxp matrix of estimated regression coefficient (beta) of outcome on each genotype, where q is the number of SNPs and p is the number of phenotype.
#' @param fullBetaSD A qxp matrix of the standard error of matrix 'fullbetas'
#' @param N A qxp matrix of corresponding sample size for each SNP-phenotype pair. A scalar is also acceptable when all SNP-phentype pairs have the same sample size
#' @param ss_back ss_back
#'
#' @import glmtlp stats
#' @examples
#' \dontrun{
#'
#' }
estimateV_sum <- function(fullbetas, fullBetaSD, N, ss_back) {
    p <- ncol(fullbetas)
    q <- nrow(fullbetas)
    if (length(N) == 1) N <- matrix(rep(N, p * q), nrow = q, ncol = p)

    V <- matrix(0, nrow = q, ncol = p)

    for (i in 1:p) {
        # standardize y
        xty <- ss_back$xtyMat[, i] / sqrt(ss_back$ytyVec[i])

        x1 <- glmtlp::sumtlp(XX = ss_back$xtx, Xy = xty, method = "tlp-c", kappa = 1:ncol(ss_back$xtx), nobs = median(N[, i]))

        # Extracting result for pseudo AIC/BIC function
        penBetas <- x1$beta
        uniBetas <- fullbetas[, i]
        uniSD <- fullBetaSD[, i]

        # Obtaining psuedo AIC/BIC values
        x2 <- pseudoBIC.internal(penalizedBetas = penBetas, betas = uniBetas, ses = uniSD, Ni = N[, i], sds = ss_back$sds, xtx = ss_back$xtx, ytyEst = 1, xtyEst = xty)

        # Saving V from minimum BIC
        V[, i] <- penBetas[, which.min(x2$bic)] # * sqrt(ss_back$ytyVec[i])
    }
    return(Vest = V)
}


#' Internal function: pseudoBIC
#'
#' This function computes pseudoBIC.
#'
#' @param penalizedBetas penalizedBetas
#' @param betas betas
#' @param ses ses
#' @param Ni Ni
#' @param sds sds
#' @param xtx xtx
#' @param ytyEst ytyEst
#' @param xtyEst xtyEst
#'
#' @import stats
#' @examples
#' \dontrun{
#'
#' }
pseudoBIC.internal <- function(penalizedBetas, betas, ses, Ni, sds, xtx, ytyEst, xtyEst) {
    SSEvec <- NULL
    qVec <- NULL
    aicVec <- NULL
    bicVec <- NULL
    bxxb <- NULL
    bxy <- NULL
    SSEvec <- NULL

    # because it is low dimension, so all the SNPs are used for estimating the sigma^2
    xtxinv <- solve(xtx)

    # obtain a consistent estimator for sigma^2 tilde
    qTemp <- nrow(xtx)
    sigSqTilde <- (ytyEst - t(xtyEst) %*% xtxinv %*% xtyEst) * median(Ni) / (median(Ni) - qTemp)

    # calculate SSE
    for (k in 1:ncol(penalizedBetas)) {
        penalizedBetasTemp <- penalizedBetas[, k]

        # q here is the number of parameters not equal to 0
        qVec <- sum(penalizedBetasTemp != 0)

        # bxxb for SSE (correct a mistake here on Oct 23)
        bxxb <- t(penalizedBetasTemp) %*% xtx %*% penalizedBetasTemp
        # bxy for SSE
        bxyTemp <- t(penalizedBetasTemp) %*% xtyEst

        # sse
        SSEest <- (ytyEst - (2 * bxyTemp) + bxxb) * median(Ni)
        SSEvec <- c(SSEvec, SSEest[1, 1])
        logLik <- -SSEest / (2 * sigSqTilde)

        aicTemp <- (2 * qVec) - (2 * logLik)
        bicTemp <- (log(median(Ni)) * qVec) - (2 * logLik)
        aicVec <- c(aicVec, aicTemp[1, 1])
        bicVec <- c(bicVec, bicTemp[1, 1])
    }

    toReturn <- structure(list(
        aic = aicVec, bic = bicVec, SSE = SSEvec,
        q = qVec, logLik = logLik, bxxb = bxxb, sigSqTilde = sigSqTilde
    ))
    return(toReturn)
}

#' Internal function: peeling
#'
#' This function implements peeling algorithm.
#'
#' @param V V matrix produced by peeling algorithm
#' @param thresh threshold
#'
#' @examples
#' \dontrun{
#'
#' }
peeling <- function(V, thresh) {
    peel.res <- topological_order(v = V, thresh = thresh)
    PiHat <- as.matrix(peel.res$an_mat) # ancestry
    colnames(PiHat) <- colnames(V)
    rownames(PiHat) <- colnames(V)
    PhiHat <- as.matrix(peel.res$iv_mat) # IV
    colnames(PhiHat) <- colnames(V)
    rownames(PhiHat) <- rownames(V)

    # calculate the effective number of tests
    ntest <- 0
    for (j in 1:nrow(PiHat)) {
        ancestor <- which(PiHat[, j] != 0)
        intervention <- which(PhiHat[, j] != 0)
        ntest <- ntest + sum(length(ancestor) != 0)
    }
    p <- ncol(PiHat)
    # effective number of test
    effect.test <- p^2 - p - ntest

    return(list(PiHat = PiHat, PhiHat = PhiHat, effect.test = effect.test))
}


#' Estimation of superDAG from GWAS summary statistics
#'
#' @param fullBeta A qxp matrix of estimated regression coefficient (beta) of outcome on each genotype, where q is the number of SNPs and p is the number of phenotype.
#' @param fullBetaSD A qxp matrix of the standard error of matrix 'fullbetas'
#' @param N A qxp matrix of corresponding sample size for each SNP-phenotype pair. A scalar is also acceptable when all SNP-phentype pairs have the same sample size
#' @param xtxdiag A constant number used to regularize the genetic correlation matrix; by default, no penalization is used.
#' @param ref ref = T if users provide a reference panel; ref = F if all SNPs are pruned and standardized and no reference panel is needed
#' @param geno.ref.text The location of the user-provided reference panel with the genotype of the q SNPs; only txt files are acceptable; The order SNPs needs to be matched with fullbetas and fullses.
#' @param thresh The threshold of V used in the peeling algorithm.
#' @param Xdata Xdata
#'
#' @return A list of three elements: PiHat, PhiHat, and effect.test
#' @export
#'
#' @examples
#' \dontrun{
#'
#' }
superDAG <- function(fullBeta, fullBetaSD, N, xtxdiag = 0, ref = F, geno.ref.text = NULL, thresh = 0.05, Xdata) {
    ss_back <- summary_stats_preprocess(fullBeta, fullBetaSD, N, xtxdiag, ref, geno.ref.text)
    Vest <- estimateV_sum(fullBeta, fullBetaSD, N, ss_back)
    rownames(Vest) <- Xdata$snp_names
    colnames(Vest) <- Xdata$prots
    peel_result <- peeling(V = Vest, thresh = thresh)
    return(list(PiHat = peel_result$PiHat, PhiHat = peel_result$PhiHat, effect.test = peel_result$effect.test))
}


#' TLP after peeling
#'
#' @param PiHat A pxp superset of the directed relations among phenotypes; output from superDAG function
#' @param PhiHat A qxp superset of the interventional relation; output from superDAG function
#' @param ss_back output from the summary_stats_preprocess function
#' @param cormat pxp correlation matrix of phenotypes; this matrix can be estimated by empricial correlation of null SNP.
#' @param N A qxp matrix of corresponding sample size for each SNP-phenotype pair. A scalar is also acceptable when all SNP-phentype pairs have the same sample size
#' @param penalizeXX penalizeXX = T if users prefer penalzing the IV matrix after peeling; penalizeXX = F if users do not.
#'
#' @import glmtlp stats
#' @return A list of two elements: W.new, U.new
#' @export
#'
#' @examples
#' \dontrun{
#'
#' }
TLP_afterpeel_U <- function(PiHat, PhiHat, ss_back, cormat, N, penalizeXX = T) {
    PiHat.new <- matrix(0, nrow = nrow(PiHat), ncol = ncol(PiHat)) # 15x15
    PhiHat.new <- matrix(0, nrow = nrow(PhiHat), ncol = ncol(PhiHat))

    p <- ncol(PiHat)
    q <- nrow(PhiHat)
    if (nrow(cormat) != nrow(PiHat)) stop("Number of phenotypes do not match.")
    if (length(N) == 1) N <- matrix(rep(N, p * q), nrow = q, ncol = p)

    for (j in 1:ncol(PiHat)) {
        PiHatj <- PiHat[, j]
        PhiHatj <- PhiHat[, j]
        Jp <- which(PhiHatj != 0) # IV set
        Ap <- which(PiHatj != 0) # ancestry set
        # standardize y
        xtyMat <- t(t(ss_back$xtyMat) / sqrt(ss_back$ytyVec))
        if (length(Jp) > 0) {
            XjpXjp <- ss_back$xtx[Jp, Jp, drop = F]
            XjpYp <- xtyMat[Jp, j]
            Nvec <- N[Jp, j]
            # when there is ancestor

            if (length(Ap) > 0) {
                XjpYap <- xtyMat[Jp, Ap, drop = F]
                YapYap <- cormat[Ap, Ap, drop = F]
                YapYp <- cormat[Ap, j]

                XTX <- rbind(cbind(XjpXjp, XjpYap), cbind(t(XjpYap), YapYap))
                XTY <- c(XjpYp, YapYp)

                # adding an option not to penalize XTX
                if (penalizeXX == T) {
                    x1 <- glmtlp::sumtlp(XX = XTX, Xy = XTY, method = "tlp-c", kappa = 1:ncol(XTX), nobs = median(Nvec))
                } else {
                    x1 <- glmtlp::sumtlp(XX = XTX, Xy = XTY, method = "tlp-c", kappa = 1:ncol(XTX), nobs = median(Nvec), penalty_factor = c(rep(0, length(Jp)), rep(1, length(Ap))))
                }

                penBetas <- x1$beta
                # Obtaining psuedo AIC/BIC values
                x2 <- pseudoBIC.internal(penalizedBetas = penBetas, Ni = median(Nvec), xtx = XTX, ytyEst = 1, xtyEst = XTY)

                PhiHat.new[Jp, j] <- penBetas[1:nrow(XjpXjp), which.min(x2$bic)]
                PiHat.new[Ap, j] <- penBetas[(nrow(XjpXjp) + 1):ncol(XTX), which.min(x2$bic)]
            }
        } else if (length(Jp) == 0 & length(Ap) > 0) {
            YapYap <- cormat[Ap, Ap, drop = F]
            YapYp <- cormat[Ap, j]
            BETAS <- solve(YapYap) %*% YapYp
            # not pe
            x1 <- glmtlp::sumtlp(XX = YapYap, Xy = YapYp, method = "tlp-c", kappa = 1:ncol(YapYap), nobs = median(Nvec))
            penBetas <- x1$beta
            # Obtaining psuedo AIC/BIC values
            x2 <- pseudoBIC.internal(penalizedBetas = penBetas, Ni = median(Nvec), xtx = YapYap, ytyEst = 1, xtyEst = XTY)

            PhiHat.new[Jp, j] <- penBetas[, which.min(x2$bic)]
            PiHat.new[Ap, j] <- 0
        }
    }
    rownames(PhiHat.new) <- rownames(PhiHat)
    colnames(PhiHat.new) <- colnames(PhiHat)
    rownames(PiHat.new) <- rownames(PiHat)
    colnames(PiHat.new) <- colnames(PiHat)
    return(list(W.est = PhiHat.new, U.est = PiHat.new))
}
