
#' Estimating key quantities from GWAS summary data
#'
#' @param beta_mat A q by p matrix.
#' Regression coefficient of q phenptypes on p genotypes.
#' @param se_mat A q by p matrix.
#' The standard errors of beta_mat.
#' @param n A q by p matrix or a scalar.
#' Sample sizes for corresponding genotype-phenotype pairs.
#' A scalar indicates all genotype-phentype pairs have the same sample size.
#' @param s_ A nonnegative scalar.
#' Regularize the genetic correlation matrix. The default is 0.
#' @param maf_vec A q vector of MAF for the genotypes.
#' Can be skipped if a genotype covariance matrix is given.
#' @param geno_ref A q by q genotype covariance matrix.
#'
#' @import utils stats
#' @return A list of quantities:
#' q by q matrix XX, p vector yy_vec, q by p matrix Xy_mat, q vector maf_vec
#' @export
#'
#' @examples
#' beta_mat <- matrix(rnorm(33*23), 33, 23)
#' se_mat <- matrix(rnorm(33*23)^2, 33, 23)
#' n <- 3393
#' preprocess(beta_mat, se_mat, n)
#'
preprocess <- function(beta_mat, se_mat, n,
                       s_ = 0, maf_vec = NULL, geno_ref = NULL) {
    if (all.equal(dim(beta_mat), dim(se_mat)) == FALSE) {
        stop("Dimensions of beta_mat and se_mat do not match.")
    }
    p <- ncol(beta_mat)
    q <- nrow(beta_mat)

    if (length(n) == 1) {
        message(paste("All pairs have the same sample size ", n))
        n <- matrix(n, nrow = q, ncol = p)
    }

    if (is.null(geno_ref) && is.null(maf_vec)) {
        stop("Need to specify geno_ref or maf_vec.")
    }
    if (is.null(geno_ref) && !is.null(maf_vec)) {
        message("Genotypes are assumed to be independent.")
    }
    if (!is.null(geno_ref) && !is.null(maf_vec)) {
        warning("maf_vec are ignored.")
    }
    if (!is.null(geno_ref)) {
        XX <- geno_ref + diag(s_, nrow(geno_ref))
        maf_vec <- diag(XX)
    } else {
        XX <- diag(maf_vec) + diag(s_, nrow(maf_vec))
    }

    yy_vec <- rep(NA, p)
    for (j in seq_len(p)) {
        yy_vec[j] <- median(
            n[, j] * maf_vec * se_mat[, j]^2 + maf_vec * beta_mat[, j]^2
        )
    }

    Xy_mat <- NULL
    for (j in 1:ncol(beta_mat)) {
        betas <- beta_mat[, j]
        xy <- maf_vec[j] * betas
        Xy_mat <- cbind(Xy_mat, xy)
    }

    return(list(XX = XX, yy_vec = yy_vec, Xy_mat = Xy_mat, maf_vec = maf_vec))
}

#' (Internal function) Estimating V matrix based on summary data
#'
#' @param beta_mat A q by p matrix.
#' Regression coefficient of q phenptypes on p genotypes.
#' @param se_mat A q by p matrix.
#' The standard errors of beta_mat.
#' @param n A q by p matrix or a scalar.
#' Sample sizes for corresponding genotype-phenotype pairs.
#' A scalar indicates all genotype-phentype pairs have the same sample size.
#' @param stats_list A list of XX, yy_vec, Xy_mat, maf_vec.
#' This is the output of preprocess() function.
#'
#' @import glmtlp stats
#' @return A q by p matrix V.
#' @examples
#' \dontrun{
#' beta_mat <- matrix(rnorm(33*23), 33, 23)
#' se_mat <- matrix(rnorm(33*23)^2, 33, 23)
#' n <- 3393
#' estimate_V(beta_mat, se_mat, preprocess(beta_mat, se_mat, n))
#' }
estimate_V <- function(beta_mat, se_mat, n, stats_list) {
    p <- ncol(beta_mat)
    q <- nrow(beta_mat)
    if (length(n) == 1) n <- matrix(n, nrow = q, ncol = p)

    V <- matrix(0, nrow = q, ncol = p)

    for (j in seq_len(p)) {
        Xy <- stats_list$Xy_mat[, j] / sqrt(stats_list$yy_vec[j])

        m1 <- glmtlp::sumtlp(
            XX = stats_list$XX, Xy = Xy, method = "tlp-c",
            kappa = seq_len(ncol(stats_list$XX)), nobs = median(n[, j])
        )

        m2 <- pseudoBIC(
            beta_vec = m1$beta,
            n = n[, j],
            XX = stats_list$XX,
            yy = 1,
            Xy = Xy
        )

        V[, j] <- m1$beta[, which.min(m2$BIC)]
    }

    return(V = V)
}

#' (Internal function) Computing pseudoBIC
#'
#' @param beta_vec A q vector.
#' Coefficients of the joint regression of a phenotype over genotypes.
#' @param n A q vector of sample sizes.
#' @param XX A q by q genotype covariance matrix.
#' @param yy A scalar.
#' @param Xy A q vector.
#'
#' @import stats
#' @examples
#' \dontrun{
#' X <- matrix(rnorm(100*10), 100, 10)
#' y <- rnorm(100)
#' XX <- t(X) %*% X
#' Xy <- t(X) %*% y
#' yy <- t(y) %*% y
#' n <- 100
#' beta_vec <- matrix(rnorm(10*5), 10, 5)
#' pseudoBIC(beta_vec, n, XX, yy, Xy)
#' }
pseudoBIC <- function(beta_vec, n, XX, yy, Xy) {
    SSE_vec <- rep(NA, ncol(beta_vec))
    AIC_vec <- rep(NA, ncol(beta_vec))
    BIC_vec <- rep(NA, ncol(beta_vec))

    err_var <- (yy - t(Xy) %*% solve(XX) %*% Xy) *
        median(n) / (median(n) - nrow(XX))

    for (k in seq_len(ncol(beta_vec))) {
        b <- beta_vec[, k]
        df <- sum(b != 0) # number of nonzero coefficients, degree of freedom
        bXXb <- t(b) %*% XX %*% b
        bXy <- t(b) %*% Xy

        SSE_vec[k] <- (yy - 2 * bXy + bXXb) * median(n)
        loglik <- -SSE_vec[k] / (2 * err_var)
        AIC_vec[k] <- 2 * df - 2 * loglik
        BIC_vec[k] <- log(median(n)) * df - 2 * loglik
    }

    return(list(
        AIC = AIC_vec,
        BIC = BIC_vec,
        SSE = SSE_vec,
        err_var = err_var
    ))
}

#' (Internal function) Peeling procedure implementation.
#'
#' @param V A q by p matrix. V matrix.
#' @param thresh A positive scalar. Threshold.
#'
#' @examples
#' \dontrun{
#' v <- matrix(rnorm(15*10), 15, 10)
#' thresh <- 0.8
#' topological_order(v, thresh)
#' }
peeling <- function(V, thresh) {
    peel_result <- topological_order(v = V, thresh = thresh)

    an_mat <- peel_result$an_mat
    colnames(an_mat) <- colnames(V)
    rownames(an_mat) <- colnames(V)

    iv_mat <- peel_result$iv_mat
    colnames(iv_mat) <- colnames(V)
    rownames(iv_mat) <- rownames(V)

    # calculate the effective number of tests
    n_test <- 0
    for (j in seq_len(nrow(an_mat))) {
        ancestor <- which(an_mat[, j] != 0)
        n_test <- n_test + sum(length(ancestor) != 0)
    }
    n_test <- ncol(V)^2 - ncol(V) - n_test

    return(list(
        an_mat = an_mat,
        iv_mat = iv_mat,
        n_test = n_test
    ))
}


#' Estimating super-DAG from GWAS summary data
#'
#' @param beta_mat A q by p matrix.
#' Regression coefficient of q phenptypes on p genotypes.
#' @param se_mat A q by p matrix.
#' The standard errors of beta_mat.
#' @param n A q by p matrix or a scalar.
#' Sample sizes for corresponding genotype-phenotype pairs.
#' A scalar indicates all genotype-phentype pairs have the same sample size.
#' @param s_ A nonnegative scalar.
#' Regularize the genetic correlation matrix. The default is 0.
#' @param maf_vec A q vector of MAF for the genotypes.
#' Can be skipped if a genotype covariance matrix is given.
#' @param geno_ref A q by q genotype covariance matrix.
#' @param thresh A scalar. The threshold used in peeling procedure.
#' @param Xdata #TODO Genotype matrix in reference panel.
#'
#' @return A list of three quantities:
#' A p by p matrix an_mat, a q by p matrix iv_mat, and a scalar n_test
#' @export
#'
#' @examples
#' beta_mat <- matrix(rnorm(33*23), 33, 23)
#' se_mat <- matrix(rnorm(33*23)^2, 33, 23)
#' n <- 3393
#' Xdata <- matrix(rnorm(500*33), 500, 33)
#' estimate_superDAG(beta_mat, se_mat, n, Xdata)
estimate_superDAG <- function(beta_mat, se_mat, n,
                              s_ = 0, maf_vec = NULL, geno_ref = NULL,
                              thresh = 0.05, Xdata) {
    stats_list <- preprocess(
        beta_mat = beta_mat,
        se_mat = se_mat,
        n = n,
        s_ = s_,
        maf_vec = maf_vec,
        geno_ref = geno_ref
    )
    V <- estimate_V(
        beta_mat = beta_mat,
        se_mat = se_mat,
        n = n,
        stats_list = stats_list
    )

    rownames(V) <- Xdata$snp_names
    colnames(V) <- Xdata$prots

    peel_result <- peeling(V = V, thresh = thresh)

    return(list(
        an_mat = peel_result$an_mat,
        iv_mat = peel_result$iv_mat,
        n_test = peel_result$n_test
    ))
}


#' Truncated Lasso penalized estimation after peeling
#'
#' @param an_mat A p by p ancestral relation matrix
#' @param iv_mat A q by p interventional relation matrix
#' @param stats_list A list output by preprocess() function
#' @param corr_Y A p by p correlation matrix of phenotypes;
#' Can be estimated by empricial correlation of null SNPs.
#' @param n A q by p matrix or a scalar.
#' Sample sizes for corresponding genotype-phenotype pairs.
#' A scalar indicates all genotype-phentype pairs have the same sample size.
#' @param if_penalize_iv Binary.
#' TRUE if users prefer penalzing the IVs; FALSE otherwise.
#'
#' @import glmtlp stats
#' @return A list of two elements:
#' A p by p adjacency matrix U, a q by p interventional matrix W
#' @export
#'
#' @examples
#' \dontrun{
#' beta_mat <- matrix(rnorm(33*23), 33, 23)
#' se_mat <- matrix(rnorm(33*23)^2, 33, 23)
#' n <- 3393
#' Xdata <- matrix(rnorm(500*33), 500, 33)
#' res <- estimate_superDAG(beta_mat, se_mat, n, Xdata)
#' an_mat <- res$an_mat
#' iv_mat <- res$iv_mat
#' y <- matrix(rnorm(100,23),100,23)
#' corr_Y <- cor(y)
#' stats_list <- preprocess(beta_mat, se_mat, n)
#' TLP_U(an_mat, iv_mat, stats_list, corr_Y, n)
#' }
TLP_U <- function(an_mat, iv_mat, stats_list,
                  corr_Y, n, if_penalize_iv = TRUE) {
    p <- ncol(an_mat)
    q <- nrow(iv_mat)
    U <- matrix(0, nrow = p, ncol = p)
    W <- matrix(0, nrow = q, ncol = p)

    if (all.equal(dim(corr_Y), c(p, p)) == FALSE) {
        stop("Dimensions of corr_Y do not match.")
    }
    if (length(n) == 1) n <- matrix(n, nrow = q, ncol = p)

    for (j in seq_len(p)) {
        AN <- which(an_mat[, j] != 0) # ancestor set
        IV <- which(iv_mat[, j] != 0) # IV set

        Xy_mat <- sweep(stats_list$Xy_mat,
            MARGIN = 2,
            sqrt(stats_list$yy_vec), FUN = "/"
        )

        if (length(IV) > 0) {
            XX_ <- stats_list$XX[IV, IV, drop = FALSE]
            Xy_ <- Xy_mat[IV, j]
            n_ <- n[IV, j]

            if (length(AN) > 0) {
                XY_ <- Xy_mat[IV, AN, drop = FALSE]
                YY_ <- corr_Y[AN, AN, drop = FALSE]
                Yy_ <- corr_Y[AN, j]

                ZZ <- rbind(cbind(XX_, XY_), cbind(t(XY_), YY_))
                Zy <- c(Xy_, Yy_)

                if (if_penalize_iv) {
                    m1 <- glmtlp::sumtlp(
                        XX = ZZ, Xy = Zy, method = "tlp-c",
                        kappa = 1:ncol(ZZ), nobs = median(n_)
                    )
                } else {
                    m1 <- glmtlp::sumtlp(
                        XX = ZZ, Xy = Zy, method = "tlp-c",
                        kappa = 1:ncol(ZZ), nobs = median(n_),
                        penalty_factor = c(
                            rep(0, length(IV)),
                            rep(1, length(AN))
                        )
                    )
                }

                b <- m1$beta
                m2 <- pseudoBIC(
                    beta_vec = b, n = median(n_),
                    XX = ZZ,
                    yy = 1,
                    Xy = Zy
                )

                W[IV, j] <- b[1:length(IV), which.min(m2$BIC)]
                U[AN, j] <- b[(length(IV) + 1):ncol(ZZ), which.min(m2$BIC)]
            } else {
                XX_ <- stats_list$XX[IV, IV, drop = FALSE]
                Xy_ <- Xy_mat[IV, j]
                n_ <- n[IV, j]

                if (if_penalize_iv) {
                    m1 <- glmtlp::sumtlp(
                        XX = XX_, Xy = Xy_, method = "tlp-c",
                        kappa = 1:ncol(XX_), nobs = median(n_)
                    )
                } else {
                    m1 <- glmtlp::sumtlp(
                        XX = XX_, Xy = Xy_, method = "tlp-c",
                        kappa = 1:ncol(XX_), nobs = median(n_),
                        penalty_factor = c(rep(0, length(IV)))
                    )
                }

                b <- m1$beta
                m2 <- pseudoBIC(
                    beta_vec = b, n = median(n_),
                    XX = XX_,
                    yy = 1,
                    Xy = Xy_
                )

                W[IV, j] <- b[, which.min(m2$BIC)]
            }
        } else if (length(IV) == 0 && length(AN) > 0) {
            YY_ <- corr_Y[AN, AN, drop = FALSE]
            Yy_ <- corr_Y[AN, j]

            m1 <- glmtlp::sumtlp(
                XX = YY_, Xy = Yy_, method = "tlp-c",
                kappa = 1:ncol(YY_), nobs = median(n[, j])
            )
            b <- m1$beta

            m2 <- pseudoBIC(
                beta_vec = b,
                n = median(n[, j]),
                XX = YY_,
                yy = 1,
                Xy = Yy_
            )
            U[AN, j] <- b[, which.min(m2$BIC)]
        }
    }

    rownames(U) <- colnames(an_mat)
    colnames(U) <- colnames(an_mat)
    rownames(W) <- rownames(iv_mat)
    colnames(W) <- colnames(an_mat)
    return(list(U = U, W = W))
}


#' (Internal function) Peeling procedure implementation.
#'
#' @param V A q by p matrix. V matrix.
#' @param thresh A positive scalar. Threshold.
#'
#' @examples
#' \dontrun{
#' v <- matrix(rnorm(15*10), 15, 10)
#' thresh <- 0.8
#' topological_order(v, thresh)
#' }
topological_order <- function(v, thresh = 0.15) {
    p <- ncol(v)
    q <- nrow(v)
    if (q < p) stop("No sufficient interventions: q < p.")

    v[abs(v) < thresh] <- 0

    an_mat <- matrix(0, p, p)
    in_mat <- matrix(0, q, p)
    iv_mat <- matrix(0, q, p)

    removed_x <- rep(FALSE, q)
    removed_y <- rep(FALSE, p)

    v_abs <- abs(v)
    v_nz <- (v_abs != 0)

    # check if there is a primary variable without intervention.
    while (!all(removed_y)) {

        # leaf-instrument pairs
        iv_targets <- rowSums(as.matrix(v_nz[, !removed_y]))
        one <- min(iv_targets[iv_targets > 0 & !removed_x])
        # Warning message: no non-missing arguments to min;
        # returning Inf (there is trait with no IV assigned;
        # happens when there is invalid IV or no IV assigned in the first place)

        leaf_iv <- which(!removed_x & iv_targets == one)
        # leaf is the starting point of the DAG
        if (length(leaf_iv) == 0) break
        leaf <- rep(NA, length(leaf_iv))
        leaf_iter <- 1
        for (l in leaf_iv) {
            j <- which(v_abs[l, ] == max(v_abs[l, !removed_y]) & !removed_y)[1]
            iv_mat[l, j] <- 1
            in_mat[l, j] <- 1
            leaf[leaf_iter] <- j
            leaf_iter <- leaf_iter + 1
        }
        leaf <- unique(leaf)

        # leaf-noninstrument pairs
        for (j in leaf) {
            leaf_interventions <- which(!removed_x & v_abs[, j] != 0)
            leaf_noniv <- setdiff(leaf_interventions, leaf_iv)
            in_mat[leaf_noniv, j] <- 1
        }

        # ancestral relations
        for (j in leaf) {
            j_instrument <- which(iv_mat[, j] != 0)
            j_descendant <- vector("list", length = length(j_instrument))
            for (l in seq_len(length(j_instrument))) {
                l2 <- j_instrument[l]
                j_descendant[[l]] <- which(removed_y & v_abs[l2, ] != 0)
            }

            j_descendant <- Reduce(union, j_descendant)
            an_mat[j, j_descendant] <- 1
        }

        # peeling-off
        removed_y[leaf] <- TRUE
        removed_x[leaf_iv] <- TRUE
        # one IV can only correspond to one trait (invalid IV if two traits)
    }

    # reconstruction of topological order
    an_mat <- (solve(diag(p) - an_mat) != 0) - diag(p)
    in_mat <- 1 * (in_mat %*% (diag(p) + an_mat) > 0)

    list(an_mat = an_mat, in_mat = in_mat, iv_mat = iv_mat)
}
