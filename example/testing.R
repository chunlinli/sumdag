# testing


estimate_sigma <- function(stats_list, an_mat, iv_mat, corr_Y, n, j) {
    AN <- which(an_mat[, j] != 0)
    IV <- which(iv_mat[, j] != 0)

    if (length(IV) > 0 && length(AN) > 0) {
        XX_ <- stats_list$XX[IV, IV, drop = FALSE]
        Xy_ <- stats_list$Xy[IV, j]
        XY_ <- stats_list$Xy[IV, AN, drop = FALSE]
        YY_ <- corr_Y[AN, AN, drop = FALSE]
        Yy_ <- corr_Y[AN, j]

        b <- solve(rbind(cbind(XX_, XY_), cbind(t(XY_), YY_))) %*% c(Xy_, Yy_)
        n_ <- n[IV, j]
        RSS <- (stats_list$yy_vec[j] - t(c(Xy_, Yy_)) %*% b) * median(n_)
        sigma2 <- RSS / (median(n_) - length(AN) - length(IV) + 1)
    } else if (length(IV) == 0 && length(AN) > 0) {
        YY_ <- corr_Y[AN, AN, drop = FALSE]
        Yy_ <- corr_Y[AN, j]
        b <- solve(YY_) %*% Yy_
        n_ <- n[, j]
        RSS <- (stats_list$yy_vec[j] - t(Yy_) %*% b) * median(n_)
        sigma2 <- RSS / (median(n_) - length(AN) + 1)
    } else if (length(IV) > 0 && length(AN) == 0) {
        XX_ <- stats_list$XX[IV, IV, drop = FALSE]
        Xy_ <- stats_list$Xy[IV, j]
        b <- solve(XX_) %*% Xy_
        n_ <- n[IV, j]
        RSS <- (stats_list$yy_vec[j] - t(Xy_) %*% b) * median(n_)
        sigma2 <- RSS / (median(n_) - length(IV) + 1)
    } else if (length(IV) == 0 && length(AN) == 0) {
        n_ <- n[, j]
        RSS <- stats_list$yy_vec[j] * median(n_)
        sigma2 <- RSS / (median(n_) - 1)
    }

    return(list(sigma2 = sigma2, RSS = RSS))
}

causal_inference <- function(stats_list, an_mat, iv_mat,
                             pairs, corr_Y, n, test_type = c("edge", "path")) {
    if (test_type == "edge") {
        likelihood_ratios <- asymptotic_inference_internal(
            stats_list,
            an_mat,
            iv_mat,
            pairs, corr_Y, n
        )
        likelihood_ratio <- sum(likelihood_ratios$likelihood_ratios)
        df <- sum(likelihood_ratios$df)
        p_value <- pchisq(likelihood_ratio, df = df, lower.tail = FALSE)

        list(likelihood_ratio = likelihood_ratio, df = df, p_value = p_value)
    } else if (test_type == "path") {
        likelihood_ratios <- asymptotic_inference_internal(
            stats_list,
            an_mat,
            iv_mat,
            pairs, corr_Y, n
        )
        likelihood_ratio <- likelihood_ratios$likelihood_ratios
        p_value <- max(pchisq(likelihood_ratio, df = 1, lower.tail = FALSE))

        list(likelihood_ratio = likelihood_ratios, p_value = p_value)
    } else {
        stop("Invalid test type.")
    }
}

asymptotic_inference_internal <- function(stats_list, an_mat, iv_mat,
                                          pairs, corr_Y, n) {
    d_mat <- an_mat
    alt_mat <- an_mat

    jj <- NULL
    for (j in seq_len(length(pairs))) {
        loc <- pairs[[j]]
        if (an_mat[loc[2], loc[1]] == 0) {
            d_mat[loc[1], loc[2]] <- 0
            alt_mat[loc[1], loc[2]] <- 1
            jj <- c(jj, loc[2])
        } else {
            message(paste(
                "Edge from", loc[1], "to", loc[2],
                "is degenerate"
            ))
        }
    }
    jj <- unique(jj)

    likelihood_ratios <- NULL
    df <- NULL

    if (length(jj) > 0) {
        for (j in jj) {
            iv_mat_j <- iv_mat[, j]
            test_j <- d_mat[, j]
            test_alt_j <- alt_mat[, j]
            ancestor <- setdiff(
                which(test_alt_j != 0),
                which(test_j != 0)
            )

            df_j <- 0

            if (median(n[, j]) >
                sum(test_alt_j != 0) + sum(iv_mat_j != 0) + 1) {
                if (length(ancestor) > 0) {
                    sigma_h1 <- estimate_sigma(
                        stats_list,
                        iv_mat = iv_mat,
                        an_mat = alt_mat,
                        corr_Y = corr_Y,
                        n = n,
                        j = j
                    )

                    sigma_h0 <- estimate_sigma(
                        stats_list,
                        iv_mat = iv_mat,
                        an_mat = d_mat,
                        corr_Y = corr_Y,
                        n = n,
                        j = j
                    )

                    likelihood_ratio <- (sigma_h0$RSS - sigma_h1$RSS) /
                        sigma_h1$sigma2
                    likelihood_ratios <- c(
                        likelihood_ratios,
                        likelihood_ratio
                    )

                    df_j <- df_j + 1
                    df <- c(df, df_j)

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
