source("./R/testing_path2.R")

# test each edge
covMat <- diag(sqrt(ss_back$ytyVec)) %*% as.matrix(corMat) %*% diag(sqrt(ss_back$ytyVec))
p_value <- matrix(1, nrow = nrow(U.TLP), ncol = nrow(U.TLP))
pair <- list()
for (i in 1:nrow(U.TLP)) {
    for (j in 1:ncol(U.TLP)) {
        if (PiHat[i, j] != 0) {
            pair[[1]] <- c(i, j)
            likelihood_ratios <- asymptotic_inference_internal(ss_back, PiHat, PhiHat, pair[1], covMat, fullBetaSS)
            likelihood_ratio <- sum(likelihood_ratios$likelihood_ratios)
            df <- sum(likelihood_ratios$df)
            p_value[i, j] <- pchisq(likelihood_ratio, df = df, lower.tail = FALSE)
        }
    }
}
