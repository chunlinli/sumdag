#######################################################################
# Step 0: load the package
#######################################################################
library(sumdag)

#######################################################################
# Step 1: read in the data
# SNP info including MAF is in "IV_info.RData", but not needed here
#######################################################################
beta_mat <- read.table("./example/fullBeta_revised_ukbb.txt")
se_mat <- read.table("./example/fullSE_revised_ukbb.txt")
pval_mat <- read.table("./example/fullP_revised_ukbb.txt")
beta_mat <- as.matrix(beta_mat)
se_mat <- as.matrix(se_mat)
pval_mat <- as.matrix(pval_mat)

# sample size
n <- matrix(3393, nrow = nrow(beta_mat), ncol = ncol(beta_mat))

dim(pval_mat)
# [1] 33 23
count <- apply(pval_mat, 2, function(x) sum(x < 5e-8))
sum(count > 0)
# [1] 23

snp <- rownames(beta_mat)
protein <- colnames(beta_mat)
p <- length(protein)

# read in reference panel
load("./example/ukbb_xtx.RData")
XX <- xtx

# read in the phenotype correlation matrix
cor_mat <- read.table("./example/23protein_phe_corr_allchr.txt",
    check.names = FALSE
)
cor_mat <- cor_mat[protein, protein]

#######################################################################
# Step 2: estimate the quantities estimated from the summary stats
#######################################################################
stats_list <- preprocess(
    beta_mat = beta_mat,
    se_mat = se_mat,
    n = n,
    s_ = 0,
    maf_vec = NULL,
    geno_ref = XX
)

#######################################################################
# Step 3: estimate the V matrix
#######################################################################
V <- estimate_V(beta_mat, se_mat, n, stats_list)
rownames(V) <- snp
colnames(V) <- protein

#######################################################################
# Step 4: peeling
#######################################################################
peel_result <- peeling(V = V, thresh = 0)
an_mat <- as.matrix(peel_result$an_mat)
iv_mat <- as.matrix(peel_result$iv_mat)

#######################################################################
# Step 5: TLP estimation after peeling
#######################################################################
result <- TLP_U(
    an_mat = an_mat,
    iv_mat = iv_mat,
    stats_list = stats_list,
    as.matrix(cor_mat), n
)
# U matrix
U <- result$U

# W matrix
W <- result$W

summary(abs(U[U != 0]))
summary(abs(W[W != 0]))

#######################################################################
# Step 6: test each edge
#######################################################################
cov_mat <- diag(sqrt(stats_list$yy_vec)) %*%
    as.matrix(cor_mat) %*% diag(sqrt(stats_list$yy_vec))
p_value <- matrix(1, nrow = nrow(U), ncol = nrow(U))
pair <- list()
for (i in 1:nrow(U)) {
    for (j in 1:ncol(U)) {
        if (an_mat[i, j] != 0) {
            pair[[1]] <- c(i, j)
            likelihood_ratios <- asymptotic_inference_internal(
                stats_list = stats_list,
                an_mat = an_mat,
                iv_mat = iv_mat,
                pairs = pair[1],
                corr_Y = cov_mat,
                n = n
            )
            likelihood_ratio <- sum(likelihood_ratios$likelihood_ratios)
            df <- sum(likelihood_ratios$df)
            p_value[i, j] <- pchisq(
                likelihood_ratio,
                df = df,
                lower.tail = FALSE
            )
        }
    }
}

#######################################################################
# Step 7: plot the graph
#######################################################################
U <- an_mat
colnames(U) <- gsub("[.]", "-", colnames(U))

p_value[p_value >= 0.05 / (23 * 22 - sum(an_mat != 0))] <- 0
pval <- p_value
res <- t(pval)[t(U) == 1]

res <- ifelse(res == 0, 0, 1)

source("./example/draw_dag_code.R")
plot_graph(U, res,
    highlightnode = c(
        "ADM", "PGF", "TNFRSF11B", "TEK",
        "MMP3", "CHI3L1", "IL1RL1", "HAVCR1",
        "MMP10", "CXCL6", "CXCL16", "LGALS3",
        "IL6R", "CTSD"
    ),
    underlightnode = "RETN",
    graph_name = "./example/protein23_ukbb.pdf"
)
