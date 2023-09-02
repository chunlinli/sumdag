# Estimating and Inferring a directed acyclic graph from summary statistics

This repository contains the implementation of the paper

- R Zilinskas, C Li, X Shen, W Pan, T Yang (2023). [Inferring a Directed Acyclic Graph of Phenotypes from GWAS Summary Statistics](https://www.biorxiv.org/content/10.1101/2023.02.10.528092v1.abstract). bioRxiv. Under review.


## Contents

The directory `./sumdag` contains the R package implementing the proposed approach.

The directory `./example` contains the R code for replicating the real data application.

The directory `./simulations` contains the R code for simulation studies.

The directory `./shiny/shinyapp` contains the R code for Shiny App. 

## Preliminaries

The R environment used for development is R version 4.2.1. 

The package `sumdag` requires installation of the following R packages:
```r
packages <- c("ncvreg","lassosum","MASS")
install.packages(packages)
devtools::install_github("chunlinli/glmtlp")
```
Then install `sumdag` from this GitHub repository:
```r
devtools::install_github("chunlinli/sumdag/sumdag")
```


## Shiny

To start the Shiny App, run the following in R
```R
library(shiny)
runApp("./shiny/shinyapp")
```
## Usage and example

To replicate the real data analysis, run 
```bash
Rscript ./example/demo.R
```
This will replicate the plot of the real data analysis. Here, we illustrate the usage of package by running `./example/demo` step by step.

Assume the working directory is the cloned repository [https://github.com/chunlinli/sumdag](https://github.com/chunlinli/sumdag). In R, load the package `sumdag`:
```r
library(sumdag)
```

Load the summary data. NOTE: SNP information, including MAF, is in `./example/IV_info.RData`, but is not needed here.

```R
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

# read in the UK biobank reference panel
load("./example/ukbb_xtx.RData")
XX <- xtx

# read in the phenotype correlation matrix
cor_mat <- read.table("./example/23protein_phe_corr_allchr.txt",
    check.names = FALSE
)
cor_mat <- cor_mat[protein, protein]
```

Obtain the quantities estimated from the summary statistics using `preprocess()` function.

```R
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
```

Estimate $V$ matrix.

```R
#######################################################################
# Step 3: estimate the V matrix
#######################################################################
V <- estimate_V(beta_mat, se_mat, n, stats_list)
rownames(V) <- snp
colnames(V) <- protein
```

Implement the peeling algorithm in Li et al., (2023): obtain the ancestral relation graph `an_mat` and interventional relation matrix `iv_mat`.

```R
#######################################################################
# Step 4: peeling
#######################################################################
peel_result <- peeling(V = V, thresh = 0)
an_mat <- as.matrix(peel_result$an_mat)
iv_mat <- as.matrix(peel_result$iv_mat)

```

Run truncated-Lasso estimation based on the peeling results.

```R
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
U

# W matrix
W <- result$W
W
```

Now, conduct likelihood ratio test for each edge.
```R
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
```

Finally, plot the graph of significant edges.

```R
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
```




## Citing information

If you find the code useful, please consider citing 

```tex
@article{zilinskas2023inferring,
    author = {Rachel Zilinskas, Chunlin Li, Xiaotong Shen, Wei Pan, Tianzhong Yang},
    title = {Inferring a directed acyclic graph of phenotypes from {GWAS} summary statistics},
    year = {2023},
    journal = {Under review}
}
```

The code is maintained on GitHub and the R package will be uploaded to CRAN soon. This project is in active development.

Implementing the structure learning algorithms is error-prone. If you spot any error, please file an issue [here](https://github.com/chunlinli/sumdag/issues) or contact me via email -- I will be grateful to be informed.