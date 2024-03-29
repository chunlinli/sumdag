% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/testing.R
\name{estimate_sigma}
\alias{estimate_sigma}
\title{Estimating variance of errors.}
\usage{
estimate_sigma(stats_list, an_mat, iv_mat, corr_Y, n, j)
}
\arguments{
\item{stats_list}{A list of XX, yy_vec, Xy_mat, maf_vec.
This is the output of preprocess() function.}

\item{an_mat}{A p by p ancestral relation matrix.}

\item{iv_mat}{A p by p interventional relation matrix.}

\item{corr_Y}{A p by p correlation matrix of Y.}

\item{n}{A p by p sample size matrix.}

\item{j}{The index of the target node.}
}
\value{
A list of sigma2 and RSS.
}
\description{
Estimating variance of errors.
}
\examples{
\dontrun{
beta_mat <- matrix(rnorm(33*23), 33, 23)
se_mat <- matrix(rnorm(33*23)^2, 33, 23)
n <- 3393
Xdata <- matrix(rnorm(500*33), 500, 33)
res <- estimate_superDAG(beta_mat, se_mat, n, Xdata)
an_mat <- res$an_mat
iv_mat <- res$iv_mat
y <- matrix(rnorm(100,23),100,23)
corr_Y <- cor(y)
stats_list <- preprocess(beta_mat, se_mat, n)
estimate_sigma(stats_list, an_mat, iv_mat, corr_Y, n)
}
}
