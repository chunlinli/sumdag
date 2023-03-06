# setting parameters for the network
# cd /home/yang3704/yang3704/network/

# setting 1: original structure (corrected on Nov 27th)
sim1 <- function(ES, Utxt, Wtxt) {
    U <- as.matrix(read.table(Utxt)) * ES # pxp
    W <- as.matrix(read.table(Wtxt)) * ES # qxp
    p <- nrow(U)
    I <- diag(p)
    S <- solve(I - U)
    Vsim <- W %*% S # qxp
    Vsim[which(abs(Vsim) < 5e-08)] <- 0
    return(list(U = U, W = W, Vsim = Vsim))
}


# setting 2: one IV per protein with max strength
sim2 <- function(ES, Utxt, Wtxt) {
    U <- as.matrix(read.table(Utxt)) * ES
    W <- as.matrix(read.table(Wtxt)) * ES
    gene <- which(apply(W, 2, function(x) sum(x != 0)) > 1)
    for (i in gene) {
        W[, i][abs(W[, i]) != max(abs(W[, i]))] <- 0
    }
    p <- nrow(U)
    I <- diag(p)
    S <- solve(I - U)
    Vsim <- W %*% S
    Vsim[which(abs(Vsim) < 5e-08)] <- 0
    return(list(U = U, W = W, Vsim = Vsim))
}

# input necessary data files
# K is number of replications
# sim is the simulation setting
# samplesize.sim is the sample size to generate summary stats
input_data <- function(samplesize.sim, sim, sigma2txt) {

    # read in genotype and standardized the genotype
    X <- read.table(paste0("../simulations/chr22.26SNPs.sum.", samplesize.sim, "_3.txt"))
    snp_names <- colnames(X)
    X <- t(X)
    X.std <- t(apply(X, 1, function(x) scale(x, scale = T, center = T)))

    # Reading in residual errors (stayed the same across simulation)
    sigma2 <- read.table(sigma2txt)
    sigma2 <- sigma2$V1

    # load ref and alternative allele
    # sim_snps_bim <- read.table("../simulations/chr22.26SNPs.ref.3000_3.bim")
    # names(sim_snps_bim) <- c("CHR", "ID", "POSG", "BP", "A1", "A2")

    # proteins name
    prots <- sort(c("FABP4", "HSP_27", "CXCL16", "MMP.10", "U.PAR", "KIM.1", "CHI3L1", "IL.1ra", "MMP.3", "FS", "NT.pro_BNP", "hK11", "MMP.12", "TM", "HGF"))

    return(list(X.std = X.std, snp_names = snp_names, sigma2 = sigma2, prots = prots))
}

# simulate individual level data
simulate_Y <- function(X.std, sigma2, sim, prots_name) {
    U <- sim$U
    Vsim <- sim$Vsim
    p <- nrow(U)
    I <- diag(p)
    # simulating Y
    eps <- t(MASS::mvrnorm(ncol(X.std), rep(0, length(sigma2)), diag(sigma2)))
    eps_v <- (solve(I - t(U))) %*% eps
    # use the standardized X here for varaible selection purpose
    Y <- (t(Vsim)) %*% X.std + eps_v
    rownames(Y) <- prots_name
    return(Y = Y)
}

# generate summary level data for the generated individual level data
# data is an output of simulate_Y
cal_sum_stats <- function(X.std, snp_names, sim_snps_bim, Y) {
    # system('mkdir ./sum_stats_temp/')

    # data= list(Y=Y, X.std=X.std, snp_names=snp_names, prots_name=prots)
    p <- nrow(Y)
    q <- nrow(X.std)
    sim_snps <- snp_names
    prots <- rownames(Y)

    # beta estimates
    fullBetas <- data.frame(matrix(nrow = q, ncol = p))
    row.names(fullBetas) <- sim_snps
    names(fullBetas) <- prots

    # SD of beta
    fullBetaSD <- data.frame(matrix(nrow = q, ncol = p))
    row.names(fullBetaSD) <- sim_snps
    names(fullBetaSD) <- prots

    # sample size
    fullBetaSS <- data.frame(matrix(nrow = q, ncol = p))
    row.names(fullBetaSS) <- sim_snps
    names(fullBetaSS) <- prots

    # p-value from GWAS
    fullP <- data.frame(matrix(nrow = q, ncol = p))
    row.names(fullP) <- sim_snps

    # calculate the summary stats
    for (i in 1:nrow(Y)) {
        sum.stats <- data.frame(matrix(nrow = nrow(X.std), ncol = 5))
        names(sum.stats) <- c("ID", "Beta", "SE", "P", "N")
        y <- as.numeric(Y[i, ])
        for (j in 1:nrow(X.std)) {
            x <- as.numeric(X.std[j, ])
            lm <- lm(y ~ x)
            sum.stats[j, 1] <- row.names(X.std)[j]
            sum.stats[j, 2] <- coefficients(lm)[2]
            sum.stats[j, 3] <- summary(lm)$coefficients[2, 2]
            sum.stats[j, 4] <- summary(lm)$coefficients[2, 4]
            sum.stats[j, 5] <- ncol(X.std)
        }
        # merged <- merge(sum.stats, sim_snps_bim, by = "ID") #merge change the order of SNPs
        fullBetas[, i] <- sum.stats$Beta
        fullBetaSD[, i] <- sum.stats$SE
        fullBetaSS[, i] <- sum.stats$N
        fullP[, i] <- sum.stats$P

        # write.table(merged, paste0("./sum_stats_temp/", row.names(Y)[i], "_intSS.txt"), quote = F, row.names = F, col.names = T)
    }
    return(list(fullBetas = fullBetas, fullBetaSD = fullBetaSD, fullBetaSS = fullBetaSS, fullP = fullP))
}

# estimate cormat using a set of null SNPs with p-value > 0.05
cal_corMat <- function(samplesize, Y, Xnull.std) {
    # read in chr5 SNPs as null SNPs
    # X.null <- read.table(paste0("./chr5.2200SNPs.sum.",samplesize,"_3.txt"))
    # Xnull.std <- t(apply(X.null,2,function(x) scale(x,scale=T, center=T)))
    # rm(X.null)

    extractZP <- function(x) {
        lm <- lm(y ~ x)
        Z.stats <- summary(lm)$coefficients[2, 3]
        P.stats <- summary(lm)$coefficients[2, 4]
        return(c(Z.stats, P.stats))
    }

    Z.stats <- data.frame(matrix(nrow = nrow(Xnull.std), ncol = nrow(Y)))
    P.stats <- data.frame(matrix(nrow = nrow(Xnull.std), ncol = nrow(Y)))
    for (i in 1:nrow(Y)) {
        y <- as.numeric(Y[i, ])
        result <- apply(Xnull.std, 1, function(x) extractZP(x))
        Z.stats[, i] <- result[1, ]
        P.stats[, i] <- result[2, ]
    }

    corMat.sim <- matrix(nrow = nrow(Y), ncol = nrow(Y))
    for (i in 1:nrow(Y)) {
        for (j in i:nrow(Y)) {
            IND <- P.stats[, i] > 0.05 & P.stats[, j] > 0.05
            # IND <- P.stats[,i] > 0 & P.stats[,j] > 0
            corMat.sim[i, j] <- cor(Z.stats[IND, i], Z.stats[IND, j])
            corMat.sim[j, i] <- corMat.sim[i, j]
        }
    }

    return(corMat.sim)
}

# scp /Users/yang/Documents/CCGG/network/network_Rachel/code/module/simulate_network.R yang3704@mesabi.msi.umn.edu:/home/yang3704/shared/network/
