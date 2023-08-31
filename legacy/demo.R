packages <- c("ncvreg","lassosum", "glmtlp","MASS")
invisible(lapply(packages, library, character.only = TRUE))

#setwd('/Users/yang/Documents/CCGG/network/network_Rachel/package/')
#setwd('/home/yang3704/yang3704/network/')
source('./test1/estimation.R')

#step 1: read in the data (SNP info including MAF is in IV_info.RData, but not needed here)
#fullBeta, fullBetaSD, and fullP should be correspond to the same set of SNPs
#beta
fullBeta <- read.table('./test/fullBeta_revised_ukbb.txt')
fullBetaSD <- read.table('./test/fullSE_revised_ukbb.txt')
#p-value
fullP <- read.table('./test/fullP_revised_ukbb.txt')
#sample size
fullBetaSS <- matrix(3393, nrow=nrow(fullBeta), ncol=ncol(fullBeta))

dim(fullP)
#[1] 33 23
count <- apply(fullP, 2, function(x) sum(x < 5e-8))
sum(count > 0)
#[1] 23

snp <- rownames(fullBeta)
prot <- colnames(fullBeta)
p <- length(prot)

#read in reference panel
load('./test/ukbb_xtx.RData')

#read in the phenotype correlation matrix
corMat1 = read.table('./test/23protein_phe_corr_allchr.txt', check.names=F)
#corMat1 = read.table('./test/23protein_phe_corr_allchr_revision.txt', check.names=F)
corMat = corMat1[prot,prot]


#step 2: estimate the quantities estimated from the summary stats
ss_back <- summary_stats_preprocess(fullbetas= fullBeta,
                                    fullses=fullBetaSD,
                                    N=fullBetaSS,
                                    xtxdiag=0,
                                    sds=NULL,
                                    geno.ref=xtx)

#step 3: estimate the V matrix
Vest <- estimateV_sum2(fullBeta, fullBetaSD, fullBetaSS, ss_back)
rownames(Vest) <- snp
colnames(Vest) <- prot

#step 4: peeling
peel_result <- peeling(V=Vest, thresh=0)
PiHat <- as.matrix(peel_result$PiHat)
PhiHat <- as.matrix(peel_result$PhiHat)

#step 5: TLP after peeling
result <- TLP_afterpeel_U2(PiHat, PhiHat,ss_back, as.matrix(corMat), fullBetaSS)
#U matrix
U.TLP <- result$U.est

#W matrix
W.TLP <- result$W.est


summary(abs(U.TLP[U.TLP!=0]))
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#0.04940 0.09294 0.12551 0.15448 0.20291 0.49412
 summary(abs(W.TLP[W.TLP!=0]))
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.1270  0.1841  0.3395  0.6600  0.6606  2.2289




