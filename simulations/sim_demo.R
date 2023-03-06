#setwd('/home/yang3704/shared/network')
#source('simulate_network.R')

source("../simulations/estimation.R")
source("../simulations/testing.R")
source("../simulations/simulate_network.R")

#simulate data
sim <-  sim1(ES=5, Utxt='../simulations/Uhat.txt',Wtxt='../simulations/What1.txt')

bfile <- "../simulations/chr22.26SNPs.ref.3000_3.txt"
samplesize = 3000

Xdata <- input_data(samplesize.sim = samplesize, sim = sim, sigma2txt = '../simulations/sigma2.txt')
Y <- simulate_Y(Xdata$X.std, Xdata$sigma2, sim, prots_name = Xdata$prots)
cormat <-  cor(t(Y))

ss <- cal_sum_stats(Xdata$X.std, Xdata$snp_names,  Xdata$sim_snps_bim, Y)
fullBetas <-  ss$fullBetas
fullBetaSD <-  ss$fullBetaSD
N <- ss$fullBetaSS

ss_back <- summary_stats_preprocess(fullbetas = fullBetas, fullses =  fullBetaSD, N = N, xtxdiag=0, ref=TRUE, geno.ref.txt=bfile)

#constructing DAG
superDAG <- superDAG(fullBetas, fullBetaSD, N, xtxdiag=0,  ref=T, geno.ref.text=bfile, thresh=0.05)

PiHat <- superDAG$PiHat
PhiHat <- superDAG$PhiHat
result <- TLP_afterpeel_U(PiHat, PhiHat, ss_back, cormat, N, penalizeXX=T)

#LRT

pairs = list(c(1,14))
#test a single null edge
causal_inference(PiHat, PhiHat, ss_back, pairs, cormat, N, test='edge')$p_value

#test a single non-null edge
pairs = list(c(1,6))
causal_inference(PiHat, PhiHat, ss_back, pairs, cormat, ss$fullBetaSS, test='edge')$p_value


#testing a path (1)
pairs = list(c(1,14), c(6,12))
causal_inference(PiHat, PhiHat, ss_back,  pairs, cormat,N, test='path')$p_value

pairs = list(c(6,12), c(1,6))
causal_inference(PiHat, PhiHat, ss_back,  pairs, cormat, N, test='path')$p_value



