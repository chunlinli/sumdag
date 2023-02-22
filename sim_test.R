packages <- c("ncvreg","glmtlp","MASS")
invisible(lapply(packages, library, character.only = TRUE))

setwd('/home/yang3704/shared/network')
source('keyfunctions_testing.R')


#load simulation setting
Strength <- c(1,5,15)
for (kk in 1:3){
sim <-  sim1(ES=Strength[kk], Utxt='Uhat.txt',Wtxt='What1.txt')

set.seed(1994)
for (sss in 1:4){
samplesize = c(3000,6000,9000,12000)[sss]

#simulate data
Xdata <- input_data(samplesize.sim = samplesize, sim = sim, sigma2txt = 'sigma2.txt')
X.std = Xdata$X.std
sigma2 = Xdata$sigma2

#reference panel that used to estimate XtX
bfile <- "chr22.26SNPs.ref.3000_3.txt"

pvalue.null <- c()
pvalue.alter <- c()
pvalue.path.alter <- c()
pvalue.path.null <- c()
pvalue.edges.alter <- c()
for(k in 1:1000){	
	print(samplesize)
	#fix X and simulate Y
	Y <- simulate_Y(X.std, sigma2, sim, prots_name = Xdata$prots)
	
	#calculate the summary level data from the simulated individual level snp
	ss <- cal_sum_stats(X.std, Xdata$snp_names,  Xdata$sim_snps_bim, Y)
		
	corMat = cor(t(Y))

	ss_back <- summary_stats_preprocess(ss$fullBetas, ss$fullBetaSD, fullP=ss$fullP, N=ss$fullBetaSS, xtxdiag=0, ref=T, geno.ref.txt=bfile)
		
	#estimate V using the summary stats
	Vest <- estimateV_sum(ss, ss$fullBetaSS, ss_back)
	rownames(Vest) <- Xdata$snp_names
	colnames(Vest) <- Xdata$prots
	  
    #peeling algorithm 
    peel_result <- peeling(V=Vest, thresh=0.05)
    PiHat <- as.matrix(peel_result$PiHat)
    PhiHat <- as.matrix(peel_result$PhiHat)
    
    #covmat = diag(sqrt(ss_back$ytyVec)) %*% corMat %*% diag(sqrt(ss_back$ytyVec)) 
    covmat = corMat
    
    pairs = list(c(1,14))
    #test a single null edge
    pvalue.null[k] <- causal_inference(ss_back, PiHat, PhiHat, pairs, covmat, ss$fullBetaSS, test='edge')$p_value
    
	pairs = list(c(1,6))
    #test a single non-null edge
    pvalue.alter[k] <- causal_inference(ss_back, PiHat, PhiHat, pairs, covmat, ss$fullBetaSS, test='edge')$p_value
    
    #test two edges
    f_mat = PiHat
    pairs = list(c(7,15), c(1,6))
    pvalue.edges.alter[k] <- causal_inference(ss_back, PiHat, PhiHat, pairs, covmat, ss$fullBetaSS, test='edge')$p_value
    
    #testing a path (1)
    pairs = list(c(1,14), c(6,12))
    pvalue.path.null[k] <- causal_inference(ss_back, PiHat, PhiHat, pairs, covmat, ss$fullBetaSS, test='path')$p_value
    
    #print(pvalue.path.null)
  
 	#testing a path (2)
 	pairs = list(c(6,12), c(1,6))
    pvalue.path.alter[k] <- causal_inference(ss_back, PiHat, PhiHat, pairs, covmat, ss$fullBetaSS, test='path')$p_value
    
}
 
}
}

