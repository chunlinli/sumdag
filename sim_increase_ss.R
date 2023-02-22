packages <- c("ncvreg","glmtlp","MASS")
invisible(lapply(packages, library, character.only = TRUE))

setwd('/home/yang3704/shared/network')
source('keyfunctions.R')

#load simulation setting
sim <-  sim1.v(ES.U=1, ES.W=15, Utxt='Uhat.txt',Wtxt='What1.txt')

#reference panel that used to estimate XtX
bfile <- "chr22.26SNPs.ref.3000_3.txt"

replications <- 200
	V_falsePositives <- matrix(NA, nrow=4, ncol=replications)
	V_falseNegatives <- matrix(NA, nrow=4, ncol=replications)
	U_falsePositives <- matrix(NA, nrow=4, ncol=replications)
	U_falseNegatives <- matrix(NA, nrow=4, ncol=replications)

set.seed(1998)
for (sss in 1:4){
	samplesize = c(3000,6000,9000,12000)[sss]
	#simulate data
	Xdata <- input_data(samplesize.sim = samplesize, sim = sim, sigma2txt = 'sigma2.txt')
	X.std = Xdata$X.std
	sigma2 = Xdata$sigma2
	
	for (k in 1:replications){		
	print(samplesize)
	#fix X and simulate Y
	Y <- simulate_Y(X.std, sigma2, sim, prots_name = Xdata$prots)
	
	#calculate the summary level data from the simulated individual level snp
	ss <- cal_sum_stats(X.std, Xdata$snp_names,  Xdata$sim_snps_bim, Y)
	
	#use Y empricially here
	corMat = cor(t(Y))
	
	#estimate ssq back from summary stats 
	ss_back <- summary_stats_preprocess(ss$fullBetas, ss$fullBetaSD, fullP=ss$fullP, N=ss$fullBetaSS, xtxdiag=0, ref=T, geno.ref.txt=bfile)
	
	#estimate V using the summary stats
	Vest <- estimateV_sum(ss, ss$fullBetaSS, ss_back)
	Vest[abs(Vest) < 0] <- 0		
	rownames(Vest) <- Xdata$snp_names
	colnames(Vest) <- Xdata$prots
		
	#calculate the FP and FN for V
	trueV2 <- sim$Vsim
	trueV2[trueV2!=0] <- 1
	V2 <- Vest
	V2[V2!=0] <- 1
 	res <- trueV2 - V2
  
	V_falsePositives[sss,k] <- sum(res==-1)/sum(sim$Vsim==0)
	V_falseNegatives[sss,k] <- sum(res==1)/sum(sim$Vsim!=0)
	print(V_falsePositives[sss,k] )
	print(V_falseNegatives[sss,k])
	
	#update the summary statistics to exclude variables with no IVs
	trueU <- sim$U
	
    #peeling algorithm 
    peel_result <- peeling(V=Vest, thresh=0)
    PiHat <- peel_result$PiHat
    PhiHat <- peel_result$PhiHat
    
    U.TLP <- TLP_afterpeel_U(PiHat, PhiHat,ss_back, corMat, ss$fullBetaSS)$U.est
    colnames(U.TLP) <- rownames(U.TLP) <- Xdata$prots
   
    #estimate FP and FN for U
    U.TLP[U.TLP!=0] <- 1
    trueU[trueU!=0] <- 1
	waldRes <- trueU - U.TLP
    U_falsePositives[sss,k] <- sum(waldRes==-1)/sum(trueU==0)
  	U_falseNegatives[sss,k] <- sum(waldRes==1)/sum(trueU!=0)
  	print(U_falsePositives[sss,k])
  	print(U_falseNegatives[sss,k])
}

