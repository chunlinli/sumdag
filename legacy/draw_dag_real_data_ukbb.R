#setwd('/Users/yang/Documents/CCGG/network/network_Rachel/code/realdata/')
#load('Protein23_revised_ukbb.RData')
U <- PiHat
colnames(U) <- gsub('[.]','-',colnames(U))

p_value[p_value >= 0.05/(23*22-sum(PiHat!=0))] = 0
pval <- p_value
res <- t(pval)[t(U)==1]

res <- ifelse(res==0,0,1)

source('./test1/draw_dag_code.r')
plot_graph(U, res,
           highlightnode=c("ADM","PGF","TNFRSF11B","TEK","MMP3","CHI3L1","IL1RL1","HAVCR1",
                           "MMP10","CXCL6","CXCL16","LGALS3","IL6R","CTSD"),
           underlightnode="RETN",
           graph_name='protein23_revised_ukbb.pdf')
