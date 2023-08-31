library(dagitty)
library(pcalg)
library(ggdag)
library(ggplot2)
library(stringr)
library(tidyr)
library(dplyr)
library(gt)

# Sourcing necessary code
#setwd("/Users/Rachel/Documents/statGen/GRN/DAG/DAG_WithSS/ssIntDAG_ShinyApp/")

setwd("/Users/chunlinli/Library/CloudStorage/Dropbox/sumdag/shiny/shinyapp")

PiNames <- read.table("PiHat_001.txt")


prots <- row.names(PiNames)
yty_dat <- read.table("yty_med.txt", h = T, row.names = 1)
rownames(yty_dat) <- str_replace(rownames(yty_dat), "-", ".")
yty_dat2 <- data.frame(yty_dat[colnames(PiNames),])
names(yty_dat2) <- "yty_dat"
row.names(yty_dat2) <- colnames(PiNames)
s_dat <- read.table("s_med.txt", h = T, row.names = 1)
rownames(s_dat) <- str_replace(rownames(s_dat), "-", ".")
s_dat2 <- data.frame(s_dat[colnames(PiNames),])
names(s_dat2) <- "S"
row.names(s_dat2) <- colnames(PiNames)
cormat <- read.table("corMat.Names.txt")
rownames(cormat) <- str_replace(rownames(cormat), "-", ".")
cormat <- cormat[colnames(PiNames), colnames(PiNames)]
sigma <- read.table("ref_cov_001.txt", h = T)
fullBetas <- read.table("fullBetas_001.txt", h = T, row.names = 1)
fullBetaVar <- read.table("fullBetaVar_001.txt", h = T)
fullBetaSS <- read.table("fullBetaSS_001.txt", h = T)



yty_median <- c()
for(j in 1:ncol(fullBetas)){
  ytys <- c()
  betas <- fullBetas[,prots[j]]
  var <- fullBetaVar[,prots[j]]
  n <- fullBetaSS[,prots[j]]
  s <- diag(as.matrix(sigma))
  for(k in 1:length(betas)){
    ytys[k] <- (n[k]-1)*n[k]*s[k]*var[k] + n[k]*s[k]*(betas[k])^2
  }
  yty_median[j] <- median(ytys)
}
#yty_dat <- as.matrix(yty_median)
#row.names(yty_dat) <- prots
#colnames(yty_dat) <- "yty"

s_median <-c()
for(j in 1:ncol(fullBetas)){
  s <- c()
  betas <- fullBetas[,prots[j]]
  var <- fullBetaVar[,prots[j]]
  n <- fullBetaSS[,prots[j]]
  sig <- diag(as.matrix(sigma))
  for(k in 1:length(betas)){
    s[k] <- (n[k]-1)*var[k] + (betas[k])^2
  }
  s_median[j] <- median(s)
}
#s_dat <- as.matrix(s_median)
#row.names(s_dat) <- prots
#colnames(s_dat) <- "S"

dag_wald_graph <- function(wald, p_thresh, v_est){
  if(v_est == "Penalized Regression"){
    PhiHat <- read.table("PhiHat_001.txt")
    PiHat <- read.table("PiHat_001.txt")
  }
  if(v_est == "Joint Model; p = 0.05"){
    PhiHat <- read.table("PhiHat_001_jb1.txt")
    PiHat <- read.table("PiHat_001_jb1.txt")
  }
  if(v_est == "Joint Model; p = 0.05/125"){
    PhiHat <- read.table("PhiHat_001_jb2.txt")
    PiHat <- read.table("PiHat_001_jb2.txt")
  }
  if(v_est == "Joint Model; p = 0.05/(125*63)"){
    PhiHat <- read.table("PhiHat_001_jb3.txt")
    PiHat <- read.table("PiHat_001_jb3.txt")
  }
  if(wald == TRUE){
    fin <- DAG_solution_wald(PiHat, PhiHat, yty_dat2, s_dat2, cormat, sigma, fullBetas, fullBetaVar, fullBetaSS, p_thresh=p_thresh, n=max(fullBetaSS))
  }
  if(wald == FALSE){
    fin <- DAG_solution(PiHat, PhiHat, yty_dat2, s_dat2, cormat, sigma, fullBetas, fullBetaVar, fullBetaSS)
  }
  U <- fin$U
  UDag <- pcalg2dagitty(t(U), labels = names(U),type = "dag")
  ggdag(UDag, node = T, text_col = "lightblue", stylized = F)
}

dag_wald_graph_tlp <- function(U, P.mat, p_thresh){
  U.wald.tlp  = as.data.frame.matrix(ifelse(P.mat > p_thresh, 0, U))
  UDag <- pcalg2dagitty(t(U.wald.tlp), labels = names(U.wald.tlp), type = "dag")
  ggdag(UDag, node = T, text_col = "lightblue", stylized = F)
}

num_dag_links <- function(v_est){
  if(v_est == "Penalized Regression"){
    PhiHat <- read.table("PhiHat_001.txt")
    PiHat <- read.table("PiHat_001.txt")
  }
  if(v_est == "Joint Model; p = 0.05"){
    PhiHat <- read.table("PhiHat_001_jb1.txt")
    PiHat <- read.table("PiHat_001_jb1.txt")
  }
  if(v_est == "Joint Model; p = 0.05/125"){
    PhiHat <- read.table("PhiHat_001_jb2.txt")
    PiHat <- read.table("PiHat_001_jb2.txt")
  }
  if(v_est == "Joint Model; p = 0.05/(125*63)"){
    PhiHat <- read.table("PhiHat_001_jb3.txt")
    PiHat <- read.table("PiHat_001_jb3.txt")
  }
  fin <- DAG_solution(PiHat, PhiHat, yty_dat2, s_dat2, cormat, sigma, fullBetas, fullBetaVar, fullBetaSS)
  U <- fin$U
  nlinks <- sum(U!=0)
  nprots <- ncol(U)
  return(paste0("There are ", nprots, " proteins in the DAG. The number of estimated links is: " , nlinks))
}

num_dag_links_tlp <- function(U_tlp){
  nlinks <- sum(U_tlp!=0)
  nprots <- ncol(U_tlp)
  return(paste0("There are ", nprots, " proteins in the DAG. The number of estimated links is: " , nlinks))
}


dag_graph <- function( protein, relation, p_thresh){

  U.mr  = as.data.frame.matrix(ifelse(p_value_mr > p_thresh, 0, U_MR))

  # Changing adjacency matrix to correct type
  daggity <- pcalg2dagitty(t(U.mr), labels = names(U.mr), type = "dag")
  if(relation == "descendants"){
    nodes <- descendants(daggity, protein)
  }
  if(relation == "ancestors"){
    nodes <- ancestors(daggity, protein)
  }
  if(relation == "children"){
    nodes <- c(protein, children(daggity, protein))
  }
  if(relation == "neighbors"){
    nodes <- c(protein, neighbours(daggity, protein))
  }
  if(relation == "spouses"){
    nodes <- c(protein, spouses(daggity, protein))
  }
  if(relation == "adjacent nodes"){
    nodes <- c(protein, adjacentNodes(daggity, protein))
  }
  if(relation == "markov blanket"){
    nodes <- c(protein, markovBlanket(daggity, protein))
  }
  if(length(nodes)<=1){
    df <- data.frame()
    df <- data.frame()
    DAG <- ggplot(df) + geom_point() + xlim(0, 10) + ylim(0, 10) + annotate("text", x=5, y=5, label= "This ancestral relationship does not exist for this protein")

  }
  if(length(nodes)>1){
    sub <- U.mr[nodes, nodes]
    DAG <- pcalg2dagitty(t(sub), labels = names(sub), type = "dag") %>%
      tidy_dagitty() %>%
      dplyr::mutate(legend = ifelse(name == protein, "selected node", relation))


    DAG <- DAG %>%
      ggplot(aes(
        x = x,
        y = y,
        xend = xend,
        yend = yend
      )) +
      geom_dag_point(aes(colour = legend)) +
      geom_dag_edges() +
      geom_dag_text() +
      theme_dag()
  }
  DAG
}

dag_graph_tlp <- function(mr_analysis, protein, relation){
  if(mr_analysis == "Full GWAS data"){dag <- read.table("Uhat_tlp_AF_MR.txt")}
  if(mr_analysis == "pQTL meta-analysis"){dag <- read.table("Uhat_tlp_AF_META.txt")}

  # Changing adjacency matrix to correct type
  daggity <- pcalg2dagitty(t(dag), labels = names(dag), type = "dag")
  if(relation == "descendants"){
    nodes <- descendants(daggity, protein)
  }
  if(relation == "ancestors"){
    nodes <- ancestors(daggity, protein)
  }
  if(relation == "children"){
    nodes <- c(protein, children(daggity, protein))
  }
  if(relation == "neighbors"){
    nodes <- c(protein, neighbours(daggity, protein))
  }
  if(relation == "spouses"){
    nodes <- c(protein, spouses(daggity, protein))
  }
  if(relation == "adjacent nodes"){
    nodes <- c(protein, adjacentNodes(daggity, protein))
  }
  if(relation == "markov blanket"){
    nodes <- c(protein, markovBlanket(daggity, protein))
  }
  if(length(nodes)<=1){
    df <- data.frame()
    df <- data.frame()
    DAG <- ggplot(df) + geom_point() + xlim(0, 10) + ylim(0, 10) + annotate("text", x=5, y=5, label= "This ancestral relationship does not exist for this protein")

  }
  if(length(nodes)>1){
    sub <- dag[nodes, nodes]
    DAG <- pcalg2dagitty(t(sub), labels = names(sub), type = "dag") %>%
      tidy_dagitty() %>%
      dplyr::mutate(legend = ifelse(name == protein, "selected node", relation))


    DAG <- DAG %>%
      ggplot(aes(
        x = x,
        y = y,
        xend = xend,
        yend = yend
      )) +
      geom_dag_point(aes(colour = legend)) +
      geom_dag_edges() +
      geom_dag_text() +
      theme_dag()
  }
  DAG
}
