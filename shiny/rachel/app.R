# Author: Rachel Zilinskas #
# Date: 05/24/2021 #
# Updated - 4/18/2022 to incorporate new filtering and choosing p-value threshold and method for V estimation #
# Updated 5/24/2022 to incorporate new MR analysis results
# Purpose: To create a shiny app to concisely display the results of my initial estimation of a protein
# DAG using summary statistics from Fokersen et al.

# Setting working directory - comment out when ready to publish
# setwd("/Users/chunlinli/Library/CloudStorage/Dropbox/sumdag/shiny/shinyapp")


# Sourcing helper functions
# source("helper.R")
# source("jointBetasHelper.R")
# source("DAGhelper.R")

# loading required packages
library(shiny)
library(stringr)
library(dagitty)
library(pcalg)
library(ggdag)
library(ggplot2)
library(stringr)
library(tidyr)
library(dplyr)
library(gt)
library(BiocManager)


# Definining the ancestry relationships to choose from in the graph
rel_choices <- c("descendants", "ancestors", "children", "neighbors", "spouses", "adjacent nodes", "markov blanket")

# Reading in a DAG to have access to the column names to choose from in the graph
U_MR <- as.matrix(read.table("p_value_mr.txt"))
U.TLP <- as.matrix(read.table("U_TLP.txt"))
p_value <- as.matrix(read.table("p_value.txt"))
p_value_mr <- as.matrix(read.table("p_value_mr.txt"))

# Define UI for application
ui <- fluidPage(
  tabsetPanel(
    tabPanel(
      "Real Data Results - TLP Regression",
      sidebarLayout(
        sidebarPanel(
          h3("To test the links in the DAG, we can conduct a Wald Test"),
          textOutput("num_links_tlp"),
          textInput("wald_p_thresh_tlp", "P-value threshold", value = 0.05)
        ),
        mainPanel(plotOutput("dag_wald_plot_tlp"))
      )
    ),
    tabPanel(
      "MR Results",
      sidebarLayout(
        sidebarPanel(
          h3("Visualizing the MR Results with AF as a node in the DAG"),
          textInput("mr_p_thresh", "P-value threshold", value = 0.05),
          selectInput("graph_prot", "Select a node", choices = colnames(U_MR), selected = "AD"),
          h4("We can further subset the graphs to visualize the estimated relationships in the DAG"),
          selectInput("relationship", "Select an ancestral relationship", choices = rel_choices, selected = "adjacent nodes")
        ),
        mainPanel(plotOutput("dag_plot_mr"))
      )
    )
  )
)

# Define server logic required to create output
server <- function(input, output) {
  # Making DAG plots using functions in helper file
  output$dag_wald_plot <- renderPlot({
    dag_wald_graph(input$wald_test, as.numeric(input$wald_p_thresh), input$v_est)
  })
  output$dag_wald_plot_tlp <- renderPlot({
    dag_wald_graph_tlp(U.TLP, p_value, input$wald_p_thresh_tlp)
  })
  output$num_links <- renderText({
    num_dag_links(input$v_est)
  })
  output$num_links_tlp <- renderText({
    num_dag_links_tlp(U.TLP)
  })
  output$dag_plot_mr <- renderPlot({
    dag_graph(input$graph_prot, input$relationship, input$mr_p_thresh)
  })
}


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

dag_wald_graph_tlp <- function(U, P.mat, p_thresh) {
  U.wald.tlp <- as.data.frame.matrix(ifelse(P.mat > p_thresh, 0, U))
  UDag <- pcalg2dagitty(t(U.wald.tlp), labels = names(U.wald.tlp), type = "dag")
  ggdag(UDag, node = T, text_col = "lightblue", stylized = F)
}

num_dag_links_tlp <- function(U_tlp) {
  nlinks <- sum(U_tlp != 0)
  nprots <- ncol(U_tlp)
  return(paste0("There are ", nprots, " proteins in the DAG. The number of estimated links is: ", nlinks))
}

dag_graph <- function(protein, relation, p_thresh) {
  U.mr <- as.data.frame.matrix(ifelse(p_value_mr > p_thresh, 0, U_MR))

  # Changing adjacency matrix to correct type
  daggity <- pcalg2dagitty(t(U.mr), labels = names(U.mr), type = "dag")
  if (relation == "descendants") {
    nodes <- descendants(daggity, protein)
  }
  if (relation == "ancestors") {
    nodes <- ancestors(daggity, protein)
  }
  if (relation == "children") {
    nodes <- c(protein, children(daggity, protein))
  }
  if (relation == "neighbors") {
    nodes <- c(protein, neighbours(daggity, protein))
  }
  if (relation == "spouses") {
    nodes <- c(protein, spouses(daggity, protein))
  }
  if (relation == "adjacent nodes") {
    nodes <- c(protein, adjacentNodes(daggity, protein))
  }
  if (relation == "markov blanket") {
    nodes <- c(protein, markovBlanket(daggity, protein))
  }
  if (length(nodes) <= 1) {
    df <- data.frame()
    df <- data.frame()
    DAG <- ggplot(df) +
      geom_point() +
      xlim(0, 10) +
      ylim(0, 10) +
      annotate("text", x = 5, y = 5, label = "This ancestral relationship does not exist for this protein")
  }
  if (length(nodes) > 1) {
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

# Run the application
shinyApp(ui = ui, server = server)
