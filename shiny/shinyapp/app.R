
library(devtools)
devtools::install_github("chunlinli/glmtlp")
devtools::install_github("chunlinli/sumdag/sumdag")

library(glmtlp)
library(sumdag)

library(shiny)
library(igraph)

load("ppi.RData")
p <- beta_mat %>% ncol() # 23
protein <- colnames(beta_mat)
snp <- rownames(beta_mat)
vertex.color <- rep("#7a0019", p)
highlightcolor <- "#ffcc33"
underlightcolor <- "grey"
highlightnode <- c("ADM", "PGF", "TNFRSF11B", "TEK","MMP3", "CHI3L1", "IL1RL1", "HAVCR1","MMP10", "CXCL6", "CXCL16", "LGALS3","IL6R", "CTSD")
underlightnode <- "RETN"
vertex.color <- ifelse(protein %in% highlightnode, highlightcolor, vertex.color)
vertex.color <- ifelse(protein %in% underlightnode, underlightcolor, vertex.color)


update_plot <- function(g) {
    la <- layout.circle(g)

    plot(g,
        layout = layout.circle, vertex.label = "", main = "",
        vertex.color = vertex.color,
        vertex.frame.color = vertex.color,
        vertex.shape = "circle", # One of "none", "circle", "square", "csquare", "rectangle" "crectangle", "vrectangle", "pie", "raster", or "sphere"
        vertex.size = 10, # Size of the node (default is 15)

        # edge.color = "#000000", # Edge color
        edge.width = 0.5, # Edge width, defaults to 1
        edge.arrow.size = 0.5, # Arrow size, defaults to 1
        edge.arrow.width = 0.5, # Arrow width, defaults to 1
        # edge.lty = 2 # Line type, could be 0 or "blank", 1 or "solid", 2 or "dashed", 3 or "dotted", 4 or "dotdash", 5 or "longdash", 6 or "twodash"
    )

    # Apply labels manually
    # Specify x and y coordinates of labels, adjust outward as desired
    x <- la[, 1] * 1.2
    y <- la[, 2] * 1.2

    # Create vector of angles for text based on number of nodes
    # (flipping the orientation of the words half way around so none appear upside down)
    angle <- ifelse(atan(-(la[, 1] / la[, 2])) * (180 / pi) < 0,
        90 + atan(-(la[, 1] / la[, 2])) * (180 / pi),
        270 + atan(-la[, 1] / la[, 2]) * (180 / pi)
    )

    # Apply the text labels with a loop with angle as srt
    for (i in 1:length(x)) {
        text(
            x = x[i], y = y[i], labels = V(g)$name[i],
            adj = NULL, pos = NULL, cex = 1, col = "black",
            srt = angle[i], xpd = T
        )
    }
}


# Define UI for application
ui <- fluidPage(

    # Application title
    titlePanel("Inferring cardiovascular-related protein-protein interaction network"),

    # Sidebar layout with input and output definitions
    sidebarLayout(

        # Inputs
        sidebarPanel(
            helpText("Select the edges to test for significance."),

            # Edge selection
            selectInput("from_node", "From Node:", choices = protein),
            selectInput("to_node", "To Node:", choices = protein),

            # Add and Delete edge buttons
            actionButton("add_edge", "Add Edge"),
            actionButton("delete_edge", "Delete Edge"),
            br(),
            br(),
            helpText("Test the selected edge(s)."),

            # Calculate buttons
            actionButton("calculate", "Test selected edge(s)"),
            br(),
            br(),
            helpText("Test all edge."),
            actionButton("calculate_all", "Test all edges")
        ),

        # Outputs
        mainPanel(
            # Network plot
            plotOutput("network_plot"),

            # Output Panel
            uiOutput("significant_edges")
        )
    )
)






server <- function(input, output, session) {

    # Initialize an empty adjacency matrix and graph
    d_mat <- matrix(0, p, p)
    colnames(d_mat) <- rownames(d_mat) <- protein
    g <- graph_from_adjacency_matrix(d_mat, mode = "directed")

    output$network_plot <- renderPlot({
        # plot(g, edge.arrow.size = 0.5, edge.arrow.width = 0.5, edge.lty = 2) # edge.lty = 2 for dashed lines
        E(g)$lty <- 2
        update_plot(g)
    })


    observeEvent(input$add_edge, {
        # need to pass the d_mat into handle_edge function



        handle_edge(input$from_node, input$to_node, action = "add")
    })

    # Delete Edge
    observeEvent(input$delete_edge, {
        handle_edge(input$from_node, input$to_node, action = "delete")
    })

    handle_edge <- function(from, to, action) {
        if (action == "add") {
            # Update adjacency matrix
            if (from != to) {
                d_mat[from, to] <<- 1

                # Update graph
                g <<- add_edges(g, c(from, to))
            }
        } else if (action == "delete") {
            # Update adjacency matrix
            d_mat[from, to] <<- 0

            # Get edge id and delete the edge from the graph
            edge_id <- get.edge.ids(g, c(from, to))
            if (edge_id != -1) { # Edge exists
                g <<- delete_edges(g, edge_id)
            }
        }

        # Plot the updated graph
        output$network_plot <- renderPlot({
            E(g)$lty <- 2
            update_plot(g)
        })
    }

    # Format the significant edges for display
    format_significant_edges <- function(significant_edges) {
        edge_list <- seq(1, length(significant_edges), by = 2)
        num_edges <- length(edge_list)

        formatted_edges <- sapply(edge_list, function(i) {
            paste(significant_edges[i], "->", significant_edges[i + 1])
        })

        # Combine the number of edges and the list of edges into a single string


        output_text <- paste(
            ifelse(num_edges > 1,
                paste("There are", num_edges, "significant edges (displayed in solid and red arrows)<br/>"),
                paste("There is", num_edges, "significant edge (displayed in solid and red arrow)<br/>")
            ),
            paste(formatted_edges, collapse = "<br/>")
        )
        return(output_text)
    }


    # Perform Hypothesis Testing
    observeEvent(input$calculate, {
        g <- graph_from_adjacency_matrix(d_mat, mode = "directed")
        if (sum(d_mat != 0) == 0) {

            output$network_plot <- renderPlot({
                # plot(g, edge.arrow.size = 0.5, edge.arrow.width = 0.5)
                update_plot(g)
            })
            output$significant_edges <- renderUI({
                HTML("Please select at least one edge to test.")
            })
            return()
        }

        p_value <- matrix(1, nrow = nrow(d_mat), ncol = nrow(d_mat))
        pair <- list()
        for (i in 1:nrow(d_mat)) {
            for (j in 1:ncol(d_mat)) {
                if (d_mat[i, j] != 0) {
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

        p_value[p_value >= 0.05 / sum(d_mat != 0)] <- -1
        significant_edges_matrix <- which(p_value > 0, arr.ind = TRUE)

        vertex_names <- rownames(d_mat) # Assuming row and column names are the same and represent vertex names
        significant_edges <- mapply(function(r, c) {
            c(vertex_names[r], vertex_names[c])
        }, significant_edges_matrix[, "row"], significant_edges_matrix[, "col"])

        # Now, significant_edges should be a vector like c('A', 'B', 'C', 'D') for edges A->B and C->D

        all_edges <- get.edgelist(g)
        tested_edge_ids <- get.edge.ids(g, as.vector(t(all_edges)))

        significant_edge_ids <- get.edge.ids(g, significant_edges)



        E(g)$lty <- ifelse(tested_edge_ids %in% significant_edge_ids, 1, 2) # 1 for solid, 2 for dashed
        E(g)$color <- ifelse(tested_edge_ids %in% significant_edge_ids, "red", "black")



        # Plot the updated graph
        output$network_plot <- renderPlot({
            # plot(g, edge.arrow.size = 0.5, edge.arrow.width = 0.5)
            update_plot(g)
        })

        formatted_significant_edges <- format_significant_edges(significant_edges)

        output$significant_edges <- renderUI({
            HTML(formatted_significant_edges)
        })
    })

    observeEvent(input$calculate_all, {
        g <- graph_from_adjacency_matrix(an_mat, mode = "directed")
        p_value <- matrix(1, nrow = nrow(d_mat), ncol = nrow(d_mat))
        pair <- list()
        for (i in 1:nrow(an_mat)) {
            for (j in 1:ncol(an_mat)) {
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

        p_value[p_value >= 0.05 / (23 * 22 - sum(an_mat != 0))] <- -1
        significant_edges_matrix <- which(p_value > 0, arr.ind = TRUE)

        vertex_names <- rownames(d_mat) # Assuming row and column names are the same and represent vertex names
        significant_edges <- mapply(function(r, c) {
            c(vertex_names[r], vertex_names[c])
        }, significant_edges_matrix[, "row"], significant_edges_matrix[, "col"])


        all_edges <- get.edgelist(g)
        tested_edge_ids <- get.edge.ids(g, as.vector(t(all_edges)))

        significant_edge_ids <- get.edge.ids(g, significant_edges)



        E(g)$lty <- ifelse(tested_edge_ids %in% significant_edge_ids, 1, 2) # 1 for solid, 2 for dashed
        E(g)$color <- ifelse(tested_edge_ids %in% significant_edge_ids, "red", "black")


        # Plot the updated graph
        output$network_plot <- renderPlot({
            # plot(g, edge.arrow.size = 0.5, edge.arrow.width = 0.5)
            update_plot(g)
        })

        formatted_significant_edges <- format_significant_edges(significant_edges)

        output$significant_edges <- renderUI({
            HTML(formatted_significant_edges)
        })
    })
}


# Run the app ----
shinyApp(ui = ui, server = server)

