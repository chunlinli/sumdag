setwd("~/Documents/projects/sumdag/shiny/shinyapp")

library(shiny)
library(igraph)
library(sumdag)

load("ppi.RData")
p <- beta_mat %>% ncol() # 23
protein <- colnames(beta_mat)
snp <- rownames(beta_mat)


vertex.color <- rep("#7a0019", p)



# Define UI for application
ui <- fluidPage(

    # Application title
    titlePanel("Inferring cardiovascular-related protein-protein interaction network"),

    # Sidebar layout with input and output definitions
    sidebarLayout(

        # Inputs
        sidebarPanel(

            # Edge selection
            # The choices can be changed to any list of nodes
            selectInput("from_node", "From Node:", choices = protein),
            selectInput("to_node", "To Node:", choices = protein),

            # Add edge button
            actionButton("add_edge", "Add Edge"),

            # Type of testing
            radioButtons("test_type", "Type of Testing:",
                choices = c("edge", "path")
            ),

            # Calculate button
            actionButton("calculate", "Calculate")
        ),

        # Outputs
        mainPanel(
            # Network plot
            plotOutput("network_plot"),

            # Output Panel
            h3("Results"),
            textOutput("likelihood_ratio"),
            textOutput("degrees_freedom"),
            textOutput("p_value")
        )
    )
)



server <- function(input, output, session) {

    # Initialize an empty adjacency matrix and graph
    d_mat <- matrix(0, p, p)
    colnames(d_mat) <- rownames(d_mat) <- protein
    g <- graph_from_adjacency_matrix(d_mat, mode = "directed")
    la <- layout.circle(g)


    output$network_plot <- renderPlot({
        # plot(g, edge.arrow.size = 0.5, edge.arrow.width = 0.5, edge.lty = 2) # edge.lty = 2 for dashed lines

        plot(g,
            layout = layout.circle, vertex.label = "", main = "",

            # === vertex
            # vertex.color = "darkred", # Node color
            # vertex.frame.color = "darkred", # Node border color
            vertex.color = vertex.color,
            vertex.frame.color = vertex.color,
            vertex.shape = "circle", # One of "none", "circle", "square", "csquare", "rectangle" "crectangle", "vrectangle", "pie", "raster", or "sphere"
            vertex.size = 10, # Size of the node (default is 15)

            # === Edge
            edge.color = "#000000", # Edge color
            edge.width = 0.5, # Edge width, defaults to 1
            edge.arrow.size = 0.5, # Arrow size, defaults to 1
            edge.arrow.width = 0.5, # Arrow width, defaults to 1
            edge.lty = 2 # Line type, could be 0 or "blank", 1 or "solid", 2 or "dashed", 3 or "dotted", 4 or "dotdash", 5 or "longdash", 6 or "twodash"
        )


        ## Apply labels manually
        # Specify x and y coordinates of labels, adjust outward as desired
        x <- la[, 1] * 1.2
        y <- la[, 2] * 1.2

        # Create vector of angles for text based on number of nodes (flipping the orientation of the words half way around so none appear upside down)
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
    })

    # Add Edge
    observeEvent(input$add_edge, {
        from <- input$from_node
        to <- input$to_node

        # Update adjacency matrix
        if (from != to) {
            d_mat[from, to] <- 1

            # Update graph
            g <<- add_edges(g, c(from, to))
        }

        # Plot the graph
        output$network_plot <- renderPlot({
            # plot(g, edge.arrow.size = 0.5, edge.arrow.width = 0.5, edge.lty = 2) # edge.lty = 2 for dashed lines

            plot(g,
                layout = layout.circle, vertex.label = "", main = "",

                # === vertex
                # vertex.color = "darkred", # Node color
                # vertex.frame.color = "darkred", # Node border color
                vertex.color = vertex.color,
                vertex.frame.color = vertex.color,
                vertex.shape = "circle", # One of "none", "circle", "square", "csquare", "rectangle" "crectangle", "vrectangle", "pie", "raster", or "sphere"
                vertex.size = 10, # Size of the node (default is 15)

                # === Edge
                edge.color = "#000000", # Edge color
                edge.width = 0.5, # Edge width, defaults to 1
                edge.arrow.size = 0.5, # Arrow size, defaults to 1
                edge.arrow.width = 0.5, # Arrow width, defaults to 1
                edge.lty = 2 # Line type, could be 0 or "blank", 1 or "solid", 2 or "dashed", 3 or "dotted", 4 or "dotdash", 5 or "longdash", 6 or "twodash"
            )


            ## Apply labels manually
            # Specify x and y coordinates of labels, adjust outward as desired
            x <- la[, 1] * 1.2
            y <- la[, 2] * 1.2

            # Create vector of angles for text based on number of nodes (flipping the orientation of the words half way around so none appear upside down)
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
        })
    })

    # Perform Hypothesis Testing
    observeEvent(input$calculate, {

        # convert the nonzero entries in d_mat to pairs
        edges <- which(d_mat != 0, arr.ind = TRUE)

        # convert pairs to a list of edges
        pairs <- list()
        for (i in 1:nrow(edges)) {
            pairs[[i]] <- c(edges[i, 1], edges[i, 2])
        }



        # Assume causal_inference function returns a list of significant edges
        results <- sumdag::causal_inference(d_mat, test_type = input$test_type)

        # Update graph to highlight significant edges
        # Assuming results$significant_edges is a list of significant edges in the format c("from", "to")
        E(g)$lty <- ifelse(get.edge.ids(g, results$significant_edges), 1, 2) # edge.lty = 1 for solid lines
        E(g)$color <- ifelse(get.edge.ids(g, results$significant_edges), "red", "black")

        # Plot the updated graph
        output$network_plot <- renderPlot({
            #plot(g, edge.arrow.size = 0.5, edge.arrow.width = 0.5)




        })

        # Display the results
        output$likelihood_ratio <- renderText({
            paste("Likelihood Ratio:", results$likelihood_ratio)
        })
        output$degrees_freedom <- renderText({
            paste("Degrees of Freedom:", results$df)
        })
        output$p_value <- renderText({
            paste("P-value:", results$p_value)
        })
    })
}


# Run the app ----
shinyApp(ui = ui, server = server)
