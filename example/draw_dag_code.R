################# PLOT

library(igraph)

'#F8766D'
'#00BFC4'

plot_graph <- function(E, res, highlightcolor='#ffcc33', highlightnode=NULL, underlightnode=NULL, graph_name = "graph-ad.pdf") {
  p <- ncol(E)
  par(mfrow = c(1, 2))
  g <- graph_from_adjacency_matrix(E, mode = "directed", diag = F)
  V(g)$name <- colnames(E)
  nedge = sum(E!=0)

  ## laid out as a circle to begin with
  la <- layout.circle(g)

  pdf(
    file = graph_name, # The file name you want to save the plot in
    width = 7, # The width & height of the plot in inches
    height = 7
  )

  vertex.color = rep('#7a0019', p)
  if (!is.null(highlightnode)) {
    		vertex.color = ifelse(colnames(E) %in% highlightnode, highlightcolor, vertex.color)
    }
    if (!is.null(underlightnode)) {
    		vertex.color = ifelse(colnames(E) %in% underlightnode, 'grey', vertex.color)
    }

  plot(g,
    layout = layout.circle, vertex.label = "", main = " ",

    # === vertex
    #vertex.color = "darkred", # Node color
    #vertex.frame.color = "darkred", # Node border color
    vertex.color = vertex.color,
    vertex.frame.color=vertex.color,
    vertex.shape = "circle", # One of "none", "circle", "square", "csquare", "rectangle" "crectangle", "vrectangle", "pie", "raster", or "sphere"
    vertex.size = 10, # Size of the node (default is 15)

    # === Edge
    edge.color="#000000",                           # Edge color
    edge.width = rep(0.7, nedge) + res * 0.3, # Edge width, defaults to 1
    edge.arrow.size = rep(0.7, nedge) + res * 0.3, # Arrow size, defaults to 1
    edge.arrow.width = rep(0.7, nedge) + res * 0.3, # Arrow width, defaults to 1
    edge.lty = rep(2, nedge) - res # Line type, could be 0 or "blank", 1 or "solid", 2 or "dashed", 3 or "dotted", 4 or "dotdash", 5 or "longdash", 6 or "twodash"
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

  dev.off()
}
