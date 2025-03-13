setMethod("plot", "ConsortiumMetabolism", function(x) {
  plot(x@Graphs[[1]])
})


setMethod("plot", "ConsortiumMetabolismAlignment", function(x, type = NULL) {
  # Plot network graph showing alignment of metabolisms
  if (is.null(type)) {
    # Get levels matrix and filter by max weight
    levels_mat <- assays(x)$Levels
    max_weight <- length(x@Graphs)
    levels_mat[levels_mat < max_weight] <- 0
    
    # Create graph from adjacency matrix
    g <- igraph::graph_from_adjacency_matrix(
      levels_mat,
      mode = "directed", 
      weighted = TRUE
    )
    
    # Remove isolated vertices
    g <- igraph::delete_vertices(
      g,
      igraph::V(g)[igraph::degree(g) == 0]
    )
    
    # Set edge visual properties
    edge_width <- igraph::E(g)$weight / max_weight * 4
    igraph::E(g)$width <- edge_width
    igraph::E(g)$arrow.size <- 0.5

    # Plot the graph
    plot(
      g,
      layout = igraph::layout_with_fr(g),
      edge.curved = 0.5,
      vertex.label.color = "black", 
      vertex.color = "lightblue",
      main = "Alignment of Consortia Metabolisms"
    )
  }
  # Plot Venn diagram of vertices
  else if (type == "venn-vertex") {
    vertex_list <- lapply(x@Graphs, igraph::V)
    ggVennDiagram::ggVennDiagram(vertex_list)
  }
  # Plot Venn diagram of edges  
  else if (type == "venn-edge") {
    edge_list <- lapply(x@Graphs, igraph::E)
    ggVennDiagram::ggVennDiagram(edge_list)
  }
})
