setMethod(
  "plot",
  "ConsortiumMetabolism",
  function(
    x,
    type = c(
      "Binary",
      "nEdges",
      "Consumption",
      "Production",
      "EffectiveConsumption",
      "EffectiveProduction"
    )
  ) {
    type <- match.arg(type)
    g <- igraph::graph_from_adjacency_matrix(
      assays(x)[[type]],
      mode = "directed",
      weighted = TRUE
    )

    plotDirectedFlow(
      g,
      color_edges_by_weight = TRUE,
      edge_width_range = c(0.5, 1),
      main = paste(type, "in", x@Name)
    )
  }
)

setMethod("plot", "ConsortiumMetabolismSet", function(x) {
  if (length(x@Dendrogram) == 0) {
    stop("Not yet clustered!")
  }
  dend <- x@Dendrogram[[1]]
  node_data <- x@NodeData

  gg_dend <- dendextend::as.ggdend(dend, labels = TRUE, type = "rectangle")
  ggplot2::ggplot(gg_dend) +
    ggplot2::geom_point(
      data = node_data,
      ggplot2::aes(x = .data$x, y = .data$y),
      color = "red",
      size = 7
    ) +
    ggplot2::geom_text(
      data = node_data,
      ggplot2::aes(x = .data$x, y = .data$y, label = .data$node_id),
      color = "white",
      size = 4,
      fontface = "bold"
    ) +
    ggplot2::scale_y_continuous(
      expand = c(0, 2)
    ) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(
        size = 2
      )
    )
})

setMethod("plot", "ConsortiumMetabolismAlignment", function(x, type = NULL) {
  # Plot network graph showing alignment of metabolisms
  if (is.null(type)) {
    # Get levels matrix and filter by max weight
    levels_mat <- assays(x)$Levels
    max_weight <- length(x@Consortia)

    # levels_mat[levels_mat < max_weight] <- 0
    colnames(levels_mat) <- SummarizedExperiment::colData(x)$met
    rownames(levels_mat) <- SummarizedExperiment::colData(x)$met

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
    edge_width <- igraph::E(g)$weight / max_weight * 2
    igraph::E(g)$width <- edge_width
    igraph::E(g)$arrow.size <- 0.2

    # Plot the graph
    plot(
      g,
      layout = igraph::layout_with_kk(g),
      # edge.curved = 0.5,
      vertex.label.color = "black",
      vertex.color = "lightblue",
      vertex.size = 4
      # main = "Alignment of Consortia Metabolisms"
    )
  } else if (type == "venn-vertex") {
    # Plot Venn diagram of vertices
    vertex_list <- lapply(x@Graphs, igraph::V)
    ggVennDiagram::ggVennDiagram(vertex_list)
  } else if (type == "venn-edge") {
    # Plot Venn diagram of edges
    edge_list <- lapply(x@Graphs, igraph::E)
    ggVennDiagram::ggVennDiagram(edge_list)
  }
})
