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

setMethod(
  "plot",
  "ConsortiumMetabolismSet",
  function(x, label_colours = NULL, max_nodes = 20, label_size = 4) {
    if (length(x@Dendrogram) == 0) {
      stop("Not yet clustered!")
    }
    dend <- x@Dendrogram[[1]]
    node_data <- x@NodeData |>
      dplyr::filter(.data$node_id <= max_nodes)

    gg_dend <- dendextend::as.ggdend(dend, labels = TRUE, type = "rectangle")

    if (!is.null(label_colours)) {
      label_tb <- gg_dend$labels |>
        dplyr::left_join(label_colours, by = "label")
    } else {
      label_tb <- gg_dend$labels |>
        dplyr::mutate(colour = "black")
    }

    ggplot2::ggplot() +
      # Dendrogram branches
      ggplot2::geom_segment(
        data = gg_dend$segments,
        ggplot2::aes(x = x, y = y, xend = xend, yend = yend)
      ) +
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
        expand = ggplot2::expansion(mult = c(0, 0), add = c(2, 1))
      ) +
      ggplot2::theme(
        axis.text = ggplot2::element_blank(),
        axis.ticks = ggplot2::element_blank(),
        axis.title = ggplot2::element_blank(),
        panel.background = ggplot2::element_blank(),
        plot.background = ggplot2::element_blank(),
        panel.grid = ggplot2::element_blank()
      ) +
      # ggplot2::ggtitle(x@Name) +
      ggplot2::geom_text(
        # manually draw the axis labels (i.e., leaf labels)
        data = label_tb,
        ggplot2::aes(x = x, y = y - 0.1, label = label), # y - 2 to push below branches
        angle = 90,
        hjust = 1,
        size = label_size,
        color = label_tb$colour # vector of colors
      )
  }
)

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
