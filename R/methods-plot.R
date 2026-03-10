#' @exportMethod plot
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

#' @exportMethod plot
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
        ggplot2::aes(
          x = .data$x,
          y = .data$y,
          xend = .data$xend,
          yend = .data$yend
        )
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
        # y - 2 to push below branches
        ggplot2::aes(x = .data$x, y = .data$y - 0.1, label = .data$label),
        angle = 90,
        hjust = 1,
        size = label_size,
        color = label_tb$colour # vector of colors
      )
  }
)

#' @exportMethod plot
setMethod(
    "plot",
    "ConsortiumMetabolismAlignment",
    function(x, type = NULL) {
        cli::cli_abort(
            "plot() for CMA will be reimplemented in Phase 3."
        )
    }
)
