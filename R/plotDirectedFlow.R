#' Plot a Directed Graph Emphasizing Flow
#'
#' Arranges nodes of a directed graph into three distinct columns representing
#' sources (nodes with only outgoing edges), sinks (nodes with only incoming
#' edges), and intermediate nodes (with both incoming and outgoing edges).
#' The layout within the intermediate column is determined using the
#' Fruchterman-Reingold algorithm. This visualization helps understand the
#' overall flow structure within the graph.
#'
#' @param g An igraph object. Must be a directed graph.
#' @param source_x Numeric scalar. The x-coordinate for source nodes.
#'   Defaults to 0.
#' @param mixed_x Numeric scalar. The central x-coordinate for intermediate
#'   nodes. The actual layout spans \code{mixed_x} +/- 0.4. Defaults to 1.
#' @param sink_x Numeric scalar. The x-coordinate for sink nodes. Defaults to 2.
#' @param vertical_spacing Numeric scalar. The maximum y-coordinate, controlling
#'   the vertical spread of the layout. Defaults to 1.
#' @param vertex_size Numeric scalar. The size of the vertices (nodes) in the
#'   plot. Defaults to 10. 
#' @param vertex_label_cex Numeric scalar. The character expansion factor for
#'   vertex labels. Defaults to 0.8.
#' @param edge_arrow_size Numeric scalar. The size of the arrows on the edges.
#'   Defaults to 0.5.
#' @param edge_width_range Numeric vector of length 2. The minimum and maximum
#'   width for edges when scaled by weight. If the graph is unweighted or all
#'   weights are identical, the mean of this range is used. Defaults to
#'   \code{c(0.1, 0.1)}.
#' @param color_edges_by_weight Logical scalar.
#'   If \code{TRUE} and the graph has a
#'   'weight' edge attribute, edges will be colored based on their weight,
#'   interpolating between \code{edge_color_low} and \code{edge_color_high}.
#' @param edge_color_low Character string. The color for the lowest edge weight
#'   when \code{color_edges_by_weight} is \code{TRUE}. Defaults to "gray80".
#' @param edge_color_high Character string. The color for the highest edge
#' weight when \code{color_edges_by_weight} is \code{TRUE}. Defaults to "black".
#' @param ... Additional arguments passed to igraph.
#'
#' @return Invisibly returns \code{NULL}. This function is called for its
#'   side effect of generating a plot.
#'
#' @export
#'
#' @importFrom igraph is_directed degree vcount induced_subgraph
#' @importFrom igraph layout_with_fr E ecount is_weighted plot.igraph
#' @importFrom scales rescale gradient_n_pal
#' @importFrom graphics plot
plotDirectedFlow <- function(
  g,
  source_x = 0,
  mixed_x = 1,
  sink_x = 2,
  vertical_spacing = 1,
  vertex_size = 10,
  vertex_label_cex = 0.8,
  edge_arrow_size = 0.5,
  edge_width_range = c(0.1, 0.1),
  color_edges_by_weight = FALSE,
  edge_color_low = "gray80",
  edge_color_high = "black",
  ...
) {
  if (!igraph::is_directed(g)) {
    stop("Graph must be directed.")
  }

  # Compute degrees
  in_deg <- igraph::degree(g, mode = "in")
  out_deg <- igraph::degree(g, mode = "out")

  # Classify nodes
  sources <- which(in_deg == 0 & out_deg > 0)
  sinks <- which(in_deg > 0 & out_deg == 0)
  mixed <- which(in_deg > 0 & out_deg > 0)

  layout <- matrix(NA, nrow = igraph::vcount(g), ncol = 2)

  ## Layout for mixed nodes
  if (length(mixed) > 1) {
    subg <- igraph::induced_subgraph(g, vids = mixed)
    mix_layout <- igraph::layout_with_kk(subg)

    mix_layout[, 1] <- scales::rescale(
      mix_layout[, 1],
      to = c(mixed_x - 0.4, mixed_x + 0.4)
    )
    mix_layout[, 2] <- scales::rescale(
      mix_layout[, 2],
      to = c(0, vertical_spacing)
    )

    layout[mixed, ] <- mix_layout
  } else if (length(mixed) == 1) {
    layout[mixed, ] <- c(mixed_x, vertical_spacing / 2)
  }

  ## Sources
  y_sources <- seq(
    from = vertical_spacing,
    to = 0,
    length.out = max(1, length(sources))
  )
  layout[sources, ] <- cbind(rep(source_x, length(sources)), y_sources)

  ## Sinks
  y_sinks <- seq(
    from = vertical_spacing,
    to = 0,
    length.out = max(1, length(sinks))
  )
  layout[sinks, ] <- cbind(rep(sink_x, length(sinks)), y_sinks)

  ## Edge weights → width
  if (igraph::is_weighted(g)) {
    w <- igraph::E(g)$weight
    if (length(unique(w)) == 1) {
      edge_widths <- rep(mean(edge_width_range), length(w))
    } else {
      edge_widths <- scales::rescale(w, to = edge_width_range)
    }
  } else {
    edge_widths <- rep(mean(edge_width_range), igraph::ecount(g))
  }

  ## Edge color mapping
  if (color_edges_by_weight && igraph::is_weighted(g)) {
    w <- igraph::E(g)$weight
    if (length(unique(w)) > 1) {
      edge_color_scaled <- scales::rescale(w, to = c(0, 1))
    } else {
      edge_color_scaled <- rep(0.5, length(w))
    }
    edge_colors <- scales::gradient_n_pal(c(edge_color_low, edge_color_high))(
      edge_color_scaled
    )
  } else {
    edge_colors <- rep("gray50", igraph::ecount(g))
  }

  ## Node colors by role
  colors <- rep("gray", igraph::vcount(g))
  colors[sources] <- "lightblue"
  colors[mixed] <- "lightgoldenrod"
  colors[sinks] <- "salmon"

  ## Plot
  graphics::plot(
    g,
    layout = layout,
    vertex.color = colors,
    vertex.size = vertex_size,
    vertex.label.cex = vertex_label_cex,
    edge.arrow.size = edge_arrow_size,
    edge.width = edge_widths,
    edge.color = edge_colors,
    ...
  )
}
