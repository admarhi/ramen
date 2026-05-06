#' Plot a Directed Graph Emphasising Flow
#'
#' Arranges nodes of a directed graph into three distinct columns representing
#' sources (nodes with only outgoing edges), sinks (nodes with only incoming
#' edges), and intermediate nodes (with both incoming and outgoing edges).
#' The layout within the intermediate column is determined using the
#' Kamada-Kawai algorithm. This visualisation helps understand the overall
#' flow structure within the graph and is rendered with \pkg{ggraph}, so the
#' returned object composes with the usual \pkg{ggplot2} operators.
#'
#' Edges may be coloured in three ways: by a continuous \code{weight}
#' attribute (\code{color_edges_by_weight = TRUE}); by a categorical edge
#' attribute (\code{edge_color_attr = "<name>"}), in which case
#' \code{edge_color_values} and \code{edge_color_labels} customise the
#' palette and legend; or with a single fixed colour (default).
#'
#' @param g An igraph object. Must be a directed graph.
#' @param source_x Numeric scalar. The x-coordinate for source nodes.
#'   Defaults to 0.
#' @param mixed_x Numeric scalar. The central x-coordinate for intermediate
#'   nodes. The actual layout spans \code{mixed_x} +/- 0.4. Defaults to 1.
#' @param sink_x Numeric scalar. The x-coordinate for sink nodes. Defaults to 2.
#' @param vertical_spacing Numeric scalar. The maximum y-coordinate, controlling
#'   the vertical spread of the layout. Defaults to 1.
#' @param vertex_size Numeric scalar. The size of node points. Defaults to 4.
#' @param vertex_label_cex Numeric scalar. The size of node labels. Defaults
#'   to 3.
#' @param edge_arrow_size Numeric scalar. The arrow length on edges, in
#'   millimetres. Defaults to 2.
#' @param edge_width_range Numeric vector of length 2. The minimum and maximum
#'   width for edges when scaled by weight. If the graph is unweighted or all
#'   weights are identical, the mean of this range is used. Defaults to
#'   \code{c(0.3, 1.2)}.
#' @param color_edges_by_weight Logical scalar. If \code{TRUE} and the graph
#'   has a \code{weight} edge attribute, edges are coloured on a continuous
#'   gradient from \code{edge_color_low} to \code{edge_color_high}.
#' @param edge_color_attr Optional character scalar naming a categorical edge
#'   attribute to colour edges by (e.g. \code{"source"}). When set, takes
#'   precedence over \code{color_edges_by_weight}.
#' @param edge_color_values Optional named character vector mapping levels of
#'   \code{edge_color_attr} to colours. Used only when \code{edge_color_attr}
#'   is non-NULL.
#' @param edge_color_labels Optional named character vector mapping levels of
#'   \code{edge_color_attr} to legend labels. Used only when
#'   \code{edge_color_attr} is non-NULL.
#' @param edge_color_legend_title Character scalar. Legend title for the
#'   categorical edge colour scale. Defaults to \code{"Type"}.
#' @param edge_color_low Character string. The colour for the lowest edge
#'   weight when \code{color_edges_by_weight} is \code{TRUE}. Defaults to
#'   \code{"gray80"}.
#' @param edge_color_high Character string. The colour for the highest edge
#'   weight when \code{color_edges_by_weight} is \code{TRUE}. Defaults to
#'   \code{"black"}.
#' @param main Character string. Optional plot title.
#' @param ... Additional arguments. Currently unused; retained for backwards
#'   compatibility with the previous igraph-based implementation.
#'
#' @return A \code{ggplot} object.
#'
#' @export
#'
#' @examples
#' cm <- synCM("test", n_species = 3, max_met = 5)
#' g <- igraph::graph_from_adjacency_matrix(
#'     SummarizedExperiment::assays(cm)[["Binary"]],
#'     mode = "directed",
#'     weighted = TRUE
#' )
#' plotDirectedFlow(g)
#'
#' @importFrom igraph is_directed degree vcount induced_subgraph
#' @importFrom igraph layout_with_kk E ecount is_weighted V delete_edge_attr
#' @importFrom scales rescale
#' @importFrom ggraph ggraph create_layout geom_edge_link geom_node_point
#' @importFrom ggraph geom_node_text scale_edge_colour_gradient
#' @importFrom ggraph scale_edge_colour_identity scale_edge_colour_manual
#' @importFrom ggraph scale_edge_width_continuous circle
#' @importFrom ggraph guide_edge_colourbar
#' @importFrom ggplot2 aes arrow unit theme element_blank labs
#' @importFrom ggplot2 scale_colour_manual
#' @importFrom rlang .data
plotDirectedFlow <- function(
    g,
    source_x = 0,
    mixed_x = 1,
    sink_x = 2,
    vertical_spacing = 1,
    vertex_size = 4,
    vertex_label_cex = 3,
    edge_arrow_size = 2,
    edge_width_range = c(0.3, 1.2),
    color_edges_by_weight = FALSE,
    edge_color_attr = NULL,
    edge_color_values = NULL,
    edge_color_labels = NULL,
    edge_color_legend_title = "Type",
    edge_color_low = "gray80",
    edge_color_high = "black",
    main = NULL,
    ...
) {
    if (!igraph::is_directed(g)) {
        cli::cli_abort("{.arg g} must be a directed graph.")
    }

    ## ---- Layout: 3-column source / intermediate / sink ---------------------
    coords <- .threeColumnFlowLayout(
        g,
        source_x = source_x,
        mixed_x = mixed_x,
        sink_x = sink_x,
        vertical_spacing = vertical_spacing
    )

    ## ---- Node-role classification (also used for node colours) -------------
    in_deg <- igraph::degree(g, mode = "in")
    out_deg <- igraph::degree(g, mode = "out")
    role <- rep("intermediate", igraph::vcount(g))
    role[in_deg == 0 & out_deg > 0] <- "source"
    role[in_deg > 0 & out_deg == 0] <- "sink"
    role <- factor(role, levels = c("source", "intermediate", "sink"))

    ## ---- Build a manual ggraph layout and attach node aesthetics -----------
    lay <- ggraph::create_layout(
        g,
        layout = "manual",
        x = coords[, 1L],
        y = coords[, 2L]
    )
    lay$role <- role
    if (!"name" %in% names(lay)) {
        nm <- igraph::V(g)$name
        if (is.null(nm)) {
            nm <- as.character(seq_len(igraph::vcount(g)))
        }
        lay$name <- nm
    }

    ## ---- Edge weight handling ---------------------------------------------
    has_weight <- igraph::is_weighted(g)

    arrow_obj <- ggplot2::arrow(
        length = ggplot2::unit(edge_arrow_size, "mm"),
        type = "closed"
    )
    end_cap_obj <- ggraph::circle(vertex_size, "pt")
    start_cap_obj <- ggraph::circle(vertex_size, "pt")

    p <- ggraph::ggraph(lay)

    if (!is.null(edge_color_attr)) {
        ## ---- Categorical edge colours --------------------------------------
        if (!edge_color_attr %in% igraph::edge_attr_names(g)) {
            cli::cli_abort(
                "Edge attribute {.val {edge_color_attr}} not found on {.arg g}."
            )
        }
        ## Build an ad-hoc per-edge factor from the attribute values, ordered
        ## to match the names of the supplied palette where possible.
        attr_vals <- igraph::edge_attr(g, edge_color_attr)
        edge_levels <- if (!is.null(edge_color_values)) {
            names(edge_color_values)
        } else {
            sort(unique(attr_vals))
        }
        attr_fac <- factor(attr_vals, levels = edge_levels)

        ## Stash the factor as a column ggraph will see at edge-extraction
        ## time. The cleanest way is to set it as an edge attribute under a
        ## stable name and reference it in `aes()`.
        igraph::edge_attr(g, "._edge_colour_") <- attr_fac
        ## Re-create layout against the updated graph so the edge data carry
        ## the new attribute.
        lay <- ggraph::create_layout(
            g,
            layout = "manual",
            x = coords[, 1L],
            y = coords[, 2L]
        )
        lay$role <- role
        if (is.null(lay$name)) {
            nm <- igraph::V(g)$name
            if (is.null(nm)) {
                nm <- as.character(seq_len(igraph::vcount(g)))
            }
            lay$name <- nm
        }
        p <- ggraph::ggraph(lay) +
            ggraph::geom_edge_link(
                ggplot2::aes(edge_colour = .data$`._edge_colour_`),
                edge_width = mean(edge_width_range),
                arrow = arrow_obj,
                end_cap = end_cap_obj,
                start_cap = start_cap_obj
            )
        if (!is.null(edge_color_values)) {
            p <- p +
                ggraph::scale_edge_colour_manual(
                    values = edge_color_values,
                    labels = if (!is.null(edge_color_labels)) {
                        edge_color_labels
                    } else {
                        ggplot2::waiver()
                    },
                    name = edge_color_legend_title,
                    drop = FALSE
                )
        }
    } else if (color_edges_by_weight && has_weight) {
        eweight <- igraph::E(g)$weight
        if (length(unique(eweight)) == 1L) {
            p <- p +
                ggraph::geom_edge_link(
                    edge_colour = edge_color_high,
                    edge_width = mean(edge_width_range),
                    arrow = arrow_obj,
                    end_cap = end_cap_obj,
                    start_cap = start_cap_obj
                )
        } else {
            p <- p +
                ggraph::geom_edge_link(
                    ggplot2::aes(
                        edge_colour = .data$weight,
                        edge_width = .data$weight
                    ),
                    arrow = arrow_obj,
                    end_cap = end_cap_obj,
                    start_cap = start_cap_obj
                ) +
                ggraph::scale_edge_colour_gradient(
                    low = edge_color_low,
                    high = edge_color_high,
                    name = "Weight",
                    guide = ggraph::guide_edge_colourbar()
                ) +
                ggraph::scale_edge_width_continuous(
                    range = edge_width_range,
                    guide = "none"
                )
        }
    } else if (!is.null(igraph::E(g)$color)) {
        p <- p +
            ggraph::geom_edge_link(
                ggplot2::aes(edge_colour = .data$color),
                edge_width = mean(edge_width_range),
                arrow = arrow_obj,
                end_cap = end_cap_obj,
                start_cap = start_cap_obj
            ) +
            ggraph::scale_edge_colour_identity()
    } else {
        p <- p +
            ggraph::geom_edge_link(
                edge_colour = "gray50",
                edge_width = mean(edge_width_range),
                arrow = arrow_obj,
                end_cap = end_cap_obj,
                start_cap = start_cap_obj
            )
    }

    p <- p +
        ggraph::geom_node_point(
            ggplot2::aes(colour = .data$role),
            size = vertex_size,
            show.legend = FALSE
        ) +
        ggplot2::scale_colour_manual(
            values = c(
                source = "lightblue",
                intermediate = "lightgoldenrod",
                sink = "salmon"
            ),
            drop = FALSE
        ) +
        ggraph::geom_node_text(
            ggplot2::aes(label = .data$name),
            size = vertex_label_cex,
            repel = TRUE,
            max.overlaps = Inf
        ) +
        ggplot2::theme(
            panel.background = ggplot2::element_blank(),
            panel.grid = ggplot2::element_blank(),
            axis.text = ggplot2::element_blank(),
            axis.ticks = ggplot2::element_blank(),
            axis.title = ggplot2::element_blank()
        )

    if (!is.null(main)) {
        p <- p + ggplot2::labs(title = main)
    }

    p
}

#' Compute the 3-column source / intermediate / sink layout
#'
#' Internal helper. Returns a 2-column matrix of (x, y) coordinates, one row
#' per vertex of \code{g}, in vertex order. Sources occupy \code{source_x},
#' sinks \code{sink_x}, and intermediates are spread across
#' \code{mixed_x +/- 0.4} via Kamada-Kawai on the induced subgraph.
#'
#' @param g An igraph object.
#' @param source_x Numeric scalar.
#' @param mixed_x Numeric scalar.
#' @param sink_x Numeric scalar.
#' @param vertical_spacing Numeric scalar.
#'
#' @return Numeric matrix with \code{vcount(g)} rows and 2 columns.
#' @keywords internal
#' @noRd
.threeColumnFlowLayout <- function(
    g,
    source_x = 0,
    mixed_x = 1,
    sink_x = 2,
    vertical_spacing = 1
) {
    in_deg <- igraph::degree(g, mode = "in")
    out_deg <- igraph::degree(g, mode = "out")

    sources <- which(in_deg == 0 & out_deg > 0)
    sinks <- which(in_deg > 0 & out_deg == 0)
    mixed <- which(in_deg > 0 & out_deg > 0)
    isolated <- which(in_deg == 0 & out_deg == 0)

    layout <- matrix(NA_real_, nrow = igraph::vcount(g), ncol = 2L)

    if (length(mixed) > 1L) {
        subg <- igraph::induced_subgraph(g, vids = mixed)
        ## Kamada-Kawai requires strictly positive edge weights. If any
        ## edge has a non-positive weight, drop the weight attribute on the
        ## induced subgraph so layout_with_kk falls back to topological
        ## (unweighted) placement.
        if (igraph::is_weighted(subg)) {
            sub_weights <- igraph::E(subg)$weight
            if (any(!is.finite(sub_weights)) || any(sub_weights <= 0)) {
                subg <- igraph::delete_edge_attr(subg, "weight")
            }
        }
        mix_layout <- igraph::layout_with_kk(subg)
        mix_layout[, 1L] <- scales::rescale(
            mix_layout[, 1L],
            to = c(mixed_x - 0.4, mixed_x + 0.4)
        )
        mix_layout[, 2L] <- scales::rescale(
            mix_layout[, 2L],
            to = c(0, vertical_spacing)
        )
        layout[mixed, ] <- mix_layout
    } else if (length(mixed) == 1L) {
        layout[mixed, ] <- c(mixed_x, vertical_spacing / 2)
    }

    if (length(sources) > 0L) {
        y_sources <- seq(
            from = vertical_spacing,
            to = 0,
            length.out = max(1L, length(sources))
        )
        layout[sources, ] <- cbind(
            rep(source_x, length(sources)),
            y_sources
        )
    }

    if (length(sinks) > 0L) {
        y_sinks <- seq(
            from = vertical_spacing,
            to = 0,
            length.out = max(1L, length(sinks))
        )
        layout[sinks, ] <- cbind(
            rep(sink_x, length(sinks)),
            y_sinks
        )
    }

    if (length(isolated) > 0L) {
        y_iso <- seq(
            from = vertical_spacing * 1.2,
            to = vertical_spacing * 1.05,
            length.out = max(1L, length(isolated))
        )
        layout[isolated, ] <- cbind(
            rep(source_x, length(isolated)),
            y_iso
        )
    }

    layout
}
