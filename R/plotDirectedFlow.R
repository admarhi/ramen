#' @include theme-ramen.R
NULL

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
#' attribute (\code{colourEdgesByWeight = TRUE}); by a categorical edge
#' attribute (\code{edgeColourAttr = "<name>"}), in which case
#' \code{edgeColourValues} and \code{edgeColourLabels} customise the
#' palette and legend; or with a single fixed colour (default).
#'
#' @details
#' The Kamada-Kawai layout used for the intermediate column is always
#' run unweighted. Within-column position has no semantic meaning
#' (edge weights are encoded via colour and width on the edges), and
#' running Kamada-Kawai with heterogeneous weights would collapse
#' strongly-connected intermediates into the centre on weighted assays
#' such as \code{Consumption} or \code{Production}. The source / sink /
#' intermediate column assignment is unaffected.
#'
#' The size-bearing arguments (\code{nodeSize}, \code{nodeLabelSize},
#' \code{edgeArrowSize}) use \pkg{ggraph} millimetre units rather than
#' the \pkg{igraph} \code{cex} factors of the previous implementation.
#'
#' @param g An igraph object. Must be a directed graph.
#' @param sourceX Numeric scalar. The x-coordinate for source nodes.
#'   Defaults to 0.
#' @param mixedX Numeric scalar. The central x-coordinate for intermediate
#'   nodes. The actual layout spans \code{mixedX} +/- 0.8. Defaults to 1.
#' @param sinkX Numeric scalar. The x-coordinate for sink nodes. Defaults
#'   to 2.
#' @param verticalSpacing Numeric scalar. The maximum y-coordinate,
#'   controlling the vertical spread of the layout. Defaults to 1.
#' @param nodeSize Numeric scalar. The size of node points in millimetres
#'   (\pkg{ggraph} unit). Defaults to 6.
#' @param nodeLabelSize Numeric scalar. The size of node labels in
#'   millimetres. Defaults to 3.
#' @param edgeArrowSize Numeric scalar. The arrow length on edges, in
#'   millimetres. Defaults to 2.
#' @param edgeWidthRange Numeric vector of length 2. The minimum and
#'   maximum width for edges when scaled by weight. If the graph is
#'   unweighted or all weights are identical, the mean of this range is
#'   used. Defaults to \code{c(0.3, 1.2)}.
#' @param colourEdgesByWeight Logical scalar. If \code{TRUE} and the graph
#'   has a \code{weight} edge attribute, edges are coloured on a continuous
#'   gradient from \code{edgeColourLow} to \code{edgeColourHigh}.
#' @param edgeColourAttr Optional character scalar naming a categorical
#'   edge attribute to colour edges by (e.g. \code{"source"}). When set,
#'   takes precedence over \code{colourEdgesByWeight}.
#' @param edgeColourValues Optional named character vector mapping levels
#'   of \code{edgeColourAttr} to colours. Used only when
#'   \code{edgeColourAttr} is non-NULL.
#' @param edgeColourLabels Optional named character vector mapping levels
#'   of \code{edgeColourAttr} to legend labels. Used only when
#'   \code{edgeColourAttr} is non-NULL.
#' @param edgeColourLegendTitle Character scalar. Legend title for the
#'   categorical edge colour scale. Defaults to \code{"Type"}.
#' @param edgeColourLow Character string. The colour for the lowest edge
#'   weight when \code{colourEdgesByWeight} is \code{TRUE}. Defaults to
#'   \code{ramenPalette$edgeWeight[["low"]]} (ColorBrewer Blues, light
#'   end).
#' @param edgeColourHigh Character string. The colour for the highest
#'   edge weight when \code{colourEdgesByWeight} is \code{TRUE}. Defaults
#'   to \code{ramenPalette$edgeWeight[["high"]]} (ColorBrewer Blues,
#'   dark end).
#' @param nodeColourValues Named character vector of length 3 mapping
#'   the node roles \code{"source"}, \code{"intermediate"}, and
#'   \code{"sink"} to colours. Defaults to \code{ramenPalette$nodeRole}.
#' @param nodeColourLabels Optional named character vector mapping node
#'   roles to legend labels. Defaults to \code{c(source = "Source",
#'   intermediate = "Intermediate", sink = "Sink")}.
#' @param nodeColourLegendTitle Character scalar. Legend title for the
#'   node-role colour scale. Defaults to \code{"Node role"}.
#' @param title Character string. Optional plot title. Falls back to
#'   the deprecated \code{main} alias when \code{title} is NULL.
#' @param subtitle Character string. Optional plot subtitle.
#' @param main Deprecated alias for \code{title}; will be removed in a
#'   future release.
#' @param ... Additional arguments. Currently unused; retained for forward
#'   compatibility.
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
#' @importFrom ggplot2 aes arrow unit theme theme_void element_text labs
#' @importFrom ggplot2 element_rect margin scale_colour_manual
#' @importFrom rlang .data
plotDirectedFlow <- function(
    g,
    sourceX = 0,
    mixedX = 1,
    sinkX = 2,
    verticalSpacing = 1,
    nodeSize = 6,
    nodeLabelSize = 3,
    edgeArrowSize = 2,
    edgeWidthRange = c(0.3, 1.2),
    colourEdgesByWeight = FALSE,
    edgeColourAttr = NULL,
    edgeColourValues = NULL,
    edgeColourLabels = NULL,
    edgeColourLegendTitle = "Type",
    # nolint start: object_usage_linter.
    edgeColourLow = ramenPalette$edgeWeight[["low"]],
    edgeColourHigh = ramenPalette$edgeWeight[["high"]],
    nodeColourValues = ramenPalette$nodeRole,
    # nolint end
    nodeColourLabels = c(
        source = "Source",
        intermediate = "Intermediate",
        sink = "Sink"
    ),
    nodeColourLegendTitle = "Node role",
    title = NULL,
    subtitle = NULL,
    main = NULL,
    ...
) {
    if (!igraph::is_directed(g)) {
        cli::cli_abort("{.arg g} must be a directed graph.")
    }
    if (is.null(title) && !is.null(main)) {
        title <- main
    }

    ## ---- Layout: 3-column source / intermediate / sink ---------------------
    coords <- .threeColumnFlowLayout(
        g,
        sourceX = sourceX,
        mixedX = mixedX,
        sinkX = sinkX,
        verticalSpacing = verticalSpacing
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
        length = ggplot2::unit(edgeArrowSize, "mm"),
        type = "closed"
    )
    end_cap_obj <- ggraph::circle(nodeSize, "pt")
    start_cap_obj <- ggraph::circle(nodeSize, "pt")

    p <- ggraph::ggraph(lay)

    if (!is.null(edgeColourAttr)) {
        ## ---- Categorical edge colours --------------------------------------
        if (!edgeColourAttr %in% igraph::edge_attr_names(g)) {
            cli::cli_abort(
                "Edge attribute {.val {edgeColourAttr}} not found on {.arg g}."
            )
        }
        ## Build an ad-hoc per-edge factor from the attribute values, ordered
        ## to match the names of the supplied palette where possible.
        attr_vals <- igraph::edge_attr(g, edgeColourAttr)
        edge_levels <- if (!is.null(edgeColourValues)) {
            names(edgeColourValues)
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
                edge_width = mean(edgeWidthRange),
                arrow = arrow_obj,
                end_cap = end_cap_obj,
                start_cap = start_cap_obj
            )
        if (!is.null(edgeColourValues)) {
            p <- p +
                ggraph::scale_edge_colour_manual(
                    values = edgeColourValues,
                    labels = if (!is.null(edgeColourLabels)) {
                        edgeColourLabels
                    } else {
                        ggplot2::waiver()
                    },
                    name = edgeColourLegendTitle,
                    drop = FALSE
                )
        }
    } else if (colourEdgesByWeight && has_weight) {
        eweight <- igraph::E(g)$weight
        if (length(unique(eweight)) == 1L) {
            p <- p +
                ggraph::geom_edge_link(
                    edge_colour = edgeColourHigh,
                    edge_width = mean(edgeWidthRange),
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
                    low = edgeColourLow,
                    high = edgeColourHigh,
                    name = "Weight",
                    guide = ggraph::guide_edge_colourbar()
                ) +
                ggraph::scale_edge_width_continuous(
                    range = edgeWidthRange,
                    guide = "none"
                )
        }
    } else if (!is.null(igraph::E(g)$color)) {
        p <- p +
            ggraph::geom_edge_link(
                ggplot2::aes(edge_colour = .data$color),
                edge_width = mean(edgeWidthRange),
                arrow = arrow_obj,
                end_cap = end_cap_obj,
                start_cap = start_cap_obj
            ) +
            ggraph::scale_edge_colour_identity()
    } else {
        p <- p +
            ggraph::geom_edge_link(
                edge_colour = "gray50",
                edge_width = mean(edgeWidthRange),
                arrow = arrow_obj,
                end_cap = end_cap_obj,
                start_cap = start_cap_obj
            )
    }

    ## Drop unused role rows from the legend so categories never
    ## represented in the data don't clutter it.
    presentRoles <- levels(droplevels(role))
    nodeColourValuesUsed <- nodeColourValues[
        names(nodeColourValues) %in% presentRoles
    ]
    nodeColourLabelsUsed <- nodeColourLabels[
        names(nodeColourLabels) %in% presentRoles
    ]

    p <- p +
        ggraph::geom_node_point(
            ggplot2::aes(colour = .data$role),
            size = nodeSize
        ) +
        ggplot2::scale_colour_manual(
            name = nodeColourLegendTitle,
            values = nodeColourValuesUsed,
            labels = nodeColourLabelsUsed,
            drop = TRUE,
            limits = presentRoles
        ) +
        ggraph::geom_node_text(
            ggplot2::aes(label = .data$name),
            size = nodeLabelSize,
            repel = TRUE,
            max.overlaps = Inf
        ) +
        ggplot2::coord_cartesian(clip = "off") +
        theme_ramen(network = TRUE) # nolint: object_usage_linter.

    if (!is.null(title) || !is.null(subtitle)) {
        p <- p + ggplot2::labs(title = title, subtitle = subtitle)
    }

    p
}

#' Compute the 3-column source / intermediate / sink layout
#'
#' Internal helper. Returns a 2-column matrix of (x, y) coordinates, one row
#' per vertex of \code{g}, in vertex order. Sources occupy \code{sourceX},
#' sinks \code{sinkX}, and intermediates are spread across
#' \code{mixedX +/- 0.8} via Kamada-Kawai on the induced subgraph.
#'
#' @param g An igraph object.
#' @param sourceX Numeric scalar.
#' @param mixedX Numeric scalar.
#' @param sinkX Numeric scalar.
#' @param verticalSpacing Numeric scalar.
#'
#' @return Numeric matrix with \code{vcount(g)} rows and 2 columns.
#' @keywords internal
#' @noRd
.threeColumnFlowLayout <- function(
    g,
    sourceX = 0,
    mixedX = 1,
    sinkX = 2,
    verticalSpacing = 1
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
        ## Lay out the intermediate column unweighted: within-column
        ## position is purely visual (edge weights are already encoded
        ## via colour and width). Heterogeneous weights would otherwise
        ## make Kamada-Kawai collapse strongly-connected nodes into the
        ## centre, producing an unreadable star on weighted assays
        ## (e.g. Consumption, Production).
        if (igraph::is_weighted(subg)) {
            subg <- igraph::delete_edge_attr(subg, "weight")
        }
        mix_layout <- igraph::layout_with_kk(subg)
        mix_layout[, 1L] <- scales::rescale(
            mix_layout[, 1L],
            to = c(mixedX - 0.8, mixedX + 0.8)
        )
        mix_layout[, 2L] <- scales::rescale(
            mix_layout[, 2L],
            to = c(0, verticalSpacing)
        )
        layout[mixed, ] <- mix_layout
    } else if (length(mixed) == 1L) {
        layout[mixed, ] <- c(mixedX, verticalSpacing / 2)
    }

    if (length(sources) > 0L) {
        y_sources <- seq(
            from = verticalSpacing,
            to = 0,
            length.out = max(1L, length(sources))
        )
        layout[sources, ] <- cbind(
            rep(sourceX, length(sources)),
            y_sources
        )
    }

    if (length(sinks) > 0L) {
        y_sinks <- seq(
            from = verticalSpacing,
            to = 0,
            length.out = max(1L, length(sinks))
        )
        layout[sinks, ] <- cbind(
            rep(sinkX, length(sinks)),
            y_sinks
        )
    }

    if (length(isolated) > 0L) {
        y_iso <- seq(
            from = verticalSpacing * 1.2,
            to = verticalSpacing * 1.05,
            length.out = max(1L, length(isolated))
        )
        layout[isolated, ] <- cbind(
            rep(sourceX, length(isolated)),
            y_iso
        )
    }

    layout
}
