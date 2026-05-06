#' @include AllClasses.R AllGenerics.R theme-ramen.R plot-helpers.R
NULL

#' Plot a ConsortiumMetabolism object
#'
#' @description
#' Render a single assay of a \code{ConsortiumMetabolism} as a
#' directed metabolite-to-metabolite flow network. Each edge is a
#' pathway (consumed metabolite \eqn{\rightarrow} produced
#' metabolite); the chosen \code{type} determines what the edge
#' weight encodes.
#'
#' @param x A \code{ConsortiumMetabolism} object.
#' @param type Character. Which of the eight CM assays to render.
#'   One of:
#'   \describe{
#'     \item{\code{"Binary"}}{Presence / absence of each pathway
#'       (consumed \eqn{\rightarrow} produced).}
#'     \item{\code{"nSpecies"}}{Count of species that participate
#'       in each pathway.}
#'     \item{\code{"Consumption"}}{Total consumption flux per
#'       pathway (sum of uptake fluxes).}
#'     \item{\code{"Production"}}{Total production flux per pathway
#'       (sum of secretion fluxes).}
#'     \item{\code{"EffectiveConsumption"}}{Consumption flux scaled
#'       by Hill-1 perplexity, \eqn{F \cdot 2^{H(p)}}; same units
#'       as \code{Consumption} but combines magnitude with species
#'       evenness across the pathway.}
#'     \item{\code{"EffectiveProduction"}}{As above on the
#'       production side.}
#'     \item{\code{"nEffectiveSpeciesConsumption"}}{Hill-1
#'       effective number of consuming species. Unitless, in
#'       \eqn{[1, S]}.}
#'     \item{\code{"nEffectiveSpeciesProduction"}}{As above on the
#'       production side.}
#'   }
#' @return A \code{ggplot} object rendered with \pkg{ggraph}.
#' @examples
#' cm <- synCM("test", n_species = 3, max_met = 5)
#' plot(cm)
#' plot(cm, type = "EffectiveConsumption")
#' @exportMethod plot
setMethod(
    "plot",
    "ConsortiumMetabolism",
    function(
        x,
        type = c(
            "Binary",
            "nSpecies",
            "Consumption",
            "Production",
            "EffectiveConsumption",
            "EffectiveProduction",
            "nEffectiveSpeciesConsumption",
            "nEffectiveSpeciesProduction"
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
            colourEdgesByWeight = TRUE,
            edgeWidthRange = c(0.5, 1),
            title = x@Name,
            subtitle = paste(type, "network")
        )
    }
)

#' Plot a ConsortiumMetabolismSet object
#'
#' @description
#' Render the dendrogram of a \code{ConsortiumMetabolismSet}'s
#' hierarchical clustering of consortia, optionally annotated with
#' the internal cluster identifiers used by
#' \code{\link{extractCluster}}.
#'
#' @details
#' The dendrogram is the one computed at CMS construction time
#' from the FOS overlap matrix; tip labels are consortium names.
#'
#' When \code{showClusterIds = TRUE} (the default), small numbered
#' circles overlay the internal nodes of the dendrogram (up to
#' \code{max_nodes} of them). The numbers are the cluster IDs that
#' \code{\link{extractCluster}} accepts to retrieve a sub-CMS for
#' the subtree rooted at that node.
#'
#' For figure export, set \code{showClusterIds = FALSE} to obtain
#' a clean dendrogram without the cluster-ID overlays. Use
#' \code{label_colours} to recolour individual tip labels by
#' supplying a tibble with \code{label} and \code{colour} columns
#' (one row per consortium to recolour).
#'
#' @param x A \code{ConsortiumMetabolismSet} object.
#' @param label_colours Optional tibble with \code{label} and
#'   \code{colour} columns mapping consortium names to custom tip
#'   colours.
#' @param max_nodes Integer. Maximum number of internal cluster
#'   nodes to overlay with cluster-ID badges (only used when
#'   \code{showClusterIds = TRUE}). Defaults to \code{20}.
#' @param label_size Numeric tip-label size.
#' @param showClusterIds Logical. If \code{TRUE} (default), draw
#'   the internal cluster identifiers used by
#'   \code{\link{extractCluster}} as small filled circles on top
#'   of the dendrogram. Set to \code{FALSE} for a clean dendrogram
#'   suitable for figure export.
#' @return A \code{ggplot} object (returned invisibly).
#' @seealso \code{\link{extractCluster}}
#' @examples
#' cm1 <- synCM("comm_1", n_species = 3, max_met = 5)
#' cm2 <- synCM("comm_2", n_species = 4, max_met = 6)
#' cms <- ConsortiumMetabolismSet(cm1, cm2, name = "test")
#' plot(cms)
#' plot(cms, showClusterIds = FALSE)
#' @exportMethod plot
setMethod(
    "plot",
    "ConsortiumMetabolismSet",
    function(
        x,
        label_colours = NULL,
        max_nodes = 20,
        label_size = 3,
        showClusterIds = TRUE
    ) {
        if (length(x@Dendrogram) == 0) {
            cli::cli_abort(
                "The {.cls ConsortiumMetabolismSet} has not been clustered yet."
            )
        }
        dend <- x@Dendrogram[[1]]

        # nolint next: object_usage_linter.
        p <- .dendrogramGgplot(
            dend,
            labelSize = label_size,
            labelColours = label_colours,
            showHeightAxis = TRUE,
            title = if (length(name(x))) name(x) else NULL,
            subtitle = sprintf(
                "Hierarchical clustering · %d consortia",
                length(labels(dend))
            )
        )

        if (isTRUE(showClusterIds) && nrow(x@NodeData) > 0L) {
            nodeData <- x@NodeData |>
                dplyr::filter(.data$node_id <= max_nodes)
            p <- p +
                ggplot2::geom_point(
                    data = nodeData,
                    ggplot2::aes(x = .data$x, y = .data$y),
                    colour = "grey25",
                    fill = "grey25",
                    shape = 21,
                    size = 5
                ) +
                ggplot2::geom_text(
                    data = nodeData,
                    ggplot2::aes(
                        x = .data$x,
                        y = .data$y,
                        label = .data$node_id
                    ),
                    colour = "white",
                    size = 2.8,
                    fontface = "bold"
                )
        }
        p
    }
)

#' Plot a ConsortiumMetabolismAlignment object
#'
#' @param x A \code{ConsortiumMetabolismAlignment} object.
#' @param type Character specifying the plot type:
#'   \code{"heatmap"}, \code{"network"}, or
#'   \code{"scores"}.
#' @param ... Extra arguments forwarded to the underlying plot
#'   helper. For \code{type = "network"}, the most useful are
#'   \code{edgeColourValues} (named character vector mapping
#'   \code{shared}, \code{query}, \code{reference} to colours;
#'   defaults emphasise the shared edges with solid black against
#'   light/dark grey) and \code{nodeColourValues} (override the node
#'   role palette).
#' @return A \code{ggplot} object.
#' @examples
#' cm1 <- synCM("comm_1", n_species = 3, max_met = 5)
#' cm2 <- synCM("comm_2", n_species = 4, max_met = 6)
#' cma <- align(cm1, cm2)
#' plot(cma)
#' @exportMethod plot
#' @importFrom rlang .data
setMethod(
    "plot",
    "ConsortiumMetabolismAlignment",
    function(x, type = NULL, ...) {
        ## Default type based on alignment type
        if (is.null(type)) {
            type <- if (x@Type == "multiple") {
                "heatmap"
            } else {
                "network"
            }
        }
        type <- match.arg(
            type,
            c("heatmap", "network", "scores")
        )

        switch(
            type,
            heatmap = .plotHeatmap(x),
            network = .plotNetwork(x, ...),
            scores = .plotScores(x)
        )
    }
)

#' Plot similarity heatmap for multiple alignment
#' @noRd
.plotHeatmap <- function(cma) {
    if (cma@Type != "multiple") {
        cli::cli_abort(
            "Heatmap plot requires a multiple alignment."
        )
    }

    sim <- cma@SimilarityMatrix
    dend <- cma@Dendrogram[[1L]]

    ## Reorder by dendrogram leaf order
    leaf_order <- labels(dend)
    sim <- sim[leaf_order, leaf_order]

    ## Melt to long form
    n <- nrow(sim)
    idx <- expand.grid(
        row = seq_len(n),
        col = seq_len(n)
    )
    plot_df <- data.frame(
        consortium1 = factor(
            rownames(sim)[idx$row],
            levels = leaf_order
        ),
        consortium2 = factor(
            colnames(sim)[idx$col],
            levels = leaf_order
        ),
        similarity = as.vector(sim),
        stringsAsFactors = FALSE
    )

    p <- ggplot2::ggplot(
        plot_df,
        ggplot2::aes(
            x = .data$consortium1,
            y = .data$consortium2,
            fill = .data$similarity
        )
    ) +
        ggplot2::geom_tile(color = "white")

    if (n <= 30L) {
        ## Flip text colour for legibility on dark tiles.
        p <- p +
            ggplot2::geom_text(
                ggplot2::aes(
                    label = round(.data$similarity, 2),
                    colour = .data$similarity > 0.6
                ),
                size = 3,
                show.legend = FALSE
            ) +
            ggplot2::scale_colour_manual(
                values = c(`TRUE` = "white", `FALSE` = "grey15")
            )
    }

    p +
        ggplot2::scale_fill_gradient2(
            # nolint start: object_usage_linter.
            low = ramenPalette$heatmapFill[["low"]],
            mid = ramenPalette$heatmapFill[["mid"]],
            high = ramenPalette$heatmapFill[["high"]],
            # nolint end
            midpoint = 0.5,
            limits = c(0, 1),
            name = cma@Metric
        ) +
        ggplot2::coord_fixed() +
        theme_ramen() + # nolint: object_usage_linter.
        ggplot2::theme(
            axis.text.x = ggplot2::element_text(
                angle = 45,
                hjust = 1
            ),
            axis.title = ggplot2::element_blank(),
            panel.grid = ggplot2::element_blank()
        ) +
        ggplot2::labs(
            title = "Pairwise similarity",
            subtitle = sprintf(
                "%s · %d consortia",
                cma@Metric,
                n
            )
        )
}

#' Plot network for pairwise alignment
#' @noRd
.plotNetwork <- function(
    cma,
    # nolint next: object_usage_linter.
    edgeColourValues = ramenPalette$edgeCategorical,
    edgeColourLabels = c(
        shared = "Shared",
        query = "Query-only",
        reference = "Reference-only"
    ),
    edgeColourLegendTitle = "Pathway type",
    ...
) {
    if (cma@Type != "pairwise") {
        cli::cli_abort(
            "Network plot requires a pairwise alignment."
        )
    }

    ## Combine pathways with source labels
    pw_df <- data.frame(
        consumed = character(0L),
        produced = character(0L),
        source = character(0L),
        stringsAsFactors = FALSE
    )

    if (nrow(cma@SharedPathways) > 0L) {
        pw_df <- rbind(
            pw_df,
            data.frame(
                consumed = cma@SharedPathways$consumed,
                produced = cma@SharedPathways$produced,
                source = "shared",
                stringsAsFactors = FALSE
            )
        )
    }
    if (nrow(cma@UniqueQuery) > 0L) {
        pw_df <- rbind(
            pw_df,
            data.frame(
                consumed = cma@UniqueQuery$consumed,
                produced = cma@UniqueQuery$produced,
                source = "query",
                stringsAsFactors = FALSE
            )
        )
    }
    if (nrow(cma@UniqueReference) > 0L) {
        pw_df <- rbind(
            pw_df,
            data.frame(
                consumed = cma@UniqueReference$consumed,
                produced = cma@UniqueReference$produced,
                source = "reference",
                stringsAsFactors = FALSE
            )
        )
    }

    if (nrow(pw_df) == 0L) {
        cli::cli_warn("No pathways to plot.")
        return(invisible(NULL))
    }

    ## Build igraph; carry `source` as a categorical edge attribute so
    ## plotDirectedFlow can render a colour-blind-safe legend.
    edge_df <- pw_df[, c("consumed", "produced", "source")]
    g <- igraph::graph_from_data_frame(edge_df, directed = TRUE)

    fos <- if (!is.null(cma@Scores$FOS)) {
        sprintf("%s = %.3f", cma@Metric, cma@Scores$FOS)
    } else {
        cma@Metric
    }
    plotDirectedFlow(
        g,
        colourEdgesByWeight = FALSE,
        edgeColourAttr = "source",
        edgeColourValues = edgeColourValues,
        edgeColourLabels = edgeColourLabels,
        edgeColourLegendTitle = edgeColourLegendTitle,
        edgeWidthRange = c(1, 1),
        title = paste(cma@QueryName, "vs", cma@ReferenceName),
        subtitle = fos,
        ...
    )
}

#' Plot score comparison bar chart
#' @noRd
.plotScores <- function(cma) {
    s <- cma@Scores
    ## Filter out non-numeric and NA values
    vals <- vapply(
        s,
        function(v) {
            if (is.numeric(v) && length(v) == 1L) {
                v
            } else {
                NA_real_
            }
        },
        numeric(1L)
    )
    vals <- vals[!is.na(vals)]
    vals <- vals[!names(vals) %in% c("nPairs", "sd")]

    if (length(vals) == 0L) {
        cli::cli_abort("No scores available to plot.")
    }

    plot_df <- data.frame(
        metric = factor(
            names(vals),
            levels = rev(names(vals))
        ),
        value = unname(vals),
        stringsAsFactors = FALSE
    )

    ## Place value labels inside the bar when the bar is long enough,
    ## outside otherwise; flip text colour to match.
    insideBar <- plot_df$value > 0.7
    plot_df$labelHjust <- ifelse(insideBar, 1.1, -0.1)
    plot_df$labelColour <- ifelse(insideBar, "white", "grey15")

    ctx <- if (cma@Type == "pairwise") {
        sprintf("%s · %s vs %s", cma@Metric, cma@QueryName, cma@ReferenceName)
    } else {
        sprintf("%s · %d consortia", cma@Metric, nrow(cma@SimilarityMatrix))
    }

    ggplot2::ggplot(
        plot_df,
        ggplot2::aes(
            x = .data$metric,
            y = .data$value
        )
    ) +
        ggplot2::geom_col(
            fill = ramenPalette$bar, # nolint: object_usage_linter.
            width = 0.6
        ) +
        ggplot2::geom_text(
            ggplot2::aes(
                label = round(.data$value, 3),
                hjust = .data$labelHjust,
                colour = .data$labelColour
            ),
            size = 3,
            show.legend = FALSE
        ) +
        ggplot2::scale_colour_identity() +
        ggplot2::scale_y_continuous(
            limits = c(0, 1),
            breaks = seq(0, 1, 0.25),
            expand = ggplot2::expansion(mult = c(0, 0.05))
        ) +
        ggplot2::coord_flip() +
        theme_ramen() + # nolint: object_usage_linter.
        ggplot2::theme(
            panel.grid.major.y = ggplot2::element_blank()
        ) +
        ggplot2::labs(
            x = NULL,
            y = "Score",
            title = "Alignment scores",
            subtitle = ctx
        )
}
