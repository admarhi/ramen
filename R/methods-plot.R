#' @include AllClasses.R AllGenerics.R
NULL

#' Plot a ConsortiumMetabolism object
#'
#' @param x A \code{ConsortiumMetabolism} object.
#' @param type Character specifying the assay to plot.
#' @return A \code{ggplot} object (returned invisibly).
#' @examples
#' \donttest{
#' cm <- synCM("test", n_species = 3, max_met = 5)
#' plot(cm)
#' }
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

#' Plot a ConsortiumMetabolismSet object
#'
#' @param x A \code{ConsortiumMetabolismSet} object.
#' @param label_colours Optional tibble with label and
#'   colour columns.
#' @param max_nodes Maximum number of dendrogram nodes.
#' @param label_size Numeric label size.
#' @return A \code{ggplot} object (returned invisibly).
#' @examples
#' \donttest{
#' cm1 <- synCM("comm_1", n_species = 3, max_met = 5)
#' cm2 <- synCM("comm_2", n_species = 4, max_met = 6)
#' cms <- ConsortiumMetabolismSet(cm1, cm2, name = "test")
#' plot(cms)
#' }
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

        gg_dend <- dendextend::as.ggdend(
            dend,
            labels = TRUE,
            type = "rectangle"
        )

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
                ggplot2::aes(
                    x = .data$x,
                    y = .data$y - 0.1,
                    label = .data$label
                ),
                angle = 90,
                hjust = 1,
                size = label_size,
                color = label_tb$colour # vector of colors
            )
    }
)

#' Plot a ConsortiumMetabolismAlignment object
#'
#' @param x A \code{ConsortiumMetabolismAlignment} object.
#' @param type Character specifying the plot type:
#'   \code{"heatmap"}, \code{"network"}, or
#'   \code{"scores"}.
#' @return A \code{ggplot} object (returned invisibly).
#' @examples
#' \donttest{
#' cm1 <- synCM("comm_1", n_species = 3, max_met = 5)
#' cm2 <- synCM("comm_2", n_species = 4, max_met = 6)
#' cma <- align(cm1, cm2)
#' plot(cma)
#' }
#' @exportMethod plot
#' @importFrom rlang .data
setMethod(
    "plot",
    "ConsortiumMetabolismAlignment",
    function(x, type = NULL) {
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
            network = .plotNetwork(x),
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

    ggplot2::ggplot(
        plot_df,
        ggplot2::aes(
            x = .data$consortium1,
            y = .data$consortium2,
            fill = .data$similarity
        )
    ) +
        ggplot2::geom_tile(color = "white") +
        ggplot2::geom_text(
            ggplot2::aes(
                label = round(.data$similarity, 2)
            ),
            size = 3
        ) +
        ggplot2::scale_fill_gradient2(
            low = "#2166AC",
            mid = "white",
            high = "#B2182B",
            midpoint = 0.5,
            limits = c(0, 1),
            name = cma@Metric
        ) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
            axis.text.x = ggplot2::element_text(
                angle = 45,
                hjust = 1
            ),
            axis.title = ggplot2::element_blank(),
            panel.grid = ggplot2::element_blank()
        ) +
        ggplot2::ggtitle(
            paste("Similarity Heatmap -", cma@Metric)
        )
}

#' Plot network for pairwise alignment
#' @noRd
.plotNetwork <- function(cma) {
    if (cma@Type != "pairwise") {
        cli::cli_abort(
            "Network plot requires a pairwise alignment."
        )
    }

    ## Combine edges with source labels
    edges <- data.frame(
        consumed = character(0L),
        produced = character(0L),
        source = character(0L),
        stringsAsFactors = FALSE
    )

    if (nrow(cma@SharedPathways) > 0L) {
        edges <- rbind(
            edges,
            data.frame(
                consumed = cma@SharedPathways$consumed,
                produced = cma@SharedPathways$produced,
                source = "shared",
                stringsAsFactors = FALSE
            )
        )
    }
    if (nrow(cma@UniqueQuery) > 0L) {
        edges <- rbind(
            edges,
            data.frame(
                consumed = cma@UniqueQuery$consumed,
                produced = cma@UniqueQuery$produced,
                source = "query",
                stringsAsFactors = FALSE
            )
        )
    }
    if (nrow(cma@UniqueReference) > 0L) {
        edges <- rbind(
            edges,
            data.frame(
                consumed = cma@UniqueReference$consumed,
                produced = cma@UniqueReference$produced,
                source = "reference",
                stringsAsFactors = FALSE
            )
        )
    }

    if (nrow(edges) == 0L) {
        cli::cli_warn("No edges to plot.")
        return(invisible(NULL))
    }

    ## Build igraph
    g <- igraph::graph_from_data_frame(
        edges[, c("consumed", "produced")],
        directed = TRUE
    )

    ## Set edge colors by source
    color_map <- c(
        shared = "#4DAF4A",
        query = "#377EB8",
        reference = "#E41A1C"
    )
    igraph::E(g)$color <- color_map[edges$source]

    plotDirectedFlow(
        g,
        color_edges_by_weight = FALSE,
        edge_width_range = c(1, 1),
        main = paste(
            cma@QueryName,
            "vs",
            cma@ReferenceName
        )
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

    if (length(vals) == 0L) {
        cli::cli_abort("No scores available to plot.")
    }

    plot_df <- data.frame(
        metric = factor(
            names(vals),
            levels = names(vals)
        ),
        value = unname(vals),
        stringsAsFactors = FALSE
    )

    ggplot2::ggplot(
        plot_df,
        ggplot2::aes(
            x = .data$metric,
            y = .data$value
        )
    ) +
        ggplot2::geom_col(
            fill = "#377EB8",
            width = 0.6
        ) +
        ggplot2::geom_text(
            ggplot2::aes(
                label = round(.data$value, 3)
            ),
            hjust = -0.1,
            size = 3
        ) +
        ggplot2::coord_flip() +
        ggplot2::ylim(0, 1.1) +
        ggplot2::theme_minimal() +
        ggplot2::labs(
            x = NULL,
            y = "Score",
            title = paste("Scores -", cma@Metric)
        )
}
