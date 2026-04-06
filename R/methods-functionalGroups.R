#' @include AllClasses.R AllGenerics.R
NULL

#' @rdname functionalGroups
#' @aliases functionalGroups,ConsortiumMetabolismSet-method
#' @export
setMethod(
    "functionalGroups",
    "ConsortiumMetabolismSet",
    function(
        object,
        k = 4,
        label_size = 6,
        label_colours = NULL,
        linkage = "complete"
    ) {
        linkage <- match.arg(
            linkage,
            c("complete", "average", "single", "ward.D2")
        )
        # Extract pathways from the ConsortiumMetabolismSet
        tb <- object@Pathways

        # Build unique species-pathway pairs
        rxns_per_species <- tb |>
            dplyr::select(
                "consumed",
                "produced",
                "species"
            ) |>
            unique() |>
            dplyr::mutate(
                pathway = paste0(
                    .data$consumed,
                    "-",
                    .data$produced
                )
            ) |>
            dplyr::select("species", "pathway")

        # Build binary species-by-pathway incidence matrix
        sp_levels <- sort(unique(rxns_per_species$species))
        pw_levels <- unique(rxns_per_species$pathway)
        sp_idx <- match(
            rxns_per_species$species,
            sp_levels
        )
        pw_idx <- match(
            rxns_per_species$pathway,
            pw_levels
        )
        incidence <- Matrix::sparseMatrix(
            i = sp_idx,
            j = pw_idx,
            x = 1,
            dims = c(
                length(sp_levels),
                length(pw_levels)
            ),
            dimnames = list(sp_levels, pw_levels)
        )

        # Jaccard via crossprod:
        #   intersection = A %*% t(A)
        #   union = |A_i| + |A_j| - intersection
        intersection <- Matrix::tcrossprod(incidence)
        pathway_counts <- Matrix::rowSums(incidence)
        union_mat <- outer(
            pathway_counts,
            pathway_counts,
            "+"
        ) - as.matrix(intersection)
        similarity_matrix <- as.matrix(intersection) /
            union_mat
        # Fix 0/0 = NaN on diagonal or empty species
        similarity_matrix[is.nan(similarity_matrix)] <- 0
        diag(similarity_matrix) <- 1

        # Hierarchical clustering on distance matrix
        dend <- similarity_matrix |>
            stats::dist() |>
            stats::hclust(method = linkage) |>
            stats::as.dendrogram()

        # Color branches of the dendrogram
        gg_dend <- dendextend::color_branches(
            dend,
            k = k
        ) |>
            dendextend::as.ggdend()
        if (!is.null(label_colours)) {
            label_tb <- gg_dend$labels |>
                dplyr::left_join(
                    label_colours,
                    by = "label"
                )
        } else {
            label_tb <- gg_dend$labels |>
                dplyr::mutate(colour = "black")
        }
        # Remove the default label layer
        gg_dend$labels <- gg_dend$labels[0, ]

        # Make ggplot
        plot_obj <- ggplot2::ggplot(
            gg_dend,
            horiz = FALSE
        ) +
            ggplot2::theme_minimal() +
            ggplot2::theme(
                axis.title = ggplot2::element_blank(),
                axis.text = ggplot2::element_blank(),
                axis.ticks = ggplot2::element_blank(),
                axis.line = ggplot2::element_blank(),
                panel.background = ggplot2::element_blank(),
                panel.grid = ggplot2::element_blank()
            ) +
            ggplot2::scale_y_continuous(
                expand = ggplot2::expansion(
                    mult = c(0, 0),
                    add = c(2, 1)
                )
            ) +
            ggplot2::geom_text(
                data = label_tb,
                ggplot2::aes(
                    x = .data$x,
                    y = .data$y - 0.1,
                    label = .data$label
                ),
                angle = 90,
                hjust = 1,
                size = label_size,
                color = label_tb$colour
            )

        # Return list with plot, dendrogram, and data
        result <- list(
            plot = plot_obj,
            dendrogram = dend,
            similarity_matrix = similarity_matrix,
            incidence_matrix = incidence,
            reactions_per_species = rxns_per_species
        )

        # Return the result invisibly
        invisible(result)
    }
)
