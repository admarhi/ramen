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
        label_colours = NULL
    ) {
        # Extract edges from the ConsortiumMetabolismSet
        tb <- object@Edges

        # Create a tibble with all unique reactions per species
        rxns_per_species <- tb |>
            dplyr::select(1:3) |>
            unique() |>
            dplyr::mutate(
                edge = paste0(
                    .data$consumed,
                    "-",
                    .data$produced
                )
            ) |>
            dplyr::select("species", "edge")

        # Create a tibble of unique species with a row ID
        unique_species <- rxns_per_species |>
            dplyr::distinct(.data$species) |>
            tibble::rowid_to_column("ind")

        # Generate all unique combinations of species pairs
        species_combinations <- tidyr::expand_grid(
            x = unique_species$ind,
            y = unique_species$ind
        ) |>
            # Keep only unique pairs (x <= y)
            dplyr::filter(.data$x <= .data$y) |>
            dplyr::left_join(
                unique_species,
                by = c("x" = "ind")
            ) |>
            dplyr::left_join(
                unique_species,
                by = c("y" = "ind")
            ) |>
            dplyr::rename(
                species_x = "species.x",
                species_y = "species.y"
            ) |>
            dplyr::mutate(similarity = 0)

        # Calculate Jaccard similarity for each pair
        species_combinations$similarity <- mapply(
            function(sp_x, sp_y) {
                rxns_set_x <- rxns_per_species$edge[
                    rxns_per_species$species == sp_x
                ]
                rxns_set_y <- rxns_per_species$edge[
                    rxns_per_species$species == sp_y
                ]
                intersection <- length(
                    intersect(rxns_set_x, rxns_set_y)
                )
                union <- length(
                    union(rxns_set_x, rxns_set_y)
                )
                intersection / union
            },
            species_combinations$species_x,
            species_combinations$species_y
        )

        # Construct a sparse matrix from combinations
        similarity_matrix <- Matrix::sparseMatrix(
            species_combinations$x,
            species_combinations$y,
            x = species_combinations$similarity,
            dimnames = list(
                unique_species$species,
                unique_species$species
            ),
            symmetric = TRUE
        ) |>
            as.matrix()

        # Hierarchical clustering on distance matrix
        dend <- similarity_matrix |>
            stats::dist() |>
            stats::hclust(method = "complete") |>
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
                plot.background = ggplot2::element_blank(),
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
            species_combinations = species_combinations,
            reactions_per_species = rxns_per_species
        )

        # Print the plot
        print(plot_obj)

        # Return the result invisibly
        invisible(result)
    }
)
