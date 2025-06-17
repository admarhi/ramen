#' @rdname getFunctionalGroups
#' @aliases getFunctionalGroups,ConsortiumMetabolismSet-method
#' @export
setMethod(
  "getFunctionalGroups",
  "ConsortiumMetabolismSet",
  function(object, k = 4, label_size = 6, label_colours = NULL) {
    # Extract edges from the ConsortiumMetabolismSet object
    tb <- object@Edges

    # Create a tibble with all unique reactions for each species
    # This identifies which species performs which metabolic reactions
    rxns_per_species <- tb |>
      dplyr::select(1:3) |>
      unique() |>
      dplyr::mutate(edge = paste0(.data$consumed, "-", .data$produced)) |>
      dplyr::select("species", "edge")

    # Create a tibble of unique species with a row ID for easier indexing
    unique_species <- rxns_per_species |>
      dplyr::distinct(.data$species) |>
      tibble::rowid_to_column("ind")

    # Generate all unique combinations of species pairs
    # This is to calculate similarity between each pair of species
    species_combinations <- tidyr::expand_grid(
      x = unique_species$ind,
      y = unique_species$ind
    ) |>
      # Keep only unique pairs (x <= y) to avoid redundant calculations
      dplyr::filter(.data$x <= .data$y) |>
      # Join with unique_species to get species names for x and y
      dplyr::left_join(unique_species, by = c("x" = "ind")) |>
      dplyr::left_join(unique_species, by = c("y" = "ind")) |>
      janitor::clean_names() |>
      # Initialize a column to store similarity scores
      dplyr::mutate(similarity = 0)

    # Calculate Jaccard similarity for each pair of species
    # Jaccard similarity = (intersection of reactions) / (union of reactions)
    for (i in seq_len(nrow(species_combinations))) {
      species_x <- species_combinations$species_x[i]
      species_y <- species_combinations$species_y[i]

      # Get the set of reactions for species x and species y
      rxns_set_x <- rxns_per_species$edge[rxns_per_species$species == species_x]
      rxns_set_y <- rxns_per_species$edge[rxns_per_species$species == species_y]

      # Calculate intersection and union of reaction sets
      intersection <- length(intersect(rxns_set_x, rxns_set_y))
      union <- length(union(rxns_set_x, rxns_set_y))

      # Calculate and store Jaccard similarity
      species_combinations$similarity[i] <- intersection / union
    }

    # Construct a sparse matrix from the species combinations and similarities
    # This matrix represents the similarity between all pairs of species
    similarity_matrix <- Matrix::sparseMatrix(
      species_combinations$x,
      species_combinations$y,
      x = species_combinations$similarity,
      dimnames = list(unique_species$species, unique_species$species),
      symmetric = TRUE # Add symmetric = TRUE as the matrix is symmetric
    ) |>
      as.matrix()

    # Perform hierarchical clustering on the distance matrix derived from similarities
    # The result is a dendrogram representing functional groups
    dend <- similarity_matrix |>
      # Convert similarity to distance (1 - similarity could be an option if needed)
      stats::dist() |> # Use stats::dist for clarity
      stats::hclust(method = "complete") |> # Use stats::hclust for clarity
      stats::as.dendrogram() # Use stats::as.dendrogram for clarity

    # Color branches of the dendrogram for visualization (e.g., k=4 clusters)
    # The plot is displayed to the user
    # dend |>
    #   dendextend::color_branches(k = k) |>
    #   # ggplot2::ggplot() |>
    #   plot()
    gg_dend <- dendextend::color_branches(dend, k = k) |>
      dendextend::as.ggdend()
    if (!is.null(label_colours)) {
      label_tb <- gg_dend$labels |>
        dplyr::left_join(label_colours, by = "label")
    } else {
      # Make it such that otherwise it takes the leaves colour
      label_tb <- gg_dend$labels |>
        dplyr::mutate(colour = "black")
    }
    # Remove the default label layer
    gg_dend$labels$label <- NA

    # Make ggplot
    ggplot2::ggplot(gg_dend, horiz = FALSE) +
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
        expand = ggplot2::expansion(mult = c(0, 0), add = c(2, 1))
      ) +
      ggplot2::geom_text(
        # manually draw the axis labels (i.e., leaf labels)
        data = label_tb,
        ggplot2::aes(x = x, y = y - 0.1, label = label), # y - 2 to push below branches
        angle = 90,
        hjust = 1,
        size = label_size,
        color = label_tb$colour # vector of colors
      )
    ### Need to output also the tibbles with the data long term
    # Return the dendrogram object invisibly
    # This allows the user to assign it to a variable if desired
    # return(invisible(dend))
  }
)
