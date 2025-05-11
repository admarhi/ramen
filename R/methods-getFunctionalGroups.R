### Make this into a method and make applicable for CMA
#' @export
getFunctionalGroups <- function(object) {
  tb <- object@Edges

  # To calculate the functional groups
  # Get a tibble with all rxns for each species
  rxns_per_species <- tb |>
    dplyr::select(1:3) |>
    unique() |>
    dplyr::mutate(edge = paste0(.data$consumed, "-", .data$produced)) |>
    dplyr::select("species", "edge")

  unique_species <-
    tibble::tibble(species = unique(rxns_per_species$species)) |>
    tibble::rowid_to_column("ind")

  # Get only unique combinations to save on computation time
  species_combinations <- tidyr::expand_grid(
    x = unique_species$ind,
    y = unique_species$ind
  ) |>
    dplyr::filter(x <= y) |>
    dplyr::left_join(unique_species, by = c("x" = "ind")) |>
    dplyr::left_join(unique_species, by = c("y" = "ind")) |>
    janitor::clean_names() |>
    dplyr::mutate(similarity = 0)

  for (i in seq_len(nrow(species_combinations))) {
    species_x <- species_combinations$species_x[i]
    species_y <- species_combinations$species_y[i]
    rxns_set_x <- rxns_per_species$edge[rxns_per_species$species == species_x]
    rxns_set_y <- rxns_per_species$edge[rxns_per_species$species == species_y]

    intersection <- length(intersect(rxns_set_x, rxns_set_y))
    union <- length(union(rxns_set_x, rxns_set_y))

    species_combinations$similarity[i] <- intersection / union
  }

  species_combinations

  Matrix::sparseMatrix(
    species_combinations$x,
    species_combinations$y,
    x = species_combinations$similarity,
    dimnames = list(unique_species$species, unique_species$species)
  ) |>
    as.matrix() |>
    ### Give these arguments as options?
    dist() |>
    hclust(method = "complete") |>
    as.dendrogram() |>
    dendextend::color_branches(k = 4) |>
    plot()
}
