#' @title Functional Microbiome Representation based on
#' \code{TreeSummarizedExperiment}
#'
#' @param data a DataFrame-like object that includes columns specfiying
#' the species, metabolites and fluxes in the microbiome. The fluxes can
#' either be weighted or unweighted (all of magnitude 1).
#' @param name a character scalar specifying the name of the Microbiome
#' @param species_col Character scalar specfiying the name of the species
#' column, defaults to 'spec'.
#' @param metabolite_col Character scalar specifying the name of the metabolite
#' column, defaults to 'met'.
#' @param flux_col Character scalar specifying the name of the flux column,
#' defaults to 'flux'.
#' @param ... Additional arguments to be passed to the constructor.
#'
#' @export
#'
#' @importFrom SummarizedExperiment metadata<-
#' @importFrom TreeSummarizedExperiment TreeSummarizedExperiment
#' @importFrom tibble tibble
#' @importFrom dplyr filter
ConsortiumMetabolism <- function(
  data,
  name = NA_character_,
  species_col = "species",
  metabolite_col = "met",
  flux_col = "flux",
  ...
) {
  ### Include checkfor all columns and use str_detect
  ### Build the renaming of the columns

  data <- data |>
    dplyr::rename(
      species = {{ species_col }},
      met = {{ metabolite_col }},
      flux = {{ flux_col }}
    )

  stopifnot(exprs = {
    all(c("species", "met", "flux") %in% names(data))
  })

  # Determine whether or not the fluxes are weighted or not
  weighted <- !all(data$flux**2 == 1)

  # if (!weighted) assays <- list(Binary = bin_mat, nEdges = n_edges)

  mets <- tibble::tibble(
    met = sort(unique(data$met)),
    index = seq_len(length(unique(data$met)))
  )

  tb <- data |>
    # Ensure that no zero values are included
    dplyr::filter(.data$flux != 0) # |>
  # recode the species names
  # dplyr::left_join(species, by = "species") |>
  # dplyr::select(-"species", species = "species_code")

  ### Check if any rows have flux == 0
  cons <- tb |>
    dplyr::filter(.data$flux < 0) |>
    dplyr::mutate(flux = .data$flux * -1) |>
    dplyr::reframe(flux = sum(.data$flux), .by = c("species", "met")) |>
    dplyr::rename(consumed = "met")

  prod <- tb |>
    dplyr::filter(.data$flux > 0) |>
    dplyr::reframe(flux = sum(.data$flux), .by = c("species", "met")) |>
    dplyr::rename(produced = "met")

  out <- cons |>
    ### Need to double check the logic with alberto
    dplyr::inner_join(
      prod,
      by = "species",
      suffix = c("_cons", "_prod"),
      relationship = "many-to-many"
    ) |>
    tidyr::nest(data = c("species", "flux_cons", "flux_prod")) |>
    dplyr::mutate(
      n_species = purrr::map_dbl(.data$data, \(x) nrow(x)),

      # Get the sum of the fluxes
      c_sum = purrr::map_dbl(.data$data, \(x) sum(x$flux_cons)),
      p_sum = purrr::map_dbl(.data$data, \(x) sum(x$flux_prod)),

      c_prob = purrr::map(.data$data, \(x) x$flux_cons / sum(x$flux_cons)),
      p_prob = purrr::map(.data$data, \(x) x$flux_prod / sum(x$flux_prod)),

      c_eff = purrr::map_dbl(
        .data$c_prob,
        \(x) round(2**(-sum(x * log2(x))), 2)
      ),
      p_eff = purrr::map_dbl(
        .data$p_prob,
        \(x) round(2**(-sum(x * log2(x))), 2)
      )
    ) |>
    # Replace with the indeces for the metabolites
    dplyr::left_join(mets, by = c(consumed = "met")) |>
    dplyr::rename(c_ind = "index") |>
    dplyr::left_join(mets, by = c(produced = "met")) |>
    dplyr::rename(p_ind = "index")

  assays <- list(
    Binary = sparseMatrix(out$c_ind, out$p_ind, x = 1),
    nEdges = sparseMatrix(out$c_ind, out$p_ind, x = out$n_species),
    Consumption = sparseMatrix(out$c_ind, out$p_ind, x = out$c_sum),
    Production = sparseMatrix(out$c_ind, out$p_ind, x = out$p_sum),
    EffectiveConsumption = sparseMatrix(out$c_ind, out$p_ind, x = out$c_eff),
    EffectiveProduction = sparseMatrix(out$c_ind, out$p_ind, x = out$p_eff)
  )
  ### How to properly deal with col data?
  tse <- TreeSummarizedExperiment(
    assays = assays,
    rowData = mets,
    colData = mets
  )

  graphs <- list(
    igraph::graph_from_adjacency_matrix(
      adjmatrix = assays$Binary,
      mode = "directed"
    )
  )

  newConsortiumMetabolism(
    tse,
    Name = name,
    Edges = list(),
    Weighted = weighted,
    InputData = data,
    Metabolites = unique(data$met),
    Graphs = graphs
  )
}
