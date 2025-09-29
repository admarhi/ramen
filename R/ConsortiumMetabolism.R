#' @title Functional Microbiome Representation based on TreeSummarizedExperiment
#'
#' @description Creates a ConsortiumMetabolism object representing metabolic
#' interactions in a microbial community. The object contains information about
#' metabolite consumption
#' and production by different species, along with various metrics like flux
#' sums and
#' effective fluxes.
#'
#' @param data a DataFrame-like object that includes columns specifying
#' the species, metabolites and fluxes in the microbiome. The fluxes can
#' either be weighted or unweighted (all of magnitude 1).
#' @param name a character scalar specifying the name of the Microbiome
#' @param species_col Character scalar specifying the name of the species
#' column, defaults to 'species'.
#' @param metabolite_col Character scalar specifying the name of the metabolite
#' column, defaults to 'met'.
#' @param flux_col Character scalar specifying the name of the flux column,
#' defaults to 'flux'.
#' @param ... Additional arguments to be passed to the constructor.
#'
#' @return A ConsortiumMetabolism object containing:
#' \itemize{
#'   \item Assays for binary interactions, edge counts, consumption/production
#' metrics
#'   \item Row and column data about metabolites
#'   \item Graph representation of the metabolic network
#'   \item Original input data and computed edge information
#' }
#'
#' @export
#'
#' @importFrom SummarizedExperiment metadata<-
#' @importFrom TreeSummarizedExperiment TreeSummarizedExperiment
ConsortiumMetabolism <- function(
  data,
  name = NA_character_,
  species_col = "species",
  metabolite_col = "met",
  flux_col = "flux",
  ...
) {
  # Validate and prepare input data
  data <- prepare_input_data(data, species_col, metabolite_col, flux_col)
  # Create metabolite index mapping
  mets <- create_metabolite_index(data)

  # Process flux data into consumption and production
  tb <- filter_nonzero_flux(data)
  cons <- calculate_consumption(tb)
  prod <- calculate_production(tb)

  # Get only producing or consuming species
  only_cons <- setdiff(unique(cons$species), unique(prod$species))
  only_prod <- setdiff(unique(prod$species), unique(cons$species))
  if (length(only_cons) == 0 && length(only_prod) == 0) {
    mets <- dplyr::filter(mets, .data$met != "media")
  } else {
    if (length(only_cons) > 0) {
      cli::cli_alert_info(
        "{.val {name}} {.code {only_cons}} only consume,
        production set to 'media'."
      )
      prod <- prod |>
        dplyr::bind_rows(
          tibble::tibble(species = only_cons, produced = "media", flux = 1)
        )
    }

    if (length(only_prod) > 0) {
      cli::cli_alert_info(
        "{.val {name}} {.code {only_prod}} only produce,
        consumption set to 'media'."
      )
      # message(only_prod, " only produce, consumption set to 'media'")
      cons <- cons |>
        dplyr::bind_rows(
          tibble::tibble(species = only_prod, consumed = "media", flux = 1)
        )
    }
  }

  # Create edge data with metrics
  out <- create_edge_data(cons, prod, mets)

  # Generate assay matrices
  assays <- create_assay_matrices(out, mets)

  # Create TreeSummarizedExperiment object
  tse <- TreeSummarizedExperiment(
    assays = assays,
    rowData = mets,
    colData = mets
  )

  # Create graph representation
  graphs <- list(
    igraph::graph_from_adjacency_matrix(
      adjmatrix = assays$Binary,
      mode = "directed"
    )
  )

  # Return final ConsortiumMetabolism object
  newConsortiumMetabolism(
    tse,
    Name = name,
    Edges = out,
    Weighted = !all(data$flux**2 == 1),
    InputData = as.data.frame(data),
    Metabolites = unique(data$met),
    Graphs = graphs
  )
}

#' Prepare and validate input data
#' @noRd
prepare_input_data <- function(data, species_col, metabolite_col, flux_col) {
  stopifnot(exprs = {
    all(c(species_col, metabolite_col, flux_col) %in% names(data))
  })

  data |>
    rename(
      species = {{ species_col }},
      met = {{ metabolite_col }},
      flux = {{ flux_col }}
    )
}

#' Create metabolite index mapping
#' @noRd
create_metabolite_index <- function(data) {
  # tibble(met = sort(unique(data$met))) |>
  #   tibble::rowid_to_column(var = "index")

  # Testing whether this works
  tibble(met = c(sort(unique(data$met)), "media")) |>
    tibble::rowid_to_column(var = "index")
}

#' Filter out zero flux values
#' @noRd
filter_nonzero_flux <- function(data) {
  filter(data, .data$flux != 0)
}

#' Calculate consumption metrics
#' @noRd
calculate_consumption <- function(tb) {
  tb |>
    filter(.data$flux < 0) |>
    mutate(flux = .data$flux * -1) |>
    # Sum bc in eg cooc data each edge is given with its own flux so then
    # there can be multiple prod and cons values per met
    reframe(flux = sum(.data$flux), .by = c("species", "met")) |>
    rename(consumed = "met")
}

#' Calculate production metrics
#' @noRd
calculate_production <- function(tb) {
  tb |>
    filter(.data$flux > 0) |>
    # Sum bc in eg cooc data each edge is given with its own flux so then
    # there can be multiple prod and cons values per met
    reframe(flux = sum(.data$flux), .by = c("species", "met")) |>
    rename(produced = "met")
}

#' Create edge data with all metrics
#' @noRd
create_edge_data <- function(cons, prod, mets) {
  cons |>
    dplyr::inner_join(
      prod,
      by = "species",
      suffix = c("_cons", "_prod"),
      relationship = "many-to-many"
    ) |>
    nest(data = c("species", "flux_cons", "flux_prod")) |>
    mutate(
      n_species = map_dbl(.data$data, \(x) nrow(x)),
      c_sum = map_dbl(.data$data, \(x) sum(x$flux_cons)),
      p_sum = map_dbl(.data$data, \(x) sum(x$flux_prod)),
      c_prob = map(.data$data, \(x) x$flux_cons / sum(x$flux_cons)),
      p_prob = map(.data$data, \(x) x$flux_prod / sum(x$flux_prod)),
      c_eff = map_dbl(.data$c_prob, \(x) round(2**(-sum(x * log2(x))), 2)),
      p_eff = map_dbl(.data$p_prob, \(x) round(2**(-sum(x * log2(x))), 2))
    ) |>
    left_join(mets, by = c(consumed = "met")) |>
    rename(c_ind = "index") |>
    left_join(mets, by = c(produced = "met")) |>
    rename(p_ind = "index")
}

#' Create assay matrices
#' @noRd
create_assay_matrices <- function(out, mets) {
  # Create the dimnames
  dimnames <- list(mets$met, mets$met)
  n <- nrow(mets)
  list(
    Binary = sparseMatrix(
      out$c_ind,
      out$p_ind,
      x = 1,
      dims = c(n, n),
      dimnames = dimnames
    ),
    nEdges = sparseMatrix(
      out$c_ind,
      out$p_ind,
      x = out$n_species,
      dims = c(n, n),
      dimnames = dimnames
    ),
    Consumption = sparseMatrix(
      out$c_ind,
      out$p_ind,
      x = out$c_sum,
      dims = c(n, n),
      dimnames = dimnames
    ),
    Production = sparseMatrix(
      out$c_ind,
      out$p_ind,
      x = out$p_sum,
      dims = c(n, n),
      dimnames = dimnames
    ),
    EffectiveConsumption = sparseMatrix(
      out$c_ind,
      out$p_ind,
      x = out$c_eff,
      dims = c(n, n),
      dimnames = dimnames
    ),
    EffectiveProduction = sparseMatrix(
      out$c_ind,
      out$p_ind,
      x = out$p_eff,
      dims = c(n, n),
      dimnames = dimnames
    )
  )
}
