#' @title Functional Microbiome Representation
#'
#' @description Creates a \code{ConsortiumMetabolism} object
#'   representing metabolic interactions in a microbial community.
#'   The object stores metabolite consumption and production by
#'   different species, along with flux sums and effective fluxes.
#'
#' @slot Name character. Display name for the consortium.
#' @slot Description character. Optional short description.
#' @slot Edges data.frame. Edge list of metabolic interactions
#'   with per-edge metrics (species, flux sums, effective
#'   diversity).
#' @slot Weighted logical. Whether flux magnitudes are used.
#' @slot InputData data.frame. Original input data (species,
#'   metabolite, flux columns).
#' @slot Metabolites character. Unique metabolite identifiers.
#' @slot Graphs list. List containing an igraph object of the
#'   metabolic network.
#'
#' @param data a data.frame with columns for species, metabolites
#'   and fluxes. Fluxes can be weighted or unweighted (magnitude
#'   1).
#' @param name Character scalar giving the consortium name.
#' @param species_col Character scalar for the species column
#'   name, defaults to \code{"species"}.
#' @param metabolite_col Character scalar for the metabolite
#'   column name, defaults to \code{"met"}.
#' @param flux_col Character scalar for the flux column name,
#'   defaults to \code{"flux"}.
#' @param ... Additional arguments passed to the constructor.
#'
#' @return A \code{ConsortiumMetabolism} object.
#'
#' @seealso \link[TreeSummarizedExperiment:TreeSummarizedExperiment-class]{TreeSummarizedExperiment}
#'
#' @examples
#' cm <- synCM("example", n_species = 3, max_met = 5)
#' cm
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
    data <- .prepareInputData(data, species_col, metabolite_col, flux_col)
    # Create metabolite index mapping
    mets <- .createMetaboliteIndex(data)

    # Process flux data into consumption and production
    tb <- .filterNonzeroFlux(data)
    cons <- .calculateConsumption(tb)
    prod <- .calculateProduction(tb)

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
                    tibble::tibble(
                        species = only_cons,
                        produced = "media",
                        flux = 1
                    )
                )
        }

        if (length(only_prod) > 0) {
            cli::cli_alert_info(
                "{.val {name}} {.code {only_prod}} only produce,
        consumption set to 'media'."
            )
            cons <- cons |>
                dplyr::bind_rows(
                    tibble::tibble(
                        species = only_prod,
                        consumed = "media",
                        flux = 1
                    )
                )
        }
    }

    # Create edge data with metrics
    out <- .createEdgeData(cons, prod, mets)

    # Generate assay matrices
    assays <- .createAssayMatrices(out, mets)

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
.prepareInputData <- function(data, species_col, metabolite_col, flux_col) {
    required <- c(species_col, metabolite_col, flux_col)
    missing_cols <- setdiff(required, names(data))
    if (length(missing_cols) > 0) {
        cli::cli_abort(
            "Column{?s} {.val {missing_cols}} not found
            in {.arg data}."
        )
    }

    data |>
        rename(
            species = {{ species_col }},
            met = {{ metabolite_col }},
            flux = {{ flux_col }}
        )
}

#' Create metabolite index mapping
#' @noRd
.createMetaboliteIndex <- function(data) {
    # tibble(met = sort(unique(data$met))) |>
    #   tibble::rowid_to_column(var = "index")

    # Testing whether this works
    tibble(met = c(sort(unique(data$met)), "media")) |>
        tibble::rowid_to_column(var = "index")
}

#' Filter out zero flux values
#' @noRd
.filterNonzeroFlux <- function(data) {
    dplyr::filter(data, .data$flux != 0)
}

#' Calculate consumption metrics
#' @noRd
.calculateConsumption <- function(tb) {
    tb |>
        dplyr::filter(.data$flux < 0) |>
        mutate(flux = .data$flux * -1) |>
        # Sum bc in eg cooc data each edge is given with its own flux so then
        # there can be multiple prod and cons values per met
        reframe(flux = sum(.data$flux), .by = c("species", "met")) |>
        rename(consumed = "met")
}

#' Calculate production metrics
#' @noRd
.calculateProduction <- function(tb) {
    tb |>
        dplyr::filter(.data$flux > 0) |>
        # Sum bc in eg cooc data each edge is given with its own flux so then
        # there can be multiple prod and cons values per met
        reframe(flux = sum(.data$flux), .by = c("species", "met")) |>
        rename(produced = "met")
}

#' Create edge data with all metrics
#' @noRd
.createEdgeData <- function(cons, prod, mets) {
    cons |>
        dplyr::inner_join(
            prod,
            by = "species",
            suffix = c("_cons", "_prod"),
            relationship = "many-to-many"
        ) |>
        nest(data = c("species", "flux_cons", "flux_prod")) |>
        mutate(
            n_species = vapply(
                .data$data, \(x) nrow(x), numeric(1)
            ),
            c_sum = vapply(
                .data$data, \(x) sum(x$flux_cons),
                numeric(1)
            ),
            p_sum = vapply(
                .data$data, \(x) sum(x$flux_prod),
                numeric(1)
            ),
            c_prob = lapply(
                .data$data,
                \(x) x$flux_cons / sum(x$flux_cons)
            ),
            p_prob = lapply(
                .data$data,
                \(x) x$flux_prod / sum(x$flux_prod)
            ),
            c_eff = vapply(.data$c_prob, \(x) {
                round(2**(-sum(x * log2(x))), 2)
            }, numeric(1)),
            p_eff = vapply(
                .data$p_prob,
                \(x) round(
                    2**(-sum(x * log2(x))), 2
                ),
                numeric(1)
            )
        ) |>
        left_join(mets, by = c(consumed = "met")) |>
        rename(c_ind = "index") |>
        left_join(mets, by = c(produced = "met")) |>
        rename(p_ind = "index")
}

#' Create assay matrices
#' @noRd
.createAssayMatrices <- function(out, mets) {
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
