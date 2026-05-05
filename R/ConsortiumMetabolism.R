#' @title Functional Microbiome Representation
#'
#' @description Creates a \code{ConsortiumMetabolism} object
#'   representing metabolic interactions in a microbial community.
#'   The object stores metabolite consumption and production by
#'   different species, along with flux sums and effective fluxes.
#'
#' @slot Name character. Display name for the consortium.
#' @slot Description character. Optional short description.
#' @slot Pathways data.frame. Pathway list of metabolic
#'   interactions with per-pathway metrics (species, flux
#'   sums, effective diversity).
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
#' @param growth Optional named numeric vector of per-species
#'   growth rates (e.g. FBA objective values). Names must match
#'   species in \code{data}. Stored in
#'   \code{metadata(cm)$growth} and retrievable via the
#'   \code{\link{growth}} accessor.
#' @param species_col Character scalar for the species column
#'   name, defaults to \code{"species"}.
#' @param metabolite_col Character scalar for the metabolite
#'   column name, defaults to \code{"metabolite"}.
#' @param flux_col Character scalar for the flux column name,
#'   defaults to \code{"flux"}.
#' @param verbose Logical scalar. If \code{TRUE}, prints progress
#'   messages during construction. Defaults to \code{FALSE}.
#'
#' @return A \code{ConsortiumMetabolism} object.
#'
#' @seealso \link[TreeSummarizedExperiment]{TreeSummarizedExperiment-class}
#'
#' @examples
#' cm <- synCM("example", n_species = 3, max_met = 5)
#' cm
#'
#' @export
#'
#' @importFrom S4Vectors metadata
#' @importFrom SummarizedExperiment metadata<-
#' @importFrom TreeSummarizedExperiment TreeSummarizedExperiment
ConsortiumMetabolism <- function(
    data,
    name = NA_character_,
    growth = NULL,
    species_col = "species",
    metabolite_col = "metabolite",
    flux_col = "flux",
    verbose = FALSE
) {
    # Validate and prepare input data
    data <- .prepareInputData(data, species_col, metabolite_col, flux_col)

    # Validate growth parameter
    if (!is.null(growth)) {
        growth <- .validateGrowth(growth, data$species)
    }

    # Create metabolite index mapping
    mets <- .createMetaboliteIndex(data)

    # Process flux data into consumption and production
    tb <- .filterNonzeroFlux(data)
    if (nrow(tb) == 0L) {
        cli::cli_abort(
            c(
                "All flux values in {.arg data} are zero.",
                "i" = "Provide at least one row with a non-zero flux \\
                (negative for consumption, positive for production)."
            )
        )
    }
    cons <- .calculateConsumption(tb)
    prod <- .calculateProduction(tb)

    # Get only producing or consuming species
    only_cons <- setdiff(unique(cons$species), unique(prod$species))
    only_prod <- setdiff(unique(prod$species), unique(cons$species))
    if (length(only_cons) == 0 && length(only_prod) == 0) {
        mets <- dplyr::filter(mets, .data$met != "media")
    } else {
        if (length(only_cons) > 0) {
            if (verbose) {
                cli::cli_inform(
                    "Species {.code {only_cons}} in {.val {name}} \\
                have no production flux; \\
                synthetic {.val media} production added."
                )
            }
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
            if (verbose) {
                cli::cli_inform(
                    "Species {.code {only_prod}} in {.val {name}} \\
                have no consumption flux; \\
                synthetic {.val media} consumption added."
                )
            }
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

    # Create pathway data with metrics
    out <- .createPathwayData(cons, prod, mets)

    # Generate assay matrices
    assays <- .createAssayMatrices(out, mets)

    n_selfloops <- sum(Matrix::diag(assays$Binary) != 0)
    if (n_selfloops > 0L) {
        cli::cli_warn(
            "{n_selfloops} self-loop pathway{?s} \\
            found in {.val {name}} \\
            (metabolite consumed and produced by the \\
            same pathway). Check input data."
        )
    }

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
    cm <- newConsortiumMetabolism(
        tse,
        Name = name,
        Pathways = out,
        Weighted = !all(abs(data$flux) == 1),
        InputData = as.data.frame(data),
        Metabolites = unique(data$met),
        Graphs = graphs
    )

    # Store growth rates in metadata
    if (!is.null(growth)) {
        metadata(cm)$growth <- growth
    }

    cm
}

#' Prepare and validate input data
#' @noRd
.prepareInputData <- function(data, species_col, metabolite_col, flux_col) {
    if (!is.data.frame(data)) {
        cli::cli_abort(
            "{.arg data} must be a data.frame, \\
            not {.cls {class(data)}}."
        )
    }
    required <- c(species_col, metabolite_col, flux_col)
    missing_cols <- setdiff(required, names(data))
    if (length(missing_cols) > 0) {
        cli::cli_abort(
            "Column{?s} {.val {missing_cols}} not found
            in {.arg data}."
        )
    }

    if (nrow(data) == 0L) {
        cli::cli_abort(
            c(
                "{.arg data} has 0 rows.",
                "i" = "Provide at least one row of \\
                species/metabolite/flux data."
            )
        )
    }

    bad_type <- character()
    if (!is.character(data[[species_col]])) {
        bad_type <- c(bad_type, species_col)
    }
    if (!is.character(data[[metabolite_col]])) {
        bad_type <- c(bad_type, metabolite_col)
    }
    if (length(bad_type) > 0L) {
        bad_first <- bad_type[1] # nolint: object_usage_linter.
        cli::cli_abort(
            c(
                "Column{?s} {.val {bad_type}} in {.arg data} \\
                must be {.cls character}.",
                "i" = "Cast with {.code \\
                data[[\"{bad_first}\"]] <- \\
                as.character(data[[\"{bad_first}\"]])}."
            )
        )
    }

    na_rows <- which(
        is.na(data[[species_col]]) |
            is.na(data[[metabolite_col]]) |
            is.na(data[[flux_col]])
    )
    if (length(na_rows) > 0L) {
        preview <- if (length(na_rows) > 10L) { # nolint: object_usage_linter.
            c(as.character(na_rows[seq_len(10L)]), "...")
        } else {
            as.character(na_rows)
        }
        # nolint next: object_usage_linter.
        cols <- c(species_col, metabolite_col, flux_col)
        cli::cli_abort(
            c(
                "{.arg data} contains NA values in the \\
                species/metabolite/flux columns \\
                ({.val {cols}}) at row(s) {.val {preview}}.",
                "i" = "Drop incomplete rows with \\
                {.code tidyr::drop_na(data, \\
                {species_col}, {metabolite_col}, {flux_col})}."
            )
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

#' Create pathway data with all metrics
#' @noRd
.createPathwayData <- function(cons, prod, mets) {
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
                .data$data,
                \(x) nrow(x),
                numeric(1)
            ),
            c_sum = vapply(
                .data$data,
                \(x) sum(x$flux_cons),
                numeric(1)
            ),
            p_sum = vapply(
                .data$data,
                \(x) sum(x$flux_prod),
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
            c_eff = vapply(
                .data$c_prob,
                \(x) {
                    round(2**(-sum(x * log2(x))), 2)
                },
                numeric(1)
            ),
            p_eff = vapply(
                .data$p_prob,
                \(x) {
                    round(
                        2**(-sum(x * log2(x))),
                        2
                    )
                },
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
        nSpecies = sparseMatrix(
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

#' Validate growth rate vector
#' @noRd
.validateGrowth <- function(growth, species_in_data) {
    if (!is.numeric(growth)) {
        cli::cli_abort(
            "{.arg growth} must be a numeric vector, \\
            not {.cls {class(growth)}}."
        )
    }
    if (is.null(names(growth))) {
        cli::cli_abort(
            "{.arg growth} must be a named numeric \\
            vector (names = species)."
        )
    }
    if (anyDuplicated(names(growth))) {
        # nolint next: object_usage_linter.
        dupes <- names(growth)[duplicated(names(growth))]
        cli::cli_abort(
            "{.arg growth} has duplicate name{?s}: \\
            {.val {unique(dupes)}}."
        )
    }
    if (any(growth < 0, na.rm = TRUE)) {
        cli::cli_abort(
            "{.arg growth} values must be non-negative."
        )
    }
    known_species <- unique(species_in_data)
    unknown <- setdiff(names(growth), known_species)
    if (length(unknown) > 0L) {
        cli::cli_abort(
            "{.arg growth} name{?s} not found in \\
            {.arg data}: {.val {unknown}}."
        )
    }
    growth
}
