#' Import SMETANA Detailed Output
#'
#' Imports raw SMETANA `--detailed` output (one TSV file per consortium) into
#' a ConsortiumMetabolism (single file/data.frame) or ConsortiumMetabolismSet
#' (directory of files). The raw format has columns `community`, `medium`,
#' `receiver`, `donor`, `compound`, `scs`, `mus`, `mps`, `smetana`, where each
#' row describes a cross-feeding interaction: `donor` produces `compound`,
#' `receiver` consumes it. This function collapses those directed interactions
#' into ramen's (species, metabolite, flux) edge-list representation.
#'
#' @param data Either a file path to a single SMETANA TSV, a directory path
#'   containing multiple SMETANA TSVs, or a pre-loaded data.frame with the
#'   required columns (`receiver`, `donor`, `compound`, and `smetana` if
#'   `use_scores = TRUE`).
#' @param name Consortium (or set) name as a length-1 character scalar. If
#'   `NULL` (default), names are derived from the input filename(s) by stripping
#'   the `.tsv_detailed.tsv` suffix. Required when `data` is a data.frame. If
#'   `data` is a directory, `name` is passed through to
#'   [ConsortiumMetabolismSet()].
#' @param use_scores If `FALSE`, flux is binary (+1 for production,
#'   -1 for consumption). If `TRUE` (default), flux magnitude is the
#'   SMETANA score,
#'   aggregated by `max()` when the same (species, compound) pair appears
#'   across multiple interaction partners.
#' @param normalize_ids If `TRUE` (default), normalize compound IDs via
#'   [.normalizeBiggIds()] (strips `M_` prefix and `_e` compartment suffix).
#' @param verbose If `TRUE` (default), emit progress messages via
#'   cli progress messages.
#'
#' @return A `ConsortiumMetabolism` object if `data` is a single file or
#'   data.frame; a `ConsortiumMetabolismSet` if `data` is a directory.
#'
#' @seealso [ConsortiumMetabolism()], [ConsortiumMetabolismSet()],
#'   [importMisosoup()]
#'
#' @export
#'
#' @examples
#' if (FALSE) {
#' # Single file -> ConsortiumMetabolism
#' cm <- importSmetana("path/to/bq_0.tsv_detailed.tsv")
#'
#' # Directory of files -> ConsortiumMetabolismSet
#' cms <- importSmetana("path/to/bq_subsample/")
#'
#' # Weighted fluxes
#' cm <- importSmetana("path/to/bq_0.tsv_detailed.tsv", use_scores = TRUE)
#' }
importSmetana <- function(
    data,
    name = NULL,
    use_scores = TRUE,
    normalize_ids = TRUE,
    verbose = TRUE
) {
    if (is.character(data) && dir.exists(data)) {
        file_paths <- list.files(
            data,
            pattern = "\\.tsv_detailed\\.tsv$",
            full.names = TRUE
        )
        if (length(file_paths) == 0) {
            cli::cli_abort("No {.val .tsv} files found in {.path {data}}.")
        }
        indices <- if (verbose) {
            cli::cli_progress_along(
                file_paths,
                name = "Importing SMETANA"
            )
        } else {
            seq_along(file_paths)
        }
        cm_list <- lapply(indices, function(i) {
            importSmetana(
                file_paths[[i]],
                use_scores = use_scores,
                normalize_ids = normalize_ids,
                verbose = FALSE
            )
        })
        cms_name <- if (is.null(name)) basename(normalizePath(data)) else name
        return(ConsortiumMetabolismSet(cm_list, name = cms_name))
    } else if (is.character(data) && file.exists(data)) {
        if (is.null(name)) {
            name <- sub("(\\.tsv_detailed)?\\.tsv$", "", basename(data))
        }
        data <- utils::read.delim(data)
    } else if (is.data.frame(data)) {
        if (is.null(name)) {
            cli::cli_abort(
                "{.arg name} is required when {.arg data} is a\\
                {.cls data.frame}."
            )
        }
    } else {
        cli::cli_abort(c(
            "{.arg data} must be a file path, directory, or data.frame.",
            "i" = "Got {.cls {class(data)}}."
        ))
    }
    required <- c("receiver", "donor", "compound")
    if (use_scores) {
        required <- c(required, "smetana")
    }
    missing_cols <- setdiff(required, names(data))
    if (length(missing_cols) > 0) {
        cli::cli_abort(
            "Column{?s} {.val {missing_cols}} not found in {.arg data}."
        )
    }
    for (col in required) {
        if (anyNA(data[[col]])) {
            cli::cli_abort("{.field {col}} contains missing values.")
        }
    }
    if (use_scores && any(data$smetana <= 0)) {
        cli::cli_abort("{.field smetana} contains zero or negative values")
    }
    if (!use_scores) {
        production <- data |>
            dplyr::distinct(.data$donor, .data$compound) |>
            dplyr::rename(
                species = "donor",
                metabolite = "compound"
            ) |>
            dplyr::mutate(flux = 1L)
        consumption <- data |>
            dplyr::distinct(.data$receiver, .data$compound) |>
            dplyr::rename(
                species = "receiver",
                metabolite = "compound",
            ) |>
            dplyr::mutate(flux = -1L)
    } else {
        production <- data |>
            dplyr::group_by(.data$donor, .data$compound) |>
            dplyr::summarise(flux = sum(.data$smetana), .groups = "drop") |>
            dplyr::rename(species = "donor", metabolite = "compound")
        consumption <- data |>
            dplyr::group_by(.data$receiver, .data$compound) |>
            dplyr::summarise(flux = -sum(.data$smetana), .groups = "drop") |>
            dplyr::rename(species = "receiver", metabolite = "compound")
    }
    edges <- dplyr::bind_rows(production, consumption)
    if (normalize_ids) {
        edges <- dplyr::mutate(
            edges,
            metabolite = .normalizeBiggIds(.data$metabolite)
        )
    }
    if (verbose) {
        cli::cli_inform(c("i" = "Importing {.val {name}}: {nrow(edges)} edges"))
    }
    ConsortiumMetabolism(edges, name = name)
}
