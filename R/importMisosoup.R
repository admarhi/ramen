#' Import MiSoSoup YAML Output
#'
#' Imports raw MiSoSoup YAML output into a
#' [ConsortiumMetabolismSet-class] object. MiSoSoup YAMLs describe many
#' consortia per file (one per `substrate` / second-level-key /
#' `solution` triple), so a single call always produces a CMS regardless
#' of whether the input is a single file, a directory of files, or a
#' pre-loaded nested list.
#'
#' The function auto-detects format differences between MiSoSoup runs:
#' the growth-row key (`Growth_<species>` in older runs,
#' `R_Biomass_<species>` in newer runs), COBRApy-style `__NN__` escape
#' sequences in reaction IDs (e.g. `__40__` for `(`), and the presence
#' or absence of a focal-strain layer below each substrate.
#'
#' Media-level exchange fluxes (reactions without a species suffix) are
#' stashed in `metadata(cm)$media` for each consortium, preserving the
#' original bounds information if downstream code needs it.
#'
#' @param data Either a file path to a single MiSoSoup YAML, a
#'   directory path containing multiple YAML files, or a pre-loaded
#'   nested list from `yaml::read_yaml()`.
#' @param name Name for the returned `ConsortiumMetabolismSet`. If
#'   `NULL` (default), the name is derived from the input filename (for
#'   a file) or directory basename (for a directory). **Required** when
#'   `data` is a pre-loaded list.
#' @param normalize_ids If `TRUE` (default), normalize metabolite IDs
#'   via `.normalizeBiggIds()` (strip `R_EX_` prefix, `_e` suffix) and
#'   decode COBRApy `__NN__` escape sequences via `.decodeBiggEscapes()`.
#'   If `FALSE`, keep metabolite IDs in their raw form.
#' @param verbose If `TRUE` (default), emit progress messages via
#'   [cli::cli_inform()] and show a progress bar when importing a
#'   directory.
#'
#' @return A [ConsortiumMetabolismSet-class] object containing one
#'   [ConsortiumMetabolism-class] per viable consortium found in the
#'   input. Zero-growth solutions are silently skipped.
#'
#' @seealso [overviewMisosoup()] for a pre-import summary of a MiSoSoup
#'   data structure; [importSmetana()] for the sibling SMETANA import.
#'
#' @export
#'
#' @examples
#' \donttest{
#' # Single YAML file -> CMS
#' # cms <- importMisosoup("path/to/misosoup.yaml")
#'
#' # Directory of YAML files -> CMS with consortia from all files
#' # cms <- importMisosoup("path/to/misosoup_dir/", name = "experiment1")
#'
#' # Pre-loaded list -> CMS (name is required here)
#' # raw <- yaml::read_yaml("path/to/misosoup.yaml")
#' # cms <- importMisosoup(raw, name = "experiment1")
#' }
importMisosoup <- function(
    data,
    name = NULL,
    normalize_ids = TRUE,
    verbose = TRUE
) {
    if (is.character(data) && length(data) == 1 && dir.exists(data)) {
        files <- list.files(
            data,
            pattern = "\\.(yaml|yml)$",
            full.names = TRUE
        )
        if (length(files) == 0) {
            cli::cli_abort(
                "No {.val .yaml} files found in {.path {data}}."
            )
        }
        indices <- if (verbose) {
            cli::cli_progress_along(
                files,
                name = "Importing MiSoSoup"
            )
        } else {
            seq_along(files)
        }
        cm_list <- list()
        for (i in indices) {
            raw <- yaml::read_yaml(files[[i]])
            file_cms <- .buildCMsFromMisosoup(
                raw,
                normalize_ids = normalize_ids
            )
            cm_list <- c(cm_list, file_cms)
        }
        cms_name <- if (is.null(name)) {
            basename(normalizePath(data))
        } else {
            name
        }
        return(ConsortiumMetabolismSet(cm_list, name = cms_name))
    } else if (is.character(data) && length(data) == 1 && file.exists(data)) {
        cms_name <- if (is.null(name)) {
            sub("\\.(yaml|yml)$", "", basename(data))
        } else {
            name
        }
        raw <- yaml::read_yaml(data)
        cm_list <- .buildCMsFromMisosoup(
            raw,
            normalize_ids = normalize_ids
        )
        if (verbose) {
            cli::cli_inform(c(
                "i" = "Imported {length(cm_list)} consort{?ium/ia} \\
                       from {.path {basename(data)}}."
            ))
        }
        return(ConsortiumMetabolismSet(cm_list, name = cms_name))
    } else if (is.list(data) && !is.data.frame(data)) {
        if (is.null(name)) {
            cli::cli_abort(
                "{.arg name} is required when {.arg data} is a \\
                 pre-loaded list."
            )
        }
        cm_list <- .buildCMsFromMisosoup(
            data,
            normalize_ids = normalize_ids
        )
        if (verbose) {
            cli::cli_inform(c(
                "i" = "Imported {length(cm_list)} consort{?ium/ia}."
            ))
        }
        return(ConsortiumMetabolismSet(cm_list, name = name))
    }
    cli::cli_abort(c(
        "{.arg data} must be a file path, directory, or pre-loaded list.",
        "i" = "Got {.cls {class(data)}}."
    ))
}


#' Overview of MiSoSoup Data Structure
#'
#' Provides a summary of a MiSoSoup nested list, including the number of
#' solutions and zero-growth solutions for each substrate /
#' second-level-key combination. Useful for inspecting a data structure
#' before importing with [importMisosoup()].
#'
#' @param data A nested list from `yaml::read_yaml()`. The structure
#'   should be `data[[substrate]][[sec_level]][[solution]]` where
#'   `sec_level` is typically either a focal strain name or `"min"`.
#'
#' @return A tibble with columns `substrate`, `sec_level`, `n_sol`, and
#'   `n_zero_growth`.
#'
#' @seealso [importMisosoup()] for the full import.
#'
#' @export
#'
#' @examples
#' \donttest{
#' # raw <- yaml::read_yaml("path/to/misosoup.yaml")
#' # overviewMisosoup(raw)
#' }
overviewMisosoup <- function(data) {
    rows <- list()
    for (sub in names(data)) {
        for (sec in names(data[[sub]])) {
            sols <- data[[sub]][[sec]]
            n_zero <- sum(vapply(
                sols,
                \(s) is.null(s$community) || length(s$community) == 0,
                logical(1)
            ))
            rows[[length(rows) + 1]] <- tibble::tibble(
                substrate = sub,
                sec_level = sec,
                n_sol = length(sols),
                n_zero_growth = n_zero
            )
        }
    }
    dplyr::bind_rows(rows)
}


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

#' Build CMs from a raw MiSoSoup list
#'
#' Iterates over `raw[[substrate]][[sec]][[sol]]`, parses each viable
#' solution, and constructs a [ConsortiumMetabolism-class] object.
#' Zero-growth solutions (empty community) are silently skipped.
#'
#' @param raw Nested list from `yaml::read_yaml()`.
#' @param normalize_ids Whether to normalize metabolite IDs.
#' @return Named list of CM objects.
#' @keywords internal
.buildCMsFromMisosoup <- function(raw, normalize_ids = TRUE) {
    cm_list <- list()
    for (sub in names(raw)) {
        for (sec in names(raw[[sub]])) {
            sols <- raw[[sub]][[sec]]
            for (i in seq_along(sols)) {
                sol <- sols[[i]]
                if (is.null(sol$community) || length(sol$community) == 0) {
                    next
                }
                if (is.null(sol$solution) || length(sol$solution) == 0) {
                    next
                }
                cons_id <- paste(sub, sec, i, sep = "_")
                parsed <- .parseMisosoupSolution(
                    sol,
                    normalize_ids = normalize_ids
                )
                valid_species <- unique(parsed$consortia$species)
                growth_vec <- parsed$growth[
                    names(parsed$growth) %in% valid_species
                ]
                cm <- ConsortiumMetabolism(
                    parsed$consortia,
                    name = cons_id,
                    growth = growth_vec
                )
                S4Vectors::metadata(cm)$media <- parsed$media
                cm_list[[cons_id]] <- cm
            }
        }
    }
    cm_list
}


#' Parse a single MiSoSoup solution entry
#'
#' Extracts the species list from the `community` block, classifies
#' each reaction in the `solution` block as growth / species-scoped
#' exchange / media-level exchange, and returns the three components as
#' a list.
#'
#' Species-scoped exchanges are matched by stripping `_<species>_i`
#' suffixes against the known species list (from `community`), which is
#' robust to both escaped metabolite names and species IDs containing
#' underscores (the old `_e_` delimiter split failed on both).
#'
#' The growth row key is matched against both `Growth_<species>` (older
#' MiSoSoup runs) and `R_Biomass_<species>` (newer runs).
#'
#' @param sol A single solution entry with `community` and `solution`
#'   fields.
#' @param normalize_ids Whether to decode `__NN__` escapes and strip
#'   `R_EX_`/`_e` from metabolite IDs.
#' @return A list with `consortia` (tibble of species, metabolite, flux
#'   for species-scoped exchanges), `media` (tibble of metabolite, flux
#'   for media-level bounds), and `growth` (named numeric vector).
#' @keywords internal
.parseMisosoupSolution <- function(sol, normalize_ids = TRUE) {
    species_list <- sub("^y_", "", names(sol$community))
    rxns <- names(sol$solution)
    fluxes <- unlist(sol$solution, use.names = FALSE)

    is_growth <- grepl("^(R_Biomass|Growth)_", rxns)
    growth_species <- sub("^(R_Biomass|Growth)_", "", rxns[is_growth])
    growth_vec <- stats::setNames(fluxes[is_growth], growth_species)

    non_growth <- !is_growth
    n <- sum(non_growth)
    species_col <- rep("", n)
    met_col <- rep("", n)
    ng_rxns <- rxns[non_growth]
    ng_flux <- fluxes[non_growth]

    for (sp in species_list) {
        pattern <- paste0("_", sp, "_i$")
        hits <- grepl(pattern, ng_rxns) & species_col == ""
        species_col[hits] <- sp
        met_col[hits] <- sub(pattern, "", ng_rxns[hits])
    }
    is_media <- species_col == ""
    met_col[is_media] <- ng_rxns[is_media]

    if (normalize_ids) {
        met_col <- .normalizeBiggIds(.decodeBiggEscapes(met_col))
    }

    consortia <- tibble::tibble(
        species = species_col[!is_media],
        metabolite = met_col[!is_media],
        flux = ng_flux[!is_media]
    )
    media <- tibble::tibble(
        metabolite = met_col[is_media],
        flux = ng_flux[is_media]
    )
    list(
        consortia = consortia,
        media = media,
        growth = growth_vec
    )
}
