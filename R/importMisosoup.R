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
#' the biomass-reaction name (via a configurable regex that defaults to
#' matching `Growth_<sp>`, `R_Biomass_<sp>`, and `R_R_BIOMASS_<sp>`
#' case-insensitively — MiSoSoup derives this name from the underlying
#' model, so it varies), COBRApy-style `__NN__` escape sequences in
#' reaction IDs (e.g. `__40__` for `(`), and the second-level key of
#' each substrate. That second-level key is either `min` — indicating a
#' **CMSC** (complete minimal supplying community: no viable focal
#' strain because no single member grows alone) — or a focal-strain ID,
#' indicating an **MSC** (minimal supplying community for that strain).
#'
#' Per-consortium metadata set on each returned CM:
#' - `metadata(cm)$media` — media-level exchange fluxes (reactions
#'   without a species suffix), preserving the original bounds.
#' - `metadata(cm)$communityGrowth` — the scalar `community_growth`
#'   summary row from the MiSoSoup solution, or `NA_real_` if absent.
#' - `metadata(cm)$misosoupMode` — `"CMSC"` or `"MSC"`.
#' - `metadata(cm)$focalStrain` — focal-strain ID (character) for MSC
#'   consortia, `NA_character_` for CMSC.
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
#' @param biomassPattern Perl-compatible regex identifying biomass
#'   reaction rows. The pattern must match *and consume* the entire
#'   prefix from the start of the reaction ID up through the separator
#'   before the species ID, because it is used for both detection and
#'   species extraction (via `sub()`). Defaults to
#'   `"(?i)^.*?(?:^|_)(?:biomass|growth)_"`, which covers the three
#'   observed MiSoSoup variants (`Growth_<sp>`, `R_Biomass_<sp>`,
#'   `R_R_BIOMASS_<sp>`) case-insensitively and correctly rejects the
#'   `community_growth` summary row. Override this if your MiSoSoup
#'   output uses a non-standard biomass reaction name, e.g.
#'   `biomassPattern = "(?i)^.*?MYBIO_"`.
#' @param verbose If `TRUE` (default), emit progress messages via
#'   [cli::cli_inform()] and show a progress bar when importing a
#'   directory.
#'
#' @return A [ConsortiumMetabolismSet-class] object containing one
#'   [ConsortiumMetabolism-class] per viable consortium found in the
#'   input. Zero-growth solutions are silently skipped.
#'
#' @seealso [importSmetana()] for the sibling SMETANA import.
#'
#' @export
#'
#' @examples
#' # Build a minimal in-memory MiSoSoup-shaped list and import it. In
#' # practice you would point at a real YAML file or a directory of
#' # YAMLs from a MiSoSoup run.
#' raw <- list(
#'     glucose = list(
#'         min = list(
#'             list(
#'                 community = list(y_sp1 = 1, y_sp2 = 1),
#'                 solution = list(
#'                     R_Biomass_sp1 = 0.4,
#'                     R_Biomass_sp2 = 0.3,
#'                     R_EX_glc__D_e_sp1_i = -10,
#'                     R_EX_ac_e_sp1_i = 5,
#'                     R_EX_ac_e_sp2_i = -4,
#'                     R_EX_co2_e_sp2_i = 3,
#'                     community_growth = 0.7
#'                 )
#'             )
#'         )
#'     )
#' )
#' cms <- importMisosoup(raw, name = "demo", verbose = FALSE)
#' cms
#'
#' # Real-world usage (commented; needs an actual MiSoSoup YAML on disk):
#' # cms <- importMisosoup("path/to/misosoup.yaml")
#' # cms <- importMisosoup("path/to/misosoup_dir/", name = "experiment1")
importMisosoup <- function(
    data,
    name = NULL,
    normalize_ids = TRUE,
    biomassPattern = "(?i)^.*?(?:^|_)(?:biomass|growth)_",
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
                normalize_ids = normalize_ids,
                biomassPattern = biomassPattern
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
            normalize_ids = normalize_ids,
            biomassPattern = biomassPattern
        )
        if (verbose) {
            # nolint start: object_usage_linter.
            n_cm <- length(cm_list)
            f_in <- basename(data)
            # nolint end
            cli::cli_inform(c(
                "i" = "Imported {n_cm} consort{?ium/ia} from {.path {f_in}}."
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
            normalize_ids = normalize_ids,
            biomassPattern = biomassPattern
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
#' @param biomassPattern Perl regex identifying biomass reaction rows.
#' @return Named list of CM objects.
#' @keywords internal
.buildCMsFromMisosoup <- function(
    raw,
    normalize_ids = TRUE,
    biomassPattern = "(?i)^.*?(?:^|_)(?:biomass|growth)_"
) {
    cm_list <- list()
    for (sub in names(raw)) {
        for (sec in names(raw[[sub]])) {
            mode <- if (identical(sec, "min")) "CMSC" else "MSC"
            focal <- if (identical(sec, "min")) NA_character_ else sec
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
                    normalize_ids = normalize_ids,
                    biomassPattern = biomassPattern
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
                S4Vectors::metadata(cm)$communityGrowth <-
                    parsed$communityGrowth
                S4Vectors::metadata(cm)$misosoupMode <- mode
                S4Vectors::metadata(cm)$focalStrain <- focal
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
#' The biomass-reaction name is model-dependent (MiSoSoup carries it
#' through from the underlying GEM), so matching is done with a
#' configurable Perl regex that consumes the full prefix through the
#' species-ID separator. The default handles the three observed
#' variants: `Growth_<sp>`, `R_Biomass_<sp>`, and `R_R_BIOMASS_<sp>`.
#'
#' A `community_growth` summary row (MiSoSoup's community-level growth
#' rate, distinct from per-species biomass fluxes) is captured as a
#' scalar in the returned list rather than treated as a media flux.
#' Non-suffixed `R_EX_<metabolite>` rows remain in media — they are
#' legitimate community-level exchange totals.
#'
#' @param sol A single solution entry with `community` and `solution`
#'   fields.
#' @param normalize_ids Whether to decode `__NN__` escapes and strip
#'   `R_EX_`/`_e` from metabolite IDs.
#' @param biomassPattern Perl regex identifying biomass rows; must
#'   consume the full prefix through the species-ID separator (used
#'   for both detection and species extraction via `sub()`).
#' @return A list with `consortia` (tibble of species, metabolite, flux
#'   for species-scoped exchanges), `media` (tibble of metabolite, flux
#'   for media-level bounds), `growth` (named numeric vector), and
#'   `communityGrowth` (scalar `community_growth` value or `NA_real_`).
#' @keywords internal
.parseMisosoupSolution <- function(
    sol,
    normalize_ids = TRUE,
    biomassPattern = "(?i)^.*?(?:^|_)(?:biomass|growth)_"
) {
    species_list <- sub("^y_", "", names(sol$community))
    rxns <- names(sol$solution)
    fluxes <- unlist(sol$solution, use.names = FALSE)

    # Community-level growth summary row: stash as scalar, drop from
    # further classification. Distinct from per-species biomass rows
    # (which the biomass regex catches) and from R_EX_* community
    # exchange totals (which remain in media).
    is_comm <- rxns == "community_growth"
    community_growth <- if (any(is_comm)) fluxes[is_comm][[1]] else NA_real_
    rxns <- rxns[!is_comm]
    fluxes <- fluxes[!is_comm]

    is_growth <- grepl(biomassPattern, rxns, perl = TRUE)
    growth_species <- sub(biomassPattern, "", rxns[is_growth], perl = TRUE)
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
        growth = growth_vec,
        communityGrowth = community_growth
    )
}
