#' @include AllClasses.R AllGenerics.R
NULL

#' @describeIn pathways Get Pathways From a
#'   \code{ConsortiumMetabolism} Object
#' @param verbose Logical scalar. If \code{FALSE}
#'   (default), returns a concise summary with columns
#'   \code{consumed}, \code{produced}, and
#'   \code{n_species}. If \code{TRUE}, returns the full
#'   pathway data including flux statistics, effective
#'   values, and per-species detail.
#' @export
setMethod(
    "pathways",
    "ConsortiumMetabolism",
    function(object, verbose = FALSE) {
        if (verbose) return(object@Pathways)
        object@Pathways[
            , c("consumed", "produced", "n_species")
        ]
    }
)

#' @describeIn pathways Get Pathways From a
#'   \code{ConsortiumMetabolismSet} Object
#' @param type Character scalar giving the type of
#'   pathways to output.
#' @param quantileCutoff Numeric scalar between 0 and 1
#'   giving the quantile threshold to use for filtering
#'   pathways. For \code{"pan-cons"} and \code{"core"}
#'   types, pathways above \code{1 - quantileCutoff}
#'   are returned. For \code{"niche"} and \code{"aux"}
#'   types, pathways below \code{quantileCutoff} are
#'   returned. Defaults to 0.1 (i.e., top/bottom 10
#'   percent).
#' @export
setMethod(
    "pathways",
    "ConsortiumMetabolismSet",
    function(
        object,
        type = c(
            "all", "pan-cons", "niche", "core", "aux"
        ),
        quantileCutoff = 0.1,
        verbose = FALSE
    ) {
        type <- match.arg(type)

        # Validate quantileCutoff parameter
        if (quantileCutoff <= 0 || quantileCutoff >= 1) {
            cli::cli_abort(
                "{.arg quantileCutoff} must be between \\
                 0 and 1 (exclusive), not \\
                 {.val {quantileCutoff}}."
            )
        }

        # Build unified summary: one row per pathway
        # with both n_cons and n_species
        pathway_summary <- object@Pathways |>
            dplyr::reframe(
                n_species = dplyr::n_distinct(
                    .data$species
                ),
                n_cons = dplyr::n_distinct(
                    .data$cm_name
                ),
                .by = c("consumed", "produced")
            ) |>
            dplyr::arrange(dplyr::desc(.data$n_cons))

        total_cons <- length(object@Consortia)
        total_species <- nrow(species(object))

        result <- if (type == "all") {
            pathway_summary
        } else if (type == "pan-cons") {
            quant <- stats::quantile(
                seq_len(total_cons),
                p = 1 - quantileCutoff
            )
            dplyr::filter(
                pathway_summary,
                .data$n_cons > quant
            ) |>
                dplyr::arrange(
                    dplyr::desc(.data$n_cons)
                )
        } else if (type == "niche") {
            quant <- stats::quantile(
                seq_len(total_cons),
                p = quantileCutoff
            )
            dplyr::filter(
                pathway_summary,
                .data$n_cons < quant
            ) |>
                dplyr::arrange(.data$n_cons)
        } else if (type == "core") {
            quant <- stats::quantile(
                2:total_species,
                p = 1 - quantileCutoff
            )
            dplyr::filter(
                pathway_summary,
                .data$n_species > quant
            ) |>
                dplyr::arrange(
                    dplyr::desc(.data$n_species)
                )
        } else if (type == "aux") {
            quant <- stats::quantile(
                seq_len(total_species),
                p = quantileCutoff
            )
            dplyr::filter(
                pathway_summary,
                .data$n_species < quant
            ) |>
                dplyr::arrange(.data$n_species)
        }

        if (verbose) {
            # Join n_cons back onto full slot data
            n_cons_lookup <- pathway_summary[
                , c("consumed", "produced", "n_cons")
            ]
            dplyr::left_join(
                object@Pathways,
                n_cons_lookup,
                by = c("consumed", "produced")
            ) |>
                dplyr::semi_join(
                    result,
                    by = c("consumed", "produced")
                )
        } else {
            result[
                , c(
                    "consumed", "produced",
                    "n_species", "n_cons"
                )
            ]
        }
    }
)

#' @describeIn pathways Get Pathways From a
#'   \code{ConsortiumMetabolismAlignment} Object
#' @export
setMethod(
    "pathways",
    "ConsortiumMetabolismAlignment",
    function(
        object,
        type = c(
            "all", "shared", "unique", "consensus"
        ),
        verbose = FALSE
    ) {
        type <- match.arg(type)
        alnType <- object@Type

        if (type == "shared") {
            if (alnType != "pairwise") {
                cli::cli_abort(
                    "{.arg type} = {.val shared} is \\
                     only available for pairwise \\
                     alignments, not \\
                     {.val {alnType}}."
                )
            }
            pw <- object@SharedPathways
            if (verbose) return(pw)
            pw[, c("consumed", "produced")]
        } else if (type == "unique") {
            if (alnType != "pairwise") {
                cli::cli_abort(
                    "{.arg type} = {.val unique} is \\
                     only available for pairwise \\
                     alignments, not \\
                     {.val {alnType}}."
                )
            }
            if (verbose) {
                list(
                    query = object@UniqueQuery,
                    reference =
                        object@UniqueReference
                )
            } else {
                list(
                    query = object@UniqueQuery[
                        , c("consumed", "produced")
                    ],
                    reference =
                        object@UniqueReference[
                            , c(
                                "consumed",
                                "produced"
                            )
                        ]
                )
            }
        } else if (type == "consensus") {
            if (alnType != "multiple") {
                cli::cli_abort(
                    "{.arg type} = {.val consensus} \\
                     is only available for multiple \\
                     alignments, not \\
                     {.val {alnType}}."
                )
            }
            pw <- object@ConsensusPathways
            if (verbose) return(pw)
            pw[
                , c(
                    "consumed", "produced",
                    "nConsortia", "proportion"
                )
            ]
        } else {
            # type == "all"
            if (verbose) return(object@Pathways)
            object@Pathways[
                , c("consumed", "produced")
            ]
        }
    }
)
