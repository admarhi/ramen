#' @include AllClasses.R AllGenerics.R
NULL

## ---- CMS accessor methods --------------------------------------------------

#' @describeIn pathways Get Pathways From a
#'   \code{ConsortiumMetabolismSet} Object
#' @param type Character scalar giving the type of
#'   pathways to output. Note that the four filtered
#'   types use two different rulers:
#'   \code{"pan-cons"} and \code{"niche"} rank by the
#'   integer per-pathway consortium count and quantile
#'   over a uniform integer ruler
#'   (\code{seq_len(total_cons)}); \code{"core"} and
#'   \code{"aux"} rank by per-pathway species count and
#'   quantile over the empirical \code{n_species}
#'   distribution across pathways.
#' @param quantileCutoff Numeric scalar between 0 and 1
#'   giving the quantile threshold to use for filtering
#'   pathways. For \code{"pan-cons"} and \code{"core"}
#'   types, pathways strictly above \code{1 - quantileCutoff}
#'   are returned. For \code{"niche"} and \code{"aux"}
#'   types, pathways at or below \code{quantileCutoff}
#'   are returned (the boundary is included on the lower
#'   tail so that ties at the quantile floor are not
#'   silently dropped). Defaults to 0.1 (i.e., top/bottom
#'   10 percent).
#' @export
setMethod(
    "pathways",
    "ConsortiumMetabolismSet",
    function(
        object,
        type = c(
            "all",
            "pan-cons",
            "niche",
            "core",
            "aux"
        ),
        quantileCutoff = 0.1,
        verbose = FALSE
    ) {
        type <- match.arg(type)

        # Validate quantileCutoff parameter
        if (quantileCutoff <= 0 || quantileCutoff >= 1) {
            cli::cli_abort(
                "{.arg quantileCutoff} must be \\
                between 0 and 1 (exclusive), \\
                not {.val {quantileCutoff}}."
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
            # Use <= so ties at the floor (e.g. n_cons = 1)
            # are included; without this, lower-tail filters
            # collapse to empty when many pathways tie at
            # the minimum.
            dplyr::filter(
                pathway_summary,
                .data$n_cons <= quant
            ) |>
                dplyr::arrange(.data$n_cons)
        } else if (type == "core") {
            quant <- stats::quantile(
                pathway_summary$n_species,
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
                pathway_summary$n_species,
                p = quantileCutoff
            )
            # Use <= so ties at the floor (very common when
            # many pathways have n_species = 1) are
            # included; otherwise the filter is empty.
            dplyr::filter(
                pathway_summary,
                .data$n_species <= quant
            ) |>
                dplyr::arrange(.data$n_species)
        }

        if (verbose) {
            # Join n_cons back onto full slot data
            n_cons_lookup <- pathway_summary[,
                c("consumed", "produced", "n_cons")
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
            result[,
                c(
                    "consumed",
                    "produced",
                    "n_species",
                    "n_cons"
                )
            ]
        }
    }
)

#' @rdname overlapMatrix
#' @exportMethod overlapMatrix
setMethod(
    "overlapMatrix",
    "ConsortiumMetabolismSet",
    function(object) {
        object@OverlapMatrix
    }
)

#' @rdname metabolites
setMethod("metabolites", "ConsortiumMetabolismSet", function(object) {
    sort(object@Metabolites$met)
})

#' @describeIn consortia Get the list of
#'   \code{ConsortiumMetabolism} objects from a
#'   \code{ConsortiumMetabolismSet}.
#' @param object A \code{ConsortiumMetabolismSet} object.
#' @export
setMethod(
    "consortia",
    "ConsortiumMetabolismSet",
    function(object) object@Consortia
)

#' @title Coerce a ConsortiumMetabolismSet to a data.frame
#'
#' @description
#' Row-binds the per-consortium edge lists of every
#' \code{ConsortiumMetabolism} in the set, prefixing each
#' row with a \code{consortium} column.
#'
#' @param x A \code{\linkS4class{ConsortiumMetabolismSet}}
#'   object.
#' @param row.names Ignored.
#' @param optional Ignored.
#' @param ... Additional arguments (currently unused).
#'
#' @return A \code{data.frame} with columns
#'   \code{consortium}, \code{met}, \code{species}, and
#'   \code{flux}. Empty sets return a 0-row
#'   \code{data.frame} with the same column names.
#'
#' @examples
#' cm1 <- synCM("a", n_species = 3, max_met = 5)
#' cm2 <- synCM("b", n_species = 3, max_met = 5)
#' cms <- ConsortiumMetabolismSet(list(cm1, cm2), name = "demo")
#' head(as.data.frame(cms))
#'
#' @rdname as.data.frame-ConsortiumMetabolismSet
#' @export
setMethod(
    "as.data.frame",
    "ConsortiumMetabolismSet",
    function(x, ...) {
        cms_list <- x@Consortia
        if (length(cms_list) == 0L) {
            return(data.frame(
                consortium = character(),
                met = character(),
                species = character(),
                flux = numeric(),
                stringsAsFactors = FALSE
            ))
        }
        nm <- vapply(
            cms_list,
            \(cm) cm@Name,
            character(1L)
        )
        per_cm <- lapply(seq_along(cms_list), function(i) {
            df <- as.data.frame(cms_list[[i]])
            cbind(
                consortium = nm[[i]],
                df,
                stringsAsFactors = FALSE
            )
        })
        do.call(rbind, per_cm)
    }
)

#' @describeIn growth Get Growth Rates From a
#'   \code{ConsortiumMetabolismSet}
#' @export
setMethod("growth", "ConsortiumMetabolismSet", function(object) {
    stats::setNames(
        lapply(object@Consortia, growth),
        vapply(
            object@Consortia,
            \(x) x@Name,
            character(1L)
        )
    )
})

#' @describeIn species Return Species in a
#'   \code{ConsortiumMetabolismSet}, optionally filtered to
#'   metabolic generalists or specialists by pathway count.
#' @param object A \code{ConsortiumMetabolismSet} object.
#' @param type Character scalar. One of \code{"all"}
#'   (default), \code{"generalists"} (top fraction by
#'   pathway count), or \code{"specialists"} (bottom
#'   fraction).
#' @param quantileCutoff Numeric scalar in (0, 1) giving
#'   the fraction of species to label as generalists or
#'   specialists. Defaults to \code{0.15}.
setMethod(
    "species",
    "ConsortiumMetabolismSet",
    function(
        object,
        type = c("all", "generalists", "specialists"),
        quantileCutoff = 0.15,
        ...
    ) {
        type <- match.arg(type)
        all_species <- sort(unique(object@Pathways$species))

        if (type == "all") {
            return(all_species)
        }

        if (quantileCutoff <= 0 || quantileCutoff >= 1) {
            cli::cli_abort(
                "{.arg quantileCutoff} must be between 0 \\
                and 1 (exclusive), not \\
                {.val {quantileCutoff}}."
            )
        }

        n_pathways_per_species <- object@Pathways |>
            dplyr::mutate(
                pathway_name = paste0(
                    .data$consumed,
                    "-",
                    .data$produced
                )
            ) |>
            dplyr::reframe(
                n_pathways = dplyr::n_distinct(
                    .data$pathway_name
                ),
                .by = "species"
            )

        counts <- n_pathways_per_species$n_pathways
        names(counts) <- n_pathways_per_species$species

        if (type == "generalists") {
            quant <- stats::quantile(
                counts,
                p = 1 - quantileCutoff
            )
            sort(names(counts)[counts >= quant])
        } else {
            quant <- stats::quantile(
                counts,
                p = quantileCutoff
            )
            sort(names(counts)[counts <= quant])
        }
    }
)

#' @describeIn speciesSummary Species summary for a
#'   \code{ConsortiumMetabolismSet}
#' @param object A \code{ConsortiumMetabolismSet} object.
#' @export
setMethod(
    "speciesSummary",
    "ConsortiumMetabolismSet",
    function(object, ...) {
        object@Pathways |>
            dplyr::mutate(
                pathway_name = paste0(
                    .data$consumed,
                    "-",
                    .data$produced
                )
            ) |>
            dplyr::reframe(
                n_consortia = dplyr::n_distinct(
                    .data$cm_name
                ),
                n_pathways = dplyr::n_distinct(
                    .data$pathway_name
                ),
                .by = "species"
            ) |>
            dplyr::arrange(dplyr::desc(.data$n_pathways))
    }
)

#' @describeIn filterConsortia Filter consortia from a
#'   \code{ConsortiumMetabolismSet}
#' @param object A \code{ConsortiumMetabolismSet} object.
#' @param i Integer vector, character vector of consortium
#'   names, or logical vector.
#' @importFrom S4Vectors metadata
#' @export
setMethod(
    "filterConsortia",
    "ConsortiumMetabolismSet",
    function(object, i) {
        n <- length(object@Consortia)
        cm_names <- vapply(
            object@Consortia,
            \(x) x@Name,
            character(1L)
        )

        idx <- if (is.logical(i)) {
            if (length(i) != n) {
                cli::cli_abort(
                    "Logical {.arg i} must have length \\
                    {n}, not {length(i)}."
                )
            }
            which(i)
        } else if (is.character(i)) {
            not_found <- setdiff(i, cm_names)
            if (length(not_found) > 0L) {
                cli::cli_abort(
                    "Consortium name{?s} not found: \\
                    {.val {not_found}}."
                )
            }
            match(i, cm_names)
        } else if (is.numeric(i)) {
            idx_int <- as.integer(i)
            oob <- idx_int[idx_int < 1L | idx_int > n]
            if (length(oob) > 0L) {
                cli::cli_abort(
                    "Index {.val {oob}} out of bounds \\
                    (1 to {n})."
                )
            }
            idx_int
        } else {
            cli::cli_abort(
                "{.arg i} must be integer, character, \\
                or logical, not {.cls {class(i)}}."
            )
        }

        selected <- object@Consortia[idx]
        lnk <- metadata(object)$linkage
        if (is.null(lnk)) {
            lnk <- "complete"
        }
        ConsortiumMetabolismSet(
            selected,
            name = object@Name,
            desc = object@Description,
            linkage = lnk,
            verbose = FALSE
        )
    }
)
