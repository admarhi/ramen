#' @include AllClasses.R AllGenerics.R
NULL

## ---- internal helpers ------------------------------------------------------

## Counts rows in a slot that may be NULL or an empty data.frame.
.nrowOr0 <- function(x) {
    if (is.null(x) || !is.data.frame(x)) {
        return(0L)
    }
    nrow(x)
}

## Format a level histogram as a single inline string. For wide alignments
## (k > maxInline), only the first three and last bucket are printed.
.formatLevelCounts <- function(counts, maxInline = 8L) {
    k <- length(counts)
    if (k == 0L) {
        return("")
    }
    parts <- sprintf("n=%d: %d", seq_along(counts), counts)
    if (k <= maxInline) {
        return(paste(parts, collapse = "  "))
    }
    paste0(
        paste(parts[seq_len(3L)], collapse = "  "),
        "  ...  ",
        parts[k]
    )
}

#' Show Method for \code{ConsortiumMetabolism} Object
#'
#' @param object An object of class \code{ConsortiumMetabolism}
#' @return The object, invisibly.
#' @examples
#' cm <- synCM("test", n_species = 3, max_met = 5)
#' show(cm)
#' @export
setMethod("show", "ConsortiumMetabolism", function(object) {
    # nolint start: object_usage_linter.
    weight_label <- if (object@Weighted) "Weighted" else "Unweighted"
    n_met <- length(object@Metabolites)
    n_sp <- length(unique(object@InputData$species))
    n_pw <- nrow(object@Pathways)
    # nolint end
    cli::cli_h3("ConsortiumMetabolism")
    cli::cli_text("Name: {.val {object@Name}}")
    cli::cli_text(
        "{weight_label} metabolic network: ",
        "{.val {n_sp}} species, ",
        "{.val {n_met}} metabolites, ",
        "{.val {n_pw}} pathways."
    )

    ## Per-species pathway range -- only meaningful with multiple species.
    if (n_sp >= 2L) {
        pw_per_sp <- tryCatch(
            {
                ss <- speciesSummary(object)
                ss$n_pathways
            },
            error = function(e) NULL
        )
        if (!is.null(pw_per_sp) && length(pw_per_sp) > 0L) {
            cli::cli_text(
                "Pathways per species: ",
                "min {.val {min(pw_per_sp)}}, ",
                "mean {.val {round(mean(pw_per_sp), 1)}}, ",
                "max {.val {max(pw_per_sp)}}."
            )
        }
    }
    invisible(object)
})

#' Show Method for \code{ConsortiumMetabolismSet} Object
#'
#' @param object An object of class \code{ConsortiumMetabolismSet}
#' @return The object, invisibly.
#' @examples
#' cm1 <- synCM("comm_1", n_species = 3, max_met = 5)
#' cm2 <- synCM("comm_2", n_species = 4, max_met = 6)
#' cms <- ConsortiumMetabolismSet(cm1, cm2, name = "test")
#' show(cms)
#' @export
setMethod("show", "ConsortiumMetabolismSet", function(object) {
    # nolint start: object_usage_linter.
    n_cons <- length(object@Consortia)
    n_sp <- dplyr::n_distinct(object@Pathways$species)
    n_met <- length(unique(
        c(object@Pathways$consumed, object@Pathways$produced)
    ))
    sp_counts <- vapply(
        object@Consortia,
        function(cm) length(species(cm)),
        integer(1L)
    )
    met_counts <- vapply(
        object@Consortia,
        function(cm) length(setdiff(metabolites(cm), "media")),
        integer(1L)
    )
    # nolint end

    cli::cli_h3("ConsortiumMetabolismSet")
    cli::cli_text("Name: {.val {object@Name}}")
    cli::cli_text(
        "{.val {n_cons}} consortia, ",
        "{.val {n_sp}} species, ",
        "{.val {n_met}} metabolites."
    )
    cli::cli_text(
        "Community size (species): ",
        "min {.val {min(sp_counts)}}, ",
        "mean {.val {round(mean(sp_counts), 1)}}, ",
        "max {.val {max(sp_counts)}}."
    )
    cli::cli_text(
        "Community size (metabolites): ",
        "min {.val {min(met_counts)}}, ",
        "mean {.val {round(mean(met_counts), 1)}}, ",
        "max {.val {max(met_counts)}}."
    )

    ## Per-class pathway breakdown. Quantile annotated so the printed
    ## counts are reproducible (defaults: 0.1 for pathways, 0.15 for
    ## species). Skipped when there is only one consortium since
    ## pan-cons / niche collapse trivially.
    if (n_cons >= 2L) {
        pw_q <- 0.1
        pw_classes <- tryCatch(
            list(
                pan = nrow(pathways(
                    object,
                    type = "pan-cons",
                    quantileCutoff = pw_q
                )),
                niche = nrow(pathways(
                    object,
                    type = "niche",
                    quantileCutoff = pw_q
                )),
                core = nrow(pathways(
                    object,
                    type = "core",
                    quantileCutoff = pw_q
                )),
                aux = nrow(pathways(
                    object,
                    type = "aux",
                    quantileCutoff = pw_q
                ))
            ),
            error = function(e) NULL
        )
        if (!is.null(pw_classes)) {
            cli::cli_text(
                "Pathways: ",
                "{.val {pw_classes$pan}} pan-cons, ",
                "{.val {pw_classes$niche}} niche, ",
                "{.val {pw_classes$core}} core, ",
                "{.val {pw_classes$aux}} aux ",
                "(quantile = {.val {pw_q}})."
            )
        }
    }

    ## Generalist / specialist split needs at least a handful of species
    ## for the 15% quantile to mean anything.
    if (n_sp >= 3L) {
        sp_q <- 0.15
        sp_classes <- tryCatch(
            list(
                generalists = length(species(
                    object,
                    type = "generalists",
                    quantileCutoff = sp_q
                )),
                specialists = length(species(
                    object,
                    type = "specialists",
                    quantileCutoff = sp_q
                ))
            ),
            error = function(e) NULL
        )
        if (!is.null(sp_classes)) {
            cli::cli_text(
                "Species: ",
                "{.val {sp_classes$generalists}} generalists, ",
                "{.val {sp_classes$specialists}} specialists ",
                "(quantile = {.val {sp_q}})."
            )
        }
    }

    if (!is.na(object@Description)) {
        cli::cli_text(
            "Description: {.val {object@Description}}"
        )
    }
    invisible(object)
})

#' Show method for \code{ConsortiumMetabolismAlignment} Objects
#'
#' @param object A \code{ConsortiumMetabolismAlignment} object.
#' @return The object, invisibly.
#' @examples
#' cm1 <- synCM("comm_1", n_species = 3, max_met = 5)
#' cm2 <- synCM("comm_2", n_species = 4, max_met = 6)
#' cma <- align(cm1, cm2)
#' show(cma)
#' @export
setMethod("show", "ConsortiumMetabolismAlignment", function(object) {
    ## Helper: TRUE if slot is a length-1 non-NA value
    .hasValue <- function(x) length(x) == 1L && !is.na(x)

    # nolint next: object_usage_linter.
    type_label <- if (.hasValue(object@Type)) object@Type else "uninitialized"
    cli::cli_h3("ConsortiumMetabolismAlignment")
    cli::cli_text("Name: {.val {object@Name}}")
    cli::cli_text("Type: {.val {type_label}}")
    cli::cli_text("Metric: {.val {object@Metric}}")
    if (.hasValue(object@PrimaryScore)) {
        cli::cli_text(
            "Score: {.val {round(object@PrimaryScore, 4)}}"
        )
    }
    if (.hasValue(object@Pvalue)) {
        cli::cli_text(
            "P-value: {.val {format(object@Pvalue, digits = 3)}}"
        )
    }
    if (.hasValue(object@Type) && object@Type == "pairwise") {
        cli::cli_text(
            "Query: {.val {object@QueryName}}, ",
            "Reference: {.val {object@ReferenceName}}"
        )
        cq <- object@Scores$coverageQuery
        cr <- object@Scores$coverageReference
        if (!is.null(cq) && !is.null(cr)) {
            cli::cli_text(
                "Coverage: query ",
                "{.val {round(cq, 3)}}, ",
                "reference ",
                "{.val {round(cr, 3)}}"
            )
        }
        ## Pairwise pathway intersection: shared / unique-each-side.
        ## This is the n=1 / n=2 case of the multi-level hierarchy.
        n_shared <- .nrowOr0(object@SharedPathways)
        n_uq <- .nrowOr0(object@UniqueQuery)
        n_ur <- .nrowOr0(object@UniqueReference)
        if (n_shared + n_uq + n_ur > 0L) {
            cli::cli_text(
                "Pathways: ",
                "{.val {n_shared}} shared, ",
                "{.val {n_uq}} query-only, ",
                "{.val {n_ur}} reference-only."
            )
        }
    }
    if (.hasValue(object@Type) && object@Type == "multiple") {
        n <- nrow(object@SimilarityMatrix) # nolint: object_usage_linter.
        cli::cli_text("Consortia: {.val {n}}")

        ## Multi-level hierarchy: pathways at n=1..n=k consortia.
        prev <- object@Prevalence
        if (
            is.data.frame(prev) &&
                nrow(prev) > 0L &&
                "nConsortia" %in% names(prev) &&
                n > 0L
        ) {
            counts <- tabulate(
                as.integer(prev$nConsortia),
                nbins = n
            )
            n_core <- counts[n] # nolint: object_usage_linter.
            cli::cli_text(
                "Core (shared by all {.val {n}}): ",
                "{.val {n_core}} pathways."
            )
            cli::cli_text(
                "Levels: ",
                .formatLevelCounts(counts) # nolint
            )
        }
    }
    if (.hasValue(object@Type) && object@Type == "search") {
        nDb <- ncol(object@SimilarityMatrix) # nolint: object_usage_linter.
        cli::cli_text(
            "Query: {.val {object@QueryName}}, ",
            "Top hit: {.val {object@ReferenceName}} ",
            "(of {.val {nDb}} consortia)"
        )
    }
    invisible(object)
})
