#' @include AllClasses.R AllGenerics.R
NULL

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
    }
    if (.hasValue(object@Type) && object@Type == "multiple") {
        n <- nrow(object@SimilarityMatrix) # nolint: object_usage_linter.
        cli::cli_text("Consortia: {.val {n}}")
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
