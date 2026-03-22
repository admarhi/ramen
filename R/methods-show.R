#' Show Method for \code{ConsortiumMetabolism} Object
#'
#' @param object An object of class \code{ConsortiumMetabolism}
#' @export
setMethod("show", "ConsortiumMetabolism", function(object) {
    weight_label <- if (object@Weighted) {
        "Weighted"
    } else {
        "Unweighted"
    }
    n_met <- length(object@Metabolites)
    cli::cli_h3("ConsortiumMetabolism")
    cli::cli_text("Name: {.val {object@Name}}")
    cli::cli_text(
        "{weight_label} metabolic network with ",
        "{.val {n_met}} metabolites."
    )
    invisible(object)
})

#' Show Method for \code{ConsortiumMetabolismSet} Object
#'
#' @param object An object of class \code{ConsortiumMetabolismSet}
#' @export
setMethod("show", "ConsortiumMetabolismSet", function(object) {
    n_cons <- length(object@Consortia)
    cli::cli_h3("ConsortiumMetabolismSet")
    cli::cli_text("Name: {.val {object@Name}}")
    cli::cli_text(
        "Containing {.val {n_cons}} consortia."
    )
    cli::cli_text(
        "Description: {.val {object@Description}}"
    )
    invisible(object)
})

#' Show method for \code{ConsortiumMetabolismAlignment} Objects
#'
#' @param object A \code{ConsortiumMetabolismAlignment} object.
#' @export
setMethod("show", "ConsortiumMetabolismAlignment", function(object) {
    ## Helper: TRUE if slot is a length-1 non-NA value
    .hasValue <- function(x) length(x) == 1L && !is.na(x)

    type_label <- if (.hasValue(object@Type)) {
        object@Type
    } else {
        "uninitialized"
    }
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
    }
    if (.hasValue(object@Type) && object@Type == "multiple") {
        n <- nrow(object@SimilarityMatrix)
        cli::cli_text("Consortia: {.val {n}}")
    }
    invisible(object)
})
