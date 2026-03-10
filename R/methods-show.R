#' Show Method for \code{ConsortiumMetabolism} Object
#'
#' @param object An object of class \code{ConsortiumMetabolism}
#' @export
setMethod("show", "ConsortiumMetabolism", function(object) {
    stringr::str_glue(
        "\n{object@Name}: ConsortiumMetabolism Object\n",
        "{ifelse(object@Weighted, 'Weighted', 'Unweighted')} ",
        "metabolic network with {length(object@Metabolites)} metabolites.\n\n"
    ) |>
        cat()
    invisible(object)
})

#' Show Method for \code{ConsortiumMetabolismSet} Object
#'
#' @param object An object of class \code{ConsortiumMetabolismSet}
#' @export
setMethod("show", "ConsortiumMetabolismSet", function(object) {
    stringr::str_glue(
        "\n{object@Name}: ConsortiumMetabolismSet Object\n",
        "containing {length(object@Consortia)} consortia.\n",
        "Description: {object@Description}\n\n"
    ) |>
        cat()
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
