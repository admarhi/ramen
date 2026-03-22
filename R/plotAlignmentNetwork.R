#' Plot Alignment Network (Deprecated)
#'
#' This function is deprecated and will be reimplemented in
#' Phase 3.
#'
#' @param object A \code{ConsortiumMetabolismAlignment} object.
#' @param frac Numeric; minimum fraction to plot.
#'
#' @return None. Raises an error.
#' @export
plotAlignmentNetwork <- function(object, frac) {
    cli::cli_abort(
        paste0(
            "{.fn plotAlignmentNetwork} is deprecated. ",
            "Use {.code plot(object, type = \"network\")} ",
            "instead."
        )
    )
}
