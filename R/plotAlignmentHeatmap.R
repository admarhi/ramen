#' Plot Alignment Heatmap (Deprecated)
#'
#' This function is deprecated and will be reimplemented in
#' Phase 3.
#'
#' @param object A \code{ConsortiumMetabolismAlignment} object.
#' @param top Numeric; top fraction filter.
#' @param bottom Numeric; bottom fraction filter.
#'
#' @return \code{invisible(NULL)} with a deprecation warning.
#' @export
#'
#' @examples
#' \donttest{
#' # This function is deprecated. Use plot() instead:
#' # plot(cma, type = "heatmap")
#' }
plotAlignmentHeatmap <- function(object, top = NULL, bottom = NULL) {
    cli::cli_warn(
        paste0(
            "{.fn plotAlignmentHeatmap} is deprecated. ",
            "Use {.code plot(object, type = \"heatmap\")} ",
            "instead."
        )
    )
    invisible(NULL)
}
