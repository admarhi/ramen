#' Compare Alignments (Deprecated)
#'
#' This function is deprecated and will be reimplemented in a
#' future version.
#'
#' @param ... Alignment objects.
#' @param names Character vector of names.
#' @param smooth Logical; add smooth line.
#' @param se Logical; add standard error bands.
#' @param min_frac Numeric; minimum fraction.
#' @param max_frac Numeric; maximum fraction.
#'
#' @return None. Raises an error.
#' @export
compareAlignments <- function(
    ...,
    names,
    smooth = FALSE,
    se = FALSE,
    min_frac = NULL,
    max_frac = NULL
) {
    cli::cli_abort(
        paste0(
            "{.fn compareAlignments} is deprecated and will ",
            "be reimplemented in a future version."
        )
    )
}
