#' @title Pivot \code{ConsortiumMetabolism} Input Data
#'
#' @description
#' Wrapper function around tidyr's `pivot_longer()` function to
#' transform wide-format community data (one column per
#' direction) into the long format required by
#' \code{\link{ConsortiumMetabolism}}.
#'
#' @param tb A data.frame or tibble with one row per
#'   species-metabolite pair in wide format (separate
#'   columns for consumed and produced metabolites).
#' @param species Column name of the species column.
#' @param from Column name specifying the metabolite
#'   consumed.
#' @param to Column name specifying the metabolite
#'   produced.
#' @param flux Column name of the flux column.
#'
#' @return A \code{\link[tibble]{tibble}} with three columns:
#'   \describe{
#'     \item{species}{Character, the species identifier.}
#'     \item{met}{Character, the metabolite name.}
#'     \item{flux}{Numeric, the metabolic flux (negative for
#'       consumption, positive for production).}
#'   }
#'
#' @export
#'
#' @examples
#' tb <- data.frame(
#'   uptake = c("m1", "m2", "m3", "m4"),
#'   secretion = c("m2", "m3", "m4", "m1"),
#'   flux = c(1, 1, 1, 1),
#'   species = c("s1", "s2", "s3", "s4")
#' )
#'
#' pivotCM(
#'   tb = tb,
#'   species = "species",
#'   from = "uptake",
#'   to = "secretion",
#'   flux = "flux"
#' )
pivotCM <- function(tb, species, from, to, flux) {
    if (!is.data.frame(tb)) {
        cli::cli_abort(
            "{.arg tb} must be a data.frame, not {.cls {class(tb)}}."
        )
    }
    required <- c(species, from, to, flux)
    missing_cols <- setdiff(required, names(tb))
    if (length(missing_cols) > 0) {
        cli::cli_abort(
            "Column{?s} {.val {missing_cols}} not found
            in {.arg tb}."
        )
    }
    tb |>
        tidyr::pivot_longer(cols = dplyr::all_of(c(from, to))) |>
        dplyr::rename(
            met = "value",
            species = {{ species }},
            flux = {{ flux }}
        ) |>
        dplyr::mutate(
            flux = dplyr::if_else(
                .data$name == from,
                .data$flux * -1,
                .data$flux
            )
        ) |>
        dplyr::select("species", "met", "flux")
}
