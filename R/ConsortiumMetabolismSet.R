#' @title Set of \code{ConsortiumMetabolism} Objects
#'
#' @param ... Lists or individual \code{ConsortiumMetabolism} objects.
#' @param name Character scalar giving the name of the set.
#'
#' @export
#'
ConsortiumMetabolismSet <- function(
  ...,
  name = NA_character_
) {
  coms <- list(...)
  # Check that all list entries are MF objects
  stopifnot(exprs = {
    all(lapply(coms, class) == "ConsortiumMetabolism")
  })
}
