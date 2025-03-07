#' @title Set of \code{ConsortiumMetabolism} Objects
#'
#' @param ... Lists or individual \code{ConsortiumMetabolism} objects.
#' @param name Character scalar giving the name of the set.
#' @param desc Optional, short description of the sets purpose.
#'
#' @export
ConsortiumMetabolismSet <- function(
  ...,
  name = NA_character_,
  desc = NA_character_
) {
  ### Include a by = NULL argument that allows the user to read a tibble with
  ### minimum 4 columns in which one gives the name of the consortia. These will
  ### then be used as the name of the consortia.

  cons <- list(...)

  # Check that all list entries are CM objects
  stopifnot(exprs = {
    all(lapply(cons, class) == "ConsortiumMetabolism")
  })

  newConsortiumMetabolismSet(
    Name = name,
    Consortia = cons,
    Description = desc
  )
}

#' @noRd
cms <- ConsortiumMetabolismSet
