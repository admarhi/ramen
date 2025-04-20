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
  args <- list(...)
  cons <- unlist(args, recursive = FALSE, use.names = FALSE)

  # Check that all list entries are CM objects
  stopifnot(exprs = {
    all(lapply(cons, class) == "ConsortiumMetabolism")
  })

  # Get binary matrices in a tibble to create the combinations
  bm_tb <- purrr::set_names(
    purrr::map(cons, \(x) assays(x)$Binary),
    purrr::map_chr(cons, \(x) x@Name)
  ) |>
    tibble::enframe() |>
    tibble::rowid_to_column("ind")

  # Create unique combinations of the consortia to reduce computation time
  tb <- tidyr::expand_grid(x = bm_tb$ind, y = bm_tb$ind) |>
    dplyr::filter(x <= y) |>
    dplyr::left_join(bm_tb, by = c("x" = "ind")) |>
    dplyr::left_join(bm_tb, by = c("y" = "ind")) |>
    dplyr::select(
      x,
      name_x = "name.x",
      cm_x = "value.x",
      y,
      name_y = "name.y",
      cm_y = "value.y"
    ) |>
    dplyr::rowwise() |>
    # Calculate the overlap score for all combinations
    dplyr::mutate(overlap_score = bin_mat_overlap(.data$cm_x, .data$cm_y))

  # Create sparse matrix because not all combinations exist
  overlap_matrix <- Matrix::sparseMatrix(
    tb$x,
    tb$y,
    x = tb$overlap_score,
    dimnames = list(bm_tb$name, bm_tb$name)
  ) |>
    # Convert to dense matrix so we can use it for dendrogram
    Matrix::as.matrix()

  newConsortiumMetabolismSet(
    Name = name,
    Consortia = cons,
    Description = desc,
    OverlapMatrix = overlap_matrix
  )
}

#' @noRd
cms <- ConsortiumMetabolismSet
