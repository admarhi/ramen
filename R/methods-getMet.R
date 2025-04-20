#' @rdname getMet
setMethod("getMet", "ConsortiumMetabolism", function(object) {
  object@Metabolites
})

#' @rdname getMet
setMethod("getMet", "ConsortiumMetabolismSet", function(object) {
  purrr::map2(
    .x = purrr::map(object@Consortia, \(x) tibble::as_tibble(x@colData)),
    .y = purrr::map_chr(object@Consortia, \(x) x@Name),
    .f = \(x, y) dplyr::mutate(x, consortium = y)
  ) |>
    purrr::reduce(\(x, y) dplyr::bind_rows(x, y)) |>
    dplyr::pull("met") |>
    unique() |>
    sort()

  ### Should go in another method bc the output format differs
  # all_met |>
  #   # Change the value to TRUE to allow for easy filtering
  #   dplyr::mutate(index = TRUE) |>
  #   tidyr::pivot_wider(names_from = "consortium", values_from = "index") |>
  #   dplyr::arrange(.data$met, .locale = "C")
})

#' @rdname getMet
setMethod("getMet", "ConsortiumMetabolismAlignment", function(object) {
})
