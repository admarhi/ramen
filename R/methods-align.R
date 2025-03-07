#' @describeIn align Align a \code{ConsortiumMetabolismSet} Object
#' @param object A \code{ConsortiumMetabolismSet} object.
#' @param name Character scalar giving name of the alignment, if `NULL` inherits
#' from the \code{ConsortiumMetabolismSet} object.
#' @return A \code{ConsortiumMetabolismAlignment} object.
#' @export
setMethod("align", "ConsortiumMetabolismSet", function(object) {
  ### Method to retrieve the Binary matrices from all CM in a CMS?

  # Get the consortium indeces in the list
  cons_ind <- tibble::tibble(
    name = purrr::map_chr(object@Consortia, \(x) x@Name)
  ) |>
    tibble::rowid_to_column("index")

  all_met <- purrr::map(object@Consortia, \(x) tibble::as_tibble(x@colData)) |>
    purrr::reduce(\(x, y) dplyr::full_join(x, y, by = "met")) |>
    dplyr::rename_with(
      .fn = \(x) paste0("cm_", cons_ind$index),
      .cols = dplyr::starts_with("index")
    ) |>
    dplyr::arrange(.data$met, .locale = "C") |>
    tibble::rowid_to_column("met_ind")

  # Get all of the tibbles from the consortia and use to calculate a new sparse
  # matrix for the alignment, based on the new indeces of the metabolites.
  all_edges <-
    purrr::map(object@Consortia, \(x) dplyr::mutate(x@Edges, name = x@Name)) |>
    purrr::map_dfr(\(x) rbind(x)) |>
    dplyr::left_join(dplyr::select(all_met, 1:2), by = c(consumed = "met")) |>
    dplyr::rename(c_ind_alig = "met_ind") |>
    dplyr::left_join(dplyr::select(all_met, 1:2), by = c(produced = "met")) |>
    dplyr::rename(p_ind_alig = "met_ind")

  n_edges <- all_edges |>
    dplyr::select(13:15) |>
    dplyr::reframe(n = dplyr::n(), .by = c("c_ind_alig", "p_ind_alig"))

  n_edges_mat <-
    sparseMatrix(n_edges$c_ind_alig, n_edges$p_ind_alig, x = n_edges$n)

  # This currently creates a simple
  graph_list <- purrr::map(
    .x = unique(all_edges$name),
    .f = \(.n)
      x |>
        dplyr::filter(.data$name == .n) |>
        dplyr::select(14:15) |>
        igraph::graph_from_data_frame(directed = TRUE)
  ) |>
    purrr::set_names(nm = unique(all_edges$name))

  newConsortiumMetabolismAlignment(
    tse = TreeSummarizedExperiment(
      assays = list(Levels = n_edges_mat)
    ),
    Name = object@Name,
    Graphs = graph_list
  )
})
