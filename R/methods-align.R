#' @describeIn align Align a \code{ConsortiumMetabolismSet} Object
#' @param object A \code{ConsortiumMetabolismSet} object.
#' @param name Character scalar giving name of the alignment, if `NULL` inherits
#' from the \code{ConsortiumMetabolismSet} object.
#' @return A \code{ConsortiumMetabolismAlignment} object.
#' @export
setMethod("align", "ConsortiumMetabolismSet", function(object) {
  # Make a tibble with new indeces for all metabolites (might not be needed)
  all_met <- purrr::map(object@Consortia, \(x) tibble::as_tibble(x@colData)) |>
    purrr::reduce(\(x, y) dplyr::full_join(x, y, by = "met")) |>
    dplyr::rename_with(
      .fn = \(x) paste0("cm_", cons_ind$index),
      .cols = dplyr::starts_with("index")
    ) |>
    dplyr::arrange(.data$met, .locale = "C") |>
    dplyr::relocate("met") |>
    tibble::rowid_to_column("met_ind")

  # Retrieve the edges tibbles from all objects and bind the new indeces for met
  all_edges <-
    purrr::map(object@Consortia, \(x) dplyr::mutate(x@Edges, name = x@Name)) |>
    purrr::map_dfr(\(x) rbind(x)) |>
    dplyr::left_join(dplyr::select(all_met, 1:2), by = c(consumed = "met")) |>
    dplyr::rename(c_ind_alig = "met_ind") |>
    dplyr::left_join(dplyr::select(all_met, 1:2), by = c(produced = "met")) |>
    dplyr::rename(p_ind_alig = "met_ind")

  # Count the number of consortia for each edge
  n_edges <- all_edges |>
    dplyr::select(13:15) |>
    dplyr::reframe(n = dplyr::n(), .by = c("c_ind_alig", "p_ind_alig"))

  # Make a matrix to create graphs
  n_edges_mat <-
    sparseMatrix(n_edges$c_ind_alig, n_edges$p_ind_alig, x = n_edges$n)

  ### Review this
  graph_list <- purrr::map(
    .x = unique(all_edges$name),
    .f = \(.n)
      all_edges |>
        dplyr::filter(.data$name == .n) |>
        dplyr::select(1:2) |>
        igraph::graph_from_data_frame(directed = TRUE)
  ) |>
    purrr::set_names(nm = unique(all_edges$name))

  # Get the consortia's binary matrices, could be method
  mat_list <- purrr::map(comp_set@Consortia, \(x) assays(x)$Binary)

  # Get vector from 1 to length of the list
  ll <- seq_len(length(mat_list))

  tb <- dplyr::bind_cols(
    tidyr::expand_grid(cm1 = mat_list, cm2 = mat_list),
    tidyr::expand_grid(x = ll, y = ll)
  ) |>
    dplyr::filter(x < y) |>
    dplyr::rowwise() |>
    dplyr::mutate(overlap_score = bin_mat_overlap(.data$cm1, .data$cm2))

  mat <- sparseMatrix(tb$x, tb$y, x = tb$overlap_score) |>
    Matrix::as.matrix()

  dend <- mat |>
    dist(method = "euclidean") |>
    hclust(method = "complete") |>
    as.dendrogram()

  newConsortiumMetabolismAlignment(
    TreeSummarizedExperiment(
      assays = list(Levels = n_edges_mat)
    ),
    Name = object@Name,
    Graphs = graph_list
  )
})


#' @export
bin_mat_overlap <- function(bm1, bm2) {
  # Get the intersection of the metabolites
  met <- intersect(rownames(bm1), rownames(bm2))

  # Calculate the intersection
  int <- bm1[met, met] * bm2[met, met]

  # Calculate the overlap
  sum(int) / min(sum(bm1), sum(bm2))
}
