#' @describeIn align Align a \code{ConsortiumMetabolismSet} Object
#' @param object A \code{ConsortiumMetabolismSet} object.
#' @param name Character scalar giving name of the alignment, if `NULL` inherits
#' from the \code{ConsortiumMetabolismSet} object.
#' @return A \code{ConsortiumMetabolismAlignment} object.
#' @export
setMethod("align", "ConsortiumMetabolismSet", function(object) {
  cli::cli_h1(paste0("Aligning ", object@Name))

  cli::cli_status("Getting metabolites")
  # Get indeces of met in each consortium for re-indexing
  all_met <- purrr::map2(
    .x = purrr::map(object@Consortia, \(x) tibble::as_tibble(x@colData)),
    .y = purrr::map_chr(object@Consortia, \(x) x@Name),
    .f = \(x, y) dplyr::mutate(x, consortium = y)
  ) |>
    purrr::reduce(\(x, y) dplyr::bind_rows(x, y)) |>
    dplyr::rename(consortium_ind = "index")
  cli::cli_process_done()
  cli::cli_status("Re-indexing metabolites")

  # Re-index all metabolites
  new_met_ind <- tibble::tibble(met = unique(all_met$met)) |>
    dplyr::arrange(.data$met, .locale = "C") |>
    tibble::rowid_to_column("met_ind")

  # Join the new indeces
  all_met <- all_met |>
    dplyr::left_join(new_met_ind, by = "met") |>
    dplyr::relocate("met_ind", "met", "consortium") |>
    # Still need the wide format here to join into the all_edges tibble
    tidyr::pivot_wider(
      names_from = "consortium",
      values_from = "consortium_ind"
    ) |>
    dplyr::arrange(.data$met_ind)
  cli::cli_process_done()
  cli::cli_status("Getting all edges")
  # Retrieve the edges tibbles from all objects and bind the new indeces for met
  all_edges <-
    purrr::map(object@Consortia, \(x) dplyr::mutate(x@Edges, name = x@Name)) |>
    purrr::map_dfr(\(x) rbind(x)) |>
    dplyr::left_join(dplyr::select(all_met, 1:2), by = c(consumed = "met")) |>
    dplyr::rename(c_ind_alig = "met_ind") |>
    dplyr::left_join(dplyr::select(all_met, 1:2), by = c(produced = "met")) |>
    dplyr::rename(p_ind_alig = "met_ind")
  cli::cli_process_done()
  cli::cli_status("Counting consortia for each edge")
  # Count the number of consortia for each edge
  n_cons <- all_edges |>
    dplyr::select(13:15) |>
    dplyr::reframe(n = dplyr::n(), .by = c("c_ind_alig", "p_ind_alig"))

  # Make a matrix to create graphs
  levels_mat <-
    sparseMatrix(n_cons$c_ind_alig, n_cons$p_ind_alig, x = n_cons$n)
  cli::cli_process_done()
  cli::cli_status("Getting all graphs")
  # Store the graphs in a named list to enable the graph based alignment
  graph_list <- purrr::set_names(
    purrr::map(object@Consortia, \(x) x@Graphs[[1]]),
    nm = purrr::map_chr(object@Consortia, \(x) x@Name)
  )
  cli::cli_process_done()

  invisible(
    return(
      newConsortiumMetabolismAlignment(
        TreeSummarizedExperiment(
          assays = list(Levels = levels_mat)
        ),
        Name = object@Name,
        Graphs = graph_list,
        Metabolites = all_met,
        Dendrogram = list(dend)
      )
    )
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
