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
  cli::cli_h1(paste0("Creating CMS: ", name))

  # ---- Checking args ---------------------------------------------------------
  ### Include a by = NULL argument that allows the user to read a tibble with
  ### minimum 4 columns in which one gives the name of the consortia. These will
  ### then be used as the name of the consortia.
  cli::cli_status("Checking arguments")
  args <- list(...)
  cons <- unlist(args, recursive = FALSE, use.names = FALSE)
  stopifnot(exprs = {
    all(lapply(cons, class) == "ConsortiumMetabolism")
  })
  cli::cli_process_done()

  # ---- Get all mets ----------------------------------------------------------
  cli::cli_status("Getting metabolites")
  # Get indeces of met in each consortium for re-indexing
  all_met <- purrr::map2(
    .x = purrr::map(cons, \(x) tibble::as_tibble(x@colData)),
    .y = purrr::map_chr(cons, \(x) x@Name),
    .f = \(x, y) dplyr::mutate(x, consortium = y)
  ) |>
    purrr::reduce(\(x, y) dplyr::bind_rows(x, y)) |>
    dplyr::rename(consortium_ind = "index")
  cli::cli_process_done()

  # ---- Create new indeces for metabolites ------------------------------------
  # Create an ordered tibble to re-index metabolites. This tibble can then be
  # joined to the all_met tibble to get the new inds for c and p in each edge.
  cli::cli_status("Re-indexing metabolites")
  new_met_ind <- tibble::tibble(met = unique(all_met$met)) |>
    dplyr::arrange(.data$met, .locale = "C") |>
    tibble::rowid_to_column("met_ind")

  # Join the new indeces to
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

  # ---- Get all edges ---------------------------------------------------------
  # Retrieve the edges tibbles from all objects and bind the new indeces for met
  cli::cli_status("Getting all edges")
  all_edges <-
    purrr::map(
      cons,
      \(x) dplyr::mutate(x@Edges, cm_name = x@Name)
    ) |>
    ### Use reduce and dplyr::bind_rows?
    purrr::map_dfr(\(x) rbind(x)) |>
    dplyr::left_join(dplyr::select(all_met, 1:2), by = c(consumed = "met")) |>
    dplyr::rename(c_ind_alig = "met_ind") |>
    dplyr::left_join(dplyr::select(all_met, 1:2), by = c(produced = "met")) |>
    dplyr::rename(p_ind_alig = "met_ind") |>
    tidyr::unnest("data") |>
    tidyr::unnest("c_prob") |>
    tidyr::unnest("p_prob")
  cli::cli_process_done()

  # ---- Levels Matrix ---------------------------------------------------------
  # Count the number of consortia for each edge to allow easy creation of the
  # levels matrix. This serves as the basis for visualisation and evaluation of
  # the different levels of the alignment object later.
  cli::cli_status("Counting consortia for each edge")
  n_cons <- all_edges |>
    dplyr::reframe(
      n = dplyr::n_distinct(.data$cm_name),
      .by = c("c_ind_alig", "p_ind_alig")
    )
  # Make a sparse matrix as not all edges exist, convert to dense for graph use.
  levels_mat <- sparseMatrix(
    n_cons$c_ind_alig,
    n_cons$p_ind_alig,
    x = n_cons$n,
    dims = c(nrow(new_met_ind), nrow(new_met_ind)),
    dimnames = list(new_met_ind$met, new_met_ind$met)
  ) |>
    Matrix::as.matrix()
  cli::cli_process_done()

  # ---- Getting binary matrices -----------------------------------------------
  # Store the binary matrices in a tibble with the names and indeces. This is
  # used to create only unique combinations between consortia in order to min
  # computation time.
  cli::cli_status("Getting binary matrices")
  bm_tb <- purrr::set_names(
    purrr::map(cons, \(x) assays(x)$Binary),
    purrr::map_chr(cons, \(x) x@Name)
  ) |>
    ### Here I should have the reindexing step
    tibble::enframe() |>
    tibble::rowid_to_column("ind")
  cli::cli_process_done()

  # ---- Overlap score ---------------------------------------------------------
  cli::cli_status("Calculating overlap scores")
  tb <-
    # Create unique combinations of the consortia to reduce computation time
    tidyr::expand_grid(x = bm_tb$ind, y = bm_tb$ind) |>
    dplyr::filter(.data$x <= .data$y) |>
    # Left join the binary matrices by new index
    dplyr::left_join(bm_tb, by = c("x" = "ind")) |>
    dplyr::left_join(bm_tb, by = c("y" = "ind")) |>
    dplyr::select(
      "x",
      name_x = "name.x",
      cm_x = "value.x",
      "y",
      name_y = "name.y",
      cm_y = "value.y"
    ) |>
    ### So if I had the metabolites matrices re-indexed here it would mean that
    ### I can simply multiply them and don't need to rely on the helper function
    ### Likely performance increase.
    dplyr::mutate(
      overlap_score = purrr::map2_dbl(
        .data$cm_x,
        .data$cm_y,
        \(x, y) .binMatOverlap(x, y),
        .progress = TRUE
      )
    )
  cli::cli_process_done()

  # Create sparse matrix because not all combinations exist
  cli::cli_status("Create new edge matrix")
  overlap_matrix <- Matrix::sparseMatrix(
    tb$x,
    tb$y,
    x = 1 - tb$overlap_score,
    dimnames = list(bm_tb$name, bm_tb$name)
  ) |>
    # Convert to dense matrix so we can use it for dendrogram
    Matrix::as.matrix()
  cli::cli_process_done()

  # ---- Dendrogram ------------------------------------------------------------
  # Make the dendrogram
  cli::cli_status("Creating dendrogram")
  dend <-
    stats::dist(overlap_matrix) |>
    stats::hclust() |>
    stats::as.dendrogram()
  cli::cli_process_done()

  # ---- Dendrogram Node Data --------------------------------------------------
  # Get node positions in the dendrogram and index without leaves s.t. they can
  # be used later to retrieve individual clusters.
  cli::cli_status("Calculating node data")
  node_data <- dendextend::get_nodes_xy(dend) |>
    as.data.frame() |>
    tibble::as_tibble() |>
    dplyr::rename(x = "V1", y = "V2") |>
    dplyr::mutate(original_node_id = dplyr::row_number()) |>
    # Filter out leaves (which have y=0)
    dplyr::filter(.data$y != 0) |>
    dplyr::arrange(dplyr::desc(.data$y)) |>
    dplyr::mutate(node_id = dplyr::row_number())
  cli::cli_process_done()

  # ---- Graphs ----------------------------------------------------------------
  # Store the graphs in a named list to enable the graph based alignment
  cli::cli_status("Getting all graphs")
  graph_list <- purrr::set_names(
    purrr::map(cons, \(x) x@Graphs[[1]]),
    nm = purrr::map_chr(cons, \(x) x@Name)
  )
  cli::cli_process_done()

  # ---- Pathways --------------------------------------------------------------

  # ---- cli -------------------------------------------------------------------
  msg <- paste0("ConsortiumMetabolismSet '", name, "' successfully created.")
  cli::cli_alert_success(msg)
  d <- cli::cli_div(theme = list(rule = list(color = "cyan")))
  cli::cli_rule()
  cli::cli_end(d)

  newConsortiumMetabolismSet(
    TreeSummarizedExperiment(
      assays = list(Levels = levels_mat),
      colData = new_met_ind
    ),
    Name = name,
    Consortia = cons,
    Description = desc,
    OverlapMatrix = overlap_matrix,
    Dendrogram = list(dend),
    NodeData = node_data,
    Graphs = graph_list,
    Edges = all_edges,
    Metabolites = all_met
  )
}

#' @noRd
cms <- ConsortiumMetabolismSet

#' @noRd
.binMatOverlap <- function(bm1, bm2) {
  # Get the intersection of the metabolites
  met <- intersect(rownames(bm1), rownames(bm2))

  # Calculate the intersection
  int <- bm1[met, met] * bm2[met, met]

  # Calculate the overlap
  sum(int) / min(sum(bm1), sum(bm2))
}
