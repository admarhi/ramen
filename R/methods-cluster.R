#' @rdname cluster
#' @export
setMethod("cluster", "ConsortiumMetabolismSet", function(object) {
  dend <- object@OverlapMatrix |>
    ### Give these arguments as options?
    dist(method = "euclidean") |>
    hclust(method = "complete") |>
    as.dendrogram()

  # Get node positions in the dendrogram plot
  node_data <- dendextend::get_nodes_xy(dend) |>
    as.data.frame() |>
    tibble::as_tibble() |>
    dplyr::rename(x = "V1", y = "V2") |>
    dplyr::mutate(original_node_id = dplyr::row_number()) |>
    # Filter out leaves (which have y=0)
    dplyr::filter(y != 0) |>
    dplyr::arrange(dplyr::desc(.data$y)) |>
    dplyr::mutate(node_id = dplyr::row_number())

  object@Dendrogram <- list(dend)
  object@NodeData <- node_data
  object
})
