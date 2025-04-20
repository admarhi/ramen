#' @rdname getCluster
#' @export
setMethod(
  "getCluster",
  "ConsortiumMetabolismSet",
  function(object, node_id, name = NA_character_, description = NA_character_) {
    # Get the tibble with the node data
    tb <- object@NodeData

    # Get the internal node id
    int_node_id <- tb$original_node_id[tb$node_id == node_id]

    cons_names <- partition_leaves(object@Dendrogram[[1]])[[int_node_id]]

    selected_consortia <- purrr::set_names(
      object@Consortia,
      nm = purrr::map_chr(object@Consortia, \(x) x@Name)
    )[cons_names]

    if (is.na(name)) name <- paste0("Cluster ", node_id, " from ", object@Name)
    newConsortiumMetabolismSet(
      Name = name,
      Consortia = selected_consortia,
      OverlapMatrix = object@OverlapMatrix[cons_names, cons_names]
    )
  }
)
