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

    # Partition retrieves all leaves under a given node
    cons_names <- dendextend::partition_leaves(object@Dendrogram[[1]])[[
      int_node_id
    ]]

    # Create named list from the CMS
    selected_consortia <- purrr::set_names(
      object@Consortia,
      nm = purrr::map_chr(object@Consortia, \(x) x@Name)
    )[cons_names]

    # In case no names were supplied in the CMs, apply based on index
    if (is.na(name)) name <- paste0("Cluster ", node_id, " from ", object@Name)

    # Create a new set with the selected CM
    ConsortiumMetabolismSet(selected_consortia, name = name)
  }
)
