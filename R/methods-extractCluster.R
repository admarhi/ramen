#' @include AllClasses.R AllGenerics.R
NULL

#' @rdname extractCluster
#' @export
setMethod(
    "extractCluster",
    "ConsortiumMetabolismSet",
    function(
        object,
        node_id,
        name = NA_character_,
        description = NA_character_
    ) {
        ## ---- Input validation -----------------------------------
        ## The empty-NodeData early return must precede the range
        ## check below, otherwise the range message would read
        ## "between 1 and 0" which is unactionable.
        n_nodes <- nrow(object@NodeData)
        if (n_nodes == 0L) {
            cli::cli_abort(
                c(
                    "{.fn extractCluster} requires a \\
                    {.cls ConsortiumMetabolismSet} with at least \\
                    2 consortia (so a dendrogram exists).",
                    "i" = "{.code object@NodeData} has 0 rows."
                )
            )
        }
        if (
            !is.numeric(node_id) ||
                length(node_id) != 1L ||
                !is.finite(node_id) ||
                node_id != as.integer(node_id) ||
                node_id < 1L ||
                node_id > n_nodes
        ) {
            cli::cli_abort(
                c(
                    "{.arg node_id} must be a single integer \\
                    between 1 and {n_nodes}.",
                    "i" = "Inspect valid IDs with \\
                    {.code object@NodeData}."
                )
            )
        }
        node_id <- as.integer(node_id)

        if (!is.character(name) || length(name) != 1L) {
            cli::cli_abort(
                "{.arg name} must be a single character \\
                string (NA allowed)."
            )
        }
        if (
            !is.character(description) ||
                length(description) != 1L
        ) {
            cli::cli_abort(
                "{.arg description} must be a single character \\
                string (NA allowed)."
            )
        }

        # Get the tibble with the node data
        tb <- object@NodeData

        # Get the internal node id
        int_node_id <- tb$original_node_id[
            tb$node_id == node_id
        ]

        # Partition retrieves all leaves under a given node
        cons_names <- dendextend::partition_leaves(
            object@Dendrogram[[1]]
        )[[int_node_id]]

        # Create named list from the CMS
        selected_consortia <- stats::setNames(
            object@Consortia,
            vapply(
                object@Consortia,
                \(x) x@Name,
                character(1)
            )
        )[cons_names]

        # In case no name was specified, create from selection
        if (is.na(name)) {
            name <- paste0(
                "Cluster ",
                node_id,
                " from ",
                object@Name
            )
        }

        # Create a new set with the selected CM
        ConsortiumMetabolismSet(
            selected_consortia,
            name = name
        )
    }
)
