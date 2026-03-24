#' @title Set of \code{ConsortiumMetabolism} Objects
#'
#' @description Creates a \code{ConsortiumMetabolismSet} combining
#'   multiple \code{ConsortiumMetabolism} objects into a unified
#'   metabolite space. Computes pairwise overlap scores and builds
#'   a dendrogram for clustering.
#'
#' @slot Name character. Display name for the set.
#' @slot Consortia list. List of \code{ConsortiumMetabolism}
#'   objects.
#' @slot Description character. Optional short description.
#' @slot OverlapMatrix matrix. Pairwise dissimilarity matrix
#'   (1 - overlap) between consortia.
#' @slot Dendrogram list. Hierarchical clustering dendrogram.
#' @slot NodeData data.frame. Internal node positions from the
#'   dendrogram.
#' @slot Graphs list. Named list of igraph objects, one per
#'   consortium.
#' @slot BinaryMatrices list. Named list of binary matrices
#'   expanded to universal metabolite space.
#' @slot Edges data.frame. Combined edge list from all consortia
#'   with re-indexed metabolite positions.
#' @slot Metabolites data.frame. Metabolite mapping between
#'   per-consortium and universal indices.
#'
#' @param ... Lists or individual \code{ConsortiumMetabolism}
#'   objects.
#' @param name Character scalar giving the name of the set.
#' @param desc Optional short description of the set.
#'
#' @return A \code{ConsortiumMetabolismSet} object.
#'
#' @seealso \linkS4class{ConsortiumMetabolism},
#'   \link[TreeSummarizedExperiment:TreeSummarizedExperiment-class]{TreeSummarizedExperiment}
#'
#' @examples
#' cm1 <- synCM("comm_1", n_species = 3, max_met = 5)
#' cm2 <- synCM("comm_2", n_species = 4, max_met = 6)
#' cms <- ConsortiumMetabolismSet(cm1, cm2, name = "example")
#' cms
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
    if (!all(vapply(cons, is, logical(1),
                     "ConsortiumMetabolism"))) {
        cli::cli_abort(
            "All elements in {.arg ...} must be
            {.cls ConsortiumMetabolism} objects."
        )
    }
    cli::cli_process_done()

    # ---- Get all mets ----------------------------------------------------------
    cli::cli_status("Getting metabolites")
    # Get indeces of met in each consortium for re-indexing
    all_met <- Map(
        \(x, y) dplyr::mutate(x, consortium = y),
        lapply(cons, \(x) {
            tibble::as_tibble(
                SummarizedExperiment::colData(x)
            )
        }),
        vapply(cons, \(x) x@Name, character(1))
    ) |>
        dplyr::bind_rows() |>
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
        lapply(
            cons,
            \(x) dplyr::mutate(x@Edges, cm_name = x@Name)
        ) |>
        dplyr::bind_rows() |>
        dplyr::left_join(
            dplyr::select(all_met, 1:2),
            by = c(consumed = "met")
        ) |>
        dplyr::rename(c_ind_alig = "met_ind") |>
        dplyr::left_join(
            dplyr::select(all_met, 1:2),
            by = c(produced = "met")
        ) |>
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
    bm_tb <- stats::setNames(
        lapply(cons, \(x) assays(x)$Binary),
        vapply(cons, \(x) x@Name, character(1))
    ) |>
        tibble::enframe() |>
        tibble::rowid_to_column("ind")
    cli::cli_process_done()

    # ---- Pre-expand binary matrices to universal space -------------------------
    cli::cli_status("Expanding binary matrices to universal space")
    universal_mets <- sort(unique(all_met$met))
    expanded_bm <- stats::setNames(
        lapply(
            bm_tb$value,
            \(bm) .expandToUniversalSpace(
                bm, universal_mets
            )
        ),
        bm_tb$name
    )
    cli::cli_process_done()

    # ---- Overlap score ---------------------------------------------------------
    cli::cli_status("Calculating overlap scores")
    ## Build tibble of pre-expanded matrices for pairwise comparison
    exp_tb <- tibble::tibble(
        ind = seq_along(names(expanded_bm)),
        name = names(expanded_bm),
        value = unname(expanded_bm)
    )
    tb <-
        tidyr::expand_grid(x = exp_tb$ind, y = exp_tb$ind) |>
        dplyr::filter(.data$x <= .data$y) |>
        dplyr::left_join(exp_tb, by = c("x" = "ind")) |>
        dplyr::left_join(exp_tb, by = c("y" = "ind")) |>
        dplyr::select(
            "x",
            name_x = "name.x",
            cm_x = "value.x",
            "y",
            name_y = "name.y",
            cm_y = "value.y"
        ) |>
        dplyr::mutate(
            overlap_score = mapply(
                \(x, y) .binMatOverlap(x, y),
                .data$cm_x,
                .data$cm_y
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
    # Make the dendrogram (hclust requires n >= 2)
    if (nrow(overlap_matrix) >= 2L) {
        cli::cli_status("Creating dendrogram")
        dend <-
            stats::dist(overlap_matrix) |>
            stats::hclust() |>
            stats::as.dendrogram()
        cli::cli_process_done()

        # ---- Dendrogram Node Data ------------------------------------------------
        # Get node positions in the dendrogram and index without leaves s.t.
        # they can be used later to retrieve individual clusters.
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
    } else {
        dend <- stats::as.dendrogram(
            stats::hclust(stats::dist(rbind(overlap_matrix, overlap_matrix)))
        )
        dend <- dend[[1L]]
        node_data <- data.frame(
            x = numeric(0),
            y = numeric(0),
            original_node_id = integer(0),
            node_id = integer(0)
        )
    }

    # ---- Graphs ----------------------------------------------------------------
    # Store the graphs in a named list to enable the graph based alignment
    cli::cli_status("Getting all graphs")
    graph_list <- stats::setNames(
        lapply(cons, \(x) x@Graphs[[1]]),
        vapply(cons, \(x) x@Name, character(1))
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
        BinaryMatrices = expanded_bm,
        Edges = all_edges,
        Metabolites = all_met
    )
}

#' @noRd
cms <- ConsortiumMetabolismSet

#' Expand a sparse binary matrix to universal metabolite space
#'
#' Takes a binary matrix from a single CM (with its local metabolite
#' set as row/colnames) and expands it to the full metabolite
#' universe. Missing metabolites get zero rows/columns.
#'
#' @param bm A sparse binary matrix from a CM's Binary assay.
#' @param all_mets Character vector of all metabolites in universal
#'   space, sorted.
#'
#' @return A sparse binary matrix with dimensions
#'   `length(all_mets) x length(all_mets)`.
#'
#' @noRd
.expandToUniversalSpace <- function(bm, all_mets) {
    local_mets <- rownames(bm)
    n <- length(all_mets)
    ## Map local indices to universal indices
    idx <- match(local_mets, all_mets)
    ## Extract triplet form from the local matrix
    bm_t <- Matrix::summary(bm)
    if (nrow(bm_t) == 0L) {
        return(
            Matrix::sparseMatrix(
                i = integer(0L),
                j = integer(0L),
                x = numeric(0L),
                dims = c(n, n),
                dimnames = list(all_mets, all_mets)
            )
        )
    }
    Matrix::sparseMatrix(
        i = idx[bm_t$i],
        j = idx[bm_t$j],
        x = bm_t$x,
        dims = c(n, n),
        dimnames = list(all_mets, all_mets)
    )
}

#' @noRd
.binMatOverlap <- function(bm1, bm2) {
    ## Both matrices must be in the same universal space
    ## (same dimensions and dimnames)
    int <- bm1 * bm2
    denom <- min(sum(bm1), sum(bm2))
    if (denom == 0L) {
        return(0)
    }
    sum(int) / denom
}
