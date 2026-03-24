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
#' @include helpers-align.R
#' @export
ConsortiumMetabolismSet <- function(
    ...,
    name = NA_character_,
    desc = NA_character_
) {
    t_start <- proc.time()[["elapsed"]]
    cli::cli_h1("Creating CMS {.val {name}}")

    ## ---- 1. Validate arguments --------------------------------
    args <- list(...)
    cons <- unlist(args, recursive = FALSE, use.names = FALSE)
    n_cons <- length(cons)
    cli::cli_progress_step(
        "Validating {.val {n_cons}} \\
        {.cls ConsortiumMetabolism} object{?s}"
    )
    if (!all(vapply(
        cons, is, logical(1L), "ConsortiumMetabolism"
    ))) {
        cli::cli_abort(
            "All elements in {.arg ...} must be
            {.cls ConsortiumMetabolism} objects."
        )
    }
    cm_names <- vapply(
        cons, \(x) x@Name, character(1L)
    )

    ## ---- 2. Collect metabolites -------------------------------
    cli::cli_progress_step(
        "Collecting metabolites from \\
        {.val {n_cons}} consortia"
    )
    all_met <- Map(
        \(x, y) dplyr::mutate(x, consortium = y),
        lapply(cons, \(x) {
            tibble::as_tibble(x@colData)
        }),
        cm_names
    ) |>
        dplyr::bind_rows() |>
        dplyr::rename(consortium_ind = "index")

    ## ---- 3. Re-index metabolites ------------------------------
    new_met_ind <- tibble::tibble(
        met = unique(all_met$met)
    ) |>
        dplyr::arrange(.data$met, .locale = "C") |>
        tibble::rowid_to_column("met_ind")
    universal_mets <- new_met_ind$met
    n_mets <- length(universal_mets)

    cli::cli_progress_step(
        "Re-indexing {.val {n_mets}} unique metabolites"
    )
    all_met <- all_met |>
        dplyr::left_join(new_met_ind, by = "met") |>
        dplyr::relocate("met_ind", "met", "consortium") |>
        tidyr::pivot_wider(
            names_from = "consortium",
            values_from = "consortium_ind"
        ) |>
        dplyr::arrange(.data$met_ind)

    ## ---- 4. Expand binary matrices to universal space ---------
    cli::cli_progress_step(
        "Expanding {.val {n_cons}} binary matrices \\
        to {.val {n_mets}}-dimensional space"
    )
    expanded_bm <- stats::setNames(
        lapply(seq_len(n_cons), \(i) {
            .expandMatrix(
                cons[[i]]@assays@data$Binary,
                universal_mets
            )
        }),
        cm_names
    )

    ## ---- 5. Levels matrix from binary matrices ----------------
    cli::cli_progress_step(
        "Computing {.val {n_mets}} x {.val {n_mets}} \\
        levels matrix"
    )
    levels_mat <- as.matrix(Reduce(`+`, expanded_bm))

    ## ---- 6. Pairwise overlap via crossprod --------------------
    n_pairs <- n_cons * (n_cons - 1L) / 2L
    cli::cli_progress_step(
        "Computing pairwise overlap \\
        ({.val {n_pairs}} pairs via crossprod)"
    )
    overlap_matrix <- .computeFOSMatrix(
        expanded_bm, cm_names
    )

    ## ---- 7. Assemble edges ------------------------------------
    cli::cli_progress_step(
        "Assembling edge data from \\
        {.val {n_cons}} consortia"
    )
    all_edges <-
        lapply(seq_len(n_cons), \(i) {
            dplyr::mutate(
                cons[[i]]@Edges, cm_name = cm_names[[i]]
            )
        }) |>
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
        dplyr::select(-"c_prob", -"p_prob")

    ## ---- 8. Dendrogram ----------------------------------------
    if (n_cons >= 2L) {
        cli::cli_progress_step(
            "Building dendrogram from \\
            {.val {n_cons}} x {.val {n_cons}} \\
            dissimilarity matrix"
        )
        dend <-
            stats::dist(overlap_matrix) |>
            stats::hclust() |>
            stats::as.dendrogram()

        cli::cli_progress_step(
            "Extracting dendrogram node positions"
        )
        node_data <- dendextend::get_nodes_xy(dend) |>
            as.data.frame() |>
            tibble::as_tibble() |>
            dplyr::rename(x = "V1", y = "V2") |>
            dplyr::mutate(
                original_node_id = dplyr::row_number()
            ) |>
            dplyr::filter(.data$y != 0) |>
            dplyr::arrange(dplyr::desc(.data$y)) |>
            dplyr::mutate(node_id = dplyr::row_number())
    } else {
        dend <- stats::as.dendrogram(
            stats::hclust(stats::dist(
                rbind(overlap_matrix, overlap_matrix)
            ))
        )
        dend <- dend[[1L]]
        node_data <- data.frame(
            x = numeric(0),
            y = numeric(0),
            original_node_id = integer(0),
            node_id = integer(0)
        )
    }

    ## ---- 9. Graphs --------------------------------------------
    cli::cli_progress_step(
        "Collecting {.val {n_cons}} consortium graphs"
    )
    graph_list <- stats::setNames(
        lapply(cons, \(x) x@Graphs[[1L]]),
        cm_names
    )

    ## ---- Done -------------------------------------------------
    elapsed <- round(
        proc.time()[["elapsed"]] - t_start, 1L
    )
    cli::cli_alert_success(
        "CMS {.val {name}} created: \\
        {.val {n_cons}} consortia, \\
        {.val {n_mets}} metabolites ({elapsed}s)"
    )

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

#' Compute pairwise FOS dissimilarity matrix via crossprod
#'
#' Flattens each m x m expanded binary matrix into a length-m^2
#' sparse column vector, stacks them into one m^2 x n matrix, and
#' uses a single \code{Matrix::crossprod()} to compute all pairwise
#' intersections at once.
#'
#' @param expanded_bm Named list of expanded sparse binary matrices
#'   (all same dimensions).
#' @param cm_names Character vector of consortium names.
#'
#' @return A dense n x n symmetric dissimilarity matrix
#'   (1 - FOS) with \code{cm_names} as dimnames.
#'
#' @noRd
.computeFOSMatrix <- function(expanded_bm, cm_names) {
    n <- length(expanded_bm)
    m <- nrow(expanded_bm[[1L]])

    ## Flatten each m x m sparse matrix to a length-m^2
    ## sparse column, then cbind into one m^2 x n matrix
    flat_list <- lapply(expanded_bm, function(bm) {
        trip <- Matrix::summary(bm)
        if (nrow(trip) == 0L) {
            return(Matrix::sparseMatrix(
                i = integer(0L),
                j = integer(0L),
                x = numeric(0L),
                dims = c(m * m, 1L)
            ))
        }
        ## Column-major flat index: (j-1)*m + i
        flat_idx <- (trip$j - 1L) * m + trip$i
        Matrix::sparseMatrix(
            i = flat_idx,
            j = rep(1L, length(flat_idx)),
            x = trip$x,
            dims = c(m * m, 1L)
        )
    })
    big_mat <- do.call(cbind, flat_list)

    ## Single crossprod = all n*(n-1)/2 intersections
    intersection <- as.matrix(
        Matrix::crossprod(big_mat)
    )

    ## Pre-computed column sums = sum of each binary matrix
    col_sums <- Matrix::colSums(big_mat)

    ## FOS = intersection / min(sum_i, sum_j)
    denom <- outer(col_sums, col_sums, pmin)
    fos_mat <- matrix(
        0, n, n,
        dimnames = list(cm_names, cm_names)
    )
    nonzero <- denom > 0
    fos_mat[nonzero] <- intersection[nonzero] /
        denom[nonzero]

    ## Return full symmetric dissimilarity matrix
    1 - fos_mat
}


