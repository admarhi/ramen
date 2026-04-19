#' @include AllClasses.R AllGenerics.R ConsortiumMetabolismSet.R
NULL

#' @title Subset ramen Objects by Metabolite
#'
#' @description
#' Subsetting methods for ramen S4 classes. For
#' \code{ConsortiumMetabolismSet}, \code{[} follows the
#' \code{TreeSummarizedExperiment} contract: the assays are
#' \emph{metabolite x metabolite} (m x m), so \code{i} and
#' \code{j} index the metabolite space, not the consortium
#' list. After subsetting, all custom slots are updated
#' eagerly: \code{BinaryMatrices} is subsetted to the
#' remaining metabolite space and \code{OverlapMatrix},
#' \code{Dendrogram}, \code{NodeData}, \code{Pathways},
#' \code{Metabolites}, and \code{Graphs} are recomputed
#' accordingly. \code{@Consortia} is not modified — each
#' \code{ConsortiumMetabolism} object retains its own
#' metabolite space.
#'
#' \strong{Note on consortium-level selection:} to extract a
#' subset of consortia rather than a subset of the metabolite
#' space, use \code{\link{extractCluster}} (dendrogram-based)
#' or the planned \code{filterConsortia()} method. Metabolite
#' subsetting is appropriate when the analysis should focus
#' on a specific part of the metabolic network — for example,
#' subsetting to the core-pathway metabolites
#' (\code{pathways(cms, type = "core")}) to compare consortia
#' on their shared metabolic backbone only.
#'
#' \code{ConsortiumMetabolism} and
#' \code{ConsortiumMetabolismAlignment} do not support
#' \code{[} and return an informative error.
#'
#' @param x A ramen S4 object.
#' @param i Row (metabolite) indices — integer, logical, or
#'   character names from \code{metabolites(cms)}.
#' @param j Column indices. Must be identical to \code{i}
#'   (or omitted): since assays are m x m, asymmetric
#'   subsetting is not meaningful.
#' @param ... Additional arguments passed to the TSE parent
#'   method.
#' @param drop Ignored; retained for S4 signature
#'   compatibility.
#'
#' @return For \code{ConsortiumMetabolismSet}: a subsetted
#'   object with all custom slots updated to reflect the
#'   remaining metabolite space. For
#'   \code{ConsortiumMetabolism} and
#'   \code{ConsortiumMetabolismAlignment}: an error.
#'
#' @seealso \code{\link{extractCluster}} for dendrogram-based
#'   consortium selection.
#'
#' @examples
#' cm1 <- synCM("a", n_species = 3, max_met = 6, seed = 1)
#' cm2 <- synCM("b", n_species = 3, max_met = 6, seed = 2)
#' cms <- ConsortiumMetabolismSet(
#'     cm1, cm2, name = "test", verbose = FALSE
#' )
#' ## Subset to first 3 metabolites (m x m assay semantics)
#' sub <- cms[seq_len(3), seq_len(3)]
#' nrow(sub@Metabolites)
#'
#' @name subset-ramen
NULL

## ---- Private helper ---------------------------------------------------------

.rebuildDendrogram <- function(overlap_matrix, linkage = "complete") {
    n <- nrow(overlap_matrix)
    if (n >= 2L) {
        dend <- stats::dist(overlap_matrix) |>
            stats::hclust(method = linkage) |>
            stats::as.dendrogram()
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
            stats::hclust(
                stats::dist(
                    rbind(overlap_matrix, overlap_matrix)
                ),
                method = linkage
            )
        )
        dend <- dend[[1L]]
        node_data <- data.frame(
            x = numeric(0),
            y = numeric(0),
            original_node_id = integer(0),
            node_id = integer(0)
        )
    }
    list(dend = dend, node_data = node_data)
}

## ---- CM [ -------------------------------------------------------------------

#' @describeIn subset-ramen Subsetting a
#'   \code{ConsortiumMetabolism} is not supported.
#' @export
setMethod(
    "[",
    "ConsortiumMetabolism",
    function(x, i, j, ..., drop = FALSE) {
        cli::cli_abort(
            c(
                "{.cls ConsortiumMetabolism} does not support \\
                {.code [} subsetting.",
                "i" = "Subsetting a metabolic network by \\
                individual metabolites is not biologically \\
                meaningful.",
                "i" = "To filter consortia from a set, subset \\
                a {.cls ConsortiumMetabolismSet} instead."
            )
        )
    }
)

## ---- CMS [ ------------------------------------------------------------------

#' @describeIn subset-ramen Subset a
#'   \code{ConsortiumMetabolismSet} by metabolite index with
#'   full slot synchronisation.
#' @importFrom S4Vectors metadata
#' @importFrom igraph induced_subgraph V
#' @importFrom dplyr rename mutate filter arrange row_number desc
#' @importFrom dendextend get_nodes_xy
#' @importFrom stats dist hclust as.dendrogram
#' @importFrom tibble as_tibble
#' @importFrom rlang %||%
#' @export
setMethod(
    "[",
    "ConsortiumMetabolismSet",
    function(x, i, j, ..., drop = FALSE) {
        if (missing(i)) {
            return(x)
        }
        if (!missing(j) && !identical(i, j)) {
            cli::cli_abort(
                c(
                    "{.cls ConsortiumMetabolismSet} assays are \\
                    metabolite x metabolite.",
                    "i" = "{.arg i} and {.arg j} must be \\
                    identical. Supply only {.arg i} to subset \\
                    metabolites symmetrically."
                )
            )
        }

        ## TSE handles assay + rowData/colData subsetting
        x <- callNextMethod(x, i, i, ..., drop = drop)

        remaining_mets <- rownames(x)

        ## Metabolites data.frame
        met_df <- x@Metabolites
        met_df <- met_df[met_df$met %in% remaining_mets, ]
        met_df$met_ind <- seq_len(nrow(met_df))
        x@Metabolites <- met_df

        ## Pathways
        x@Pathways <- x@Pathways[
            x@Pathways$consumed %in% remaining_mets &
                x@Pathways$produced %in% remaining_mets,
        ]

        ## BinaryMatrices
        x@BinaryMatrices <- lapply(
            x@BinaryMatrices,
            function(bm) bm[remaining_mets, remaining_mets]
        )

        ## OverlapMatrix
        x@OverlapMatrix <- .computeFOSMatrix(
            x@BinaryMatrices,
            names(x@BinaryMatrices)
        )

        ## Dendrogram + NodeData
        linkage_used <- metadata(x)[["linkage"]] %||% "complete"
        dr <- .rebuildDendrogram(x@OverlapMatrix, linkage_used)
        x@Dendrogram <- list(dr$dend)
        x@NodeData <- dr$node_data

        ## Graphs
        x@Graphs <- lapply(
            x@Graphs,
            function(g) {
                igraph::induced_subgraph(
                    g,
                    intersect(
                        igraph::V(g)$name,
                        remaining_mets
                    )
                )
            }
        )

        validObject(x)
        x
    }
)

## ---- CMA [ ------------------------------------------------------------------

#' @describeIn subset-ramen Subsetting a
#'   \code{ConsortiumMetabolismAlignment} is not supported.
#' @export
setMethod(
    "[",
    "ConsortiumMetabolismAlignment",
    function(x, i, j, ..., drop = FALSE) {
        cli::cli_abort(
            c(
                "{.cls ConsortiumMetabolismAlignment} does \\
                not support {.code [} subsetting.",
                "i" = "Use {.fun pathways}, {.fun scores}, \\
                or {.fun similarityMatrix} to access results."
            )
        )
    }
)
