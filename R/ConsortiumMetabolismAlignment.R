#' @include AllClasses.R AllGenerics.R
#'
#' @title Constructor for
#'   \code{ConsortiumMetabolismAlignment} Objects
#'
#' @description Builds a valid
#'   \code{ConsortiumMetabolismAlignment} by constructing a TSE
#'   base and setting CMA-specific slots with sensible defaults.
#'
#' @slot Name character. Display name for the alignment.
#' @slot Description character. Optional short description.
#' @slot Type character. One of \code{"pairwise"},
#'   \code{"multiple"}, or \code{"search"}.
#' @slot Metric character. Similarity metric used (e.g.,
#'   \code{"FOS"}).
#' @slot Params list. Additional parameters passed to the
#'   alignment method.
#' @slot QueryName character. Name of the query consortium
#'   (pairwise).
#' @slot ReferenceName character. Name of the reference
#'   consortium (pairwise).
#' @slot Scores list. Named list of all computed score
#'   components.
#' @slot PrimaryScore numeric. Primary similarity score,
#'   between 0 and 1.
#' @slot Pvalue numeric. Permutation p-value, between 0
#'   and 1.
#' @slot SharedPathways data.frame. Pathways present in both
#'   query and reference.
#' @slot UniqueQuery data.frame. Pathways unique to the query.
#' @slot UniqueReference data.frame. Pathways unique to the
#'   reference.
#' @slot SimilarityMatrix matrix. Pairwise similarity matrix
#'   (multiple alignment).
#' @slot ConsensusEdges data.frame. Consensus network edges
#'   (multiple alignment).
#' @slot Prevalence data.frame. Edge prevalence across
#'   consortia (multiple alignment).
#' @slot Dendrogram list. Hierarchical clustering dendrogram
#'   (multiple alignment).
#' @slot Edges data.frame. Combined edge list.
#' @slot Graphs list. List of igraph objects.
#' @slot Metabolites data.frame. Metabolite mapping table.
#'
#' @param ... Named CMA slot values to set (e.g.,
#'   \code{Type = "pairwise"}, \code{PrimaryScore = 0.8}).
#' @param tse An optional
#'   \code{\link[TreeSummarizedExperiment]{TreeSummarizedExperiment}}
#'   to use as the base. Default creates an empty TSE.
#'
#' @return A validated
#'   \code{\linkS4class{ConsortiumMetabolismAlignment}} object.
#'
#' @seealso \link[TreeSummarizedExperiment:TreeSummarizedExperiment-class]{TreeSummarizedExperiment},
#'   \code{\link{align}}
#'
#' @examples
#' # Empty alignment
#' cma <- ConsortiumMetabolismAlignment()
#'
#' # Pairwise alignment via align()
#' cm1 <- synCM("comm_1", n_species = 3, max_met = 5)
#' cm2 <- synCM("comm_2", n_species = 4, max_met = 6)
#' cma <- align(cm1, cm2)
#' cma
#'
#' @export
#'
#' @importFrom TreeSummarizedExperiment TreeSummarizedExperiment
ConsortiumMetabolismAlignment <- function(..., tse = NULL) {
    if (is.null(tse)) {
        tse <- TreeSummarizedExperiment::TreeSummarizedExperiment()
    }
    args <- list(...)

    ## Defaults for all CMA-specific slots
    defaults <- list(
        Name = NA_character_,
        Description = NA_character_,
        Type = NA_character_,
        Metric = "FOS",
        Params = list(),
        QueryName = NA_character_,
        ReferenceName = NA_character_,
        Scores = list(),
        PrimaryScore = NA_real_,
        Pvalue = NA_real_,
        SharedPathways = data.frame(),
        UniqueQuery = data.frame(),
        UniqueReference = data.frame(),
        SimilarityMatrix = matrix(0, 0, 0),
        ConsensusEdges = data.frame(),
        Prevalence = data.frame(),
        Dendrogram = list(),
        Edges = data.frame(),
        Graphs = list(),
        Metabolites = data.frame()
    )

    ## Merge user values with defaults
    slot_vals <- defaults
    for (nm in names(args)) {
        slot_vals[[nm]] <- args[[nm]]
    }

    ## Build CMA: pass TSE + all slot values to the generator,
    ## matching the pattern used by ConsortiumMetabolism() and
    ## ConsortiumMetabolismSet().
    do.call(
        newConsortiumMetabolismAlignment,
        c(list(tse), slot_vals)
    )
}
