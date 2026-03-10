#' @include AllClasses.R AllGenerics.R
#'
#' @title Constructor for
#'   \code{ConsortiumMetabolismAlignment} objects
#'
#' @description Builds a valid CMA by constructing a TSE base
#'   and setting CMA-specific slots with sensible defaults.
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
