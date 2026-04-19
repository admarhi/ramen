#' @include AllClasses.R AllGenerics.R
#'
#' @title Constructor for
#'   \code{ConsortiumMetabolismAlignment} Objects
#'
#' @description Builds a valid
#'   \code{ConsortiumMetabolismAlignment} by setting CMA-specific
#'   slots with sensible defaults.
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
#' @slot ConsensusPathways data.frame. Consensus network
#'   pathways (multiple alignment).
#' @slot Prevalence data.frame. Pathway prevalence across
#'   consortia (multiple alignment).
#' @slot Dendrogram list. Hierarchical clustering dendrogram
#'   (multiple alignment).
#' @slot Pathways data.frame. Combined pathway list.
#' @slot Graphs list. List of igraph objects.
#' @slot Metabolites data.frame. Metabolite mapping table.
#'
#' @param ... Named CMA slot values to set (e.g.,
#'   \code{Type = "pairwise"}, \code{PrimaryScore = 0.8}).
#'
#' @return A validated
#'   \code{\linkS4class{ConsortiumMetabolismAlignment}} object.
#'
#' @seealso \code{\link{align}}
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
ConsortiumMetabolismAlignment <- function(...) {
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
        ConsensusPathways = data.frame(),
        Prevalence = data.frame(),
        Dendrogram = list(),
        Pathways = data.frame(),
        Graphs = list(),
        Metabolites = data.frame()
    )

    ## Merge user values with defaults
    slot_vals <- defaults
    for (nm in names(args)) {
        slot_vals[[nm]] <- args[[nm]]
    }

    do.call(newConsortiumMetabolismAlignment, slot_vals)
}

#' @rdname ConsortiumMetabolismAlignment
#'
#' @title Coerce a ConsortiumMetabolismAlignment to a data.frame
#'
#' @description Exports alignment results to a plain
#'   \code{data.frame} suitable for downstream analysis.
#'   For \code{"pairwise"} alignments the three pathway sets
#'   (\code{SharedPathways}, \code{UniqueQuery},
#'   \code{UniqueReference}) are row-bound with a
#'   \code{pathway_type} column added.
#'   For \code{"multiple"} alignments \code{ConsensusPathways}
#'   is returned. For all other types (or empty objects)
#'   \code{Pathways} is returned, falling back to an empty
#'   \code{data.frame()}.
#'
#' @param x A
#'   \code{\linkS4class{ConsortiumMetabolismAlignment}} object.
#' @param ... Additional arguments (currently unused).
#'
#' @return A \code{data.frame}. For pairwise alignments the
#'   result contains at minimum \code{consumed}, \code{produced},
#'   and \code{pathway_type} columns.
#'
#' @examples
#' cm1 <- synCM("comm_1", n_species = 3, max_met = 5)
#' cm2 <- synCM("comm_2", n_species = 4, max_met = 6)
#' cma <- align(cm1, cm2)
#' df <- as.data.frame(cma)
#' head(df)
#'
#' @export
setMethod(
    "as.data.frame",
    "ConsortiumMetabolismAlignment",
    function(x, ...) {
        type <- x@Type
        if (isTRUE(type == "pairwise")) {
            shared <- x@SharedPathways
            uq <- x@UniqueQuery
            ur <- x@UniqueReference
            if (nrow(shared) > 0L)
                shared$pathway_type <- "shared"
            if (nrow(uq) > 0L)
                uq$pathway_type <- "unique_query"
            if (nrow(ur) > 0L)
                ur$pathway_type <- "unique_reference"
            return(as.data.frame(dplyr::bind_rows(shared, uq, ur)))
        }
        if (isTRUE(type == "multiple")) {
            cp <- x@ConsensusPathways
            if (nrow(cp) > 0L)
                return(as.data.frame(cp))
            return(as.data.frame(x@Pathways))
        }
        pways <- x@Pathways
        if (nrow(pways) > 0L)
            return(as.data.frame(pways))
        data.frame()
    }
)
