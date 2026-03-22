#' @title Return Species in a Consortium
#'
#' @param object A \code{ConsortiumMetabolism} object
#' @param ... Object specific arguments.
#'
#' @return A character vector containing the names of species in the consortium
#' @export
setGeneric("getSpecies", function(object, ...) standardGeneric("getSpecies"))

#' @title Get Metabolites
#'
#' @description
#' Retrieves the metabolites involved in the metabolic network.
#'
#' @param object A \code{ConsortiumMetabolism} or
#' \code{ConsortiumMetabolismAlignment} object
#' @return A character vector containing the names of metabolites in the network
#' @export
setGeneric("getMet", function(object) standardGeneric("getMet"))

#' @title Get Edges From a \code{ConsortiumMetabolism} Object
#'
#' @description
#' Retrieves the edges representing metabolic interactions between species.
#' The argument `type` can be used to return only specfic types of Edges from
#' a `ConsortiumMetabolismSet` type object.
#' \itemize{
#'   \item `all` will return all edges
#'   \item `pan-cons` will return only edges that exist in all consortia
#'   \item `niche` will return niche edges specific to individual consortia
#'   \item `core` will return core metabolic edges
#'   \item `aux` will return auxiliary edges
#' }
#'
#' @param object A \code{ConsortiumMetabolism} object
#' @param ... Object specific arguments.
#'
#' @return A tibble containing edge information including:
#' \itemize{
#'   \item consumed/produced metabolites
#'   \item number of species involved
#'   \item consumption/production sums
#'   \item effective consumption/production metrics
#' }
#' @export
setGeneric("getEdges", function(object, ...) standardGeneric("getEdges"))

#' @title Get the Consortia
#'
#' @description
#' Returns the consortia data in a tabular format.
#'
#' @param object A \code{ConsortiumMetabolism}, \code{ConsortiumMetabolismSet},
#' or \code{ConsortiumMetabolismAlignment} object
#' @return For \code{ConsortiumMetabolism} objects, returns a tb with species,
#' metabolite and flux information. For \code{ConsortiumMetabolismSet} and
#' \code{ConsortiumMetabolismAlignment} objects, returns a list of such tibbles.
#' @export
setGeneric("getCo", function(object) standardGeneric("getCo"))

#' @title Align Consortium Metabolisms
#'
#' @description
#' Computes functional alignment between consortium metabolisms.
#' Dispatches on the combination of `x` and `y`:
#' \itemize{
#'   \item `align(CM, CM)`: Pairwise alignment of two consortia
#'   \item `align(CMS)`: Multiple alignment across all consortia
#'     in the set
#'   \item `align(CM, CMS)`: Database search -- align one
#'     consortium against all in a set
#' }
#'
#' @param x A [ConsortiumMetabolism] or
#'   [ConsortiumMetabolismSet] object.
#' @param y A [ConsortiumMetabolism],
#'   [ConsortiumMetabolismSet], or `missing`.
#' @param method Character scalar specifying the similarity metric.
#'   One of `"FOS"` (default), `"jaccard"`,
#'   `"brayCurtis"`, `"redundancyOverlap"`, or `"MAAS"`.
#' @param ... Additional arguments passed to methods. Common
#'   arguments include:
#'   \describe{
#'     \item{`computePvalue`}{Logical; compute permutation
#'       p-value? Default `FALSE`.}
#'     \item{`nPermutations`}{Integer; number of permutations
#'       for null model. Default `999L`.}
#'     \item{`BPPARAM`}{A
#'       [BiocParallel::BiocParallelParam] object for
#'       parallel execution.}
#'   }
#'
#' @return A [ConsortiumMetabolismAlignment] object.
#'
#' @export
setGeneric(
    "align",
    function(x, y, method = "FOS", ...) {
        standardGeneric("align")
    }
)

#' @title Modify a \code{ConsortiaMetabolismSet} Object
#'
#' @description
#' Modifies a \code{ConsortiumMetabolismSet} by adding or removing
#' \code{ConsortiumMetabolism} objects. This allows for dynamic updating of
#' consortium sets for comparative analyses.
#'
#' @param object A \code{ConsortiumMetabolismSet} object to modify
#' @return A modified \code{ConsortiumMetabolismSet} object
#' @export
setGeneric("modify", function(object) standardGeneric("modify"))


#' @title Get a cluster from a \code{ConsortiaMetabolismSet} Object
#'
#' @description
#' Gets a cluster from a \code{ConsortiumMetabolismSet} object.
#'
#' @param object A \code{ConsortiumMetabolismSet} object
#' @param node_id Numeric scalar giving the node to be extracted.
#' @param name Character scalar specifying a name for the selection.
#' @param description Character scalar describing the selection.
#'
#' @return A cluster from a \code{ConsortiumMetabolismSet} object
#' @export
setGeneric(
  "getCluster",
  function(object, node_id, name = NA_character_, description = NA_character_) {
    standardGeneric("getCluster")
  }
)


#' @title Set Names or description of ramen Objects
#'
#' @description
#' Sets a new name or description of \code{ConsortiumMetabolism},
#' \code{ConsortiumMetabolismSet}, or \code{ConsortiumMetabolismAlignment}
#' objects
#'
#' @param object A \code{ConsortiumMetabolism},
#' \code{ConsortiumMetabolismSet}, or \code{ConsortiumMetabolismAlignment}
#' object
#' @param value Character scalar specifying a name.
#'
#' @return The modified object with updated name
#' @export
setGeneric(
  "setName",
  function(object, value = NA_character_) standardGeneric("setName")
)

#' @title Set Names or description of ramen Objects
#'
#' @description
#' Sets a new name or description of \code{ConsortiumMetabolism},
#' \code{ConsortiumMetabolismSet}, or \code{ConsortiumMetabolismAlignment}
#' objects
#'
#' @param object A \code{ConsortiumMetabolism},
#' \code{ConsortiumMetabolismSet}, or \code{ConsortiumMetabolismAlignment}
#' object
#' @param value Character scalar specifying a description.
#'
#' @return The modified object with updated description
#' @export
setGeneric(
  "setDesc",
  function(object, value = NA_character_) standardGeneric("setDesc")
)

#' @title Get Functional Groups
#'
#' @description
#' Calculates and returns functional groups based on metabolic reactions.
#' For \code{ConsortiumMetabolismSet} objects, this involves analyzing shared
#' reactions across species to identify clusters of species with similar
#' metabolic capabilities.
#'
#' @details
#' This method is currently implemented for \code{ConsortiumMetabolismSet}
#' objects. Future versions will extend functionality to
#' \code{ConsortiumMetabolismAlignment} objects to allow for comparative
#' functional group analysis across different alignments.
#'
#' @param object A \code{ConsortiumMetabolismSet} object.
#' @param k An integer scalar specifying the number of clusters to color in the
#'   dendrogram.
#' @param label_size Numeric scalar specifying label size in the output plot.
#' @param label_colours If not `NULL`, tb with two columns, label and colour.
#' @param ... Additional arguments to be passed to specific methods.
#'
#' @return A list (returned invisibly) containing:
#' \itemize{
#'   \item plot: The ggplot2 dendrogram visualization
#'   \item dendrogram: The dendrogram object
#'   \item similarity_matrix: Matrix of Jaccard similarities between species
#'   \item species_combinations: Tibble with pairwise species comparisons
#'   \item reactions_per_species: Tibble mapping species to their reactions
#' }
#' The plot is automatically displayed.
#'
#' @export
setGeneric(
  "getFunctionalGroups",
  function(object, k = 4, ...) standardGeneric("getFunctionalGroups")
)


#' @title Compare Species
#'
#' @description
#' Compares the metabolisms of two species and outputs a list of tibbles,
#' containing the tibbles intersection, unique, and consistency.
#' As the name suggests, intersection, unique contain the intersection and
#' unique pathways per species compared, while the consistency tibble contains
#' information on whether or not a specie's set of pathways is consistent in
#' all of the consortia in which it is present or not. If all species contain
#' the same edges in all consortia in which they appear, this tibble will be
#' returned with 0 rows.
#' For \code{ConsortiumMetabolismSet} objects.
#'
#' @details
#' This method is currently implemented for \code{ConsortiumMetabolismSet}
#' objects. Future versions will extend functionality to
#' \code{ConsortiumMetabolism} objects to allow for the analysis of species
#' within a single consortium different alignments.
#'
#' @param object A \code{ConsortiumMetabolismSet} object.
#' @param species A character vector of species names to compare
#'
#' @return A list of tibbles.
#'
#' @export
setGeneric(
  "compareSpecies",
  function(object, species) standardGeneric("compareSpecies")
)

## ---- CMA accessors -------------------------------------------------------

#' @title Get Alignment Scores
#'
#' @description
#' Returns the scores from a
#' [ConsortiumMetabolismAlignment] object.
#'
#' @param object A [ConsortiumMetabolismAlignment] object.
#'
#' @return A named list of scores.
#'
#' @export
setGeneric(
    "scores",
    function(object) standardGeneric("scores")
)

#' @title Get Shared Pathways
#'
#' @description
#' Returns pathways shared between query and reference
#' in a pairwise alignment.
#'
#' @param object A [ConsortiumMetabolismAlignment] object
#'   of type `"pairwise"`.
#'
#' @return A data.frame of shared pathway edges.
#'
#' @export
setGeneric(
    "sharedPathways",
    function(object) standardGeneric("sharedPathways")
)

#' @title Get Unique Pathways
#'
#' @description
#' Returns pathways unique to query and reference in a
#' pairwise alignment.
#'
#' @param object A [ConsortiumMetabolismAlignment] object
#'   of type `"pairwise"`.
#'
#' @return A list with `query` and `reference` data.frames.
#'
#' @export
setGeneric(
    "uniquePathways",
    function(object) standardGeneric("uniquePathways")
)

#' @title Get Similarity Matrix
#'
#' @description
#' Returns the pairwise similarity matrix from a
#' multiple alignment.
#'
#' @param object A [ConsortiumMetabolismAlignment] object
#'   of type `"multiple"`.
#'
#' @return A numeric n x n matrix.
#'
#' @export
setGeneric(
    "similarityMatrix",
    function(object) standardGeneric("similarityMatrix")
)

#' @title Get Edge Prevalence
#'
#' @description
#' Returns edge prevalence across consortia from a
#' multiple alignment.
#'
#' @param object A [ConsortiumMetabolismAlignment] object
#'   of type `"multiple"`.
#'
#' @return A data.frame with columns `consumed`, `produced`,
#'   `nConsortia`, and `proportion`.
#'
#' @export
setGeneric(
    "prevalence",
    function(object) standardGeneric("prevalence")
)

#' @title Get Consensus Edges
#'
#' @description
#' Returns the consensus network edges from a multiple
#' alignment.
#'
#' @param object A [ConsortiumMetabolismAlignment] object
#'   of type `"multiple"`.
#'
#' @return A data.frame of consensus edges with prevalence.
#'
#' @export
setGeneric(
    "consensusEdges",
    function(object) standardGeneric("consensusEdges")
)
