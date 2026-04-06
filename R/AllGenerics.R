#' @title Return Species in a Consortium
#'
#' @description
#' Returns the species present in a consortium or set of
#' consortia.
#'
#' @param object A [ConsortiumMetabolism],
#'   [ConsortiumMetabolismSet], or
#'   [ConsortiumMetabolismAlignment] object.
#' @param ... Additional arguments passed to methods.
#'
#' @return A character vector of species names, or a
#'   tibble with species and pathway counts (for
#'   [ConsortiumMetabolismSet]).
#'
#' @examples
#' cm <- synCM("test", n_species = 3, max_met = 5)
#' species(cm)
#'
#' @export
setGeneric(
    "species",
    function(object, ...) standardGeneric("species")
)

#' @title Get Metabolites
#'
#' @description
#' Retrieves the metabolites involved in the metabolic
#' network.
#'
#' @param object A \code{ConsortiumMetabolism} or
#'   \code{ConsortiumMetabolismAlignment} object.
#'
#' @return A character vector containing the names of
#'   metabolites in the network.
#'
#' @examples
#' cm <- synCM("test", n_species = 3, max_met = 5)
#' metabolites(cm)
#'
#' @export
setGeneric(
    "metabolites",
    function(object) standardGeneric("metabolites")
)

#' @title Retrieve Metabolic Pathways
#'
#' @description
#' Retrieves the pathways representing metabolic
#' interactions between species.
#'
#' By default, returns a concise summary with columns
#' \code{consumed}, \code{produced}, and
#' \code{n_species}. For
#' \code{ConsortiumMetabolismSet} objects, \code{n_cons}
#' is included as well. Set \code{verbose = TRUE} to
#' return the full pathway data including flux
#' statistics, indices, and per-species detail.
#'
#' The argument \code{type} can be used to return only
#' specific types of pathways from a
#' \code{ConsortiumMetabolismSet} object:
#' \itemize{
#'   \item \code{"all"} returns all pathways
#'   \item \code{"pan-cons"} returns pathways present in
#'     most consortia
#'   \item \code{"niche"} returns niche pathways specific
#'     to few consortia
#'   \item \code{"core"} returns core metabolic pathways
#'     shared across most species
#'   \item \code{"aux"} returns auxiliary pathways found
#'     in few species
#' }
#'
#' For \code{ConsortiumMetabolismAlignment} objects,
#' \code{type} selects the pathway subset:
#' \itemize{
#'   \item \code{"all"} returns the union of all pathways
#'   \item \code{"shared"} (pairwise only) returns
#'     pathways shared between query and reference
#'   \item \code{"unique"} (pairwise only) returns
#'     pathways unique to query and reference as a list
#'   \item \code{"consensus"} (multiple only) returns
#'     consensus network pathways with prevalence
#' }
#'
#' @param object A \code{ConsortiumMetabolism},
#'   \code{ConsortiumMetabolismSet}, or
#'   \code{ConsortiumMetabolismAlignment} object.
#' @param ... Object specific arguments. See methods for
#'   details.
#'
#' @return A data.frame of pathway information. With
#'   \code{verbose = FALSE} (default): \code{consumed},
#'   \code{produced}, \code{n_species} (and
#'   \code{n_cons} for CMS objects). With
#'   \code{verbose = TRUE}: all available columns
#'   including flux statistics and indices. For CMA
#'   objects, the return depends on \code{type}: a
#'   data.frame for \code{"all"}, \code{"shared"}, and
#'   \code{"consensus"}; a list of data.frames for
#'   \code{"unique"}.
#'
#' @examples
#' cm <- synCM("test", n_species = 3, max_met = 5)
#' pathways(cm)
#' pathways(cm, verbose = TRUE)
#'
#' cm1 <- synCM("comm_1", n_species = 3, max_met = 5)
#' cm2 <- synCM("comm_2", n_species = 4, max_met = 6)
#' cma <- align(cm1, cm2)
#' pathways(cma)
#' pathways(cma, type = "shared")
#' pathways(cma, type = "unique")
#'
#' @export
setGeneric(
    "pathways",
    function(object, ...) standardGeneric("pathways")
)

#' @title Get the Consortia
#'
#' @description
#' Returns the consortia data in a tabular format.
#'
#' @param object A \code{ConsortiumMetabolism},
#'   \code{ConsortiumMetabolismSet}, or
#'   \code{ConsortiumMetabolismAlignment} object.
#'
#' @return For \code{ConsortiumMetabolism} objects, returns
#'   a tibble with species, metabolite and flux
#'   information. For \code{ConsortiumMetabolismSet} and
#'   \code{ConsortiumMetabolismAlignment} objects, returns
#'   a list of such tibbles.
#'
#' @examples
#' cm <- synCM("test", n_species = 3, max_met = 5)
#' consortia(cm)
#'
#' @export
setGeneric(
    "consortia",
    function(object) standardGeneric("consortia")
)

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
#' @examples
#' cm1 <- synCM("comm_1", n_species = 3, max_met = 5)
#' cm2 <- synCM("comm_2", n_species = 4, max_met = 6)
#' cma <- align(cm1, cm2)
#' cma
#'
#' @export
setGeneric(
    "align",
    function(x, y, method = "FOS", ...) {
        standardGeneric("align")
    }
)

#' @title Extract a Cluster
#'
#' @description
#' Extracts a cluster from a
#' \code{ConsortiumMetabolismSet} object.
#'
#' @param object A \code{ConsortiumMetabolismSet} object.
#' @param node_id Numeric scalar giving the node to be
#'   extracted.
#' @param name Character scalar specifying a name for the
#'   selection.
#' @param description Character scalar describing the
#'   selection.
#'
#' @return A \code{ConsortiumMetabolismSet} object
#'   containing the extracted cluster.
#'
#' @examples
#' \donttest{
#' cm1 <- synCM("comm_1", n_species = 3, max_met = 5)
#' cm2 <- synCM("comm_2", n_species = 4, max_met = 6)
#' cms <- ConsortiumMetabolismSet(
#'     cm1, cm2, name = "test"
#' )
#' extractCluster(cms, node_id = 1)
#' }
#'
#' @export
setGeneric(
    "extractCluster",
    function(
        object,
        node_id,
        name = NA_character_,
        description = NA_character_
    ) {
        standardGeneric("extractCluster")
    }
)


#' @title Get Functional Groups
#'
#' @description
#' Calculates and returns functional groups based on
#' metabolic reactions. For
#' \code{ConsortiumMetabolismSet} objects, this involves
#' analyzing shared reactions across species to identify
#' clusters of species with similar metabolic
#' capabilities.
#'
#' @details
#' This method is currently implemented for
#' \code{ConsortiumMetabolismSet} objects. Future versions
#' will extend functionality to
#' \code{ConsortiumMetabolismAlignment} objects to allow
#' for comparative functional group analysis across
#' different alignments.
#'
#' @param object A \code{ConsortiumMetabolismSet} object.
#' @param k An integer scalar specifying the number of
#'   clusters to color in the dendrogram.
#' @param label_size Numeric scalar specifying label size
#'   in the output plot.
#' @param label_colours If not \code{NULL}, a tibble with
#'   columns \code{label} and \code{colour}.
#' @param ... Additional arguments to be passed to
#'   specific methods.
#'
#' @return A list (returned invisibly) containing:
#' \itemize{
#'   \item plot: The ggplot2 dendrogram visualization
#'   \item dendrogram: The dendrogram object
#'   \item similarity_matrix: Matrix of Jaccard
#'     similarities between species
#'   \item species_combinations: Tibble with pairwise
#'     species comparisons
#'   \item reactions_per_species: Tibble mapping species
#'     to their reactions
#' }
#' The plot is automatically displayed.
#'
#' @examples
#' \donttest{
#' cm1 <- synCM("comm_1", n_species = 3, max_met = 5)
#' cm2 <- synCM("comm_2", n_species = 4, max_met = 6)
#' cms <- ConsortiumMetabolismSet(
#'     cm1, cm2, name = "test"
#' )
#' functionalGroups(cms, k = 2)
#' }
#'
#' @export
setGeneric(
    "functionalGroups",
    function(object, k = 4, ...) {
        standardGeneric("functionalGroups")
    }
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
#' @examples
#' cm1 <- synCM("comm_1", n_species = 3, max_met = 5)
#' cm2 <- synCM("comm_2", n_species = 4, max_met = 6)
#' cma <- align(cm1, cm2)
#' scores(cma)
#'
#' @export
setGeneric(
    "scores",
    function(object) standardGeneric("scores")
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
#' @examples
#' \donttest{
#' cm1 <- synCM("comm_1", n_species = 3, max_met = 5)
#' cm2 <- synCM("comm_2", n_species = 4, max_met = 6)
#' cms <- ConsortiumMetabolismSet(cm1, cm2, name = "test")
#' cma <- align(cms)
#' similarityMatrix(cma)
#' }
#'
#' @export
setGeneric(
    "similarityMatrix",
    function(object) standardGeneric("similarityMatrix")
)

#' @title Get Pathway Prevalence
#'
#' @description
#' Returns pathway prevalence across consortia from a
#' multiple alignment.
#'
#' @param object A [ConsortiumMetabolismAlignment] object
#'   of type `"multiple"`.
#'
#' @return A data.frame with columns `consumed`, `produced`,
#'   `nConsortia`, and `proportion`.
#'
#' @examples
#' \donttest{
#' cm1 <- synCM("comm_1", n_species = 3, max_met = 5)
#' cm2 <- synCM("comm_2", n_species = 4, max_met = 6)
#' cms <- ConsortiumMetabolismSet(cm1, cm2, name = "test")
#' cma <- align(cms)
#' prevalence(cma)
#' }
#'
#' @export
setGeneric(
    "prevalence",
    function(object) standardGeneric("prevalence")
)


## ---- Read accessors --------------------------------------------------------

#' @title Get or Set Object Name
#'
#' @description
#' Returns or sets the name of a ramen object.
#'
#' @param object A [ConsortiumMetabolism],
#'   [ConsortiumMetabolismSet], or
#'   [ConsortiumMetabolismAlignment] object.
#' @param value Character scalar specifying the new name.
#'
#' @return A character scalar (getter), or the modified
#'   object (setter).
#'
#' @examples
#' cm <- synCM("test", n_species = 3, max_met = 5)
#' name(cm)
#'
#' @export
setGeneric("name", function(object) standardGeneric("name"))

#' @title Get or Set Object Description
#'
#' @description
#' Returns or sets the description of a ramen object.
#'
#' @param object A [ConsortiumMetabolism],
#'   [ConsortiumMetabolismSet], or
#'   [ConsortiumMetabolismAlignment] object.
#' @param value Character scalar specifying the new
#'   description.
#'
#' @return A character scalar (getter), or the modified
#'   object (setter).
#'
#' @examples
#' cm <- synCM("test", n_species = 3, max_met = 5)
#' description(cm)
#'
#' @export
setGeneric(
    "description",
    function(object) standardGeneric("description")
)

## ---- Replacement methods ---------------------------------------------------

#' @rdname name
#' @export
setGeneric(
    "name<-",
    function(object, value) standardGeneric("name<-")
)

#' @rdname description
#' @export
setGeneric(
    "description<-",
    function(object, value) standardGeneric("description<-")
)
