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
#' @return A character vector of species names.
#'
#' @note `BiocGenerics::species` (deprecated) sometimes
#'   resolves first when typing `?species` interactively.
#'   Use `?ramen::species` to land on this page.
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

#' @title Species Summary
#'
#' @description
#' Returns an enriched per-species summary as a tibble.
#' Provides more detail than \code{species()}, which
#' returns only a character vector.
#'
#' @param object A \code{ConsortiumMetabolism},
#'   \code{ConsortiumMetabolismSet}, or
#'   \code{ConsortiumMetabolismAlignment} object.
#' @param ... Additional arguments passed to methods.
#'
#' @return A tibble with per-species metrics. Columns
#'   depend on the class:
#'   \itemize{
#'     \item \code{ConsortiumMetabolism}: \code{species},
#'       \code{n_pathways}, \code{n_consumed},
#'       \code{n_produced}.
#'     \item \code{ConsortiumMetabolismSet}:
#'       \code{species}, \code{n_consortia},
#'       \code{n_pathways}.
#'     \item \code{ConsortiumMetabolismAlignment}
#'       (pairwise): \code{species}, \code{role}
#'       (\code{"shared"}, \code{"unique_query"}, or
#'       \code{"unique_reference"}).
#'   }
#'
#' @examples
#' cm <- synCM("test", n_species = 3, max_met = 5)
#' speciesSummary(cm)
#'
#' @export
setGeneric(
    "speciesSummary",
    function(object, ...) standardGeneric("speciesSummary")
)

#' @title Filter Consortia from a Set
#'
#' @description
#' Selects a subset of consortia from a
#' \code{ConsortiumMetabolismSet}, returning a fully
#' recomputed \code{ConsortiumMetabolismSet} containing
#' only the selected consortia.
#'
#' @param object A \code{ConsortiumMetabolismSet} object.
#' @param i Integer vector of indices, character vector
#'   of consortium names, or logical vector of length
#'   equal to the number of consortia.
#'
#' @return A \code{ConsortiumMetabolismSet} containing
#'   only the selected consortia.
#'
#' @examples
#' cm1 <- synCM("a", n_species = 3, max_met = 5)
#' cm2 <- synCM("b", n_species = 3, max_met = 5)
#' cms <- ConsortiumMetabolismSet(cm1, cm2, name = "test")
#' filterConsortia(cms, 1L)
#'
#' @export
setGeneric(
    "filterConsortia",
    function(object, i) standardGeneric("filterConsortia")
)

#' @title Compare Two Species by Pathway Set
#'
#' @description
#' Computes similarity metrics between the pathway sets
#' of two species. Two dispatch modes are supported:
#' \itemize{
#'   \item \code{compareSpecies(cm, sp1, sp2)}: compare
#'     two species within the same
#'     \code{ConsortiumMetabolism} object.
#'   \item \code{compareSpecies(cm1, cm2, sp1, sp2)}:
#'     compare one species from each of two
#'     \code{ConsortiumMetabolism} objects (e.g. the
#'     same species under different growth conditions).
#' }
#'
#' The pathway set for a species is the set of
#' (consumed, produced) pairs in which that species
#' participates.
#'
#' @param x A \code{ConsortiumMetabolism} object, or
#'   the first consortium in a cross-CM comparison.
#' @param y Character scalar naming species 1 (for
#'   same-CM comparison), or a second
#'   \code{ConsortiumMetabolism} object (for cross-CM).
#' @param ... For same-CM: \code{sp2} (character scalar
#'   naming the second species). For cross-CM:
#'   \code{sp1} and \code{sp2} (character scalars naming
#'   one species from each consortium).
#'
#' @return A named list with elements:
#'   \describe{
#'     \item{\code{fos}}{Szymkiewicz-Simpson overlap
#'       score (intersection over min set size).}
#'     \item{\code{jaccard}}{Jaccard similarity
#'       (intersection over union).}
#'     \item{\code{n_shared}}{Number of shared
#'       pathways.}
#'     \item{\code{n_unique_sp1}}{Pathways only in
#'       sp1.}
#'     \item{\code{n_unique_sp2}}{Pathways only in
#'       sp2.}
#'   }
#'
#' @examples
#' cm <- synCM("test", n_species = 3, max_met = 5)
#' sp <- species(cm)
#' compareSpecies(cm, sp[1], sp[2])
#'
#' @export
setGeneric(
    "compareSpecies",
    function(x, y, ...) standardGeneric("compareSpecies")
)

#' @title Get Metabolites
#'
#' @description
#' Retrieves the metabolites involved in the metabolic
#' network. For \code{ConsortiumMetabolism} objects, the
#' result can optionally be restricted to a specific
#' species and/or direction (\code{"consumed"} or
#' \code{"produced"}).
#'
#' @param object A \code{ConsortiumMetabolism},
#'   \code{ConsortiumMetabolismSet}, or
#'   \code{ConsortiumMetabolismAlignment} object.
#' @param ... Additional arguments. For
#'   \code{ConsortiumMetabolism}: \code{species}
#'   (character scalar; restrict to metabolites involved
#'   with this species) and \code{direction} (one of
#'   \code{"all"}, \code{"consumed"}, or
#'   \code{"produced"}; defaults to \code{"all"}).
#'
#' @return A character vector containing the names of
#'   metabolites in the network.
#'
#' @examples
#' cm <- synCM("test", n_species = 3, max_met = 5)
#' metabolites(cm)
#' ## Metabolites consumed by a specific species:
#' sp <- species(cm)[1]
#' metabolites(cm, species = sp, direction = "consumed")
#'
#' @export
setGeneric(
    "metabolites",
    function(object, ...) standardGeneric("metabolites")
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

#' @title Get the Constituent Consortia
#'
#' @description
#' Returns the \code{ConsortiumMetabolism} objects that a
#' container is built from. Defined for
#' \code{ConsortiumMetabolismSet} (the consortia in the set)
#' and \code{ConsortiumMetabolismAlignment} (the consortia
#' that produced the alignment). The plural noun reflects
#' that the result is always a collection; for a single
#' \code{ConsortiumMetabolism}, use
#' \code{\link[=as.data.frame]{as.data.frame()}} to obtain
#' the underlying edge list.
#'
#' @param object A \code{ConsortiumMetabolismSet} or
#'   \code{ConsortiumMetabolismAlignment} object.
#'
#' @return A named list of \code{ConsortiumMetabolism}
#'   objects.
#'
#' @examples
#' cm1 <- synCM("a", n_species = 3, max_met = 5)
#' cm2 <- synCM("b", n_species = 3, max_met = 5)
#' cms <- ConsortiumMetabolismSet(list(cm1, cm2), name = "demo")
#' consortia(cms)
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
#'     \item{`linkage`}{Character; agglomeration method for
#'       hierarchical clustering (multiple alignment only).
#'       One of `"complete"` (default), `"average"`,
#'       `"single"`, or `"ward.D2"`. Passed to
#'       \code{\link[stats]{hclust}}.}
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
#' @note `method = "brayCurtis"` is defined only for
#'   weighted CMs (i.e. constructed from non-unit fluxes).
#'   On unweighted inputs the score returns `NA`. To
#'   suppress this, build CMs from a weighted edge list or
#'   pick a different metric.
#'
#' @note `align(CM, CMS)` (database search) currently
#'   raises a not-yet-implemented error; the dispatch is
#'   reserved for the upcoming MinHash-prefiltered search
#'   feature.
#'
#' @note If `tibble` (or a package that re-exports
#'   `pillar::align`) is on the search path, plain
#'   `?align` may resolve to `pillar::align`. Use
#'   `?ramen::align` to land on this page.
#'
#' @note For the formal definitions of FOS (Szymkiewicz-Simpson),
#'   Jaccard, Bray-Curtis, RedundancyOverlap, MAAS, coverage
#'   ratios, the Hill-1 perplexity construction underlying the
#'   EffectiveConsumption / EffectiveProduction (flux-corrected)
#'   and nEffectiveSpeciesConsumption / nEffectiveSpeciesProduction
#'   (effective species counts) assays, and the
#'   "Mathematical formulation" section of
#'   \code{vignette("alignment", package = "ramen")}. The same
#'   section discusses cross-product inflation, Hill-1 saturation
#'   on small consortia, flux reversibility, and alternate-optima
#'   caveats relevant to interpreting the scores returned here.
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
#' cm1 <- synCM("comm_1", n_species = 3, max_met = 5)
#' cm2 <- synCM("comm_2", n_species = 4, max_met = 6)
#' cms <- ConsortiumMetabolismSet(
#'     cm1, cm2, name = "test"
#' )
#' extractCluster(cms, node_id = 1)
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
#' metabolic pathways. For
#' \code{ConsortiumMetabolism} objects, species are
#' clustered within a single consortium. For
#' \code{ConsortiumMetabolismSet} objects, the analysis
#' pools species across all consortia in the set,
#' identifying clusters of species with similar metabolic
#' capabilities regardless of which consortium they
#' belong to.
#'
#' @details
#' This method computes a Jaccard similarity matrix
#' between species based on shared pathways, then
#' performs hierarchical clustering. A pathway is
#' represented as the unique \code{(consumed, produced)}
#' metabolite pair. To visualise the resulting
#' dendrogram, pass the output to
#' \code{\link{plotFunctionalGroups}}.
#'
#' If a \code{ConsortiumMetabolism} contains fewer than
#' two species, a warning is emitted and the returned
#' list has \code{dendrogram = NULL}; the incidence
#' matrix and (trivial) similarity matrix are still
#' returned so downstream code can inspect them.
#'
#' @param object A \code{ConsortiumMetabolism} or
#'   \code{ConsortiumMetabolismSet} object.
#' @param ... Additional arguments passed to methods.
#'   Supported arguments include:
#'   \describe{
#'     \item{\code{linkage}}{Character scalar specifying
#'       the agglomeration method for hierarchical
#'       clustering. Passed to
#'       \code{\link[stats]{hclust}} as the
#'       \code{method} argument. One of
#'       \code{"complete"} (default),
#'       \code{"average"}, \code{"single"}, or
#'       \code{"ward.D2"}.}
#'   }
#'
#' @return A list (returned invisibly) containing:
#' \itemize{
#'   \item \code{dendrogram}: The dendrogram object, or
#'     \code{NULL} when fewer than two species are
#'     available.
#'   \item \code{similarity_matrix}: Matrix of Jaccard
#'     similarities between species.
#'   \item \code{incidence_matrix}: Sparse binary
#'     species-by-pathway incidence matrix.
#'   \item \code{reactions_per_species}: Data frame
#'     mapping species to their pathways.
#' }
#'
#' @examples
#' ## Single consortium
#' cm <- synCM("test", n_species = 4, max_met = 8)
#' fg_cm <- functionalGroups(cm)
#' plotFunctionalGroups(fg_cm, k = 2)
#'
#' ## Set of consortia
#' cm1 <- synCM("comm_1", n_species = 3, max_met = 5)
#' cm2 <- synCM("comm_2", n_species = 4, max_met = 6)
#' cms <- ConsortiumMetabolismSet(
#'     cm1, cm2, name = "test"
#' )
#' fg <- functionalGroups(cms)
#' plotFunctionalGroups(fg, k = 2)
#'
#' @seealso \code{\link{plotFunctionalGroups}} for
#'   visualising the dendrogram.
#'
#' @export
setGeneric(
    "functionalGroups",
    function(object, ...) {
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
#' cm1 <- synCM("comm_1", n_species = 3, max_met = 5)
#' cm2 <- synCM("comm_2", n_species = 4, max_met = 6)
#' cms <- ConsortiumMetabolismSet(cm1, cm2, name = "test")
#' cma <- align(cms)
#' similarityMatrix(cma)
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
#' cm1 <- synCM("comm_1", n_species = 3, max_met = 5)
#' cm2 <- synCM("comm_2", n_species = 4, max_met = 6)
#' cms <- ConsortiumMetabolismSet(cm1, cm2, name = "test")
#' cma <- align(cms)
#' prevalence(cma)
#'
#' @export
setGeneric(
    "prevalence",
    function(object) standardGeneric("prevalence")
)


## ---- CMS accessors ---------------------------------------------------------

#' @title Get Overlap Matrix
#'
#' @description
#' Returns the pairwise dissimilarity matrix from a
#' \code{ConsortiumMetabolismSet} object. Values are
#' \code{1 - FOS} (Functional Overlap Score), so 0
#' indicates identical consortia and 1 indicates no
#' shared pathways.
#'
#' @param object A \code{ConsortiumMetabolismSet} object.
#'
#' @return A numeric \eqn{n \times n} matrix of pairwise
#'   dissimilarities, where \eqn{n} is the number of
#'   consortia. Row and column names are consortium
#'   names.
#'
#' @examples
#' cm1 <- synCM("comm_1", n_species = 3, max_met = 5)
#' cm2 <- synCM("comm_2", n_species = 4, max_met = 6)
#' cms <- ConsortiumMetabolismSet(
#'     cm1, cm2, name = "test"
#' )
#' overlapMatrix(cms)
#'
#' @export
setGeneric(
    "overlapMatrix",
    function(object) standardGeneric("overlapMatrix")
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

#' @title Get Growth Rates
#'
#' @description
#' Returns per-species growth rates (e.g. FBA objective
#' values) stored in a ramen object. For
#' [ConsortiumMetabolism] objects, returns the named
#' numeric vector supplied at construction. For
#' [ConsortiumMetabolismSet] objects, returns a named
#' list with one entry per consortium.
#'
#' @param object A [ConsortiumMetabolism] or
#'   [ConsortiumMetabolismSet] object.
#'
#' @return For [ConsortiumMetabolism]: a named numeric
#'   vector of growth rates (names = species), or
#'   \code{NULL} if no growth data was supplied. For
#'   [ConsortiumMetabolismSet]: a named list of such
#'   vectors.
#'
#' @examples
#' data <- data.frame(
#'     species = c("s1", "s1", "s2", "s2"),
#'     metabolite = c("m1", "m2", "m1", "m3"),
#'     flux = c(-1, 1, -1, 1)
#' )
#' cm <- ConsortiumMetabolism(
#'     data, name = "test",
#'     growth = c(s1 = 0.5, s2 = 0.3)
#' )
#' growth(cm)
#'
#' @export
setGeneric(
    "growth",
    function(object) standardGeneric("growth")
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
