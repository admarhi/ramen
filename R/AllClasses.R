#' @rdname ConsortiumMetabolism
#' @exportClass ConsortiumMetabolism
#' @import methods
#' @importClassesFrom TreeSummarizedExperiment TreeSummarizedExperiment
newConsortiumMetabolism <- setClass(
    Class = "ConsortiumMetabolism",
    contains = "TreeSummarizedExperiment",
    slots = list(
        Name = "character",
        Description = "character",
        Pathways = "data.frame",
        Weighted = "logical",
        InputData = "data.frame",
        Metabolites = "character",
        Graphs = "list"
    )
)

#' @rdname ConsortiumMetabolismSet
#' @exportClass ConsortiumMetabolismSet
#' @import methods
#' @importClassesFrom TreeSummarizedExperiment TreeSummarizedExperiment
newConsortiumMetabolismSet <- setClass(
    Class = "ConsortiumMetabolismSet",
    contains = "TreeSummarizedExperiment",
    slots = list(
        Name = "character",
        Consortia = "list",
        Description = "character",
        OverlapMatrix = "matrix",
        Dendrogram = "list",
        NodeData = "data.frame",
        Graphs = "list",
        BinaryMatrices = "list",
        Pathways = "data.frame",
        Metabolites = "data.frame"
    )
)

#' @rdname ConsortiumMetabolismAlignment
#' @exportClass ConsortiumMetabolismAlignment
#' @import methods
#' @importClassesFrom TreeSummarizedExperiment TreeSummarizedExperiment
newConsortiumMetabolismAlignment <- setClass(
    Class = "ConsortiumMetabolismAlignment",
    contains = "TreeSummarizedExperiment",
    slots = list(
        ## -- Metadata --
        Name = "character",
        Description = "character",
        Type = "character",
        Metric = "character",
        Params = "list",
        ## -- References --
        QueryName = "character",
        ReferenceName = "character",
        ## -- Scores --
        Scores = "list",
        PrimaryScore = "numeric",
        Pvalue = "numeric",
        ## -- Pairwise pathway detail --
        SharedPathways = "data.frame",
        UniqueQuery = "data.frame",
        UniqueReference = "data.frame",
        ## -- Multiple alignment --
        SimilarityMatrix = "matrix",
        ConsensusPathways = "data.frame",
        Prevalence = "data.frame",
        Dendrogram = "list",
        ## -- Carried over --
        Pathways = "data.frame",
        Graphs = "list",
        Metabolites = "data.frame"
    )
)

## ---- Validity: ConsortiumMetabolism ----
setValidity("ConsortiumMetabolism", function(object) {
    errors <- character()

    ## Name must be a length-1 character (NA is allowed for unnamed objects)
    if (length(object@Name) != 1L) {
        errors <- c(
            errors,
            "'Name' must be a length-1 character (NA allowed)"
        )
    }

    ## Weighted must be length-1 non-NA logical
    if (
        length(object@Weighted) != 1L ||
            is.na(object@Weighted)
    ) {
        errors <- c(
            errors,
            "'Weighted' must be a length-1 non-NA logical"
        )
    }

    ## Metabolites must be character (can be length 0)
    if (!is.character(object@Metabolites)) {
        errors <- c(
            errors,
            "'Metabolites' must be a character vector"
        )
    }

    ## Pathways must be a data.frame
    if (!is.data.frame(object@Pathways)) {
        errors <- c(
            errors,
            "'Pathways' must be a data.frame"
        )
    }

    ## Graphs must be a list
    if (!is.list(object@Graphs)) {
        errors <- c(
            errors,
            "'Graphs' must be a list"
        )
    }

    if (length(errors) == 0L) TRUE else errors
})

## ---- Validity: ConsortiumMetabolismSet ----
setValidity("ConsortiumMetabolismSet", function(object) {
    errors <- character()

    ## Name must be a length-1 character (NA is allowed for unnamed objects)
    if (length(object@Name) != 1L) {
        errors <- c(
            errors,
            "'Name' must be a length-1 character (NA allowed)"
        )
    }

    ## Consortia must be a list
    if (!is.list(object@Consortia)) {
        errors <- c(
            errors,
            "'Consortia' must be a list"
        )
    } else if (length(object@Consortia) > 0L) {
        ## If non-empty, all elements must be CM objects
        is_cm <- vapply(
            object@Consortia,
            function(x) is(x, "ConsortiumMetabolism"),
            logical(1L)
        )
        if (!all(is_cm)) {
            errors <- c(
                errors,
                paste0(
                    "All elements of 'Consortia' must be ",
                    "'ConsortiumMetabolism' objects"
                )
            )
        }
    }

    ## OverlapMatrix must be a matrix
    if (!is.matrix(object@OverlapMatrix)) {
        errors <- c(
            errors,
            "'OverlapMatrix' must be a matrix"
        )
    } else if (
        length(object@OverlapMatrix) > 0L &&
            nrow(object@OverlapMatrix) != ncol(object@OverlapMatrix)
    ) {
        errors <- c(
            errors,
            "'OverlapMatrix' must be square"
        )
    }

    ## Dendrogram must be a list
    if (!is.list(object@Dendrogram)) {
        errors <- c(
            errors,
            "'Dendrogram' must be a list"
        )
    }

    ## NodeData must be a data.frame
    if (!is.data.frame(object@NodeData)) {
        errors <- c(
            errors,
            "'NodeData' must be a data.frame"
        )
    }

    ## Pathways must be a data.frame
    if (!is.data.frame(object@Pathways)) {
        errors <- c(
            errors,
            "'Pathways' must be a data.frame"
        )
    }

    if (length(errors) == 0L) TRUE else errors
})

## ---- Validity: ConsortiumMetabolismAlignment ----
setValidity("ConsortiumMetabolismAlignment", function(object) {
    errors <- character()

    ## Helper: check if slot has a non-NA scalar value
    .hasValue <- function(x) {
        length(x) == 1L && !is.na(x)
    }

    ## Type must be one of the allowed values
    valid_types <- c("pairwise", "multiple", "search")
    if (
        .hasValue(object@Type) &&
            !object@Type %in% valid_types
    ) {
        errors <- c(
            errors,
            paste0(
                "Type must be one of: ",
                paste(valid_types, collapse = ", ")
            )
        )
    }

    ## PrimaryScore must be in [0, 1] or NA
    if (
        .hasValue(object@PrimaryScore) &&
            (object@PrimaryScore < 0 || object@PrimaryScore > 1)
    ) {
        errors <- c(
            errors,
            "PrimaryScore must be between 0 and 1"
        )
    }

    ## Pvalue must be in [0, 1] or NA
    if (
        .hasValue(object@Pvalue) &&
            (object@Pvalue < 0 || object@Pvalue > 1)
    ) {
        errors <- c(
            errors,
            "Pvalue must be between 0 and 1"
        )
    }

    ## Pairwise: QueryName and ReferenceName required
    if (.hasValue(object@Type) && object@Type == "pairwise") {
        if (!.hasValue(object@QueryName)) {
            errors <- c(
                errors,
                "QueryName required for pairwise alignment"
            )
        }
        if (!.hasValue(object@ReferenceName)) {
            errors <- c(
                errors,
                "ReferenceName required for pairwise alignment"
            )
        }
    }

    ## SimilarityMatrix shape depends on Type
    sm <- object@SimilarityMatrix
    if (length(sm) > 0L) {
        if (
            .hasValue(object@Type) &&
                object@Type == "multiple" &&
                nrow(sm) != ncol(sm)
        ) {
            errors <- c(
                errors,
                "SimilarityMatrix must be square for \\
                multiple alignment"
            )
        }
        if (
            .hasValue(object@Type) &&
                object@Type == "search" &&
                nrow(sm) != 1L
        ) {
            errors <- c(
                errors,
                "SimilarityMatrix must have exactly 1 row \\
                for search (query x database)"
            )
        }
    }

    if (length(errors) == 0L) TRUE else errors
})
