#' @noRd
#' @import methods
#' @importClassesFrom TreeSummarizedExperiment TreeSummarizedExperiment
newConsortiumMetabolism <- setClass(
  Class = "ConsortiumMetabolism",
  contains = "TreeSummarizedExperiment",
  slots = list(
    Name = "character",
    Description = "character",
    Edges = "data.frame",
    Weighted = "logical",
    InputData = "data.frame",
    Metabolites = "character",
    Graphs = "list"
  )
)

#' @noRd
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
        Edges = "data.frame",
        Metabolites = "data.frame"
    )
)

#' @noRd
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
        ConsensusEdges = "data.frame",
        Prevalence = "data.frame",
        Dendrogram = "list",
        ## -- Carried over --
        Edges = "data.frame",
        Graphs = "list",
        Metabolites = "data.frame"
    )
)

setValidity("ConsortiumMetabolismAlignment", function(object) {
    errors <- character()

    ## Helper: check if slot has a non-NA scalar value
    .hasValue <- function(x) {
        length(x) == 1L && !is.na(x)
    }

    ## Type must be one of the allowed values
    valid_types <- c("pairwise", "multiple", "search")
    if (.hasValue(object@Type) &&
        !object@Type %in% valid_types) {
        errors <- c(
            errors,
            paste0(
                "Type must be one of: ",
                paste(valid_types, collapse = ", ")
            )
        )
    }

    ## PrimaryScore must be in [0, 1] or NA
    if (.hasValue(object@PrimaryScore) &&
        (object@PrimaryScore < 0 || object@PrimaryScore > 1)) {
        errors <- c(
            errors,
            "PrimaryScore must be between 0 and 1"
        )
    }

    ## Pvalue must be in [0, 1] or NA
    if (.hasValue(object@Pvalue) &&
        (object@Pvalue < 0 || object@Pvalue > 1)) {
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

    ## Multiple: SimilarityMatrix must be square if non-empty
    sm <- object@SimilarityMatrix
    if (length(sm) > 0L && nrow(sm) != ncol(sm)) {
        errors <- c(
            errors,
            "SimilarityMatrix must be square"
        )
    }

    if (length(errors) == 0L) TRUE else errors
})
