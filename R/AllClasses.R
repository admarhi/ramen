#' @noRd
#' @import methods
#' @importClassesFrom TreeSummarizedExperiment TreeSummarizedExperiment
newConsortiumMetabolism <- setClass(
  Class = "ConsortiumMetabolism",
  contains = "TreeSummarizedExperiment",
  slots = list(
    Name = "character",
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
    Edges = "data.frame",
    Metabolites = "data.frame"
  ) #,
  # prototype = list(
  #   Name = NA_character_,
  #   Consortia = list(),
  #   Description = NA_character_
  # )
)

#' @noRd
#' @import methods
#' @importClassesFrom TreeSummarizedExperiment TreeSummarizedExperiment
newConsortiumMetabolismAlignment <- setClass(
  Class = "ConsortiumMetabolismAlignment",
  contains = "TreeSummarizedExperiment",
  slots = list(
    Name = "character",
    Alignment = "list",
    Edges = "data.frame",
    Consortia = "list",
    Score = "data.frame",
    Graphs = "list",
    Metabolites = "data.frame",
    Dendrogram = "list"
  )
)
