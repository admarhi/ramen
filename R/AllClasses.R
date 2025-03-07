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
newConsortiumMetabolismSet <- setClass(
  Class = "ConsortiumMetabolismSet",
  slots = list(
    Name = "character",
    Consortia = "list",
    Description = "character"
  ),
  prototype = list(
    Name = NA_character_,
    Consortia = list(),
    Description = NA_character_
  )
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
    Communities = "list",
    Score = "data.frame",
    Graphs = "list"
  )
)
