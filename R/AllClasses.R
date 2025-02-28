#' @noRd
#' @import methods
#' @importClassesFrom TreeSummarizedExperiment TreeSummarizedExperiment
newConsortiumMetabolism <- setClass(
  Class = "ConsortiumMetabolism",
  contains = "TreeSummarizedExperiment",
  slots = list(
    Name = "character",
    Edges = "list",
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
  ### don't think I need to import TSE here
  # contains = "TreeSummarizedExperiment",
  slots = list(
    Name = "character",
    Communities = "list",
    Description = "character"
  ),
  prototype = list(
    Name = NA_character_,
    Communities = list(),
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
    Score = "data.frame"
  )
)
