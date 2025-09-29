#' @describeIn align Align a \code{ConsortiumMetabolismSet} Object
#' @export
setMethod("align", "ConsortiumMetabolismSet", function(object, name) {
  cli::cli_h1(paste0("Aligning ", object@Name))

  # ---- Construct Alignment Object --------------------------------------------
  newConsortiumMetabolismAlignment(
    TreeSummarizedExperiment(
      # assays = list(Levels = levels_mat),
      # colData = new_met_ind
    ),
    # Name = object@Name,
    # Edges = all_edges,
    # # Graphs = graph_list,
    # Metabolites = all_met,
    # ### This should really not store consortia again but only a list of the
    # ### names.
    # Consortia = object@Consortia
    # ### Need
    # # - Number of consortia
    # # -
  )
})
