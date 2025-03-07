setMethod("plot", "ConsortiumMetabolism", function(x) {
  plot(x@Graphs[[1]])
})


setMethod("plot", "ConsortiumMetabolismAlignment", function(x, type) {
  if (type == "venn-vertex") {
    vertex_list <- lapply(x@Graphs, function(g) igraph::V(g)$name)
    ggVennDiagram::ggVennDiagram(vertex_list)
  } else if (type == "venn-edge") {
    edge_list <- lapply(x@Graphs, function(g) igraph::E(g))
    ggVennDiagram::ggVennDiagram(edge_list)
  }
})
