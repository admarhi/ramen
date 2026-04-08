# Plot a Directed Graph Emphasizing Flow

Arranges nodes of a directed graph into three distinct columns
representing sources (nodes with only outgoing edges), sinks (nodes with
only incoming edges), and intermediate nodes (with both incoming and
outgoing edges). The layout within the intermediate column is determined
using the Fruchterman-Reingold algorithm. This visualization helps
understand the overall flow structure within the graph.

## Usage

``` r
plotDirectedFlow(
  g,
  source_x = 0,
  mixed_x = 1,
  sink_x = 2,
  vertical_spacing = 1,
  vertex_size = 10,
  vertex_label_cex = 0.8,
  edge_arrow_size = 0.5,
  edge_width_range = c(0.1, 0.1),
  color_edges_by_weight = FALSE,
  edge_color_low = "gray80",
  edge_color_high = "black",
  ...
)
```

## Arguments

- g:

  An igraph object. Must be a directed graph.

- source_x:

  Numeric scalar. The x-coordinate for source nodes. Defaults to 0.

- mixed_x:

  Numeric scalar. The central x-coordinate for intermediate nodes. The
  actual layout spans `mixed_x` +/- 0.4. Defaults to 1.

- sink_x:

  Numeric scalar. The x-coordinate for sink nodes. Defaults to 2.

- vertical_spacing:

  Numeric scalar. The maximum y-coordinate, controlling the vertical
  spread of the layout. Defaults to 1.

- vertex_size:

  Numeric scalar. The size of the vertices (nodes) in the plot. Defaults
  to 10.

- vertex_label_cex:

  Numeric scalar. The character expansion factor for vertex labels.
  Defaults to 0.8.

- edge_arrow_size:

  Numeric scalar. The size of the arrows on the edges. Defaults to 0.5.

- edge_width_range:

  Numeric vector of length 2. The minimum and maximum width for edges
  when scaled by weight. If the graph is unweighted or all weights are
  identical, the mean of this range is used. Defaults to `c(0.1, 0.1)`.

- color_edges_by_weight:

  Logical scalar. If `TRUE` and the graph has a 'weight' edge attribute,
  edges will be colored based on their weight, interpolating between
  `edge_color_low` and `edge_color_high`.

- edge_color_low:

  Character string. The color for the lowest edge weight when
  `color_edges_by_weight` is `TRUE`. Defaults to "gray80".

- edge_color_high:

  Character string. The color for the highest edge weight when
  `color_edges_by_weight` is `TRUE`. Defaults to "black".

- ...:

  Additional arguments passed to igraph.

## Value

Invisibly returns `NULL`. This function is called for its side effect of
generating a plot.

## Examples

``` r
cm <- synCM("test", n_species = 3, max_met = 5)
g <- igraph::graph_from_adjacency_matrix(
    SummarizedExperiment::assays(cm)[["Binary"]],
    mode = "directed",
    weighted = TRUE
)
plotDirectedFlow(g)

```
