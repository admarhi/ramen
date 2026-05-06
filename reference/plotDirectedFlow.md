# Plot a Directed Graph Emphasising Flow

Arranges nodes of a directed graph into three distinct columns
representing sources (nodes with only outgoing edges), sinks (nodes with
only incoming edges), and intermediate nodes (with both incoming and
outgoing edges). The layout within the intermediate column is determined
using the Kamada-Kawai algorithm. This visualisation helps understand
the overall flow structure within the graph and is rendered with ggraph,
so the returned object composes with the usual ggplot2 operators.

## Usage

``` r
plotDirectedFlow(
  g,
  sourceX = 0,
  mixedX = 1,
  sinkX = 2,
  verticalSpacing = 1,
  nodeSize = 6,
  nodeLabelSize = 3,
  edgeArrowSize = 2,
  edgeWidthRange = c(0.3, 1.2),
  colourEdgesByWeight = FALSE,
  edgeColourAttr = NULL,
  edgeColourValues = NULL,
  edgeColourLabels = NULL,
  edgeColourLegendTitle = "Type",
  edgeColourLow = "#DEEBF7",
  edgeColourHigh = "#08519C",
  nodeColourValues = c(source = "#0072B2", intermediate = "#F0E442", sink = "#D55E00"),
  nodeColourLabels = c(source = "Source", intermediate = "Intermediate", sink = "Sink"),
  nodeColourLegendTitle = "Node role",
  main = NULL,
  ...
)
```

## Arguments

- g:

  An igraph object. Must be a directed graph.

- sourceX:

  Numeric scalar. The x-coordinate for source nodes. Defaults to 0.

- mixedX:

  Numeric scalar. The central x-coordinate for intermediate nodes. The
  actual layout spans `mixedX` +/- 0.4. Defaults to 1.

- sinkX:

  Numeric scalar. The x-coordinate for sink nodes. Defaults to 2.

- verticalSpacing:

  Numeric scalar. The maximum y-coordinate, controlling the vertical
  spread of the layout. Defaults to 1.

- nodeSize:

  Numeric scalar. The size of node points in millimetres (ggraph unit).
  Defaults to 6.

- nodeLabelSize:

  Numeric scalar. The size of node labels in millimetres. Defaults to 3.

- edgeArrowSize:

  Numeric scalar. The arrow length on edges, in millimetres. Defaults to
  2.

- edgeWidthRange:

  Numeric vector of length 2. The minimum and maximum width for edges
  when scaled by weight. If the graph is unweighted or all weights are
  identical, the mean of this range is used. Defaults to `c(0.3, 1.2)`.

- colourEdgesByWeight:

  Logical scalar. If `TRUE` and the graph has a `weight` edge attribute,
  edges are coloured on a continuous gradient from `edgeColourLow` to
  `edgeColourHigh`.

- edgeColourAttr:

  Optional character scalar naming a categorical edge attribute to
  colour edges by (e.g. `"source"`). When set, takes precedence over
  `colourEdgesByWeight`.

- edgeColourValues:

  Optional named character vector mapping levels of `edgeColourAttr` to
  colours. Used only when `edgeColourAttr` is non-NULL.

- edgeColourLabels:

  Optional named character vector mapping levels of `edgeColourAttr` to
  legend labels. Used only when `edgeColourAttr` is non-NULL.

- edgeColourLegendTitle:

  Character scalar. Legend title for the categorical edge colour scale.
  Defaults to `"Type"`.

- edgeColourLow:

  Character string. The colour for the lowest edge weight when
  `colourEdgesByWeight` is `TRUE`. Defaults to `"#DEEBF7"` (ColorBrewer
  Blues, light end).

- edgeColourHigh:

  Character string. The colour for the highest edge weight when
  `colourEdgesByWeight` is `TRUE`. Defaults to `"#08519C"` (ColorBrewer
  Blues, dark end).

- nodeColourValues:

  Named character vector of length 3 mapping the node roles `"source"`,
  `"intermediate"`, and `"sink"` to colours. Defaults to a maximally
  distinct Okabe-Ito triple (blue / yellow / vermilion) chosen for
  colour-blind safety and high mutual contrast.

- nodeColourLabels:

  Optional named character vector mapping node roles to legend labels.
  Defaults to
  `c(source = "Source", intermediate = "Intermediate", sink = "Sink")`.

- nodeColourLegendTitle:

  Character scalar. Legend title for the node-role colour scale.
  Defaults to `"Node role"`.

- main:

  Character string. Optional plot title.

- ...:

  Additional arguments. Currently unused; retained for forward
  compatibility.

## Value

A `ggplot` object.

## Details

Edges may be coloured in three ways: by a continuous `weight` attribute
(`colourEdgesByWeight = TRUE`); by a categorical edge attribute
(`edgeColourAttr = "<name>"`), in which case `edgeColourValues` and
`edgeColourLabels` customise the palette and legend; or with a single
fixed colour (default).

When the graph carries non-positive edge weights, the Kamada-Kawai
layout used for intermediate nodes falls back to topological placement
(the `weight` attribute is dropped on the induced intermediate subgraph,
since
[`igraph::layout_with_kk`](https://r.igraph.org/reference/layout_with_kk.html)
requires strictly positive weights, and `EffectiveConsumption`-style
assays carry zeros). This affects only the within-column layout; the
source / sink / intermediate column assignment is unaffected.

The size-bearing arguments (`nodeSize`, `nodeLabelSize`,
`edgeArrowSize`) use ggraph millimetre units rather than the igraph
`cex` factors of the previous implementation.

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
