# Ramen plot theme

Shared ggplot2 theme used across every plot returned by the package.
Provides a consistent typography, white background, bottom legend, and
tight title spacing. Network plots pass `network = TRUE` to drop axes,
ticks, and grid lines while keeping the rest of the look-and-feel.

## Usage

``` r
theme_ramen(baseSize = 11, baseFamily = "", network = FALSE)
```

## Arguments

- baseSize:

  Numeric. Base font size in points. Defaults to 11.

- baseFamily:

  Character. Base font family. Defaults to `""` (device default).

- network:

  Logical. If `TRUE`, build on top of
  [`theme_void()`](https://ggplot2.tidyverse.org/reference/ggtheme.html)
  instead of
  [`theme_minimal()`](https://ggplot2.tidyverse.org/reference/ggtheme.html);
  used by
  [`plotDirectedFlow`](https://admarhi.github.io/ramen/reference/plotDirectedFlow.md)
  and `plot` methods that render networks.

## Value

A `ggplot2` theme object.

## Examples

``` r
library(ggplot2)
ggplot(mtcars, aes(mpg, hp)) +
    geom_point() +
    theme_ramen()
```
