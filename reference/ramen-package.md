# ramen: Reconstruction and Alignment of Microbial Exchange Networks

Reconstructs metabolic exchange networks from microbial community data
and aligns them using functional overlap scores to quantify the
similarity and dissimilarity of consortium-level metabolic function.
Provides tools for visualization and comparison of microbial consortia
from a network perspective.

## See also

Useful links:

- <https://admarhi.github.io/ramen/>

- Report bugs at <https://github.com/admarhi/ramen/issues>

## Author

**Maintainer**: Adrian Hirt <am.hirt@pm.me>
([ORCID](https://orcid.org/0009-0008-1929-6629))

## Examples

``` r
# Create a synthetic consortium
cm <- synCM("example", n_species = 3, max_met = 5)
cm
#> 
#> ── ConsortiumMetabolism 
#> Name: "example"
#> Weighted metabolic network: 3 species, 4 metabolites, 6 pathways.
```
