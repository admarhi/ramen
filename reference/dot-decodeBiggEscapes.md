# Decode COBRApy/SBML `__NN__` escape sequences

COBRApy encodes non-alphanumeric characters in reaction and metabolite
IDs as `__<decimal-ASCII-code>__` sequences when writing SBML. For
example, `(` becomes `__40__` and `)` becomes `__41__`. This helper
reverses that encoding so downstream parsing sees the literal
characters.

## Usage

``` r
.decodeBiggEscapes(ids)
```

## Arguments

- ids:

  Character vector of identifiers possibly containing `__NN__` escape
  sequences.

## Value

Character vector with escape sequences replaced by their literal
characters. IDs with no escape sequences pass through unchanged.
