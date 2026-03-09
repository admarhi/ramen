# Get the Consortia

Returns the consortia data in a tabular format.

## Usage

``` r
getCo(object)

# S4 method for class 'ConsortiumMetabolism'
getCo(object)
```

## Arguments

- object:

  An object of class `ConsortiumMetabolism`

## Value

For `ConsortiumMetabolism` objects, returns a tb with species,
metabolite and flux information. For `ConsortiumMetabolismSet` and
`ConsortiumMetabolismAlignment` objects, returns a list of such tibbles.

A list with the community data.

## Methods (by class)

- `getCo(ConsortiumMetabolism)`: Get the Community
