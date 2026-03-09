# ConsortiumMetabolismSet

``` r
library(ramen)
```

## Class

A MFS can be created by uniting any number of MF objects together by
calling ConsortiumMetabolismSet() or MFS() on either a list of MF
objects or inidividual MF objects. This process it self does not trigger
any processing of the data, the only thing that happens is that metadata
of the MFS is calculated and an UUID is generated. If upon creation, the
mode of alignment is already known it can be supplied to the constructor
function of the MF as follows: ConsortiumMetabolismSet(…,
alignment=“multiple”)

## ToDo

- Inherit the MF class and methods
- Write `removeMF()` and `addMF()` functions
