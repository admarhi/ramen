% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AllGenerics.R, R/methods-getSpecies.R
\name{getSpecies}
\alias{getSpecies}
\alias{getSpecies,ConsortiumMetabolism-method}
\alias{getSpecies,ConsortiumMetabolismSet-method}
\alias{getSpecies,ConsortiumMetabolismAlignment-method}
\title{Return Species in a Consortium}
\usage{
getSpecies(object, ...)

\S4method{getSpecies}{ConsortiumMetabolism}(object)

\S4method{getSpecies}{ConsortiumMetabolismSet}(object, type = c("all", "generalists", "specialists"))

\S4method{getSpecies}{ConsortiumMetabolismAlignment}(object)
}
\arguments{
\item{object}{a \code{ConsortiumMetabolismAlignment} Object}

\item{...}{Object specific arguments.}

\item{type}{Character scalar giving the type of species to output.}
}
\value{
A character vector containing the names of species in the consortium

A character vector representing the microorganisms.

A character vector representing the microorganisms.

A character vector representing the microorganisms.
}
\description{
Return Species in a Consortium
}
\section{Methods (by class)}{
\itemize{
\item \code{getSpecies(ConsortiumMetabolism)}: Return Species in a Microbiome

\item \code{getSpecies(ConsortiumMetabolismSet)}: Return Species in a Microbiome

\item \code{getSpecies(ConsortiumMetabolismAlignment)}: Return Species in a Microbiome

}}
