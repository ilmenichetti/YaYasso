% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/litterfluxes_AWEN.r
\name{root.litter}
\alias{root.litter}
\title{Calculates root litter}
\usage{
root.litter(Mr, spec)
}
\arguments{
\item{Mr}{Mass roots (usually in Mg ha$^{-1}$)}

\item{spec}{tree species, 1=pine 2=spruce, 3=residues (broadleaves)}
}
\value{
litter inputs (mass) from a mass of roots, to be passed to the AWEN partitioning function
}
\description{
Calculates root litter
}
\seealso{
\code{\link{branches.AWEN}}
}
