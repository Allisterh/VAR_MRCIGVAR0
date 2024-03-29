% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CIVAR_b.r
\name{CIB2B}
\alias{CIB2B}
\title{Transformation of coefficient matrix of a VECM to the coefficients matrix of the corresponding CIGVAR in level.}
\usage{
CIB2B(tst)
}
\arguments{
\item{tst}{: an output of CIGVARest}
}
\value{
\preformatted{ A list containing three components.
}

\itemize{
\item B          : the coefficient matrices of the domestic variables
\item A	        : the coefficient matrices of the foreign variables
\item C          : the coefficient matrices of the deterministic components
}
}
\description{
Transformation of coefficient matrix of a VECM to the coefficients matrix of the corresponding CIGVAR in level.
}
