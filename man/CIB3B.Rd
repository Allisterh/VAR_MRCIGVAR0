% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CIVAR_b.r
\name{CIB3B}
\alias{CIB3B}
\title{Transformation of the coefficients matrix of a VECM to the coefficients of the corresponding VAR in level}
\usage{
CIB3B(CIB)
}
\arguments{
\item{CIB}{an (m,m,L) array of the coefficients matrices of a VECM with lag L-1.}
}
\value{
B (m,m,L) array the coefficients of a VAR with lag L
}
\description{
Transformation of the coefficients matrix of a VECM to the coefficients of the corresponding VAR in level
}
