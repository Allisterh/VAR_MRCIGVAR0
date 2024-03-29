% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/VARs_functions_b.R
\name{IRF_graph}
\alias{IRF_graph}
\title{Plot impulse response functions}
\usage{
IRF_graph(
  IRF_CB = IRF_CB,
  Names = NA,
  INames = NA,
  response = c(1:n),
  impulse = c(1:n),
  ncol = n
)
}
\arguments{
\item{IRF_CB}{An (n x n x L x 3) array of impulse response function with confidence bands}

\item{Names}{An n-vector of strings of the variable names}

\item{INames}{An n-vector of string of the impulse names}

\item{response}{An vector of impulse indices}

\item{impulse}{An vector of response indices}

\item{ncol}{Number of columns of impulse response functions in the plot}
}
\value{
An ggplot object of impulse response functions
}
\description{
Plot impulse response functions
}
