% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/VARs_functions_b.R
\name{abMRCVECMest2testm}
\alias{abMRCVECMest2testm}
\title{Estimation of constrained CIVAR models}
\usage{
abMRCVECMest2testm(
  y,
  x,
  s,
  model = c("I", "II", "III", "IV", "V"),
  type = c("eigen", "trace"),
  ret = c("statistic", "test"),
  ctable = c("A3", "A1", "A5"),
  crk = crk,
  P = matrix(2, 2, 2),
  q = 0.95,
  Dxflag = 0,
  H = H,
  h = h,
  phi = phi,
  G = G,
  psi = psi
)
}
\arguments{
\item{y}{The series of endogenous variables}

\item{x}{The series of conditioning/exogeneous variables}

\item{s}{The series of the regime indicator function}

\item{model}{Types of deterministic components in CIVAR models}

\item{type}{Type of the Johansen test}

\item{ret}{Statistics}

\item{ctable}{Critical value tables}

\item{crk}{Cointegration rank}

\item{P}{Lags of MRCIVAR}

\item{q}{Significance level}

\item{Dxflag}{A flag specifying whether the conditioning variablesare within the cointegration space or not}

\item{H}{A matrix specifying restrictions on beta}

\item{h}{A vector specifying restrictions on beta}

\item{phi}{freely varying parameters in beta}

\item{G}{A matrix specifying restrictions on alpha}

\item{psi}{freely varying parameters in alpha}
}
\value{
A list of estimated constrained MRCIVAR
}
\description{
Estimation of constrained CIVAR models
}
