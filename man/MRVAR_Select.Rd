% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MRVAR_b.R
\name{MRVAR_Select}
\alias{MRVAR_Select}
\title{Calculation of information criteria for MRVAR(n,p,S) models.}
\usage{
MRVAR_Select(res, L_V = L_V, TH_V = TH_V)
}
\arguments{
\item{res}{: an MRVAR object generated from MRVARData or estimated from MRVARest.}

\item{L_V}{: a 2 components vector containing the maxima lags in the two regimes, respectively.}

\item{TH_V}{: a vector containing the possible threshold values over which the model selection criteria values will be calculated.}
}
\value{
\preformatted{ a matrix with different lag specifications and threshold values as well as the corresponding information criterion values.
}
}
\description{
Calculation of information criteria for MRVAR(n,p,S) models.
}
\examples{

res_d = MRVARData(n=4,p=matrix(c(2,1,2,2,0,0,0,0),4,2),T=800,S=2,SESVI=1)
max(res_d$Y)
colnames(res_d$Y) = c("P","Y","R","U")
res_e = MRVARest(res=res_d)
TH_v = c(0,0.0)
L_v = c(5,5)
Sel = MRVAR_Select(res=res_e,L_V=L_v,TH_V=TH_v)
MRVAR_Select_Summary(Sel)
}
