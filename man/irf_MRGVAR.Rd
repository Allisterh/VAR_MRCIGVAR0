% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MRGVAR_b.R
\name{irf_MRGVAR}
\alias{irf_MRGVAR}
\title{Regime specific impulse response functions of MRGVAR(n,p,S)}
\usage{
irf_MRGVAR(
  res = res,
  state = state,
  nstep = nstep,
  comb = comb,
  irf = c("gen", "chol", "chol1", "gen1", "genN1", "comb1", "smat"),
  G = NA,
  smat = NA,
  sigmaNPDS = NA
)
}
\arguments{
\item{res}{a list of estimated MRGVAR as output of MRGVARest}

\item{state}{an n vector specifying the specific state for each country.}

\item{nstep}{the length of impulse response function}

\item{comb}{a vector specify the concerted action in policy-simulation impulse response function}

\item{irf}{: types of the impulse response irf=c("gen","chol","chol1","gen1","comb1")
"gen" for generalized impulse response with one standard deviation impulses, "gen1" for GIRF with one unit impulses, "chol" Cholezky decomposition, "chol1" Cholezky decomposition with unit impulses, "comb1" concerted action with unit impulse.}

\item{G}{For permanent and transitory decomposition}

\item{smat}{For explicit structural decomposition of the correlated shocks}

\item{sigmaNPDS}{the state-dependent covariance matrix}
}
\value{
a list containing the impulse response functions and the accumulated impulse response function, and the boostrap parameters as well.
}
\description{
This function calculates the regime specific impulse response functions of an estimated MRGVAR(n,p,S).
Using the estimated G[,,,s] and Sigma[,,s] matrices of the MRGVAR, this function calculated the regime specific impulse response functions.
}
\examples{

## case of n = 2, m = 2, S = 2     ## m: number of variables, n: number of countries
p = rep(1,12); dim(p) = c(2,3,2)
p[1,1,2] = 2; p[2,2,2]=2; p[,3,] = 0
TH = c(1:2)*0; dim(TH) = c(1,2)
res_d <- MRGVARData(m=2,n=2,p=p,TH=TH,T=100,S=2,SESVI=c(1,3),type="const")
max(res_d$Y)

### estimation of the MRGVAR model
colnames(res_d$Y) = c("P","Q","Pa","Qa")
res_e = MRGVARest(res=res_d)
res_e$Summary

IRF  = irf_MRGVAR(res=res_e,nstep=10,comb=NA,state=c(1,1),irf="gen1")

}
