% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MRCIGVAR_b.R
\name{girf_MRCIGVAR_RM_CB}
\alias{girf_MRCIGVAR_RM_CB}
\title{Generalized impulse response functions of MRCIGVAR(m,p,n,S,W,TH) with regime migrations}
\usage{
girf_MRCIGVAR_RM_CB(
  res,
  shock,
  R,
  nstep,
  Omega_hist = NA,
  resid_method = "parametric",
  conf,
  N
)
}
\arguments{
\item{res}{: a MRCIGVAR object containing the components of an output of MRCIGVARestm.}

\item{shock}{: an mn-vector containing the shocks as impulse.}

\item{R}{: the number of runs to integrate out the random effects in order to obtain the means (see equation above).}

\item{nstep}{: the length of the responses}

\item{Omega_hist}{: the initial values from which the simulation runs of impulse and response functions start}

\item{resid_method}{: resid_method = c("resid", "parametric"), It generate random residuals either from residuals bootstrap or parametric bootstrap.}

\item{conf}{: a two component vector containing the tail probabilities of the bootstrap confidence interval.}

\item{N}{: number of bootstrapping runs}
}
\value{
an (n x n x nstep+1 x 3) array containing of impulse response functions with lower and upper confidence bonds. The rows represent response the columns represent impulses.
}
\description{
This function calculates the generalized impulse response functions of an estimated MRVAR(n,p,S) for given a shock vector and initial values.
}
\details{
\preformatted{              GIRF(shock=SHCK) = mean(Y(resid)) - mean(Y(SHCK))
}

It also generates the bootstrapped confidence intervals.
}
\examples{
m = 2
n = 3
p = c(2,2,2,2,2,2,2,2,2,2,0,0,0,0,0,2,2,2,2,2,2,2,2,2,2,0,0,0,0,0); dim(p) = c(5,3,2)
p = p[1:n,,]; p[,1,] = 3; p[,2,] = 2

TH = c(1:n)*0; dim(TH) = c(1,n)
SESVI=rep(1,3,5)
r  = rep(1,n)

## case of n = 3, m = 2, S = 2    m: number of variables, n: number of countries
res_d <- MRCIGVARData(m=2,n=3,p=p,TH=TH,T=200,S=2, SESVI=c(1,3,5),r=rep(1,3),Ncommtrend=1)
max(abs(res_d$Y))
STAT(res_d$Go[,,,2])
STAT(res_d$Go[,,,1])
res_e  = MRCIGVARestm(res=res_d)

STAT(res_e$Go[,,,2])
STAT(res_e$Go[,,,1])

plot(ts(res_d$Y))
#res_e$Summary

if (!((max(Mod(STAT(res_e$Go[,,,1])))>1)|(max(Mod(STAT(res_e$Go[,,,2])))>1)) ) {
  GIRF <- girf_MRCIGVAR_RM(res=res_e,shock=c(1,1,1,1,1,1),R=100,nstep=10,Omega_hist=NA,resid_method="parametric")
  GIRF_CB <- girf_MRCIGVAR_RM_CB(res=res_e,shock=c(1,1,1,1,1,1),R=100,nstep=10,Omega_hist=NA,resid_method="parametric",conf=c(0.05,0.95),N=100)
  IRF_g = IRF_graph(GIRF_CB)
}

}
