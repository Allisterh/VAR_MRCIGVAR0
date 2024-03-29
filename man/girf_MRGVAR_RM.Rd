% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MRGVAR_b.R
\name{girf_MRGVAR_RM}
\alias{girf_MRGVAR_RM}
\title{Generalized impulse response functions of MRGVAR(m,n,p,S) with regime migrations}
\usage{
girf_MRGVAR_RM(res, shock, R, nstep, Omega_hist, resid_method)
}
\arguments{
\item{res}{: an MRGVAR object containing the components of the output of MRGVARData or MRGVARest.}

\item{shock}{: an mn-vector containing the shocks as impulse.}

\item{R}{: the number of runs to integrate out the random effects in order to obtain the means (see equation above).}

\item{nstep}{: the length of the responses}

\item{Omega_hist}{: the initial values from which the simulation runs start.For Omega_hist=NA the most recent values are taken as the initial values. For Omega_hist=0, the initial values are zeros.}

\item{resid_method}{: resid_method = c("resid", "parametric"), It generates the random residuals from residuals bootstrap or parametric bootstrap.}
}
\value{
an (mn x mn x nstep+1) matrix of impulse response functions. The rows represent response the columns represent impulses.
}
\description{
This function calculates the generalized impulse response functions of an estimated MRGVAR(n,p,S) for given a shock vector.
}
\details{
\preformatted{              GIRF(shock=SHCK) = mean(Y(resid)) - mean(Y(SHCK))
}
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

IRF_CB  = irf_MRGVAR_CB1(res=res_e,nstep=10,comb=NA,state=c(1,1),irf="gen1",runs=20,conf=c(0.05,0.95))
IRF_g = IRF_graph(IRF_CB[[1]],Names=c("P","Q","Pa","Qa"))    #IRF
IRF_g = IRF_graph(IRF_CB[[2]])   # accumulated IRF

GIRF    = girf_MRGVAR_RM(res=res_e,shock=c(1,1,1,1),R=100,nstep=10,Omega_hist=NA,resid_method='parametric')
GIRF_CB = girf_MRGVAR_RM_CB(res=res_e,shock=c(1,1,1,1),R=100,nstep=10,Omega_hist=NA,resid_method='parametric',conf=c(0.05,0.95),N=10)
IRF_g   = IRF_graph(GIRF_CB,Names=c("P","Q","Pa","Qa"))    #IRF
}
