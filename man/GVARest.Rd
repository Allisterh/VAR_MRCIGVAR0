% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GVAR_b.R
\name{GVARest}
\alias{GVARest}
\title{Estimation of GVAR(m,n,p)}
\usage{
GVARest(res)
}
\arguments{
\item{res}{: an GVAR object that is an output of GVARData including at least: m,n,p,type,Y and optionally X.}
}
\value{
res  : an GVAR object with estimated parameter values, AIC, BIC, AIC_g, BIC_g and LH, where AIC and BIC are n-vectors of the country equations' AIC and BIC and AIG_g amd BIC_g are the GVAR information criteria respectively.
}
\description{
This function estimates the parameters of a specified GVAR(m,n,p) model based on provided data.
}
\examples{
n = 5
p = (1:15)*0; dim(p) = c(5,3)
p[,1] = 2; p[,2]=1;
res_d = GVARData(m=2,n=5,p=p,T=100,type="const")
max(res_d$Y)
dim(res_d$Y)
res_e = GVARest(res = res_d)
res_e$Summary

X1 = matrix(1,200,1)
X2 = matrix(rnorm(200),200,1)
X3 = matrix(rnorm(200),200,1)
X4 = matrix(rnorm(200),200,1)
X  = cbind(X1,X2,X3,X4)
dim(X) = c(200,1,4)
n = 4
p = (1:12)*0; dim(p) = c(4,3);p[,1] = 2; p[,2]=1;   p[,3]=1; p[2,2]=2;
p

res_d = GVARData(m=2,n=4,p=p,T=200,type="exog0",X=X)
res_e = GVARest(res = res_d)
res_e$Summary

IRF_CB = irf_GVAR_CB(res_e,nstep=10,comb=NA,irf="gen",runs=200,conf=c(0.05,0.95))
dim(IRF_CB)
IRF_g = IRF_graph(IRF_CB,Names=NA,response=c(1,4),impulse=c(1,2,3,4), ncol=4)


}
