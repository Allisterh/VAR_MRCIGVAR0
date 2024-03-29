% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/VAR_b.R
\name{VARData}
\alias{VARData}
\title{Data generating process of VAR(p)}
\usage{
VARData(n, p, T, r_np, A, B, Co, U, Sigma, type, X, mu, Yo)
}
\arguments{
\item{n}{: number of variables}

\item{p}{: lag length}

\item{T}{: number of observations}

\item{r_np}{: an n x p matrix of roots of the characteristic polynomials of n independent AR(p)-processes. If not provided, it will be generated randomly.}

\item{A}{: an n x n full rank matrix of transformation to generate correlated VAR(p) from the n independent AR(p)}

\item{B}{: (n,n,p) array of the coefficients of the VAR(p) process. If B is not given, it will be calculated out of r_np and A.}

\item{Co}{: (n,k+1) matrix of the coefficients of deterministic components in a VAR(p) process. For type="none" Co = 0*(1:n), for type="const" Co is an n vector, for type="exog0" Co is a (n,k) matrix, andfor type="exog1" Co is an (n,1+k) matrix.}

\item{U}{: residuals, if it is not NA it will be used as input to generate the VAR(p) process.}

\item{Sigma}{: an n x n covariance matrix of the VAR(p) residuals.}

\item{type}{: types of deterministic components. "none", "const" "exog0" and "exog1" are four options}

\item{X}{: a (T x k) matrix of exogenous variables.}

\item{mu}{: an n vector of the expected mean of the VAR(p) process}

\item{Yo}{: a p x n matrix of initial values of the VAR(p) process}
}
\value{
\preformatted{ An object of VAR(p) containing the generated data, the used parameters and the exogenous variables. res = list(n,p,type,r_np,Phi,A,B,Co,Sigma,Y,X,resid,U,Y1,Yo,check)
}
}
\description{
This function will generate data from a stationary VAR(p) process and return a list containing data and parameters used in the VAR(p) process.
}
\examples{
res_d = VARData(n=2,p=2,T=100,type="const")
res_d = VARData(n=2,p=2,T=10,Co=c(1:2)*0,type="none")
res_d = VARData(n=2,p=2,T=10,Co=c(1:2),  type="const")
res_d = VARData(n=3,p=2,T=200,type="exog1",X=matrix(rnorm(400),200,2))
res_d = VARData(n=3,p=2,T=200,Co=matrix(c(0,0,0,1,2,3,3,2,1),3,3), type="exog0",X=matrix(rnorm(400),200,2))

res_d = VARData(n=2,p=2,T=100,type="const")
res_d = VARData(n=3,p=2,T=200,type="exog1",X=matrix(rnorm(400),200,2))
}
