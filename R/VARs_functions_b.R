### Root2coef
### This is function which maps roots (in L) of the characteristic function of an AR process as inout
### to the AR coefficients
###
### input
###
### p:   nummber of lags
### r_p  (optional) A p-vector of roots outside the unit circle
###
### output
### a_p  p-vector of lag-coefficients
###
### the recursion follows from the equation:    (1-a[3,1]L^1-a[3,2]L^2-a[3,3]L^3)(L-r_4)
###                                           = (1-a[4,1]L^1-a[4,2]L^2-a[4,3]L^3-a[4,4]L^4)r_4
###
###
###
#' Roots2Coefficients
#'
#' Given the roots of a characteristic polynom in lags, output the coefficients of the corresponding AR process.
#' @param p The lag length
#' @param r_p A p-vector of roots outside the unit circle
#'
#' @return A vector of AR coefficients
#' @export
#'
#' @examples
#'
#' Roots2coef(3,c(1.1,1.2,1.3))
Roots2coef = function(p,r_p) {
   if (missing(r_p)) r_p <- 0.5/(stats::runif(p)-0.5)  #random number outside unit circle
   if (min(abs(r_p)) < 1) {
   	return("r_p is within the unit circle")
   	#print("r_p is within the unit circle")
   }
   a = matrix(0,p,p)
   a[1,1] = 1/r_p[1]
   if ( p>1 ) {
   for (i in 2:p) {
      for ( j in 1:i ) {
         if (j == 1)           a[i,j] = a[i-1,j] + 1/r_p[i]
         if ((j > 1) & (j<i))  a[i,j] = a[i-1,j] - a[i-1,j-1]/r_p[i]
         if (j == i)           a[i,j] = -a[i-1,i-1]/r_p[i]
      }
   }
   }

   #R = matrix(0,p,p)
   #for (i in 1: p)     {
   #   for (j in 1: p ) {
   #        R[i,j] = r_p[j]^i
   #   }
   #}
   #a[p,]%*%R[,]
   return(a[p,])
}


#' Multivariate normal random series
#'
#' This function will generate iid multivariate normal random time series.
#'
#' @param T Length of the generated time series
#' @param sigma An (n x n) covariance matrix of the normal series
#' @return T x n matrix of iid normal time series
#' @export
rnormSIGMA = function(T,sigma) {
    # generate random numbers from iid multivariate normal distribution with covariance matrix Sigma
    n = dim(sigma)[1]
    U = stats::rnorm(T*n)
    dim(U) = c(T,n)
    U = U%*%chol(sigma)
    return(U)
}



#' Conditional normal random numbers
#'
#' This function generates random numbers from iid multivariate conditional normal distribution with covariance matrix Sigma, given i-th component has the value of v, this will be an (n-1) dimensional random number
#'
#' @param T Length of generated time series
#' @param sigma The (n x n) covariance matrix of the normal series
#' @param I Index of conditioning component
#' @param v The value of the conditioning component
#' @param switch A switch variable: switch = 1 gives the conditional random series and switch = 0 gives the expected values.
#'
#' @return A (T x (n-1)) matrix of iid conditional normal time series or the conditional expected values
#' @export
rnormSIGMA_cond = function(T,sigma,I,v,switch) {
    # generate random numbers from iid multivariate conditional normal distribution with covariance matrix Sigma, given
    # i-th component has the value of v, this will be an (n-1) dimensional random number
      sigma_cond = as.matrix(sigma[-I,-I])-sigma[-I,I]%*%solve(sigma[I,I])%*%sigma[I,-I]
      mu_cond = sigma[-I,I]%*%solve(sigma[I,I])%*%v
      U = rnormSIGMA(T,sigma_cond)*switch+as.vector(c(1:T)/c(1:T))%*%t(mu_cond)
      return(U)
}


#'  This function selects randomly N elements out of a set of T elements
#'
#' @param N The number of elements to be selected
#' @param T The total number of elements
#' @export
NoutofT = function(N,T) {
  unique(round(stats::runif(3*N)*(T-1)))[1:N]+1
}

#' Impulse response function of a vector autoregressive model
#'
#' This function generates impulse response functions for VAR,CIVAR,MRVAR MRCIVAR, also for GVAR, CIGVAR, MRGVAR and MRCIGVAR.
#' For the later four classes of models it also provides the functionalities to calculate the global, regional and country-specific shocks.
#' It also calculates global and regional responses and coordinated policy actions.
#'
#' @param B An (nxnxp) coefficients array of an n-dimensional VAR(p) model
#' @param sigma The covariance matrix of the VAR(p) model
#' @param nstep Number of steps of the impulse response function
#' @param comb An n-vector of weights of the coordinated policy actions
#' @param irf Type of impulse response function
#' @param G The matrix used in the permanent and transitory decomposition
#' @param smat An explicit decomposition matrix that defines a structural shock.
#'
#' @export
irf_B_sigma = function (B, sigma, nstep, comb, irf = c("gen", "chol", "chol1","gen1","genN1", "comb1","smat","concerts1"),G=NA,smat=NA)
{
    neq <- dim(B)[1]
    nvar <- dim(B)[2]
    lags <- dim(B)[3]
    n = dim(sigma)[1]

    if (irf == "smat") {
        smat = (smat)
    }

    if (irf == "chol") 	   {        smat = chol(sigma)    					    }
    if (irf == "PTdecomp") {        smat = chol(G%*%sigma%*%t(G))          		    }
    if (irf == "gen")      {        smat = t(sigma %*% diag(1/(sqrt(diag(sigma)))))     }
    if (irf == "chol1")    {        smat = chol(sigma) %*% diag(1/diag(chol(sigma)))    }
    if (irf == "gen1")     {        smat = t(sigma %*% diag(1/(diag(sigma))))           }
    if (irf == "genN1")    {        smat = t(sigma %*% diag(-1/(diag(sigma))))          }
    if ( irf == "concerts1") { smat  = t((sigma%*%diag(1/(diag(sigma))))%*%comb)        };
    if ( irf == "concerts0") { smat  = t((sigma%*%diag(1/(sqrt(diag(sigma)))))%*%comb)  };

    if ( irf == "concertc") {
	   c     = as.numeric(!(comb[,1]==0));
         smat  = matrix(0,n,n);
         for (i in 1: n) { smat[i,] = sigma[i,]%*%diag(c)%*%INVI(sigma,c,i)%*%comb; }
         smat = t(smat);
    };

    if (irf == "comb") {
        DD = diag(t(comb) %*% sigma %*% comb)
        for (i in 1:length(DD)) {
            if (DD[i] > 0)
                DD[i] = sqrt(1/DD[i])
        }
        DD = diag(DD)
        smat = t(sigma %*% comb %*% DD)
    }
    if (irf == "comb1") {
        DD = diag(t(comb) %*% sigma %*% comb)
        for (i in 1:length(DD)) {
            if (DD[i] > 0)
                DD[i] = (1/DD[i])
        }
        DD = diag(DD)
        smat = t(sigma %*% comb %*% DD)
    }
    if (dim(smat)[2] != dim(B)[2])
        stop("B and smat conflict on # of variables")
    response <- array(0, dim = c(neq, nvar, nstep))
    response[, , 1] <- t(smat)
    for (it in 2:nstep) {
        for (ilag in 1:min(lags, it - 1)) response[, , it] <- response[,
            , it] + B[, , ilag] %*% response[, , it - ilag]
    }
    dimnames(response) <- list(dimnames(B)[[2]], dimnames(smat)[[1]],
        as.character(0:(nstep - 1)))
    return(response)
}




#' Plot impulse response functions
#'
#' @param IRF_CB An (n x n x L x 3) array of impulse response function with confidence bands
#' @param Names An n-vector of strings of the variable names
#' @param INames An n-vector of string of the impulse names
#' @param response An vector of impulse indices
#' @param impulse An vector of response indices
#' @param ncol Number of columns of impulse response functions in the plot
#'
#' @return An ggplot object of impulse response functions
#' @export
IRF_graph <- function(IRF_CB=IRF_CB,Names = NA,INames=NA,response=c(1:n),impulse=c(1:n),ncol=n) {
  ### This function create a list of ggplot objects for the impulse response functions
  IRF_list = list()
  n = dim(IRF_CB)[1]
  if (anyNA(Names))   Names  <-  colnames(IRF_CB)
  if (is.null(Names)) Names  <-  paste0(rep("Y",n),c(1:n))
  if (anyNA(INames))  INames <-  Names

  k = 0
  for (i in response) for (j in impulse) {
    k =  k +1
    myData <- as.data.frame(cbind(c(0:(dim(IRF_CB)[3]-1)),IRF_CB[i,j,,]))

    IRF_list[[k]] <-  ggplot2::ggplot(myData,
                                      ggplot2::aes(x=V1, y=V2, ymin=V3, ymax=V4)) +
      ggplot2::geom_hline(yintercept = 0, color="red") +
      ggplot2::geom_ribbon(fill="grey", alpha=0.5) +
      ggplot2::geom_line() +
      ggplot2::theme_light() +
      ggplot2::ggtitle(paste(Names[i],"response to", INames[j], "shock"))+
      ggplot2::ylab("")+
      ggplot2::xlab("") +
      ggplot2::theme(plot.title = ggplot2::element_text(size = 11, hjust=0.5), axis.title.y = ggplot2::element_text(size=11))
  }

  do.call(gridExtra::grid.arrange,c(IRF_list,ncol=ncol))
  return(IRF_list)
}


#' Embedding a time series
#'
#' Embeds the time series y into a low-dimensional Euclidean space with column names.
#'
#' @param y The time series to be embedded
#' @param p Number of lags for embedding
#' @param prefix Prefix of the column names
#'
#' @return A matrix containing the embedded time series y.
#' @export
Embed <- function(y=tseries::ts(c(1:20)),p=3,prefix="") {
  YY = stats::embed(y,p)
  if (!is.null(colnames(y))) {
    strc = colnames(y)
    strc = paste(prefix,strc,sep="")
    str0 = strc
    if (p > 1) {
      for (i in 1:(p-1)) {
        stri <- paste(str0, ".", sep="")
        stri <- paste(stri,as.character(i),sep="")
        strc <- cbind(strc,stri)
      }
    }
    colnames(YY) <- as.vector(strc)
  }
  return(YY)
}



#' Help function for nonlinear optimization in constrained estimation of CIVAR
#'
#' @param x Optimization variables
#' @param beta Cointegration vectors
#' @param alpha Adjustment vectors
#' @param G A matrix specifying restrictions on alpha
#' @param H A matrix specifying restrictions on beta
#' @param phi Freely varying parameters in beta
#' @param psi Freely varying parameters in alpha
#' @param h A vector specifying restrictions in beta
#' @param Z1 I(1) data matrix
#' @param St Regime 1 indicator series
#' @param NSt Regime 2 indicator series
#' @param Y0 I(0) data matrix
#' @param Z2 I(0) data matrix
#'
#' @return Sum of squared residuals
#' @export
f_constrained <- function(x, beta=beta,alpha=alpha,G=G,H=H,phi=phi,psi=psi,h=h, Z1, St, NSt, Y0, Z2) {
  ## this is a help function to incoporate restrictions on beta and regime specific alpha
  ## vec(alpha'_1) = G_1 psi_1,  vec(alpha'_2) = G_2 psi_2   , vec(beta) = H phi + h
  ##
  x1 = x[1:length(phi)]
  x2 = x[(1+length(phi)):(length(phi)+length(psi[[1]]))]
  x3 = x[(1+(length(phi)+length(psi[[1]]))):(length(phi)+length(psi[[1]])+length(psi[[2]]))]

  ### restriction on beta
  phi  = x1
  beta_v = H%*% phi + h
  dim(beta_v) = dim(beta)
  #?dim(x) = dim(beta)

  CI = Z1 %*% beta_v
  ### restrictions on alpha_1,alpha_2
  psi_1 = x2
  psi_2 = x3
  alpha_1 = G[[1]]%*%psi_1
  alpha_2 = G[[2]]%*%psi_2
  dim(alpha_1) = dim(alpha)
  dim(alpha_2) = dim(alpha)
  CI1 = (CI * St)%*%t(alpha_1)
  CI2 = (CI * NSt)%*%t(alpha_2)
  residuals <- stats::lm(Y0-CI1-CI2 ~ 0 + Z2)$residuals
  #SSR = sum(diag(t(residuals) %*% residuals))
  SSR <- det(t(residuals) %*% residuals)
  return(SSR)
}


#' Test restrictions within the cointegration space
#'
#' This function runs a likelihood ratio test of linear restrictions on \eqn{\alpha} and \eqn{\beta} in a CIVAR model
#' @param res An object of the output of CIVARest
#' @param H A matrix specifying the restrictions on beta
#' @param h A vector specifying the restrictions on beta
#' @param phi freely varying parameters in beta
#' @param G A matrix specifying the restrictions on alpha
#' @param psi freely varying parameters in alpha
#' @param method Method used to run the test
#'
#' @section Details:
#' This function runs a likelihood ratio test of linear restrictions on alpha and beta in a CIVAR model in the following form:
#'		\deqn{vec(\alpha') = G \psi,  vec(\beta) = H\phi + h}
#'
#' @return A list containing constrained estimated CIVAR model
#' @examples
#' Sigma = 1:(4*4*2)
#' dim(Sigma) = c(4,4,2)
#' Sigma[,,1] = diag(4)
#' Sigma[,,2] = diag(4)
#' p=matrix(0,2,2)
#' p[,1] = c(2,2)
#'
#' res_d = MRCIVARDatam(n=4,p=p,T=261,S=2,SESVI=1,TH=0,Sigmao=Sigma,type="const",r=1)
#' colnames(res_d$Y) = c("w","p","y","r")
#' res_e = MRCIVARestm1(res=res_d)
#' res_e$Summary
#'
#' ### case 1 Compare MRCIVARestm1 with AB_MEVARTestm
#' n = 4; crk = 3
#' G = diag(12); G2 = list(G,G);  psi = matrix(1,12,1); psi2=list(psi,psi); # No restrictions on alpha
#' H = diag(12); H2 = H[,-c(1,5,9)]; h = matrix(0,12,1); h[c(1,5,9),1] = c(1,1,1);  phi = matrix(1,9,1)
#'
#' G2[[1]]%*%psi2[[1]]; G2[[2]]%*%psi2[[2]]; H2%*%phi+h
#'
#' res_t  = MRCIVARTestm(res=res_d,H=H2,h=h,phi=phi,G=G2,psi=psi2,method=c("parametric"))
#' res_eR = ABC_MRCIVARestm(res=res_d,H=H2,h=h,phi=phi,G=G2,psi=psi2)
#'
#' res_eR$LR
#' res_eR$p_value                                 #### p-value = 1 because thre is no effective restrictions.
#' res_eR$code
#' summaryCIVAR(res_eR$VECMR,sname="Z2")
#'
#'
#' ### case 2 test of weak exogeneity of 1st variable
#' n = 4; crk = 3
#'
#' G = diag(12); G2 = list(G[,4:12],G[,4:12]);  psi = matrix(1,9,1); psi2=list(psi,psi);
#'
#'
#' H = diag(12); H2 = H[,-c(1,5,9)]; h = matrix(0,12,1); h[c(1,5,9),1] = c(1,1,1);  phi = matrix(1,9,1)
#'
#'
#'  G2[[1]]%*%psi2[[1]]; G2[[2]]%*%psi2[[2]]; H2%*%phi+h
#'
#'
#' res_t = MRCIVARTestm(res=res_d,H=H2,h=h,phi=phi,G=G2,psi=psi2,method=c("parametric"))
#' res_eR = ABC_MRCIVARestm(res=res_d,H=H2,h=h,phi=phi,G=G2,psi=psi2)
#'
#' res_eR$LR
#' res_eR$p_value                                 #### p-value = 1 because thre is no effective restrictions.
#' res_eR$code
#' summaryCIVAR(res_eR$VECMR,sname="Z2")
#' @export
MRCIVARTestm <- function(res=res,H=H,h=h,phi=phi,G=G,psi=psi,method=c("parametric","bootstrap")) {
  ###
  ### This function runs a likelihood ratio test of linear restrictions on alpha and beta in a CIVAR model
  ###
  ###		vec(alpha') = G psi , vec(beta) = H phi + h
  ###
  ###        example 1 (restrictions on alpha) test of exogeneity
  ###			   vec(alpha) is 45 x 1  vector  ( N = 9 crk = 5, defined by the model )
  ###                  G          is 45 x 40 matrix  ( the first variable is exogeneous, i.e. the 5 adjustment coefficient of the first variable are zero )
  ###		         psi        40 x 1 vector (free variaring parameters not appearing in the specification but implied by G)
  ###                  vec(beta)  is 45 x 1 vector   ( N = 0 crk = 5 )
  ###                  H          is 45 x 45 identity matrix
  ###                  phi        45 x 1 vector (free variaring parameters not appearing in the specification but implied by h =0
  ###                  h          45 x 1 zero matrix implying ver(beta) = phi  >> no restrictions on beta.
  ###			   (H is identity and h is zero vector implies only restrictions on alpha)
  ###
  ###        example 2 (restrictions on beta ) test of PPP
  ###  			   vec(alpha) is 40 x 1  vector  ( N = 8 crk = 5 ) conditioal VECM  or VECMX model
  ###                  G          is 40 x 40 identity matrix, implying there is no restriction on alpha
  ###		         psi        40 x 1 vector (free variaring parameters not appearing in the specification but implied by the identity matrix G
  ###                  vec(beta)  is 45 x 1 vector   ( N = 0 crk = 5, defined by the model )
  ###                  H          is 45 x 2 matrix that picks out the elements under restricitons ( two colunms out of the identity matrix ) a zero row in H and the corresponding h implies zero-restrictions on beta.
  ###                             ones in a row of H and zero in the corresponding h implies non-restricted beta.
  ###                  phi        2  x 1  vector (free variaring parameters not appearing in the specification but implied by h =0
  ###                  h          45 x 1  non zero elements in this vector together with the zero elements in the corresponding row in H are the normalization conditions.
  ###			   (H is identity and h is zero vector implies only restrictions on alpha)
  ###
  ###
  ### (y,x,model = c("I","II","III","IV","V"), bnorm = c("1","2"),type = c("eigen", "trace"),p = 1, r = 2, q = 0.95, H=H,h=h,G=G)
  y = res$Y
  if (is.na(res$X)) x = 0 else x = res$X

  if (res$type == "const") model = "III"
  if (res$type == "none")  model = "I"
  bnorm = "1"
  type  = "eigen"
  p     =  res$p
  r     =  res$crk
  St    = res$St
  St    = St - 1
  crk   = res$crk


  #y = res$res$Y
  #if (is.na(res$res$X)) x = 0 else x = res$res$X
  #if (res$res$type == "const") model = "III"
  #if (res$res$type == "none")  model = "I"
  #bnorm = "1"
  #type  = "eigen"
  #p     =  res$res$p
  #r     =  res$res$crk
  #St    = res$res$St
  #St    = St - 1
  #crk   = res$res$crk
  ### N = 7, crk = 5, alphe beta are N  x crk  = 7 x 5   H = 35x35 G 35x35
  ###  G = diag(35); h = matrix(0,35,1); h[c(1,8,15,22,29)] = c(1,1,1,1,1)  ; H = diag(35); H2 = H[,-c(1,8,15,22,29)]
  test = abMRCVECMest2testm(y=y, x=0, s=St, model = model, type = type, P=p, crk = crk, q = 0.95, Dxflag = 0,H=H,h=h,phi=phi,G=G,psi=psi)

  return(test)
}



#' Estimation of constrained CIVAR models
#'
#' @param y The series of endogenous variables
#' @param x The series of conditioning/exogeneous variables
#' @param s The series of the regime indicator function
#' @param model Types of deterministic components in CIVAR models
#' @param type Type of the Johansen test
#' @param ret  Statistics
#' @param ctable Critical value tables
#' @param crk Cointegration rank
#' @param P Lags of MRCIVAR
#' @param q Significance level
#' @param Dxflag A flag specifying whether the conditioning variablesare within the cointegration space or not
#' @param H A matrix specifying restrictions on beta
#' @param h A vector specifying restrictions on beta
#' @param phi freely varying parameters in beta
#' @param G A matrix specifying restrictions on alpha
#' @param psi freely varying parameters in alpha
#'
#' @return A list of estimated constrained MRCIVAR
#' @export
abMRCVECMest2testm <- function (y, x, s, model = c("I", "II", "III","IV", "V"), type = c("eigen", "trace"), ret = c("statistic", "test"), ctable = c("A3", "A1", "A5"), crk = crk,P = matrix(2, 2, 2), q = 0.95, Dxflag = 0, H=H,h=h,phi=phi,G=G,psi=psi)
{
  y <- as.matrix(y)
  #model <- match.arg(model)
  #type <- match.arg(type)
  #ret <- match.arg(ret)
  #ctable <- match.arg(ctable)
  if (q != 0.9 && q != 0.95 && q != 0.99) {
    print("please correct significance level")
    (break)()
  }
  S = 2
  p <- as.integer(max(P))
  pmin <- as.integer(min(P))
  N1 <- ncol(as.matrix(y))
  NN1 <- ncol(as.matrix(x))
  N <- ncol(as.matrix(y))
  n = N
  if (N1 < crk) {
    print("y's dimension must be larger than crk")
    (break)()
  }
  if (missing(s)) {
    s = NA
  }
  if (missing(Dxflag)) {
    Dxflag = 0
  }
  if (!anyNA(s)) {
    St = s[(p + 1):length(s)]
    NSt = 1 - s[(p + 1):length(s)]
  }
  if (!sum(abs(x)) == 0) {
    z <- cbind(y, x)
    Zy <- stats::embed(diff(y), p)
    Zx <- stats::embed(diff(x), p)
    if (P[1, 1] < p) {
      Aa = P[1, 1] * N1 + (1:((p - P[1, 1]) * N1))
      ZyI <- Zy[, -Aa]
    }
    else ZyI = Zy
    if (P[1, 2] < p) {
      Ba = P[1, 2] * NN1 + (1:((p - P[1, 2]) * NN1))
      ZxI = Zx[, -Ba]
    }
    else ZxI = Zx
    if (P[2, 1] < p) {
      Aa = P[2, 1] * N1 + (1:((p - P[2, 1]) * N1))
      ZyII = Zy[, -Aa]
    }
    else ZyII = Zy
    if (P[2, 2] < p) {
      Bb = P[2, 2] * NN1 + (1:((p - P[2, 2]) * NN1))
      ZxII = Zx[, -Bb]
    }
    else ZxII = Zx
    Z = cbind(ZyI, ZxI)
    Z_2 = cbind(ZyII, ZxII)
    if (Dxflag == 0) {
      Z2 <- Z[, -c(1:N1, P[1, 1] * N1 + 1:NN1)]
      ZS_2 <- Z_2[, -c(1:N1, P[2, 1] * N1 + 1:NN1)]
    }
    else {
      Z2 <- Z[, -c(1:N1)]
      ZS_2 <- Z_2[, -c(1:N1)]
    }
  }
  else {
    z = y
    Z = stats::embed(diff(y), p)
    if (Dxflag == 0) {
      Z2 <- Z[, -c(1:N1)]
      ZS_2 <- Z2
      Z2 <- Z2[, 1:((P[1, 1] - 1) * n)]
      ZS_2 <- ZS_2[, 1:((P[2, 1] - 1) * n)]
    }
    else {
      Z2 <- Z[, -c(1:N1)]
      Zs_2 <- Z2
      Z2 <- Z2[, 1:((P[1, 1] - 1) * n)]
      ZS_2 <- ZS_2[, 1:((P[2, 1] - 1) * n)]
    }
  }
  M1 <- ncol(as.matrix(z))
  T <- nrow(as.matrix(z))
  MM1 = ncol(Z)
  Y0 <- Z[, c(1:N1)]
  Z1 <- z[-T, ][p:(T - 1), ]
  if (!anyNA(s))
    Z2 = cbind(St * Z2, NSt * ZS_2)
  T1 <- nrow(as.matrix(Y0))
  lT = (1:T1)/(1:T1)
  Trend <- matrix(1:T1, T1, 1)
  if (model == "I") {
    Y0 = Y0
    Z1 = Z1
    Z2 = Z2
  }
  if (model == "II") {
    Y0 = Y0
    Z1 = cbind(lT, Z1)
    Z2 = Z2
  }
  if (model == "III") {
    Y0 = Y0
    Z1 = Z1
    if (!anyNA(s))
      Z2 = cbind(Z2, St, NSt)
    else Z2 = cbind(Z2, lT)
  }
  if (model == "IV") {
    Y0 = Y0
    Z1 = cbind(Trend, Z1)
    Z2 = cbind(Z2, lT)
  }
  if (model == "V") {
    Y0 = Y0
    Z1 = cbind(Z1)
    Z2 = cbind(Z2, lT, Trend)
  }
  M00 <- crossprod(Y0)/T1
  M11 <- crossprod(Z1)/T1
  M22 <- crossprod(Z2)/T1
  M01 <- crossprod(Y0, Z1)/T1
  M02 <- crossprod(Y0, Z2)/T1
  M10 <- crossprod(Z1, Y0)/T1
  M20 <- crossprod(Z2, Y0)/T1
  M12 <- crossprod(Z1, Z2)/T1
  M21 <- crossprod(Z2, Z1)/T1
  M22inv <- solve(M22)
  R0 <- Y0 - t(M02 %*% M22inv %*% t(Z2))
  R1 <- Z1 - t(M12 %*% M22inv %*% t(Z2))
  S00 <- crossprod(R0)/T1
  S01 <- crossprod(R0, R1)/T1
  S10 <- crossprod(R1, R0)/T1
  S11 <- crossprod(R1)/T1
  Ctemp <- chol(S11, pivot = TRUE)
  pivot <- attr(Ctemp, "pivot")
  oo <- order(pivot)
  C <- t(Ctemp[, oo])
  Cinv <- solve(C)
  S00inv <- solve(S00)
  valeigen <- eigen(Cinv %*% S10 %*% S00inv %*% S01 %*% t(Cinv))
  lambda <- valeigen$values
  e <- valeigen$vector
  V <- t(Cinv) %*% e
  Vorg <- V
  V <- sapply(1:M1, function(j) V[, j]/V[1, j])
  W <- S01 %*% V %*% solve(t(V) %*% S11 %*% V)
  PI <- S01 %*% solve(S11)
  DELTA <- S00 - S01 %*% V %*% solve(t(V) %*% S11 %*% V) %*%
    t(V) %*% S10
  GAMMA <- M02 %*% M22inv - PI %*% M12 %*% M22inv
  beta <- as.matrix(V[, 1:crk])
  beta0 <- as.matrix(V[, 1:crk])

  if ((crk>0)&(length(Z2)>0))  {
    CI = Z1 %*% beta
    VECM <- stats::lm(Y0 ~ 0 + CI + Z2)
    if (crk==1 )  { alpha = as.matrix(VECM$coefficients[1:1,]) }
    if (crk > 1)  { alpha = t(as.matrix(VECM$coefficients[1:crk,])) }
    LSKOEF = stats::lm(Y0-Z1%*%beta%*%t(alpha)~0+Z2)$coefficients
  }


  if (crk > 0) {
    betaS = stats::nlm(f, as.vector(beta0), beta, Z1, St, NSt, Y0, Z2)$estimate
    dim(betaS) = dim(as.matrix(beta))
    for (i in 1:ncol(as.matrix(betaS))) betaS[, i] = betaS[,i]/(betaS[1, i])
    CI = Z1 %*% betaS
    estimation <- stats::lm(Y0 ~ 0 + CI + Z2)
    CI1 = CI * St
    CI2 = CI * NSt
    VECM1 <- stats::lm(Y0 ~ 0 + CI1 + CI2 + Z2)
    alphaS_1 = t(VECM1$coefficients[1:crk,])
    alphaS_2 = t(VECM1$coefficients[(1+crk):(2*crk),])
    if (dim(alphaS_1)[1]==1) alphaS_1 <-t(alphaS_1)
    if (dim(alphaS_2)[1]==1) alphaS_2 <-t(alphaS_2)

  }
  if (crk == 0) {
    CI = 0
    estimation <- stats::lm(Y0 ~ 0 + Z2)
  }


  E = -T1 * log(1 - lambda)
  E = E[1:N]
  resultsvecm <- summary(VECM)
  if (model == "I") {
    Tab <- Tab1
  }
  if (model == "II") {
    Tab <- Tab2
  }
  if (model == "III") {
    Tab <- Tab3
  }
  if (model == "IV") {
    Tab <- Tab4
  }
  if (model == "V") {
    Tab <- Tab5
  }
  b = c(1:12)
  for (i in 1:12) {
    b[i] = Tab[2 * (i - 1) + 1, 2 + M1 - N1]
  }
  a = c(1:12)
  for (i in 1:12) {
    a[i] = Tab[2 * i, 2 + M1 - N1]
  }



  if (type == "eigen") {
    critical_vals = b
    M = matrix(0, N, 1)
    j = 1
    rank = 0
    while (j <= N && E[j] > critical_vals[N + 1 - j]) {
      M[j, ] = M[j, ] + 1
      j = j + 1
      rank = rank + 1
    }
    #if (ret == "test") {
    #    return(M)
    #}
    erg <- cbind(E, critical_vals[N:1])
    colnames(erg) <- c("teststatistic", "critical_value")
    if (N > 1) {
      rownames(erg) <- c("crk <= 0 |", paste("crk <= ",
                                             1:(N - 1), " |", sep = ""))
    }
    if (N == 1) {
      rownames(erg) <- c("crk <= 0 |")
    }
    coint_rank <- paste("Johansen-Test (with maximum-eigenvalue-teststatistic) indicates",
                        rank, "cointegrating equation(s) at the", 1 -
                          q, "level")

    LOGL = -T1/2.0*(log(det(S00))+sum(log((1-lambda)[1:crk])))-T1*N1/2.0*log(2*pi)-T1*N1/2.0

    Omega = 1.0/(nrow(VECM$residuals)-(p-1)*N1-crk)*t(VECM$residuals)%*%VECM$residuals



    ######### Calculation of restricted MLE for non-regime-specific cointegrated system as initial value
    alphai = alpha
    betai = beta
    Omegai = Omega
    VS10 = S10
    dim(VS10) <- c(length(S10),1)
    error = 1
    tol   = 0.000001
    vecpip = solve(S11)%*%S10;dim(vecpip) = c(length(vecpip),1)

    phii = solve(t(H)%*%kronecker(t(alphai)%*%solve(Omegai)%*%alphai,S11)%*%H)%*%
      (t(H)%*%kronecker(t(alphai)%*%solve(Omegai),diag(nrow(beta)))%*%VS10
       - t(H)%*%kronecker(t(alphai)%*%solve(Omegai)%*%alphai,S11)%*%h)

    betai <-H%*%phii+h ; dim(betai) <- c(nrow(beta),crk)

    gammai = solve(t(G[[1]])%*%kronecker(solve(Omegai),t(betai)%*%S11%*%betai)%*%G[[1]])%*%
      t(G[[1]])%*%kronecker(solve(Omegai),t(betai))%*%VS10
    talphai <- G[[1]]%*%gammai;dim(talphai) <- c(crk,N1)
    alphai = t(talphai)
    Omegai = S00 - alphai%*%t(betai)%*%S10-S01%*%betai%*%t(alphai)+alphai%*%t(betai)%*%S11%*%betai%*%t(alphai)

    ##################################
    ##################################


    ################################
    ################################

    while  ( error > tol ) {

      phir   = phii
      gammar = gammai
      Omegar = Omegai

      phii = solve(t(H)%*%kronecker(t(alphai)%*%solve(Omegai)%*%alphai,S11)%*%H)%*%
        (t(H)%*%kronecker(t(alphai)%*%solve(Omegai),diag(nrow(beta)))%*%VS10
         - t(H)%*%kronecker(t(alphai)%*%solve(Omegai)%*%alphai,S11)%*%h)

      betai <-H%*%phii+h ; dim(betai) <- c(nrow(beta),crk)

      gammai = solve(t(G[[1]])%*%kronecker(solve(Omegai),t(betai)%*%S11%*%betai)%*%G[[1]])%*%t(G[[1]])%*%kronecker(solve(Omegai),t(betai))%*%VS10

      talphai <- G[[1]]%*%gammai;dim(talphai) <- c(crk,N1)
      alphai = t(talphai)
      Omegai = S00 - alphai%*%t(betai)%*%S10-S01%*%betai%*%t(alphai)+alphai%*%t(betai)%*%S11%*%betai%*%t(alphai)
      error = max( abs(log(det(Omegai))-log(det(Omegar))) )
      log(det(Omegai))
      ### literature:  Identifying, Estimating and Testing Restricted CointegratedSystems: An Overview H. Peter Boswijk Jurgen A. Doornik 2003

    }

    alphar = alphai
    betar  = betai
    Omegar = Omegai

    tst = AB_MRCIVARTest(R0,R1,G=G,H=H,h=h,alphaR1=alphaS_1,alphaR2=alphaS_2,betaR=betaS,alpha1=alphaS_1,alpha2=alphaS_2,beta1=betaS,beta2=betaS,OmegaR=Omega,Omega1=Omega,Omega2=Omega,IC=1,T1=T1,St=St,NSt=NSt)

    #### regime-specific constrained  ML

    #### vec(alpha'_1) = G_1 psi_1,  vec(alpha'_2) = G_2 psi_2   , vec(beta) = H phi + h

    #phi = as.vector(betaS)[! (as.vector(betai) - h)==0]
    #Psi0 = as.vector(t(alphai))[ !as.vector(t(alphai))==0]
    #psi[[1]] = as.vector(t(alphai))[ !as.vector(t(alphai))==0]
    #psi[[2]] = (as.vector(t(alphai))[ !as.vector(t(alphai))==0])[1:30]

    phio = as.vector(betaS)[(H%*%phi-h)==1]
    psi_1 = as.vector(t(alphaS_1))[!G[[1]]%*%psi[[1]]==0]
    psi_2 = as.vector(t(alphaS_2))[!G[[2]]%*%psi[[2]]==0]



    x = c(t(as.vector(phio)),t(as.vector(psi_1)),t(as.vector(psi_2)))
    #f_constrained(x=c(phi,psi[[1]],psi[[2]]), beta=beta,alpha=alpha,G=G,H=H,phi=phi,psi=psi,h=h, Z1, St, NSt, Y0, Z2)
    if (crk > 0) {
      x       = c(phio,psi_1,psi_2)
      XX      = stats::nlm(f_constrained,x, beta, alpha, G, H, phi, psi,h, Z1, St, NSt, Y0, Z2,iterlim=100)
      xR      = XX$estimate
      #### not converge for bad initial value but converge for good initial value
      phir    = xR[1:length(phi)]
      psi_1r  = xR[(1+length(phi)):(length(phi)+length(psi[[1]]))]
      psi_2r  = xR[(1+(length(phi)+length(psi[[1]]))):(length(phi)+length(psi[[1]])+length(psi[[1]]))]
      betaR   = H%*%phir + h
      alpha_1 = G[[1]]%*%psi_1r
      alpha_2 = G[[2]]%*%psi_2r

      dim(betaR) = dim(beta)
      CI = Z1 %*% betaR
      ### restrictions on alpha_1,alpha_2
      dim(alpha_1) = dim(t(alpha))
      alpha_1 = t(alpha_1)
      dim(alpha_2) = dim(t(alpha))
      alpha_2 = t(alpha_2)

      alphaR  = list(alpha_1,alpha_2)
      CI1 = (CI * St) %*%t(alpha_1)
      CI2 = (CI * NSt)%*%t(alpha_2)
      LM  <- stats::lm(Y0-CI1-CI2 ~ 0 + Z2)
      residuals <- LM$residuals
      CI10 = CI * St
      CI20 = CI * NSt
      VECMR <- stats::lm(Y0 ~ 0 + CI10 + CI20 + Z2)
      VECMR$residuals <-  residuals
      #dim(VECMR$coefficients)
      VECMR$coefficients[1:crk,]  = t(alpha_1)
      VECMR$coefficients[1:(2*crk),] = t(alpha_2)
      VECMR$coefficients[(2*crk+1):nrow(VECMR$coefficients),] = LM$coefficients

    }

    ##LR, likelihood Ratio between unrestricted model crk = n, vs crk = crk
    ##LR1 likelihood ratio between unrestricted model crk = n, vs crk = crk with restricions beta and alpha
    ##LR2 likelihood ratio between unrestricted alpha and beta vs restrictions on alpha and beta
    Omega1 = t(VECM1$residuals)%*%(VECM1$residuals)/(nrow(VECM1$residuals)-(p-1)*N1-crk)
    OmegaR = t(VECMR$residuals)%*%(VECMR$residuals)/(nrow(VECMR$residuals)-(p-1)*N1-crk)

    LR  =  T1*(log(det(Omega1))-log(det(S00))-sum(log(1-lambda)[1:crk]))
    LR1 =  T1*(log(det(OmegaR)) -log(det(S00))-sum(log(1-lambda)[1:crk]))
    LR2 =  T1*(log(det(OmegaR))-log(det(Omega1)))

    result <- new.env()
    result$erg = erg
    result$estimation = VECMR
    result$lambda     = E
    result$z = z
    result$Z2 = Z2
    result$beta = betaS
    result$alpha = list(alphaS_1,alphaS_2)
    result$PI = PI
    result$GAMMA = GAMMA
    result$model = model
    result$P = P
    result$NN1 = NN1
    result$s = s
    result$LR = LR
    result$betar  = betaR
    result$alphar = alphaR
    result$tst    = tst
    result$LR1    = LR1
    result$LR2    = LR2
    rr <- as.list(result)
    return(rr)
  }
  else {
    type = "trace"
    critical_vals = a
    stat = matrix(0, N, 1)
    for (i in 1:N) {
      sum = 0
      for (j in i:N) {
        sum = sum + E[j]
      }
      stat[i] = sum
    }
    M = matrix(0, N, 1)
    j = 1
    rank = 0
    while (stat[j] > critical_vals[N + 1 - j] && j <= N) {
      M[j, ] = M[j, ] + 1
      j = j + 1
      rank = rank + 1
    }
    if (ret == "test") {
      return(M)
    }
    erg <- cbind(stat, critical_vals[N:1])
    colnames(erg) <- c("teststatistic", "critical_value")
    if (N > 1) {
      rownames(erg) <- c("crk <= 0 |", paste("crk <= ",
                                             1:(N - 1), " |", sep = ""))
    }
    if (N == 1) {
      rownames(erg) <- c("crk <= 0 |")
    }
    coint_rank <- paste("Johansen-Test (with trace-teststatistic) indicates",
                        rank, "cointegrating equation(s) at the", 1 -
                          q, "level")

    ############ added on 12/12/21

    LOGL = -T1/2.0*(log(det(S00))+sum(log((1-lambda)[1:crk])))-T1*N1/2.0*log(2*pi)-T1*N1/2.0

    Omega = 1.0/(nrow(VECM$residuals)-(p-1)*N1-crk)*t(VECM$residuals)%*%VECM$residuals



    ######### Calculation of restricted MLE for cointegrated system
    alphai = alpha
    betai = beta
    Omegai = Omega
    VS10 = S10
    dim(VS10) <- c(length(S10),1)
    error = 1
    tol   = 0.0001
    vecpip = solve(S11)%*%S10;dim(vecpip) = c(length(vecpip),1)

    phii = solve(t(H)%*%kronecker(t(alphai)%*%solve(Omegai)%*%alphai,S11)%*%H)%*%
      (t(H)%*%kronecker(t(alphai)%*%solve(Omegai),diag(nrow(beta)))%*%VS10
       - t(H)%*%kronecker(t(alphai)%*%solve(Omegai)%*%alphai,S11)%*%h)

    betai <-H%*%phii+h ; dim(betai) <- c(nrow(beta),crk)

    gammai = solve(t(G)%*%kronecker(solve(Omegai),t(betai)%*%S11%*%betai)%*%G)%*%
      t(G)%*%kronecker(solve(Omegai),t(betai))%*%VS10

    talphai <- G%*%gammai;dim(talphai) <- c(crk,N1)
    alphai = t(talphai)

    Omegai = S00 - alphai%*%t(betai)%*%S10-S01%*%betai%*%t(alphai)+alphai%*%t(betai)%*%S11%*%betai%*%t(alphai)


    while  ( error > tol ) {

      phir   = phii
      gammar = gammai
      Omegar = Omegai

      phii = solve(t(H)%*%kronecker(t(alphai)%*%solve(Omegai)%*%alphai,S11)%*%H)%*%
        (t(H)%*%kronecker(t(alphai)%*%solve(Omegai),diag(nrow(beta)))%*%VS10
         - t(H)%*%kronecker(t(alphai)%*%solve(Omegai)%*%alphai,S11)%*%h)

      betai <-H%*%phii+h ; dim(betai) <- c(nrow(beta),crk)

      gammai = solve(t(G)%*%kronecker(solve(Omegai),t(betai)%*%S11%*%betai)%*%G)%*%
        t(G)%*%kronecker(solve(Omegai),t(betai))%*%VS10

      talphai <- G%*%gammai;dim(talphai) <- c(crk,N1)
      alphai = t(talphai)

      Omegai = S00 - alphai%*%t(betai)%*%S10-S01%*%betai%*%t(alphai)+alphai%*%t(betai)%*%S11%*%betai%*%t(alphai)
      error = max( abs(log(det(Omegai))-log(det(Omegar))) )


      log(det(Omegai))
    }
    LR = T1*(log(det(Omegar))-log(det(S00))-sum(log(1-lambda)[1:crk]))

    result <- new.env()
    result$erg = erg
    result$estimation = VECMR
    result$lambda     = E
    result$z = z
    result$Z2 = Z2
    result$beta = betaS
    result$alpha = list(alphaS_1,alphaS_2)
    result$PI = PI
    result$GAMMA = GAMMA
    result$model = model
    result$P = P
    result$NN1 = NN1
    result$s = s
    result$LR = LR
    result$betar  = betaR
    result$alphar = alphaR
    result$tst    = tst
    result$LR1    = LR1
    result$LR2    = LR2
    rr <- as.list(result)
    return(rr)

  }
}



#' Likelihood ratio test of restrictions in CIVAR models
#'
#' This function runs Doornik and Boswijk Iteration to estimate restricted alpha and beta in a CIVAR model.
#'
#' @param R0 Contorled I(0) data matrix
#' @param R1 Controled I(1) data matrix
#' @param G A matrix specifying restrictions on alpha
#' @param H A matrix specifying restrictions on beta
#' @param h A vector specifying restrictions on alpha
#' @param alphaR1 Initial value for constrined estimation
#' @param alphaR2 Initial value for constrined estimation
#' @param betaR Initial value for constrined estimation
#' @param alpha1 Initial value for unconstrined estimation
#' @param alpha2 Initial value for unconstrined estimation
#' @param beta1 Initial value for unconstrined estimation
#' @param beta2 Initial value for unconstrined estimation
#' @param OmegaR Initial value for constrined estimation
#' @param Omega1 Initial value for unconstrined estimation
#' @param Omega2 Initial value for unconstrined estimation
#' @param IC Indicator of BB dimension
#' @param T1 Number of  used observations
#' @param St Regime 1 indicator series
#' @param NSt Regime 2 indicator series
#'
#' @return A list containing the constrained estimates and likelihood test results
#' @export
AB_MRCIVARTest <- function(R0,R1,G,H,h,alphaR1,alphaR2,betaR,alpha1,alpha2,beta1,beta2,OmegaR,Omega1,Omega2,IC=1,T1,St,NSt) {
  ### this function runs DoornikBoswijk Iteration to estimate restricted alpha and beta
  ### It also estimates the separate MRCIVAR and compare it with the restricted including the model with identical beta.
  ### Add output of of separate Models.

  n   = dim(alphaR1)[1]
  crk = dim(alphaR1)[2]

  T1 = nrow(R1)

  #R0 <- Y0 - t(M02 %*% M22inv %*% t(Z2))
  #R1 <- Z1 - t(M12 %*% M22inv %*% t(Z2))

  R0_1 <- R0*St
  R0_2 <- R0*NSt
  R1_1 <- R1*St
  R1_2 <- R1*NSt


  R0 = cbind(R0_1,R0_2)
  R1 = cbind(R1_1,R1_2)

  S00 <- crossprod(R0)%*%diag(c(rep(1/sum(St),n),rep(1/sum(NSt),n)))
  S01 <- crossprod(R0, R1)%*%diag(c(rep(1/sum(St),n),rep(1/sum(NSt),n)))
  S10 <- crossprod(R1, R0)%*%diag(c(rep(1/sum(St),n),rep(1/sum(NSt),n)))
  S11 <- crossprod(R1)%*%diag(c(rep(1/sum(St),n),rep(1/sum(NSt),n)))


  alphai = kronecker(diag(2),alphaR1); alphai[(n+1):(n+n),(crk+1):(crk+crk)] = alphaR2
  betai  = kronecker(diag(2),betaR)
  Omegai = kronecker(diag(2),OmegaR)

  dim(alphai)
  dim(betai)
  dim(Omegai)
  ##
  A = diag(length(alphai))
  G0 =  A[,which(!as.vector(t(alphai))==0)]
  dim(G0)
  #### from G 2 G
  AA = matrix(0,2*dim(G[[1]])[1]+2*dim(G[[2]])[1],dim(G[[1]])[2]+dim(G[[2]])[2])
  GG = AA
  for (i in 1:n ) {
    GG[((i-1)*2*crk+1):((i-1)*2*crk+2*crk),1:dim(G[[1]])[2]]     = kronecker(c(1,0),G[[1]][ ((i-1)*crk+1): ((i-1)*crk+crk),] )
    GG[((i-1+n)*2*crk+1):((i-1+n)*2*crk+2*crk),(1+dim(G[[1]])[2]):(dim(G[[1]])[2]+dim(G[[2]])[2])] = kronecker(c(0,1),G[[2]][ ((i-1)*crk+1): ((i-1)*crk+crk),] )
  }

  if (!IC==1) BB = matrix(0,length(betai),dim(H)[2]+dim(H)[2])
  if ( IC==1) BB = matrix(0,length(betai),dim(H)[2])
  HH = BB
  hh = matrix(0,dim(HH)[1],1)

  for (i in 1:crk ) {
    if (!IC==1) {
      HH[((i-1)*2*n+1):((i-1)*2*n+2*n),1:dim(H)[2]]     = kronecker(c(1,0),H[ ((i-1)*n+1): ((i-1)*n+n),] )
      HH[((i-1+crk)*2*n+1):((i-1+crk)*2*n+2*n),(1+dim(H)[2]):(dim(H)[2]+dim(H)[2])] = kronecker(c(0,1),H[ ((i-1)*n+1): ((i-1)*n+n),] )

      hh[((i-1)*2*n+1):((i-1)*2*n+2*n),1]             = kronecker(c(1,0),h[ ((i-1)*n+1): ((i-1)*n+n),1] )
      hh[((i-1+crk)*2*n+1):((i-1+crk)*2*n+2*n),1]     = kronecker(c(0,1),h[ ((i-1)*n+1): ((i-1)*n+n),1] )
    }
    if (IC==1) {
      HH[((i-1)*2*n+1):((i-1)*2*n+2*n),1:dim(H)[2]]             = kronecker(c(1,0),H[ ((i-1)*n+1): ((i-1)*n+n),] )
      HH[((i-1+crk)*2*n+1):((i-1+crk)*2*n+2*n),1:dim(H)[2]] = kronecker(c(0,1),H[ ((i-1)*n+1): ((i-1)*n+n),] )

      hh[((i-1)*2*n+1):((i-1)*2*n+2*n),1]                     = kronecker(c(1,0),h[ ((i-1)*n+1): ((i-1)*n+n),1] )
      hh[((i-1+crk)*2*n+1):((i-1+crk)*2*n+2*n),1]             = kronecker(c(0,1),h[ ((i-1)*n+1): ((i-1)*n+n),1] )
    }
  }

  HH0 = A[,-which((as.vector(betai)==0)|(as.vector(betai)==1)  )]
  dim(HH0)[2]
  H0 <- HH0[,1: (dim(HH0)[2]/2)]
  H0[(dim(HH0)[1]/2+1):dim(HH0)[1],] = HH0[(dim(HH0)[1]/2+1):dim(HH0)[1],(dim(HH0)[2]/2+1):dim(HH0)[2]]
  #H0 - HH
  h0 = matrix(0,length(betai),1)
  h0[ which(as.vector(betai)==1) ,] = 1
  # h0-hh
  #H = A[,-which((as.vector(betai)==0)|(as.vector(betai)==1)  )] - HH
  #hhh = matrix(0,length(betai),1)
  #hhh[ which(as.vector(betai)==1) ,] = 1

  VS10 = S10
  dim(VS10) <- c(length(S10),1)
  error = 1
  tol   = 0.000000001
  vecpip = solve(S11)%*%S10;dim(vecpip) = c(length(vecpip),1)

  phii = solve(t(HH)%*%kronecker(t(alphai)%*%solve(Omegai)%*%alphai,S11)%*%HH)%*%(t(HH)%*%kronecker(t(alphai)%*%solve(Omegai),diag(nrow(betai)))%*%VS10
                                                                                  - t(HH)%*%kronecker(t(alphai)%*%solve(Omegai)%*%alphai,S11)%*%hh)

  betai <-HH%*%phii+hh ; dim(betai) <- dim(kronecker(diag(2),betaR))


  gammai = solve(t(GG)%*%kronecker(solve(Omegai),t(betai)%*%S11%*%betai)%*%GG)%*% t(GG)%*%kronecker(solve(Omegai),t(betai))%*%VS10
  talphai <- GG%*%gammai;  dim(talphai) <- dim(t(kronecker(diag(2),alphaR1)))
  alphai = t(talphai)
  Omegai = S00 - alphai%*%t(betai)%*%S10-S01%*%betai%*%t(alphai)+alphai%*%t(betai)%*%S11%*%betai%*%t(alphai)


  error = 1
  ######### restriction on identical beta in two regimes: beta_1 = betat_2
  while  ( error > tol ) {
    betar  = betai
    alphar = alphai
    phir   = phii
    gammar = gammai
    Omegar = Omegai

    phii = solve(t(HH)%*%kronecker(t(alphai)%*%solve(Omegai)%*%alphai,S11)%*%HH)%*%
      (t(HH)%*%kronecker(t(alphai)%*%solve(Omegai),diag(nrow(betai)))%*%VS10
       - t(HH)%*%kronecker(t(alphai)%*%solve(Omegai)%*%alphai,S11)%*%hh)

    betai <-HH%*%phii+hh ; dim(betai) <- dim(kronecker(diag(2),betaR))


    gammai = solve(t(GG)%*%kronecker(solve(Omegai),t(betai)%*%S11%*%betai)%*%GG)%*%t(GG)%*%kronecker(solve(Omegai),t(betai))%*%VS10

    talphai <- GG%*%gammai;  dim(talphai) <- dim(t(kronecker(diag(2),alphaR1)))

    alphai = t(talphai)
    Omegai = S00 - alphai%*%t(betai)%*%S10-S01%*%betai%*%t(alphai)+alphai%*%t(betai)%*%S11%*%betai%*%t(alphai)
    error = max( abs(log(det(Omegai))-log(det(Omegar))) )
    #log(det(Omegai))
    ### literature:  Identifying, Estimating and Testing Restricted CointegratedSystems: An Overview H. Peter Boswijk Jurgen A. Doornik 2003
  }

  #### output: uncontrained estimation under identical beta
  error = 1

  #### initial value for uncontrained iteration
  betai[1:n,1:crk]  = beta1;  betai[(n+1):(2*n),(crk+1):(2*crk)]   = beta2;
  alphai[1:n,1:crk] = alpha1; alphai[(n+1):(2*n),(crk+1):(2*crk)]  = alpha2;
  Omegai[1:n,1:n]   = Omega1; Omegai[(n+1):(n+n),(n+1):(n+n)]      = Omega2;
  phii   = solve(t(HH0)%*%kronecker(t(alphai)%*%solve(Omegai)%*%alphai,S11)%*%HH0)%*%
    (t(HH0)%*%kronecker(t(alphai)%*%solve(Omegai),diag(nrow(betai)))%*%VS10
     - t(HH0)%*%kronecker(t(alphai)%*%solve(Omegai)%*%alphai,S11)%*%h0)

  gammai = solve(t(G0)%*%kronecker(solve(Omegai),t(betai)%*%S11%*%betai)%*%G0)%*%t(G0)%*%kronecker(solve(Omegai),t(betai))%*%VS10


  error = 1

  while  ( error > tol ) {
    beta0  = betai
    alpha0 = alphai
    phi0   = phii
    gamma0 = gammai
    Omega0 = Omegai

    phii = solve(t(HH0)%*%kronecker(t(alphai)%*%solve(Omegai)%*%alphai,S11)%*%HH0)%*%
      (t(HH0)%*%kronecker(t(alphai)%*%solve(Omegai),diag(nrow(betai)))%*%VS10
       - t(HH0)%*%kronecker(t(alphai)%*%solve(Omegai)%*%alphai,S11)%*%h0)

    betai <-HH0%*%phii+h0 ; dim(betai) <- dim(kronecker(diag(2),betaR))


    gammai = solve(t(G0)%*%kronecker(solve(Omegai),t(betai)%*%S11%*%betai)%*%G0)%*%t(G0)%*%kronecker(solve(Omegai),t(betai))%*%VS10

    talphai <- G0%*%gammai;  dim(talphai) <- dim(t(kronecker(diag(2),alphaR1)))

    alphai = t(talphai)
    Omegai = S00 - alphai%*%t(betai)%*%S10-S01%*%betai%*%t(alphai)+alphai%*%t(betai)%*%S11%*%betai%*%t(alphai)
    error = max( abs(log(det(Omegai))-log(det(Omega0))) )
    #O0 = diag(c(rep(1,4),rep(0,4)))%*%Omega0%*%diag(c(rep(1,4),rep(0,4)))+diag(c(rep(0,4),rep(1,4)))%*%Omega0%*%diag(c(rep(0,4),rep(1,4)))
    #Oi = diag(c(rep(1,4),rep(0,4)))%*%Omegai%*%diag(c(rep(1,4),rep(0,4)))+diag(c(rep(0,4),rep(1,4)))%*%Omegai%*%diag(c(rep(0,4),rep(1,4)))
    #error = max( abs(log(det(Oi))-log(det(O0))) )

    #log(det(Omegai))
    ### literature:  Identifying, Estimating and Testing Restricted CointegratedSystems: An Overview H. Peter Boswijk Jurgen A. Doornik 2003
  }

  LR = T*(log(det(Omegar))-log(det(Omega0)))
  p_value <- 1- stats::pchisq(LR,dim(G0)[2]-dim(GG)[2]+dim(HH0)[2]-dim(HH)[2])
  ret = list(alphar,betar,Omegar,alpha0,beta0,Omega0,LR,p_value)
  names(ret) = c("alphar","betar","Omegar","alpha0","beta0","Omega0","LR","p_value")
  return(ret)
}



#' Inverse of a partial covariance matrix
#'
#' @param sigma The input covariance matrix
#' @param c The weighting vector of a concerted policy action
#' @param i The index of responding variables
#'
#' @return A response matrix
#' @export
INVI = function(sigma,c,i) {
  n = dim(sigma)[1]
  sigmaout = matrix(0,n,n)
  cc = as.numeric(c>0)*(1:n)
  sigmai = sigma
  for (j in 1:n)         {
    for (k in 1:n)  {
      if ((j==i)|(k==i)) sigmai[j,k] = 0
    }
  }
  sigmai[i,i] = sigma[i,i]
  invsigmai = solve(sigmai[cc,cc])
  sigmaout[cc,cc] = invsigmai
  return(sigmaout)
}

