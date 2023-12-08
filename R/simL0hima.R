#' Simulation data generator for L0-regularized high-dimensional mediation analysis
#'
#'\code{simL0hima} is used to generate simulation data for L0-regularized high-dimensional mediation analysis.
#'
#' @param n An integer specifying sample size.
#' @param p An integer specifying the dimension of mediators.
#' @param q An integer specifying the dimension of covariates.
#' @param family One of the following models:
#' \code{"gaussian"} (continuous outcome),
#' \code{"binomial"} (binomial outcome),
#' \code{"poisson"} (count outcome),
#' \code{"cox"} (right-censored outcome),
#' \code{"AFT"} (right-censored outcome).
#' @param A An vector of the regression coefficients A \code{X ~ M}.
#' @param B An vector of the regression coefficients B \code{M ~ Y}.
#' @param rhoM An coefficient between 0 and 1, which determines the correlation matrix of random errors for M.
#' @param seed An integer specifying a seed for random number generation.
#'
#' @return A list containing the information about the generated data
#' \itemize{
#'     \item{X: }{A vector of predictor.}
#'     \item{M: }{A matrix of mediators. Rows represent samples, columns represent variables.}
#'     \item{Y: }{The outcome variable.}
#'     \item{Z: }{A matrix of covariates.}
#'     \item{A: }{An vector of the regression coefficients A \code{X ~ M}.}
#'     \item{B: }{An vector of the regression coefficients B \code{M ~ Y}.}
#' }
#' @importFrom stats binomial coef complete.cases ks.test glm lm lsfit p.adjust pnorm rnorm runif sd
#' @examples
#' n <- 400  # sample size
#' p <- 1000 # the dimension of mediators
#' q <- 2 # the number of covariates
#' family <- "gaussian"  # the type of model
#'
#' # the regression coefficients A (X --> M)
#' A<-rep(0,p)
#'
#' # the regression coefficients B (M --> Y)
#' B<-rep(0,p)
#'
#' # the first four markers are true mediators.
#' A[1:6] <- c(0.45, 0.4, -0.35, -0.3, 0.0, 0.5)
#' B[1:6] <- c(0.5, -0.45, 0.4, -0.35, 0.5,0.0)
#' rhoM<-0.0 #the correlation coefficient of random errors for M
#' # Generate simulation data
#' dat = simL0hima(n,p,q,family,A,B,rhoM,seed=1234)
#' @export
simL0hima<-function(n,p,q,family=c("gaussian","binomial","poisson", "cox","AFT"),A,B,rhoM=0.0,seed=1e3){
  set.seed(seed)
  C<-0.5

  if(family=="gaussian"){
    X<-rnorm(n,0,sqrt(1))
    #starting value of 0 for each variable: M, Y
    M<-matrix(0,n,p);Y<-rep(0,n)
    var_M=2;var_Y=2

    if(q==0) { Z=NULL;D=NULL;E=NULL
    } else{
      Z<-matrix(rnorm(n*q,0,sqrt(1.0)),n,q)
      D<-matrix(0,q,p);E<-rep(0,q)
      for(i in 1:q){ E[i]<-0.5
      for(j in 1:p){ D[i,j]<-0.5}}
    }

    #Generate error
    mean_M_error <- rep(0,p)
    sigma_M_error <- matrix(0,nrow=p,ncol=p)
    for(i in 1:p){
      for(j in 1:p){
        sigma_M_error[i,j]<-rhoM^(abs(i-j))
      }
    }
    sigma_M_error<-var_M*sigma_M_error;E_M<-mvtnorm::rmvnorm(n,mean_M_error,sigma_M_error)
    E_Y<-rnorm(n,0,sqrt(var_Y)) #for Y

    if(q==0){
      M<-X%*%t(A)+E_M
      Y<-M%*%B+X*C+E_Y
    } else{
      M<-X%*%t(A)+Z%*%D+E_M
      Y<-M%*%B+Z%*%E+X*C+E_Y
    }
    Y<-as.vector(Y)
  } else if(family=="binomial"){
    X<-rnorm(n,0,sqrt(1))
    #starting value of 0 for each variable: M, Y
    M<-matrix(0,n,p);Y<-rep(0,n)
    var_M=2

    if(q==0) { Z=NULL;D=NULL;E=NULL
    } else{
      Z<-matrix(rnorm(n*q,0,sqrt(1.0)),n,q)
      D<-matrix(0,q,p);E<-rep(0,q)
      for(i in 1:q){ E[i]<-0.5
      for(j in 1:p){ D[i,j]<-0.5}}
    }

    #Generate error
    mean_M_error <- rep(0,p)
    sigma_M_error <- matrix(0,nrow=p,ncol=p)
    for(i in 1:p){
      for(j in 1:p){
        sigma_M_error[i,j]<-rhoM^(abs(i-j))
      }
    }
    sigma_M_error<-var_M*sigma_M_error;E_M<-mvtnorm::rmvnorm(n,mean_M_error,sigma_M_error)

    if(q==0){
      M<-X%*%t(A)+E_M
      P_Y<-apply(M%*%B+X*C,1,sigmoid)
      Y<-stats::rbinom(n, 1, P_Y)
    } else{
      M<-X%*%t(A)+Z%*%D+E_M
      P_Y<-apply(M%*%B+Z%*%E+X*C,1,sigmoid)
      Y<-stats::rbinom(n, 1, P_Y)
    }
    Y<-as.vector(Y)
  } else if(family=="poisson"){
    X<-rnorm(n,0,sqrt(1))
    #starting value of 0 for each variable: M, Y
    M<-matrix(0,n,p);Y<-rep(0,n)
    var_M=2

    if(q==0) { Z=NULL;D=NULL;E=NULL
    } else{
      Z<-matrix(rnorm(n*q,0,sqrt(1.0)),n,q)
      D<-matrix(0,q,p);E<-rep(0,q)
      for(i in 1:q){ E[i]<-0.5
      for(j in 1:p){ D[i,j]<-0.5}}
    }

    #Generate error
    mean_M_error <- rep(0,p)
    sigma_M_error <- matrix(0,nrow=p,ncol=p)
    for(i in 1:p){
      for(j in 1:p){
        sigma_M_error[i,j]<-rhoM^(abs(i-j))
      }
    }
    sigma_M_error<-var_M*sigma_M_error;E_M<-mvtnorm::rmvnorm(n,mean_M_error,sigma_M_error)

    if(q==0){
      M<-X%*%t(A)+E_M
      P_Y<-apply(M%*%B+X*C,1,exp)
      Y<-stats::rpois(n, P_Y)
    } else{
      M<-X%*%t(A)+Z%*%D+E_M
      P_Y<-apply(M%*%B+Z%*%E+X*C,1,exp)
      Y<-stats::rpois(n, P_Y)
    }
    Y<-as.vector(Y)
  } else if(family=="cox"){
    var_M=2
    sigma_e <- matrix(0.0, p, p)
    for(i in 1:p){
      for(j in 1:p){
        sigma_e[i,j]<-rhoM^(abs(i-j))
      }
    }
    sigma_e<-var_M*sigma_e

    if(q==0) { Z=NULL;D=NULL;E=NULL
    } else{
      Z<-matrix(rnorm(n*q,0,sqrt(1.0)),n,q)
      D<-matrix(0,q,p);E<-rep(0,q)
      for(i in 1:q){
        E[i]<-0.5
        for(j in 1:p) D[i,j]<-0.5}
    }

    C <- matrix(0.5, 1, 1)
    X <- matrix(rnorm(n, mean = 0, sd = sqrt(1)), n, 1) # expoure
    e <- mvtnorm::rmvnorm(n, matrix(0, p, 1), sigma_e)
    M <- X%*%t(A) +Z%*%D + e
    MX <- cbind(M,  X)
    BC <- c(B, C)
    ## generate the failure time T
    u <- runif(n, 0, 1)
    T <- matrix(0, n, 1)
    for (i in 1:n)
      T[i] <- -log(1 - u[i])*exp(-(sum(BC*MX[i,])+sum(E*Z[i,])))
    ## generate censoring time 0.45 censoring rate
    CT <- runif(n, min = 0, max =3)
    status <- as.integer(T <= CT)
    ## the observed failure time
    OT <- apply(cbind(CT,T), 1, min)
    Y<-cbind(OT,status)
  } else if(family=="AFT"){
    var_M=2
    sigma_e <- matrix(0.0, p, p)
    for(i in 1:p){
      for(j in 1:p){
        sigma_e[i,j]<-rhoM^(abs(i-j))
      }
    }
    sigma_e<-var_M*sigma_e
    if(q==0) { Z=NULL;D=NULL;E=NULL
    } else{
      Z<-matrix(rnorm(n*q,0,sqrt(1.0)),n,q)
      D<-matrix(0,q,p);E<-rep(0,q)
      for(i in 1:q){
        E[i]<-0.5
        for(j in 1:p) D[i,j]<-0.5
      }
    }

    C <- matrix(0.5, 1, 1)
    X <- matrix(rnorm(n, mean = 0, sd = sqrt(1)), n, 1) # expoure
    mu <- matrix(0, p, 1)
    e_m <- mvtnorm::rmvnorm(n, mu, sigma_e)
    M <- X%*%t(A) +Z%*%D + e_m
    MXZ <- cbind(M, X, Z)
    BCE <- c(B, C,E)
    e_y<-rnorm(n, 0,0.5)
    # e_y<-rweibull(n, 1, scale=1)
    ## generate the failure time T
    T <- rep(0, n)
    for (i in 1:n)
      T[i] <- exp(sum(BCE*MXZ[i,])+e_y[i])
    CT <- runif(n, min = 0, max =5)
    status <- as.integer(T <= CT)
    ## the observed failure time
    OT <- apply(cbind(CT,T), 1, min)
    Y<-cbind(OT,status)
  } else {
    stop(paste0("Family ", family, " is not supported."))
  }
  data<-list(X=as.vector(X),M=as.matrix(M),Y=Y,Z=Z,A=A,B=B)
  return(data)
}
