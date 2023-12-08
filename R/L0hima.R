#' L0-regularized high-dimensional mediation analysis
#'
#'\code{L0hima} is used to estimate and test high-dimensional mediation effects.
#'
#' @param X A vector of exposure.
#' @param Y The outcome variable, of \code{n} observations.
#' For \code{family = "gaussian"} \code{Y} should be a vector of continious outcome.
#' For \code{family = "binomial"} \code{Y} should be a vector of binomial outcome.
#' For \code{family = "poisson"} \code{Y} should be a vector of count outcome.
#' For \code{family = "cox"} or \code{"AFT"}, \code{Y} should be a two-column matrix with columns
#' named \code{"time"} and \code{"status"}.
#' @param M a \code{matrix} of mediators. Rows represent samples, columns represent variables.
#' @param COV.XM a \code{matrix} of covariates dataset for testing the association \code{X ~ M}.
#' Covariates specified here will not participate penalization.
#' @param COV.MY a \code{matrix} of covariates dataset for testing the association \code{M ~ Y}.
#' Covariates specified here will not participate penalization.
#' @param family One of the following models:
#' \code{"gaussian"} (continuous outcome),
#' \code{"binomial"} (binomial outcome),
#' \code{"poisson"} (count outcome),
#' \code{"cox"} (right-censored outcome),
#' \code{"AFT"} (right-censored outcome).
#' @param test.type The type of significant test. Available options are \code{"JS-uniform"}, and \code{"JS-mixture"}.
#' @param test.control The objective of significant test. Available options are \code{"FWER"}, and \code{"FDR"}.
#' @param sig.level The significance level of test. Default is 0.05.
#' @param screening Whether to conduct the sure independent screening. Default = \code{TRUE}.
#' @param topN An integer specifying the number of top markers from sure independent screening. Default = \code{NULL}.
#'
#' @return A data.frame containing results of L0-regularized high-dimensional mediation analysis.
#' \itemize{
#'     \item{SM: }{The index set of selected mediators.}
#'     \item{A_est: }{estimates of coefficient A (predictor (X) --> mediators (M)).}
#'     \item{B_est: }{estimates of coefficient B (mediators (M) --> outcome (Y)).}
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
#' L0hima.fit<-L0hima(dat$X, dat$M, dat$Y, dat$Z, dat$Z, family,
#' test.type="JS-uniform",test.control="FDR", sig.level=0.05 )
#' L0hima.fit
#' @export
L0hima<- function(X, M, Y,COV.XM = NULL, COV.MY = NULL, family = c("gaussian", "binomial","poisson","cox","AFT"), test.type=c("JS-uniform","JS-mixture"),
                  test.control=c("FWER","FDR"), sig.level=0.05,screening=FALSE,topN = NULL) {
  n <- nrow(M);p <- ncol(M);q_MY<-ncol(COV.MY);q_XM<-ncol(COV.XM)
  ID<-1:p;M_o<-M
  A<-rep(0,p);B<-rep(0,p);C<-0

  #STEP 1 Mediators screening
  if(screening==TRUE){   if (is.null(topN)) d <- 6*n  else d <- topN
  d <- min(p, d) # if d > p select all mediators
  SIS_results<-L0himams(X, M, Y, COV.XM, COV.MY, family)
  SIS_est<-abs(SIS_results[,1]*SIS_results[,2])
  SIS_sort <- sort(SIS_est,decreasing = TRUE)
  ID <- which(SIS_est>= SIS_sort[d])  # the index of top mediators
  M <- M[, ID];p <- ncol(M)
  }

  #STEP 2 L0-regularized estimates
  if(family=="gaussian"){
    MXZ<-cbind(M,X,COV.MY)
    abess_fit <- abess::abess(MXZ,Y, family = "gaussian", tune.type = "gic",
                              always.include=(p+1):(p+q_MY+1))
    B[ID]<-abess_fit$beta[1:p,abess_fit$best.size-q_MY]
    C<-abess_fit$beta[p+1,abess_fit$best.size-q_MY]
    SB<-which(B!=0)
    M <- matrix(M_o[, SB],n,length(SB))
  } else if(family=="binomial"){
    MXZ<-cbind(M,X,COV.MY)
    abess_fit <- abess::abess(MXZ,Y, family = "binomial", tune.type = "gic",
                              always.include=(p+1):(p+q_MY+1),num.threads=1)
    B[ID]<-abess_fit$beta[1:p,abess_fit$best.size-q_MY]
    C<-abess_fit$beta[p+1,abess_fit$best.size-q_MY]
    SB<-which(B!=0)
    M <- matrix(M_o[, SB],n,length(SB))
  } else if(family=="poisson"){
    MXZ<-cbind(M,X,COV.MY)
    abess_fit <- abess::abess(MXZ,Y, family = "poisson", tune.type = "gic",
                              always.include=(p+1):(p+q_MY+1),num.threads=1)
    B[ID]<-abess_fit$beta[1:p,abess_fit$best.size-q_MY]
    C<-abess_fit$beta[p+1,abess_fit$best.size-q_MY]
    SB<-which(B!=0)
    M <- matrix(M_o[, SB],n,length(SB))
  } else if(family=="cox"){
    MXZ<-cbind(M,X,COV.MY)
    abess_fit <- abess::abess(MXZ,Y, family = "cox", tune.type = "gic",
                              always.include=(p+1):(p+q_MY+1),num.threads=1)
    B[ID]<-abess_fit$beta[1:p,abess_fit$best.size-q_MY]
    C<-abess_fit$beta[p+1,abess_fit$best.size-q_MY]
    SB<-which(B!=0)
    M <- matrix(M_o[, SB],n,length(SB))
  } else if(family=="AFT"){
    MXZ<-cbind(M,X,COV.MY,rep(1,n))
    Y_order<-log(sort(Y[,1]))
    delta<-Y[order(Y[,1]),2]
    MXZ_order<-MXZ[order(Y[,1]),]
    w<-rep(0,n)
    w[1]<-delta[1]/n
    for(i in 2:n){
      w[i]<-delta[i]/(n-i+1)
      for(j in 1:(i-1)) w[i]<-w[i]*((n-j)/(n-j+1))^delta[j]
    }
    for(i in 1:n)   w[i]<-w[i]*n
    W<-diag(sqrt(w))
    MXZ_order_n<-W%*%MXZ_order
    Y_order_n<-W%*%Y_order
    abess_fit <- abess::abess(MXZ_order_n,Y_order_n, family = "gaussian", tune.type = "gic",
                              always.include=(p+1):(p+q_MY+2),num.threads=1)
    BC<-as.vector(abess_fit$beta[,abess_fit$best.size-q_MY-1])
    B[ID]<-BC[1:p];C<-BC[p+1]
    SB<-which(B!=0)
    M <- matrix(M_o[, SB],n,length(SB))
    MXZ_order_n<- matrix(MXZ_order_n[, which(BC!=0)],n)
  }  else {
    stop(paste0("Family ", family, " is not supported."))
  }

  #STEP3 Significant test
  if(length(SB)==0){
    A_f<-rep(0,ncol(M_o));B_f<-rep(0,p);C_f<-C
    SM<-NULL
  } else{
    p_A<-rep(0,length(SB));A<-rep(0,ncol(M_o))
    for(i in 1:length(SB)){
      if(q_XM==0){
        MX1 <- data.frame(M = M[, i], X = X)
      }else{
        MX1 <- data.frame(M = M[, i], X = X, COV =COV.XM)
      }
      fit <- glm(M ~., data = MX1)
      A[SB[i]] <- summary(fit)$coef[2,1]           #coefficients for a
      p_A[i]<- summary(fit)$coef[2,4]         #p-value for a
    }
    if(family=="gaussian"){
      data1<-data.frame(Y=Y,M=M,X=X,COV=COV.MY)
      fit <- glm(Y ~., data = data1)
      p_B<-rep(0,length(SB))
      p_B<- summary(fit)$coef[2:(1+length(SB)),4]         #p-value for b
    } else if(family=="binomial"){
      data1<-data.frame(Y=Y,M=M,X=X,COV=COV.MY)
      fit <- glm(Y ~.,family=binomial(link='logit'), data = data1)
      p_B<-rep(0,length(SB))
      p_B<- summary(fit)$coef[2:(1+length(SB)),4]          #p-value for b
    } else if(family=="poisson"){
      data1<-data.frame(Y=Y,M=M,X=X,COV=COV.MY)
      fit <- glm(Y ~.,family = 'poisson', data = data1)
      p_B<-rep(0,length(SB))
      p_B<- summary(fit)$coef[2:(1+length(SB)),4]          #p-value for b
    }else if(family=="cox"){
      data1<-data.frame(time=Y[,1],status=Y[,2],M=M,X=X,COV=COV.MY)
      fit <-survival::coxph(survival::Surv(time, status) ~ .,data=data1, ties="exact")
      p_B<-rep(0,length(SB))
      p_B<- summary(fit)$coef[1:length(SB),5]         #p-value for b
    } else if(family=="AFT"){
      est_bootsrap<-replicate(n=500,fit_AFT_bootstrap(X,M,Y,COV.XM,COV.MY))
      if(ncol(M)==1) sd_bootsrap<-sd(est_bootsrap) else sd_bootsrap<-apply(est_bootsrap,1,sd)
      if(ncol(M)==1) mean_bootsrap<-mean(est_bootsrap) else mean_bootsrap<-apply(est_bootsrap,1,mean)
      p_B<-rep(0,length(SB))
      p_B<-2*(1-pnorm(abs(B[SB]/sd_bootsrap)))  #p-value for b
    }  else {
      stop(paste0("Family ", family, " is not supported."))
    }
    if(test.type=="JS-uniform"){
      if(test.control=="FDR") {
        p_fdr_A <- p.adjust(p_A, "fdr");p_fdr_B <- p.adjust(p_B, "fdr") } else if(test.control=="FWER"){
          p_fdr_A <- p.adjust(p_A, "bonferroni");p_fdr_B <- p.adjust(p_B, "bonferroni")
        } else {
          stop(paste0("Test.control: ", test.control, " is not supported."))
        }
      p_fdr_A[p_fdr_A> 1] <- 1;p_fdr_B[p_fdr_B> 1] <- 1
      PA<- cbind(p_fdr_A, p_fdr_B)
      p_value <- apply(PA,1,max)
      SM<-SB[which(p_value<sig.level)]#The index set of selected mediators
    } else if(test.type=="JS-mixture"){
      PA <- cbind(p_A,p_B)
      N0 <- dim(PA)[1]*dim(PA)[2]
      input_pvalues <- PA + matrix(runif(N0,0,10^{-10}),dim(PA)[1],2)
      # Estimate the proportions of the three component nulls
      nullprop <- null_estimation_adapt(input_pvalues)
      if(test.control=="FDR") {
        # Compute the estimated pointwise FDR for every observed p-max
        p_value<-  fdr_est_adapt(nullprop$alpha00,nullprop$alpha01,nullprop$alpha10, nullprop$alpha1,nullprop$alpha2,input_pvalues,exact=0)
        SM<-SB[which(p_value<sig.level)]#The index set of selected mediators
      } else if(test.control=="FWER"){
        p_value <- apply(PA,1,max)
        # Compute FWER cutoff for p-max
        FWER_cutoff<-  fwer_est_adapt(nullprop$alpha10,nullprop$alpha01,nullprop$alpha00, nullprop$alpha1,nullprop$alpha2,input_pvalues,alpha=sig.level,exact=0)
        SM<-SB[which(p_value<FWER_cutoff)]#The index set of selected mediators
      } else {
        stop(paste0("Test.control: ", test.control, " is not supported."))
      }
    } else {
      stop(paste0("Test type: ", test.type, " is not supported."))
    }
  }
  A_est=A[SM];B_est=B[SM];p_value=as.vector(p_value[which(p_value<sig.level)])
  if(length(SM)==0) {SM<-NULL;A_est=NULL;B_est=NULL;p_value=NULL}
  results<-data.frame(SM=SM,A_est=A_est,B_est=B_est,p_value=p_value)
  results


}
