#' @importFrom stats binomial coef complete.cases ks.test glm lm lsfit p.adjust pnorm rnorm runif sd
L0himams<-function(X, M, Y, COV.XM = NULL, COV.MY = NULL, family = c("gaussian", "binomial","poisson","cox","AFT")){
  n <- nrow(M);p <- ncol(M);q_MY<-ncol(COV.MY);q_XM<-ncol(COV.XM)
  MXZ<-cbind(M,X,COV.MY)

  # Estimate the regression coefficients b (mediators --> outcome)
  if(family=="gaussian"){
    b_SIS <- rep(0,p)
    for (i in 1:p){
      ID_S <- c(i, (p+1):(p+q_MY+1))
      MXZ_SIS <- MXZ[,ID_S]
      data1<-data.frame(MXZ_SIS,Y)
      fit <- glm(Y ~., data = data1)
      b_SIS[i] <- fit$coefficients[2]
    }
  } else if(family=="binomial"){
    b_SIS <- rep(0,p)
    for (i in 1:p){
      ID_S <- c(i, (p+1):(p+q_MY+1))
      MXZ_SIS <- MXZ[,ID_S]
      data1<-data.frame(MXZ_SIS,Y)
      fit <- glm(Y ~.,family=binomial(link='logit'), data = data1)
      b_SIS[i] <- fit$coefficients[2]
    }
  } else if(family=="poisson"){
    b_SIS <- rep(0,p)
    for (i in 1:p){
      ID_S <- c(i, (p+1):(p+q_MY+1))
      MXZ_SIS <- MXZ[,ID_S]
      data1<-data.frame(MXZ_SIS,Y)
      fit <- glm(Y ~.,family="poisson", data = data1)
      b_SIS[i] <- fit$coefficients[2]
    }
  }else if(family=="cox"){
    b_SIS <- rep(0,p)
    for(i in 1:p){
      data1<-data.frame(time=Y[,1],status=Y[,2],M[,i],X,COV.MY)
      fit <-survival::coxph(survival::Surv(time, status) ~ .,data=data1, ties="exact")
      b_SIS[i] <- fit$coefficients[1]
    }
  } else if(family=="AFT"){
    b_SIS <- rep(0,p)
    Y_order<-log(sort(Y[,1]))
    delta<-Y[order(Y[,1]),2]
    w<-rep(0,n)
    w[1]<-delta[1]/n
    for(i in 2:n){
      w[i]<-delta[i]/(n-i+1)
      for(j in 1:(i-1)) w[i]<-w[i]*((n-j)/(n-j+1))^delta[j]
    }
    for(i in 1:n)   w[i]<-w[i]*n
    W<-diag(sqrt(w))
    Y_order_n<-W%*%Y_order
    for(i in 1:p){
      MXZ<-cbind(M[,i],X,COV.MY,rep(1,n))
      MXZ_order<-MXZ[order(Y[,1]),]
      MXZ_order_n<-W%*%MXZ_order
      fit <- lsfit(MXZ_order_n,Y_order_n,intercept = FALSE)
      b_SIS[i] <- fit$coefficients[1]
    }
  } else {
    stop(paste0("Family ", family, " is not supported."))
  }


  # Estimate the regression coefficients a (predictor --> mediators)
  a_SIS <- rep(0,p)
  XZ <- cbind(X,COV.XM)
  for (i in 1:p){
    fit_a  <- lsfit(XZ,M[,i],intercept = TRUE)
    est_a <- matrix(coef(fit_a))[2]
    a_SIS[i] <- est_a
  }
  result=data.frame(A_est=a_SIS,B_est=b_SIS)
  return(result)
}

sigmoid <- function(x) {
  1 / ( 1 + exp(-x) )
}

null_estimation_adapt <- function(input_pvalues)
{
  ## updated function that automatically choose best lambda that result in better behave q-q plot

  if (is.null(ncol(input_pvalues)))
    stop("input_pvalues should be a matrix or data frame")
  if (ncol(input_pvalues) != 2)
    stop("inpute_pvalues should have 2 column")
  input_pvalues <- matrix(as.numeric(input_pvalues), nrow = nrow(input_pvalues))
  if (sum(complete.cases(input_pvalues)) < nrow(input_pvalues))
    warning("input_pvalues contains NAs to be removed from analysis")
  input_pvalues <- matrix(input_pvalues[complete.cases(input_pvalues),],nrow = nrow(input_pvalues))
  if (!is.null(nrow(input_pvalues)) & nrow(input_pvalues) <1)
    stop("input_pvalues doesn't have valid p-values")

  ### first identify features that may come from H11 (alternatives for both hypotheses) ###
  #library(qvalue)
  #ish11 <- qvalue(input_pvalues[,1])$qvalue<0.25 & qvalue(input_pvalues[,2])$qvalue<0.25

  pcut <- seq(0.1, 0.8, 0.1)
  frac1 <- rep(0, 8)
  frac2 <- rep(0, 8)
  frac12 <- rep(0, 8)
  for (i in 1:8) {
    frac1[i] <- mean(input_pvalues[, 1] >= pcut[i])/(1 - pcut[i])
    frac2[i] <- mean(input_pvalues[, 2] >= pcut[i])/(1 - pcut[i])
    frac12[i] <- mean(input_pvalues[, 2] >= pcut[i] & input_pvalues[,1] >= pcut[i])/(1 - pcut[i])^2
  }

  alphaout <- matrix(0,4,5)
  ll <- 1
  qqslope <- rep(0,4)
  for (lambda in c(0.5,0.6,0.7,0.8)) {
    alpha00 <- min(frac12[pcut >= lambda][1], 1)
    if (ks.test(input_pvalues[, 1], "punif", 0, 1, alternative = "greater")$p > 0.05)
      alpha1 <- 1 else alpha1 <- min(frac1[pcut >= lambda][1], 1)
      if (ks.test(input_pvalues[, 2], "punif", 0, 1, alternative = "greater")$p > 0.05)
        alpha2 <- 1 else alpha2 <- min(frac2[pcut >= lambda][1], 1)
        if (alpha00 == 1) {
          alpha01 <- 0
          alpha10 <- 0
          alpha11 <- 0
        } else {
          if (alpha1 == 1 & alpha2 == 1) {
            alpha01 <- 0
            alpha10 <- 0
            alpha11 <- 0
            alpha00 <- 1
          }
          if (alpha1 == 1 & alpha2 != 1) {
            alpha10 <- 0
            alpha11 <- 0
            alpha01 <- alpha1 - alpha00
            alpha01 <- max(0, alpha01)
            alpha00 <- 1 - alpha01
          }
          if (alpha1 != 1 & alpha2 == 1) {
            alpha01 <- 0
            alpha11 <- 0
            alpha10 <- alpha2 - alpha00
            alpha10 <- max(0, alpha10)
            alpha00 <- 1 - alpha10
          }
          if (alpha1 != 1 & alpha2 != 1) {
            alpha10 <- alpha2 - alpha00
            alpha10 <- max(0, alpha10)
            alpha01 <- alpha1 - alpha00
            alpha01 <- max(0, alpha01)
            if ((1 - alpha00 - alpha01 - alpha10) < 0) {
              alpha11 <- 0
              alpha10 <- 1 - alpha1
              alpha01 <- 1 - alpha2
              alpha00 <- 1 - alpha10 - alpha01
            }
            else {
              alpha11 <- 1 - alpha00 - alpha01 - alpha10
            }
          }
        }

        pmax <- apply(input_pvalues,1,max)
        pmax <- pmax[order(pmax)]
        nnulls <- sum(pmax>0.8)
        if(nnulls<2) break
        nmed <- nrow(input_pvalues)
        pexp <- rep(0,nnulls)
        for (i in 1:nmed) {
          c <- (-i/nmed)
          b <- alpha01+alpha10
          a <- 1-b
          if (alpha00==0) pexp[i] <--c/b else pexp[i] <- (-b+sqrt(b^2-4*a*c))/(2*a)
        }
        xx <- -log(pexp[(nmed-nnulls+1):nmed],base=10)
        yy <- -log(pmax[(nmed-nnulls+1):nmed],base=10)
        fit1 <- lm(yy~xx-1)

        qqslope[ll]<- fit1$coef[1]
        alphaout[ll,1] <- alpha10
        alphaout[ll,2] <- alpha01
        alphaout[ll,3] <- alpha00
        alphaout[ll,4] <- alpha1
        alphaout[ll,5] <- alpha2

        ll <- ll+1

  }

  if(nnulls>=2){
    bestslope <- which.min(qqslope)
    alpha.null <- list(alpha10 = alphaout[bestslope,1], alpha01 = alphaout[bestslope,2],alpha00 = alphaout[bestslope,3], alpha1 = alphaout[bestslope,4], alpha2 = alphaout[bestslope,5])
  } else alpha.null <- list(alpha10= alpha10, alpha01=alpha01,alpha00=alpha00 , alpha1=alpha1 , alpha2= alpha2)
  if(sum(as.vector(unlist(alpha.null)))==0) alpha.null=list(alpha10=0.5,alpha01=0.5,alpha00=0,alpha1=0.5,alpha2=0.5)
  return(alpha.null)
}

fdr_est_adapt<-function(alpha00,alpha01,alpha10,alpha1,alpha2,input_pvalues,exact=0){

  ## alpha10,alpha01,alpha00 are estimated three types of null proportions
  ## alpha1 is the marginal null proportion for first p-value
  ## alpha2 is the marginal null proportion for second p-value
  ## input pvalues are two columns of p-values
  ## alpha is the level of FWER to be control at
  ## exact=0 corresponding to the approximation used in section 2.2-2.3 in the paper, the default value for exact is 0
  ## exact=1 corresponding to the exact used in section 2.4 in the paper
  ## check input

  if (is.null(ncol(input_pvalues)))
    stop("input_pvalues should be a matrix or data frame")
  if (ncol(input_pvalues) !=2)
    stop("inpute_pvalues should have 2 column")
  input_pvalues <- matrix(as.numeric(input_pvalues),nrow=nrow(input_pvalues))
  if (sum(complete.cases(input_pvalues))<nrow(input_pvalues))
    warning("input_pvalues contains NAs to be removed from analysis")
  input_pvalues <- matrix(input_pvalues[complete.cases(input_pvalues),],nrow = nrow(input_pvalues))
  if (!is.null(nrow(input_pvalues)) & nrow(input_pvalues)<1)
    stop("input_pvalues doesn't have valid p-values")

  pmax <- apply(input_pvalues,1,max)
  nmed <- length(pmax)
  efdr1 <- rep(0,nmed)

  if (exact==0) {
    for (i in 1:nmed) {
      fdr11 <-  (pmax[i]*alpha01)/mean(pmax<=pmax[i])
      fdr12 <-  (pmax[i]*alpha10)/mean(pmax<=pmax[i])
      fdr2  <-  (pmax[i]*pmax[i]*alpha00)/mean(pmax<=pmax[i])
      efdr1[i] <- fdr11+fdr12+fdr2
    }
  }
  if (exact==1) {
    #library(fdrtool)
    ish11 <- qvalue::qvalue(input_pvalues[,1])$qvalue<0.25 & qvalue::qvalue(input_pvalues[,2])$qvalue<0.25


    out1 <- input_pvalues[!ish11,]
    nmed1 <- nrow(out1)

    nmed  <- nrow(input_pvalues)
    cdf12 <- input_pvalues

    xx1 <- c(0,out1[order(out1[,1]),1])
    yy1 <- c(0,seq(1,nmed1,by=1)/nmed1)
    gfit1<- fdrtool::gcmlcm(xx1,yy1,type="lcm")
    xknots1 <- gfit1$x.knots[-1]
    Fknots1 <- cumsum(diff(gfit1$x.knots)*gfit1$slope.knots)

    xx2 <- c(0,out1[order(out1[,2]),2])
    yy2 <- c(0,seq(1,nmed1,by=1)/nmed1)
    gfit2<- fdrtool::gcmlcm(xx2,yy2,type="lcm")
    xknots2 <- gfit2$x.knots[-1]
    Fknots2 <- cumsum(diff(gfit2$x.knots)*gfit2$slope.knots)

    alpha1 <- (alpha00+alpha01)/(alpha00+alpha01+alpha10)
    alpha2 <- (alpha00+alpha10)/(alpha00+alpha01+alpha10)

    if (alpha1!=1) Fknots1 <- (Fknots1 - alpha1*xknots1)/(1-alpha1) else Fknots1 <- rep(0,length(xknots1))
    if (alpha2!=1) Fknots2 <- (Fknots2 - alpha2*xknots2)/(1-alpha2) else Fknots2 <- rep(0,length(xknots2))

    orderq1 <- pmax
    orderq2 <- pmax

    gcdf1 <- pmax
    gcdf2 <- pmax
    for (i in 1:length(xknots1)) {
      if (i==1) {
        gcdf1[orderq1<=xknots1[i]] <- (Fknots1[i]/xknots1[i])*orderq1[orderq1<=xknots1[i]]
      } else {
        if (sum(orderq1>xknots1[i-1] & orderq1<=xknots1[i])>0){
          temp <- orderq1[orderq1>xknots1[i-1] & orderq1<=xknots1[i]]
          gcdf1[orderq1>xknots1[i-1] & orderq1<=xknots1[i]] <- Fknots1[i-1] + (Fknots1[i]-Fknots1[i-1])/(xknots1[i]-xknots1[i-1])*(temp-xknots1[i-1])
        }
      }
    }

    for (i in 1:length(xknots2)) {
      if (i==1) {
        gcdf2[orderq2<=xknots2[i]] <- (Fknots2[i]/xknots2[i])*orderq2[orderq2<=xknots2[i]]
      } else {
        if (sum(orderq2>xknots2[i-1] & orderq2<=xknots2[i])>0){
          temp <- orderq2[orderq2>xknots2[i-1] & orderq2<=xknots2[i]]
          gcdf2[orderq2>xknots2[i-1] & orderq2<=xknots2[i]] <- Fknots2[i-1] + (Fknots2[i]-Fknots2[i-1])/(xknots2[i]-xknots2[i-1])*(temp-xknots2[i-1])
        }
      }
    }


    gcdf1 <- ifelse(gcdf1>1,1,gcdf1)
    gcdf2 <- ifelse(gcdf2>1,1,gcdf2)

    cdf12[,1] <- gcdf1
    cdf12[,2] <- gcdf2

    for (i in 1:nmed) {
      fdr11 <-  (pmax[i]*cdf12[i,2]*alpha01)/mean(pmax<=pmax[i])
      fdr12 <-  (pmax[i]*cdf12[i,1]*alpha10)/mean(pmax<=pmax[i])
      fdr2  <-  (pmax[i]*pmax[i]*alpha00)/mean(pmax<=pmax[i])
      efdr1[i] <- fdr11+fdr12+fdr2
    }
  }

  efdr1.order <- efdr1[order(pmax,decreasing=T)]

  for (i in 2:nmed)  {
    efdr1.order[i] <- min(efdr1.order[i],efdr1.order[i-1])
  }

  efdr1 <- efdr1.order[rank(-pmax)]
  return(efdr=efdr1)
}

fwer_est_adapt <-function(alpha10,alpha01,alpha00,alpha1,alpha2,input_pvalues,alpha=0.05,exact=0) {

  ## alpha10,alpha01,alpha00 are estimated null proportions
  ## alpha1 is the estimated marginal null proportion for first p-value
  ## alpha2 is the estimated marginal null proportion for second p-value

  ## input_pvalues is data matrix with two columns of p-values
  ## alpha is the level of FWER to be controlled at, default 0.05
  ## exact=0 corresponding to the approximation used in section 2.2-2.3, the default value is 0
  ## exact=1 corresponding to the exact used in section 2.4
  ## check input

  if (is.null(ncol(input_pvalues)))
    stop("input_pvalues should be a matrix or data frame")
  if (ncol(input_pvalues) !=2)
    stop("inpute_pvalues should have 2 column")
  input_pvalues <- matrix(as.numeric(input_pvalues),nrow=nrow(input_pvalues))
  if (sum(complete.cases(input_pvalues))<nrow(input_pvalues))
    warning("input_pvalues contains NAs to be removed from analysis")
  input_pvalues <- input_pvalues[complete.cases(input_pvalues),]
  if (!is.null(nrow(input_pvalues)) & nrow(input_pvalues)<1)
    stop("input_pvalues doesn't have valid p-values")

  pmax <- apply(input_pvalues,1,max)
  nmed <- length(pmax)

  ## first compute the approximation fwer cut-off, using it as starting value if approximation method 2 is called
  c <- (-alpha)/nmed
  b <- alpha01+alpha10
  a <- alpha00
  if (exact==0) { if (a==0) fwer_alpha<- -c/b else  fwer_alpha<- (-b+sqrt(b^2-4*a*c))/(2*a) }


  if (exact==1) {
    ish11 <- qvalue::qvalue(input_pvalues[,1])$qvalue<0.25 & qvalue::qvalue(input_pvalues[,2])$qvalue<0.25

    out1 <- input_pvalues[!ish11,]
    nmed1 <- nrow(out1)

    xx1 <- c(0,out1[order(out1[,1]),1])
    yy1 <- c(0,seq(1,nmed1,by=1)/nmed1)

    gfit1<- fdrtool::gcmlcm(xx1,yy1,type="lcm")
    xknots1 <- gfit1$x.knots[-1]
    Fknots1 <- cumsum(diff(gfit1$x.knots)*gfit1$slope.knots)

    xx2 <- c(0,out1[order(out1[,2]),2])
    yy2 <- c(0,seq(1,nmed1,by=1)/nmed1)

    gfit2<- fdrtool::gcmlcm(xx2,yy2,type="lcm")
    xknots2 <- gfit2$x.knots[-1]
    Fknots2 <- cumsum(diff(gfit2$x.knots)*gfit2$slope.knots)

    alpha1 <- (alpha00+alpha01)/(alpha00+alpha01+alpha10)
    alpha2 <- (alpha00+alpha10)/(alpha00+alpha01+alpha10)

    if (alpha1!=1) Fknots1 <- (Fknots1 - alpha1*xknots1)/(1-alpha1) else Fknots1 <- rep(0,length(xknots1))
    if (alpha2!=1) Fknots2 <- (Fknots2 - alpha2*xknots2)/(1-alpha2) else Fknots2 <- rep(0,length(xknots2))

    fwer_alpha<- (-b+sqrt(b^2-4*a*c))/(2*a)
    qfwer <- fwer_alpha


    ite <- 1
    difff <- 1
    while (abs(difff)>1e-6 & ite<10){
      cat(ite,"..")
      if (sum(input_pvalues[,1]<fwer_alpha)<60 & sum(input_pvalues[,1]<fwer_alpha)>15) {
        if (alpha1==1) cdf1 <- 0 else {
          cdf1 <- max(0,(mean(input_pvalues[,1]<fwer_alpha) - alpha1*fwer_alpha)/(1-alpha1))
          cdf1 <- min(cdf1,1)
        }
      }  else{
        if (sum(input_pvalues[,1]<fwer_alpha)<=15) cdf1 <- 1
        if (sum(input_pvalues[,1]<fwer_alpha)>=60) {
          if (fwer_alpha<=xknots1[1]) cdf1 <- Fknots1[1] else {
            for (i in 2:length(xknots1)) {
              if (sum(fwer_alpha>xknots1[i-1] & fwer_alpha<=xknots1[i])>0){
                cdf1 <- Fknots1[i-1] + (Fknots1[i]-Fknots1[i-1])/(xknots1[i]-xknots1[i-1])*(fwer_alpha-xknots1[i-1])
              }
            }
            if (fwer_alpha>xknots1[length(xknots1)]) cdf1 <- 1
          }
        }
      }

      if (sum(input_pvalues[,2]<fwer_alpha)<60 & sum(input_pvalues[,2]<fwer_alpha)>15) {
        if (alpha2==1) cdf2 <- 0 else {
          cdf2 <- max(0,(mean(input_pvalues[,2]<fwer_alpha) - alpha2*fwer_alpha)/(1-alpha2))
          cdf2 <- min(cdf2,1)
        }
      }  else{
        if (sum(input_pvalues[,2]<fwer_alpha)<=15) cdf2 <- 1
        if (sum(input_pvalues[,2]<fwer_alpha)>=60) {
          if (fwer_alpha<=xknots2[1]) cdf2 <- Fknots2[1] else {
            for (i in 2:length(xknots2)) {
              if (sum(fwer_alpha>xknots2[i-1] & fwer_alpha<=xknots2[i])>0){
                cdf2 <- Fknots2[i-1] + (Fknots2[i]-Fknots2[i-1])/(xknots2[i]-xknots2[i-1])*(fwer_alpha-xknots2[i-1])
              }
            }
            if (fwer_alpha>xknots2[length(xknots2)]) cdf2 <- 1
          }
        }
      }
      c <- (-alpha/nmed)
      if (cdf1>1) cdf1 <- 1
      if (cdf2>1) cdf2 <- 1

      b <- alpha10*cdf1+alpha01*cdf2
      a <- alpha00
      fwer_alpha <- (-b+sqrt(b^2-4*a*c))/(2*a)
      difff <- max(qfwer-fwer_alpha)
      qfwer <- fwer_alpha
      ite <- ite+1
    }
  }
  return(fwer_alpha)
}

fit_AFT<- function(X, M, Y,COV.XM = NULL, COV.MY = NULL) {
  n <- nrow(M);p <- ncol(M);q_MY<-ncol(COV.MY);q_XM<-ncol(COV.XM)
  ID<-1:p;M_o<-M
  A<-rep(0,p);B<-rep(0,p);C<-0

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
  data1<-data.frame(Y=Y_order_n)
  data1<-data.frame(data1,X=MXZ_order_n)
  fit <- glm(Y ~.-1, data = data1)
  results<-fit$coefficients[1:p]
  return(results)
}

fit_AFT_bootstrap<- function(X, M, Y, COV.XM, COV.MY) {
  n<-nrow(M)
  nid<-sample(1:n, size = ceiling(0.632*n), replace = FALSE)
  est<-fit_AFT(X[nid], as.matrix(M[nid,]), Y[nid,], COV.XM[nid,], COV.MY[nid,])
  return(est)
}
