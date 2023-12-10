#Simulation code for HIMA2 with 1000 replications
library(hdi);library(HDMT)
filterSel<- function(X, M, Y, dir = TRUE) {
  y <- Y
  x <- X
  
  n <- length(x)                                  # number of samples
  cpx <- c(crossprod(x))                          # cross product of X
  cpxM <- apply(M, 2, crossprod, y = x)           # cross product of x and M
  alpha <- cpxM/cpx                               # alpha paths
  res_M <- M - tcrossprod(x, alpha)               # residual of m~x+0
  var_M <- apply(res_M, 2, crossprod) / (n - 1)   # rss variance
  var_a <- var_M / cpx                            # variance of alpha
  
  beta <- var_b <- numeric(ncol(M))
  for (i in 1:ncol(M)) {
    m <- M[,i]
    # then the beta path
    if (dir) {
      mm <- cbind(x, m)                                 # model matrix
    } else {
      mm <- cbind(m)
    }
    cpm <- crossprod(mm)                                # cross product of mm
    b <- solve(cpm, crossprod(mm, y))                   # beta
    res_y <- y - mm %*% c(b)                            # residual of y~m+x+0
    var_y <- as.numeric(crossprod(res_y) / (n - 1))     # rss variance
    vb <- diag(var_y * chol2inv(chol(cpm)))             # variance of beta
    beta[i] <- b[2]
    var_b[i] <- vb[2]
  }
  
  stat <- alpha * beta # product of coefficients
  se <- sqrt(alpha^2 * var_b + beta^2 * var_a) # estimated standard error Sobel test - var_a * var_b
  
  return(stat/se)
}
HIMA2<-function(X,Y,M,Z)
{
  n <- dim(X)[1]  # number of samples
  p <- dim(M)[2]  # number of mediators
  d <- dim(X)[2]  # number of exposures
  q <- dim(Z)[2]  # number of covariates
  
  MZX<-cbind(M,Z,X)
  
  #########################################################################
  ########################### (Step 1) SIS step ###########################
  #########################################################################
  message("Step 1: Sure Independent Screening ...", "  (", Sys.time(), ")")
  
  d_0 <- min(2*round(n/log(n)),p)
  beta_SIS <- matrix(0,1,p)
  
  # Estimate the regression coefficients beta (mediators --> outcome)
  for (i in 1:p){
    ID_S <- c(i, (p+1):(p+q+1))
    MZX_SIS <- MZX[,ID_S]
    fit <- lsfit(MZX_SIS,Y,intercept = TRUE)
    beta_SIS[i] <- fit$coefficients[2]
  }
  
  # Estimate the regression coefficients alpha (exposure --> mediators)
  alpha_SIS <- matrix(0,1,p)
  XZ <- cbind(X,Z)
  for (i in 1:p){
    fit_a  <- lsfit(XZ,M[,i],intercept = TRUE)
    est_a <- matrix(coef(fit_a))[2]
    alpha_SIS[i] <- est_a
  }
  
  # Select the d_0 number of mediators with top largest effect
  ab_SIS <- alpha_SIS*beta_SIS
  ID_SIS  <- which(-abs(ab_SIS) <= sort(-abs(ab_SIS))[d_0])
  d <- length(ID_SIS)
  
  #########################################################################
  ################### (Step 2) De-biased Lasso Estimates ##################
  #########################################################################
  message("Step 2: De-biased Lasso Estimates ...", "   (", Sys.time(), ")")
  
  P_beta_SIS <- matrix(0,1,d)
  beta_DLASSO_SIS_est <- matrix(0,1,d)
  beta_DLASSO_SIS_SE <- matrix(0,1,d)
  MZX_SIS <- MZX[,c(ID_SIS, (p+1):(p+q+1))]
  
  DLASSO_fit <- lasso.proj(x=MZX_SIS, y=Y, family = "gaussian",Z = NULL)
  beta_DLASSO_SIS_est <- DLASSO_fit$bhat[1:d]
  beta_DLASSO_SIS_SE <- DLASSO_fit$se
  P_beta_SIS <- t(DLASSO_fit$pval[1:d])
  
  ################### Estimate alpha ################
  alpha_SIS_est <- matrix(0,1,d)
  alpha_SIS_SE <- matrix(0,1,d)
  P_alpha_SIS <- matrix(0,1,d)
  
  XZ <- cbind(X,Z)
  for (i in 1:d){
    fit_a  <- lsfit(XZ,M[,ID_SIS[i]],intercept = TRUE)
    est_a <- matrix(coef(fit_a))[2]
    se_a <- ls.diag(fit_a)$std.err[2]
    sd_1 <- abs(est_a)/se_a
    P_alpha_SIS[i] <- 2*(1-pnorm(sd_1,0,1))  ## the SIS for alpha
    alpha_SIS_est[i] <- est_a
    alpha_SIS_SE[i] <- se_a
  }
  
  #########################################################################
  ################ (step 3) The multiple-testing  procedure ###############
  #########################################################################
  message("Step 3: Joint significance test ...", "     (", Sys.time(), ")")
  
  PA <- cbind(t(P_alpha_SIS),(t(P_beta_SIS)))
  P_value <- apply(PA,1,max)  #The joint p-values for SIS variable
  
  N0 <- dim(PA)[1]*dim(PA)[2]
  input_pvalues <- PA + matrix(runif(N0,0,10^{-10}),dim(PA)[1],2)
  
  # Estimate the proportions of the three component nulls
  nullprop <- null_estimation(input_pvalues)
  
  # Compute the estimated pointwise FDR for every observed p-max
  fdrcut  <- fdr_est(nullprop$alpha00,nullprop$alpha01,nullprop$alpha10, nullprop$alpha1,nullprop$alpha2,input_pvalues,exact=0)
  
  ID_fdr <- which(fdrcut <= 0.05)
  
  # Following codes extract the estimates for mediators with fdrcut<=0.05
  beta_hat_est <- beta_DLASSO_SIS_est[ID_fdr]
  beta_hat_SE  <- beta_DLASSO_SIS_SE[ID_fdr]
  
  alpha_hat_est <-  alpha_SIS_est[ID_fdr]
  alpha_hat_SE  <-  alpha_SIS_SE[ID_fdr]
  
  P.value_raw <- P_value[ID_fdr]
  
  # Indirect effect
  IDE <- beta_hat_est*alpha_hat_est # mediation(indirect) effect
  
  # Here we name the mediators as M1-Mp and extract the names of significant ones.
  M<-(sprintf("M%d", 1:p))[ID_SIS[ID_fdr]]
  
  # create a data frame with output values
  output<-data.frame(cbind(M, alpha=alpha_hat_est,alpha_SE=alpha_hat_SE,beta=beta_hat_est,beta_SE=beta_hat_SE,"alpha*beta"=IDE,
                           p_val=P.value_raw))
  
  message("Done!", "     (", Sys.time(), ")")
  
  return(output)
}
null_estimation_adaptation <- function(input_pvalues)
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


Sim_ms<-function(Q){
  set.seed(1e3+Q)
  fdr=0.05
  number <- 100  # sample size
  pnumber <- 100 # the dimension of mediators
  qnumber <- 2 # the number of covariates
  family <- "gaussian"  # the type of model
  # the regression coefficients A (X --> M)
  A<-rep(0,pnumber)
  # the regression coefficients B (M --> Y)
  B<-rep(0,pnumber)
  # the first four markers are true mediators.
  A[1:6] <- c(0.45, 0.4, -0.35, -0.3, 0.0, 0.5)
  B[1:6] <- c(0.5, -0.45, 0.4, -0.35, 0.5,0.0)
  rhoM<-0.0 #the correlation coefficient of random errors for M
  # Generate simulation data
  data<-simhima(number,pnumber,qnumber,family,A,B,rhoM,seed=1e3+Q)
  A_t<-data$A;B_t<-data$B;C_t<-0
  X<-data$X;M<-data$M;Y<-data$Y;Z=data$Z
  
  hima2.fit <-tryCatch({HIMA2(X=as.matrix(X),Y=Y,M=M,Z = Z)},  error = function(e) { rep("Error",80)  })
  hima2.fit
  
  
  if(as.vector(hima2.fit)[1]=="Error") {
    ID_true<-as.numeric(A_t*B_t!=0)
    results<-c(number,pnumber,ID_true,A_t*B_t,C_t,rep(0,pnumber),rep(0,pnumber),0) }else{
      #FDR
      id<-NULL
      if(sum(as.numeric(hima2.fit[,7])<fdr)>0){
        id<-as.numeric(as.numeric(gsub('M','',rownames(hima2.fit)))[as.numeric(hima2.fit[,7])<fdr])
        A_h<-rep(0,pnumber);B_h<-rep(0,pnumber);C_h<-0
        A_h[id]=as.numeric(hima2.fit[as.numeric(hima2.fit[,7])<fdr,2]);B_h[id]=as.numeric(hima2.fit[as.numeric(hima2.fit[,7])<fdr,4])
        data0<-data.frame(Y=Y)
        data0<-data.frame(data0,X)
        data0<-data.frame(data0,Z)
        fit <- lm(Y ~.-1, data = data0)
        C_h<-fit$coefficients[1]-sum(A_h*B_h)
        
      } else{
        A_h<-rep(0,pnumber);B_h<-rep(0,pnumber);C_h<-0
      }
      ID_hima2<-rep(0,pnumber)
      ID_hima2[id]<-1
      id
      
      ID_true<-as.numeric(A_t*B_t!=0)
      results<-c(number,pnumber,ID_true,A_t*B_t,C_t,ID_hima2,A_h*B_h,C_h)
    }
  
  return(results)
  
}


NUM<-100
t1<-Sys.time()
library(parallel)
clnum<-detectCores()
cl <- makeCluster(getOption("cl.cores",20))
clusterExport(cl, c('simhima','Sim_ms','SIM_ms','HIMA2','filterSel',
                    'lasso.proj','null_estimation','fdr_est','tryCatch'))
system.time(result1<-parLapply(cl,1:NUM, Sim_ms))
stopCluster(cl)
t2<-Sys.time()
t_hima2=t2-t1
t_hima2
result<-result1
N1<-length(result)
j=0
for(i in 1:N1){
  if(result[[i-j]][1]=="Error") {result[[i-j]] <- NULL;j<-j+1}
  
}
length(result)
number=result[[1]][1];pnumber=result[[1]][2] #n,p setting-2
L<-length(result[[1]]);NUM=length(result)
res1<-matrix(unlist(result),L,length(result))
id_true<-result[[1]][(3:(2+pnumber))]
me_true<-result[[1]][((3+pnumber):(2+2*pnumber))]
de_true<-result[[1]][3+2*pnumber]
res1<-res1[-(1:(3+2*pnumber)),]

FWER_it<-rep(0,NUM);TPR_it<-rep(0,NUM);FNR_it<-rep(0,NUM);FPR_it<-rep(0,NUM);TNR_it<-rep(0,NUM)
PPV_it<-rep(0,NUM);FOR_it<-rep(0,NUM);FDR_it<-rep(0,NUM);NPV_it<-rep(0,NUM)
F1_score_it<-rep(0,NUM)
for(i in 1:NUM){
  FWER_it[i]<-sum((res1[1:pnumber,i]-id_true)==1)>0
  TPR_it[i]<-sum(((res1[1:pnumber,i]-id_true)==0)&(id_true==1))/sum(id_true==1)
  FNR_it[i]<-sum(((res1[1:pnumber,i]-id_true)==-1)&(id_true==1))/sum(id_true==1)
  FPR_it[i]<-sum((res1[1:pnumber,i]-id_true)==1)/sum(id_true==0)
  TNR_it[i]<-sum(((res1[1:pnumber,i]-id_true)==0)&(id_true==0))/sum(id_true==0)
  PPV_it[i]<-sum((res1[1:pnumber,i]==1)&(id_true==1))/sum(res1[1:pnumber,i]==1)
  if(is.nan(PPV_it[i])) PPV_it[i]=1
  FOR_it[i]<-sum(((res1[1:pnumber,i]-id_true)==-1)&(id_true==1))/sum(res1[1:pnumber,i]==0)
  FDR_it[i]<-sum((res1[1:pnumber,i]-id_true)==1)/sum(res1[1:pnumber,i]==1)
  if(is.nan(FDR_it[i])) FDR_it[i]=0
  NPV_it[i]<-sum(((res1[1:pnumber,i]-id_true)==0)&(id_true==0))/sum(res1[1:pnumber,i]==0)
  if(TPR_it[i]==0) {F1_score_it[i]=0}
  else {F1_score_it[i]<-2*(PPV_it[i]*TPR_it[i])/(PPV_it[i]+TPR_it[i])}
}
FWER<-mean(FWER_it,na.rm=T)
TPR<-mean(TPR_it,na.rm=T)
FNR<-mean(FNR_it,na.rm=T)
FPR<-mean(FPR_it,na.rm=T)
TNR<-mean(TNR_it,na.rm=T)
PPV<-mean(PPV_it,na.rm=T)
FOR<-mean(FOR_it,na.rm=T)
FDR<-mean(FDR_it,na.rm=T)
NPV<-mean(NPV_it,na.rm=T)
F1_score<-mean(F1_score_it,na.rm=T)
R_hima2_1<-round(100*c(FWER,TPR,FNR,FPR,TNR,PPV,FOR,FDR,NPV,F1_score),2)
names(R_hima2_1)<-c("FWER","TPR","FNR","FPR","TNR","PPV","FOR","FDR","NPV","F1_score")
# R_hima2_1
sel<-cbind((apply(res1[1:pnumber,],1,mean)))
R_hima2_2<-sel
R_hima2_2<-round(100*R_hima2_2,3)
RES1<-matrix(0,pnumber+1,4);colnames(RES1)<-c("truevalue","mean","mse","sd")
RES1[,1]=c(me_true,de_true)
for(i in 1:(pnumber+1)){
  RES1[i,2]=mean(res1[pnumber+i,])
  RES1[i,3]=mean((res1[pnumber+i,]-RES1[i,1])^2)
  RES1[i,4]=sd(res1[pnumber+i,])^2
}

RES<-round(cbind(as.matrix(1:(pnumber+1)),RES1)[,],3)
RES[,3]<-RES[,3]-RES[,2]
RES[,4]<-round(sqrt(RES[,4]),3);RES[,5]<-round(sqrt(RES[,5]),3)
colnames(RES)<-c("","truevalue","bias","root_mse","root_sd")

R_hima2_3<-RES
#Results of FWER, PR, FNR, FPR, TNR, PPV, FOR, FDR, NPV, and F1_score
R_hima2_1
#Estimation results of the mediation effects
R_hima2_3[c(1:10,pnumber+1),c(2,3,4)]