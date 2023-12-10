#Simulation code for HIMAS with 1000 replications
library(survival);library(HIMA)
SurvHIMA<-function(X,Z,M,OT,status,FDRcut,verbose=TRUE){
  tryCatch({survHIMA(X=X, Z=Z, M=M, OT=OT, status=status,FDRcut = FDRcut,verbose=verbose)}
           ,  error = function(e) {
             NULL  #n,T setting-2
           }
  )
}

Sim_ms_surv<-function(Q){
  set.seed(1e3+Q)
  number <- 100  # sample size
  pnumber <- 100 # the dimension of mediators
  qnumber <- 2 # the number of covariates
  family <- "cox"  # the type of model
  # the regression coefficients A (X --> M)
  A<-rep(0,pnumber)
  # the regression coefficients B (M --> Y)
  B<-rep(0,pnumber)
  # the first four markers are true mediators.
  A[1:6] <- c(0.45, 0.4, -0.35, -0.3, 0.0, 0.5)
  B[1:6] <- c(0.5, -0.45, 0.4, -0.35, 0.5,0.0)
  rhoM<-0.0 #the correlation coefficient of random errors for M
  # Generate simulation data
  dat = simhima(number,pnumber,qnumber,family,A,B,rhoM,seed=1e3+Q)
  X<-dat$X;M<-dat$M;Y<-dat$Y;Z<-dat$Z
  A_t<-dat$A;B_t<-dat$B;C_t<-0
  colnames(Y)<-c("time","status")
  colnames(M) <- paste0("M", 1:ncol(M))
  survHIMA.fit <- SurvHIMA(X, Z, M, Y[,1], Y[,2],FDRcut =0.05)
  
  if(is.null(survHIMA.fit)){
    B_f<-rep(0,pnumber)
    A_f<-rep(0,pnumber)
    C_f<-0
    ID_hima<-rep(0,pnumber)
  } else{
    survHIMA.fit$ID
    B_f<-rep(0,pnumber)
    B_f[survHIMA.fit$ID]<-survHIMA.fit$beta
    A_f<-rep(0,pnumber)
    A_f[survHIMA.fit$ID]<-survHIMA.fit$alpha
    C_f<-0
    ID_hima<-rep(0,pnumber)
    ID_hima[survHIMA.fit$ID]<-1
  }
  
  data1<-data.frame(time=Y[,1],status=Y[,2])
  data1<-data.frame(data1,X)
  data1<-data.frame(data1,Z)
  fit <-coxph(Surv(time, status) ~ .,data=data1)
  C_f<-fit$coefficients[1]-sum(A_f*B_f)
  
  ID_true<-as.numeric(A_t*B_t!=0)
  c(number,pnumber,ID_true,A_t*B_t,C_t,ID_hima,A_f*B_f,C_f)
  
}

SIM_ms_surv<-function(Q){
  tryCatch({Sim_ms_surv(Q)}
           ,  error = function(e) {
             rep("Error",84)  #n,T setting-2
           }
  )
}
NUM<-1000
t1<-Sys.time()
library(parallel)
clnum<-detectCores()
cl <- makeCluster(getOption("cl.cores",20))
clusterExport(cl, c('simhima','Sim_ms_surv','SIM_ms_surv','survHIMA','SurvHIMA','Surv','coxph'))
system.time(result1<-parLapply(cl,1:NUM, SIM_ms_surv))
stopCluster(cl)
t2<-Sys.time()
t2-t1
result<-result1
j=0
for(i in 1:length(result)){
  if( result[[i-j]][1]=="Error") {result[[i-j]] <- NULL;j<-j+1}
  
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
R_HIMAS_1<-round(100*c(FWER,TPR,FNR,FPR,TNR,PPV,FOR,FDR,NPV,F1_score),2)
names(R_HIMAS_1)<-c("FWER","TPR","FNR","FPR","TNR","PPV","FOR","FDR","NPV","F1_score")
# R_HIMAS_1
sel<-cbind((apply(res1[1:pnumber,],1,mean)))
R_HIMAS_2<-sel
R_HIMAS_2<-round(100*R_HIMAS_2,3)
RES1<-matrix(0,pnumber+1,3);colnames(RES1)<-c("truevalue","mean","sd")
RES1[,1]=c(me_true,de_true)
for(i in 1:(pnumber+1)){
  RES1[i,2]=mean(res1[pnumber+i,])
  RES1[i,3]=mean((res1[pnumber+i,]-RES1[i,1])^2)
}

RES<-round(cbind(as.matrix(1:(pnumber+1)),RES1)[,],4)
RES[,3]<-RES[,3]-RES[,2]
RES[,4]<-round(sqrt(RES[,4]),4)
colnames(RES)<-c("","truevalue","bias","root_mse")
R_HIMAS_3<-RES
#Results of FWER, PR, FNR, FPR, TNR, PPV, FOR, FDR, NPV, and F1_score
R_HIMAS_1
#Estimation results of the mediation effects
R_HIMAS_3[c(1:10,pnumber+1),]
t_hima_surv