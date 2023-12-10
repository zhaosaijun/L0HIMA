library(HIMA)
Sim_ms<-function(Q){
  set.seed(1e3+Q)
  fdr=0.05
  n <- 100  # sample size
  p <- 100 # the dimension of mediators
  q <- 2 # the number of covariates
  family <- "gaussian"  #or "binomial" # the type of model
  # the regression coefficients A (X --> M)
  A<-rep(0,p)
  # the regression coefficients B (M --> Y)
  B<-rep(0,p)
  # the first four markers are true mediators.
  A[1:6] <- c(0.45, 0.4, -0.35, -0.3, 0.0, 0.5)
  B[1:6] <- c(0.5, -0.45, 0.4, -0.35, 0.5,0.0)
  rhoM<-0.0 #the correlation coefficient of random errors for M
  # Generate simulation data
  data<-simhima(n,p,q,family,A,B,rhoM,seed=1e3+Q)
  A_t<-data$A;B_t<-data$B;C_t<-0
  X<-data$X;M<-data$M;Y<-data$Y;Z=data$Z
  colnames(M) <- paste0("M", 1:ncol(M))
  hima.fit <-hima(X=as.matrix(X),Y=Y,M=M,COV.XM = Z, COV.MY = Z,
                  Y.family = "gaussian",M.family = "gaussian", penalty = "MCP",verbose = TRUE)
  hima.fit
  
  # #FWER
  # id<-NULL
  # if(sum(hima.fit[,6]<fdr)>0){
  #   id<-as.numeric(as.numeric(gsub('M','',rownames(hima.fit)))[hima.fit[,6]<fdr])
  #   A_h<-rep(0,p);B_h<-rep(0,p);C_h<-0
  #   A_h[id]=hima.fit[hima.fit[,6]<fdr,1];B_h[id]=hima.fit[hima.fit[,6]<fdr,2]
  #   C_h<-hima.fit[1,3]-sum(A_h*B_h)
  #
  # } else{
  #   A_h<-rep(0,p);B_h<-rep(0,p);C_h<-0
  # }
  # ID_hima<-rep(0,p)
  # ID_hima[id]<-1
  # id
  
  #FDR
  id<-NULL
  if(sum(hima.fit[,7]<fdr)>0){
    id<-as.numeric(as.numeric(gsub('M','',rownames(hima.fit)))[hima.fit[,7]<fdr])
    A_h<-rep(0,p);B_h<-rep(0,p);C_h<-0
    A_h[id]=hima.fit[hima.fit[,7]<fdr,1];B_h[id]=hima.fit[hima.fit[,7]<fdr,2]
    C_h<-hima.fit[1,3]-sum(A_h*B_h)
    
  } else{
    A_h<-rep(0,p);B_h<-rep(0,p);C_h<-0
  }
  ID_hima<-rep(0,p)
  ID_hima[id]<-1
  id
  
  ID_true<-as.numeric(A_t*B_t!=0)
  c(n,p,ID_true,A_t*B_t,C_t,ID_hima,A_h*B_h,C_h)
}

NUM<-1000
t1<-Sys.time()
library(parallel)
clnum<-detectCores()
cl <- makeCluster(getOption("cl.cores",20))
clusterExport(cl, c('simhima','hima'))
system.time(result1<-parLapply(cl,1:NUM, Sim_ms))
stopCluster(cl)
t2<-Sys.time()
t_hima=t2-t1
t_hima

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
R_hima_1<-round(100*c(FWER,TPR,FNR,FPR,TNR,PPV,FOR,FDR,NPV,F1_score),2)
names(R_hima_1)<-c("FWER","TPR","FNR","FPR","TNR","PPV","FOR","FDR","NPV","F1_score")
# R_hima_1
sel<-cbind((apply(res1[1:pnumber,],1,mean)))
R_hima_2<-sel
R_hima_2<-round(100*R_hima_2,3)
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
R_hima_3<-RES
#Results of FWER, PR, FNR, FPR, TNR, PPV, FOR, FDR, NPV, and F1_score
R_hima_1
##Estimation results of the mediation effects
R_hima_3[c(1:10,pnumber+1),c(2,3,4)]