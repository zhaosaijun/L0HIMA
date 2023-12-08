#Code for the real data analysis

#real_data_1_linear
M_1<-read.xlsx("/Users/zhaosaijun/Desktop/DNA_real_data/M_1.xlsx")
M_2<-read.xlsx("/Users/zhaosaijun/Desktop/DNA_real_data/M_2.xlsx")
XYZ<-read.xlsx("/Users/zhaosaijun/Desktop/DNA_real_data/XY.xlsx")

#M_1,M_2, and XYZ are the raw data
#data preprocessing
X<-XYZ[,5];Y<-XYZ[,4];Z1<-XYZ[,2];Z2<-XYZ[,3]
for(i in 1:85){
  if(Z2[i]=="male") Z2[i]=0
  else Z2[i]=1
}
Z<-cbind(as.numeric(Z1),as.numeric(Z2))
M<-cbind(M_1[,-1],M_2[,-1]);M<-t(M)
colnames(M) <- paste0("M", 1:ncol(M));colnames(Z) <- paste0("Z", 1:ncol(Z));DNAm<-M_1[,1]
X<-scale(X,center = TRUE, scale = TRUE)
Z<-scale(Z,center = TRUE, scale = TRUE)
Y<-scale(Y,center = TRUE, scale = TRUE)

L0hima.fit_linear<-L0hima(X, M, Y, Z, Z, family="gaussian", test.type="JS-uniform",test.control="FDR",
                  sig.level=0.05,screening = TRUE, topN = 6*nrow(M))
L0hima.fit_linear


#real_data2
data_raw_1<-read.table("/Users/zhaosaijun/Desktop/real data/HumanMethylation450")
data_raw_2<-read.csv("/Users/zhaosaijun/Desktop/real data/clinical_754.csv")

#data_raw_1 and data_raw_2 are the raw data
#data preprocessing
Tes<-function(P){
  sum(is.na(P))
}
D<-NULL;D<-apply(data_raw_1,1,Tes)
CPGnmae<-data_raw_1[,1]
data_raw_1<-data_raw_1[-which(D>0),];colnames(data_raw_1)<-data_raw_1[1,];rownames(data_raw_1)<-data_raw_1[,1]
data_raw_1<-data_raw_1[-1,-1];data_raw_1<-t(data_raw_1)
id<-data_raw_2[,1]
M<-matrix(as.numeric(data_raw_1[id,]),nrow=length(id))
Y<-matrix(0,nrow(M_s),2);Y[,1]=data_raw_2[,2];Y[,2]=data_raw_2[,3];colnames(Y)<-c("time","status")
Z<-matrix(0,nrow(M_s),4)
Z[,1]=data_raw_2[,7];Z[,2]=as.numeric(data_raw_2[,10]=="MALE");Z[,4]=as.numeric(data_raw_2[,9]=="YES")
for(i in 1:nrow(M_s)){
  if(data_raw_2[i,8]=="Stage IA"|data_raw_2[i,8]=="Stage IB"|data_raw_2[i,8]=="Stage I") Z[i,3]=1
  if(data_raw_2[i,8]=="Stage IIA"|data_raw_2[i,8]=="Stage IIB"|data_raw_2[i,8]=="Stage II") Z[i,3]=2
  if(data_raw_2[i,8]=="Stage IIIA"|data_raw_2[i,8]=="Stage IIIB"|data_raw_2[i,8]=="Stage III") Z[i,3]=3
  if(data_raw_2[i,8]=="Stage IVA"|data_raw_2[i,8]=="Stage IVB"|data_raw_2[i,8]=="Stage IV") Z[i,3]=4
}
colnames(Z)<-c("age","sex","cancerstage","radiotherapyindicator")
Z_s[,1]<-scale(Z[,1],center = TRUE, scale = TRUE);Z_s[,2]<-scale(Z[,2],center = TRUE, scale = TRUE)
Z_s[,3]<-scale(Z[,3],center = TRUE, scale = TRUE);Z_s[,4]<-scale(Z[,4],center = TRUE, scale = TRUE)
X<-data_raw_2[,6];X<-scale(X,center = TRUE, scale = TRUE)
Y[,1]<-scale(Y[,1],center =FALSE, scale = TRUE)

L0hima.fit_cox<-L0hima(X, M, Y, Z, Z, family="cox", test.type="JS-uniform",test.control="FDR",
                      sig.level=0.05,screening = TRUE, topN = 6*nrow(M))
L0hima.fit_cox

L0hima.fit_AFT<-L0hima(X, M, Y, Z, Z, family="AFT", test.type="JS-uniform",test.control="FDR",
                 sig.level=0.05,screening = TRUE, topN = 6*nrow(M))
L0hima.fit_AFT

