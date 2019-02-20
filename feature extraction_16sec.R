# peaks fetaure extraction
# 
install.packages("moments")
install.packages("pracma")
install.packages("e1071")

require(e1071); require(pracma); require(moments)
##defining the functions used in features extraction
cv<-function(vector){
  return(sd(vector, na.rm = TRUE)/mean(vector, na.rm = TRUE))*100
}
pow.625 <- function(vm,M)
{
  mods <- Mod(fft(vm))
  mods <- mods[-1]
  n <- length(mods)
  n <- floor(n/2)
  freq <- (M*(1:n))/(2*n)
  mods <- mods[1:n]
  inds <- (1:n)[(freq>0.6)&(freq<2.5)]
  pow625 <- sum(mods[inds])/sum(mods)
  mods[is.na(mods)] <- 0
  if (sd(vm)==0)
    pow625 <- 0
  return(pow625)
}

dom.freq <- function(vm,M)
{
  if(length(vm)==1)
    return(NA)
  mods <- Mod(fft(vm))
  mods <- mods[-1]
  n <- length(mods)
  n <- floor(n/2)
  freq <- (M*(1:n))/(2*n)
  mods <- mods[1:n]
  dom.ind <- which.max(mods)
  d.f <- as.vector(freq[which.max(mods)])
  return(d.f)
}

frac.pow.dom.freq <- function(vm,M)
{
  mods <- Mod(fft(vm))
  mods <- mods[-1]
  n <- length(mods)
  n <- floor(n/2)
  freq <- (M*(1:n))/(2*n)
  mods <- mods[1:n]
  rat <- max(mods)/sum(mods)
  mods[is.na(mods)] <- 0
  if (sd(vm)==0)
    rat <- 0
  return(rat)
}
## loading raw timmeseries data 
load("D:/peaks_project/ankle_data.rda")
features_ankle_16sec<-array(dim=c(47034,20))
features_ankle_16sec<-as.data.frame(features_ankle_16sec)
deleting<-which(is.na(ankle_data$X_axis)==TRUE)
ankle_data$VM<-(ankle_data$X_axis^2+ankle_data$Y_axis^2+ankle_data$Z_axis^2)^0.5
ankle_data<-cbind(ankle_data[,c(1:3,7,4:6)])
#ankle_data<-ankle_data[-deleting,]
t=0

## building the features
for(i in 1:length(unique(ankle_data$ID))){  
  print(i);
  subset_1<-ankle_data[which(ankle_data$ID == unique(ankle_data$ID)[i]),]
  for(j in 1:length(unique(subset_1$label))){
    subset_2<-subset_1[which(subset_1$label == unique(subset_1$label)[j]),]
    if(nrow(subset_2)>=400){
      for(s in 1:(floor(nrow(subset_2)/25/16))){ # 25Hz is resolution of data, 16 sec is window-length
        features_ankle_16sec[s+t,1]<-subset_2$ID[1]
        features_ankle_16sec[s+t,2]<-subset_2$label[1]
        features_ankle_16sec[s+t,3]<-dom.freq(subset_2$VM[((s-1)*400+1):(s*400)], M=25)# 25 is the resolution (25Hz data)
        features_ankle_16sec[s+t,4]<-pow.625(subset_2$VM[((s-1)*400+1):(s*400)], 25)
        features_ankle_16sec[s+t,5]<-frac.pow.dom.freq(subset_2$VM[((s-1)*400+1):(s*400)], 25)
        features_ankle_16sec[s+t,6]<-mean(subset_2$VM[((s-1)*400+1):(s*400)])
        features_ankle_16sec[s+t,7]<-sd(subset_2$VM[((s-1)*400+1):(s*400)])
        features_ankle_16sec[s+t,8]<-mean(90*asin((subset_2$X_axis[((s-1)*400+1):(s*400)]/subset_2$VM[((s-1)*400+1):(s*400)])/(pi/2)), na.rm = TRUE)
        features_ankle_16sec[s+t,9]<-sd(90*asin((subset_2$X_axis[((s-1)*400+1):(s*400)]/subset_2$VM[((s-1)*400+1):(s*400)])/(pi/2)), na.rm = TRUE)
        features_ankle_16sec[s+t,10]<-var(subset_2$VM[((s-1)*400+1):(s*400)], na.rm = TRUE)
        del<-acf(subset_2$VM[((s-1)*400+1):(s*400)], lag.max = 1, type = "covariance", plot = FALSE)
        features_ankle_16sec[s+t,11]<-del$acf[2]
        features_ankle_16sec[s+t,12]<-skewness(subset_2$VM[((s-1)*400+1):(s*400)], na.rm = TRUE)
        features_ankle_16sec[s+t,13]<-kurtosis(subset_2$VM[((s-1)*400+1):(s*400)], na.rm = TRUE)
        features_ankle_16sec[s+t,14]<-geary(subset_2$VM[((s-1)*400+1):(s*400)], na.rm = TRUE)
        features_ankle_16sec[s+t,15]<-approx_entropy(subset_2$VM[((s-1)*400+1):(s*400)])
        #features_ankle_16sec[s+t,14]<-spectrum(subset_2$VM[((s-1)*400+1):(s*400)])
        features_ankle_16sec[s+t,16]<-cv(subset_2$VM[((s-1)*400+1):(s*400)])
        features_ankle_16sec[s+t,17]<-cor(subset_2$Y_axis[((s-1)*400+1):(s*400)],
                                                     subset_2$X_axis[((s-1)*400+1):(s*400)])
        features_ankle_16sec[s+t,18]<-cor(subset_2$Y_axis[((s-1)*400+1):(s*400)],
                                                     subset_2$Z_axis[((s-1)*400+1):(s*400)])
        features_ankle_16sec[s+t,19]<-cor(subset_2$X_axis[((s-1)*400+1):(s*400)],
                                                     subset_2$Z_axis[((s-1)*400+1):(s*400)])
        features_ankle_16sec[s+t,20]<-mean(90*asin((subset_2$Y_axis[((s-1)*400+1):(s*400)]/subset_2$VM[((s-1)*400+1):(s*400)])/(pi/2)), na.rm = TRUE)
        features_ankle_16sec[s+t,21]<-sd(90*asin((subset_2$Y_axis[((s-1)*400+1):(s*400)]/subset_2$VM[((s-1)*400+1):(s*400)])/(pi/2)), na.rm = TRUE)
        features_ankle_16sec[s+t,22]<-mean(90*asin((subset_2$Z_axis[((s-1)*400+1):(s*400)]/subset_2$VM[((s-1)*400+1):(s*400)])/(pi/2)), na.rm = TRUE)
        features_ankle_16sec[s+t,23]<-sd(90*asin((subset_2$Z_axis[((s-1)*400+1):(s*400)]/subset_2$VM[((s-1)*400+1):(s*400)])/(pi/2)), na.rm = TRUE)
        
      }
      t=t+s
    }
  }
}

save(file = "D:/features_ankle_16sec_new2.rda", features_ankle_16sec)
## started 9/4/2018 5:25 pm
print(Sys.time())

