

#This script contains the functions used in order to fit the Imprecise-response regression model
#proposed in the 2011 paper "A determination coefficient for a linear
#regression model with imprecise response" written by Ferraro, Ana Colubi, González-Rodríguez and Coppi.

#extension of the Yang-Ko metric, used in order to embed the F_lr space into R^3.
#In this version, the two parameters lambda and rho are both equal to 1/2

Dlp<-function(X,Y,lambda=1/2,rho=1/2){
  return ((X[1]-Y[1])^2 + ((X[1]-lambda*X[2])-(Y[1]-lambda*Y[2]))^2 + ((X[1]+rho*X[3])-(Y[1]+rho*Y[3]))^2)
}

#The log-linear model used in the imprecise-response model.
logFit<-function(Y,X){
  return(lm(log(abs(Y))~X))
}

#linear model.
myFit<-function(Y,X){
  return(lm(Y~X))
}

#function that computes the Tn statistic which is the estimate of the 
#determination coefficient scaled by n
#X is the explanatory variable, Yl and Yr are respectively the left and right spread
#and Ym is the center

Tn<-function(X,Ym,Yl,Yr){
  n<-length(X)
  
  #fitting the regression model
  YmFit<-predict(myFit(Ym,X))
  YlFit<-predict(logFit(abs(Yl),X))
  YrFit<-predict(logFit(abs(Yr),X))
  
  #here we put the abs() in order to prevent problems with negative values 
  #(log function)
  Y<-cbind(Ym,log(abs(Yl)),log(abs(Yr)))
  Yh<-cbind(YmFit,YlFit,YrFit) #fitted values
  Ybar<-c(mean(Ym),mean(log(abs(Yl))),mean(log(abs(Yr))))
  
  y1<-0
  y2<-0
  for (i in 1:n){
  y1<-y1+Dlp(Yh[i,],Ybar)
  y2<-y2+Dlp(Y[i,],Ybar)
  }
  Tn=n*(y1/y2)
  return (Tn)
  }


#function that computes the bootstrap replicates of the n*R^2 statistic
#ix is the vector indices the bootstrap will be based on
Tnstar<-function(X,Ym,Yl,Yr,ix){
  n<-length(X)

  #we fit again the regression model
  YmFit<-predict(myFit(Ym,X))
  YlFit<-predict(logFit(abs(Yl),X))
  YrFit<-predict(logFit(abs(Yr),X))

  #compute the residuals (Z)
  Zm<-Ym-YmFit
  Zl<-Yl-YlFit
  Zr<-Yr-YrFit
  
  #Here comes the bootstrap: we resample the observed data
  X<-X[ix]
  Zm<-Zm[ix]
  Zl<-Zl[ix]
  Zr<-Zr[ix]

  Z<-cbind(Zm,Zl,Zr) #bootstrapped residuals
  
  #model fitting of the residuals
  Zmh<-predict(myFit(Zm,X))
  Zlh<-predict(myFit(Zl,X))
  Zrh<-predict(myFit(Zr,X))
  Zh<-cbind(Zmh,Zlh,Zrh) #Zhat
  
  Zem<-c(mean(Zm),mean(Zl),mean(Zr))
  
  y1<-0
  y2<-0
  for (i in 1:n){
    y1<-y1+Dlp(Zh[i,],Zem)
    y2<-y2+Dlp(Z[i,],Zem)
  }
  Tstar<-n*(y1/y2)
  return (Tstar)
  }
  







