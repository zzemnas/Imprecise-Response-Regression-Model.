
###################################################MAIN
###########I: APPLICATION TO REAL DATA USING THE TEST BASED ON THE BOOTSTRAP (trees quality data)
library("readxl")
myData<-as.data.frame(read_excel("Quality.xls")) #importing excel data
Ym<-myData[,1]
Yl<-myData[,2]
Yr<-myData[,3]
X<-myData[,4]

alpha<-0.01 #level of significance
n<-length(X)
B<-500

#compute the Tn statistics using the observed data
TnData<-as.numeric(Tn(X,Ym,Yl,Yr))
TnData

#initialization of the bootstrap replicates vector
TnB<-numeric(B)

#we could haves also used the boot function: proc.time() computes the time
#spent running this method to our data
#ptm <- proc.time()
#TnB<-boot(data=X,statistic=Tnstar,R=B,Ym=Ym,Yl=Yl,Yr=Yr)$t
#proc.time() - ptm

for (j in 1:B)
{
  TnB[j]<-Tnstar(X,Ym,Yl,Yr,ix=sample(1:n,replace=TRUE))
}

#The approximate p-value equals 0, thus we reject the hypothesis of linear independence without a doubt
p<-mean(TnB>TnData)
p<alpha
#######################################################


#############II: TEST BASED ON THE ASYMPTOTICAL ANALYSIS OF THE PAPER
#############THIS SIMULATION USES RANDOM VA (BOTH ERRORS AND EXPLANATORY)
rm(list=ls())
alpha<-0.1 #significance level
n<-50 #number of feature vectors
delta = seq(0,1,length.out=6)/sqrt(n) #sequence of pitman-like alternative hypothesis
a<-as.data.frame(rbind(delta,delta,delta),nrow=3,ncol=3) #matrix of delta

#m<-10000 #number of tests made
m<-100
#B<-1000 #number of bootstrap replicates inside each test
B<-500

TestB<-numeric(m) #results of the tests
TnB<-numeric(B) #bootstrap replicates

sim <- matrix(0, 6, 4) #matrix used to keep in memory the power of the tests made

for(k in 1:1) {
  #delta <- a[,k]  
  delta<-c(1,1,1)
  for (i in 1:m) {
    X<-rnorm(n)
    Ym<-delta[1]*X+rnorm(n)
    Yl<-delta[2]*X+rnorm(n)
    Yr<-delta[3]*X+rnorm(n)
    TnData<-Tn(X,Ym,Yl,Yr)
    for (j in 1:B){
      ix=sample(1:n,replace=TRUE)
      TnB[j]<-Tnstar(X,Ym,Yl,Yr,ix=ix)
    }
    p<-mean(TnB>TnData) #probability of observing a value that is more "extreme" than the one observed
    TestB[i] <- as.integer(p <= alpha)
  }
  print(c(delta, mean(TestB))) #alternative hypothesis + power of the test based on the delta-alt-hypothesis
  sim[k, ] <- c(delta, mean(TestB))
} 
#plot(delta,sim[,4])
#Results: delta=1 -> pow=0.99, delta=0.1 -> pow=0.12, per delta=0 -> pow=0.08

#https://nceas.github.io/oss-lessons/parallel-computing-in-r/parallel-computing-in-r.html
