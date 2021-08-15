##Performing a Hedge expirement, using the Polynomial regression. The degree of the polynomial can be changed in line 19
source("functions.r")
S0<-1
T<-1
K<-1.0
vol=0.2

TruePrice<-ZeroRateBachelierCall(S0,T,K,vol)$Price

Npaths<-1000
Nhedge <-52
dt<-T/Nhedge 

Poly<-TRUE

if (Poly){
  # Polynomial Regression Estimation
  set.seed(1)
  N<-1000
  Powers<-0:8
  M<-length(Powers)
  CoefMatrix<-matrix(0,nrow=length(Powers),ncol=Nhedge)
  Sim1<-rnorm(Npaths,0,1)
  Sim2<-rnorm(N,0,1)
  
  for (i in 1:Nhedge){
    
    S1<-S0+vol*sqrt(T)*Sim1
    S2<-S1+vol*sqrt(T-(i-1)*dt)*Sim2
    CallPayoff<-pmax(S2-K,0)
    CallDelta<-S2>K
    X1<-S1^0; 
    X2<-rep(0, length(S1))
    
    dummy<-sd(CallPayoff)/sd(CallDelta)
    w<-0.5
    print(c(i,round(w,2)))
    
    for (k in 2:length(Powers)) {
      X1<-cbind(X1,S1^Powers[k])
      X2<-cbind(X2,Powers[k]*S1^(Powers[k]-1))
    }
    
    OLSCoef<-solve((w*t(X1)%*%X1+(1-w)*t(X2)%*%X2),(w*t(X1)%*%CallPayoff + (1-w)*t(X2)%*%CallDelta))
    
    Criterion<-function(a) {
      M<-length(a)
      N<-length(S1)
      eps1<-eps2<-rep(0,N)
      for (i in 1:N){
        eps1[i]<-CallPayoff[i]-a%*%(S1[i]^Powers)
        eps2[i]<-CallDelta[i]-a[2:M]%*%(Powers[2:M]*S1[i]^((Powers[2:M])-1))
      }
      
      Criterion<-w*sum(eps1^2)+(1-w)*sum(eps2^2)
    }
    #NumCoef<-nlm(f=Criterion,p=OLSCoef)$estimate
    CoefMatrix[,i]<-NumCoef<-OLSCoef
  }
}



start<-ZeroRateBachelierCall(S0,T,K,vol)
if (Poly){
  start$Price<-CoefMatrix[,1]%*%(S0^Powers)
  start$Delta<-CoefMatrix[2:M,1]%*%(Powers[2:M]*S0^(Powers[2:M]-1))
} 

InitialPrice<-TruePrice
temp<-S<-rep(S0,Npaths)
Vpf<-rep(InitialPrice,Npaths)
a<-rep(start$Delta,Npaths)
b<-Vpf-a*S

set.seed(11)

for (i in 2:(Nhedge-1)) {
  S<-S+vol*sqrt(dt)*rnorm(Npaths)
  tau<-T-dt*(i-1)
  dummy<-ZeroRateBachelierCall(S,tau,K,vol)
  if (Poly){
    for (k in 1:Npaths)dummy$Delta[k]<-CoefMatrix[2:M,i]%*%(Powers[2:M]*S[k]^(Powers[2:M]-1))
    
  }
  Vpf<-a*S+b
  a<-dummy$Delta
  b<-Vpf-a*S
}

S<-S+vol*sqrt(dt)*rnorm(Npaths)
CallPayoff<-pmax(S-K,0)
Vpf<-a*S+b

#par(mfrow=c(1,1))

#plot(S,Vpf,col="light grey")
#points(S,CallPayoff,lty=2,pch=1, col = "blue")
#legend("topleft" , legend=c("Hedge portfolio", "PayOff Call Option"),
#       col=c("light grey", "blue"),lty=1,  cex = 0.75)

z<-sd(CallPayoff-Vpf)/InitialPrice

print(round(c(z,1.96*sqrt(2)*z/sqrt(Npaths)),8))





##plot
l1 = c(0.32517609,0.30200278,0.31745755, 0.29913710, 0.30960417, 0.29914881,0.30304256)
m1 = c(0.31577579,0.30828383, 0.30637829, 0.3048026,0.30677868, 0.30489957,0.30464292)
l2 = c(0.25955224,0.24598245, 0.23935104, 0.21646387,0.20855105 ,0.20919292, 0.21168448)
m2 = c(0.23398942 ,0.21705811, 0.21398887, 0.21147153, 0.21130118,0.21078038, 0.21056656)
l3 = c(0.25581098,0.2487103, 0.2532273,0.21749940, 0.20959071,  0.20896209,0.21151641)
m3 = c(0.22446275, 0.2152528 ,0.21467602,0.21127953, 0.21145284, 0.21065851, 0.21022381)
l4 = c(0.25222845,0.21594038,  0.18860597,0.17115013, 0.16764668,0.16645541,0.16579001)
m4 = c(0.18735005,0.17607270,0.17332253,0.17185867,0.17213401, 0.17228225, 0.17235648)
l5 = c(0.26515548, 0.23905409,0.18516271, 0.19142325, 0.17093196, 0.16626577, 0.16708303)
m5 = c(0.18448952, 0.17515806,0.17249577, 0.17077526,  0.17156632, 0.17224800,0.17231740)
l6 = c(0.29081515,0.22281876, 0.1894684, 0.17559399,0.15988587, 0.14811182, 0.14911261)
m6 = c(0.16974758,0.15741007, 0.15291336, 0.14902174 ,0.15036149, 0.14966914, 0.15035207)
l7 = c(0.23518632,0.21458137,0.16606963, 0.15091882, 0.15092098 ,0.14880402, 0.14958354)
m7 = c(0.21038779, 0.1971347,0.15917348,0.14875176,  0.15041104, 0.14990426, 0.15064045)


simula <- c(500,1000,2000,3000,4000,5000,6000)

plot(simula,l1,ylim=c(0,0.5),type = "l",col="red",xlab=c("Number of simulations"),ylab=c("Standard deviation of relative hedge error"))
lines(simula,m1,col="blue",type = "b")
lines(simula,l2,col="red",type = "l")
lines(simula,m2,col="blue",type = "b")
lines(simula,l3,col="red",type="l")
lines(simula,m3,col="blue",type = "b")
lines(simula,l4,col="red",type="l")
lines(simula,m4,col="blue",type = "b")
lines(simula,l5,col="red",type="l")
lines(simula,m5,col="blue",type = "b")
lines(simula,l6,col="red",type="l")
lines(simula,m6,col="blue",type = "b")
lines(simula,l7,col="red",type="l")
lines(simula,m7,col="blue",type = "b")
abline(h=0.1208,col="orange")
legend("topleft" , legend=c("Price only-regression", "Equal variance","True Delta"),
       col=c("red", "blue","orange"),lty=1,  cex = 0.65)
text(5975,0.33,"3rd",cex=0.6)
text(5950,0.23,"4th,5th",cex=0.6)
text(5950,0.1852,"6th,7th",cex=0.6)
text(5950,0.1362,"8th,9th",cex=0.6)
