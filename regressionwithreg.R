source("functions.R")
set.seed(1)
S0<-1
T<-1
K<-1
vol=0.2

TruePrice<-ZeroRateBachelierCall(S0,T,K,vol)$Price

N<-10000
Powers<-0:6
M<-length(Powers)
CoefMatrix<-matrix(0,nrow=length(Powers),1)
Sim1<-rnorm(N,0,1)
Sim2<-rnorm(N,0,1)

S1<-S0+vol*sqrt(T)*Sim1
S2<-S1+vol*sqrt(T)*Sim2
CallPayoff<-pmax(S2-K,0)
CallDelta<-S2>K
X1<-S1^0; 
X2<-rep(0, length(S1))
  
dummy<-sd(CallPayoff)/sd(CallDelta)
w<-1/(1+dummy)

for (k in 2:length(Powers)) {
  X1<-cbind(X1,S1^Powers[k])
  X2<-cbind(X2,Powers[k]*S1^(Powers[k]-1))
}
loss <- function(beta) {
  b = beta
  return(w*(t(X1%*%b-CallPayoff))%*%(X1%*%b-CallPayoff)+(1-w)*(t(X2%*%b-CallDelta))%*%(X2%*%b-CallDelta))
}
Grad <- function(beta) {
  b = beta
  return(2*w*((t(X1)%*%X1)%*%b)-2*w*(t(X1)%*%CallPayoff)+2*(1-w)*((t(X2)%*%X2)%*%b)-2*w*(t(X2)%*%CallDelta))
}

func1 <- function(beta){
  out <- loss(beta)
  attr(out, 'gradient') <- Grad(beta)
  return(out)
}

Coef = optim(c(0,0,0,0,0), func1, method = "BFGS")

Srange<-EstPricesimple<-EstDeltasimple<-SimEstDelta<-10:200/100

for (i in 1:length(Srange)){
  EstPricesimple[i]<-t(Coef$par)%*%Srange[i]^Powers
  EstDeltasimple[i]<-t(Coef$par[2:length(Powers)])%*% ( Powers[2:length(Powers)]*(Srange[i]^(Powers[2:length(Powers)]-1)))
}

tst<-ZeroRateBachelierCall(Srange,T,K,vol)

##First plot
xrange<-c(0.8,1.2)
par(mfrow=c(1,1))
#jpeg("Fig1.jpg", units="in", width=5, height=5, res=300)
plot(S1,CallPayoff,col="light grey", ylim=c(-0.05,0.5),xlim=xrange,xlab="Stock price", ylab="Call  value")
points(Srange,EstPricesimple,type='l',col="red", lwd = 2.2)
points(Srange,tst$Price,type='l', lwd = 2.2)
legend("topleft" , legend=c("Simulated payoff", "3th degree poly", "True"),
       col=c("light grey", "red","black" ),lty=1,  cex = 0.75)
dev.off()

xrange<-c(min(S1),max(S1))
par(mfrow=c(1,1))
#jpeg("Fig1.jpg", units="in", width=5, height=5, res=300)
plot(S1,CallPayoff,col="light grey", ylim=c(-0.05,1.0),xlim=xrange,xlab="Stock price", ylab="Call  value")
points(Srange,EstPricesimple,type='l',col="red", lwd = 2.2)
points(Srange,tst$Price,type='l', lwd = 2.2)
legend("topleft" , legend=c("Simulated payoff", "6th degree poly", "True"),
       col=c("light grey", "red","black" ),lty=1,  cex = 0.75)



plot(Srange,tst$Delta,type='l',col="black",lwd = 2.2, xlim=xrange,ylim=c(0,1.2), xlab="Stock price", ylab="Call Delta")
points(S1,CallDelta,col="light grey")
points(Srange,EstDeltasimple,type='l',col="red", lwd = 2.2)
legend("topleft" , legend=c("Simulated deltas", "6th degree poly", "True"),
       col=c("light grey", "red","black" ),lty=1,  cex = 0.65)



##Plotting with different weights in 6th degree poly
set.seed(1)
S0<-1
T<-1
K<-1
vol=0.2

TruePrice<-ZeroRateBachelierCall(S0,T,K,vol)$Price

N<-10000
Powers<-0:5
M<-length(Powers)
CoefMatrix<-matrix(0,nrow=length(Powers),1)
Sim1<-rnorm(N,0,1)
Sim2<-rnorm(N,0,1)

S1<-S0+vol*sqrt(T)*Sim1
S2<-S1+vol*sqrt(T)*Sim2
CallPayoff<-pmax(S2-K,0)
CallDelta<-S2>K
X1<-S1^0; 
X2<-rep(0, length(S1))

dummy<-sd(CallPayoff)/sd(CallDelta)
w<-0.5

for (k in 2:length(Powers)) {
  X1<-cbind(X1,S1^Powers[k])
  X2<-cbind(X2,Powers[k]*S1^(Powers[k]-1))
}
loss <- function(beta) {
  b = beta
  return(w*(t(X1%*%b-CallPayoff))%*%(X1%*%b-CallPayoff)+(1-w)*(t(X2%*%b-CallDelta))%*%(X2%*%b-CallDelta))
}
Grad <- function(beta) {
  b = beta
  return(2*w*((t(X1)%*%X1)%*%b)-2*w*(t(X1)%*%CallPayoff)+2*(1-w)*((t(X2)%*%X2)%*%b)-2*w*(t(X2)%*%CallDelta))
}

func1 <- function(beta){
  out <- loss(beta)
  attr(out, 'gradient') <- Grad(beta)
  return(out)
}

Coef = optim(c(0,0,0,0,0,0), func1, method = "BFGS")

Srange<-EstPricesimple50<-EstDeltasimple50<-10:200/100

for (i in 1:length(Srange)){
  EstPricesimple50[i]<-t(Coef$par)%*%Srange[i]^Powers
  EstDeltasimple50[i]<-t(Coef$par[2:length(Powers)])%*% ( Powers[2:length(Powers)]*(Srange[i]^(Powers[2:length(Powers)]-1)))
}


##Delta only
set.seed(1)
S0<-1
T<-1
K<-1
vol=0.2

TruePrice<-ZeroRateBachelierCall(S0,T,K,vol)$Price

N<-10000
Powers<-0:5
M<-length(Powers)
CoefMatrix<-matrix(0,nrow=length(Powers),1)
Sim1<-rnorm(N,0,1)
Sim2<-rnorm(N,0,1)

S1<-S0+vol*sqrt(T)*Sim1
S2<-S1+vol*sqrt(T)*Sim2
CallPayoff<-pmax(S2-K,0)
CallDelta<-S2>K
X1<-S1^0; 
X2<-rep(0, length(S1))

dummy<-sd(CallPayoff)/sd(CallDelta)
w<-0.01

for (k in 2:length(Powers)) {
  X1<-cbind(X1,S1^Powers[k])
  X2<-cbind(X2,Powers[k]*S1^(Powers[k]-1))
}
loss <- function(beta) {
  b = beta
  return(w*(t(X1%*%b-CallPayoff))%*%(X1%*%b-CallPayoff)+(1-w)*(t(X2%*%b-CallDelta))%*%(X2%*%b-CallDelta))
}
Grad <- function(beta) {
  b = beta
  return(2*w*((t(X1)%*%X1)%*%b)-2*w*(t(X1)%*%CallPayoff)+2*(1-w)*((t(X2)%*%X2)%*%b)-2*w*(t(X2)%*%CallDelta))
}

func1 <- function(beta){
  out <- loss(beta)
  attr(out, 'gradient') <- Grad(beta)
  return(out)
}

Coef = optim(c(0,0,0,0,0,0), func1, method = "BFGS")

Srange<-EstPricesimpleonly<-EstDeltasimpleonly<-10:200/100

for (i in 1:length(Srange)){
  EstPricesimpleonly[i]<-t(Coef$par)%*%Srange[i]^Powers
  EstDeltasimpleonly[i]<-t(Coef$par[2:length(Powers)])%*% ( Powers[2:length(Powers)]*(Srange[i]^(Powers[2:length(Powers)]-1)))
}


#plot
xrange<-c(min(S1),max(S1))
par(mfrow=c(1,1))
#jpeg("Fig1.jpg", units="in", width=5, height=5, res=300)
plot(S1,CallPayoff,col="light grey", ylim=c(-0.05,1.0),xlim=xrange,xlab="Stock price", ylab="Call  value")
points(Srange,EstPricesimple,type='l',col="red", lwd = 2.2)
points(Srange,tst$Price,type='l', lwd = 2.2)
points(Srange,EstPricesimple50,type='l',col="green", lwd = 2.2)
points(Srange,EstPricesimpleonly,type='l',col="blue", lwd = 2.2)
legend("topleft" , legend=c("Simulated payoff", "Price only","true", "50/50", "Delta(almost)"),
       col=c("light grey", "red","black","yellow","blue" ),lty=1,  cex = 0.75)


plot(Srange,tst$Delta,type='l',col="black",lwd = 2.2, xlim=xrange,ylim=c(0,1.2), xlab="Stock price", ylab="Call Delta")
points(S1,CallDelta,col="light grey")
points(Srange,EstDeltasimple,type='l',col="red", lwd = 2.2)
points(Srange,EstDeltasimple50,type='l',col="green", lwd = 2.2)
points(Srange,EstDeltasimpleonly,type='l',col="blue", lwd = 2.2)
legend("topleft" , legend=c("Simulated payoff", "Price only","true", "50/50", "Delta(almost)"),
       col=c("light grey", "red","black","yellow","blue" ),lty=1,  cex = 0.65)



1/(1+(1/2))
2/(2+1)
