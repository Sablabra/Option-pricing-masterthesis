##Hedge expirement using the SoftPlus function. First by setting S = seq(0.10,1.9999,c) giving us N=19000. Change the order is done in line 77. 
#We also used regular simulation for the hedge expirement, where one can change order in line 233, and the number of simulation in line 171
source("functions")
set.seed(1)
S0 = 1
T = 1
K = 1
vol = 0.2
c = 0.0001
S = seq(0.10,1.9999,c)
TruePrice <- ZeroRateBachelierCall(S,T,K,vol)$Price

Nhedge = 52
dt = T/Nhedge

N<-10000
Npaths <- length(S)
Sim1<-rnorm(Npaths,0,1)
Sim2<-rnorm(Npaths,0,1)

Criterion<-function(a) {
  Dummy = a[1]
  dif_dummy = 0
  for (i in 1:((length(a)-1)/3)){
    Dummy<-Dummy+a[3*i-1]*SoftPlus(a[3*i]*S + a[3*i+1])
    dif_dummy = dif_dummy + a[3*i-1]*dXSoftPlus(a[3*i]*S+a[3*i+1])*a[3*i]
  }
  
  return(w*sum((CallPrice-Dummy)^2)+(1-w)*sum((CallDelta-dif_dummy)^2))
}

Gradient <- function(a){
  Dummy = a[1]
  dif_dummy = 0
  for (j in 1:((length(a)-1)/3)){
    Dummy = Dummy+a[3*j-2]*SoftPlus(a[3*j-1]*S + a[3*j])
    dif_dummy = dif_dummy+a[3*j-2]*dXSoftPlus(a[3*j-1]*S + a[3*j])*a[3*j-1]
    Dummy = Dummy-CallPrice
    dif_dummy = dif_dummy-CallDelta
    Grad = matrix(0,length(a),1)
    Grad[1] = w*2*sum(Dummy)
  }
  for (j in 1:((length(a)-1)/3)){
    d1 = SoftPlus(a[3*j-1]*S + a[3*j])
    dif_d1 = dXSoftPlus(a[3*j-1]*S+a[3*j])*a[3*j-1]
    Grad[3*j-1] = w*2*sum(Dummy*d1)+(1-w)*2*sum(dif_dummy*dif_d1)
    
    d2 = a[3*j-2]*dXSoftPlus(a[3*j-1]*S + a[3*j])*S
    dif_d2 = a[3*j-2]*dXSoftPlus2(a[3*j-1]*S+a[3*j])*S*a[3*j-1]+a[3*j-2]*dXSoftPlus(a[3*j-1]*S+a[3*j])
    Grad[3*j] = w*2*sum(Dummy*d2)+(1-w)*2*sum(dif_dummy*dif_d2)
    
    d3 = a[3*j-2]*dXSoftPlus(a[3*j-1]*S + a[3*j])
    dif_d3 = a[3*j-2]*dXSoftPlus2(a[3*j-1]*S+a[3*j])*a[3*j-1]
    Grad[3*j+1] = w*2*sum(Dummy*d3)+(1-w)*2*sum(dif_dummy*dif_d3)
  }
  return(Grad)
}


if (length(x) == 4){
  y = Regression_SP(Coef$par,S)
  delta = Coef$par[2]*dXSoftPlus(Coef$par[3]*S+Coef$par[4])*Coef$par[3]
} else if (length(x)==7){
  y = Regression_SP2(Coef$par,S)
  delta = Coef$par[2]*dXSoftPlus(Coef$par[3]*S+Coef$par[4])*Coef$par[3]+Coef$par[5]*dXSoftPlus(Coef$par[6]*S+Coef$par[7])*Coef$par[6]
} else if (length(x)==10){
  y = Regression_SP3(Coef$par,S)
  delta = Coef$par[2]*dXSoftPlus(Coef$par[3]*S+Coef$par[4])*Coef$par[3]+Coef$par[5]*dXSoftPlus(Coef$par[6]*S+Coef$par[7])*Coef$par[6]+Coef$par[8]*dXSoftPlus(Coef$par[9]*S+Coef$par[10])*Coef$par[9]
} else if (length(x)==13){
  y = Regression_SP4(Coef$par,S)
  delta = Coef$par[2]*dXSoftPlus(Coef$par[3]*S+Coef$par[4])*Coef$par[3]+Coef$par[5]*dXSoftPlus(Coef$par[6]*S+Coef$par[7])*Coef$par[6]+Coef$par[8]*dXSoftPlus(Coef$par[9]*S+Coef$par[10])*Coef$par[9]+Coef$par[11]*dXSoftPlus(Coef$par[12]*S+Coef$par[13])*Coef$par[12]
} else if (length(x)==16){
  y = Regression_SP5(Coef$par,S)
  delta = Coef$par[2]*dXSoftPlus(Coef$par[3]*S+Coef$par[4])*Coef$par[3]+Coef$par[5]*dXSoftPlus(Coef$par[6]*S+Coef$par[7])*Coef$par[6]+Coef$par[8]*dXSoftPlus(Coef$par[9]*S+Coef$par[10])*Coef$par[9]+Coef$par[11]*dXSoftPlus(Coef$par[12]*S+Coef$par[13])*Coef$par[12]+Coef$par[14]*dXSoftPlus(Coef$par[15]*S+Coef$par[16])*Coef$par[15]
}


order = 1
errorHedge = matrix(0,order,1)
for (h in 1:order){
  if (h==1){
    x = c(0,0,0,0)
  } else {
    z = runif(3)
    x <- c(Coef$par, z)
  }
  Coef = matrix(0,length(x),52)
  S0 = 1
  for (i in 1:52){
    S = seq(0.25,1.9999,c)
    S2 = S+vol*sqrt(T-(i-1)*dt)*Sim2
    CallPrice = pmax(S2-K,0)
    CallDelta = matrix(0,length(S),1)
    for (y in 1:length(S)){
      if (S2[y]>=K){
        CallDelta[y]=1
      } else {
        CallDelta[y]=0
      }
    }
    dummy = sd(CallPrice)/sd(CallDelta)
    w = 1/(1+dummy)
    
    func <- function(a){
      out <- Criterion(a)
      attr(out, 'gradient') <- Gradient(a)
      return(out)
    }
    Coef = optim(x, func, method = "BFGS")
    x = Coef$par
    Coef[,i] = x 
  }
  Initialprice = ZeroRateBachelierCall(S0,T,K,vol)$Price
  Initialdelta = ZeroRateBachelierCall(S0,T,K,vol)$Delta
  start = matrix(0,2,1)
  start[1] = Regression_SP(Coef$par,S0)
  start[2] = Coef$par[2]*dXSoftPlus(Coef$par[3]*S0+Coef$par[4])*Coef$par[3]
  Stock = matrix(1,N,1)
  Vpf = rep(Initialprice,N)
  a = rep(start[2],N)
  b = Vpf-a*Stock
  
  error = matrix(1,N,Nhedge)
  error[,1] = a*Stock+b-Initialprice
  S = Stock
  for (i in 2:(Nhedge-1)) {
    S<-S+vol*sqrt(dt)*rnorm(N)
    dummy<-T-dt*(i-1)
    dummy<-ZeroRateBachelierCall(S,dummy,K,vol)
    for (k in 1:N){
      dummy$Delta[k] = 0
      for (j in 1:((length(Coef$par)-1)/3)){
        dummy$Delta[k]<-dummy$Delta[k] + Coef[3*j-1,i]*dXSoftPlus(Coef[3*j,i]*S[k]+Coef[3*j+1,i])*Coef[3*j,i]
      }
    }
    Vpf<-a*S+b
    a<-dummy$Delta
    b<-Vpf-a*S
  }
}

S<-S+vol*sqrt(dt)*rnorm(N)
CallPayoff<-pmax(S-K,0)
Vpf<-a*S+b

par(mfrow=c(1,1))

plot(S,Vpf,col="light grey", ylim = c(0,1))
points(S,CallPayoff,lty=2, col = "blue")
legend("topleft" , legend=c("Hedge portfolio", "PayOff Call Option"),
       col=c("light grey", "blue"),lty=1,  cex = 0.75)

z<-sd(CallPayoff-Vpf)/Initialprice

print(round(c(z,1.96*sqrt(2)*z/sqrt(Npaths)),8))


##Hedge for softplus simulation
source("functions")
S0 = 1
T = 1
K = 1
vol = 0.2

TruePrice <- ZeroRateBachelierCall(S0,T,K,vol)$Price

Nhedge = 52
dt = T/Nhedge

set.seed(1)
N<-10000
Npaths <- N
Sim1<-rnorm(Npaths,0,1)
Sim2<-rnorm(Npaths,0,1)

Criterion<-function(a) {
  Dummy = a[1]
  dif_dummy = 0
  for (i in 1:((length(a)-1)/3)){
    Dummy<-Dummy+a[3*i-1]*SoftPlus(a[3*i]*S + a[3*i+1])
    dif_dummy = dif_dummy + a[3*i-1]*dXSoftPlus(a[3*i]*S+a[3*i+1])*a[3*i]
  }
  
  return(w*sum((CallPrice-Dummy)^2)+(1-w)*sum((CallDelta-dif_dummy)^2))
}

Gradient <- function(a){
  Dummy = a[1]
  dif_dummy = 0
  for (j in 1:((length(a)-1)/3)){
    Dummy = Dummy+a[3*j-2]*SoftPlus(a[3*j-1]*S + a[3*j])
    dif_dummy = dif_dummy+a[3*j-2]*dXSoftPlus(a[3*j-1]*S + a[3*j])*a[3*j-1]
    Dummy = Dummy-CallPrice
    dif_dummy = dif_dummy-CallDelta
    Grad = matrix(0,length(a),1)
    Grad[1] = w*2*sum(Dummy)
  }
  for (j in 1:((length(a)-1)/3)){
    d1 = SoftPlus(a[3*j-1]*S + a[3*j])
    dif_d1 = dXSoftPlus(a[3*j-1]*S+a[3*j])*a[3*j-1]
    Grad[3*j-1] = w*2*sum(Dummy*d1)+(1-w)*2*sum(dif_dummy*dif_d1)
    
    d2 = a[3*j-2]*dXSoftPlus(a[3*j-1]*S + a[3*j])*S
    dif_d2 = a[3*j-2]*dXSoftPlus2(a[3*j-1]*S+a[3*j])*S*a[3*j-1]+a[3*j-2]*dXSoftPlus(a[3*j-1]*S+a[3*j])
    Grad[3*j] = w*2*sum(Dummy*d2)+(1-w)*2*sum(dif_dummy*dif_d2)
    
    d3 = a[3*j-2]*dXSoftPlus(a[3*j-1]*S + a[3*j])
    dif_d3 = a[3*j-2]*dXSoftPlus2(a[3*j-1]*S+a[3*j])*a[3*j-1]
    Grad[3*j+1] = w*2*sum(Dummy*d3)+(1-w)*2*sum(dif_dummy*dif_d3)
  }
  return(Grad)
}


if (length(x) == 4){
  y = Regression_SP(Coef$par,S)
  delta = Coef$par[2]*dXSoftPlus(Coef$par[3]*S+Coef$par[4])*Coef$par[3]
} else if (length(x)==7){
  y = Regression_SP2(Coef$par,S)
  delta = Coef$par[2]*dXSoftPlus(Coef$par[3]*S+Coef$par[4])*Coef$par[3]+Coef$par[5]*dXSoftPlus(Coef$par[6]*S+Coef$par[7])*Coef$par[6]
} else if (length(x)==10){
  y = Regression_SP3(Coef$par,S)
  delta = Coef$par[2]*dXSoftPlus(Coef$par[3]*S+Coef$par[4])*Coef$par[3]+Coef$par[5]*dXSoftPlus(Coef$par[6]*S+Coef$par[7])*Coef$par[6]+Coef$par[8]*dXSoftPlus(Coef$par[9]*S+Coef$par[10])*Coef$par[9]
} else if (length(x)==13){
  y = Regression_SP4(Coef$par,S)
  delta = Coef$par[2]*dXSoftPlus(Coef$par[3]*S+Coef$par[4])*Coef$par[3]+Coef$par[5]*dXSoftPlus(Coef$par[6]*S+Coef$par[7])*Coef$par[6]+Coef$par[8]*dXSoftPlus(Coef$par[9]*S+Coef$par[10])*Coef$par[9]+Coef$par[11]*dXSoftPlus(Coef$par[12]*S+Coef$par[13])*Coef$par[12]
} else if (length(x)==16){
  y = Regression_SP5(Coef$par,S)
  delta = Coef$par[2]*dXSoftPlus(Coef$par[3]*S+Coef$par[4])*Coef$par[3]+Coef$par[5]*dXSoftPlus(Coef$par[6]*S+Coef$par[7])*Coef$par[6]+Coef$par[8]*dXSoftPlus(Coef$par[9]*S+Coef$par[10])*Coef$par[9]+Coef$par[11]*dXSoftPlus(Coef$par[12]*S+Coef$par[13])*Coef$par[12]+Coef$par[14]*dXSoftPlus(Coef$par[15]*S+Coef$par[16])*Coef$par[15]
}


order = 2
errorHedge = matrix(0,order,1)
for (h in 1:order){
  if (h==1){
    x = c(0,0,0,0)
  } else {
    z = runif(3)
    x <- c(Coef$par, z)
  }
  Coef = matrix(0,length(x),52)
  S0 = 1
  for (i in 1:52){
    S1 = S0+vol*sqrt(T)*Sim1
    S = S1
    S2 = S+vol*sqrt(T-(i-1)*dt)*Sim2
    CallPrice = pmax(S2-K,0)
    CallDelta = matrix(0,length(S),1)
    for (y in 1:length(S)){
      if (S2[y]>=K){
        CallDelta[y]=1
      } else {
        CallDelta[y]=0
      }
    }
    dummy = sd(CallPrice)/sd(CallDelta)
    w = 1/(1+dummy)
    
    func <- function(a){
      out <- Criterion(a)
      attr(out, 'gradient') <- Gradient(a)
      return(out)
    }
    Coef = optim(x, func, method = "BFGS")
    x = Coef$par
    print(x)
    Coef[,i] = x 
  }
  Initialprice = ZeroRateBachelierCall(S0,T,K,vol)$Price
  Initialdelta = ZeroRateBachelierCall(S0,T,K,vol)$Delta
  start = matrix(0,2,1)
  start[1] = Regression_SP2(Coef$par,S0)
  start[2] = Coef$par[2]*dXSoftPlus(Coef$par[3]*S0+Coef$par[4])*Coef$par[3]+Coef$par[5]*dXSoftPlus(Coef$par[6]*S0+Coef$par[7])*Coef$par[6]
  Stock = matrix(1,N,1)
  Vpf = rep(Initialprice,N)
  a = rep(start[2],N)
  b = Vpf-a*Stock
  
  error = matrix(1,N,Nhedge)
  error[,1] = a*Stock+b-Initialprice
  S = Stock
  for (i in 2:(Nhedge-1)) {
    S<-S+vol*sqrt(dt)*rnorm(N)
    dummy<-T-dt*(i-1)
    dummy<-ZeroRateBachelierCall(S,dummy,K,vol)
    for (k in 1:N){
      dummy$Delta[k] = 0
      for (j in 1:((length(Coef$par)-1)/3)){
        dummy$Delta[k]<-dummy$Delta[k] + Coef[3*j-1,i]*dXSoftPlus(Coef[3*j,i]*S[k]+Coef[3*j+1,i])*Coef[3*j,i]
      }
    }
    Vpf<-a*S+b
    a<-dummy$Delta
    b<-Vpf-a*S
  }
}

S<-S+vol*sqrt(dt)*rnorm(N)
CallPayoff<-pmax(S-K,0)
Vpf<-a*S+b

par(mfrow=c(1,1))

plot(S,Vpf,col="light grey", ylim = c(0,1))
points(S,CallPayoff,lty=2, col = "blue")
legend("topleft" , legend=c("Hedge portfolio", "PayOff Call Option"),
       col=c("light grey", "blue"),lty=1,  cex = 0.75)

z<-sd(CallPayoff-Vpf)/Initialprice

print(round(c(z,1.96*sqrt(2)*z/sqrt(Npaths)),8))

#Softplus hedge
simula <- c(500,1000,2000,3000,4000,5000,6000)

ll1 = c(0.21132416,0.21232809,0.14378216,0.12596808,0.12307295,0.12162969,0.12330648 )
mm1 = c(0.13727797,0.12617370,0.12648002,0.12164828, 0.12127743,0.12101818, 0.12450184)
ll2 = c(0.36478379,0.20930125,0.23645605,0.15068388,0.14006135, 0.13535888, 0.15387726)
mm2 = c(0.16436148,0.13207614, 0.1289644,0.12443212,0.12254308,0.12553841,0.12712284)
ll3 = c(0.37869920,0.23021317,0.33990492,0.23110745,0.13505870,0.14301045,0.18472897)
mm3 = c(0.1604251,0.13072361 ,0.13506467 ,0.124202468 ,0.12262579 ,0.12586261,0.12766549)

plot(simula,ll1,ylim=c(0,0.5),type = "l",col="red",xlab=c("Number of simulations"),ylab=c("SDE of relative hedge error in SoftPlus"))
lines(simula,mm1,col="blue",type = "b")
lines(simula,ll2,col="red",type = "l")
lines(simula,mm2,col="blue",type = "b")
lines(simula,ll3,col="red",type="l")
lines(simula,mm3,col="blue",type = "b")
abline(h=0.1208,col="orange")
legend("topleft" , legend=c("Price only-regression", "Equal variance","True Delta"),
       col=c("red", "blue","orange"),lty=1,  cex = 0.65)
text(5975,0.135,"1st",cex=0.6)
text(5950,0.163,"2nd",cex=0.6)
text(5950,0.212,"3rd",cex=0.6)
