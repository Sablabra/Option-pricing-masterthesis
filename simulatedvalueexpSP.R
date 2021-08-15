##Simulation study
set.seed(1)
S0 = 1
T = 1
K = 1
vol = 0.2

set.seed(1)
N<-10000
Npaths <- N
Sim1<-rnorm(Npaths,0,1)
Sim2<-rnorm(Npaths,0,1)

S1 = S0+vol*sqrt(T)*Sim1
S = S1
S2 = S+vol*sqrt(T)*Sim2
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

TruePrice <- ZeroRateBachelierCall(S,T,K,vol)$Price
TrueDelta <- ZeroRateBachelierCall(S,T,K,vol)$Delta

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
for (h in 1:order){
  if (h==1){
    x = c(0,0,0,0)
  } else {
    z = runif(3)
    x <- c(Coef$par, z)
  }
  Srange<-10:200/100
  Coef <- optim(x, func, method = "BFGS")
  print(Coef$par)
  y = Regression_SP(Coef$par,Srange)
  delta = 0
  for (j in 1:((length(Coef$par)-1)/3)){
    delta = delta + Coef$par[3*j-1]*dXSoftPlus(Coef$par[3*j]*Srange+Coef$par[3*j+1])*Coef$par[3*j]
  }
}


###Regression
set.seed(1)
S0<-1
T<-1
K<-1
vol=0.2

TruePrice<-ZeroRateBachelierCall(S0,T,K,vol)$Price

N<-10000
Powers<-0:14
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

Coef = optim(c(0,0,0,0,0,0,0), func1, method = "BFGS")

Srange<-EstPricesimple<-EstDeltasimple<-SimEstDelta<-10:200/100

for (i in 1:length(Srange)){
  EstPricesimple[i]<-t(Coef$par)%*%Srange[i]^Powers
  EstDeltasimple[i]<-t(Coef$par[2:length(Powers)])%*% ( Powers[2:length(Powers)]*(Srange[i]^(Powers[2:length(Powers)]-1)))
}

tst<-ZeroRateBachelierCall(Srange,T,K,vol)
#regression
#plot of option price and m = 1
xrange<-c(min(S1),max(S1))
plot(Srange,tst$Price,type = "l",col = 'black',xlim=xrange,lwd = 2.2, ylab = "Call Price", xlab = "S")
lines(Srange,y,col = 'chartreuse', lwd = 2.2)
lines(Srange, EstPricesimple, col = "blue4", lwd = 2.2)
legend("topleft" , legend=c("True Bachelier", "Softplus", "Poly 9th degree"),
       col=c("black", "chartreuse","blue4" ),lty=1,  cex = 0.75)

#plotting delta 3 degree and m = 1
plot(Srange,tst$Delta,type = "l",col = 'black',xlim=xrange,ylim =c(0,1.5) ,lwd = 2.2,ylab = "Call Price", xlab = "S")
lines(Srange,delta,col = 'chartreuse',lwd = 2.2)
lines(Srange, EstDeltasimple, col = "blue4", lwd = 2.2)
legend("topleft" , legend=c("True Bachelier", "Softplus", "Poly 6th degree"),
       col=c("black", "chartreuse","blue4" ),lty=1,  cex = 0.75)


##Error
##Simulation study
set.seed(1)
S0 = 1
T = 1
K = 1
vol = 0.2

set.seed(1)
N<-10000
Npaths <- N
Sim1<-rnorm(Npaths,0,1)
Sim2<-rnorm(Npaths,0,1)

S1 = S0+vol*sqrt(T)*Sim1
S = S1
S2 = S+vol*sqrt(T)*Sim2
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

TruePrice <- ZeroRateBachelierCall(S,T,K,vol)$Price
TrueDelta <- ZeroRateBachelierCall(S,T,K,vol)$Delta

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

order = 12
Price_error = matrix(0,order,1)
Delta_error = matrix(0,order,1)
for (h in 1:order){
  if (h==1){
    x = c(0,0,0,0)
  } else {
    z = runif(3)
    x <- c(Coef$par, z)
  }
  Srange<-10:200/100
  Coef <- optim(x, func, method = "BFGS")
  print(Coef$par)
  y = Regression_SP(Coef$par,Srange)
  delta = 0
  tst<-ZeroRateBachelierCall(Srange,T,K,vol)
  for (j in 1:((length(Coef$par)-1)/3)){
    delta = delta + Coef$par[3*j-1]*dXSoftPlus(Coef$par[3*j]*Srange+Coef$par[3*j+1])*Coef$par[3*j]
  }
  Price_error[h] = round(mean(abs(y-tst$Price)),8)
  Delta_error[h] = round(mean(abs(delta-tst$Delta)),8)
}

##Reg
set.seed(1)
CallPayoff = CallPrice
NN = 12
Price_errorreg = matrix(0,NN-2,1)
Delta_errorreg = matrix(0,NN-2,1)
for(j in 3:NN){
  if (j == 3){
    x = c(0,0,0,0)
  } else {
    z = runif(1)
    x = c(res2$par, z)
  }
  p = length(x)
  #print(p)
  Powers = 1:p-1
  Powers = c(Powers)
  #print(Powers)
  N = 10000
  
  X = matrix(1,N,1)
  for (i in 2:p){
    X = cbind(X,S1^(Powers[i]))
  }
  X
  ######### NUMERICAL OPTIMIZATION
  x0 = x
  x0 = matrix(0,length(x),1)
  Y = matrix(0,N,1)
  for (i in 2:p){
    Y = cbind(Y,Powers[i]*S1^(Powers[i]-1))
  }
  
  tau = sd(CallPrice)/sd(CallDelta)
  w = 1/(1+tau)
  
  loss <- function(beta) {
    b = beta
    return(w*(t(X%*%b-CallPayoff))%*%(X%*%b-CallPayoff)+(1-w)*(t(Y%*%b-CallDelta))%*%(Y%*%b-CallDelta))
  } 
  
  Grad <- function(beta) {
    b = beta
    return(2*w*((t(X)%*%X)%*%b)-2*w*(t(X)%*%CallPayoff)+2*(1-w)*((t(Y)%*%Y)%*%b)-2*w*(t(Y)%*%CallDelta))
  }
  
  func2 <- function(beta){
    out <- loss(beta)
    attr(out, 'gradient') <- Grad(beta)
    return(out)
  }
  
  res2 = optim(x, func2, method = "BFGS")
  Srange<-EstPricereg<-EstDeltareg<-10:200/100
  tst<-ZeroRateBachelierCall(Srange,T,K,vol)
  for (i in 1:length(Srange)){
    EstPricereg[i]<-res2$par%*%(Srange[i]^Powers)
    EstDeltareg[i] = res2$par[2:length(Powers)]%*%(Powers[2:length(Powers)]*Srange[i]^(Powers[2:length(Powers)]-1))
    #print(Powers[2:length(Powers)])
  }
  
  Price_errorreg[j-2] = round(mean(abs(EstPricereg-tst$Price)),8)
  Delta_errorreg[j-2] = round(mean(abs(EstDeltareg-tst$Delta)),8)
  
}

#Price
x = 4:13
k = 1:order
plot(x,Price_errorreg,type = "l", col = "red",ylim=c(0,0.02),ylab = "error",xlab="order of polynomial")
lines(k, Price_error, col = "blue")
legend("topright" , legend=c("Polynomial regression", "Softplus regression"),
       col=c("red","blue" ),lty=1,  cex = 0.75)
#Delta
x = Powers[4:length(Powers)]
k = Powers[2:length(Powers)]
plot(x,Delta_errorreg,type = "l", col = "red",ylim=c(0,0.1),ylab = "error",xlab="order of polynomial")
lines(k, Delta_error, col = "blue")
legend("topright" , legend=c("Polynomial regression", "Softplus regression"),
       col=c("red","blue" ),lty=1,  cex = 0.75)
