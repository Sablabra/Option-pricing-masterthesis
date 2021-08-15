##Softplus vs poly without reg
source("functions.R")
set.seed(1)
T = 1
K = 1
vol = 0.2

c = 0.01
S = seq(0.10,1.99,c)

CallPrice = TruePrice = ZeroRateBachelierCall(S,1,1,0.2)$Price
CallDelta = TrueDelta = ZeroRateBachelierCall(S,1,1,0.2)$Delta

dummy = sd(CallPrice)/sd(CallDelta)
w = 1/(1+dummy)

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

func <- function(a){
  out <- Criterion(a)
  attr(out, 'gradient') <- Gradient(a)
  return(out)
}

order = 15
Price_error = matrix(0,order,1)
Delta_error = matrix(0,order,1)
for (h in 1:order){
  if (h==1){
    x = c(0,0,0,0)
  } else {
    z = runif(3)
    x <- c(Coef$par, z)
  }
  
  Coef <- optim(x, func, method = "BFGS")
  
  y = Regression_SP(Coef$par,S)
  delta = 0
  for (j in 1:((length(Coef$par)-1)/3)){
    delta = delta + Coef$par[3*j-1]*dXSoftPlus(Coef$par[3*j]*S+Coef$par[3*j+1])*Coef$par[3*j]
  }
  
  Price_error[h] = round(mean(abs(y-TruePrice)),8)
  Delta_error[h] = round(mean(abs(delta-TrueDelta)),8)
}

##regression error with same setup
T = 1
K = 1
vol = 0.2

c = 0.01
S = seq(0.10,1.99,c)

CallPayoff = TruePrice = ZeroRateBachelierCall(S,1,1,0.2)$Price
CallDelta = TrueDelta = ZeroRateBachelierCall(S,1,1,0.2)$Delta

PolyDeg = 15
Price_errorreg = matrix(0,PolyDeg-2,1)
Delta_errorreg = matrix(0,PolyDeg-2,1)
for(j in 1:PolyDeg){
  if (j == 3){
    x = c(0,0,0,0)
  } else {
    z = runif(1)
    x = c(Coef1$par, z)
  }
  p = length(x)
  Powers = 1:p-1
  Powers = c(Powers)
  N = length(S)
  
  X = matrix(1,N,1)
  for (i in 2:p){
    X = cbind(X,S^(Powers[i]))
  }
  x0 = x
  x0 = matrix(0,length(x),1)
  Y = matrix(0,N,1)
  for (i in 2:p){
    Y = cbind(Y,Powers[i]*S^(Powers[i]-1))
  }
  
  dummy = sd(CallPrice)/sd(CallDelta)
  w = 1/(1+dummy)
  
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
  
  Coef1 = optim(x, func2, method = "BFGS")
  
  EstPricereg = matrix(0,length(S),1)
  EstDeltareg = matrix(0,length(S),1)
  for (i in 1:length(S)){
    EstPricereg[i] =Coef1$par%*%(S[i]^Powers)
    EstDeltareg[i] = Coef1$par[2:length(Powers)]%*%(Powers[2:length(Powers)]*S[i]^(Powers[2:length(Powers)]-1))
    #print(Powers[2:length(Powers)])
  }
  
  Price_errorreg[j-2] = round(mean(abs(EstPricereg-TruePrice)),8)
  Delta_errorreg[j-2] = round(mean(abs(EstDeltareg-TrueDelta)),8)
  
}

x = Powers[4:length(Powers)]
k = Powers[2:length(Powers)]
plot(x,Price_errorreg,type = "l", col = "red",ylim=c(0,0.004),ylab = "error",xlab="order of polynomial")
lines(k, Price_error, col = "blue")
legend("topright" , legend=c("Polynomial regression", "Softplus regression"),
       col=c("red","blue" ),lty=1,  cex = 0.75)

#Delta
x = Powers[4:length(Powers)]
k = Powers[2:length(Powers)]
plot(x,Delta_errorreg,type = "l", col = "red",ylim=c(0,0.04),ylab = "error",xlab="order of polynomial")
lines(k, Delta_error, col = "blue")
legend("topright" , legend=c("Polynomial regression", "Softplus regression"),
       col=c("red","blue" ),lty=1,  cex = 0.75)
