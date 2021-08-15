##SoftPlus with no regulization 3 degree
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

order = 1
for (h in 1:order){
  if (h==1){
    x = c(0,0,0,0)
  } else {
    z = runif(3)
    x <- c(CoefSP$par, z)
  }

  CoefSP <- optim(x, func, method = "BFGS")

  if (length(x) == 4){
    y = Regression_SP(CoefSP$par,S)
    delta = CoefSP$par[2]*dXSoftPlus(CoefSP$par[3]*S+CoefSP$par[4])*CoefSP$par[3]
  } else if (length(x)==7){
    y = Regression_SP2(CoefSP$par,S)
    delta = CoefSP$par[2]*dXSoftPlus(CoefSP$par[3]*S+CoefSP$par[4])*CoefSP$par[3]+CoefSP$par[5]*dXSoftPlus(CoefSP$par[6]*S+CoefSP$par[7])*CoefSP$par[6]
  } else if (length(x)==10){
    y = Regression_SP3(CoefSP$par,S)
    delta = CoefSP$par[2]*dXSoftPlus(CoefSP$par[3]*S+CoefSP$par[4])*CoefSP$par[3]+CoefSP$par[5]*dXSoftPlus(CoefSP$par[6]*S+CoefSP$par[7])*CoefSP$par[6]+CoefSP$par[8]*dXSoftPlus(CoefSP$par[9]*S+CoefSP$par[10])*CoefSP$par[9]
  } else if (length(x)==13){
    y = Regression_SP4(CoefSP$par,S)
    delta = CoefSP$par[2]*dXSoftPlus(CoefSP$par[3]*S+CoefSP$par[4])*CoefSP$par[3]+CoefSP$par[5]*dXSoftPlus(CoefSP$par[6]*S+CoefSP$par[7])*CoefSP$par[6]+CoefSP$par[8]*dXSoftPlus(CoefSP$par[9]*S+CoefSP$par[10])*CoefSP$par[9]+CoefSP$par[11]*dXSoftPlus(CoefSP$par[12]*S+CoefSP$par[13])*CoefSP$par[12]
  } else if (length(x)==16){
    y = Regression_SP5(CoefSP$par,S)
    delta = CoefSP$par[2]*dXSoftPlus(CoefSP$par[3]*S+CoefSP$par[4])*CoefSP$par[3]+CoefSP$par[5]*dXSoftPlus(CoefSP$par[6]*S+CoefSP$par[7])*CoefSP$par[6]+CoefSP$par[8]*dXSoftPlus(CoefSP$par[9]*S+CoefSP$par[10])*CoefSP$par[9]+CoefSP$par[11]*dXSoftPlus(CoefSP$par[12]*S+CoefSP$par[13])*CoefSP$par[12]+CoefSP$par[14]*dXSoftPlus(CoefSP$par[15]*S+CoefSP$par[16])*CoefSP$par[15]
  }
}

#simple regression without regularization 3,6,9 degree
set.seed(3)
Sim1 = rnorm(length(S),0,1)
Sim2 = rnorm(length(S),0,1)
S = seq(.10, 1.99, c)

CallPrice = TruePrice = ZeroRateBachelierCall(S,1,1,0.2)$Price
CallDelta = TrueDelta = ZeroRateBachelierCall(S,1,1,0.2)$Delta

S2 = S+vol*sqrt(T)*Sim2

CallPayoff = CallPrice

Powers<-0:3
M<-length(Powers)
N <- length(S)

X1<-S^0; 
X2<-rep(0, length(S))


dummy<-sd(CallPayoff)/sd(CallDelta)
w<-1/(1+dummy)

for (k in 2:length(Powers)) {
  X1<-cbind(X1,S^Powers[k])
  X2<-cbind(X2,Powers[k]*S^(Powers[k]-1))
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

Coef = optim(c(0,0,0,0), func1, method = "BFGS")

EstPrice = matrix(0,length(S),1)
EstDelta = matrix(0,length(S),1)
for (i in 1:length(S)){
  EstPrice[i] =Coef$par%*%(S[i]^Powers)
  EstDelta[i] = Coef$par[2:length(Powers)]%*%(Powers[2:length(Powers)]*S[i]^(Powers[2:length(Powers)]-1))
}

#regression
#Bachelier payoff .. degree and m = ..
plot(S,TruePrice,type = "l",col = 'black',lwd = 2.2, ylab = "Call Price", xlab = "S")
lines(S,y,col = 'chartreuse', lwd = 2.2)
lines(S, EstPrice, col = "blue4", lwd = 2.2)
legend("topleft" , legend=c("True Bachelier", "Softplus", "Poly 9th degree"),
       col=c("black", "chartreuse","blue4" ),lty=1,  cex = 0.75)

#plotting delta .. degree and m = ..
plot(S,TrueDelta,type = "l",col = 'black',lwd = 2.2, ylim=c(-0.5,2.5),ylab = "Call Price", xlab = "S")
lines(S,delta,col = 'chartreuse',lwd = 2.2)
lines(S, EstDelta, col = "blue4", lwd = 2.2)
legend("topleft" , legend=c("True Bachelier", "Softplus", "Poly 9th degree"),
       col=c("black", "chartreuse","blue4" ),lty=1,  cex = 0.75)

