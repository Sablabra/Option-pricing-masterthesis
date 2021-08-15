ZeroRateBachelierCall<-function(S,T,K,vol){
  d<-(S-K)/(vol*sqrt(T))
  CallDelta<-pnorm(d,0,1)
  CallPrice<-(S-K)*CallDelta+vol*sqrt(T)*dnorm(d,0,1)
  return(list(Price=CallPrice, Delta=CallDelta))
}


###Softplus function
## first order
SoftPlus <- function(x){
  a = log(1+exp(x))
  return(a)
}

dXSoftPlus <- function(x){
  a = exp(x)/(1+exp(x))
  return(a)
}


dXSoftPlus2 <- function(x){
  a = exp(x)/((1+exp(x))^2)
  return(a)
}

Regression_SP <- function(beta,X){
  y = beta[1]
  for (i in 1:((length(beta)-1)/3)){
    y = y+beta[3*i-1]*SoftPlus(beta[3*i]*X+beta[3*i+1])
  }
  return(y)
}



## order 2
Regression_SP2 <- function(beta,X){
  y = beta[1]+beta[2]*SoftPlus(beta[3]*X+beta[4])+beta[5]*SoftPlus(beta[6]*X+beta[7])
  return(y)
}
## order 3
Regression_SP3 <- function(beta,X){
  y = beta[1]+beta[2]*SoftPlus(beta[3]*X+beta[4])+beta[5]*SoftPlus(beta[6]*X+beta[7])+beta[8]*SoftPlus(beta[9]*X+beta[10])
  return(y)
}

## order 4
Regression_SP4 <- function(beta,X){
  y = beta[1]+beta[2]*SoftPlus(beta[3]*X+beta[4])+beta[5]*SoftPlus(beta[6]*X+beta[7])+beta[8]*SoftPlus(beta[9]*X+beta[10])+beta[11]*SoftPlus(beta[12]*X+beta[13])
  return(y)
}  
## order 5
Regression_SP5 <- function(beta,X){
  y = beta[1]+beta[2]*SoftPlus(beta[3]*X+beta[4])+beta[5]*SoftPlus(beta[6]*X+beta[7])+beta[8]*SoftPlus(beta[9]*X+beta[10])+beta[11]*SoftPlus(beta[12]*X+beta[13])+beta[14]*SoftPlus(beta[15]*X+beta[16])
  return(y)
}  

## order 6
Regression_SP6 <-function(beta,X){
  y = beta[1]+beta[2]*SoftPlus(beta[3]*X+beta[4])+beta[5]*SoftPlus(beta[6]*X+beta[7])+beta[8]*SoftPlus(beta[9]*X+beta[10])+beta[11]*SoftPlus(beta[12]*X+beta[13])+beta[14]*SoftPlus(beta[15]*X+beta[16])+beta[17]*SoftPlus(beta[18]*X+beta[19])
  return(y)
}

##order 7
Regression_SP7 <-function(beta,X){
  y = beta[1]+beta[2]*SoftPlus(beta[3]*X+beta[4])+beta[5]*SoftPlus(beta[6]*X+beta[7])+beta[8]*SoftPlus(beta[9]*X+beta[10])+beta[11]*SoftPlus(beta[12]*X+beta[13])+beta[14]*SoftPlus(beta[15]*X+beta[16])+beta[17]*SoftPlus(beta[18]*X+beta[19])+beta[20]*SoftPlus(beta[21]*X+beta[22])
  return(y)
}

