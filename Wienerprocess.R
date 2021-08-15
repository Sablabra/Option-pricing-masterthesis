set.seed(1)
x0 <- 0
n <- 1
T <- 1
m <- 1000
dt <- T/m
w <- matrix(0,m+1,n)
w[1,] <- x0
for(j in 1:n)
{for(i in 2:(m+1)){
  dw <- sqrt(dt)*rnorm(1,0,1)
  w[i,j] <- x[i-1,j] + dx
}
}
t <- seq(0, T, dt)
matplot(t, x[,1], type="l", lty=1, col=c("blue4"), ylab = "W(t)")
