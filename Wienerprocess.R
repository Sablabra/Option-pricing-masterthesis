##Plotting of the trajectory of a Wiener Process
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
  w[i,j] <- w[i-1,j] + dw
}
}
t <- seq(0, T, dt)
#jpeg("Fig56.jpg", units="in", width=5, height=5, res=300)
matplot(t, w[,1], type="l", lty=1, col=c("blue4"), ylab = "W(t)")
dev.off()
