##Hedge ekspirement ML
simula <- c(1000,2000,3000,4000,5000,6000)

ll1 = c(0.23191002,0.17117363,0.15281885,0.1511261,0.147391,0.14718048)
mm1 = c(0.17151353,0.13547317,0.12981786,0.12737579, 0.12405673,0.12624865)
ll2 = c(0.16345892,0.14368641,0.1310048,0.12441858,0.12988172,0.124568)
mm2 = c(0.12715786,0.12187419, 0.12109752,0.12036936,0.1217376,0.11964939)
ll3 = c(0.13956431,0.13345105,0.1255626,0.12427101,0.12517285,0.122297)
mm3 = c(0.12267188,0.12138273,0.12121968,0.12039326,0.12196928,0.11949841)
ll4 = c(0.14849349,0.12757571,0.12887116,0.12775931,0.13313605,0.12263616)
mm4 = c(0.12334028,0.12067686,0.12410815,0.12014395,0.12050922,0.11954408)


plot(simula,ll1,ylim=c(0.09,0.23),type = "l",col="red",xlab=c("Number of simulations"),ylab=c("SDE of relative hedge error in SoftPlus"))
lines(simula,mm1,col="blue",type = "b")
lines(simula,ll2,col="red",type = "l")
lines(simula,mm2,col="blue",type = "b")
lines(simula,ll3,col="red",type="l")
lines(simula,mm3,col="blue",type = "b")
lines(simula,ll4,col="red",type="l")
lines(simula,mm4,col="blue",type = "b")
abline(h=0.1208,col="orange")
legend("topright" , legend=c("Standard Neural Network", "Twin Neural Network","True Delta"),
       col=c("red", "blue","orange"),lty=1,  cex = 0.65)
text(5975,0.152,"L=1",cex=0.6)
text(5950,0.135,"L=2,3,4",cex=0.6)
text(1000,0.18,"L=1",cex=0.6)
text(1050,0.135,"L=2,3,4",cex=0.6)
