base.dir=file.path('C:/Users/Tharshan-PC/Desktop/Research/QuasiSeq')
out.dir=file.path(base.dir, 'out')
prg.dir=file.path(base.dir, 'prg')
dat.dir=file.path(base.dir, 'dat')
img.dir=file.path(base.dir, 'img')

if(F){
	library(parallel)
	RNGkind("L'Ecuyer-CMRG")
	cl=makeCluster(8L)
	clusterEvalQ(cl, {
		RNGkind("L'Ecuyer-CMRG"); 
		options(contrasts=c('contr.sum', 'contr.sum'))
	})
	clusterSetRNGStream(cl , 9992722L)
}
set.seed(9992722L)

prg.name=paste('Coner', Sys.Date(), sep='_')
img.name=file.path(img.dir, paste(prg.name, 'RData', sep='.'))
out.file=function(suffix, subdir) if(missing(subdir)) file.path(out.dir, paste(prg.name, suffix, sep='_')) else file.path(out.dir, subdir, paste(prg.name, suffix, sep='_')) 

setwd(out.dir)

options(contrasts=c('contr.sum', 'contr.sum'))

.sessionInfo=list(sessionInfo=sessionInfo(), Sys.info=Sys.info(),
  Sys.getenv=Sys.getenv(names=TRUE), capabilities=capabilities(),
  options=options(), RNGkind=RNGkind(), .libPaths=.libPaths(), 
  date=date(), wd=getwd(),.Random.seed=.Random.seed)


####################################### END OF HEADER ################################################
<<<<<<< HEAD
=======
		
	
QNB=function(mu, y, tau, phi)
{
	
}

BNB=function(mu, y, tau, phi)
{
	
}

#BNB and QLNB changing with different pi values.
=====================

z<-1:100
pi<-c(2,4,6,1000000)
t<-c(0.1,0.3,0.5,0.8)#tho values
#x_mean values
curve(QNB(mu=x, y=10, tau=.5, phi=1.1), 0, 50)

curve((((gamma(z[7]+x*(t[1]*pi[1]+1)/(pi[1]-1)))*(beta((x*(t[1]*pi[1]+1))/(pi[1]-1)+(t[1]*(2*pi[1]-1)+1)/(t[1]*(pi[1]-1)),z[7]+1/t[1])))/
((factorial(z[7])*gamma(x*(t[1]*pi[1]+1)/(pi[1]-1)))*(beta((t[1]*(2*pi[1]-1)+1)/(t[1]*(pi[1]-1)),1/t[1])
)))/(exp(1/pi[1])*((1/(1+t[1]*x))^(1/t[1])*(x/(x+1/t[1]))^z[7])),0.01,0.02,100,ylab="c",col="red")
curve((((gamma(z[7]+x*(t[1]*pi[2]+1)/(pi[2]-1)))*(beta((x*(t[1]*pi[2]+1))/(pi[2]-1)+(t[1]*(2*pi[2]-1)+1)/(t[1]*(pi[2]-1)),z[7]+1/t[1])))/
((factorial(z[7])*gamma(x*(t[1]*pi[2]+1)/(pi[2]-1)))*(beta((t[1]*(2*pi[2]-1)+1)/(t[1]*(pi[2]-1)),1/t[1])
)))/(exp(1/pi[2])*((1/(1+t[1]*x))^(1/t[1])*(x/(x+1/t[1]))^z[7])),0.01,0.02,100,ylab="c",col="blue",add=T)
curve((((gamma(z[7]+x*(t[1]*pi[3]+1)/(pi[3]-1)))*(beta((x*(t[1]*pi[3]+1))/(pi[3]-1)+(t[1]*(2*pi[3]-1)+1)/(t[1]*(pi[3]-1)),z[7]+1/t[1])))/
((factorial(z[7])*gamma(x*(t[1]*pi[3]+1)/(pi[3]-1)))*(beta((t[1]*(2*pi[3]-1)+1)/(t[1]*(pi[3]-1)),1/t[1])
)))/(exp(1/pi[3])*((1/(1+t[1]*x))^(1/t[1])*(x/(x+1/t[1]))^z[7])),0.01,0.02,100,ylab="c",col="black",add=T)
curve((((gamma(z[7]+x*(t[1]*pi[4]+1)/(pi[4]-1)))*(beta((x*(t[1]*pi[4]+1))/(pi[4]-1)+(t[1]*(2*pi[4]-1)+1)/(t[1]*(pi[4]-1)),z[7]+1/t[1])))/
((factorial(z[7])*gamma(x*(t[1]*pi[4]+1)/(pi[4]-1)))*(beta((t[1]*(2*pi[4]-1)+1)/(t[1]*(pi[4]-1)),1/t[1])
)))/(exp(1/pi[4])*((1/(1+t[1]*x))^(1/t[1])*(x/(x+1/t[1]))^z[7])),0.01,0.02,100,ylab="c",col="gray",add=T)
legend(0.016,5*10^15, c("pi=2","pi=4","pi=6","pi=1000000"), cex=0.8,col=c("red","blue","black","gray"),pch=20:23,lwd=1:3)

===========================
curve((((gamma(z[7]+x*(t[1]*pi[1]+1)/(pi[1]-1)))*(beta((x*(t[1]*pi[1]+1))/(pi[1]-1)+(t[1]*(2*pi[1]-1)+1)/(t[1]*(pi[1]-1)),z[7]+1/t[1])))/
((factorial(z[7])*gamma(x*(t[1]*pi[1]+1)/(pi[1]-1)))*(beta((t[1]*(2*pi[1]-1)+1)/(t[1]*(pi[1]-1)),1/t[1])
)))/(exp(1/pi[1])*((1/(1+t[1]*x))^(1/t[1])*(x/(x+1/t[1]))^z[7])),0.02,0.05,100,ylab="c",col="red")
curve((((gamma(z[7]+x*(t[1]*pi[2]+1)/(pi[2]-1)))*(beta((x*(t[1]*pi[2]+1))/(pi[2]-1)+(t[1]*(2*pi[2]-1)+1)/(t[1]*(pi[2]-1)),z[7]+1/t[1])))/
((factorial(z[7])*gamma(x*(t[1]*pi[2]+1)/(pi[2]-1)))*(beta((t[1]*(2*pi[2]-1)+1)/(t[1]*(pi[2]-1)),1/t[1])
)))/(exp(1/pi[2])*((1/(1+t[1]*x))^(1/t[1])*(x/(x+1/t[1]))^z[7])),0.02,0.05,100,ylab="c",col="blue",add=T)
curve((((gamma(z[7]+x*(t[1]*pi[3]+1)/(pi[3]-1)))*(beta((x*(t[1]*pi[3]+1))/(pi[3]-1)+(t[1]*(2*pi[3]-1)+1)/(t[1]*(pi[3]-1)),z[7]+1/t[1])))/
((factorial(z[7])*gamma(x*(t[1]*pi[3]+1)/(pi[3]-1)))*(beta((t[1]*(2*pi[3]-1)+1)/(t[1]*(pi[3]-1)),1/t[1])
)))/(exp(1/pi[3])*((1/(1+t[1]*x))^(1/t[1])*(x/(x+1/t[1]))^z[7])),0.02,0.05,100,ylab="c",col="black",add=T)
curve((((gamma(z[7]+x*(t[1]*pi[4]+1)/(pi[4]-1)))*(beta((x*(t[1]*pi[4]+1))/(pi[4]-1)+(t[1]*(2*pi[4]-1)+1)/(t[1]*(pi[4]-1)),z[7]+1/t[1])))/
((factorial(z[7])*gamma(x*(t[1]*pi[4]+1)/(pi[4]-1)))*(beta((t[1]*(2*pi[4]-1)+1)/(t[1]*(pi[4]-1)),1/t[1])
)))/(exp(1/pi[4])*((1/(1+t[1]*x))^(1/t[1])*(x/(x+1/t[1]))^z[7])),0.02,0.05,100,ylab="c",col="gray",add=T)
legend(0.040,8*10^13, c("pi=2","pi=4","pi=6","pi=1000000"), cex=0.8,col=c("red","blue","black","gray"),pch=20:23,lwd=1:3)

=====================================

curve((((gamma(z[7]+x*(t[1]*pi[1]+1)/(pi[1]-1)))*(beta((x*(t[1]*pi[1]+1))/(pi[1]-1)+(t[1]*(2*pi[1]-1)+1)/(t[1]*(pi[1]-1)),z[7]+1/t[1])))/
((factorial(z[7])*gamma(x*(t[1]*pi[1]+1)/(pi[1]-1)))*(beta((t[1]*(2*pi[1]-1)+1)/(t[1]*(pi[1]-1)),1/t[1])
)))/(exp(1/pi[1])*((1/(1+t[1]*x))^(1/t[1])*(x/(x+1/t[1]))^z[7])),0.1,0.2,100,ylab="y",col="red")
curve((((gamma(z[7]+x*(t[1]*pi[2]+1)/(pi[2]-1)))*(beta((x*(t[1]*pi[2]+1))/(pi[2]-1)+(t[1]*(2*pi[2]-1)+1)/(t[1]*(pi[2]-1)),z[7]+1/t[1])))/
((factorial(z[7])*gamma(x*(t[1]*pi[2]+1)/(pi[2]-1)))*(beta((t[1]*(2*pi[2]-1)+1)/(t[1]*(pi[2]-1)),1/t[1])
)))/(exp(1/pi[2])*((1/(1+t[1]*x))^(1/t[1])*(x/(x+1/t[1]))^z[7])),0.1,0.2,100,ylab="y",col="blue",add=T)
curve((((gamma(z[7]+x*(t[1]*pi[3]+1)/(pi[3]-1)))*(beta((x*(t[1]*pi[3]+1))/(pi[3]-1)+(t[1]*(2*pi[3]-1)+1)/(t[1]*(pi[3]-1)),z[7]+1/t[1])))/
((factorial(z[7])*gamma(x*(t[1]*pi[3]+1)/(pi[3]-1)))*(beta((t[1]*(2*pi[3]-1)+1)/(t[1]*(pi[3]-1)),1/t[1])
)))/(exp(1/pi[3])*((1/(1+t[1]*x))^(1/t[1])*(x/(x+1/t[1]))^z[7])),0.1,0.2,100,ylab="y",col="black",add=T)
curve((((gamma(z[7]+x*(t[1]*pi[4]+1)/(pi[4]-1)))*(beta((x*(t[1]*pi[4]+1))/(pi[4]-1)+(t[1]*(2*pi[4]-1)+1)/(t[1]*(pi[4]-1)),z[7]+1/t[1])))/
((factorial(z[7])*gamma(x*(t[1]*pi[4]+1)/(pi[4]-1)))*(beta((t[1]*(2*pi[4]-1)+1)/(t[1]*(pi[4]-1)),1/t[1])
)))/(exp(1/pi[4])*((1/(1+t[1]*x))^(1/t[1])*(x/(x+1/t[1]))^z[7])),0.1,0.2,100,ylab="y",col="gray",add=T)
legend(0.170,6*10^9, c("pi=2","pi=4","pi=6","pi=1000000"), cex=0.8,col=c("red","blue","black","gray"),pch=20:23,lwd=1:3)

====================================================

curve((((gamma(z[7]+x*(t[1]*pi[1]+1)/(pi[1]-1)))*(beta((x*(t[1]*pi[1]+1))/(pi[1]-1)+(t[1]*(2*pi[1]-1)+1)/(t[1]*(pi[1]-1)),z[7]+1/t[1])))/
((factorial(z[7])*gamma(x*(t[1]*pi[1]+1)/(pi[1]-1)))*(beta((t[1]*(2*pi[1]-1)+1)/(t[1]*(pi[1]-1)),1/t[1])
)))/(exp(1/pi[1])*((1/(1+t[1]*x))^(1/t[1])*(x/(x+1/t[1]))^z[7])),0.2,0.5,100,ylab="y",col="red")
curve((((gamma(z[7]+x*(t[1]*pi[2]+1)/(pi[2]-1)))*(beta((x*(t[1]*pi[2]+1))/(pi[2]-1)+(t[1]*(2*pi[2]-1)+1)/(t[1]*(pi[2]-1)),z[7]+1/t[1])))/
((factorial(z[7])*gamma(x*(t[1]*pi[2]+1)/(pi[2]-1)))*(beta((t[1]*(2*pi[2]-1)+1)/(t[1]*(pi[2]-1)),1/t[1])
)))/(exp(1/pi[2])*((1/(1+t[1]*x))^(1/t[1])*(x/(x+1/t[1]))^z[7])),0.2,0.5,100,ylab="y",col="blue",add=T)
curve((((gamma(z[7]+x*(t[1]*pi[3]+1)/(pi[3]-1)))*(beta((x*(t[1]*pi[3]+1))/(pi[3]-1)+(t[1]*(2*pi[3]-1)+1)/(t[1]*(pi[3]-1)),z[7]+1/t[1])))/
((factorial(z[7])*gamma(x*(t[1]*pi[3]+1)/(pi[3]-1)))*(beta((t[1]*(2*pi[3]-1)+1)/(t[1]*(pi[3]-1)),1/t[1])
)))/(exp(1/pi[3])*((1/(1+t[1]*x))^(1/t[1])*(x/(x+1/t[1]))^z[7])),0.2,0.5,100,ylab="y",col="black",add=T)
curve((((gamma(z[7]+x*(t[1]*pi[4]+1)/(pi[4]-1)))*(beta((x*(t[1]*pi[4]+1))/(pi[4]-1)+(t[1]*(2*pi[4]-1)+1)/(t[1]*(pi[4]-1)),z[7]+1/t[1])))/
((factorial(z[7])*gamma(x*(t[1]*pi[4]+1)/(pi[4]-1)))*(beta((t[1]*(2*pi[4]-1)+1)/(t[1]*(pi[4]-1)),1/t[1])
)))/(exp(1/pi[4])*((1/(1+t[1]*x))^(1/t[1])*(x/(x+1/t[1]))^z[7])),0.2,0.5,100,ylab="y",col="gray",add=T)
legend(0.40,10^8, c("pi=2","pi=4","pi=6","pi=1000000"), cex=0.8,col=c("red","blue","black","gray"),pch=20:23,lwd=1:3)

==============================================
curve((((gamma(z[7]+x*(t[1]*pi[1]+1)/(pi[1]-1)))*(beta((x*(t[1]*pi[1]+1))/(pi[1]-1)+(t[1]*(2*pi[1]-1)+1)/(t[1]*(pi[1]-1)),z[7]+1/t[1])))/
((factorial(z[7])*gamma(x*(t[1]*pi[1]+1)/(pi[1]-1)))*(beta((t[1]*(2*pi[1]-1)+1)/(t[1]*(pi[1]-1)),1/t[1])
)))/(exp(1/pi[1])*((1/(1+t[1]*x))^(1/t[1])*(x/(x+1/t[1]))^z[7])),0.5,1,100,ylab="y",col="red")
curve((((gamma(z[7]+x*(t[1]*pi[2]+1)/(pi[2]-1)))*(beta((x*(t[1]*pi[2]+1))/(pi[2]-1)+(t[1]*(2*pi[2]-1)+1)/(t[1]*(pi[2]-1)),z[7]+1/t[1])))/
((factorial(z[7])*gamma(x*(t[1]*pi[2]+1)/(pi[2]-1)))*(beta((t[1]*(2*pi[2]-1)+1)/(t[1]*(pi[2]-1)),1/t[1])
)))/(exp(1/pi[2])*((1/(1+t[1]*x))^(1/t[1])*(x/(x+1/t[1]))^z[7])),0.5,1,100,ylab="y",col="blue",add=T)
curve((((gamma(z[7]+x*(t[1]*pi[3]+1)/(pi[3]-1)))*(beta((x*(t[1]*pi[3]+1))/(pi[3]-1)+(t[1]*(2*pi[3]-1)+1)/(t[1]*(pi[3]-1)),z[7]+1/t[1])))/
((factorial(z[7])*gamma(x*(t[1]*pi[3]+1)/(pi[3]-1)))*(beta((t[1]*(2*pi[3]-1)+1)/(t[1]*(pi[3]-1)),1/t[1])
)))/(exp(1/pi[3])*((1/(1+t[1]*x))^(1/t[1])*(x/(x+1/t[1]))^z[7])),0.5,1,100,ylab="y",col="black",add=T)
curve((((gamma(z[7]+x*(t[1]*pi[4]+1)/(pi[4]-1)))*(beta((x*(t[1]*pi[4]+1))/(pi[4]-1)+(t[1]*(2*pi[4]-1)+1)/(t[1]*(pi[4]-1)),z[7]+1/t[1])))/
((factorial(z[7])*gamma(x*(t[1]*pi[4]+1)/(pi[4]-1)))*(beta((t[1]*(2*pi[4]-1)+1)/(t[1]*(pi[4]-1)),1/t[1])
)))/(exp(1/pi[4])*((1/(1+t[1]*x))^(1/t[1])*(x/(x+1/t[1]))^z[7])),0.5,1,100,ylab="y",col="gray",add=T)
legend(0.80,2*10^6, c("pi=2","pi=4","pi=6","pi=1000000"), cex=0.8,col=c("red","blue","black","gray"),pch=20:23,lwd=1:3)
============================================

curve((((gamma(z[7]+x*(t[1]*pi[1]+1)/(pi[1]-1)))*(beta((x*(t[1]*pi[1]+1))/(pi[1]-1)+(t[1]*(2*pi[1]-1)+1)/(t[1]*(pi[1]-1)),z[7]+1/t[1])))/
((factorial(z[7])*gamma(x*(t[1]*pi[1]+1)/(pi[1]-1)))*(beta((t[1]*(2*pi[1]-1)+1)/(t[1]*(pi[1]-1)),1/t[1])
)))/(exp(1/pi[1])*((1/(1+t[1]*x))^(1/t[1])*(x/(x+1/t[1]))^z[7])),1,4,1000,ylab="y",col="red")
curve((((gamma(z[7]+x*(t[1]*pi[2]+1)/(pi[2]-1)))*(beta((x*(t[1]*pi[2]+1))/(pi[2]-1)+(t[1]*(2*pi[2]-1)+1)/(t[1]*(pi[2]-1)),z[7]+1/t[1])))/
((factorial(z[7])*gamma(x*(t[1]*pi[2]+1)/(pi[2]-1)))*(beta((t[1]*(2*pi[2]-1)+1)/(t[1]*(pi[2]-1)),1/t[1])
)))/(exp(1/pi[2])*((1/(1+t[1]*x))^(1/t[1])*(x/(x+1/t[1]))^z[7])),1,4,1000,ylab="y",col="blue",add=T)
curve((((gamma(z[7]+x*(t[1]*pi[3]+1)/(pi[3]-1)))*(beta((x*(t[1]*pi[3]+1))/(pi[3]-1)+(t[1]*(2*pi[3]-1)+1)/(t[1]*(pi[3]-1)),z[7]+1/t[1])))/
((factorial(z[7])*gamma(x*(t[1]*pi[3]+1)/(pi[3]-1)))*(beta((t[1]*(2*pi[3]-1)+1)/(t[1]*(pi[3]-1)),1/t[1])
)))/(exp(1/pi[3])*((1/(1+t[1]*x))^(1/t[1])*(x/(x+1/t[1]))^z[7])),1,4,1000,ylab="y",col="black",add=T)
curve((((gamma(z[7]+x*(t[1]*pi[4]+1)/(pi[4]-1)))*(beta((x*(t[1]*pi[4]+1))/(pi[4]-1)+(t[1]*(2*pi[4]-1)+1)/(t[1]*(pi[4]-1)),z[7]+1/t[1])))/
((factorial(z[7])*gamma(x*(t[1]*pi[4]+1)/(pi[4]-1)))*(beta((t[1]*(2*pi[4]-1)+1)/(t[1]*(pi[4]-1)),1/t[1])
)))/(exp(1/pi[4])*((1/(1+t[1]*x))^(1/t[1])*(x/(x+1/t[1]))^z[7])),1,4,1000,ylab="y",col="gray",add=T)
legend(3,60000, c("pi=2","pi=4","pi=6","pi=1000000"), cex=0.8,col=c("red","blue","black","gray"),pch=20:23,lwd=1:3)

=========================

curve((((gamma(z[7]+x*(t[1]*pi[1]+1)/(pi[1]-1)))*(beta((x*(t[1]*pi[1]+1))/(pi[1]-1)+(t[1]*(2*pi[1]-1)+1)/(t[1]*(pi[1]-1)),z[7]+1/t[1])))/
((factorial(z[7])*gamma(x*(t[1]*pi[1]+1)/(pi[1]-1)))*(beta((t[1]*(2*pi[1]-1)+1)/(t[1]*(pi[1]-1)),1/t[1])
)))/(exp(1/pi[1])*((1/(1+t[1]*x))^(1/t[1])*(x/(x+1/t[1]))^z[7])),3,5,1000,ylab="y",col="red")
curve((((gamma(z[7]+x*(t[1]*pi[2]+1)/(pi[2]-1)))*(beta((x*(t[1]*pi[2]+1))/(pi[2]-1)+(t[1]*(2*pi[2]-1)+1)/(t[1]*(pi[2]-1)),z[7]+1/t[1])))/
((factorial(z[7])*gamma(x*(t[1]*pi[2]+1)/(pi[2]-1)))*(beta((t[1]*(2*pi[2]-1)+1)/(t[1]*(pi[2]-1)),1/t[1])
)))/(exp(1/pi[2])*((1/(1+t[1]*x))^(1/t[1])*(x/(x+1/t[1]))^z[7])),3,5,1000,ylab="y",col="blue",add=T)
curve((((gamma(z[7]+x*(t[1]*pi[3]+1)/(pi[3]-1)))*(beta((x*(t[1]*pi[3]+1))/(pi[3]-1)+(t[1]*(2*pi[3]-1)+1)/(t[1]*(pi[3]-1)),z[7]+1/t[1])))/
((factorial(z[7])*gamma(x*(t[1]*pi[3]+1)/(pi[3]-1)))*(beta((t[1]*(2*pi[3]-1)+1)/(t[1]*(pi[3]-1)),1/t[1])
)))/(exp(1/pi[3])*((1/(1+t[1]*x))^(1/t[1])*(x/(x+1/t[1]))^z[7])),3,5,1000,ylab="y",col="black",add=T)
curve((((gamma(z[7]+x*(t[1]*pi[4]+1)/(pi[4]-1)))*(beta((x*(t[1]*pi[4]+1))/(pi[4]-1)+(t[1]*(2*pi[4]-1)+1)/(t[1]*(pi[4]-1)),z[7]+1/t[1])))/
((factorial(z[7])*gamma(x*(t[1]*pi[4]+1)/(pi[4]-1)))*(beta((t[1]*(2*pi[4]-1)+1)/(t[1]*(pi[4]-1)),1/t[1])
)))/(exp(1/pi[4])*((1/(1+t[1]*x))^(1/t[1])*(x/(x+1/t[1]))^z[7])),3,5,1000,ylab="y",col="gray",add=T)
legend(4,7000, c("pi=2","pi=4","pi=6","pi=1000000"), cex=0.8,col=c("red","blue","black","gray"),pch=20:23,lwd=1:3)

===============================
curve((((gamma(z[7]+x*(t[1]*pi[1]+1)/(pi[1]-1)))*(beta((x*(t[1]*pi[1]+1))/(pi[1]-1)+(t[1]*(2*pi[1]-1)+1)/(t[1]*(pi[1]-1)),z[7]+1/t[1])))/
((factorial(z[7])*gamma(x*(t[1]*pi[1]+1)/(pi[1]-1)))*(beta((t[1]*(2*pi[1]-1)+1)/(t[1]*(pi[1]-1)),1/t[1])
)))/(exp(1/pi[1])*((1/(1+t[1]*x))^(1/t[1])*(x/(x+1/t[1]))^z[7])),5,10,1000,ylab="y",col="red")
curve((((gamma(z[7]+x*(t[1]*pi[2]+1)/(pi[2]-1)))*(beta((x*(t[1]*pi[2]+1))/(pi[2]-1)+(t[1]*(2*pi[2]-1)+1)/(t[1]*(pi[2]-1)),z[7]+1/t[1])))/
((factorial(z[7])*gamma(x*(t[1]*pi[2]+1)/(pi[2]-1)))*(beta((t[1]*(2*pi[2]-1)+1)/(t[1]*(pi[2]-1)),1/t[1])
)))/(exp(1/pi[2])*((1/(1+t[1]*x))^(1/t[1])*(x/(x+1/t[1]))^z[7])),5,10,1000,ylab="y",col="blue",add=T)
curve((((gamma(z[7]+x*(t[1]*pi[3]+1)/(pi[3]-1)))*(beta((x*(t[1]*pi[3]+1))/(pi[3]-1)+(t[1]*(2*pi[3]-1)+1)/(t[1]*(pi[3]-1)),z[7]+1/t[1])))/
((factorial(z[7])*gamma(x*(t[1]*pi[3]+1)/(pi[3]-1)))*(beta((t[1]*(2*pi[3]-1)+1)/(t[1]*(pi[3]-1)),1/t[1])
)))/(exp(1/pi[3])*((1/(1+t[1]*x))^(1/t[1])*(x/(x+1/t[1]))^z[7])),5,10,1000,ylab="y",col="black",add=T)
curve((((gamma(z[7]+x*(t[1]*pi[4]+1)/(pi[4]-1)))*(beta((x*(t[1]*pi[4]+1))/(pi[4]-1)+(t[1]*(2*pi[4]-1)+1)/(t[1]*(pi[4]-1)),z[7]+1/t[1])))/
((factorial(z[7])*gamma(x*(t[1]*pi[4]+1)/(pi[4]-1)))*(beta((t[1]*(2*pi[4]-1)+1)/(t[1]*(pi[4]-1)),1/t[1])
)))/(exp(1/pi[4])*((1/(1+t[1]*x))^(1/t[1])*(x/(x+1/t[1]))^z[7])),5,10,1000,ylab="y",col="gray",add=T)
legend(6,6200, c("pi=2","pi=4","pi=6","pi=1000000"), cex=0.8,col=c("red","blue","black","gray"),pch=20:23,lwd=1:3)

============================
>>>>>>> 07ba1cc93d8956b6f4005e90379e318e67cb66c1

#####Graphical comparison QLNB and LBNB

my_QLNB<-function(x,pi,t,y){
QLNB<-exp(1/pi)*((1/(1+t*x))^(1/t)*(t*x/(t*x+1))^y)
return(QLNB)
}

my_LBNB<-function(x,pi,t,y){
LBNB<-
exp(lgamma(y+x*(t*pi+1)/(pi-1))+lbeta((x*(t*pi+1))/(pi-1)+(t*(2*pi-1)+1)/(t*(pi-1)),y+1/t)
-log(factorial(y))-lgamma(x*(t*pi+1)/(pi-1))-lbeta((t*(2*pi-1)+1)/t*(pi-1),1/t))
return(LBNB)
}
		
curve(my_LBNB(x,pi=2,t=0.1,y=10),0.01,100,
main="comparison of BNB and QLNB(pi=2,y=10,t=0.1)",xlab="Mean",ylab="QLNB&BNB",pch=21,lwd=2)
curve(my_QLNB(x,pi=2,t=0.1,y=10),0,01,100,add=T,xlab="Mean",ylab="QLNB&BNB",col="red",
main="comparison of BNB and QLNB(pi=2,y=10,t=0.1)",pch=22,lwd=3)
legend("topright", c("BNB","QLNB"), cex=0.8,col=c("black","red"),pch=21:22,lwd=2:3)

##################END####################################

#This is the function to find exp(h(y,tau,pi))=H,where 
#exp(h(y,tau=t,pi))=integration(f(QLNB)-f(LBNB)) over the mean values
#f(QLNB)=exp(1/pi)*((1/(1+t*x))^(1/t)*(t*x/(t*x+1))^y)
#f(LBNB)=(exp(lgamma(y+x*(t*pi+1)/(pi-1))+lbeta((x*(t*pi+1))/(pi-1)
#             +(t*(2*pi-1)+1)/(t*(pi-1)),y+1/t)-log(factorial(y))-lgamma(x*(t*pi+1)/(pi-1))
#             -lbeta((t*(2*pi-1)+1)/t*(pi-1),1/t)))

my.fun<-function(pi,y,t){
        f1<-function(x){
             exp(1/pi)*((1/(1+t*x))^(1/t)*(t*x/(t*x+1))^y)-(exp(lgamma(y+x*(t*pi+1)/(pi-1))+lbeta((x*(t*pi+1))/(pi-1)
              +(t*(2*pi-1)+1)/(t*(pi-1)),y+1/t)-log(factorial(y))-lgamma(x*(t*pi+1)/(pi-1))
              -lbeta((t*(2*pi-1)+1)/t*(pi-1),1/t)))
}
return(f1)
}
H<-function(pi,y,t){
 integrate(my.fun(pi,y,t),0,Inf)
}

H(2,10,0.2)

####################END########################











