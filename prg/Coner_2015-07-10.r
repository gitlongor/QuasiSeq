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











