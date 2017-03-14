library(parallel)
r<-function(x) 1-1/(1+exp(x))

r0<-Vectorize(function(mu) integrate(function(x)(1-r(x))*dnorm(x,mu,1),lower=-Inf,upper=Inf)$value)
ur0<-Vectorize(function(mu) integrate(function(x)(1-r(x))*(x-mu)*dnorm(x,mu,1),lower=-Inf,upper=Inf)$value)


one.value<-function(epsilon){
#epsilon<-0.012
mu<-0
pdf.unnorm<-function(x,mu,epsilon) dnorm(x,mu,1)*exp(-epsilon*(1-1/r(x))*(x-mu)-epsilon*ur0(mu))
cdf.unnorm<-function(x,mu,epsilon) integrate(Vectorize(function(s) pdf.unnorm(s,mu,epsilon)),lower=-Inf,upper=x)$value
c.epsilon<-cdf.unnorm(Inf,mu,epsilon)
pdf<-function(x,mu,epsilon) pdf.unnorm(x,mu,epsilon)/c.epsilon
cdf<-function(x,mu,epsilon) cdf.unnorm(x,mu,epsilon)/c.epsilon
pdfratio<-function(x,mu,epsilon) exp(-epsilon*(1-1/r(x))*(x-mu)-epsilon*ur0(mu))

popnorm<-rnorm(1e7)
popratio<-pdfratio(popnorm,0,epsilon)
keepnorm<-rbinom(1e7,1,popratio/max(popratio))
popwrong<-popnorm[keepnorm==1]


bestapprox<-replicate(1000,{
	X<-sample(popwrong,10000)
	R<-rbinom(10000,1,r(X))
	X[R==0]<-NA

	ell0<-Vectorize(function(mu){
		sum(dnorm(X[R==1],mu,1,log=TRUE))+sum(R==0)*log(r0(mu))	
	})

	U0<-Vectorize(function(mu){
		sum((X-mu)[R==1])+ sum(R==0)*log(ur0(mu))
	})


	optim(0,ell0,U0,method="Brent",lower=-5,upper=5,control=list(fnscale=-1/100))$par
})

centering<-mean(bestapprox)
print(c(epsilon,centering))
r1<-function(mu) integrate(function(x) (1-r(x))*pdf(x,mu,epsilon),lower=-Inf,upper=Inf,rel.tol=1e-8 )$value


replicate(1000,{
	X<-sample(popwrong,10000)
	R<-rbinom(10000,1,r(X))
	X[R==0]<-NA

	ell0<-Vectorize(function(mu){
		sum(dnorm(X[R==1],mu,1,log=TRUE))+sum(R==0)*log(r0(mu))	
	})

	U0<-Vectorize(function(mu){
		sum((X-mu)[R==1])+ sum(R==0)*log(ur0(mu))
	})

	ell1<-Vectorize(function(mu){
		sum(log(pdf(X[R==1],mu,epsilon)))+sum(R==0)*log(r1(mu))	
	})

	c(mle=optim(0,ell0,U0,method="Brent",lower=-5,upper=5,control=list(fnscale=-1/100))$par,
	  ipw=sum(X[R==1]/r(X[R==1]))/sum(1/r(X[R==1])), 
	  ell0=ell0(centering), ell1=ell1(0)
	)
})

}


allresults<-lapply(c(0.000001,.003,0.006,0.009,0.012,0.015),one.value)

save(allresults, file="~/nearlytrue-toy-sim.rda")
