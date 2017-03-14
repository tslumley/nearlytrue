
beta<-1
expit<-function(eta) exp(eta)/(1+exp(eta))
logit<-function(p) log(p/(1-p))

POPN<-2e5
x<-rnorm(POPN)
mu<-expit(x*beta-beta^2/2-6)
pi0<- sum(mu)/(POPN-sum(mu))

y<-rbinom(POPN, 1, mu)

keep<-c(which(y==1), sample(which(y==0), sum(y==1)))

X<-cbind(1,x)
mu1<-expit(x*beta-beta^2/2-6-log(pi0))

V<- t(X)%*%(mu*(1-mu)*(X))/POPN

V1<- t(X)%*%(mu1*(1-mu1)*X*(mu+pi0))/sum(mu+pi0)

if0<- (X[keep,]*((y-mu)[keep])*ifelse(y[keep]==1,1,1/pi0))%*%solve(V)
if1<-(X[keep,]*((y-mu1)[keep]))%*%solve(V1)


sims<-function(epsilon,LOTS){
delta<-(if0-if1)[,2]
mustar<-mu[keep]*exp(delta*epsilon/sqrt(sum(mu[keep])))
mufn<-lowess(x[keep],mustar,f=1/10)
#plot(x[keep],mustar)
#lines(mufn,col="red")

nearlyzero<-1e-10
nearlyone<-1-nearlyzero


one.sim<-function(){
x<-rnorm(POPN)	
mm<-pmax(nearlyzero ,pmin(nearlyone, approxfun(mufn$x,mufn$y,rule=2)(x)))
y<-rbinom(POPN,1, mm)
cases<-which(y==1)
ctrls<-sample(which(y==0),sum(y))
yy<-y[c(cases,ctrls)]
xx<-x[c(cases,ctrls)]
c(
	coef(glm(y~x,family=binomial)),
	coef(glm(yy~xx,family=binomial)),
	coef(glm(yy~xx,family=binomial,weights=ifelse(yy==1,1,1/pi0)))
)
}


two.sim<-function(beta0,beta1){
x<-rnorm(POPN)	
mm<-pmax(nearlyzero,pmin(nearlyone,approxfun(mufn$x,mufn$y,rule=2)(x)))
y<-rbinom(POPN,1,mm)
cases<-which(y==1)
ctrls<-sample(which(y==0),sum(y))
yy<-y[c(cases,ctrls)]
xx<-x[c(cases,ctrls)]
llr<-sum(dbinom(yy,1,expit(xx*beta1+beta0),log=TRUE))-sum(dbinom(yy,1,expit(logit(mm[c(cases,ctrls)])-log(pi0)),log=TRUE))
}

rval<-replicate(LOTS,one.sim())
beta0<-mean(rval[3,])
beta1<-mean(rval[4,])
rval<-rbind(rval,replicate(LOTS,two.sim(beta0,beta1)))
rval
}

