
beta<-1
expit<-function(eta) exp(eta)/(1+exp(eta))
logit<-function(p) log(p/(1-p))

#epsilon<-0.0004
#POPN<-2e5

epsilons<-c(0,0.0001,0.00015, 0.0002,0.00023,0.00026,0.0003,0.0004)


one.sim<-function(epsilon,POPN=2e5,beta=1) {
x<-rnorm(POPN)
mu<-expit(x*beta-5)
pi0<- sum(mu)/(POPN-sum(mu))

y<-rbinom(POPN, 1, mu)

keep<-c(which(y==1), sample(which(y==0), 2*sum(y==1)))

X<-cbind(1,x)
mu1<-expit(x*beta-5-log(pi0))

V<- t(X)%*%(mu*(1-mu)*(X))/POPN

V1<- t(X)%*%(mu1*(1-mu1)*X*(mu+pi0))/sum(mu+pi0)

Xk<-X[keep,]
xk<-x[keep]
yk<-y[keep]
muk<-mu[keep]
mu1k<-mu1[keep]

if0<- (Xk*((yk-muk))*ifelse(yk==1,1,1/pi0))%*%solve(V)
if1<-(Xk*((yk-mu1k)))%*%solve(V1)

delta<-(if0-if1)[,2]

pxy<-exp(delta*epsilon)

id1<-sample(which(yk==1),sum(yk==1)/2,replace=TRUE,prob=pxy[yk==1])
id0<-sample(which(yk==0),sum(yk==1)/2,replace=TRUE,prob=pxy[yk==0])


idx<-c(id1,id0)
xx<-xk[idx]
yy<-yk[idx]

cbind(xx,yy)
}

lots.sim<-function(epsilon) do.call(rbind, replicate(100, one.sim(epsilon)))

power.sim<-function(epsilon){
	df<-as.data.frame(lots.sim(epsilon))
	 
	l<-gam(yy~s(xx),data=df,family=binomial)
	etafn<-approxfun(df$x, logit(fitted(l)),rule=2)
	coefs<-coef(glm(yy~xx,data=df,family=binomial))
    replicate(1000,{
	   	df1<-as.data.frame(one.sim(epsilon))
		mline<-glm(yy~0+offset(coefs[1]+coefs[2]*xx),family=binomial,data=df1)
		mcurve<-glm(yy~0+offset(etafn(xx)),family=binomial,data=df1)
		logLik(mcurve)-logLik(mline)
	})
	
}


rcond<-lapply(epsilons, power.sim)

condkappas<-sapply(rcond,function(x) c(mean(x)*2,var(x)))
condpowers<-pnorm(qnorm(0.05)+condkappas[1,])
save(condkappas,condpowers,epsilons,rcond,file="~/nearlytrue-ccworst-conditional.rda")


load("~/nearlytrue-ccworst.rda")
kappa<-sapply(rcc,function(d) c(2*mean(d[3,]),var(d[3,])))
errors<-sapply(rcc,function(d) 1000*c(mean(d[1,]-d[2,])^2, var(d[2,])-var(d[1,])))
mse<-sapply(rcc,function(d) 1000*c(mean(d[1,]-d[2,])^2+var(d[1,]), var(d[2,])))
 matplot(epsilons,t(mse),type="b",lty=1,pch=c(19,1),xlab=expression("% Power"~(list("joint",rho==0.5))),ylab="MSE",xaxt="n",ylim=range(mse,0) )
 powers<-pnorm(qnorm(0.05)+kappa[1,])
 axis(1,at=epsilons,labels=round(100*powers))
 axis(3,at=epsilons,labels=round(100*condpowers))
 mtext(expression("% Power"~(conditional)),side=3,line=2)

## show the worst-case example 
df<-as.data.frame(lots.sim(epsilons[6])) 
l<-gam(yy~s(xx),data=df,family=binomial,weights=ifelse(yy==1,1,2e5*100/nrow(df)))
i<-order(df$xx)
plot(df$xx[i],logit(fitted(l)[i]),type="l",lwd=2,xlab="x",ylab="P[Y=1]")

l0s<-glm(yy~xx,data=df,family=quasibinomial,weights=ifelse(yy==1,1,2e5*100/nrow(df)))
abline(coef(l0s)[1],coef(l0s)[2],lty=2,lwd=2,col="blue")

