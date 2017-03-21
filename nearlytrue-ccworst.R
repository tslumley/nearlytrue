
beta<-1
expit<-function(eta) exp(eta)/(1+exp(eta))
logit<-function(p) log(p/(1-p))

#epsilon<-0.0004
#POPN<-2e5

one.sim<-function(epsilon,POPN=2e5,beta=1) replicate(10000,{
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
ww<-ifelse(yy==1,1,(POPN-sum(yy))/sum(yy))

llr<-sum(log(pxy[idx]/mean(pxy)))

m0<-glm(yy~xx,family=binomial)
m1<-glm(yy~xx,weights=ww,family=quasibinomial)
c(coef(m0)[2],coef(m1)[2],llr)
}
)


epsilons<-c(0,0.0001,0.00015, 0.0002,0.00023,0.00026,0.0003,0.0004)
rcc<-lapply(epsilons, one.sim)
save(epsilons, rcc, beta, "~/nearlytrue-ccworst.rda")

kappa<-sapply(rcc,function(d) c(2*mean(d[3,]),var(d[3,])))
rho<-sapply(rcc, function(d) cor(d[3,],d[2,]-d[1,]))
errors<-sapply(rcc,function(d) 1000*c(mean(d[1,]-d[2,])^2, var(d[2,])-var(d[1,])))
mse<-sapply(rcc,function(d) 1000*c(mean(d[1,]-d[2,])^2+var(d[1,]), var(d[2,])))


matplot(epsilons,t(mse),type="b",lty=1,pch=c(19,1),xlab=expression("% Power"~(rho==0.5)),ylab="MSE",xaxt="n",ylim=range(mse,0) )
powers<-pnorm(qnorm(0.05)+kappa[1,]/sqrt(kappa[2,]))
powers[1]<-0.05
axis(1,at=epsilons,labels=round(100*powers))
