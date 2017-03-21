
expit<-function(eta) exp(eta)/(1+exp(eta))


one.sim<-function(beta=0.07,knot=1.8){
df<-data.frame(x=rnorm(1e5))
df$mu<-expit(df$x-5+beta*pmax(df$x-knot,0))
df$y<-rbinom(nrow(df),1,df$mu)

cases<-which(df$y==1)
controls<-sample(which(df$y==0),length(cases))

sdf<-df[c(cases,controls),]
sdf$wt<-ifelse(sdf$y==1,1,sum(df$y==0)/length(cases))

m0<-glm(y~x+I(pmin(x-knot,0)),data=sdf,family=binomial)
m1<-glm(y~x,data=sdf,family=binomial)
m2<-glm(y~x,data=sdf,family=quasibinomial,weights=wt)

c(coef(summary(m0))[3,3],coef(m1)[2],coef(m2)[2])
}


epsilons<-c(0,0.2,0.3,0.35,0.4,0.6)
rr<-lapply(epsilons, function(b) replicate(1000,one.sim(b,knot=1.8)))
errors<-sapply(rr, function(d) c(var(d[3,]),var(d[2,])+(mean(d[3,])-mean(d[2,]))^2))

matplot(epsilons,t(errors),pch=c(19,1),type="b",lty=1,ylim=range(0,errors),ylab="MSE")

kappa<-sapply(rr, function(d) -mean(d[1,])/sd(d[1,]))

power<-pnorm(qnorm(0.05)+pmax(0,kappa))


matplot(epsilons,t(errors)*1000,pch=c(1,19),col=c("red","black"),type="b",lty=1,ylim=range(0,errors*1000),ylab=expression("MSE"%*%1000),xaxt="n",xlab="Power")
axis(1,at=epsilons,labels=round(100*power))

save(epsilons, rr,errors,kappa,power,file="~/nearlytrue-linspline.rda")
