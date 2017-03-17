
expit<-function(eta) exp(eta)/(1+exp(eta))


one.sim<-function(beta=0.07){
df<-data.frame(x=rnorm(1e5))
df$x2<-df$x^2
df$mu<-expit(df$x-5+df$x2*beta)
df$y<-rbinom(nrow(df),1,df$mu)

cases<-which(df$y==1)
controls<-sample(which(df$y==0),length(cases))

sdf<-df[c(cases,controls),]
sdf$wt<-ifelse(sdf$y==1,1,sum(df$y==0)/length(cases))

m0<-glm(y~x+x2,data=df,family=binomial)
m1<-glm(y~x,data=sdf,family=binomial)
m2<-glm(y~x,data=sdf,family=quasibinomial,weights=wt)

c(coef(summary(m0))[3,3],coef(m1)[2],coef(m2)[2])
}
