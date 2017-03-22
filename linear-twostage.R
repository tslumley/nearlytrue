library(missreg3)
library(survey)
library(splines)

estfun.glm<-function (model) 
{
    xmat <- model.matrix(model)
    residuals(model, "working") * model$weights * xmat
}

rr<-replicate(500,{
df<-data.frame(x=rnorm(4000))
df$y<-df$x+rnorm(4000)
df$z<-df$x+rnorm(4000)
df$id<-1:nrow(df)
df$fullx<-df$x

df$zstrat<-cut(df$z, c(-Inf, qnorm(.05,s=sqrt(2)),qnorm(0.95,s=sqrt(2)),Inf))
df$insample<-df$id %in% c(which(df$zstrat=="(-Inf,-2.33]"), which(df$zstrat=="(2.33, Inf]"),sample(which(df$zstrat=="(-2.33,2.33]"),200))
df$x[!df$insample]<-NA

impmodel<-glm(x~ns(y,3)*zstrat+ns(z,3),data=subset(df,insample))

df$predx<-predict(impmodel,newdata=df)

calmodel<-glm(y~predx,data=df)

ef<-as.data.frame(estfun.glm(calmodel))
names(ef)<-c("ef1","ef2")

impmodel2<-glm(x~ns(y,3)*zstrat,data=subset(df,insample))

df$predx2<-predict(impmodel2,newdata=df)

calmodel2<-glm(y~predx2,data=df)

eff<-as.data.frame(estfun.glm(calmodel2))
names(eff)<-c("eff1","eff2")

df2<-cbind(df,ef,eff)

des<-twophase(id=list(~id,~id), strata=list(NULL,~zstrat), data=df2, subset=~insample)

cdes<-calibrate(des,formula=~(ef1+ef2)*zstrat,calfun="raking")
cdes2<-calibrate(des,formula=~(eff1+eff2)*zstrat,calfun="raking")


df$obstype.name<-ifelse(df$insample,"retro","strata")
yCuts<-c(-3,-2.33,-2,-1,0,1,2,2.33,3)
mz<-locsc2stg(y~x,~1,xstrata="zstrat",data=df,xs.includes=FALSE,method="ycutmeth",obstype.name="obstype.name",yCuts=yCuts,errdistn="normal")
mz2<-locsc2stg(y~x,~1,xstrata="zstrat",data=df,xs.includes=FALSE,method="direct",obstype.name="obstype.name",start=mz$coefficients,errdistn="normal")

c(
coef(svyglm(y~x,design=des))[2],
coef(svyglm(y~x,design=cdes))[2],
coef(svyglm(y~x,design=cdes2))[2],
summary(mz)$coef.table2[2,1],
summary(mz2)$coef.table2[2,1],
coef(glm(y~fullx,data=df))[2]
)
})

save(rr,file="~/lin-mle.rda")


jk<-function(data,fit){
	one.jk<-function(i){ 
		print(i)
	   refit<-locsc2stg(y~x,~1,xstrata="zstrat",data=data[-i,],
	   	xs.includes=FALSE,method="direct",obstype.name="obstype.name", print.progress=0,
	   	start=fit$coefficients,errdistn="normal", Qstart=fit$Qmat)
	   coef(refit)-coef(fit)
	}
t(sapply(1:nrow(data), one.jk))
}


jkvalues<-jk(df,mz2)

save(rr,jkvalues,file="~/lin-mle.rda")





epsilon<-0.05
one.sim<-function(epsilon){
df<-data.frame(x=rnorm(4000))
df$z<-df$x+rnorm(4000)
df$zstrat<-cut(df$z, c(-Inf, qnorm(.05,s=sqrt(2)),qnorm(0.95,s=sqrt(2)),Inf))
df$zmidx<-ifelse(as.numeric(df$zstrat)==2,df$x,0)
df$y<-df$x+rnorm(4000)-epsilon*df$zmidx
df$id<-1:nrow(df)
df$fullx<-df$x

df$insample<-df$id %in% c(which(df$zstrat=="(-Inf,-2.33]"), which(df$zstrat=="(2.33, Inf]"),sample(which(df$zstrat=="(-2.33,2.33]"),200))
df$x[!df$insample]<-NA

impmodel<-glm(x~ns(y,3)*zstrat+ns(z,3),data=subset(df,insample))

df$predx<-predict(impmodel,newdata=df)

calmodel<-glm(y~predx,data=df)

ef<-as.data.frame(estfun.glm(calmodel))
names(ef)<-c("ef1","ef2")

impmodel2<-glm(x~ns(y,3)*zstrat,data=subset(df,insample))

df$predx2<-predict(impmodel2,newdata=df)

calmodel2<-glm(y~predx2,data=df)

eff<-as.data.frame(estfun.glm(calmodel2))
names(eff)<-c("eff1","eff2")

df2<-cbind(df,ef,eff)

des<-twophase(id=list(~id,~id), strata=list(NULL,~zstrat), data=df2, subset=~insample)

#cdes<-calibrate(des,formula=~(ef1+ef2)*zstrat,calfun="raking")
cdes2<-calibrate(des,formula=~(eff1+eff2)*zstrat,calfun="raking")


df$obstype.name<-ifelse(df$insample,"retro","strata")
yCuts<-c(-3,-2.33,-2,-1,0,1,2,2.33,3)
mz<-locsc2stg(y~x,~1,xstrata="zstrat",data=df,xs.includes=FALSE,method="ycutmeth",obstype.name="obstype.name",yCuts=yCuts,errdistn="normal", print.progress=0)
mz2<-locsc2stg(y~x,~1,xstrata="zstrat",data=df,xs.includes=FALSE,method="direct",obstype.name="obstype.name",start=mz$coefficients,errdistn="normal", print.progress=0)
mzdiv<-tryCatch(
locsc2stg(y~x+zmidx,~1,xstrata="zstrat",data=df,xs.includes=FALSE,method="ycutmeth",obstype.name="obstype.name",yCuts=yCuts,errdistn="normal", print.progress=0),
error=function(e) NULL
)

cat(".")
flush.console()
c(
coef(svyglm(y~x,design=des))[2],
#coef(svyglm(y~x,design=cdes))[2],
coef(svyglm(y~x,design=cdes2))[2],
#summary(mz)$coef.table2[2,1],
summary(mz2)$coef.table2[2,1],
if (is.null(mzdiv)) NA else summary(mzdiv)$coef.table2[3,3]
#coef(glm(y~fullx,data=df))[2]
)}

epsilons<-c(0,0.025,0.05,0.06,0.07,0.08,0.1,0.15)
rr2<-lapply(epsilons, function(eps) {print(eps);replicate(1000, one.sim(eps))})

errors<-sapply(rr2, function(d) c(var(d[2,]),var(d[3,])+(mean(d[3,])-mean(d[2,]))^2))

kappa<-sapply(rr2, function(d) -mean(d[4,])/sd(d[4,]))

power<-pnorm(qnorm(0.05)+pmax(0,kappa))


matplot(epsilons,t(errors)*1000,pch=c(1,19),col=c("red","black"),type="b",lty=1,ylim=range(0,errors*1000),ylab=expression("MSE"%*%1000),xaxt="n",xlab="Power")
axis(1,at=epsilons,labels=round(100*power))

