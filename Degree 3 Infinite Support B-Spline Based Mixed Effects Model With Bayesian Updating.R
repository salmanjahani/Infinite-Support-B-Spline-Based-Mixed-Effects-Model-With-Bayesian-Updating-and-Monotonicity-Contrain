#' Code for Degree 3 Infinite SUpport B-Spline Based Mixed Effects Model With Bayesian Updating
#' Simulation
list.of.packages <- c("Matrix","nlme","rootSolve","splines","mvtnorm","gmm","tmvtnorm")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

# Load the packages
library(Matrix)
library(nlme)                    # Nonlinear Mixed Effects Model
library(rootSolve)               # Nonlinear Root Finding
library(splines)                 # Regression spline functions and classes
library(mvtnorm)
library(gmm)
library(tmvtnorm)

rm(list=ls())                    # Clean-up the memory
########################################### System Input ##############################################

setwd("D:/Job Hunt/Sample Code/B-Splines With Infinite Support")  #please set the correct directory 
source("Infinite Sup.R") 
endroot=15
threshold=30
iterations=100
##########################################
iit=0
counter=0
result_InfBS=matrix(0,iterations,3)
f=1
q=c()
while(1)
{
  iit=iit+1
  ########################################### Construct Data Groups ##############################################
  data <- data.frame(id = numeric(0), x = numeric(0), y = numeric(0))
  m=15
  for (i in 1:m) {
    x <- seq(0,10,length.out=30)
    w <-rnorm(1,1,0.25)
    y <- w*x^2+rnorm(length(x),0,0.5)               
    id <- rep(i, length(x))
    data <- rbind(data, cbind(id = id, x = x, y = y))
  }
  data <- groupedData(y ~ x | id, data)
  
  ########################################## Train Test Mean ###############################################
  n=15                              # first n-1 out of 104 available signals are selected to train the model
  # The n th signal is considered the new signal for which we are interested to do prediction
  train=data[data$id%in%(1:n-1),]
  test=data[data$id%in%n,]
  trains=lapply(1:(n-1),function(i){train$x[train$id==i]})
  trainy=lapply(1:(n-1),function(i){train$y[train$id==i]})
  for(i in 1:(n-1)){
    cr=trainy[[i]][trainy[[i]]==trainy[[i]][trainy[[i]]>=threshold][1]]
    trains[[i]]=trains[[i]][which(trainy[[i]]<=cr)]
    trainy[[i]]=trainy[[i]][which(trainy[[i]]<=cr)]
  }
  #### Truncate all data after threshold level
  cr=test$y[test$y==test$y[test$y>=threshold][1]]
  death=c()
  tr=test$x[which(test$y==cr)]
  death[n]=tr
  tstar=0.4*tr                ## Time of Bayesian Update
  for (i in 1:(n-1)) {
    x <- trains[[i]]
    y <- trainy[[i]]              
    id <- rep(i, length(x))
    data <- rbind(data, cbind(id = id, x = x, y = y))
  }
  train <- groupedData(y ~ x | id, data)
  tests=test$x[test$x<=tstar];m1=length(tests);testy=test$y[1:m1];  # Test Signal
  ########################################## Exploratory data visualization ###############################################
  dev.off()                        #Remove any previous plot
  for(i in 1:(n-1))                 #Plot the first 20 Evaporator outlet temperature signals
  {
    plot(trains[[i]],trainy[[i]],xlim=c(0,15),ylim=c(-2,50),type="l",pch=2,lwd=2,cex.axis=1.2,xlab='Time',ylab='Resistance',main='Generated Degradation Signals')
    par(new=T)
  }
  abline(h=threshold,lty=3)
  
  ##################################### B-Spline Mixed effects model ##########################
  #Fit B-spline mixed effect model
  knots=seq(2,8,length.out = 4)
  fitlme <- tryCatch(lme(y~bas(x, knots=seq(2,8,length.out = 4),degree=3)-1,random=~bas(x, knots=seq(2,8,length.out = 4),degree=3)-1|id,data=train,control=lmeControl(returnObject=TRUE,opt="optim",optimMethod = "SANN")), error = function(e) e)
  if(any(class(fitlme) == "error")==T){cat("LME fit Prob")}  #Check for any error in fitting
  
  ##################################### Specifiy Segments of Infinite SUpport B-SPlines ###########################
  seg=vector("list", length = 5)
  xx=seq(-2,13,length.out=80)
  yy=bas(xx, knots=seq(2,8,length.out = 4),degree=3)
  seg1=vector("list", length = 5)
  for(i in 1:length(xx)){
    if((xx[i]<knots[1])){seg1[[1]]=c(seg1[[1]],xx[i])}
    if((xx[i]>=knots[1])&(xx[i]<knots[2])){seg1[[2]]=c(seg1[[2]],xx[i])}
    if((xx[i]>=knots[2])&(xx[i]<knots[3])){seg1[[3]]=c(seg1[[3]],xx[i])}
    if((xx[i]>=knots[3])&(xx[i]<knots[4])){seg1[[4]]=c(seg1[[4]],xx[i])}
    # if((xx[i]>=knots[4])&(xx[i]<knots[5])){seg1[[5]]=c(seg1[[5]],xx[i])}
    # if((xx[i]>=knots[5])&(xx[i]<knots[6])){seg1[[6]]=c(seg1[[6]],xx[i])}
    if((xx[i]>=knots[4])){seg1[[5]]=c(seg1[[5]],xx[i])}
  }
  l=c();ll=list()
  for(i in 1:length(seg1)){
    ll[[i]]=length(l)+seq(1:length(seg1[[i]]))
    l=c(l,length(l)+seq(1:length(seg1[[i]])))
  }
  ##################################### Extract Coefficents of Each Segment ##############################
  yy=bas(xx, knots=seq(2,8,length.out = 4),degree=3)
  nc=ncol(yy)
  coef=lapply(1:5,function(i){lapply(1:nc,function(j){lm(yy[ll[[i]],j]~seg1[[i]]+I(seg1[[i]]^2)+I(seg1[[i]]^3))$coefficient})})
  coef2=lapply(1:5,function(i){do.call('rbind',coef[[i]])})
  
  # Extract parameters
  sigma2f = (fitlme$sigma)^2                      # Signal noise
  MUBf = fitlme$coefficients$fixed                # Mean vector
  SIGMABf = var(fitlme$coefficients$random$id)    # Covariance matrix
  SIGMABf=matrix(as.vector(SIGMABf),nrow=nrow(SIGMABf),byrow=T)
  
  ##################################### Bayesian Updating for New Unit ##############################
  Rp=testy          
  Rp = t(t(Rp))
  Zp = matrix(0,m1,8)
  Zp[,(1:8)]=bas(tests, knots=seq(2,8,length.out = 4),degree=3)[1:m1,]

  Baysupdate = function(tvec, obs, X0, P0, sig)
  {
    tp = tvec
    Rp = obs[1:length(tvec)]
    Rp = t(t(Rp))
    Zp = matrix(0,length(tp),8)
    Zp[,(1:8)]=bas(tests, knots=seq(2,8,length.out = 4),degree=3)[1:length(tp),]#############################################################################
    SIGMABp = solve(solve(P0)+t(Zp)%*%Zp/sig)
    MUBp = SIGMABp%*%(t(Zp)%*%Rp/sig+solve(P0)%*%X0)
    return(list(MUBp=MUBp,SIGMABp=SIGMABp))
  }
  result1=Baysupdate(tests,testy,MUBf,SIGMABf,sigma2f)
  MUBpf1=result1$MUBp
  
  ##################################### Bayesian Updating for New Unit Considering Monotonocity Constraints##############################
  tt=proc.time()
  N=2000
  samplepool=rmvnorm(N,result1$MUBp,result1$SIGMABp)
  eta=list()
  csamplepool=rep(NA,length(MUBpf2))
  csamplepool=list()
  for(i in 1:N){
    MU=samplepool[i,]
    eta=lapply(1:5,function(i){lapply(1:4,function(j){as.vector(t(MU)%*%coef2[[i]][,j])})})
    cs=matrix(,nrow=5,ncol=3)
    for(t in 1:5){
      if(t==1){
        if((-eta[[t]][[3]]<=3*eta[[t]][[4]]*knots[t])&(eta[[t]][[4]]>=0)){cs[t,1]=(eta[[t]][[2]]-eta[[t]][[3]]^2/eta[[t]][[4]]/3)>=0}else{cs[t,1]=TRUE}
        cs[t,2]=eta[[t]][[2]]+2*eta[[t]][[3]]*knots[t]+3*eta[[t]][[4]]*knots[t]^2>=0
        cs[t,3]=eta[[t]][[2]]+2*eta[[t]][[3]]*(0)+3*eta[[t]][[4]]*(0)^2>=0}
      else if(t==5){
        if((3*eta[[t]][[4]]*knots[t-1]<=-eta[[t]][[3]])&(eta[[t]][[4]]>0)){cs[t,1]=(eta[[t]][[2]]-eta[[t]][[3]]^2/eta[[t]][[4]]/3)>=0}else{cs[t,1]=TRUE}
        cs[t,2]=eta[[t]][[2]]+2*eta[[t]][[3]]*knots[t-1]+3*eta[[t]][[4]]*knots[t-1]^2>=0
        cs[t,3]=eta[[t]][[2]]+2*eta[[t]][[3]]*10+3*eta[[t]][[4]]*10^2>=0}
      else{
        if((3*eta[[t]][[4]]*knots[t-1]<=-eta[[t]][[3]])&(-eta[[t]][[3]]<=3*eta[[t]][[4]]*knots[t])&(eta[[t]][[4]]>=0)){cs[t,1]=(eta[[t]][[2]]-eta[[t]][[3]]^2/eta[[t]][[4]]/3)>=0}else{cs[t,1]=TRUE}
        cs[t,2]=eta[[t]][[2]]+2*eta[[t]][[3]]*knots[t-1]+3*eta[[t]][[4]]*knots[t-1]^2>=0
        cs[t,3]=eta[[4]][[2]]+2*eta[[4]][[3]]*knots[4]+3*eta[[4]][[4]]*knots[4]^2>=0
      }
    }
    if(all(cs)){csamplepool[[i]]=MU}
  }
  # nrow(do.call("rbind",csamplepool))
  MUBpf2=apply(do.call("rbind",csamplepool),2,mean)
  tt=proc.time()-tt
  ################################# Visual illustration ########################################
  X=seq(0,11,length.out = 60)
  priorfit=bas(X, knots=seq(2,8,length.out = 4),degree=3)%*%as.matrix(MUBf)
  ####Posterior Fit
  posteriorfit1=bas(X, knots=seq(2,8,length.out = 4),degree=3)%*%as.matrix(MUBpf1)
  posteriorfit2=bas(X, knots=seq(2,8,length.out = 4),degree=3)%*%as.matrix(MUBpf2)
  # dev.off()
  plot(test$x,test$y,xlim=c(0,11),ylim=c(0,100),type="l",lwd=2,xlab='Time',ylab='',main='Online Updating for the In-field Unit',cex.axis=1.2) ## True
  ##Current observations of new signal
  points(tests,testy,xlim=c(0,11),ylim=c(0,100),pch=16,cex=2,xlab="",ylab="",main="",cex.axis=1.2)
  par(new=TRUE)
  ## Posterior and Prior predictions
  plot(X,priorfit,type="l",lwd=2,lty=3,col="blue",xlab="",ylab="",main="",xlim=c(0,11),ylim=c(0,100),cex.axis=1.2)
  par(new=TRUE)
  plot(X,posteriorfit1,type="l",lwd=2,lty=2,col="red",xlab="",ylab="",main="",xlim=c(0,11),ylim=c(0,100),cex.axis=1.2)
  par(new=TRUE)
  plot(X,posteriorfit2,type="l",lwd=2,lty=2,col="green",xlab="",ylab="",main="",xlim=c(0,11),ylim=c(0,100),cex.axis=1.2)
  legend("bottomright",col=c("black","blue","red","green"),legend=c("True","Prior","Posterior-No Constraint","Posterior-With Constraint"),lty=c(1,3,2,2),lwd=3.5,cex=0.8)
  
  #################################### Find RUL ###############################################
  zt=function(t) bas(t, knots=seq(2,8,length.out = 4),degree=3)
  resME=Vectorize(function(t)  (zt(t))%*%as.matrix(MUBpf2) - threshold )
  rootME=try(uniroot.all(resME,c(1,endroot)))
  if(inherits(rootME, "try-error")){iit=iit-1;next}
  if(length(rootME)==0){iit=iit-1;next}
  result_InfBS[iit,1:3]=c(rootME[1],death[n],tests[m1])
  
  cat(" Iteration ",iit,"is reached\n")
  if (iit==iterations) {break}
  
}
q2=abs(result_InfBS[,2]-result_InfBS[,1])
boxplot(q2,main="50% Observation percentile-Model Setting II ")
save(q2,file='settingII-0.5.RData')
save(result_InfBS,file="infbs04.RData")

