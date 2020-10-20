# rm(list = ls())
library(splines)
library(nlme)

bas=function(x, degree, knots){
  np=degree+length(knots)+1
  s=rep(0, np)
  rr=lapply(1:length(x),function(i){lapply(1:np,function(j){basis(x[i], degree, (j-1), knots)})})
  r=matrix(unlist(rr),nrow=length(x),ncol=np,byrow = T)
  return(r)
}

basis=function(x, degree, i, knots){
  if(degree==0){if(i==0){if(x<knots[1]){y=1}else{y=0}}
    else if(i==length(knots)){if(x>=tail(knots,n=1)){y=1}else{y=0}}
    else {if((x<knots[i+1])&(x>=knots[i])){y=1}else{y=0}}
  }
  else{
    k=length(knots);
    if((i>=min(degree+1,k+1))&(i<=max(degree,k))){
      if((knots[i]-knots[i-degree])==0) {temp1=1} else {temp1=(x-knots[i-degree])/(knots[i]-knots[i-degree])}}
    if((i>=min(degree,k))&(i<max(degree,k))){
      if((knots[i+1]-knots[i-degree+1])==0) {temp2=0} else {temp2=(knots[i+1]-x)/(knots[i+1]-knots[i-degree+1])}}
    
    # if(i>=min(degree+2,k+2)){
    #   if((knots[i]-knots[i-degree-1])==0) {temp1=1} else {temp1=(x-knots[i-degree-1])/(knots[i-1]-knots[i-degree-1])}
    #   if((knots[i]-knots[i-degree])==0) {temp2=0} else {temp2=(knots[i]-x)/(knots[i]-knots[i-degree])}}
    if(tail(knots,n=1)>knots[1]){c=(tail(knots,n=1)-knots[1])/(k-1)}else{c=1}
    
    if(i==0){y=(knots[i+1]-x)*basis(x, (degree-1), i, knots)/c}
    else if((i>=1)&(i<min(degree,k))){y=basis(x, (degree-1),(i-1), knots)+(knots[i+1]-x)*basis(x, (degree-1), i, knots)/c}
    else if(i==min(degree,k)){if(k==degree){y=basis(x, (degree-1),(i-1), knots)+basis(x, (degree-1),i, knots)}
      else if(k>degree){y=basis(x, (degree-1),(i-1), knots)+temp2*basis(x, (degree-1), i, knots)}}
    else if((i>=min(degree+1,k+1))&(i<max(degree,k))){y=temp1*basis(x, (degree-1), (i-1), knots)+temp2*basis(x, (degree-1), i, knots)}
    else if(i==max(degree,k)){if(k==degree){y=basis(x, (degree-1),(i-1), knots)+basis(x, (degree-1),i, knots)}
      else if(k>degree){y=basis(x, (degree-1),i, knots)+temp1*basis(x, (degree-1), (i-1), knots)}}
    else if((i>=max(degree+1,k+1))&(i<(degree+k))){y=(x-knots[i-degree])*basis(x, (degree-1), (i-1), knots)/c+basis(x, (degree-1), i, knots)}
    else if(i==(k+degree)){y=(x-knots[i-degree])*basis(x, (degree-1), (i-1), knots)/c}
  }
  return(y)
}