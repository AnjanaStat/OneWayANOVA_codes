rm(list=ls())
install.packages("MASS")
# package MASS is needed for generating sample from multivariate normal distribution
library(MASS)
# by fun1 and fun2 we find crtical value of asymptotic Max-T test
fun1<-function(d1,d2,d3,d4,nu1,nu2,nu3,nu4,nu5,v2,v3,v4)
{
  # find estimate of the diagonal matrix D
  D=matrix(c(d1,0,0,0,0,d2,0,0,0,0,d3,0,0,0,0,d4),nrow=4,byrow=TRUE)
  s1=sqrt(nu1*nu3)*v2;s2=sqrt(nu2*nu4)*v3;s3=sqrt(nu3*nu5)*v4
  # find the estimate of the dispersion matrix
  S=matrix(c(d1^2,-s1,0,0,-s1,d2^2,-s2,0,0,-s2,d3^2,-s3,0,0,-s3,d4^2),nrow=4,byrow=TRUE)
  #find the estimate of covariance matrix
  P=solve(D)%*%S%*%solve(D);mu=c(0,0,0,0)
  # generate sample from multivariate normal distribution with mean vector mu and covariance matrix P
  T<-mvrnorm(1,mu,P)
  #find the max-T value
  A=max(T[1],T[2],T[3],T[4],na.rm = FALSE)
  return(A)
}
fun2<-function(d1,d2,d3,d4,nu1,nu2,nu3,nu4,nu5,v2,v3,v4)
{
  # find 500 max-T values from fun1
  x<-replicate(1000,fun1(d1,d2,d3,d4,nu1,nu2,nu3,nu4,nu5,v2,v3,v4))
  # sort them in increasing order
  y<-sort(x,decreasing=FALSE);c<-y[950]
  return(c)
}
fun3<-function(n1,n2,n3,n4,n5,mu1,mu2,mu3,mu4,mu5,v11,v12,v13,v14,v15)
{
  # generate sample from normal distribution
  g1<-rnorm(n1,mu1,sqrt(v11));g2<-rnorm(n2,mu2,sqrt(v12))
  g3<-rnorm(n3,mu3,sqrt(v13));g4<-rnorm(n4,mu4,sqrt(v14));g5<-rnorm(n5,mu5,sqrt(v15))
  # below steps are for finding max-T statistic values
  X1=mean(g1);X2=mean(g2);X3=mean(g3);X4=mean(g4); X5=mean(g5)
  s1=var(g1);s2=var(g2);s3=var(g3);s4=var(g4);s5=var(g5)
  v1=sqrt(s1/n1+s2/n2);v2=sqrt(s2/n2+s3/n3);v3=sqrt(s3/n3+s4/n4);v4=sqrt(s4/n4+s5/n5)
  T1=(X2-X1)/v1;T2=(X3-X2)/v2;T3=(X4-X3)/v3;T4=(X5-X4)/v4
  # find max-t statistic value
  A=max(T1,T2,T3,T4,na.rm=FALSE)
  # next steps are for finding critical value
  N=n1+n2+n3+n4+n5
  nu1=n1/N;nu2=n2/N;nu3=n3/N;nu4=n4/N;nu5=n5/N
  d1=sqrt(nu1*s2+nu2*s1);d2=sqrt(nu2*s3+nu3*s2);d3=sqrt(nu3*s4+nu4*s3);d4=sqrt(nu4*s5+nu5*s4)
  # find crtical value from fun1 and fun2
  out=fun2(d1,d2,d3,d4,nu1,nu2,nu3,nu4,nu5,s2,s3,s4)
  # count if the statistic value is greater than the critical value
  a=0
  if(A>out)
    a=a+1
  return(a)
}
fun4<-function(n1,n2,n3,n4,n5,mu1,mu2,mu3,mu4,mu5,v11,v12,v13,v14,v15)
{ 
  # find the number of times statistic value is greater than the crtical value among 1000 values
  out<-replicate(1000,fun3(n1,n2,n3,n4,n5,mu1,mu2,mu3,mu4,mu5,v11,v12,v13,v14,v15))
  p<-sum(out)/1000
  return(p)
}

##Size values of Table 2.5

p<-replicate(5,fun4(15,10,5,5,6,1,1,1,1,1,1,16,64,4,9))
p;alpha<-mean(p);alpha
p<-replicate(5,fun4(5,5,5,5,5,1,1,1,1,1,1,16,64,4,9))
p;alpha<-mean(p);alpha
p<-replicate(5,fun4(5,6,5,6,5,1,1,1,1,1,1,16,64,4,9))
p;alpha<-mean(p);alpha
p<-replicate(5,fun4(7,8,10,11,12,1,1,1,1,1,1,16,64,4,9))
p;alpha<-mean(p);alpha
p<-replicate(5,fun4(60,8,5,10,70,1,1,1,1,1,1,16,64,4,9))
p;alpha<-mean(p);alpha
p<-replicate(5,fun4(90,5,5,8,7,1,1,1,1,1,1,16,64,4,9))
p;alpha<-mean(p);alpha
p<-replicate(5,fun4(16,5,5,5,32,1,1,1,1,1,1,16,64,4,9))
p;alpha<-mean(p);alpha


p<-replicate(5,fun4(15,10,5,5,6,1,1,1,1,1,9,4,64,16,1))
p;alpha<-mean(p);alpha
p<-replicate(5,fun4(5,5,5,5,5,1,1,1,1,1,9,4,64,16,1))
p;alpha<-mean(p);alpha
p<-replicate(5,fun4(5,6,5,6,5,1,1,1,1,1,9,4,64,16,1))
p;alpha<-mean(p);alpha
p<-replicate(5,fun4(7,8,10,11,12,1,1,1,1,1,9,4,64,16,1))
p;alpha<-mean(p);alpha
p<-replicate(5,fun4(60,8,5,10,70,1,1,1,1,1,9,4,64,16,1))
p;alpha<-mean(p);alpha
p<-replicate(5,fun4(90,5,5,8,7,1,1,1,1,1,9,4,64,16,1))
p;alpha<-mean(p);alpha
p<-replicate(5,fun4(16,5,5,5,32,1,1,1,1,1,9,4,64,16,1))
p;alpha<-mean(p);alpha

