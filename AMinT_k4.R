rm(list=ls())
install.packages("MASS")
# package MASS is needed for generating sample from multivariate normal distribution
library(MASS)
# by fun1 and fun2 we find critical value of asymptotic Max-T test
fun1<-function(d1,d2,d3,nu1,nu2,nu3,nu4,v2,v3)
{
  # find the estimate of the diagonal matrix D
  D=matrix(c(d1,0,0,0,d2,0,0,0,d3),nrow=3,byrow=TRUE)
  s1=sqrt(nu1*nu3)*v2;s2=sqrt(nu2*nu4)*v3
  # find the estimate of the dispersion matrix
  S=matrix(c(d1^2,-s1,0,-s1,d2^2,-s2,0,-s2,d3^2),nrow=3,byrow=TRUE)
  # find the covariance matrix
  P=solve(D)%*%S%*%solve(D);mu=c(0,0,0)
  # generate sample from multivariate normal distribution with mean vector mu and covariance matrix P
  T<-mvrnorm(1,mu,P)
  # find min-T statistic value
  A=min(T[1],T[2],T[3],na.rm = FALSE)
  return(A)
}
fun2<-function(d1,d2,d3,nu1,nu2,nu3,nu4,v2,v3)
{
  # find 500 max-T statistic values
  x<-replicate(1000,fun1(d1,d2,d3,nu1,nu2,nu3,nu4,v2,v3))
  # arrange them in increasing order
  y<-sort(x,decreasing=FALSE)
  # find the critical value
  c<-y[950]
  return(c)
}

fun3<-function(n1,n2,n3,n4,mu1,mu2,mu3,mu4,v11,v12,v13,v14)
{
  # generate sample from normal distribution
  g1<-rnorm(n1,mu1,sqrt(v11)); g2<-rnorm(n2,mu2,sqrt(v12));g3<-rnorm(n3,mu3,sqrt(v13))
  g4<-rnorm(n4,mu4,sqrt(v14))
  # below steps are for finding max-T statistics values
  X1=mean(g1);X2=mean(g2);X3=mean(g3);X4=mean(g4)
  S1=var(g1);S2=var(g2);S3=var(g3);S4=var(g4)
  v1<-sqrt(S1/n1+S2/n2);v2<-sqrt(S2/n2+S3/n3);v3<-sqrt(S3/n3+S4/n4)
  T1=(X2-X1)/v1;T2=(X3-X2)/v2;T3=(X4-X3)/v3
  # find min-T statistic value
  T=min(T1,T2,T3,na.rm = FALSE)
  N=n1+n2+n3+n4
  nu1=n1/N;nu2=n2/N;nu3=n3/N;nu4=n4/N
  d1=sqrt(nu1*S2+nu2*S1);d2=sqrt(nu2*S3+nu3*S2);d3=sqrt(nu3*S4+nu4*S3)
  # find critical value
  out=fun2(d1,d2,d3,nu1,nu2,nu3,nu4,S2,S3)
  # count if the statistic value is greater than the critical value
  a=0
  if(T>out)
    a=a+1
  return(a)
}
fun4<-function(n1,n2,n3,n4,mu1,mu2,mu3,mu4,v11,v12,v13,v14)
{
  # find the number of times statistic value is greater than the crtical value among 1000 values
  out<-replicate(1000,fun3(n1,n2,n3,n4,mu1,mu2,mu3,mu4,v11,v12,v13,v14))
  p<-sum(out)/1000
  return(p)
}

#Size values of Table 2.2

p<-replicate(5,fun4(5,5,7,6,1,1,1,1,10^2,10^2,10^2,10^2))
p;alpha=mean(p);alpha
p<-replicate(5,fun4(6,8,6,8,1,1,1,1,10^2,10^2,10^2,10^2))
p;alpha=mean(p);alpha
p<-replicate(5,fun4(10,5,6,8,1,1,1,1,10^2,10^2,10^2,10^2))
p;alpha=mean(p);alpha
p<-replicate(5,fun4(15,10,10,5,1,1,1,1,10^2,10^2,10^2,10^2))
p;alpha=mean(p);alpha
p<-replicate(5,fun4(30,5,30,5,1,1,1,1,10^2,10^2,10^2,10^2))
p;alpha=mean(p);alpha
p<-replicate(5,fun4(15,5,7,6,1,1,1,1,10^2,10^2,10^2,10^2))
p;alpha=mean(p);alpha
p<-replicate(5,fun4(6,6,5,18,1,1,1,1,10^2,10^2,10^2,10^2))
p;alpha=mean(p);alpha
p<-replicate(5,fun4(16,5,6,6,1,1,1,1,10^2,10^2,10^2,10^2))
p;alpha=mean(p);alpha
p<-replicate(5,fun4(29,29,29,29,1,1,1,1,10^2,10^2,10^2,10^2))
p;alpha=mean(p);alpha

p<-replicate(5,fun4(5,5,7,6,1,1,1,1,1,4^2,6^2,8^2))
p;alpha=mean(p);alpha
p<-replicate(5,fun4(6,8,6,8,1,1,1,1,1,4^2,6^2,8^2))
p;alpha=mean(p);alpha
p<-replicate(5,fun4(10,5,6,8,1,1,1,1,1,4^2,6^2,8^2))
p;alpha=mean(p);alpha
p<-replicate(5,fun4(15,10,10,5,1,1,1,1,1,4^2,6^2,8^2))
p;alpha=mean(p);alpha
p<-replicate(5,fun4(30,5,30,5,1,1,1,1,1,4^2,6^2,8^2))
p;alpha=mean(p);alpha
p<-replicate(5,fun4(15,5,7,6,1,1,1,1,1,4^2,6^2,8^2))
p;alpha=mean(p);alpha
p<-replicate(5,fun4(6,6,5,18,1,1,1,1,1,4^2,6^2,8^2))
p;alpha=mean(p);alpha
p<-replicate(5,fun4(16,5,6,6,1,1,1,1,1,4^2,6^2,8^2))
p;alpha=mean(p);alpha
p<-replicate(5,fun4(29,29,29,29,1,1,1,1,1,4^2,6^2,8^2))
p;alpha=mean(p);alpha



p<-replicate(5,fun4(5,5,7,6,1,1,1,1,8^2,6^2,4^2,1))
p;alpha=mean(p);alpha
p<-replicate(5,fun4(6,8,6,8,1,1,1,1,8^2,6^2,4^2,1))
p;alpha=mean(p);alpha
p<-replicate(5,fun4(10,5,6,8,1,1,1,1,8^2,6^2,4^2,1))
p;alpha=mean(p);alpha
p<-replicate(5,fun4(15,10,10,5,1,1,1,1,8^2,6^2,4^2,1))
p;alpha=mean(p);alpha
p<-replicate(5,fun4(30,5,30,5,1,1,1,1,8^2,6^2,4^2,1))
p;alpha=mean(p);alpha
p<-replicate(5,fun4(15,5,7,6,1,1,1,1,8^2,6^2,4^2,1))
p;alpha=mean(p);alpha
p<-replicate(5,fun4(6,6,5,18,1,1,1,1,8^2,6^2,4^2,1))
p;alpha=mean(p);alpha
p<-replicate(5,fun4(16,5,6,6,1,1,1,1,8^2,6^2,4^2,1))
p;alpha=mean(p);alpha
p<-replicate(5,fun4(29,29,29,29,1,1,1,1,8^2,6^2,4^2,1))
p;alpha=mean(p);alpha

p<-replicate(5,fun4(5,5,7,6,1,1,1,1,3^2,4^2,5^2,6^2))
p;alpha=mean(p);alpha
p<-replicate(5,fun4(6,8,6,8,1,1,1,1,3^2,4^2,5^2,6^2))
p;alpha=mean(p);alpha
p<-replicate(5,fun4(10,5,6,8,1,1,1,1,3^2,4^2,5^2,6^2))
p;alpha=mean(p);alpha
p<-replicate(5,fun4(15,10,10,5,1,1,1,1,3^2,4^2,5^2,6^2))
p;alpha=mean(p);alpha
p<-replicate(5,fun4(30,5,30,5,1,1,1,1,3^2,4^2,5^2,6^2))
p;alpha=mean(p);alpha
p<-replicate(5,fun4(15,5,7,6,1,1,1,1,3^2,4^2,5^2,6^2))
p;alpha=mean(p);alpha
p<-replicate(5,fun4(6,6,5,18,1,1,1,1,3^2,4^2,5^2,6^2))
p;alpha=mean(p);alpha
p<-replicate(5,fun4(16,5,6,6,1,1,1,1,3^2,4^2,5^2,6^2))
p;alpha=mean(p);alpha
p<-replicate(5,fun4(29,29,29,29,1,1,1,1,3^2,4^2,5^2,6^2))
p;alpha=mean(p);alpha

p<-replicate(5,fun4(5,5,7,6,1,1,1,1,6^2,5^2,4^2,3^2))
p;alpha=mean(p);alpha
p<-replicate(5,fun4(6,8,6,8,1,1,1,1,3^2,4^2,5^2,6^2))
p;alpha=mean(p);alpha
p<-replicate(5,fun4(10,5,6,8,1,1,1,1,3^2,4^2,5^2,6^2))
p;alpha=mean(p);alpha
p<-replicate(5,fun4(15,10,10,5,1,1,1,1,3^2,4^2,5^2,6^2))
p;alpha=mean(p);alpha
p<-replicate(5,fun4(30,5,30,5,1,1,1,1,3^2,4^2,5^2,6^2))
p;alpha=mean(p);alpha
p<-replicate(5,fun4(15,5,7,6,1,1,1,1,3^2,4^2,5^2,6^2))
p;alpha=mean(p);alpha
p<-replicate(5,fun4(6,6,5,18,1,1,1,1,3^2,4^2,5^2,6^2))
p;alpha=mean(p);alpha
p<-replicate(5,fun4(16,5,6,6,1,1,1,1,3^2,4^2,5^2,6^2))
p;alpha=mean(p);alpha
p<-replicate(5,fun4(29,29,29,29,1,1,1,1,3^2,4^2,5^2,6^2))
p;alpha=mean(p);alpha


p<-replicate(5,fun4(5,5,7,6,1,1,1,1,20^2,25^2,70^2,40^2))
p;alpha=mean(p);alpha
p<-replicate(5,fun4(6,8,6,8,1,1,1,1,20^2,25^2,70^2,40^2))
p;alpha=mean(p);alpha
p<-replicate(5,fun4(10,5,6,8,1,1,1,1,20^2,25^2,70^2,40^2))
p;alpha=mean(p);alpha
p<-replicate(5,fun4(15,10,10,5,1,1,1,1,20^2,25^2,70^2,40^2))
p;alpha=mean(p);alpha
p<-replicate(5,fun4(30,5,30,5,1,1,1,1,20^2,25^2,70^2,40^2))
p;alpha=mean(p);alpha
p<-replicate(5,fun4(15,5,7,6,1,1,1,1,20^2,25^2,70^2,40^2))
p;alpha=mean(p);alpha
p<-replicate(5,fun4(6,6,5,18,1,1,1,1,20^2,25^2,70^2,40^2))
p;alpha=mean(p);alpha
p<-replicate(5,fun4(16,5,6,6,1,1,1,1,20^2,25^2,70^2,40^2))
p;alpha=mean(p);alpha
p<-replicate(5,fun4(29,29,29,29,1,1,1,1,20^2,25^2,70^2,40^2))
p;alpha=mean(p);alpha


p<-replicate(5,fun4(27,28,29,30,400,500,600,700,20^2,25^2,70^2,40^2))
p;alpha=mean(p);alpha
p<-replicate(5,fun4(29,29,29,29,400,500,600,700,20^2,25^2,70^2,40^2))
p;alpha=mean(p);alpha

##Power values of Table 2.9

p<-replicate(5,fun4(n1=20,n2=25,n3=30,n4=20,mu1=1,mu2=1,mu3=1,mu4=1,v11=1,v12=1,v13=1,v14=1))
alpha=mean(p);alpha
p<-replicate(5,fun4(n1=20,n2=25,n3=30,n4=20,mu1=1,mu2=1.1,mu3=1.2,mu4=1.3,v11=1,v12=1,v13=1,v14=1))
alpha=mean(p);alpha
p<-replicate(5,fun4(n1=20,n2=25,n3=30,n4=20,mu1=1*1.2,mu2=1.1*1.2,mu3=1.2*1.2,mu4=1.3*1.2,v11=1,v12=1,v13=1,v14=1))
alpha=mean(p);alpha
p<-replicate(5,fun4(n1=20,n2=25,n3=30,n4=20,mu1=1*1.4,mu2=1.1*1.4,mu3=1.2*1.4,mu4=1.3*1.4,v11=1,v12=1,v13=1,v14=1))
alpha=mean(p);alpha
p<-replicate(5,fun4(n1=20,n2=25,n3=30,n4=20,mu1=1*1.6,mu2=1.1*1.6,mu3=1.2*1.6,mu4=1.3*1.6,v11=1,v12=1,v13=1,v14=1))
alpha=mean(p);alpha
p<-replicate(5,fun4(n1=20,n2=25,n3=30,n4=20,mu1=1*1.8,mu2=1.1*1.8,mu3=1.2*1.8,mu4=1.3*1.8,v11=1,v12=1,v13=1,v14=1))
alpha=mean(p);alpha
p<-replicate(5,fun4(n1=20,n2=25,n3=30,n4=20,mu1=1*2,mu2=1.1*2,mu3=1.2*2,mu4=1.3*2,v11=1,v12=1,v13=1,v14=1))
alpha=mean(p);alpha
p<-replicate(5,fun4(n1=20,n2=25,n3=30,n4=20,mu1=1*2.2,mu2=1.1*2.2,mu3=1.2*2.2,mu4=1.3*2.2,v11=1,v12=1,v13=1,v14=1))
alpha=mean(p);alpha
p<-replicate(5,fun4(n1=20,n2=25,n3=30,n4=20,mu1=1*2.4,mu2=1.1*2.4,mu3=1.2*2.4,mu4=1.3*2.4,v11=1,v12=1,v13=1,v14=1))
alpha=mean(p);alpha
p<-replicate(5,fun4(n1=20,n2=25,n3=30,n4=20,mu1=1*2.6,mu2=1.1*2.6,mu3=1.2*2.6,mu4=1.3*2.6,v11=1,v12=1,v13=1,v14=1))
alpha=mean(p);alpha
p<-replicate(5,fun4(n1=20,n2=25,n3=30,n4=20,mu1=1*2.8,mu2=1.1*2.8,mu3=1.2*2.8,mu4=1.3*2.8,v11=1,v12=1,v13=1,v14=1))
alpha=mean(p);alpha
p<-replicate(5,fun4(n1=20,n2=25,n3=30,n4=20,mu1=1*3,mu2=1.1*3,mu3=1.2*3,mu4=1.3*3,v11=1,v12=1,v13=1,v14=1))
alpha=mean(p);alpha
p<-replicate(5,fun4(n1=20,n2=25,n3=30,n4=20,mu1=1*3.2,mu2=1.1*3.2,mu3=1.2*3.2,mu4=1.3*3.2,v11=1,v12=1,v13=1,v14=1))
alpha=mean(p);alpha
p<-replicate(5,fun4(n1=20,n2=25,n3=30,n4=20,mu1=1*3.4,mu2=1.1*3.4,mu3=1.2*3.4,mu4=1.3*3.4,v11=1,v12=1,v13=1,v14=1))
alpha=mean(p);alpha
p<-replicate(5,fun4(n1=20,n2=25,n3=30,n4=20,mu1=1*3.6,mu2=1.1*3.6,mu3=1.2*3.6,mu4=1.3*3.6,v11=1,v12=1,v13=1,v14=1))
alpha=mean(p);alpha
p<-replicate(5,fun4(n1=20,n2=25,n3=30,n4=20,mu1=1*3.8,mu2=1.1*3.8,mu3=1.2*3.8,mu4=1.3*3.8,v11=1,v12=1,v13=1,v14=1))
alpha=mean(p);alpha
p<-replicate(5,fun4(n1=20,n2=25,n3=30,n4=20,mu1=1*4,mu2=1.1*4,mu3=1.2*4,mu4=1.3*4,v11=1,v12=1,v13=1,v14=1))
alpha=mean(p);alpha
p<-replicate(5,fun4(n1=20,n2=25,n3=30,n4=20,mu1=1*4.2,mu2=1.1*4.2,mu3=1.2*4.2,mu4=1.3*4.2,v11=1,v12=1,v13=1,v14=1))
alpha=mean(p);alpha
p<-replicate(5,fun4(n1=20,n2=25,n3=30,n4=20,mu1=1*4.4,mu2=1.1*4.4,mu3=1.2*4.4,mu4=1.3*4.4,v11=1,v12=1,v13=1,v14=1))
alpha=mean(p);alpha
p<-replicate(5,fun4(n1=20,n2=25,n3=30,n4=20,mu1=1*4.6,mu2=1.1*4.6,mu3=1.2*4.6,mu4=1.3*4.6,v11=1,v12=1,v13=1,v14=1))
alpha=mean(p);alpha
p<-replicate(5,fun4(n1=20,n2=25,n3=30,n4=20,mu1=1*4.8,mu2=1.1*4.8,mu3=1.2*4.8,mu4=1.3*4.8,v11=1,v12=1,v13=1,v14=1))
alpha=mean(p);alpha
p<-replicate(5,fun4(n1=20,n2=25,n3=30,n4=20,mu1=1*5,mu2=1.1*5,mu3=1.2*5,mu4=1.3*5,v11=1,v12=1,v13=1,v14=1))
alpha=mean(p);alpha


p<-replicate(5,fun4(n1=40,n2=60,n3=50,n4=70,mu1=1,mu2=1,mu3=1,mu4=1,v11=1,v12=1,v13=1,v14=1))
alpha=mean(p);alpha
p<-replicate(5,fun4(n1=40,n2=60,n3=50,n4=70,mu1=1,mu2=1.1,mu3=1.2,mu4=1.3,v11=1,v12=1,v13=1,v14=1))
alpha=mean(p);alpha
p<-replicate(5,fun4(n1=40,n2=60,n3=50,n4=70,mu1=1*1.2,mu2=1.1*1.2,mu3=1.2*1.2,mu4=1.3*1.2,v11=1,v12=1,v13=1,v14=1))
alpha=mean(p);alpha
p<-replicate(5,fun4(n1=40,n2=60,n3=50,n4=70,mu1=1*1.4,mu2=1.1*1.4,mu3=1.2*1.4,mu4=1.3*1.4,v11=1,v12=1,v13=1,v14=1))
alpha=mean(p);alpha
p<-replicate(5,fun4(n1=40,n2=60,n3=50,n4=70,mu1=1*1.6,mu2=1.1*1.6,mu3=1.2*1.6,mu4=1.3*1.6,v11=1,v12=1,v13=1,v14=1))
alpha=mean(p);alpha
p<-replicate(5,fun4(n1=40,n2=60,n3=50,n4=70,mu1=1*1.8,mu2=1.1*1.8,mu3=1.2*1.8,mu4=1.3*1.8,v11=1,v12=1,v13=1,v14=1))
alpha=mean(p);alpha
p<-replicate(5,fun4(n1=40,n2=60,n3=50,n4=70,mu1=1*2,mu2=1.1*2,mu3=1.2*2,mu4=1.3*2,v11=1,v12=1,v13=1,v14=1))
alpha=mean(p);alpha
p<-replicate(5,fun4(n1=40,n2=60,n3=50,n4=70,mu1=1*2.2,mu2=1.1*2.2,mu3=1.2*2.2,mu4=1.3*2.2,v11=1,v12=1,v13=1,v14=1))
alpha=mean(p);alpha
p<-replicate(5,fun4(n1=40,n2=60,n3=50,n4=70,mu1=1*2.4,mu2=1.1*2.4,mu3=1.2*2.4,mu4=1.3*2.4,v11=1,v12=1,v13=1,v14=1))
alpha=mean(p);alpha
p<-replicate(5,fun4(n1=40,n2=60,n3=50,n4=70,mu1=1*2.6,mu2=1.1*2.6,mu3=1.2*2.6,mu4=1.3*2.6,v11=1,v12=1,v13=1,v14=1))
alpha=mean(p);alpha
p<-replicate(5,fun4(n1=40,n2=60,n3=50,n4=70,mu1=1*2.8,mu2=1.1*2.8,mu3=1.2*2.8,mu4=1.3*2.8,v11=1,v12=1,v13=1,v14=1))
alpha=mean(p);alpha
p<-replicate(5,fun4(n1=40,n2=60,n3=50,n4=70,mu1=1*3,mu2=1.1*3,mu3=1.2*3,mu4=1.3*3,v11=1,v12=1,v13=1,v14=1))
alpha=mean(p);alpha
p<-replicate(5,fun4(n1=40,n2=60,n3=50,n4=70,mu1=1*3.2,mu2=1.1*3.2,mu3=1.2*3.2,mu4=1.3*3.2,v11=1,v12=1,v13=1,v14=1))
alpha=mean(p);alpha
p<-replicate(5,fun4(n1=40,n2=60,n3=50,n4=70,mu1=1*3.4,mu2=1.1*3.4,mu3=1.2*3.4,mu4=1.3*3.4,v11=1,v12=1,v13=1,v14=1))
alpha=mean(p);alpha
p<-replicate(5,fun4(n1=40,n2=60,n3=50,n4=70,mu1=1*3.6,mu2=1.1*3.6,mu3=1.2*3.6,mu4=1.3*3.6,v11=1,v12=1,v13=1,v14=1))
alpha=mean(p);alpha
p<-replicate(5,fun4(n1=40,n2=60,n3=50,n4=70,mu1=1*3.8,mu2=1.1*3.8,mu3=1.2*3.8,mu4=1.3*3.8,v11=1,v12=1,v13=1,v14=1))
alpha=mean(p);alpha
p<-replicate(5,fun4(n1=40,n2=60,n3=50,n4=70,mu1=1*4,mu2=1.1*4,mu3=1.2*4,mu4=1.3*4,v11=1,v12=1,v13=1,v14=1))
alpha=mean(p);alpha
p<-replicate(5,fun4(n1=40,n2=60,n3=50,n4=70,mu1=1*4.2,mu2=1.1*4.2,mu3=1.2*4.2,mu4=1.3*4.2,v11=1,v12=1,v13=1,v14=1))
alpha=mean(p);alpha
p<-replicate(5,fun4(n1=40,n2=60,n3=50,n4=70,mu1=1*4.4,mu2=1.1*4.4,mu3=1.2*4.4,mu4=1.3*4.4,v11=1,v12=1,v13=1,v14=1))
alpha=mean(p);alpha
p<-replicate(5,fun4(n1=40,n2=60,n3=50,n4=70,mu1=1*4.6,mu2=1.1*4.6,mu3=1.2*4.6,mu4=1.3*4.6,v11=1,v12=1,v13=1,v14=1))
alpha=mean(p);alpha
p<-replicate(5,fun4(n1=40,n2=60,n3=50,n4=70,mu1=1*4.8,mu2=1.1*4.8,mu3=1.2*4.8,mu4=1.3*4.8,v11=1,v12=1,v13=1,v14=1))
alpha=mean(p);alpha
p<-replicate(5,fun4(n1=40,n2=60,n3=50,n4=70,mu1=1*5,mu2=1.1*5,mu3=1.2*5,mu4=1.3*5,v11=1,v12=1,v13=1,v14=1))
alpha=mean(p);alpha

##Power values of Table 2.10

p<-replicate(5,fun4(n1=20,n2=25,n3=30,n4=20,mu1=1,mu2=1,mu3=1,mu4=1,v11=1,v12=2,v13=3,v14=4))
alpha=mean(p);alpha
p<-replicate(5,fun4(n1=20,n2=25,n3=30,n4=20,mu1=1,mu2=1.1,mu3=1.2,mu4=1.3,v11=1,v12=2,v13=3,v14=4))
alpha=mean(p);alpha
p<-replicate(5,fun4(n1=20,n2=25,n3=30,n4=20,mu1=1*1.2,mu2=1.1*1.2,mu3=1.2*1.2,mu4=1.3*1.2,v11=1,v12=2,v13=3,v14=4))
alpha=mean(p);alpha
p<-replicate(5,fun4(n1=20,n2=25,n3=30,n4=20,mu1=1*1.4,mu2=1.1*1.4,mu3=1.2*1.4,mu4=1.3*1.4,v11=1,v12=2,v13=3,v14=4))
alpha=mean(p);alpha
p<-replicate(5,fun4(n1=20,n2=25,n3=30,n4=20,mu1=1*1.6,mu2=1.1*1.6,mu3=1.2*1.6,mu4=1.3*1.6,v11=1,v12=2,v13=3,v14=4))
alpha=mean(p);alpha
p<-replicate(5,fun4(n1=20,n2=25,n3=30,n4=20,mu1=1*1.8,mu2=1.1*1.8,mu3=1.2*1.8,mu4=1.3*1.8,v11=1,v12=2,v13=3,v14=4))
alpha=mean(p);alpha
p<-replicate(5,fun4(n1=20,n2=25,n3=30,n4=20,mu1=1*2,mu2=1.1*2,mu3=1.2*2,mu4=1.3*2,v11=1,v12=2,v13=3,v14=4))
alpha=mean(p);alpha
p<-replicate(5,fun4(n1=20,n2=25,n3=30,n4=20,mu1=1*2.2,mu2=1.1*2.2,mu3=1.2*2.2,mu4=1.3*2.2,v11=1,v12=2,v13=3,v14=4))
alpha=mean(p);alpha
p<-replicate(5,fun4(n1=20,n2=25,n3=30,n4=20,mu1=1*2.4,mu2=1.1*2.4,mu3=1.2*2.4,mu4=1.3*2.4,v11=1,v12=2,v13=3,v14=4))
alpha=mean(p);alpha
p<-replicate(5,fun4(n1=20,n2=25,n3=30,n4=20,mu1=1*2.6,mu2=1.1*2.6,mu3=1.2*2.6,mu4=1.3*2.6,v11=1,v12=2,v13=3,v14=4))
alpha=mean(p);alpha
p<-replicate(5,fun4(n1=20,n2=25,n3=30,n4=20,mu1=1*2.8,mu2=1.1*2.8,mu3=1.2*2.8,mu4=1.3*2.8,v11=1,v12=2,v13=3,v14=4))
alpha=mean(p);alpha
p<-replicate(5,fun4(n1=20,n2=25,n3=30,n4=20,mu1=1*3,mu2=1.1*3,mu3=1.2*3,mu4=1.3*3,v11=1,v12=2,v13=3,v14=4))
alpha=mean(p);alpha
p<-replicate(5,fun4(n1=20,n2=25,n3=30,n4=20,mu1=1*3.2,mu2=1.1*3.2,mu3=1.2*3.2,mu4=1.3*3.2,v11=1,v12=2,v13=3,v14=4))
alpha=mean(p);alpha
p<-replicate(5,fun4(n1=20,n2=25,n3=30,n4=20,mu1=1*3.4,mu2=1.1*3.4,mu3=1.2*3.4,mu4=1.3*3.4,v11=1,v12=2,v13=3,v14=4))
alpha=mean(p);alpha
p<-replicate(5,fun4(n1=20,n2=25,n3=30,n4=20,mu1=1*3.6,mu2=1.1*3.6,mu3=1.2*3.6,mu4=1.3*3.6,v11=1,v12=2,v13=3,v14=4))
alpha=mean(p);alpha
p<-replicate(5,fun4(n1=20,n2=25,n3=30,n4=20,mu1=1*3.8,mu2=1.1*3.8,mu3=1.2*3.8,mu4=1.3*3.8,v11=1,v12=2,v13=3,v14=4))
alpha=mean(p);alpha
p<-replicate(5,fun4(n1=20,n2=25,n3=30,n4=20,mu1=1*4,mu2=1.1*4,mu3=1.2*4,mu4=1.3*4,v11=1,v12=2,v13=3,v14=4))
alpha=mean(p);alpha
p<-replicate(5,fun4(n1=20,n2=25,n3=30,n4=20,mu1=1*4.2,mu2=1.1*4.2,mu3=1.2*4.2,mu4=1.3*4.2,v11=1,v12=2,v13=3,v14=4))
alpha=mean(p);alpha
p<-replicate(5,fun4(n1=20,n2=25,n3=30,n4=20,mu1=1*4.4,mu2=1.1*4.4,mu3=1.2*4.4,mu4=1.3*4.4,v11=1,v12=2,v13=3,v14=4))
alpha=mean(p);alpha
p<-replicate(5,fun4(n1=20,n2=25,n3=30,n4=20,mu1=1*4.6,mu2=1.1*4.6,mu3=1.2*4.6,mu4=1.3*4.6,v11=1,v12=2,v13=3,v14=4))
alpha=mean(p);alpha
p<-replicate(5,fun4(n1=20,n2=25,n3=30,n4=20,mu1=1*4.8,mu2=1.1*4.8,mu3=1.2*4.8,mu4=1.3*4.8,v11=1,v12=2,v13=3,v14=4))
alpha=mean(p);alpha
p<-replicate(5,fun4(n1=20,n2=25,n3=30,n4=20,mu1=1*5,mu2=1.1*5,mu3=1.2*5,mu4=1.3*5,v11=1,v12=2,v13=3,v14=4))
alpha=mean(p);alpha

p<-replicate(5,fun4(n1=40,n2=60,n3=50,n4=70,mu1=1,mu2=1,mu3=1,mu4=1,v11=1,v12=2,v13=3,v14=4))
alpha=mean(p);alpha
p<-replicate(5,fun4(n1=40,n2=60,n3=50,n4=70,mu1=1,mu2=1.1,mu3=1.2,mu4=1.3,v11=1,v12=2,v13=3,v14=4))
alpha=mean(p);alpha
p<-replicate(5,fun4(n1=40,n2=60,n3=50,n4=70,mu1=1*1.2,mu2=1.1*1.2,mu3=1.2*1.2,mu4=1.3*1.2,v11=1,v12=2,v13=3,v14=4))
alpha=mean(p);alpha
p<-replicate(5,fun4(n1=40,n2=60,n3=50,n4=70,mu1=1*1.4,mu2=1.1*1.4,mu3=1.2*1.4,mu4=1.3*1.4,v11=1,v12=2,v13=3,v14=4))
alpha=mean(p);alpha
p<-replicate(5,fun4(n1=40,n2=60,n3=50,n4=70,mu1=1*1.6,mu2=1.1*1.6,mu3=1.2*1.6,mu4=1.3*1.6,v11=1,v12=2,v13=3,v14=4))
alpha=mean(p);alpha
p<-replicate(5,fun4(n1=40,n2=60,n3=50,n4=70,mu1=1*1.8,mu2=1.1*1.8,mu3=1.2*1.8,mu4=1.3*1.8,v11=1,v12=2,v13=3,v14=4))
alpha=mean(p);alpha
p<-replicate(5,fun4(n1=40,n2=60,n3=50,n4=70,mu1=1*2,mu2=1.1*2,mu3=1.2*2,mu4=1.3*2,v11=1,v12=2,v13=3,v14=4))
alpha=mean(p);alpha
p<-replicate(5,fun4(n1=40,n2=60,n3=50,n4=70,mu1=1*2.2,mu2=1.1*2.2,mu3=1.2*2.2,mu4=1.3*2.2,v11=1,v12=2,v13=3,v14=4))
alpha=mean(p);alpha
p<-replicate(5,fun4(n1=40,n2=60,n3=50,n4=70,mu1=1*2.4,mu2=1.1*2.4,mu3=1.2*2.4,mu4=1.3*2.4,v11=1,v12=2,v13=3,v14=4))
alpha=mean(p);alpha
p<-replicate(5,fun4(n1=40,n2=60,n3=50,n4=70,mu1=1*2.6,mu2=1.1*2.6,mu3=1.2*2.6,mu4=1.3*2.6,v11=1,v12=2,v13=3,v14=4))
alpha=mean(p);alpha
p<-replicate(5,fun4(n1=40,n2=60,n3=50,n4=70,mu1=1*2.8,mu2=1.1*2.8,mu3=1.2*2.8,mu4=1.3*2.8,v11=1,v12=2,v13=3,v14=4))
alpha=mean(p);alpha
p<-replicate(5,fun4(n1=40,n2=60,n3=50,n4=70,mu1=1*3,mu2=1.1*3,mu3=1.2*3,mu4=1.3*3,v11=1,v12=2,v13=3,v14=4))
alpha=mean(p);alpha
p<-replicate(5,fun4(n1=40,n2=60,n3=50,n4=70,mu1=1*3.2,mu2=1.1*3.2,mu3=1.2*3.2,mu4=1.3*3.2,v11=1,v12=2,v13=3,v14=4))
alpha=mean(p);alpha
p<-replicate(5,fun4(n1=40,n2=60,n3=50,n4=70,mu1=1*3.4,mu2=1.1*3.4,mu3=1.2*3.4,mu4=1.3*3.4,v11=1,v12=2,v13=3,v14=4))
alpha=mean(p);alpha
p<-replicate(5,fun4(n1=40,n2=60,n3=50,n4=70,mu1=1*3.6,mu2=1.1*3.6,mu3=1.2*3.6,mu4=1.3*3.6,v11=1,v12=2,v13=3,v14=4))
alpha=mean(p);alpha
p<-replicate(5,fun4(n1=40,n2=60,n3=50,n4=70,mu1=1*3.8,mu2=1.1*3.8,mu3=1.2*3.8,mu4=1.3*3.8,v11=1,v12=2,v13=3,v14=4))
alpha=mean(p);alpha
p<-replicate(5,fun4(n1=40,n2=60,n3=50,n4=70,mu1=1*4,mu2=1.1*4,mu3=1.2*4,mu4=1.3*4,v11=1,v12=2,v13=3,v14=4))
alpha=mean(p);alpha
p<-replicate(5,fun4(n1=40,n2=60,n3=50,n4=70,mu1=1*4.2,mu2=1.1*4.2,mu3=1.2*4.2,mu4=1.3*4.2,v11=1,v12=2,v13=3,v14=4))
alpha=mean(p);alpha
p<-replicate(5,fun4(n1=40,n2=60,n3=50,n4=70,mu1=1*4.4,mu2=1.1*4.4,mu3=1.2*4.4,mu4=1.3*4.4,v11=1,v12=2,v13=3,v14=4))
alpha=mean(p);alpha
p<-replicate(5,fun4(n1=40,n2=60,n3=50,n4=70,mu1=1*4.6,mu2=1.1*4.6,mu3=1.2*4.6,mu4=1.3*4.6,v11=1,v12=2,v13=3,v14=4))
alpha=mean(p);alpha
p<-replicate(5,fun4(n1=40,n2=60,n3=50,n4=70,mu1=1*4.8,mu2=1.1*4.8,mu3=1.2*4.8,mu4=1.3*4.8,v11=1,v12=2,v13=3,v14=4))
alpha=mean(p);alpha
p<-replicate(5,fun4(n1=40,n2=60,n3=50,n4=70,mu1=1*5,mu2=1.1*5,mu3=1.2*5,mu4=1.3*5,v11=1,v12=2,v13=3,v14=4))
alpha=mean(p);alpha


## Power values of Table 2.11


p<-replicate(5,fun4(n1=20,n2=25,n3=30,n4=20,mu1=1,mu2=1,mu3=1,mu4=1,v11=1,v12=1,v13=1,v14=1))
alpha=mean(p);alpha
p<-replicate(5,fun4(n1=20,n2=25,n3=30,n4=20,mu1=1*1.6,mu2=1.2*1.6,mu3=1.2*1.6,mu4=1.5*1.6,v11=1,v12=1,v13=1,v14=1))
alpha=mean(p);alpha
p<-replicate(5,fun4(n1=20,n2=25,n3=30,n4=20,mu1=1*1.4,mu2=1.05*1.4,mu3=1.15*1.4,mu4=1.3*1.4,v11=1,v12=1,v13=1,v14=1))
alpha=mean(p);alpha
p<-replicate(5,fun4(n1=20,n2=25,n3=30,n4=20,mu1=1*1.6,mu2=1.05*1.6,mu3=1.15*1.6,mu4=1.3*1.6,v11=1,v12=1,v13=1,v14=1))
alpha=mean(p);alpha
p<-replicate(5,fun4(n1=20,n2=25,n3=30,n4=20,mu1=1*1.8,mu2=1.05*1.8,mu3=1.15*1.8,mu4=1.3*1.8,v11=1,v12=1,v13=1,v14=1))
alpha=mean(p);alpha
p<-replicate(5,fun4(n1=20,n2=25,n3=30,n4=20,mu1=1*2,mu2=1.05*2,mu3=1.15*2,mu4=1.3*2,v11=1,v12=1,v13=1,v14=1))
alpha=mean(p);alpha
p<-replicate(5,fun4(n1=20,n2=25,n3=30,n4=20,mu1=1*2.2,mu2=1.05*2.2,mu3=1.15*2.2,mu4=1.3*2.2,v11=1,v12=1,v13=1,v14=1))
alpha=mean(p);alpha
p<-replicate(5,fun4(n1=20,n2=25,n3=30,n4=20,mu1=1*2.4,mu2=1.05*2.4,mu3=1.15*2.4,mu4=1.3*2.4,v11=1,v12=1,v13=1,v14=1))
alpha=mean(p);alpha
p<-replicate(5,fun4(n1=20,n2=25,n3=30,n4=20,mu1=1*2.6,mu2=1.05*2.6,mu3=1.15*2.6,mu4=1.3*2.6,v11=1,v12=1,v13=1,v14=1))
alpha=mean(p);alpha
p<-replicate(5,fun4(n1=20,n2=25,n3=30,n4=20,mu1=1*2.8,mu2=1.05*2.8,mu3=1.15*2.8,mu4=1.3*2.8,v11=1,v12=1,v13=1,v14=1))
alpha=mean(p);alpha
p<-replicate(5,fun4(n1=20,n2=25,n3=30,n4=20,mu1=1*3,mu2=1.05*3,mu3=1.15*3,mu4=1.3*3,v11=1,v12=1,v13=1,v14=1))
alpha=mean(p);alpha
p<-replicate(5,fun4(n1=20,n2=25,n3=30,n4=20,mu1=1*3.2,mu2=1.05*3.2,mu3=1.15*3.2,mu4=1.3*3.2,v11=1,v12=1,v13=1,v14=1))
alpha=mean(p);alpha
p<-replicate(5,fun4(n1=20,n2=25,n3=30,n4=20,mu1=1*3.4,mu2=1.05*3.4,mu3=1.15*3.4,mu4=1.3*3.4,v11=1,v12=1,v13=1,v14=1))
alpha=mean(p);alpha
p<-replicate(5,fun4(n1=20,n2=25,n3=30,n4=20,mu1=1*3.6,mu2=1.05*3.6,mu3=1.15*3.6,mu4=1.3*3.6,v11=1,v12=1,v13=1,v14=1))
alpha=mean(p);alpha
p<-replicate(5,fun4(n1=20,n2=25,n3=30,n4=20,mu1=1*3.8,mu2=1.05*3.8,mu3=1.15*3.8,mu4=1.3*3.8,v11=1,v12=1,v13=1,v14=1))
alpha=mean(p);alpha
p<-replicate(5,fun4(n1=20,n2=25,n3=30,n4=20,mu1=1*4,mu2=1.05*4,mu3=1.15*4,mu4=1.3*4,v11=1,v12=1,v13=1,v14=1))
alpha=mean(p);alpha
p<-replicate(5,fun4(n1=20,n2=25,n3=30,n4=20,mu1=1*4.2,mu2=1.05*4.2,mu3=1.15*4.2,mu4=1.3*4.2,v11=1,v12=1,v13=1,v14=1))
alpha=mean(p);alpha
p<-replicate(5,fun4(n1=20,n2=25,n3=30,n4=20,mu1=1*4.4,mu2=1.05*4.4,mu3=1.15*4.4,mu4=1.3*4.4,v11=1,v12=1,v13=1,v14=1))
alpha=mean(p);alpha
p<-replicate(5,fun4(n1=20,n2=25,n3=30,n4=20,mu1=1*4.6,mu2=1.05*4.6,mu3=1.15*4.6,mu4=1.3*4.6,v11=1,v12=1,v13=1,v14=1))
alpha=mean(p);alpha
p<-replicate(5,fun4(n1=20,n2=25,n3=30,n4=20,mu1=1*4.8,mu2=1.05*4.8,mu3=1.15*4.8,mu4=1.3*4.8,v11=1,v12=1,v13=1,v14=1))
alpha=mean(p);alpha
p<-replicate(5,fun4(n1=20,n2=25,n3=30,n4=20,mu1=1*5,mu2=1.05*5,mu3=1.15*5,mu4=1.3*5,v11=1,v12=1,v13=1,v14=1))
alpha=mean(p);alpha


p<-replicate(5,fun4(n1=20,n2=25,n3=30,n4=20,mu1=1,mu2=1,mu3=1,mu4=1,v11=1,v12=1,v13=3,v14=4))
alpha=mean(p);alpha
p<-replicate(5,fun4(n1=20,n2=25,n3=30,n4=20,mu1=1*1.6,mu2=1.2*1.6,mu3=1.2*1.6,mu4=1.5*1.6,v11=1,v12=1,v13=1,v14=1))
alpha=mean(p);alpha
p<-replicate(5,fun4(n1=20,n2=25,n3=30,n4=20,mu1=1*1.4,mu2=1.05*1.4,mu3=1.15*1.4,mu4=1.3*1.4,v11=1,v12=1,v13=3,v14=4))
alpha=mean(p);alpha
p<-replicate(5,fun4(n1=20,n2=25,n3=30,n4=20,mu1=1*1.6,mu2=1.05*1.6,mu3=1.15*1.6,mu4=1.3*1.6,v11=1,v12=1,v13=3,v14=4))
alpha=mean(p);alpha
p<-replicate(5,fun4(n1=20,n2=25,n3=30,n4=20,mu1=1*1.8,mu2=1.05*1.8,mu3=1.15*1.8,mu4=1.3*1.8,v11=1,v12=1,v13=3,v14=4))
alpha=mean(p);alpha
p<-replicate(5,fun4(n1=20,n2=25,n3=30,n4=20,mu1=1*2,mu2=1.05*2,mu3=1.15*2,mu4=1.3*2,v11=1,v12=1,v13=3,v14=4))
alpha=mean(p);alpha
p<-replicate(5,fun4(n1=20,n2=25,n3=30,n4=20,mu1=1*2.2,mu2=1.05*2.2,mu3=1.15*2.2,mu4=1.3*2.2,v11=1,v12=1,v13=3,v14=4))
alpha=mean(p);alpha
p<-replicate(5,fun4(n1=20,n2=25,n3=30,n4=20,mu1=1*2.4,mu2=1.05*2.4,mu3=1.15*2.4,mu4=1.3*2.4,v11=1,v12=1,v13=3,v14=4))
alpha=mean(p);alpha
p<-replicate(5,fun4(n1=20,n2=25,n3=30,n4=20,mu1=1*2.6,mu2=1.05*2.6,mu3=1.15*2.6,mu4=1.3*2.6,v11=1,v12=1,v13=3,v14=4))
alpha=mean(p);alpha
p<-replicate(5,fun4(n1=20,n2=25,n3=30,n4=20,mu1=1*2.8,mu2=1.05*2.8,mu3=1.15*2.8,mu4=1.3*2.8,v11=1,v12=1,v13=3,v14=4))
alpha=mean(p);alpha
p<-replicate(5,fun4(n1=20,n2=25,n3=30,n4=20,mu1=1*3,mu2=1.05*3,mu3=1.15*3,mu4=1.3*3,v11=1,v12=1,v13=3,v14=4))
alpha=mean(p);alpha
p<-replicate(5,fun4(n1=20,n2=25,n3=30,n4=20,mu1=1*3.2,mu2=1.05*3.2,mu3=1.15*3.2,mu4=1.3*3.2,v11=1,v12=1,v13=3,v14=4))
alpha=mean(p);alpha
p<-replicate(5,fun4(n1=20,n2=25,n3=30,n4=20,mu1=1*3.4,mu2=1.05*3.4,mu3=1.15*3.4,mu4=1.3*3.4,v11=1,v12=1,v13=3,v14=4))
alpha=mean(p);alpha
p<-replicate(5,fun4(n1=20,n2=25,n3=30,n4=20,mu1=1*3.6,mu2=1.05*3.6,mu3=1.15*3.6,mu4=1.3*3.6,v11=1,v12=1,v13=3,v14=4))
alpha=mean(p);alpha
p<-replicate(5,fun4(n1=20,n2=25,n3=30,n4=20,mu1=1*3.8,mu2=1.05*3.8,mu3=1.15*3.8,mu4=1.3*3.8,v11=1,v12=1,v13=3,v14=4))
alpha=mean(p);alpha
p<-replicate(5,fun4(n1=20,n2=25,n3=30,n4=20,mu1=1*4,mu2=1.05*4,mu3=1.15*4,mu4=1.3*4,v11=1,v12=1,v13=3,v14=4))
alpha=mean(p);alpha
p<-replicate(5,fun4(n1=20,n2=25,n3=30,n4=20,mu1=1*4.2,mu2=1.05*4.2,mu3=1.15*4.2,mu4=1.3*4.2,v11=1,v12=1,v13=3,v14=4))
alpha=mean(p);alpha
p<-replicate(5,fun4(n1=20,n2=25,n3=30,n4=20,mu1=1*4.4,mu2=1.05*4.4,mu3=1.15*4.4,mu4=1.3*4.4,v11=1,v12=1,v13=3,v14=4))
alpha=mean(p);alpha
p<-replicate(5,fun4(n1=20,n2=25,n3=30,n4=20,mu1=1*4.6,mu2=1.05*4.6,mu3=1.15*4.6,mu4=1.3*4.6,v11=1,v12=1,v13=3,v14=4))
alpha=mean(p);alpha
p<-replicate(5,fun4(n1=20,n2=25,n3=30,n4=20,mu1=1*4.8,mu2=1.05*4.8,mu3=1.15*4.8,mu4=1.3*4.8,v11=1,v12=1,v13=3,v14=4))
alpha=mean(p);alpha
p<-replicate(5,fun4(n1=20,n2=25,n3=30,n4=20,mu1=1*5,mu2=1.05*5,mu3=1.15*5,mu4=1.3*5,v11=1,v12=1,v13=3,v14=4))
alpha=mean(p);alpha

