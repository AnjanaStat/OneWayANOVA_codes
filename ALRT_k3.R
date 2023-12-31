rm(list=ls())
# package Iso is needed for the pava function
install.packages("Iso")
# smoothest package for generating sample from Laplace distribution
install.packages("smoothmest")
library(Iso)
library(smoothmest)
fun1<-function(mu1,mu2,mu3,n1,n2,n3,v11,v12,v13)
{
  # For generating samples from Laplace distribution
  #g1=rdoublex(n1,mu1,sqrt(v11)/sqrt(2))
  #g2=rdoublex(n2,mu2,sqrt(v12)/sqrt(2))
  #g3=rdoublex(n3,mu3,sqrt(v13)/sqrt(2))
  ##For generating samples from t distribution
  #Y1=rt(n1,3,ncp=0)/sqrt(3);Y2=rt(n2,3,ncp=0)/sqrt(3)
  #Y3=rt(n3,3,ncp=0)/sqrt(3)
  #g1=sqrt(v11)*Y1+mu1;g2=sqrt(v12)*Y2+mu2;g3=sqrt(v13)*Y3+mu3
  #Y1=rexp(n1,1);Y2=rexp(n2,1);Y3=rexp(n3,1)
  #g1=mu1+sqrt(v11)*(Y1-1);g2=mu2+sqrt(v12)*(Y2-1);g3=mu3+sqrt(v13)*(Y3-1)
  ##For generating samples from lognormal distribution
  #Y1=(rlnorm(n1,0,1)-exp(.5))/sqrt(exp(2)-exp(1));Y2=(rlnorm(n2,0,1)-exp(.5))/sqrt(exp(2)-exp(1))
  #Y3=(rlnorm(n3,0,1)-exp(.5))/sqrt(exp(2)-exp(1))
  #g1=mu1+sqrt(v11)*Y1;g2=mu2+sqrt(v12)*Y2;g3=mu3+sqrt(v13)*Y3
  ##For generating samples from Weibull distribution
  #Y1=(rweibull(n1,2,1)-sqrt(pi)/2)/sqrt(1-pi/4);Y2=(rweibull(n2,2,1)-sqrt(pi)/2)/sqrt(1-pi/4)
  #Y3=(rweibull(n3,2,1)-sqrt(pi)/2)/sqrt(1-pi/4)
  #g1=mu1+sqrt(v11)*Y1;g2=mu2+sqrt(v12)*Y2;g3=mu3+sqrt(v13)*Y3
  ##For generating samples from normal distributions
  g1=rnorm(n1,mu1,sqrt(v11));g2=rnorm(n2,mu2,sqrt(v12));g3=rnorm(n3,mu3,sqrt(v13))
  # following steps are for finding coefficients of the polynomial of mu
  m1<-mean(g1);m2<-mean(g2);m3<-mean(g3)
  p1<-sum(g1^2)/n1;p2<-sum(g2^2)/n2;p3<-sum(g3^2)/n3
  #coefficients of the polynomial
  a1<-p2*p3;a2<--2*(m3*p2+m2*p3);a3<-p2+4*m2*m3+p3;a4<--2*(m2+m3)
  b1<-p1*p3;b2<--2*(m3*p1+m1*p3);b3<-p3+4*m1*m3+p1;b4<--2*(m1+m3)
  c1<-p1*p2;c2<--2*(m1*p2+m2*p1);c3<-p1+p2+4*m1*m2;c4<--2*(m1+m2)
  z1<-n1*m1*a1+n2*m2*b1+n3*m3*c1
  z2<-n1*(m1*a2-a1)+n2*(m2*b2-b1)+n3*(m3*c2-c1)
  z3<-n1*(m1*a3-a2)+n2*(m2*b3-b2)+n3*(m3*c3-c2)
  z4<-n1*(m1*a4-a3)+n2*(m2*b4-b3)+n3*(m3*c4-c3)
  z5<-n1*(m1-a4)+n2*(m2-b4)+n3*(m3-c4)
  z6<--(n1+n2+n3)
  # find the roots of the polynomial
  root<-polyroot(c(z1,z2,z3,z4,z5,z6))
  #find the real roots of the polynomial
  Rert=Re(root)[abs(Im(root))<=0.000001]
  # find the number of real roots
  n=length(Rert)
  # find the values of sigma_1^2 corresponding to the real values of mu
  m6=NULL;m7=NULL;m8=NULL;m9=NULL
  for(i in 1:n)
  {
    m6[i]=(sum((g1-Rert[i])^2))/n1
    m7[i]=(sum((g2-Rert[i])^2))/n2
    m8[i]=(sum((g3-Rert[i])^2))/n3
    m9[i]=(m6[i]^(n1/2))*(m7[i]^(n2/2))*(m8[i]^(n3/2))
    
  }
  # find MLE of sigma_1^2
  m91=min(m9)
  # following steps are for finding MLes of mu and sigma_i^2 under the alternative space
  x<-c(m1, m2, m3)
  s1<-sum((g1-m1)^2)/n1;s2<-sum((g2-m2)^2)/n2;s3<-sum((g3-m3)^2)/n3;s<-c(s1, s2, s3)
  x0=x;s0=s;x8<-x0
  repeat
  {
    w1=c(n1/s0[1],n2/s0[2],n3/s0[3])
    x6=pava(x0, w1)
    s11<-sum((g1-x6[1])^2)/n1
    s12<-sum((g2-x6[2])^2)/n2
    s13<-sum((g3-x6[3])^2)/n3
    v=c(s11, s12, s13)
    s0<-v
    x10=abs(x8-x6)
    p=max(x10[1], x10[2], x10[3], na.rm = FALSE)
    if(p<=0.0000001)
    {
      break
    }
    x8<-x6
  }
  neu<-(s0[1])^(n1/2)*(s0[2])^(n2/2)*(s0[3])^(n3/2)
  dno<-m91
  # find the likelihood ratio value
  LR=neu/dno
  value=-2*log(LR)
  #return  the value of -2log(lambda) from fun1
  return(value)
}

fun4<-function(mu1,mu2,mu3,n1,n2,n3,v11,v12,v13)
{
  # find 10000 values of -2log(lambda) from fun1
  X<-replicate(10000,fun1(mu1,mu2,mu3,n1,n2,n3,v11,v12,v13))
  # find the critical value
  q=qchisq(0.05,df=1,lower.tail=FALSE)
  # count if -2log(lambda) is greater than the critical value
  a=0
  for(i in 1:10000)
  {
    if(X[i]>q)
      a=a+1
  }
  p=a/10000
  p
  return(p)
}

##Power values of Table 2.6

p<-replicate(5,fun4(1,1,1,20,30,25,1,1,1))
p;alpha<-mean(p);alpha
p<-replicate(5,fun4(1,1.1,1.2,20,30,25,1,1,1))
p;alpha<-mean(p);alpha
p<-replicate(5,fun4(1.2,1.1*1.2,1.2*1.2,20,30,25,1,1,1))
p;alpha<-mean(p);alpha
p<-replicate(5,fun4(1.4,1.1*1.4,1.2*1.4,20,30,25,1,1,1))
p;alpha<-mean(p);alpha
p<-replicate(5,fun4(1.6,1.1*1.6,1.2*1.6,20,30,25,1,1,1))
p;alpha<-mean(p);alpha
p<-replicate(5,fun4(1.8,1.1*1.8,1.2*1.8,20,30,25,1,1,1))
p;alpha<-mean(p);alpha
p<-replicate(5,fun4(2,1.1*2,1.2*2,20,30,25,1,1,1))
p;alpha<-mean(p);alpha
p<-replicate(5,fun4(2.2,1.1*2.2,1.2*2.2,20,30,25,1,1,1))
p;alpha<-mean(p);alpha
p<-replicate(5,fun4(2.4,1.1*2.4,1.2*2.4,20,30,25,1,1,1))
p;alpha<-mean(p);alpha
p<-replicate(5,fun4(2.6,1.1*2.6,1.2*2.6,20,30,25,1,1,1))
p;alpha<-mean(p);alpha
p<-replicate(5,fun4(2.8,1.1*2.8,1.2*2.8,20,30,25,1,1,1))
p;alpha<-mean(p);alpha
p<-replicate(5,fun4(3,1.1*3,1.2*3,20,30,25,1,1,1))
p;alpha<-mean(p);alpha
p<-replicate(5,fun4(3.2,1.1*3.2,1.2*3.2,20,30,25,1,1,1))
p;alpha<-mean(p);alpha
p<-replicate(5,fun4(3.4,1.1*3.4,1.2*3.4,20,30,25,1,1,1))
p;alpha<-mean(p);alpha
p<-replicate(5,fun4(3.6,1.1*3.6,1.2*3.6,20,30,25,1,1,1))
p;alpha<-mean(p);alpha
p<-replicate(5,fun4(3.8,1.1*3.8,1.2*3.8,20,30,25,1,1,1))
p;alpha<-mean(p);alpha
p<-replicate(5,fun4(4,1.1*4,1.2*4,20,30,25,1,1,1))
p;alpha<-mean(p);alpha
p<-replicate(5,fun4(4.2,1.1*4.2,1.2*4.2,20,30,25,1,1,1))
p;alpha<-mean(p);alpha
p<-replicate(5,fun4(4.4,1.1*4.4,1.2*4.4,20,30,25,1,1,1))
p;alpha<-mean(p);alpha
p<-replicate(5,fun4(4.6,1.1*4.6,1.2*4.6,20,30,25,1,1,1))
p;alpha<-mean(p);alpha
p<-replicate(5,fun4(4.8,1.1*4.8,1.2*4.8,20,30,25,1,1,1))
p;alpha<-mean(p);alpha
p<-replicate(5,fun4(5,1.1*5,1.2*5,20,30,25,1,1,1))
p;alpha<-mean(p);alpha


p<-replicate(5,fun4(1,1,1,60,70,50,1,1,1))
p;alpha<-mean(p);alpha
p<-replicate(5,fun4(1,1.1,1.2,60,70,50,1,1,1))
p;alpha<-mean(p);alpha
p<-replicate(5,fun4(1.2,1.1*1.2,1.2*1.2,60,70,50,1,1,1))
p;alpha<-mean(p);alpha
p<-replicate(5,fun4(1.4,1.1*1.4,1.2*1.4,60,70,50,1,1,1))
p;alpha<-mean(p);alpha
p<-replicate(5,fun4(1.6,1.1*1.6,1.2*1.6,60,70,50,1,1,1))
p;alpha<-mean(p);alpha
p<-replicate(5,fun4(1.8,1.1*1.8,1.2*1.8,60,70,50,1,1,1))
p;alpha<-mean(p);alpha
p<-replicate(5,fun4(2,1.1*2,1.2*2,60,70,50,1,1,1))
p;alpha<-mean(p);alpha
p<-replicate(5,fun4(2.2,1.1*2.2,1.2*2.2,60,70,50,1,1,1))
p;alpha<-mean(p);alpha
p<-replicate(5,fun4(2.4,1.1*2.4,1.2*2.4,60,70,50,1,1,1))
p;alpha<-mean(p);alpha
p<-replicate(5,fun4(2.6,1.1*2.6,1.2*2.6,60,70,50,1,1,1))
p;alpha<-mean(p);alpha
p<-replicate(5,fun4(2.8,1.1*2.8,1.2*2.8,60,70,50,1,1,1))
p;alpha<-mean(p);alpha
p<-replicate(5,fun4(3,1.1*3,1.2*3,60,70,50,1,1,1))
p;alpha<-mean(p);alpha
p<-replicate(5,fun4(3.2,1.1*3.2,1.2*3.2,60,70,50,1,1,1))
p;alpha<-mean(p);alpha
p<-replicate(5,fun4(3.4,1.1*3.4,1.2*3.4,60,70,50,1,1,1))
p;alpha<-mean(p);alpha
p<-replicate(5,fun4(3.6,1.1*3.6,1.2*3.6,60,70,50,1,1,1))
p;alpha<-mean(p);alpha
p<-replicate(5,fun4(3.8,1.1*3.8,1.2*3.8,60,70,50,1,1,1))
p;alpha<-mean(p);alpha
p<-replicate(5,fun4(4,1.1*4,1.2*4,60,70,50,1,1,1))
p;alpha<-mean(p);alpha
p<-replicate(5,fun4(4.2,1.1*4.2,1.2*4.2,60,70,50,1,1,1))
p;alpha<-mean(p);alpha
p<-replicate(5,fun4(4.4,1.1*4.4,1.2*4.4,60,70,50,1,1,1))
p;alpha<-mean(p);alpha
p<-replicate(5,fun4(4.6,1.1*4.6,1.2*4.6,60,70,50,1,1,1))
p;alpha<-mean(p);alpha
p<-replicate(5,fun4(4.8,1.1*4.8,1.2*4.8,60,70,50,1,1,1))
p;alpha<-mean(p);alpha
p<-replicate(5,fun4(5,1.1*5,1.2*5,60,70,50,1,1,1))
p;alpha<-mean(p);alpha

###Power values of Table 2.7

p<-replicate(5,fun4(1,1,1,20,30,25,1,2,3))
p;alpha<-mean(p);alpha
p<-replicate(5,fun4(1,1.1,1.2,20,30,25,1,2,3))
p;alpha<-mean(p);alpha
p<-replicate(5,fun4(1.2,1.1*1.2,1.2*1.2,20,30,25,1,2,3))
p;alpha<-mean(p);alpha
p<-replicate(5,fun4(1.4,1.1*1.4,1.2*1.4,20,30,25,1,2,3))
p;alpha<-mean(p);alpha
p<-replicate(5,fun4(1.6,1.1*1.6,1.2*1.6,20,30,25,1,2,3))
p;alpha<-mean(p);alpha
p<-replicate(5,fun4(1.8,1.1*1.8,1.2*1.8,20,30,25,1,2,3))
p;alpha<-mean(p);alpha
p<-replicate(5,fun4(2,1.1*2,1.2*2,20,30,25,1,2,3))
p;alpha<-mean(p);alpha
p<-replicate(5,fun4(2.2,1.1*2.2,1.2*2.2,20,30,25,1,2,3))
p;alpha<-mean(p);alpha
p<-replicate(5,fun4(2.4,1.1*2.4,1.2*2.4,20,30,25,1,2,3))
p;alpha<-mean(p);alpha
p<-replicate(5,fun4(2.6,1.1*2.6,1.2*2.6,20,30,25,1,2,3))
p;alpha<-mean(p);alpha
p<-replicate(5,fun4(2.8,1.1*2.8,1.2*2.8,20,30,25,1,2,3))
p;alpha<-mean(p);alpha
p<-replicate(5,fun4(3,1.1*3,1.2*3,20,30,25,1,2,3))
p;alpha<-mean(p);alpha
p<-replicate(5,fun4(3.2,1.1*3.2,1.2*3.2,20,30,25,1,2,3))
p;alpha<-mean(p);alpha
p<-replicate(5,fun4(3.4,1.1*3.4,1.2*3.4,20,30,25,1,2,3))
p;alpha<-mean(p);alpha
p<-replicate(5,fun4(3.6,1.1*3.6,1.2*3.6,20,30,25,1,2,3))
p;alpha<-mean(p);alpha
p<-replicate(5,fun4(3.8,1.1*3.8,1.2*3.8,20,30,25,1,2,3))
p;alpha<-mean(p);alpha
p<-replicate(5,fun4(4,1.1*4,1.2*4,20,30,25,1,2,3))
p;alpha<-mean(p);alpha
p<-replicate(5,fun4(4.2,1.1*4.2,1.2*4.2,20,30,25,1,2,3))
p;alpha<-mean(p);alpha
p<-replicate(5,fun4(4.4,1.1*4.4,1.2*4.4,20,30,25,1,2,3))
p;alpha<-mean(p);alpha
p<-replicate(5,fun4(4.6,1.1*4.6,1.2*4.6,20,30,25,1,2,3))
p;alpha<-mean(p);alpha
p<-replicate(5,fun4(4.8,1.1*4.8,1.2*4.8,20,30,25,1,2,3))
p;alpha<-mean(p);alpha
p<-replicate(5,fun4(5,1.1*5,1.2*5,20,30,25,1,2,3))
p;alpha<-mean(p);alpha


p<-replicate(5,fun4(1,1,1,60,70,50,1,2,3))
p;alpha<-mean(p);alpha
p<-replicate(5,fun4(1,1.1,1.2,60,70,50,1,2,3))
p;alpha<-mean(p);alpha
p<-replicate(5,fun4(1.2,1.1*1.2,1.2*1.2,60,70,50,1,2,3))
p;alpha<-mean(p);alpha
p<-replicate(5,fun4(1.4,1.1*1.4,1.2*1.4,60,70,50,1,2,3))
p;alpha<-mean(p);alpha
p<-replicate(5,fun4(1.6,1.1*1.6,1.2*1.6,60,70,50,1,2,3))
p;alpha<-mean(p);alpha
p<-replicate(5,fun4(1.8,1.1*1.8,1.2*1.8,60,70,50,1,2,3))
p;alpha<-mean(p);alpha
p<-replicate(5,fun4(2,1.1*2,1.2*2,60,70,50,1,2,3))
p;alpha<-mean(p);alpha
p<-replicate(5,fun4(2.2,1.1*2.2,1.2*2.2,60,70,50,1,2,3))
p;alpha<-mean(p);alpha
p<-replicate(5,fun4(2.4,1.1*2.4,1.2*2.4,60,70,50,1,2,3))
p;alpha<-mean(p);alpha
p<-replicate(5,fun4(2.6,1.1*2.6,1.2*2.6,60,70,50,1,2,3))
p;alpha<-mean(p);alpha
p<-replicate(5,fun4(2.8,1.1*2.8,1.2*2.8,60,70,50,1,2,3))
p;alpha<-mean(p);alpha
p<-replicate(5,fun4(3,1.1*3,1.2*3,60,70,50,1,2,3))
p;alpha<-mean(p);alpha
p<-replicate(5,fun4(3.2,1.1*3.2,1.2*3.2,60,70,50,1,2,3))
p;alpha<-mean(p);alpha
p<-replicate(5,fun4(3.4,1.1*3.4,1.2*3.4,60,70,50,1,2,3))
p;alpha<-mean(p);alpha
p<-replicate(5,fun4(3.6,1.1*3.6,1.2*3.6,60,70,50,1,2,3))
p;alpha<-mean(p);alpha
p<-replicate(5,fun4(3.8,1.1*3.8,1.2*3.8,60,70,50,1,2,3))
p;alpha<-mean(p);alpha
p<-replicate(5,fun4(4,1.1*4,1.2*4,60,70,50,1,2,3))
p;alpha<-mean(p);alpha
p<-replicate(5,fun4(4.2,1.1*4.2,1.2*4.2,60,70,50,1,2,3))
p;alpha<-mean(p);alpha
p<-replicate(5,fun4(4.4,1.1*4.4,1.2*4.4,60,70,50,1,2,3))
p;alpha<-mean(p);alpha
p<-replicate(5,fun4(4.6,1.1*4.6,1.2*4.6,60,70,50,1,2,3))
p;alpha<-mean(p);alpha
p<-replicate(5,fun4(4.8,1.1*4.8,1.2*4.8,60,70,50,1,2,3))
p;alpha<-mean(p);alpha
p<-replicate(5,fun4(5,1.1*5,1.2*5,60,70,50,1,2,3))
p;alpha<-mean(p);alpha


###Power values of Table 2.8

p<-replicate(5,fun4(0,0,0,20,30,25,1,1,1))
p;alpha<-mean(p);alpha
p<-replicate(5,fun4(1,1,1.5,20,30,25,1,1,1))
p;alpha<-mean(p);alpha
p<-replicate(5,fun4(1.2,1.2,1.5*1.2,20,30,25,1,1,1))
p;alpha<-mean(p);alpha
p<-replicate(5,fun4(1.4,1.4,1.5*1.4,20,30,25,1,1,1))
p;alpha<-mean(p);alpha
p<-replicate(5,fun4(1.6,1.6,1.5*1.6,20,30,25,1,1,1))
p;alpha<-mean(p);alpha
p<-replicate(5,fun4(1.8,1.8,1.5*1.8,20,30,25,1,1,1))
p;alpha<-mean(p);alpha
p<-replicate(5,fun4(2,2,1.5*2,20,30,25,1,1,1))
p;alpha<-mean(p);alpha
p<-replicate(5,fun4(2.2,2.2,1.5*2.2,20,30,25,1,1,1))
p;alpha<-mean(p);alpha
p<-replicate(5,fun4(2.4,2.4,1.5*2.4,20,30,25,1,1,1))
p;alpha<-mean(p);alpha
p<-replicate(5,fun4(2.6,2.6,1.5*2.6,20,30,25,1,1,1))
p;alpha<-mean(p);alpha
p<-replicate(5,fun4(2.8,2.8,1.5*2.8,20,30,25,1,1,1))
p;alpha<-mean(p);alpha
p<-replicate(5,fun4(3,3,1.5*3,20,30,25,1,1,1))
p;alpha<-mean(p);alpha
p<-replicate(5,fun4(3.2,3.2,1.5*3.2,20,30,25,1,1,1))
p;alpha<-mean(p);alpha
p<-replicate(5,fun4(3.4,3.4,1.5*3.4,20,30,25,1,1,1))
p;alpha<-mean(p);alpha
p<-replicate(5,fun4(3.6,3.6,1.5*3.6,20,30,25,1,1,1))
p;alpha<-mean(p);alpha
p<-replicate(5,fun4(3.8,3.8,1.5*3.8,20,30,25,1,1,1))
p;alpha<-mean(p);alpha
p<-replicate(5,fun4(4,4,1.5*4,20,30,25,1,1,1))
p;alpha<-mean(p);alpha
p<-replicate(5,fun4(4.2,4.2,1.5*4.2,20,30,25,1,1,1))
p;alpha<-mean(p);alpha
p<-replicate(5,fun4(4.4,4.4,1.5*4.4,20,30,25,1,1,1))
p;alpha<-mean(p);alpha
p<-replicate(5,fun4(4.6,4.6,1.5*4.6,20,30,25,1,1,1))
p;alpha<-mean(p);alpha
p<-replicate(5,fun4(4.8,4.8,1.5*4.8,20,30,25,1,1,1))
p;alpha<-mean(p);alpha
p<-replicate(5,fun4(5,5,1.5*5,20,30,25,1,1,1))
p;alpha<-mean(p);alpha


p<-replicate(5,fun4(0,0,0,20,30,25,1,2,3))
p;alpha<-mean(p);alpha
p<-replicate(5,fun4(1,1,1.5,20,30,25,1,2,3))
p;alpha<-mean(p);alpha
p<-replicate(5,fun4(1.2,1.2,1.5*1.2,20,30,25,1,2,3))
p;alpha<-mean(p);alpha
p<-replicate(5,fun4(1.4,1.4,1.5*1.4,20,30,25,1,2,3))
p;alpha<-mean(p);alpha
p<-replicate(5,fun4(1.6,1.6,1.5*1.6,20,30,25,1,2,3))
p;alpha<-mean(p);alpha
p<-replicate(5,fun4(1.8,1.8,1.5*1.8,20,30,25,1,2,3))
p;alpha<-mean(p);alpha
p<-replicate(5,fun4(2,2,1.5*2,20,30,25,1,2,3))
p;alpha<-mean(p);alpha
p<-replicate(5,fun4(2.2,2.2,1.5*2.2,20,30,25,1,2,3))
p;alpha<-mean(p);alpha
p<-replicate(5,fun4(2.4,2.4,1.5*2.4,20,30,25,1,2,3))
p;alpha<-mean(p);alpha
p<-replicate(5,fun4(2.6,2.6,1.5*2.6,20,30,25,1,2,3))
p;alpha<-mean(p);alpha
p<-replicate(5,fun4(2.8,2.8,1.5*2.8,20,30,25,1,2,3))
p;alpha<-mean(p);alpha
p<-replicate(5,fun4(3,3,1.5*3,20,30,25,1,2,3))
p;alpha<-mean(p);alpha
p<-replicate(5,fun4(3.2,3.2,1.5*3.2,20,30,25,1,2,3))
p;alpha<-mean(p);alpha
p<-replicate(5,fun4(3.4,3.4,1.5*3.4,20,30,25,1,2,3))
p;alpha<-mean(p);alpha
p<-replicate(5,fun4(3.6,3.6,1.5*3.6,20,30,25,1,2,3))
p;alpha<-mean(p);alpha
p<-replicate(5,fun4(3.8,3.8,1.5*3.8,20,30,25,1,2,3))
p;alpha<-mean(p);alpha
p<-replicate(5,fun4(4,4,1.5*4,20,30,25,1,2,3))
p;alpha<-mean(p);alpha
p<-replicate(5,fun4(4.2,4.2,1.5*4.2,20,30,25,1,2,3))
p;alpha<-mean(p);alpha
p<-replicate(5,fun4(4.4,4.4,1.5*4.4,20,30,25,1,2,3))
p;alpha<-mean(p);alpha
p<-replicate(5,fun4(4.6,4.6,1.5*4.6,20,30,25,1,2,3))
p;alpha<-mean(p);alpha
p<-replicate(5,fun4(4.8,4.8,1.5*4.8,20,30,25,1,2,3))
p;alpha<-mean(p);alpha
p<-replicate(5,fun4(5,5,1.5*5,20,30,25,1,2,3))
p;alpha<-mean(p);alpha