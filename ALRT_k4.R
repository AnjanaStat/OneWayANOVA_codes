rm(list=ls())
#install.packages("Iso")
# package Iso is needed for the pava function
library(Iso)
fun1<-function(n1,n2,n3,n4,mu1,mu2,mu3,mu4,v11,v12,v13,v14)
{
  # generate sample 
  X1=rnorm(n1,mu1,sqrt(v11))
  X2=rnorm(n2,mu2,sqrt(v12))
  X3=rnorm(n3,mu3,sqrt(v13))
  X4=rnorm(n4,mu4,sqrt(v14))
  # follwing steps are for finding coefficients of the polynomial of mu
  m1<-mean(X1)
  m2<-mean(X2)
  m3<-mean(X3)
  m4<-mean(X4)
  p1<-sum(X1^2)/n1
  p2<-sum(X2^2)/n2
  p3<-sum(X3^2)/n3
  p4<-sum(X4^2)/n4
  f1<-p2*p3
  f2<--2*(m3*p2+m2*p3)
  f3<-p2+4*m2*m3+p3
  f4<--2*(m2+m3)
  g1<-p1*p3
  g2<--2*(m3*p1+m1*p3)
  g3<-p3+4*m1*m3+p1
  g4<--2*(m1+m3)
  h1<-p1*p2
  h2<--2*(m1*p2+m2*p1)
  h3<-p1+p2+4*m1*m2
  h4<--2*(m1+m2)
  a1<-f1*p4
  a2<-f2*p4-2*f1*m4
  a3<-f1-2*f2*m4+f3*p4
  a4<-f2-2*f3*m4+f4*p4
  a5<-f3-2*f4*m4+p4
  a6<-f4-2*m4
  b1<-g1*p4
  b2<-g2*p4-2*g1*m4
  b3<-g1-2*g2*m4+g3*p4
  b4<-g2-2*g3*m4+g4*p4
  b5<-g3-2*g4*m4+p4
  b6<-g4-2*m4
  c1<-h1*p4
  c2<-h2*p4-2*m4*h1
  c3<-h1-2*h2*m4+h3*p4
  c4<-h2-2*h3*m4+h4*p4
  c5<-h3-2*h4*m4+p4
  c6<-h4-2*m4
  d1<-h1*p3
  d2<-h2*p3-2*m3*h1
  d3<-h1-2*h2*m3+h3*p3
  d4<-h2-2*h3*m3+h4*p3
  d5<-h3-2*h4*m3+p3
  d6<-h4-2*m3
  # coefficients of the polynomial
  z1<-n1*m1*a1+n2*m2*b1+n3*m3*c1+n4*m4*d1
  z2<-n1*(m1*a2-a1)+n2*(m2*b2-b1)+n3*(m3*c2-c1)+n4*(m4*d2-d1)
  z3<-n1*(m1*a3-a2)+n2*(m2*b3-b2)+n3*(m3*c3-c2)+n4*(m4*d3-d2)
  z4<-n1*(m1*a4-a3)+n2*(m2*b4-b3)+n3*(m3*c4-c3)+n4*(m4*d4-d3)
  z5<-n1*(m1*a5-a4)+n2*(m2*b5-b4)+n3*(m3*c5-c4)+n4*(m4*d5-d4)
  z6<-n1*(m1*a6-a5)+n2*(m2*b6-b5)+n3*(m3*c6-c5)+n4*(m4*d6-d5)
  z7<-n1*(m1-a6)+n2*(m2-b6)+n3*(m3-c6)+n4*(m4-d6)
  z8<--(n1+n2+n3+n4)
  # find the roots of the polynomial
  root<-polyroot(c(z1,z2,z3,z4,z5,z6,z7,z8))
  # find the real roots of the polynomial
  Rert=Re(root)[abs(Im(root))<=0.000001]
  # find the number of real roots
  n=length(Rert)
  # find the values of sigma_1^2 corresponding to the real values of mu 
  m7=NULL;m8=NULL;m9=NULL;m10=NULL;m11=NULL
  for(i in 1:n)
  {

      m7[i]=(sum((X1-Rert[i])^2))/n1
      m8[i]=(sum((X2-Rert[i])^2))/n2
      m9[i]=(sum((X3-Rert[i])^2))/n3
      m10[i]=(sum((X4-Rert[i])^2))/n4
      m11[i]=(m7[i]^(n1/2))*(m8[i]^(n2/2))*(m9[i]^(n3/2))*(m10[i]^(n4/2))
      
  }
  # find the MLE of sigma_1^2
  m111=min(m11)  
  # bellow  steps are for finding the MLEs over the alternative space
  x<-c(m1, m2, m3, m4)
  s1<-sum((X1-m1)^2)/n1
  s2<-sum((X2-m2)^2)/n2
  s3<-sum((X3-m3)^2)/n3
  s4<-sum((X4-m4)^2)/n4
  s<-c(s1, s2, s3, s4)
  x0=x
  s0<-s
  x8<-x
  repeat
  {
    w1=c(n1/s0[1],n2/s0[2],n3/s0[3],n4/s0[4])
    x6=pava(x0, w1)
    s11<-sum((X1-x6[1])^2)/n1
    s12<-sum((X2-x6[2])^2)/n2
    s13<-sum((X3-x6[3])^2)/n3
    s14<-sum((X4-x6[4])^2)/n4
    v=c(s11, s12, s13, s14)
    s0<-v
    x10=abs(x8-x6)
    p=max(x10[1], x10[2], x10[3], x10[4], na.rm = FALSE)
    if(p<=0.0000001)
    {
      break
    }
    x8<-x6
  }
  neu<-(s0[1])^(n1/2)*(s0[2])^(n2/2)*(s0[3])^(n3/2)*(s0[4])^(n4/2)
  dno<-m111
  # find the likelihood ratio value
  LR=neu/dno
  # find the value of the statistic=-2log(lambda)
  value1=-2*log(LR)
  return(value1)
}
fun4<-function(n1,n2,n3,n4,mu1,mu2,mu3,mu4,v11,v12,v13,v14)
{
  # find 10000 values of statistic
  out<-replicate(10000,fun1(n1,n2,n3,n4,mu1,mu2,mu3,mu4,v11,v12,v13,v14))
  # find chi-square value with alpha=0.05, df=2
  q=qchisq(0.05,df=2,lower.tail=FALSE)
  # count if the statistic value is greater than critical value among 10000 values
  a=0
  for(i in 1:10000)
  {
    if(out[i]>q)
      a=a+1
  }
  p=a/10000
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
p<-replicate(5,fun4(n1=20,n2=25,n3=30,n4=20,mu1=1*1.6,mu2=1.2*1.6,mu3=1.2*1.6,mu4=1.5*1.6,v11=1,v12=1,v13=3,v14=4))
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

