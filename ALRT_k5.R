#install.packages("Iso")
rm(list=ls())
# package Iso is needed for the pava function
library(Iso)
fun1<-function(n1,n2,n3,n4,n5,mu1,mu2,mu3,mu4,mu5,v11,v12,v13,v14,v15)
{
  # generating sample
  g1<-rnorm(n1,mu1,sqrt(v11))
  g2<-rnorm(n2,mu2,sqrt(v12))
  g3<-rnorm(n3,mu3,sqrt(v13))
  g4<-rnorm(n4,mu4,sqrt(v14))
  g5<-rnorm(n5,mu5,sqrt(v15))
  # below steps are for finding coefficients of the polynomial of mu(common mean)
  m1<-mean(g1)
  m2<-mean(g2)
  m3<-mean(g3)
  m4<-mean(g4)
  m5<-mean(g5)
  p1<-sum(g1^2)/n1
  p2<-sum(g2^2)/n2
  p3<-sum(g3^2)/n3
  p4<-sum(g4^2)/n4
  p5<-sum(g5^2)/n5
  a1<-p2*p3
  a2<--2*(m3*p2+p3*m2)
  a3<-p3+4*m2*m3+p2
  a4<--2*(m2+m3)
  b1<-p4*p5
  b2<--2*(m4*p5+p4*m5)
  b3<-p4+4*m4*m5+p5
  b4<--2*(m4+m5)
  c1<-p1*p3
  c2<--2*(m1*p3+m3*p1)
  c3<-p3+4*m1*m3+p1
  c4<--2*(m1+m3)
  d1<-p1*p2
  d2<--2*(m1*p2+m2*p1)
  d3<-p2+4*m1*m2+p1
  d4<--2*(m1+m2)
  f1<-p3*p5
  f2<--2*(m3*p5+m5*p3)
  f3<-p3+4*m3*m5+p5
  f4<--2*(m3+m5)
  h1<-p3*p4
  h2<--2*(m3*p4+m4*p3)
  h3<-p3+4*m3*m4+p4
  h4<--2*(m3+m4)
  q1<-a1*b1
  q2<-a2*b1+a1*b2
  q3<-a1*b3+a2*b2+a3*b1
  q4<-a1*b4+a2*b3+a3*b2+a4*b1
  q5<-a1+a2*b4+a3*b3+a4*b2+b1
  q6<-a2+a3*b4+a4*b3+b2
  q7<-a3+a4*b4+b3
  q8<-a4+b4
  r1<-c1*b1
  r2<-c2*b1+c1*b2
  r3<-c1*b3+c2*b2+c3*b1
  r4<-c1*b4+c2*b3+c3*b2+c4*b1
  r5<-c1+c2*b4+c3*b3+c4*b2+b1
  r6<-c2+c3*b4+c4*b3+b2
  r7<-c3+c4*b4+b3
  r8<-c4+b4
  s1<-d1*b1
  s2<-d2*b1+d1*b2
  s3<-d1*b3+d2*b2+d3*b1
  s4<-d1*b4+d2*b3+d3*b2+d4*b1
  s5<-d1+d2*b4+d3*b3+d4*b2+b1
  s6<-d2+d3*b4+d4*b3+b2
  s7<-d3+d4*b4+b3
  s8<-d4+b4
  t1<-d1*f1
  t2<-d2*f1+d1*f2
  t3<-d1*f3+d2*f2+d3*f1
  t4<-d1*f4+d2*f3+d3*f2+d4*f1
  t5<-d1+d2*f4+d3*f3+d4*f2+f1
  t6<-d2+d3*f4+d4*f3+f2
  t7<-d3+d4*f4+f3
  t8<-d4+f4
  u1<-d1*h1
  u2<-d2*h1+d1*h2
  u3<-d1*h3+d2*h2+d3*h1
  u4<-d1*h4+d4*h1+d2*h3+d3*h2
  u5<-d1+d2*h4+d3*h3+d4*h2+h1
  u6<-d2+d3*h4+d4*h3+h2
  u7<-d3+d4*h4+h3
  u8<-d4+h4
  # find the coefficients of the polynomial
  z1<-n1*m1*q1+n2*m2*r1+n3*m3*s1+n4*m4*t1+n5*m5*u1
  z2<-n1*(m1*q2-q1)+n2*(m2*r2-r1)+n3*(m3*s2-s1)+n4*(m4*t2-t1)+n5*(m5*u2-u1)
  z3<-n1*(m1*q3-q2)+n2*(m2*r3-r2)+n3*(m3*s3-s2)+n4*(m4*t3-t2)+n5*(m5*u3-u2)
  z4<-n1*(m1*q4-q3)+n2*(m2*r4-r3)+n3*(m3*s4-s3)+n4*(m4*t4-t3)+n5*(m5*u4-u3)
  z5<-n1*(m1*q5-q4)+n2*(m2*r5-r4)+n3*(m3*s5-s4)+n4*(m4*t5-t4)+n5*(m5*u5-u4)
  z6<-n1*(m1*q6-q5)+n2*(m2*r6-r5)+n3*(m3*s6-s5)+n4*(m4*t6-t5)+n5*(m5*u6-u5)
  z7<-n1*(m1*q7-q6)+n2*(m2*r7-r6)+n3*(m3*s7-s6)+n4*(m4*t7-t6)+n5*(m5*u7-u6)
  z8<-n1*(m1*q8-q7)+n2*(m2*r8-r7)+n3*(m3*s8-s7)+n4*(m4*t8-t7)+n5*(m5*u8-u7)
  z9<-n1*(m1-q8)+n2*(m2-r8)+n3*(m3-s8)+n4*(m4-t8)+n5*(m5-u8)
  z10<--(n1+n2+n3+n4+n5)
  c(z1,z2,z3,z4,z5,z6,z7,z8,z9,z10)
  # find the roots of the polynomial
  root<-polyroot(c(z1,z2,z3,z4,z5,z6,z7,z8,z9,z10))
  # find the real roots 
  Rert=Re(root)[abs(Im(root))<=0.000001]
  # find the number of real roots
  n=length(Rert)
  # find the values of sigma_1^2 correponding to the real value of mu
  m8=NULL
  for(i in 1:n)
  {
      m8[i]=(sum((g1-Rert[i])^2))/n1
  }
  # find the MLE of sigma_1
  m81=min(m8)
  m9=NULL
  for(i in 1:n)
  {
      m9[i]=(sum((g2-Rert[i])^2))/n2
  }
  m91=min(m9)
  m10=NULL
  for(i in 1:n)
  {
      m10[i]=(sum((g3-Rert[i])^2))/n3
  }
  m101=min(m10)
  m11=NULL
  for(i in 1:n)
  {
      m11[i]=(sum((g4-Rert[i])^2))/n4
  }
  m111=min(m11)
  m12=NULL
  for(i in 1:n)
  {
      m12[i]=(sum((g5-Rert[i])^2))/n5
  }
  m121=min(m12)
  # below steps find the MLE of the parameters over alternative space
  x0<-c(m1, m2, m3, m4, m5)
  s21<-sum((g1-m1)^2)/n1
  s22<-sum((g2-m2)^2)/n2
  s23<-sum((g3-m3)^2)/n3
  s24<-sum((g4-m4)^2)/n4
  s25<-sum((g5-m5)^2)/n5
  s0<-c(s21, s22, s23, s24,s25)
  x8<-x0
  repeat
  {
    w1=c(n1/s0[1],n2/s0[2],n3/s0[3],n4/s0[4],n5/s0[5])
    x6=pava(x0,w1)
    s11<-sum((g1-x6[1])^2)/n1
    s12<-sum((g2-x6[2])^2)/n2
    s13<-sum((g3-x6[3])^2)/n3
    s14<-sum((g4-x6[4])^2)/n4
    s15<-sum((g5-x6[5])^2)/n5
    v=c(s11, s12, s13, s14, s15)
    s0<-v
    x10=abs(x8-x6)
    p=max(x10[1], x10[2], x10[3], x10[4], x10[5], na.rm = FALSE)
    if(p<=0.0000001)
    {
      break
    }
    x8<-x6
  }
  LR=(sqrt(s0[1]/m81))^(n1)*(sqrt(s0[2]/m91))^(n2)*(sqrt(s0[3]/m101))^(n3)*(sqrt(s0[4]/m111))^(n4)*(sqrt(s0[5]/m121))^(n5)
  value=-2*log(LR)
  return(value)
}

fun4<-function(n1,n2,n3,n4,n5,mu1,mu2,mu3,mu4,mu5,v11,v12,v13,v14,v15)
{
  # find 10000 values of LR statistic
  out=replicate(10000,fun1(n1,n2,n3,n4,n5,mu1,mu2,mu3,mu4,mu5,v11,v12,v13,v14,v15))
  # count if the -2log(lambda) value is gretaer than the critical value
  a=0
  q=qchisq(0.05,df=3,lower.tail=FALSE)
  for(i in 1:10000)
  {
    if(out[i]>q)
      a=a+1
  }
  p=a/10000
  p
  
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
