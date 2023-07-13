rm(list=ls())
# by fun1 and fun2 we find critical value

fun1<-function(n1,n2,n3,n4,v1,v2,v3,v4)
{
  # generating sample from normal distribution with mean zero and variance vi's calculated from fun3
  g1<-rnorm(n1,0,sqrt(v1));g2<-rnorm(n2,0,sqrt(v2));g3<-rnorm(n3,0,sqrt(v3));g4<-rnorm(n4,0,sqrt(v4))
  # steps 11-24 are for finding the Max-T statistic value
  X1=mean(g1);X2=mean(g2);X3=mean(g3);X4=mean(g4)
  S1=var(g1);S2=var(g2);S3=var(g3);S4=var(g4)
  v1<-sqrt(S1/n1+S2/n2);v2<-sqrt(S2/n2+S3/n3);v3<-sqrt(S3/n3+S4/n4)
  T1=(X2-X1)/v1;T2=(X3-X2)/v2;T3=(X4-X3)/v3
  T=min(T1,T2,T3,na.rm = FALSE)
  return(T)
}
fun2<-function(n1,n2,n3,n4,v1,v2,v3,v4)
{
  # find 500 max-T statistic values from fun1
  x<-replicate(1000,fun1(n1,n2,n3,n4,v1,v2,v3,v4))
  # arrange them in increasing order
  y<-sort(x,decreasing=FALSE)
  # gives the critical value
  c<-y[950]
  return(c)
}
#c<-fun2(n1=15,n2=15,n3=15,n4=15,v1=3.170854,v2=5.884281,v3=0.7796924,v4=2.109946)
#c
#c<-fun2(n1=15,n2=15,n3=15,n4=15,v1=1.3^2,v2=1.4^2,v3=1.2^2,v4=1)
#c
fun3<-function(n1,n2,n3,n4,mu1,mu2,mu3,mu4,v11,v12,v13,v14)
{
  # generate sample from normal distribution
  g1<-rnorm(n1,mu1,sqrt(v11));g2<-rnorm(n2,mu2,sqrt(v12));g3<-rnorm(n3,mu3,sqrt(v13))
  g4<-rnorm(n4,mu4,sqrt(v14))
  # steps from 44-57 are for finding max-T statistic value
  X1=mean(g1);X2=mean(g2);X3=mean(g3); X4=mean(g4)
  S1=var(g1);S2=var(g2);S3=var(g3);S4=var(g4)
  v1<-sqrt(S1/n1+S2/n2);v2<-sqrt(S2/n2+S3/n3);v3<-sqrt(S3/n3+S4/n4)
  T1=(X2-X1)/v1;T2=(X3-X2)/v2;T3=(X4-X3)/v3
  # give the statistic value
  T=min(T1,T2,T3,na.rm = FALSE)
  # find the critical value
  out<-fun2(n1,n2,n3,n4,S1,S2,S3,S4)
  # count if the statistic value is greater than the critical value
  a=0
  if(T>out)
    a=a+1
  return(a)
}
fun4<-function(n1,n2,n3,n4,mu1,mu2,mu3,mu4,v11,v12,v13,v14)
{
  # find the number of times the statistic value is greater than the critical value among 1000 values
  out<-replicate(1000,fun3(n1,n2,n3,n4,mu1,mu2,mu3,mu4,v11,v12,v13,v14))
  alpha<-sum(out)/1000
  alpha
  return(alpha)
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

