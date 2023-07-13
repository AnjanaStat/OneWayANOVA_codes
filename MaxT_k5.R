rm(list=ls())
# fun1 and fun2 are for finding critical values
fun1<-function(n1,n2,n3,n4,n5,v1,v2,v3,v4,v5)
{
  # generating sample from normal distribution with mean zero and variance vi's, calculated from fun3
  g1<-rnorm(n1,0,sqrt(v1));g2<-rnorm(n2,0,sqrt(v2));g3<-rnorm(n3,0,sqrt(v3))
  g4<-rnorm(n4,0,sqrt(v4));g5<-rnorm(n5,0,sqrt(v5))
  # steps 12-29 for finding max-T statistic values
  X1=mean(g1);X2=mean(g2);X3=mean(g3);X4=mean(g4);X5=mean(g5)
  s1=var(g1);s2=var(g2); s3=var(g3);s4=var(g4);s5=var(g5)
  v1=sqrt(s1/n1+s2/n2);v2=sqrt(s2/n2+s3/n3);v3=sqrt(s3/n3+s4/n4);v4=sqrt(s4/n4+s5/n5)
  T1=(X2-X1)/v1;T2=(X3-X2)/v2;T3=(X4-X3)/v3; T4=(X5-X4)/v4
  # gives max-t statistic value
  A=max(T1,T2,T3,T4,na.rm=FALSE)
  return(A)
}
fun2<-function(n1,n2,n3,n4,n5,v1,v2,v3,v4,v5)
{
  # calculate 500 max-t statistic values from fun1
  x<-replicate(1000,fun1(n1,n2,n3,n4,n5,v1,v2,v3,v4,v5))
  # arrange them in increasing order
  y<-sort(x,decreasing=FALSE)
  # gives the critical value
  c<-y[950]
  return(c)
}
fun3<-function(n1,n2,n3,n4,n5,mu1,mu2,mu3,mu4,mu5,v11,v12,v13,v14,v15)
{
  # generate sample from normal distribution
  g1<-rnorm(n1,mu1,sqrt(v11));g2<-rnorm(n2,mu2,sqrt(v12))
  g3<-rnorm(n3,mu3,sqrt(v13));g4<-rnorm(n4,mu4,sqrt(v14));g5<-rnorm(n5,mu5,sqrt(v15))
  # steps 47-56 for finding max-T statistic
  X1=mean(g1);X2=mean(g2);X3=mean(g3);X4=mean(g4);X5=mean(g5)
  s1=var(g1);s2=var(g2);s3=var(g3);s4=var(g4);s5=var(g5)
  v1=sqrt(s1/n1+s2/n2);v2=sqrt(s2/n2+s3/n3);v3=sqrt(s3/n3+s4/n4);v4=sqrt(s4/n4+s5/n5)
  T1=(X2-X1)/v1;T2=(X3-X2)/v2;T3=(X4-X3)/v3;T4=(X5-X4)/v4
  # give the statistic value
  A=max(T1,T2,T3,T4,na.rm=FALSE)
  # calculate the critical value
  out<-fun2(n1,n2,n3,n4,n5,s1,s2,s3,s4,s5)
  #count if the statistic value greater than the critical value
  a=0
  if(A>out)
    a=a+1
  return(a)
}
fun4<-function(n1,n2,n3,n4,n5,mu1,mu2,mu3,mu4,mu5,v11,v12,v13,v14,v15)
{
  # find the number of times the statistic value is greater than the critical value among 1000 values
  out<-replicate(1000,fun3(n1,n2,n3,n4,n5,mu1,mu2,mu3,mu4,mu5,v11,v12,v13,v14,v15))
  p<-sum(out)/1000
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

