

#p chart
p0=0.0136
q0=1-p0
p1=0.113
q1=1-p1
L=3.336
n=20
  
n=20
h=2.88
L=3.336
lambda=1/50
Delta= 0.86
E=T0=T1=5/60
T2=45/60
delta1=1
delta2=0
C0=114.24
C1=949.20
Y=977.4
W=0
a=0
b=4.22

###
f <- function(h=2.88,L=3.336,n=20,p0=0.0136,p1=0.113){
  q0 <- 1-p0
  q1 <- 1-p1
  tau <- (1-(1+lambda*h)*exp(-lambda*h))/(lambda*(1-exp(-lambda*h)))
  s <- exp(-lambda*h)/(1-exp(-lambda*h))
  IL=max(c(0,floor(n*p0-L*sqrt(n*p0*q0))+1))
  IU=min(c(n,floor(n*p0+L*sqrt(n*p0*q0))))
  temp=0
  for(i in IL:IU)
    temp=temp+choose(n,i)*p0^i*q0^(n-i)
  alpha=1-temp
  temp=0
  for(i in IL:IU)
    temp=temp+choose(n,i)*p1^i*q1^(n-i)
  beta=temp
  ARL1 <- 1/alpha
  ARL2 <- 1/(1-beta)
  ECT <- 1/lambda+(1-delta1)*s*T0/ARL1 - tau + n*E + h*ARL2 + T1 + T2
  ECC <- C0/lambda + C1*(-tau+n*E+h*ARL2+delta1*T1+delta2*T2) + s*Y/ARL1+W+(a+b*n)*(1/lambda-tau+n*E+h*ARL2+delta1*T1+delta2*T2)/h
  cost <- ECC/ECT
  return(cost)
}

h=seq(1,3,by=.1); L=seq(2,4,by=.1)
cost.min=f(h=h[1],L=L[1],n=2)
  for(n in 2:200){
    mat=outer(h,L,FUN=f,n=n)
    if(cost.min>=min(mat)){
      n.min=n
      cost.min=min(mat)
      mat.min=mat
    }
}
n
cost.min
aa <- which(mat.min==min(mat.min),arr.ind=T)
aa[1,]
h

dim(mat.min)
