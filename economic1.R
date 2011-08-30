#p chart
p0=0.0136
q0=1-p0
p1=0.113
q1=1-p1
lambda=1/50
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

f(h=3,L=3.336,n=20)
f(h=8,L=3,n=250)

h=seq(1,10,by=.1); L=seq(2,4,by=.1)
cost.min=f(h=h[1],L=L[1],n=2)
  for(n in 2:300){
    mat=outer(h,L,FUN=f,n=n)
    if(cost.min>=min(mat)){
      n.min=n
      cost.min=min(mat)
      mat.min=mat
    }
}

cost.min   ## minimum cost
n.min  ##sample size
aa <- which(mat.min==min(mat.min),arr.ind=T)
h[aa[1,][1]]  ##hours between sample
L[aa[1,][2]]  ##number of standard deviation
mat.min


dim(mat.min)


####### An Economic Model of the xbar control chart ######
f2 <- function(h=0.76,L=2.99,n=5,a1=1,a2=.1,W=25,Y=50,a4=100,lambda=.05,delta=2,g=0.0167,D=1)
  {
    alpha=2*pnorm(-L)
    beta=pnorm(L-delta*sqrt(n))-pnorm(-L-delta*sqrt(n))
    tau <- (1-(1+lambda*h)*exp(-lambda*h))/(lambda*(1-exp(-lambda*h)))
    cost <- (a1+a2*n)/h + (a4*(h/(1-beta)-tau+g*n+D) + W + Y*alpha*exp(-lambda*h)/(1-exp(-lambda*h)))/(1/lambda+h/(1-beta)-tau+g*n+D)
    return(cost)
  }

## minimum cost design
h=seq(0.1,1,by=.01); L=seq(2,4.5,by=.01)
cost.min=f2(n=1,h=h[1],L=L[1])
  for(n in 1:20){
    mat=outer(h,L,FUN=f2,n=n)
    if(cost.min>=min(mat)){
      n.min=n
      cost.min=min(mat)
      mat.min=mat
    }
}

cost.min   ## minimum cost
n.min  ##sample size
aa <- which(mat.min==min(mat.min),arr.ind=T)
h[aa[1,][1]]  ##hours between sample
L[aa[1,][2]]  ##number of standard deviation
f2()

## Figure 9-24 Page 465 of SQC Douglas
cost.frame <- NULL
h=seq(0.1,1,by=.01); L=seq(2,4.5,by=.01)
  for(n in 1:20){
    mat=outer(h,L,FUN=f2,n=n)
    aa <- which(mat==min(mat),arr.ind=T)
    h[aa[1,][1]]  ##hours between sample
#    cat(n,L[aa[1,][2]],h[aa[1,][1]],min(mat),"\n")
    cost.frame <- rbind(cost.frame,c(n,L[aa[1,][2]],h[aa[1,][1]],min(mat)))
  }
colnames(cost.frame) <- c("n","Optimum L","Optimum h","Cost")
cost.frame
