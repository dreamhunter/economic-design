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



#### economic design for xbar chart ####
eco.xbar <- function(h=seq(0.1,1,by=.01),L=seq(2,4.5,by=.01),n=1:20,a1=1,a2=.1,W=25,Y=50,a4=100,lambda=.05,delta=2,g=0.0167,D=1,...){
  # Economic design for xbar chart, minimizing the cost by finding the best sample interval h, control limit L and sample size n
  # Args:
  #   h: sample interval
  #   L: number of s.d. from control limits to center line
  #   n: sample size
  #   a1: fixed cost per sample
  #   a2: cost per unit sampled
  #   W: cost of finding an assignable cause
  #   Y: cost of investigating a false alarm
  #   a4: hourly penalty cost associated with production in the out-of-control state
  #   lambda: 1/mean time process is in control
  #   delta: number of s.d. slip when out of control
  #   g: the time required to take a sample and interpret the result
  #   D: the time required to find the assignable cause
  f2 <- function(h,L,n)
  {
    alpha=2*pnorm(-L)
    beta=pnorm(L-delta*sqrt(n))-pnorm(-L-delta*sqrt(n))
    tau <- (1-(1+lambda*h)*exp(-lambda*h))/(lambda*(1-exp(-lambda*h)))
    cost <- (a1+a2*n)/h + (a4*(h/(1-beta)-tau+g*n+D) + W + Y*alpha*exp(-lambda*h)/(1-exp(-lambda*h)))/(1/lambda+h/(1-beta)-tau+g*n+D)
    return(cost)
  }
  cost.frame <- NULL
  for(k in n){
    mat=outer(h,L,FUN=f2,n=k)
    aa <- which(mat==min(mat),arr.ind=T)
    cost.frame <- rbind(cost.frame,c(k,L[aa[1,][2]],h[aa[1,][1]],min(mat)))
  }
  colnames(cost.frame) <- c("n","Optimum L","Optimum h","Cost")
  rownames(cost.frame) <- rep("",length(n))
  optimum <- cost.frame[which(cost.frame[,4]==min(cost.frame[,4])),]
  par(mar=c(7.1,4.1,2.1,2.1))
  contour(h,L,outer(h,L,FUN=f2,n=optimum[1]),xlab="h",ylab="L",...)
  points(optimum[3],optimum[2],pch=3)
  mtext(sprintf('n=%s   Opt L=%s   Opt h=%s   Cost=%s',optimum[1], optimum[2], optimum[3],round(optimum[4],digits=4)),side=1,line=4.5)
  return(list(optimum,cost.frame))
}

eco.xbar(nlevels=50)
