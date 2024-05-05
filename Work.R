#Defining the fucntions

Create_Tapestry_M = function(){
  
  M = matrix(0,16,40)
  
  r1  = c(2,13,14,18,34,35,39)
  r2  = c(2,5,11,17,24,27,35)
  r3  = c(6,8,9,12,23,30,32)
  r4  = c(1,3,5,15,21,25,29,32)
  r5  = c(6,10,14,17,21,31,36)
  r6  = c(4,12,18,20,24,36,37)
  r7  = c(8,10,25,27,38,39,40)
  r8  = c(11,13,19,22,23,36)
  r9  = c(1,16,22,28,31,32,38)
  r10 = c(1,4,8,19,26,31,37)
  r11 = c(3,10,16,17,25,33,34,37,40)
  r12 = c(2,7,15,16,20,23,29)
  r13 = c(3,7,12,13,21,26,27,28)
  r14 = c(14,15,24,26,30,33,38)
  r15 = c(5,9,19,20,28,34,35,40)
  r16 = c(4,6,7,9,22,29,33,39)
  
  index_list = list(r1,r2,r3,r4,r5,r6,r7,r8,r9,r10,r11,r12,r13,r14,r15,r16)
  
  for (i in 1:16){
    for (j in index_list[i]){
      M[i,j] = 1
    }
  }
  
  return(M)
}

Create_Block_M = function(){
  
  M = matrix(0,6,10)
  
  r1 = c(1,2,3,4,5)
  r2 = c(1,2,6,7,8)
  r3 = c(1,3,6,9,10)
  r4 = c(2,4,7,9,10)
  r5 = c(3,5,7,8,9)
  r6 = c(4,5,6,8,10)
  
  index_list = list(r1,r2,r3,r4,r5,r6)
  
  for (i in 1:6){
    for (j in index_list[i]){
      M[i,j] = 1
    }
  }
  
  return(M)
}

row_sum = function(M){
  rs = c()
  for (i in 1:nrow(M)){
    rs = append(rs,sum(M[i,]))
  }
  return(rs)
}

col_sum = function(M){
  cs = c()
  for (i in 1:ncol(M)){
    cs = append(cs,sum(M[,i]))
  }
  return(cs)
}

library(R1magic)

MNormE = function(x,x1,p=1){
  v = (sum((x-x1)^p))^(1/p)
  return(v)
}

MSE = function(x,x1){
  v = sqrt(sum((x-x1)^2))
  return(v)
}

AE = function(x,x_){
  v = sum(abs(x-x_))
  return(v)
}

VectorEstimate = function(M,y,T1,p,lambda=0.1){
  l1 = solveL1(M,y,T1,p,lambda)
  x_ = l1$estimate
  return(x_)
}

MakeBinary = function(v){
  n = length(v)
  for (i in 1:n){
    if (v[i]>0){
      v[i] = 1
    }
  }
  return(v)
}

IsZero = function(v){
  return(sum(v==0)==length(v))
}

Zero_Indices = function(v){
  n = length(v)
  v1 = c()
  for (i in 1:n){
    if (v[i] == 0){
      v1 = append(v1,i)
    }
  }
  return(v1)
}

NonZero_Indices = function(v){
  n = length(v)
  v1 = c()
  for (i in 1:n){
    if (v[i] > 0){
      v1 = append(v1,i)
    }
  }
  return(v1)
}

pl_X = function(M,x0,T1,p,lambda=0.1){
  
  y = M%*%x0
  y = MakeBinary(y)
  x_ = VectorEstimate(M,y,T1,p,lambda)
  
  n = length(x_)
  m1 = min(min(x0),min(x_))
  m2 = max(max(x0),max(x_))
  d = (m2-m1)/n
  
  plot( 1:n, seq(m1+d,m2,d), type = "n",xlab="",ylab="")
  title(main=lambda,
        xlab="Signal Component",ylab="Spike Value")
  lines(1:n, x0 , col = "red")
  lines(1:n, x_, col = "blue", cex = 1.5)
}

Compare_X = function(M,x0,T1,p,lambda1=0.1,lambda2=0.1,Method=MSE){
  
  y = M%*%x0
  
  x1 = VectorEstimate(M,y,T1,p,lambda1)
  x2 = VectorEstimate(M,y,T1,p,lambda2)
  
  par(mfrow = c(1,2))
  pl_X(M,x0,T1,p,lambda1) 
  pl_X(M,x0,T1,p,lambda2)
  
  IsBetter = function(x1,x2,f){
    
    v1 = f(x0,x1)
    v2 = f(x0,x2)
    
    if (v1>v2){
      return(c(v2,lambda2))
    }
    else if (v1<v2){
      return(c(v1,lambda1))
    }
    else{
      print("Equal")
    }
  }
  
  IsBetter(x1,x2,Method)
} 

pl_Y = function(M,y0,T1,p,lambda=0.1){
  
  x_ = VectorEstimate(M,y0,T1,p,lambda)
  y_ = M%*%x_
  n = length(y_)
  m1 = min(min(y0),min(y_))
  m2 = max(max(y0),max(y_))
  d = (m2-m1)/n

  plot( 1:n, seq(m1+d,m2,d), type = "n",xlab="",ylab="")
  title(main=lambda,
        xlab="Signal Component",ylab="Spike Value")
  lines(1:n, y0 , col = "red")
  lines(1:n, y_, col = "blue", cex = 1.5)
}

Compare_Y = function(M,y0,T1,p,lambda1=0.1,lambda2=0.1,Method=MSE){
  
  x1 = VectorEstimate(M,y0,T1,p,lambda1)
  x2 = VectorEstimate(M,y0,T1,p,lambda2)
  y1 = M%*%x1
  y2 = M%*%x2
  
  par(mfrow = c(1,2))
  pl_Y(M,y0,T1,p,lambda1) 
  pl_Y(M,y0,T1,p,lambda2)
  
  IsBetter = function(y1,y2,f){
    
    v1 = f(y0,y1)
    v2 = f(y0,y2)
    
    if (v1>v2){
      return(c(v2,lambda2))
    }
    else if (v1<v2){
      return(c(v1,lambda1))
    }
    else{
      print("Equal")
    }
  }
  
  IsBetter(y1,y2,Method)
} 


#Testing the functions on an example, x0 is sparse
M =Create_Tapestry_M()
c = sample.int(40,5)
x0 = matrix(0,40,1)
for (i in c){
  x0[i,1] = 1
}
y = M%*%x0
y = MakeBinary(y)
T1 = diag(40)
p = matrix(0,40,1)

#Define epsilon as difference between least expected non-zero entry of x_ and 
#greatest expected zero entry of x_ . Distribution of error(MSE) vs epsilon

vec = c()
errorvec = c()
for (i in 0:1000){
  x_ = VectorEstimate(M,y,T1,p,7*i/1000)
  m1 = min(x_[NonZero_Indices(x0)])
  M1 = max(x_[Zero_Indices(x0)])
  m = abs(M1-m1)
  vec = append(vec,m)
  e = MSE(x0,x_)
  errorvec = append(errorvec,e)
}
cor(vec,errorvec,method="pearson")

par(mfrow=c(1,1))
plot(1:length(vec),vec,type="n",xlab="Index",ylab="Epsilon")
title(main="Epsilon Distribution")
lines(1:length(vec),vec,col="blue")

plot(1:length(errorvec),errorvec,type="n",xlab="Index",ylab="Error")
title(main="Error Distribution")
lines(1:length(errorvec),errorvec,col="red")

#Some simulations

c = sample.int(40,5)
x0 = matrix(0,40,1)
for (i in c){
  x0[i,1] = 1
}
y = M%*%x0
y = MakeBinary(y)

max(VectorEstimate(M,y,T1,p,1))
min(VectorEstimate(M,y,T1,p,0.5))

vec1 = c()
vec11 = c()
for (i in 0:1000){
  x_ = VectorEstimate(M,y,T1,p,7*i/1000)
  if (IsZero(x_)==FALSE){
    m1 = sort(x_[NonZero_Indices(x0)])
    M1 = sort(x_[Zero_Indices(x0)],decreasing = TRUE)[1]
    j = which(m1>=M1)[1] - 1
    vec1 = append(vec1,j)
    vec11 = append(vec11,M1)
  }
}

vec2 = c()
vec22 = c()
for (i in 0:1000){
  x_ = VectorEstimate(M,y,T1,p,7*i/1000)
  if (IsZero(x_)==FALSE){
    m1 = sort(x_[NonZero_Indices(x0)])[1]
    M1 = sort(x_[Zero_Indices(x0)],decreasing = TRUE)
    j = which(M1<=m1)[1] - 1
    vec2 = append(vec2,j)
    vec22 = append(vec22,m1)
  }
}

windows()
par(mfrow=c(2,2))

a = table(vec2)
barplot(a,main="False Positives")

#a$percent = round(100*a$Freq/sum(a$Freq), digits = 1)
#a$label = paste(a$vec2," (", a$percent,"%)", sep = "")
#pie(a$Freq,labels = a$label,col=rainbow(length(a$vec2)),cex=1,radius=1,main="False Positives")

b = table(vec1)
barplot(b,main="False Negatives")
#b$percent = round(100*b$Freq/sum(b$Freq), digits = 1)
#b$label = paste(b$vec1," (", b$percent,"%)", sep = "")
#pie(b$Freq,labels = b$label,col=rainbow(length(b$vec1)),cex=1,radius=1,main="False Negatives")

plot(1:length(vec2),vec2,type="n",xlab="Lambda",ylab="False Positives")
title(main="Distribution of False Positives and Lambda")
lines(1:length(vec2),vec2,col="red")

plot(1:length(vec1),vec1,type="n",xlab="Lambda",ylab="False Negatives")
title(main="Distribution of False Negatives and Lambda")
lines(1:length(vec1),vec1,col="blue")

windows()
par(mfrow=c(1,2))
plot(1:length(vec22),vec22,type="n",xlab="Lambda",ylab="Epsilon")
title(main="Distribution of Possible Epsilon(FP) and Lambda")
lines(1:length(vec22),vec22,col="red")

plot(1:length(vec11),vec11,type="n",xlab="Lambda",ylab="Epsilon")
title(main="Distribution of Possible Epsilon(FN) and Lambda")
lines(1:length(vec11),vec11,col="blue")

#ROC curves
library(pROC)
library(R1magic)
x_ = VectorEstimate(M,y,T1,p,0.007)
df = data.frame("Predictor"=x0,"Responses"=x_)
plot.roc(df$Predictor,df$Responses,auc.polygon=TRUE,auc.polygon.col="red",max.auc.polygon=TRUE,
         max.auc.polygon.col="green",main="ROC Curves",print.auc=TRUE,percent=TRUE,direction="<",levels=c(0,1))
auc(df$Predictor,df$Responses,percent=TRUE,direction="<",levels=c(0,1))

c = sample.int(40,5)
x0 = matrix(0,40,1)
for (i in c){
  x0[i,1] = 1
}
y = M%*%x0
y = MakeBinary(y)
auc_vec=c()
for (i in 0:1000){
  x_ = VectorEstimate(M,y,T1,p,7*i/1000)
  df = data.frame("Predictor"=x0,"Responses"=x_)
  a = auc(df$Predictor,df$Responses,percent=TRUE,direction="<",levels=c(0,1))
  auc_vec = append(auc_vec,a)
}

a = (which(auc_vec==max(auc_vec)))*7/1000
a
max(auc_vec)
plot(1:length(auc_vec),auc_vec,type="n",xlab="Lambda",ylab="AUC")
title(main="Distribution of AUC and Lambda")
lines(1:length(auc_vec),auc_vec,col="red")

x_ = VectorEstimate(M,y,T1,p,a[1])
df = data.frame("Predictor"=x0,"Responses"=x_)
plot.roc(df$Predictor,df$Responses,auc.polygon=TRUE,auc.polygon.col="red",max.auc.polygon=TRUE,
         max.auc.polygon.col="green",main="ROC Curves",print.auc=TRUE,percent=TRUE,direction="<",levels=c(0,1))

#Trying on a BIBD from Wikipedia
Mb = Create_Block_M()
xb = matrix(0,10,1)
c = sample.int(10,2)
for (i in c){
  xb[i] = 1
}
yb = Mb%*%xb
yb = MakeBinary(yb)
Tb = diag(10)
pb = matrix(0,10,1)

auc_vec=c()
for (i in 0:1000){
  x_ = VectorEstimate(Mb,yb,Tb,pb,7*i/1000)
  df = data.frame("Predictor"=xb,"Responses"=x_)
  a = auc(df$Predictor,df$Responses,percent=TRUE,direction="<",levels=c(0,1))
  auc_vec = append(auc_vec,a)
}

max(auc_vec)
max((which(auc_vec==max(auc_vec))))*7/1000
plot(1:length(auc_vec),auc_vec,type="n",xlab="Lambda",ylab="AUC")
title(main="Distribution of AUC and Lambda")
lines(1:length(auc_vec),auc_vec,col="red")


#Doing it for 100 random tests and seeing how the best AUC is distributed
max_auc_vec = c()
for (i in 1:100){
  xb = matrix(0,10,1)
  c = sample.int(10,2)
  for (i in c){
    xb[i] = 1
  }
  yb = Mb%*%xb
  yb = MakeBinary(yb)
  
  auc_vec=c()
  for (i in 0:1000){
    x_ = VectorEstimate(Mb,yb,Tb,pb,7*i/1000)
    df = data.frame("Predictor"=xb,"Responses"=x_)
    a = auc(df$Predictor,df$Responses,percent=TRUE,direction="<",levels=c(0,1))
    auc_vec = append(auc_vec,a)
  }
  
  max_auc_vec = append(max_auc_vec,max(auc_vec))
}

max_auc_vec = round(max_auc_vec,2)
a1 = table(max_auc_vec)
d1 = as.data.frame(a1)
ylim <- c(0, 1.1*max(d1$Freq))
xx <- barplot(d1$Freq, xaxt = 'n', xlab = '', width = 0.85, ylim = ylim,
              main = "Distribution of AUC for different test vectors", 
              ylab = "Frequency")
text(x = xx, y = d1$Freq, label = d1$Freq, pos = 3, cex = 0.8, col = "red")
axis(1, at=xx, labels=d1$max_auc_vec, tick=FALSE, las=2, line=-0.5, cex.axis=1)


c1 = c(1,2,3,4)
c1 = append(c1,c(5,6,7,8))
c1

library(plot.matrix)
x <- matrix(runif(35)<0.5, ncol=5) 
plot(x)

sample(10,replace=TRUE)

