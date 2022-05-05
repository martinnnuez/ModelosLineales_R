rm(list=ls());
library(nlme)
Tr=matrix(c(rep("A",4),rep("B",4),rep("C",4)),ncol=1)
X=matrix(c(1,5,10,15,1,5,10,15,1,5,10,15),ncol=1)
Y=matrix(c(39.7,42.5,44.6,47.9,45.1,47.3,50.1,53.1,39.5,42.9,46.1,46.8),ncol=1)
datos<-data.frame(Tr,X,Y)
colnames(datos)<-c("T","X","Y")

datos<-read.table("clipboard",header = T,sep="\t",dec = ".")
#colnames(datos)<-c("Tr","X","Y")
names(datos)
mean(datos$X)
mean(datos$Y)

modelo=gls(Y~T+X+T*X,data=datos)
#modelo=gls(Y~A+B+A*B+x+B*x,data=data)

b=modelo$coefficients;b

#Item 1
x=8
M=matrix(c(1,0,0,x,0,0,
           1,1,0,x,x,0,
           1,0,1,x,0,x),byrow=T,ncol=length(modelo$coefficients))

rownames(M)=c("A","B","C")

M%*%b

#Item 2
Q=matrix(c(1,0,-1),ncol=3,byrow=T)

Q%*%M%*%b

D=Q%*%M

#Item 3
varD=D%*%vcov(modelo)%*%t(D)
EED=sqrt(varD)

#Item 4

#A)Efecto del factor principal T:
x=mean(datos$X)
M%*%b
K=matrix(c(1,-1,0,
           1,0,-1),ncol=3,byrow=T)

K%*%M%*%b #Estimacion de las diferencias

H=K%*%M

W=t((H%*%b))%*%solve(H%*%vcov(modelo)%*%t(H))%*%(H%*%b)/nrow(H);W
1-pf(W,nrow(H),nrow(datos)-length(b))

anova(modelo,type="marginal")

###### diferencia de las tres ordenadas
b
K=matrix(c(1,-1,0,0,0,0,
           1,0,-1,0,0,0),ncol=6,byrow=T)

K%*%b #Estimacion de las diferencias

H=K

W=t((H%*%b))%*%solve(H%*%vcov(modelo)%*%t(H))%*%(H%*%b)/nrow(H);W
1-pf(W,nrow(H),nrow(datos)-length(b))

###### tres ordenadas IGUAL 0
b
K=matrix(c(0,1,0,0,0,0,
           0,0,1,0,0,0),ncol=6,byrow=T)

K%*%b #Estimacion de las diferencias

H=K

W=t((H%*%b))%*%solve(H%*%vcov(modelo)%*%t(H))%*%(H%*%b)/nrow(H);W
1-pf(W,nrow(H),nrow(datos)-length(b))


###### diferencia de ordenadas sin M
b
K=matrix(c(0,1,-1,0,0,0),ncol=6,byrow=T)

K%*%b #Estimacion de las diferencias

H=K

W=t((H%*%b))%*%solve(H%*%vcov(modelo)%*%t(H))%*%(H%*%b)/nrow(H);W
1-pf(W,nrow(H),nrow(datos)-length(b))

###### diferencia de las pendientes
b
K=matrix(c(0,0,0,1,0,0,
           0,0,0,1,1,0,
           0,0,0,1,0,1),ncol=6,byrow=T)

K%*%b #Estimacion de las diferencias

H=K

W=t((H%*%b))%*%solve(H%*%vcov(modelo)%*%t(H))%*%(H%*%b)/nrow(H);W
1-pf(W,nrow(H),nrow(datos)-length(b))

#C)
###### diferencia de las pendientes entre las 3
b
K=matrix(c(0,0,0,0,1,0,
           0,0,0,0,0,1),ncol=6,byrow=T)

K%*%b #Estimacion de las diferencias

H=K

W=t((H%*%b))%*%solve(H%*%vcov(modelo)%*%t(H))%*%(H%*%b)/nrow(H);W
1-pf(W,nrow(H),nrow(datos)-length(b))


#B)X
b
H=matrix(c(0,0,0,0,1,0,
           0,0,0,0,0,1),ncol=6,byrow=T)
W=t((H%*%b))%*%solve(H%*%vcov(modelo)%*%t(H))%*%(H%*%b)/nrow(H);W

1-pf(W,nrow(H),nrow(datos)-length(b))

anova(modelo)

summary(modelo)

#Item 5
b
K=matrix(c(0,0,0,0,0,2),ncol=6,byrow=T)

K%*%b #Estimacion de las diferencias

H=K

#Item 6
W=t((H%*%b))%*%solve(H%*%vcov(modelo)%*%t(H))%*%(H%*%b)/nrow(H);W
1-pf(W,nrow(H),nrow(datos)-length(b))

#Item 7
b
K=matrix(c(0,1,0,0,0,0,
           0,0,1,0,0,0,
           0,0,0,0,1,0,
           0,0,0,0,0,1),ncol=6,byrow=T)

K%*%b #Estimacion de las diferencias

H=K

W=t((H%*%b))%*%solve(H%*%vcov(modelo)%*%t(H))%*%(H%*%b)/nrow(H);W
1-pf(W,nrow(H),nrow(datos)-length(b))

#item 8
b
K=matrix(c(0,0,0,0,1,0,
           0,0,0,0,0,1),ncol=6,byrow=T)

K%*%b #Estimacion de las diferencias

H=K

W=t((H%*%b))%*%solve(H%*%vcov(modelo)%*%t(H))%*%(H%*%b)/nrow(H);W
1-pf(W,nrow(H),nrow(datos)-length(b))

#item 9
#a) 
b
P=matrix(c(1,1,0,8,8,0),ncol=6,byrow=T)

#Valor esperado:
Vesp=P%*%b

#Error estandar:
varP=P%*%vcov(modelo)%*%t(P)
EEP=sqrt(varP)

#Limites del intervalo de confianza al 95%
n=nrow(datos)
p=length(modelo$coefficients)
LI=Vesp-qt(0.975,(n-p))*EEP
LS=Vesp+qt(0.975,(n-p))*EEP
c(LI,LS)


