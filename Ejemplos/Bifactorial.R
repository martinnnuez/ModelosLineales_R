remove(list = ls())
setwd("C:/Users/Jose Maria/Desktop/MODELOS LINEALES/SEMANA 4")
data<-read.table("InteraccionyCov.txt",header = T,sep="\t",dec = ",")

modelo=lm(Y~A+B+A*B+x+B*x,data=data)

x=mean(data$x)
b=modelo$coefficients

M=matrix(c(1,0,0,0,x,0,0,0,
           1,0,0,1,x,0,0,x,
           1,1,0,0,x,0,0,0,
           1,1,0,1,x,1,0,x,
           1,0,1,0,x,0,0,0,
           1,0,1,1,x,0,1,x),byrow=T,ncol=length(modelo$coefficients))
rownames(M)=c("a1_b1","a1_b2","a2_b1","a2_b2","a3_b1","a3_b2")
M%*%b

library(nlme)
modelo2=gls(Y~A+B+A*B+x+B*x,data=data)
anova(modelo2,type="marginal")

#Efecto factor principal A
unique(data$A)

#1 calcular las medias marginales del factor
K=matrix(c(1/2,1/2,0,0,0,0, #primeras 2 marginal a1 
           0,0,1/2,1/2,0,0, #3 y 4 marginal a2
           0,0,0,0,1/2,1/2),ncol=6, byrow=T) #marginal a3

#obtengo las medias marginales de A
K%*%M%*%modelo$coefficients

#2 calcular la matriz H que permite comparar
#matriz facil de combinacion de medias
Q=matrix(c(1,-1,0,
           1,0,-1),ncol=3, byrow=T)

H=Q%*%K%*%M

H%*%b #estimacion de las diferencias de las medias de tratamiento A

#3 hago el test de hipotesis
#Obtengo el estadistico para calcular el p valor para esta prueba de hipotesis.
W=(t(H%*%b)%*%solve(H%*%vcov(modelo)%*%t(H))%*%(H%*%b))/nrow(H)

1-pf(W,nrow(H),nrow(data)-length(b))                      
#COMPARADO INFOSTAT ESTA BIEN


#Efecto factor principal B: (NO ME DA COMO INFOSTAT)
#1 calcular las medias marginales del factor
K=matrix(c(1/3,0,1/3,0,1/3,0, #b1
           0,1/3,0,1/3,0,1/3),ncol=6, byrow=T) #b2

#obtengo las medias marginales de B
K%*%M%*%modelo$coefficients

#2 calcular la matriz H que permite comparar
#matriz facil de combinacion de medias
Q=matrix(c(-1,1),ncol=2, byrow=T)

H=Q%*%K%*%M

H%*%b #estimacion de las diferencias de las medias de tratamiento A

#3 hago el test de hipotesis
#Obtengo el estadistico para calcular el p valor para esta prueba de hipotesis.
W=(t(H%*%b)%*%solve(H%*%vcov(modelo)%*%t(H))%*%(H%*%b))/nrow(H)

1-pf(W,nrow(H),nrow(data)-length(b))                      
#es nrow(data)


#Efecto de la interaccion A*B
M%*%b
#Matriz H directamente junta esos coeficientes
H=matrix(c(0,0,0,0,0,1,0,0,
           0,0,0,0,0,0,1,0),ncol=8, byrow=T)

W=(t(H%*%b)%*%solve(H%*%vcov(modelo)%*%t(H))%*%(H%*%b))/nrow(H)

1-pf(W,nrow(H),nrow(data)-length(b))    

#Intento obtener H a partir de comparacion
M%*%b
Q=matrix(c(1,-1,-1,1,0,0,
           1,-1,0,0,-1,1),ncol=6,byrow=T)

H=Q%*%M #Ahora si me da la H


#Efecto interaccion con covariable
#Matriz H que directamente junta ese coeficiente
H=matrix(c(0,0,0,0,0,0,0,1),ncol=8, byrow=T)

W=(t(H%*%b)%*%solve(H%*%vcov(modelo)%*%t(H))%*%(H%*%b))/nrow(H)

1-pf(W,nrow(H),nrow(data)-length(b))


#Efecto covariable
H=matrix(c(0,0,0,0,1,0,0,0),ncol=8, byrow=T)

W=(t(H%*%b)%*%solve(H%*%vcov(modelo)%*%t(H))%*%(H%*%b))/nrow(H)

1-pf(W,nrow(H),nrow(data)-length(b))


#Prueba de hipotesis para a1 dado x fijo en media
Q=matrix(c(1,0,-1,0,0,0,
           1,0,0,0,-1,0),ncol=6,byrow=T)

H=Q%*%M

W=(t(H%*%b)%*%solve(H%*%vcov(modelo)%*%t(H))%*%(H%*%b))/nrow(H)

1-pf(W,nrow(H),nrow(data)-length(b))


#Intervalo de confianza al 95% para los b
b
EE=sqrt(diag(vcov(modelo)))

n=nrow(data)
p=length(modelo$coefficients)

LI=b-qt(0.975,(n-p))*EE
LS=b+qt(0.975,(n-p))*EE
c(LI,LS)

#Intervalos de prediccion al 95% para los b
LIp=b-qt(0.975,(n-p))*(EE+(sigma(modelo)**2))
LSp=b+qt(0.975,(n-p))*(EE+(sigma(modelo)**2))
c(LIp,LSp)

#Prediccion particular e intervalos: a2=1, b2=1, x=8
P=matrix(c(1,1,0,1,8,1,0,8),ncol=8,byrow=T)

Est=P%*%b
Var=P%*%vcov(modelo)%*%t(P)
EE=sqrt(Var)

n=nrow(data)
p=length(modelo$coefficients)

LI=Est-qt(0.975,(n-p))*EE
LS=Est+qt(0.975,(n-p))*EE
c(LI,LS)

LIp=Est-qt(0.975,(n-p))*(EE+(sigma(modelo)**2))
LSp=Est+qt(0.975,(n-p))*(EE+(sigma(modelo)**2))
c(LIp,LSp)
