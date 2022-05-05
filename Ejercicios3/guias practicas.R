#GUIA,PRACTICA NUMERO 3:
rm(list=ls());
# actividad 1
Y=matrix(c(3,2,1,2,2,4,5,2,5,5,2,4),ncol=1)
A=factor(c(1,1,1,2,2,2,3,3,3,4,4,4))

X=matrix(c(1,1,0,0,0, 
           1,1,0,0,0, 
           1,1,0,0,0, 
           1,0,1,0,0, 
           1,0,1,0,0, 
           1,0,1,0,0,
           1,0,0,1,0, 
           1,0,0,1,0, 
           1,0,0,1,0, 
           1,0,0,0,1, 
           1,0,0,0,1, 
           1,0,0,0,1), byrow=T, ncol=5)
X

library(MASS)
G=ginv(t(X)%*%X)

bo=G%*%t(X)%*%Y
bo

q1=c(1,1,0,0,0); q1
q2=c(0,1,-1,0,0);q2
q3=c(0,1,0,0,0); q3
q4=c(0,1,0,-1,0);q4

round(q1%*%G%*%t(X)%*%X,4)==q1
round(q2%*%G%*%t(X)%*%X,4)==q2
round(q3%*%G%*%t(X)%*%X,4)==q3
round(q4%*%G%*%t(X)%*%X,4)==q4

#Parte pruebas de hipotesis por reparametrizacion
modelo=lm(Y~-1+A)
b=modelo$coefficients;b
X=model.matrix(modelo)

#Prueba hipotesis no hay efecto de tratamiento
K=matrix(c(1,-1,0,0,
           1,0,-1,0,
           1,0,0,-1),ncol=4,byrow=T); K

K%*%b #diferencias estimadas

H=K

#Hago la prueba de hipotesis:
micov=H%*%vcov(modelo)%*%t(H)

w=(t(H%*%b)%*%solve(micov)%*%H%*%b)/nrow(H); w
1-pf(w,nrow(H),nrow(Y)-length(b))

#Suma de cuadrados asociada Q
n=nrow(X)
p=ncol(X)
I=diag(1,nrow(X))
P=X%*%(solve(t(X)%*%X)%*%t(X))

Q=(t(H%*%b)%*%solve(H%*%solve(t(X)%*%X)%*%t(H))%*%H%*%b);Q
SCR=t(Y)%*%(I-P)%*%Y;SCR

#Prueba de hipotesis: 
H=matrix(c(1,-1,0,0),ncol=4)
H%*%b #estimacion diferencia
#Hago la prueba de hipotesis:
micov=H%*%vcov(modelo)%*%t(H)

w=(t(H%*%b)%*%solve(micov)%*%H%*%b)/nrow(H); w
1-pf(w,nrow(H),nrow(Y)-length(b))

#Suma de cuadrados asociada Q
n=nrow(X)
p=ncol(X)
I=diag(1,nrow(X))
P=X%*%(solve(t(X)%*%X)%*%t(X))

Q=(t(H%*%b)%*%solve(H%*%solve(t(X)%*%X)%*%t(H))%*%H%*%b);Q
SCR=t(Y)%*%(I-P)%*%Y;SCR


# actividad 2
rm(list=ls());
Y=matrix(c(6.29,6.38,6.25,6.32,6.44,6.29,5.80,5.92,5.78,5.95,6.05,5.89),ncol=1)
A=factor(c(1,1,1,1,1,1,2,2,2,2,2,2))
B=factor(c(1,1,1,2,2,2,1,1,1,2,2,2))

modelo=lm(Y~A+B)
b=modelo$coefficients;b
X=model.matrix(modelo)

M=matrix(c(1,0,0,
           1,0,1,
           1,1,0,
           1,1,1),ncol=3,byrow=T)
rownames(M)=c("a1_b1","a1_b2","a2_b1","a2_b2")
M%*%b
#No hay efecto de los factores horno y temperatura
#Factor temperatura:
#1 calcular las medias marginales para este factor
K=matrix(c(1/2,1/2,0,0,
           0,0,1/2,1/2),ncol=4,byrow=T)

K%*%M%*%b #medias estimadas para cada T

Q=matrix(c(1,-1),ncol=2,byrow=T)

Q%*%K%*%M%*%b #diferencia estimada

H=Q%*%K%*%M
#Hago la prueba de hipotesis:
micov=H%*%vcov(modelo)%*%t(H)

w=(t(H%*%b)%*%solve(micov)%*%H%*%b)/nrow(H); w
1-pf(w,nrow(H),nrow(Y)-length(b))

#Suma de cuadrados asociada Q
n=nrow(X)
p=ncol(X)
I=diag(1,nrow(X))
P=X%*%(solve(t(X)%*%X)%*%t(X))

Q=(t(H%*%b)%*%solve(H%*%solve(t(X)%*%X)%*%t(H))%*%H%*%b);Q
SCR=t(Y)%*%(I-P)%*%Y;SCR


#Factor horno:
M%*%b
#1 calcular las medias marginales para este factor
K=matrix(c(1/2,0,1/2,0,
           0,1/2,0,1/2),ncol=4,byrow=T)

K%*%M%*%b #medias estimadas para cada horno

Q=matrix(c(1,-1),ncol=2,byrow=T)

Q%*%K%*%M%*%b #diferencia estimada

H=Q%*%K%*%M
#Hago la prueba de hipotesis:
micov=H%*%vcov(modelo)%*%t(H)

w=(t(H%*%b)%*%solve(micov)%*%H%*%b)/nrow(H); w
1-pf(w,nrow(H),nrow(Y)-length(b))

#Suma de cuadrados asociada Q
n=nrow(X)
p=ncol(X)
I=diag(1,nrow(X))
P=X%*%(solve(t(X)%*%X)%*%t(X))

Q=(t(H%*%b)%*%solve(H%*%solve(t(X)%*%X)%*%t(H))%*%H%*%b);Q
SCR=t(Y)%*%(I-P)%*%Y;SCR


#Prubea modelo: no hay efecto de ningun factor 
K=matrix(c(1/2,1/2,0,0,
           0,0,1/2,1/2,
           1/2,0,1/2,0,
           0,1/2,0,1/2),ncol=4,byrow=T)

MM=K%*%M%*%b #medias temperatura y horno
rownames(MM)=c("T1","T2","H1","H2")

Q=matrix(c(1,-1,0,0,
           0,0,1,-1),ncol=4,byrow=T); H

Q%*%MM #diferencias estimadas con m1

H=Q%*%K%*%M
#Hago la prueba de hipotesis:
micov=H%*%vcov(modelo)%*%t(H)

w=(t(H%*%b)%*%solve(micov)%*%H%*%b)/nrow(H); w
1-pf(w,nrow(H),nrow(Y)-length(b))

#Suma de cuadrados asociada Q
n=nrow(X)
p=ncol(X)
I=diag(1,nrow(X))
P=X%*%(solve(t(X)%*%X)%*%t(X))

Q=(t(H%*%b)%*%solve(H%*%solve(t(X)%*%X)%*%t(H))%*%H%*%b);Q
SCR=t(Y)%*%(I-P)%*%Y;SCR


# actividad 3
rm(list=ls());
Y=matrix(c(6.29,6.38,6.25,6.32,6.44,6.29,5.80,5.92,5.78,5.95,6.05,5.89),ncol=1)
A=factor(c(1,1,1,1,1,1,2,2,2,2,2,2))
B=factor(c(1,1,1,2,2,2,1,1,1,2,2,2))

modelo=lm(Y~A+B+A*B)
b=modelo$coefficients;b
X=model.matrix(modelo)

M=matrix(c(1,0,0,0,
           1,0,1,0,
           1,1,0,0,
           1,1,1,1),ncol=4,byrow=T)
rownames(M)=c("a1_b1","a1_b2","a2_b1","a2_b2")
M%*%b #medias para cada tratamiento

#No hay efecto del factor temperatura:
K=matrix(c(1/2,1/2,0,0,
           0,0,1/2,1/2),ncol=4,byrow=T)
K%*%M%*%b#estimacion medias
Q=matrix(c(1,-1),ncol=2,byrow=T)
H=Q%*%K%*%M
#Hago la prueba de hipotesis:
micov=H%*%vcov(modelo)%*%t(H)

w=(t(H%*%b)%*%solve(micov)%*%H%*%b)/nrow(H); w
1-pf(w,nrow(H),nrow(Y)-length(b))

#Suma de cuadrados asociada Q
n=nrow(X)
p=ncol(X)
I=diag(1,nrow(X))
P=X%*%(solve(t(X)%*%X)%*%t(X))

Q=(t(H%*%b)%*%solve(H%*%solve(t(X)%*%X)%*%t(H))%*%H%*%b);Q
SCR=t(Y)%*%(I-P)%*%Y;SCR

#No hay efecto del factor horno:
K=matrix(c(1/2,0,1/2,0,
           0,1/2,0,1/2),ncol=4,byrow=T)
K%*%M%*%b#estimacion medias
Q=matrix(c(1,-1),ncol=2,byrow=T)
H=Q%*%K%*%M
#Hago la prueba de hipotesis:
micov=H%*%vcov(modelo)%*%t(H)

w=(t(H%*%b)%*%solve(micov)%*%H%*%b)/nrow(H); w
1-pf(w,nrow(H),nrow(Y)-length(b))

#Suma de cuadrados asociada Q
n=nrow(X)
p=ncol(X)
I=diag(1,nrow(X))
P=X%*%(solve(t(X)%*%X)%*%t(X))

Q=(t(H%*%b)%*%solve(H%*%solve(t(X)%*%X)%*%t(H))%*%H%*%b);Q
SCR=t(Y)%*%(I-P)%*%Y;SCR

#Efecto de la interaccion:
M%*%b

Q=matrix(c(-1,1,1,-1),ncol=4,byrow=T)

Q%*%M%*%b #diferencia estiamda

H=Q%*%M
#Hago la prueba de hipotesis:
micov=H%*%vcov(modelo)%*%t(H)

w=(t(H%*%b)%*%solve(micov)%*%H%*%b)/nrow(H); w
1-pf(w,nrow(H),nrow(Y)-length(b))

#Suma de cuadrados asociada Q
n=nrow(X)
p=ncol(X)
I=diag(1,nrow(X))
P=X%*%(solve(t(X)%*%X)%*%t(X))

Q=(t(H%*%b)%*%solve(H%*%solve(t(X)%*%X)%*%t(H))%*%H%*%b);Q
SCR=t(Y)%*%(I-P)%*%Y;SCR

#parte de las funciones estimables
#matriz de incidencia
X=matrix(c(1,1,0,1,0,1,0,0,0,
           1,1,0,1,0,1,0,0,0,
           1,1,0,1,0,1,0,0,0,
           1,1,0,0,1,0,1,0,0,
           1,1,0,0,1,0,1,0,0,
           1,1,0,0,1,0,1,0,0,
           1,0,1,1,0,0,0,1,0,
           1,0,1,1,0,0,0,1,0,
           1,0,1,1,0,0,0,1,0,
           1,0,1,0,1,0,0,0,1,
           1,0,1,0,1,0,0,0,1,
           1,0,1,0,1,0,0,0,1),ncol=9,byrow=T)
library(MASS)
G=ginv(t(X)%*%X)

bo=G%*%t(X)%*%Y
bo

#Obtengo las medias para cada tratamiento
C=matrix(c(1,1,0,1,0,1,0,0,0,
           1,1,0,0,1,0,1,0,0,
           1,0,1,1,0,0,0,1,0,
           1,0,1,0,1,0,0,0,1),ncol=9,byrow=T)
C%*%bo #medias para cada tratamiento

#Funciones estimables a partir del despeje en calculo
q1=c(0,0,0,1,-1,1,-1,0,0)
q2=c(0,1,-1,1,-1,1,0,0,-1)

round(q1%*%G%*%t(X)%*%X,4)==q1
round(q2%*%G%*%t(X)%*%X,4)==q2

#GUIA PRACTICA NUMERO 2:
rm(list=ls());
Y=matrix(c(3,2,4,9,6,7,2,6,5,8),nrow=10,byrow=TRUE);Y   
Xd=matrix(c(1,1,3,7,8,7,4,6,6,9,
           3,4,7,9,7,6,5,8,5,7), ncol=2,byrow=F)

modelo=lm(Y~X) 
b=modelo$coefficients;b
X=model.matrix(modelo)
var=vcov(modelo);var

#Intervalo confianza para prediccion:
P=matrix(c(1,5,7),ncol=3,byrow=T)
est=P%*%b
varest=P%*%var%*%t(P)
EEest=sqrt(varest)

n=nrow(X)
p=ncol(X)
LI=est-qt(0.975,(n-p))*EEest
LS=est+qt(0.975,(n-p))*EEest
c(LI,LS)

#Intervalo de confianza mas estrecho
med=apply(Xd,2,mean)
P=matrix(c(1,5.2,6.1),ncol=3,byrow=T)
est=P%*%b
varest=P%*%var%*%t(P)
EEest=sqrt(varest)

n=nrow(X)
p=ncol(X)
LI=est-qt(0.975,(n-p))*EEest
LS=est+qt(0.975,(n-p))*EEest
c(LI,LS)

#Prueba hipotesis b1=0
b
H=matrix(c(0,1,0),ncol=3,byrow=T)
H%*%b
#Hago la prueba de hipotesis:
micov=H%*%vcov(modelo)%*%t(H)

w=(t(H%*%b)%*%solve(micov)%*%H%*%b)/nrow(H); w
1-pf(w,nrow(H),nrow(Y)-length(b))

#Suma de cuadrados asociada Q
n=nrow(X)
p=ncol(X)
I=diag(1,nrow(X))
P=X%*%(solve(t(X)%*%X)%*%t(X))

Q=(t(H%*%b)%*%solve(H%*%solve(t(X)%*%X)%*%t(H))%*%H%*%b);Q
SCR=t(Y)%*%(I-P)%*%Y;SCR

#intervalo de confianza b1:
b
H=matrix(c(0,1,0),ncol=3,byrow=T)
est=H%*%b
#Hago la prueba de hipotesis:
var=H%*%vcov(modelo)%*%t(H)
EE=sqrt(var)
n=nrow(X)
p=ncol(X)
LI=est-qt(0.975,(n-p))*EE
LS=est+qt(0.975,(n-p))*EE
c(LI,LS)


#b1=b2=0
b
H=matrix(c(0,1,0,
           0,0,1),ncol=3,byrow=T)
H%*%b
#Hago la prueba de hipotesis:
micov=H%*%vcov(modelo)%*%t(H)

w=(t(H%*%b)%*%solve(micov)%*%H%*%b)/nrow(H); w
1-pf(w,nrow(H),nrow(Y)-length(b))

#Suma de cuadrados asociada Q
n=nrow(X)
p=ncol(X)
I=diag(1,nrow(X))
P=X%*%(solve(t(X)%*%X)%*%t(X))

Q=(t(H%*%b)%*%solve(H%*%solve(t(X)%*%X)%*%t(H))%*%H%*%b);Q
SCR=t(Y)%*%(I-P)%*%Y;SCR

#intervalo de confianza b1-b2:
b
H=matrix(c(0,1,-1),ncol=3,byrow=T)
est=H%*%b
#Hago la prueba de hipotesis:
var=H%*%vcov(modelo)%*%t(H)
EE=sqrt(var)
n=nrow(X)
p=ncol(X)
LI=est-qt(0.975,(n-p))*EE
LS=est+qt(0.975,(n-p))*EE
c(LI,LS)


#Ejercico 6= 
rm(list=ls());
Y=matrix(c(25,29.5,34.8,28.4,33,38.9),ncol=1,byrow=T)
X=matrix(c(0,1,2,0,1,2,
           0,0,0,1,1,1),ncol=2,byrow=F)

modelo=lm(Y~X)
b=modelo$coefficients;b
X=model.matrix(modelo)
vcov(modelo)

#Efecto del nitrogeno duplica potasio
H=matrix(c(0,1,-2),ncol=3,byrow=T)
H%*%b #diferencia estimada
#Hago la prueba de hipotesis:
micov=H%*%vcov(modelo)%*%t(H)

w=(t(H%*%b)%*%solve(micov)%*%H%*%b)/nrow(H); w
1-pf(w,nrow(H),nrow(Y)-length(b))

#Suma de cuadrados asociada Q
n=nrow(X)
p=ncol(X)
I=diag(1,nrow(X))
P=X%*%(solve(t(X)%*%X)%*%t(X))

Q=(t(H%*%b)%*%solve(H%*%solve(t(X)%*%X)%*%t(H))%*%H%*%b);Q
SCR=t(Y)%*%(I-P)%*%Y;SCR




#B1=B2+1.5
H=matrix(c(0,1,-1),nrow=1)
h=matrix(c(1.5),nrow=1)

W=t((H%*%b)-h)%*%solve(H%*%vcov(modelo)%*%t(H))%*%(H%*%b-h)/nrow(H);W

1-pf(W,nrow(H),length(Y)-length(b))
#Suma de cuadrados asociada Q
n=nrow(X)
p=ncol(X)
I=diag(1,nrow(X))
P=X%*%(solve(t(X)%*%X)%*%t(X))

Q=(t(H%*%b-h)%*%solve(H%*%solve(t(X)%*%X)%*%t(H))%*%(H%*%b-h));Q
SCR=t(Y)%*%(I-P)%*%Y;SCR


#guia practica 1:
setwd("C:/Users/Jose Maria/Desktop/MODELOS LINEALES/GUIA practica 1")
datos <- read.table("calefaccion.txt", header=TRUE);
