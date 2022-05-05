# actividad 2
Y=matrix(c(6.29,6.38,6.25,6.32,6.44,6.29,5.80,5.92,5.78,5.95,6.05,5.89),ncol=1)
A=factor(c(1,1,1,1,1,1,2,2,2,2,2,2))
B=factor(c(1,1,1,2,2,2,1,1,1,2,2,2))


# inciso b
X=matrix(c(1,1,0,1,0, 1,1,0,1,0, 1,1,0,1,0, 1,1,0,0,1, 1,1,0,0,1, 1,1,0,0,1, 1,0,1,1,0, 1,0,1,1,0, 1,0,1,1,0, 1,0,1,0,1, 1,0,1,0,1, 1,0,1,0,1), byrow=T, ncol=5)
X


# inciso c



# ALTERNATIVA 1: cálculo de una inversa generalizada para X'X

library(MASS)
G=ginv(t(X)%*%X)
G

bo=G%*%t(X)%*%Y
bo

# inciso d

q1=c(0,0,0,1,-1) ;q1
q2=c(0,1,-1,1,-1) ; q2


H=G%*%t(X)%*%X

round(q1%*%H,4)
round(q2%*%H,4)

round(q1%*%bo,4)
round(q2%*%bo,4)



# ALTERNATIVA 2: reparametrización del modelo

# OPCIÓN 1: 

# para obtener Lp a partir de U

U=X[,-c(2,4)]
U  

Lp=ginv(ginv(X)%*%U)
round(Lp,4)



# verificar que son funciones estimables

round(Lp[1,]%*%H,4)    
round(Lp[2,]%*%H,4)
round(Lp[3,]%*%H,4)


theta=solve(t(U)%*%U)%*%t(U)%*%Y
theta


modelo=lm(Y~A+B)
modelo$coefficients
model.matrix(Y~A+B)



# ALTERNATIVA 3: uso de restricciones en las soluciones del sistema

# opción 1:

C=matrix(c(0,1,0,0,0, 0,0,0,1,0),ncol=2); C
t(C)%*%H    # es una función no estimable

XpX=t(X)%*%X ;XpX
XpXrestricciones=rbind(XpX,t(C))
XpXrestricciones=cbind(XpXrestricciones,rbind(C,0,0))  
XpXrestricciones

XpYrestricciones=c(t(X)%*%Y,0,0)
XpYrestricciones
b.amp=solve(XpXrestricciones)%*%XpYrestricciones
round(b.amp,4)


# inciso d

round(q1%*%b.amp[-c(6,7)],4)
round(q2%*%b.amp[-c(6,7)],4)



# INCISO E apartado a)

K=matrix(c(0,1,-1,0,0, 0,0,0,1,-1),byrow=T, ncol=5)
K
m=matrix(c(0,0),byrow=T, ncol=1)
m

SCR=t(Y)%*%Y -t(Y)%*%X%*%G%*%t(X)%*%Y
SCR

Q=t(K%*%bo-m)%*%solve(K%*%G%*%t(K))%*%(K%*%bo-m)
Q



Y


# INCISO E apartado b)

K=matrix(c(0,0,0,1,-1),byrow=T, ncol=5)
K
m=matrix(c(0),byrow=T, ncol=1)
m

Q=t(K%*%bo-m)%*%solve(K%*%G%*%t(K))%*%(K%*%bo-m)
Q



# OBTENCIÓN DE SUMAS DE CUADRADO POR TÉCNICA DE REDUCCIÓN.

# Suma de cuadrados de un modelo con U + A + B (modelo completo)

SCR=t(Y)%*%Y -t(Y)%*%X%*%G%*%t(X)%*%Y
SCR

# Suma de cuadrados de un modelo con U 
XU=X[,-(2:5)]
XU

SCRU=t(Y)%*%Y -t(Y)%*%XU%*%ginv(t(XU)%*%XU)%*%t(XU)%*%Y
SCRU

# Suma de cuadrados de un modelo con U + A 
XUA=X[,-(4:5)]
XUA

SCRUA=t(Y)%*%Y -t(Y)%*%XUA%*%ginv(t(XUA)%*%XUA)%*%t(XUA)%*%Y
SCRUA

# Suma de cuadrados de un modelo con U + B 
XUB=X[,-(2:3)]
XUB

SCRUB=t(Y)%*%Y -t(Y)%*%XUB%*%ginv(t(XUB)%*%XUB)%*%t(XUB)%*%Y
SCRUB

# reduccion en la SCR

SCRU - SCRUA      #  SC tipo I para A   SECUENCIALES
SCRUA - SCR       #  SC tipo I para B


SCRU - SCRUB      #  SC tipo I para B   SECUENCIALES
SCRUB - SCR       #  SC tipo I para A


SCRUB - SCR       #  SC tipo III para A  PARCIALES
SCRUA - SCR       #  SC tipo III para B


