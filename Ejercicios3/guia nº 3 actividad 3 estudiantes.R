# actividad 3
Y=matrix(c(6.29,6.38,6.25,6.32,6.44,6.29,5.80,5.92,5.78,5.95,6.05,5.89),ncol=1)
A=factor(c(1,1,1,1,1,1,2,2,2,2,2,2))
B=factor(c(1,1,1,2,2,2,1,1,1,2,2,2))


# inciso b
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
           1,0,1,0,1,0,0,0,1), byrow=T, ncol=9)
X

# inciso c

# ALTERNATIVA 1: cálculo de una inversa generalizada para X'X

library(MASS)
G=ginv(t(X)%*%X)

bo=G%*%t(X)%*%Y
bo

# inciso d

q1=c(0,0,0,1,-1,1,-1,0,0)
q2=c(0,1,-1,1,-1,1,0,0,-1)
q3=c(1,0,1,0,1,0,0,0,1)

# para verificar si son funciones estimables

H=G%*%t(X)%*%X
H

round(q1%*%H,4)
round(q2%*%H,4)
round(q3%*%H,4)

round(q1%*%bo,4)
round(q2%*%bo,4)
round(q3%*%bo,4)

# ALTERNATIVA 2: reparametrización del modelo

# OPCIÓN 1: 

U=X[,-c(2,4,6,7,8)]
U 

Lp=ginv(ginv(X)%*%U)
round(Lp,4)

# verificación
round(X%*%ginv(Lp),4)


# verificar que son funciones estimables

round(Lp[1,]%*%H,4)    
round(Lp[2,]%*%H,4)
round(Lp[3,]%*%H,4)
round(Lp[4,]%*%H,4)

theta=solve(t(U)%*%U)%*%t(U)%*%Y
theta

# inciso d


modelo=lm(Y~A+B+A*B)
b=modelo$coefficients  ; b
model.matrix(Y~A+B+A*B)

M=matrix(c(1,0,0,0,
           1,0,1,0,
           1,1,0,0,
           1,1,1,1), byrow=T, ncol=4); M

M%*%modelo$coefficients

# medias marginales factor A 
D1=matrix(c(1/2,1/2,0,0,
           0,0,1/2,1/2), byrow=T, ncol=4); D1

D1%*%M%*%b


# medias marginales factor B 
D2=matrix(c(1/2,0,1/2,0,
           0,1/2,0,1/2), byrow=T, ncol=4); D2

D2%*%M%*%b



# ALTERNATIVA 3: uso de restricciones en las soluciones del sistema

# opción 1:

C=matrix(c(0,1,0,0,0,0,0,0,0,  0,0,0,1,0,0,0,0,0,  0,0,0,0,0,1,0,0,0,  0,0,0,0,0,0,1,0,0,  0,0,0,0,0,0,0,1,0),ncol=5); C3
t(C)
C%*%H    # es una función no estimable

XpX=t(X)%*%X ;XpX
XpXrestricciones=rbind(XpX,t(C))
XpXrestricciones=cbind(XpXrestricciones,rbind(C,0,0,0,0,0))  
XpXrestricciones

XpYrestricciones=c(t(X)%*%Y,0,0,0,0,0)
XpYrestricciones
b.amp=solve(XpXrestricciones)%*%XpYrestricciones
round(b.amp,4)

# inciso d

round(q1%*%b.amp[-c(10:14)],4)
round(q2%*%b.amp[-c(10:14)],4)






# INCISO E apartado a)

K=matrix(c(0,1,-1,0,0,0,0,0,0),byrow=T, ncol=9)
K
m=matrix(c(0),byrow=T, ncol=1)
m

SCR=t(Y)%*%Y -t(Y)%*%X%*%G%*%t(X)%*%Y
SCR

Q=t(K%*%bo-m)%*%solve(K%*%G%*%t(K))%*%(K%*%bo-m)
Q


# o bien ...


D1%*%M%*%b

C=matrix(c(1,-1),byrow=T, ncol=2); C

H=C%*%D1%*%M ;H

micov=H%*%vcov(modelo)%*%t(H)

w=(t(H%*%b)%*%solve(micov)%*%H%*%b)/nrow(H); w





# INCISO E apartado b)

K=matrix(c(0,0,0,1,-1,0,0,0,0),byrow=T, ncol=9)
K
m=matrix(c(0),byrow=T, ncol=1)
m

Q=t(K%*%bo-m)%*%solve(K%*%G%*%t(K))%*%(K%*%bo-m)
Q


# o bien ...



D2%*%M%*%b

C=matrix(c(1,-1),byrow=T, ncol=2); C

H=C%*%D2%*%M ;H

micov=H%*%vcov(modelo)%*%t(H)

w=(t(H%*%b)%*%solve(micov)%*%H%*%b)/nrow(H); w






# INCISO E apartado c)

K=matrix(c(0,0,0,0,0,1,-1,-1,1),byrow=T, ncol=9)
K
m=matrix(c(0),byrow=T, ncol=1)
m

Q=t(K%*%bo-m)%*%solve(K%*%G%*%t(K))%*%(K%*%bo-m)
Q



#  o bien ...


C=matrix(c(1,-1,-1,1),byrow=T, ncol=4); C

H=C%*%M ;H

micov=H%*%vcov(modelo)%*%t(H)

w=(t(H%*%b)%*%solve(micov)%*%H%*%b)/nrow(H); w




# INCISO E apartado d)

K=matrix(c(0,1,-1,0,0,0,0,0,0,  
           0,0,0,1,-1,0,0,0,0,  
           0,0,0,0,0,1,-1,0,0),byrow=T, ncol=9)
K
m=matrix(c(0,0,0),byrow=T, ncol=1)
m

Q=t(K%*%bo-m)%*%solve(K%*%G%*%t(K))%*%(K%*%bo-m)
Q





# o bien ...

C=matrix(c(1,-1,0,0,
           1,0,-1,0,
           1,0,0,-1),byrow=T, ncol=4); C

H=C%*%M ;H

micov=H%*%vcov(modelo)%*%t(H)

w=(t(H%*%b)%*%solve(micov)%*%H%*%b)/nrow(H); w









# OBTENCIÓN DE SUMAS DE CUADRADO POR TÉCNICA DE REDUCCIÓN.

# Suma de cuadrados de un modelo con U + A + B + C=A*B (modelo completo)

SCR=t(Y)%*%Y -t(Y)%*%X%*%ginv(t(X)%*%X)%*%t(X)%*%Y
SCR

# Suma de cuadrados de un modelo con U 
XU=X[,-(2:9)]
XU
SCRU=t(Y)%*%Y -t(Y)%*%XU%*%ginv(t(XU)%*%XU)%*%t(XU)%*%Y
SCRU


# Suma de cuadrados de un modelo con U + A 
XUA=X[,-(4:9)]
XUA
SCRUA=t(Y)%*%Y -t(Y)%*%XUA%*%ginv(t(XUA)%*%XUA)%*%t(XUA)%*%Y
SCRUA

# Suma de cuadrados de un modelo con U + B 
XUB=X[,-c(2,3,6:9)]
XUB
SCRUB=t(Y)%*%Y -t(Y)%*%XUB%*%ginv(t(XUB)%*%XUB)%*%t(XUB)%*%Y
SCRUB


# Suma de cuadrados de un modelo con U + A + B 
XUAB=X[,-(6:9)]
XUAB
SCRUAB=t(Y)%*%Y -t(Y)%*%XUAB%*%ginv(t(XUAB)%*%XUAB)%*%t(XUAB)%*%Y
SCRUAB



# reduccion en la SCR

SCRU - SCRUA       #  SC tipo I para A    
SCRUA - SCRUAB     #  SC tipo I para B
SCRUAB - SCR       #  SC tipo I para A*B





