# actividad 1
Y=matrix(c(3,2,1,2,2,4,5,2,5,5,2,4),ncol=1)
A=factor(c(1,1,1,2,2,2,3,3,3,4,4,4))

# inciso b
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

# inciso c
t(X)%*%X

solve(t(X)%*%X) #es una matriz singular por lo tanto no acepta inversa

det(t(X)%*%X)#determinante 0 confirma que no admite inversa
#Nos damos cuenta de que se trata de un modelo de rango incompleto, debemos buscar otra solucion

# ALTERNATIVA 1: cálculo de una inversa generalizada para X'X
# Calcula Moore-Penrose generalized inverse of a matrix

library(MASS)
G=ginv(t(X)%*%X)
G
round(t(X)%*%X%*%G%*%t(X)%*%X,4)   # verificar condición de inversa generalizada


Xinv=ginv(X)
Xinv
round(X%*%Xinv%*%X,4)   # verificar condición de inversa generalizada devuelve el t(X)%*%X

round(X%*%G%*%t(X)%*%X,4) 

#Una de las posibles soluciones
bo=G%*%t(X)%*%Y
bo

# inciso d

q1=c(1,1,0,0,0); q1
q2=c(0,1,-1,0,0);q2
q3=c(0,1,0,0,0); q3
q4=c(0,1,0,-1,0);q4

H=G%*%t(X)%*%X ;H

# prueba de estimabilidad
round(q1%*%H,4)
round(q2%*%H,4)
round(q3%*%H,4)
round(q4%*%H,4)

# estimaciones de eso que queria ver si era funcion estimable
round(q1%*%bo,4)
round(q2%*%bo,4)
round(q4%*%bo,4)



# ALTERNATIVA 2: reparametrización del modelo

# OPCIÓN 1: 

Lp1=matrix(c(1,1,0,0,0, 
             1,0,1,0,0, 
             1,0,0,1,0, 
             1,0,0,0,1),byrow=T, ncol=5)
Lp1

round(Lp1%*%H,4)#verifica que son funciones estimables


# verificar que son funciones estimables

round(Lp1[1,]%*%H,4)    
round(Lp1[2,]%*%H,4)
round(Lp1[3,]%*%H,4)
round(Lp1[4,]%*%H,4)

U1=X%*%ginv(Lp1)
round(U1,4)


theta1=solve(t(U1)%*%U1)%*%t(U1)%*%Y
theta1 #estimaciones de los parametros, tita m1,m2,m3,m4 sin ordenada al origen

#lo mismo que hacer este modelo.
modelo1=lm(Y~-1+A)
modelo1$coefficients


# inciso d




# OPCIÓN 2: 

Lp2=matrix(c(1,1,0,0,0,
             0,-1,1,0,0, 
             0,-1,0,1,0, 
             0,-1,0,0,1),byrow=T, ncol=5)
Lp2

# verificar que son funciones estimables

round(Lp2[1,]%*%H,4)    
round(Lp2[2,]%*%H,4)
round(Lp2[3,]%*%H,4)
round(Lp2[4,]%*%H,4)

U2=X%*%ginv(Lp2)
round(U2,4)


# para obtener Lp a partir de U

U=X[,-2]
U  
  
Lp=ginv(ginv(X)%*%U)
round(Lp,4)


theta2=solve(t(U2)%*%U2)%*%t(U2)%*%Y
theta2

#lo mismo que hacer este modelo
modelo2=lm(Y~A)
modelo2$coefficients
model.matrix(Y~A)

# ALTERNATIVA 3: uso de restricciones en las soluciones del sistema

# opción 1: 

C1=matrix(c(1,0,0,0,0),ncol=1); C1
t(C1)%*%H    # es una función no estimable

XpX=t(X)%*%X ;XpX
XpXrestricciones1=rbind(XpX,t(C1))
XpXrestricciones1=cbind(XpXrestricciones1,c(C1,0))  
XpXrestricciones1

inversa=round(solve(XpXrestricciones1),4); inversa

XpYrestricciones=c(t(X)%*%Y,0)
XpYrestricciones
b.amp1=solve(XpXrestricciones1)%*%XpYrestricciones
round(b.amp1,4)

# inciso d

round(q1%*%b.amp1[-6],4)
round(q2%*%b.amp1[-6],4)
round(q4%*%b.amp1[-6],4)


# opción 2:

C2=matrix(c(0,1,0,0,0),ncol=1); C2
t(C2)%*%H    # es una función no estimable

XpX=t(X)%*%X ;XpX
XpXrestricciones2=rbind(XpX,t(C2))
XpXrestricciones2=cbind(XpXrestricciones2,c(C2,0))  
XpXrestricciones2

XpYrestricciones=c(t(X)%*%Y,0)
XpYrestricciones
b.amp2=solve(XpXrestricciones2)%*%XpYrestricciones
round(b.amp2,4)

# inciso d

round(q1%*%b.amp3[-6],4)
round(q2%*%b.amp3[-6],4)
round(q4%*%b.amp3[-6],4)



# INCISO E

# para probar que no hay efecto de tratamientos

K=matrix(c(0,3,-1,-1,-1, 
           0,0,1,-1,0, 
           0,0,0,1,-1),byrow=T, ncol=5)
K

K1=matrix(c(0,1,-1,0,0, 
            0,0,1,-1,0, 
            0,0,0,1,-1),byrow=T, ncol=5)
K1

m=matrix(c(0,0,0),byrow=T, ncol=1)
m

SCR=t(Y)%*%Y -t(Y)%*%X%*%G%*%t(X)%*%Y
SCR

Q=t(K%*%bo-m)%*%solve(K%*%G%*%t(K))%*%(K%*%bo-m)
Q

## o también

modelo=lm(Y~A)
model.matrix(modelo)

b=modelo$coefficients ;b
M=matrix(c(1,0,0,0,
           1,1,0,0,
           1,0,1,0,
           1,0,0,1), byrow=T, ncol=4); M
#medias de cada uno de los tratamientos continuamos trabajando con esta
M%*%modelo$coefficients


C=matrix(c(1,-1,0,0,
           1,0,-1,0,
           1,0,0,-1), byrow=T, ncol=4);C

H=C%*%M ;H


micov=H%*%vcov(modelo)%*%t(H)


w=(t(H%*%b)%*%solve(micov)%*%H%*%b)/nrow(H); w
1-pf(w,nrow(H),nrow(Y)-length(b)) 


X1=model.matrix(modelo);X1
P=X1%*%solve(t(X1)%*%X1)%*%t(X1)
I=diag(1,nrow(X1))

SCresidual=t(Y)%*%(I-P)%*%Y



# INCISO F
#  contraste para diferencia de medias tratamiento 1 vs tratamiento 2

K=matrix(c(0,1,-1,0,0),byrow=T, ncol=5)
K

m=matrix(c(0),byrow=T, ncol=1)
m

Q=t(K%*%bo-m)%*%solve(K%*%G%*%t(K))%*%(K%*%bo-m)
Q




###  o también:

modelo=lm(Y~A)
b=modelo$coefficients ;b

#Matriz que permite recuperar las medias
M=matrix(c(1,0,0,0,
           1,1,0,0,
           1,0,1,0,
           1,0,0,1), byrow=T, ncol=4); M

#medias de mi modelo
M%*%modelo$coefficients

#Permite comparar las medias
C=matrix(c(1,-1,0,0), byrow=T, ncol=4)
H=C%*%M
micov=H%*%vcov(modelo)%*%t(H)

w=(t(H%*%b)%*%solve(micov)%*%H%*%b)/nrow(H); w
1-pf(w,nrow(H),nrow(Y)-length(b))

H%*%b
micov=H%*%vcov(modelo2)%*%t(H)
sqrt(micov)

X1=model.matrix(modelo);X1
P=X1%*%solve(t(X1)%*%X1)%*%t(X1)
I=diag(1,nrow(X1))

SCresidual=t(Y)%*%(I-P)%*%Y




# OBTENCIÓN DE SUMAS DE CUADRADO POR TÉCNICA DE REDUCCIÓN

# Suma de cuadrados de un modelo con U + A (modelo completo)

SCR=t(Y)%*%Y -t(Y)%*%X%*%G%*%t(X)%*%Y
SCR

# Suma de cuadrados de un modelo con U 
XU=X[,-(2:5)]
XU

SCRU=t(Y)%*%Y -t(Y)%*%XU%*%ginv(t(XU)%*%XU)%*%t(XU)%*%Y
SCRU

# reduccion en la SCR

SCRU - SCR



