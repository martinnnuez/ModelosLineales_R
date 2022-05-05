# actividad 3 por reparametrizacion
rm(list=ls());
Y=matrix(c(6.29,6.38,6.25,6.32,6.44,6.29,5.80,5.92,5.78,5.95,6.05,5.89),ncol=1)
A=factor(c(1,1,1,1,1,1,2,2,2,2,2,2))
B=factor(c(1,1,1,2,2,2,1,1,1,2,2,2))

modelo=lm(Y~A+B+A*B)
b=modelo$coefficients  ; b
X=model.matrix(modelo)

M=matrix(c(1,0,0,0,
           1,0,1,0,
           1,1,0,0,
           1,1,1,1), byrow=T, ncol=4); M
rownames(M)=c("a1_b1","a1_b2","a2_b1","a2_b2")
E=M%*%b

#Matriz de medias marginales:
K=matrix(c(1/2,1/2,0,0, #T1
           0,0,1/2,1/2, # t2
           1/2,0,1/2,0, #h1
           0,1/2,0,1/2),byrow = T,ncol=4)
rownames(K)=c("T1","T2","H1","H2")
MM=K%*%M%*%b

#Inciso E
#A)No hay efecto del factor temperatura
MM
#comparar las primeras 2
Q=matrix(c(1,-1,0,0),byrow = T,ncol=4)

H=Q%*%K%*%M
H

#Hago la prueba de hipotesis:
micov=H%*%vcov(modelo)%*%t(H)

w=(t(H%*%b)%*%solve(micov)%*%H%*%b)/nrow(H); w
1-pf(w,nrow(H),nrow(Y)-length(b))
#Suma de cuadrados asociada Q
n=nrow(X)
p=ncol(X)
I=diag(1,nrow(X))
P=X%*%(solve(t(X)%*%X)%*%t(X))

Q=(t(H%*%b)%*%solve(H%*%solve(t(X)%*%X)%*%t(H))%*%H%*%b)
SCR=t(Y)%*%(I-P)%*%Y;SCR


#B)No hay efecto del factor horno
MM
#comparar las primeras 2
Q=matrix(c(0,0,1,-1),byrow = T,ncol=4)

H=Q%*%K%*%M
H

#Hago la prueba de hipotesis:
micov=H%*%vcov(modelo)%*%t(H)

w=(t(H%*%b)%*%solve(micov)%*%H%*%b)/nrow(H); w
1-pf(w,nrow(H),nrow(Y)-length(b))

#C)No hay efecto de interaccion temperatura*horno
#No tengo que trabajar sobre la matriz de medias marginales, si no que la de por tratamiento
M%*%b

Q=matrix(c(1,-1,-1,1),byrow = T,ncol=4)

H=Q%*%M
H

#Hago la prueba de hipotesis:
micov=H%*%vcov(modelo)%*%t(H)

w=(t(H%*%b)%*%solve(micov)%*%H%*%b)/nrow(H); w
1-pf(w,nrow(H),nrow(Y)-length(b))



#ojo ver como calcula las sumas de cuadrados
n=nrow(X)
p=ncol(X)
I=diag(1,nrow(X))
P=X%*%(solve(t(X)%*%X)%*%t(X))

U=matrix(rep(1,nrow(X)),ncol=1) #como la corrijo
P1=U%*%solve(t(U)%*%U)%*%t(U)

SCT=t(Y)%*%Y
SCMnc=t(Y)%*%P%*%Y
SCR=t(Y)%*%(I-P)%*%Y;SCR #0.046 me da igual perfecto
SCbo=t(Y)%*%P1%*%Y
SCMc=t(Y)%*%(P-P1)%*%Y