####################parcial1
#EJERCICIO 1
#1
rm(list=ls());
#a)Especificar la distribucion de z1 y z2
#Z
A=matrix(c(2,-3,
           3,0),nrow=2,ncol=2,byrow=T)
M=matrix(c(0,0),ncol=1)
E=matrix(c(3,1,
           1,3),ncol=2)
cte=matrix(c(1,2),ncol=1)

(A%*%M)+cte #vector de esperanzas de combinacion lineal
EZ=A%*%E%*%t(A) #vector de varianzas combinacion lineal


#Z1      C*X+c = Z(1)
C=matrix(c(2,-3),nrow=1)
M=matrix(c(0,0),ncol=1)
E=matrix(c(3,1,1,3),ncol=2)

(C%*%M)+1
C%*%E%*%t(C)

#Z2     C*X+c = Z(2)
C=matrix(c(3,0),nrow=1)
M=matrix(c(0,0),ncol=1)
E=matrix(c(3,1,1,3),ncol=2)

(C%*%M)+2
C%*%E%*%t(C)

#2) distribucion conjunta de z1 y z2   C*X+c = Z(1) de interes
C=matrix(c(2,-3,
         3,0),nrow=2,ncol=2,byrow=T)

cte=matrix(c(1,2),ncol=1)

#Esperanza de la conjunta
MZ=C%*%M+cte
#Varianza de la conjunta
EZ=C%*%E%*%t(C)

#b) 
#1) Es Q una chi cuadrado

#devuelve el rango de una matriz
la.qr <- qr(A)
la.qr$rank
#verifico idemopotencia
A=matrix(c(1/24,-1/72,
           -1/72,1/24),nrow=2,ncol=2,byrow=T)

#Verificar la distribucion de la forma cuadratica ZQZ
rango <- function(x) { return(qr(solve(x))$rank) };
esSimetrica <- function(x) { all(x == t(x)) };
esIdempotente <- function(x) { all(x %*% x == x) };
n=rango(A);
esSimetrica(A);
esIdempotente(A %*% EZ);
(A%*%EZ)%*%(A%*%EZ) #verificacion es idempontente
l=1/2 * t(MZ) %*% A %*% MZ; # parametro de no centralidad

# por lo tanto Y'AY ~ chi-cuad con 2 (rango(A)) g.l. y parametro de no centr. 0.07638889

#Esperanza y varianza forma cuadratica, que es una chi cuadrado no central:
EQ=n+2*l;EQ
VQ=2*(n+4*l);VQ

#Si A no es la inversa de EZ cuales son los grados de libertad de la Q
solve(EZ)
A

#EJERCICIO 2
rm(list=ls());
Y=matrix(c(3,2,5,6,8,6,4,8,9,4),ncol=1)
A=matrix(c(rep(1,5),rep(0,5)),ncol=1)
B=matrix(c(rep(0,5),rep(1,5)),ncol=1)
datos<-data.frame(Y,A,B)
#colnames(X)<-c("trat")
modelo=lm(Y~-1+A+B)
model.matrix(modelo)
b=modelo$coefficients
anova(modelo)

#Medias para cada nivel del factor:
M=matrix(c(1,0,
           0,1),ncol=2,byrow=T)
M%*%b

#Prueba de hipotesis sobre las pendientes sean iguales:
H=matrix(c(1,-1),ncol=2,byrow=T)

W=t((H%*%b))%*%solve(H%*%vcov(modelo)%*%t(H))%*%(H%*%b)/nrow(H)

1-pf(W,nrow(H),length(Y)-length(b))

#PROBLEMA 3
remove(list = ls())
Y=matrix(c(2,3,5,6,5,5,8,9,8),ncol=1,byrow=T)
X=matrix(c("A1","A1","A1","A2","A2","A2","A2","A2","A2",
           "B1","B1","B2","B1","B1","B1","B2","B2","B2"),ncol=2,byrow=F)
colnames(X)<-c("A","B")
datos<-data.frame(Y,X)
datos$A

#Modelo lineal no aditivo para los datos
modelo=lm(Y~A+B+A*B,data=datos)
model.matrix(modelo)

#Medias para cada tratamiento
b=modelo$coefficients;b
M=matrix(c(1,0,0,0,
           1,0,1,0,
           1,1,0,0,
           1,1,1,1),ncol=4,byrow=T)
rownames(M)=c("a1b1","a1b2","a2b1","a2b2")

M%*%b

#Matriz de incidencia
model.matrix(modelo)

#Prueba hipotesis no hay efecto del factor principal A
#1 calcular las medias de A
K=matrix(c(1/2,1/2,0,0,
           0,0,1/2,1/2),ncol=4,byrow=T)
K%*%M%*%b

#2 matriz de comparacion de medias
Q=matrix(c(1,-1),ncol=2)

#Estimacion de las diferencias
Q%*%K%*%M%*%b

H=Q%*%K%*%M

#Prueba de hipotesis
W=(t(H%*%b)%*%solve(H%*%vcov(modelo)%*%t(H))%*%(H%*%b))/nrow(H)

1-pf(W,nrow(H),nrow(datos)-length(b))

#anova(modelo) #dio mal y confundio mucho, hacerlo infostat


#Prueba hipotesis no hay efecto del factor principal B
#1 calcular las medias de B
K=matrix(c(1/2,0,1/2,0,
           0,1/2,0,1/2),ncol=4,byrow=T)
K%*%M%*%b

#2 matriz de comparacion de medias
Q=matrix(c(1,-1),ncol=2)

#Estimacion de las diferencias
Q%*%K%*%M%*%b

H=Q%*%K%*%M

#Prueba de hipotesis
W=(t(H%*%b)%*%solve(H%*%vcov(modelo)%*%t(H))%*%(H%*%b))/nrow(H)

1-pf(W,nrow(H),nrow(datos)-length(b))


#Prueba hipotesis no hay efecto de interaccion
M%*%b
#Se calcula directamente sobre medias de tratamiento
#2 matriz de comparacion de medias
Q=matrix(c(1,-1,-1,1),ncol=4)

#Estimacion de las diferencias
Q%*%M%*%b

H=Q%*%M

#Prueba de hipotesis
W=(t(H%*%b)%*%solve(H%*%vcov(modelo)%*%t(H))%*%(H%*%b))/nrow(H)

1-pf(W,nrow(H),nrow(datos)-length(b))


#Prueba hipotesis las medias de los tratamientos son iguales
#condice con la suma de cuadrados corregida por la media
#ANOVA MODELO INFOSTAT
M%*%b
#Se calcula directamente sobre medias de tratamiento
#2 matriz de comparacion de medias
Q=matrix(c(1,-1,0,0,
           1,0,-1,0,
           1,0,0,-1),ncol=4,byrow=T)

#Estimacion de las diferencias
Q%*%M%*%b

H=Q%*%M

#Prueba de hipotesis
W=(t(H%*%b)%*%solve(H%*%vcov(modelo)%*%t(H))%*%(H%*%b))/nrow(H)

1-pf(W,nrow(H),nrow(datos)-length(b))


#prueba de hipotesis ma1_b1=ma2_b2
Q=matrix(c(1,0,-1,0),ncol=4)

H=Q%*%M

#Prueba de hipotesis
W=(t(H%*%b)%*%solve(H%*%vcov(modelo)%*%t(H))%*%(H%*%b))/nrow(H)

1-pf(W,nrow(H),nrow(datos)-length(b))


################################parcial 2
#1
#invariante a cambios de escala
remove(list = ls())
Y=matrix(c(3,2,5,6,8,6,4,8,9,4),ncol=1)
X=matrix(c(rep("A",5),rep("B",5)),ncol=1)
datos<-data.frame(Y,X)
modelo<-lm(Y~X,data=datos)

datos2<-data.frame(2*Y,X)
colnames(datos2)<-c("dY","X")
modelo2<-lm(dY~X,data=datos2)

anova(modelo)
anova(modelo2)

anova(modelo)==anova(modelo2)
#Se puede ver que la prueba aplicada a las observaciones originales es identica a la prueba aplicada a y por un escalar distinto de cero.

#PROBLEMA 2
remove(list = ls())
Y <- matrix(c(20,25,40,50,32,42,62,82),nrow=8,byrow=TRUE); Y
XVals <- matrix(c(8,20,33,37,7,18,32,38),ncol=1)
#X <- matrix(c(c(rep(1,4), rep(0,4)), c(rep(0,4), rep(1,4)), c(XVals[1:4], rep(0,4)), c(rep(0,4), XVals[5:8])),nrow=8,byrow=F); X

W<-matrix(c(rep("A",4),rep("B",4)),ncol=1)

datos1<-data.frame(Y,W,XVals)
colnames(datos1)<-c("Y","Tr","X")

modelo1=lm(Y~-1+Tr+X+Tr*X,data=datos1)
b=modelo1$coefficients;b

#Prueba de hipotesis sobre pendientes iguales:
H=matrix(c(0,0,0,1),ncol=4,byrow=T)
H%*%b
W=(t(H%*%b)%*%solve(H%*%vcov(modelo1)%*%t(H))%*%(H%*%b))/nrow(H);W
1-pf(W,nrow(H),nrow(datos1)-length(b))

#Ambas pendientes simultaneamente iguales a 1:
H=matrix(c(0,0,1,0,
           0,0,1,1),ncol=4,byrow=T)
H%*%b

h=matrix(c(1,1),ncol=1,byrow=T)

W=(t(H%*%b-h)%*%solve(H%*%vcov(modelo1)%*%t(H))%*%(H%*%b-h))/nrow(H);W
1-pf(W,nrow(H),nrow(datos1)-length(b))


#Probar que las rectas son iguales, prueba en dos partes:
b
H=matrix(c(-1,1,0,0,
           0,0,0,1),ncol=4,byrow=T)

H%*%b
W=(t(H%*%b)%*%solve(H%*%vcov(modelo1)%*%t(H))%*%(H%*%b))/nrow(H);W
1-pf(W,nrow(H),nrow(datos1)-length(b))


################################parcial 3
#ejercicio 1
remove(list = ls())
Y=matrix(c(3,2,5,6,8,6,4,8,9,4),ncol=1)
X=matrix(c(rep("A",5),rep("B",5)),ncol=1)
datos<-data.frame(Y,X)
modelo<-lm(Y~X,data=datos)

#invariante a traslaciones
datos2<-data.frame(4+Y,X)
colnames(datos2)<-c("cY","X")
modelo2<-lm(cY~X,data=datos2)

anova(modelo)
anova(modelo2)

anova(modelo)==anova(modelo2)
#Se puede ver que la prueba aplicada a las observaciones originales es identica a la prueba aplicada a y sumado a un escalar.

#ejercicio 2
remove(list = ls())
#Determinar si la forma lineal es independiente de la forma cuadratica
V=matrix(c(1,1,
           1,2),ncol=2,nrow=2,byrow=T)

B=matrix(c(1/2,1/2,
           -1,1),ncol=2,nrow=2,byrow=T)

A=matrix(c(2,-1,
           -1,1/2),ncol=2,nrow=2,byrow=T)

#CASO 2: independencia forma lineal y forma cuadratica
B%*%V%*%A==0

#Por lo tanto no son independientes

#ejecicio 4 
remove(list = ls())
Y=matrix(c(3.1,4.0,2.0,2.2,5.0,6.8,5.2,3.0,9.0,9.0),ncol=1,byrow=T)
X=matrix(c("A1","A1","A1","A1","A1","A2","A2","A2","A2","A2",
           "B1","B1","B2","B2","B3","B1","B1","B2","B3","B3"),ncol=2,byrow=F)
colnames(X)<-c("A","B")
datos<-data.frame(Y,X)
datos$A

#Modelo lineal no aditivo para los datos
modelo=lm(Y~A+B+A*B,data=datos)
model.matrix(modelo)
b=modelo$coefficients;b

#Medias para cada tratamiento
b=modelo$coefficients;b
M=matrix(c(1,0,0,0,0,0,
           1,0,1,0,0,0,
           1,0,0,1,0,0,
           1,1,0,0,0,0,
           1,1,1,0,1,0,
           1,1,0,1,0,1),ncol=6,byrow=T)
rownames(M)=c("a1b1","a1b2","a1b3","a2b1","a2b2","a2b3")

M%*%b

#Matriz de incidencia
model.matrix(modelo)

#Combinacion lineal para el tratamiento a2_b3:
R=matrix(c(1,1,1,0,1,0),nrow=1,byrow=T)
R%*%b

#Prueba hipotesis las medias de los tratamientos son iguales
#condice con la suma de cuadrados corregida por la media
#ANOVA MODELO INFOSTAT
M%*%b
#Se calcula directamente sobre medias de tratamiento
#2 matriz de comparacion de medias
Q=matrix(c(1,-1,0,0,0,0,
           1,0,-1,0,0,0,
           1,0,0,-1,0,0,
           1,0,0,0,-1,0,
           1,0,0,0,0,-1),ncol=6,byrow=T)

#Estimacion de las diferencias
Q%*%M%*%b

H=Q%*%M

#Prueba de hipotesis
W=(t(H%*%b)%*%solve(H%*%vcov(modelo)%*%t(H))%*%(H%*%b))/nrow(H)

1-pf(W,nrow(H),nrow(datos)-length(b))



#Prueba hipotesis no hay efecto del factor principal A
#1 calcular las medias de A
K=matrix(c(1/3,1/3,1/3,0,0,0,
           0,0,0,1/3,1/3,1/3),ncol=6,byrow=T)
K%*%M%*%b

#2 matriz de comparacion de medias
Q=matrix(c(1,-1),ncol=2)

#Estimacion de las diferencias
Q%*%K%*%M%*%b

H=Q%*%K%*%M

#Prueba de hipotesis
W=(t(H%*%b)%*%solve(H%*%vcov(modelo)%*%t(H))%*%(H%*%b))/nrow(H)

1-pf(W,nrow(H),nrow(datos)-length(b))

#anova(modelo) #dio mal y confundio mucho, hacerlo infostat


#Prueba hipotesis no hay efecto del factor principal B
#1 calcular las medias de B
K=matrix(c(1/2,0,0,1/2,0,0,
           0,1/2,0,0,1/2,0,
           0,0,1/2,0,0,1/2),ncol=6,byrow=T)
K%*%M%*%b

#2 matriz de comparacion de medias
Q=matrix(c(1,-1,0,
           1,0,-1),ncol=3,byrow=T)

#Estimacion de las diferencias
Q%*%K%*%M%*%b

H=Q%*%K%*%M

#Prueba de hipotesis
W=(t(H%*%b)%*%solve(H%*%vcov(modelo)%*%t(H))%*%(H%*%b))/nrow(H)

1-pf(W,nrow(H),nrow(datos)-length(b))


#Prueba hipotesis no hay efecto de interaccion
M%*%b
#Se calcula directamente sobre medias de tratamiento
#2 matriz de comparacion de medias
Q=matrix(c(1,-1,0,-1,1,0,
           1,0,-1,-1,0,1),ncol=6,byrow=T)

#Estimacion de las diferencias
Q%*%M%*%b

H=Q%*%M

#Prueba de hipotesis
W=(t(H%*%b)%*%solve(H%*%vcov(modelo)%*%t(H))%*%(H%*%b))/nrow(H)

1-pf(W,nrow(H),nrow(datos)-length(b))


####################### parcial 4
I=diag(4)
solve(I)

#determinante de la identidad es 1
I=diag(75)
det(I)

#Aca me di cuenta que la inversa por la matriz da identidad
X=matrix(c(4,3,5,6,
           5,2,3,4,
           2,3,3,4,
           2,1,2,3),ncol=4,byrow=T)
Xi=solve(X)
round(Xi%*%X,4)

#Voy a hacer el modelo con vectores columna ortonormales y ver si la cov de sus parametros es 0.
Y=matrix(c(2,3),ncol=1,byrow=T)
X=matrix(c(1,0,
           0,1),ncol=2,byrow=T)
datos<-data.frame(Y,X)

#Modelo lineal no aditivo para los datos
modelo=lm(Y~-1+X1+X2,data=datos)
vcov(modelo)
modelo$coefficients

#por forma de resolucion
b <- solve(t(X) %*% X) %*% t(X) %*% Y; b  # beta sombrero
sigma(modelo)
I=diag(2)
Res=t(Y)%*%(I-X%*%solve(t(X%*%X))%*%X)%*%Y
n=length(Y)
p=length(b)


