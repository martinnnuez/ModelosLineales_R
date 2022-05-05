###################### Guia 1

# si no quiere andar la funcion rango entonces probar
# library("Matrix")
# rankMatrix(x)

###### Ejercicio 1
rm(list=ls());
#cargamos los datos
V <- matrix(c(3,1,0,
              1,2,1,
              0,1,4),nrow=3,byrow=TRUE); V
M <- matrix(c(0,0,0), ncol=1); M

#calculos auxiliares: a) Para escribir la formula de la distribucion
det(V)**(-0.5)
solve(V)

#Punto b.1) Marginal de X1
#matriz C que selecciona
C <- matrix(c(1,0,0), ncol=3); C
C%*%M
C%*%V%*%t(C)

#b.2) Marginal de X3, reordenar
Vr <- matrix(c(4,1,0,
               1,2,1,
               0,1,3),nrow=3,byrow=TRUE); Vr

Vaux <- V[c(3,2,1),];
Vre <- Vaux[,c(3,2,1)];
Vre

#misma matriz C que selecciona pq reordene
C%*%M
C%*%Vr%*%t(C)

#C.1)Distribucion conjunta X1,X2 no hace falta reordenar
#matriz C que selecciona
C2 <- matrix(c(1,0,0,
               0,1,0), ncol=3, byrow=T); C2
C2%*%M
C2%*%V%*%t(C2)

#C.2)Distribucion conjunta X1,X3 hace falta reordenar
#matriz C que selecciona es la misma pq reordene
Vr <- matrix(c(3,0,1,
               0,4,1,
               1,1,2),nrow=3,byrow=T); Vr

Vaux <- V[c(1,3,2),];
Vr <- Vaux[,c(1,3,2)];
Vr
C2%*%M
C2%*%Vr%*%t(C2)

#d.1)
#Matriz reordenada para X2/x1,x3
Vaux <- V[c(2,1,3),];
Vr <- Vaux[,c(2,1,3)];
Vr

V11 <- Vr[1,1]; V11
V12 <- Vr[1,2:3]; V12
V21 <- as.matrix(Vr[2:3,1]); V21
V22 <- as.matrix(Vr[2:3,2:3]); V22

#sirve para m1.2
V12%*%solve(V22)
#V11.2
V11-V12%*%solve(V22)%*%V21

#d.2)X1,X3/X2
Vaux <- V[c(1,3,2),];
Vr <- Vaux[,c(1,3,2)];
Vr

V11 <- Vr[1:2,1:2]; V11
V12 <- as.matrix(Vr[1:2,3]); V12
V21 <- Vr[3,1:2]; V21
V22 <- Vr[3,3]; V22

#sirve para m1.2
V12%*%solve(V22)
#V11.2
V11-V12%*%solve(V22)%*%V21

#### Inciso e (matriz de correlacion)
D <- diag(diag(V)); D
Daux <- sqrt(solve(D)); Daux
R <- Daux%*%V%*%Daux ; R

#### Inciso f.1) (matriz de correlacion parcial (1,2)/3), no hace falta ordenar nada
V
V11 <- V[1:2,1:2]; V11
V12 <- as.matrix(V[1:2,3]); V12
V21 <- V[3,1:2]; V21
V22 <- V[3,3]; V22

V11.2 <- V11-(V12%*%solve(V22)%*%V21); V11.2
D11.2 <- diag(diag(V11.2)); D11.2 #diagonal de V11.2
D11.2aux <- sqrt(solve(D11.2)); D11.2aux #inversa a la -1/2
R1.2 <- D11.2aux%*%V11.2%*%D11.2aux; R1.2

# f.2) correlacion parcial (1,3)/2
# sera lo mismo pero hay que reordenar V rotando columna 2 con 3 y luego fila 2 con 3
Vaux <- V[c(1,3,2),];
Vr <- Vaux[,c(1,3,2)];
Vr
#ASI INTERCAMBIA FILAS Y COLUMNAS!!!

V11 <- Vr[1:2,1:2]; V11
V12 <- as.matrix(Vr[1:2,3]); V12
V21 <- Vr[3,1:2]; V21
V22 <- Vr[3,3]; V22

V11.2 <- V11-(V12%*%solve(V22)%*%V21); V11.2
D11.2 <- diag(diag(V11.2)); D11.2
D11.2aux <- sqrt(solve(D11.2)); D11.2aux
R1.2 <- D11.2aux%*%V11.2%*%D11.2aux; R1.2

#g)
#para encontrar var(Z)
V
C<-matrix(c(3,-2,0),nrow=1,byrow=TRUE); C
C%*%V%*%t(C)
#encontrar mu a mano

#e)
C1=matrix(c(2,-3,0),ncol=3,byrow=T)
C2=matrix(c(0,0,3),ncol=3,byrow=T)

C1%*%V%*%t(C2)


###### Ejercicio 2
rm(list=ls());

#### Inciso a
getwd()
setwd("C:/Users/Jose Maria/Desktop/MODELOS LINEALES/GUIA practica 1")
datos <- read.table("calefaccion.txt", header=TRUE);
M <- apply(datos, 2, mean); M #medias por columnas
V <- cov(datos); V # insesgada

#### Inciso b
R <- cor(datos); R

# manualmente
D <- diag(diag(V));
Daux <- sqrt(solve(D));
R <- Daux %*% V %*% Daux; R

#### Inciso c (matriz de correlacion parcial (1,2)/3)
# en forma manual #no hace falta reordenar
V11 <- V[1:2,1:2]; V11
V12 <- as.matrix(V[1:2,3]); V12
V21 <- V[3,1:2]; V21
V22 <- V[3,3]; V22

V11.2 <- V11 - (V12 %*% solve(V22) %*% V21); V11.2 # pag 15
D11.2 <- diag(diag(V11.2)); D11.2
D11.2aux <- sqrt(solve(D11.2)); D11.2aux
RP <- D11.2aux %*% V11.2 %*% D11.2aux; RP

# con software (todas: (1,2)/3;(1,3)/2 y (2,3)/1)
library(corpcor);
RP <- cor2pcor(R); RP

###### Ejercicio 3
rm(list=ls());
V <- matrix(c(2,-1,0,-1,4,0,0,0,1),nrow=3,byrow=TRUE); V
M <- matrix(c(2,1,4),nrow=3,byrow=TRUE); M

#a)Distribucion conjunta de X1,X2
#No hace falta reordenar, matriz que selecciona
C<-matrix(c(1,0,0,
            0,1,0),nrow=2,byrow=TRUE); C

V11<-C%*%V%*%t(C);V11
M1<- C%*%M;M1

C2=matrix(c(0,0,1),ncol=3,byrow=T)
M2=C2%*%M

V22<-V[3,3]
V12<-V[1:2,3]
V21<-V[3,1:2]

#b)Distribucion condicional x1,x2/x3
#para el calculo de mu
M1
(V12%*%solve(V22))

#varianza
V11.2=V11-V12%*%solve(V22)%*%V21;V11.2

#### Inciso c
#matriz de correlacion
D <- diag(diag(V))
Daux <- sqrt(solve(D))
R <- Daux %*% V %*% Daux ; R 

# para obtener la correlacion parcial (1,2)/3
#No hace falta reordenar
V11 <- V[1:2,1:2]; V11
V12 <- as.matrix(V[1:2,3]); V12
V21 <- V[3,1:2]; V21
V22 <- V[3,3]; V22

V11.2 <- V11 - (V12 %*% solve(V22) %*% V21); V11.2
D11.2 <- diag(diag(V11.2)); D11.2
D11.2aux <- sqrt(solve(D11.2));
RP <- D11.2aux %*% V11.2 %*% D11.2aux; RP


###### Ejercicio 5 resuelto con descomposicion cholesky
rm(list=ls());
Vx <- diag(3);Vx
Mx <- matrix(c(0,0,0), nrow=3);Mx

Vy <- matrix(c(7, 2, 1, 
               2, 7, -1, 
               1, -1, 4), nrow=3, byrow=T) # Sea Y la que queremos conseguir (CX + c)
gammay <- chol(Vy);
C <- t(gammay);C

###### Ejercicio 6
rm(list=ls());
Vx <- matrix(c(3, 1, 1, 3), nrow=2, byrow=T);
Mx <- matrix(c(1, 2), nrow=2, byrow=T);

Vy <- diag(2);

gammax <- chol(Vx);
Ti <- solve(t(gammax));

Ti%*%Mx

###### Ejercicio 7
rm(list=ls());
C <- diag(2);
V <- matrix(c(3,1,1,3), nrow=2, byrow=T)
Vm <- C %*% V %*% t(C)
M

A <- solve(V);

rango <- function(x) { return(qr(solve(x))$rank) };
esIdempotente <- function(x) { all(x %*% x == x) };
esSimetrica <- function(x) { all(x == t(x)) };

esSimetrica(A);
rango(A);
esIdempotente(A %*% Vm);

0.5*(t())

###### Ejercicio 8 (ver si BY y Y'AY son indep)
# Y ~ N3(My, I)
rm(list=ls());
A <- matrix(c(1/6,1/3,1/6,1/3,2/3,1/3,1/6,1/3,1/6), nrow=3, byrow=T);
B <- matrix(c(2,-1,0,-1,0,1), nrow=2, byrow=T);
rango <- function(x) { return(qr(solve(x))$rank) };
esSimetrica <- function(x) { all(x == t(x)) };

esSimetrica(A);
rango(A); # Como no es de rango completo entonces veamos el caso general Y ~ N(M,V)
V <- diag(3);
all(B %*% V %*% A == 0);
B%*%V%*%A #caso 2 pq A no es de rango K

###### Ejercicio 9
# Y ~ N2(My, Vy)
rm(list=ls());
My <- matrix(c(1,3), nrow=2);
Vy <- matrix(c(1,1,1,2), nrow=2, byrow=T);
B <- matrix(c(0.5,0.5,-1,1), nrow=2, byrow=T);
A <- matrix(c(2,-1,-1,1), nrow=2, byrow=T);

#### Inciso a (ver si BY y Y'AY son indep)
esSimetrica <- function(x) { all(x == t(x)) };
rango <- function(x) { return(qr(solve(x))$rank) };

esSimetrica(A);
rango(A)
all(B %*% Vy %*% A == 0);
B %*% Vy %*% A #por lo tanto no son independientes.

#### Inciso b (buscar distribucion de Y'AY)
rango <- function(x) { return(qr(solve(x))$rank) };
esIdempotente <- function(x) { all(x %*% x == x) };
esSimetrica(A);
rango(A);
esIdempotente(A %*% Vy);
(A%*%Vy)%*%(A%*%Vy) #verificacion es idempontente
1/2 * t(My) %*% A %*% My; # parametro de no centralidad

# por lo tanto Y'AY ~ chi-cuad con 2 (rango(A)) g.l. y parametro de no centr. 2.5

###### Ejercicio 10 (ver si Y'BY y Y'AY son indep)
# Y ~ N3(My, Vy)
rm(list=ls());
My <- matrix(c(1,3), nrow=2);
Vy <- matrix(c(1,1,1,2), nrow=2, byrow=T);
B <- matrix(c(2,-1,-1,1), nrow=2, byrow=T);
A <- matrix(c(3,1/2,1/2,2), nrow=2, byrow=T);

esSimetrica <- function(x) { all(x == t(x)) };
esSimetrica(A);
esSimetrica(B);

all(A %*% Vy %*% B == 0) #por lo tanto no son independientes


