# EJERCICIO 1 

#inciso e  (matriz de correlación)
V=matrix(c(3,1,0,1,2,1,0,1,4),nrow=3,byrow=TRUE);V
D=diag(diag(V));D
A=sqrt(solve(D));A
R=A%*%V%*%A ; R 

# inciso f (matriz de correlación parcial (1,2)/3
V11=V[1:2,1:2];V11
V12=as.matrix(V[1:2,3]);V12   
V21=V[3,1:2]; V21     
V22=V[3,3];V22

V11.2=V11-(V12%*%solve(V22)%*%V21); V11.2
E=diag(diag(V11.2)) ;E
F=sqrt(solve(E));F
RP=F%*%V11.2%*%F; RP 

# tarea: obtener la correlación parcial (1,3)/2



# EJERCICIO 2 

# inciso a y b
datos <- read.table("calefaccion.txt",header=TRUE)
M=colMeans(datos) ; M   
V=cov(datos); V      # insesgada
R=cor(datos); R

# para verificar
D=diag(diag(V));D
A=sqrt(solve(D));A
R1=A%*%V%*%A ; R 



# inciso c (matriz de correlacion parcial (1,2)/3   
# en forma manual
V11=V[1:2,1:2];V11
V12=as.matrix(V[1:2,3]);V12   
V21=V[3,1:2]; V21     
V22=V[3,3];V22

V11.2=V11-(V12%*%solve(V22)%*%V21); V11.2
E=diag(diag(V11.2)) ;E
F=sqrt(solve(E));F
RP=F%*%V11.2%*%F; RP 

# con software (todas: (1,2)/3;(1,3)/2 y (2,3)/1)
library(corpcor)
RP=cor2pcor(R); RP



# EJERCICIO 3 

# inciso c
V=matrix(c(2,-1,0,-1,4,0,0,0,1),nrow=3,byrow=TRUE);V
D=diag(diag(V))
A=sqrt(solve(D))
R=A%*%V%*%A ; R 

# para obtener la correlación parcial (1,2)/3
V11=V[1:2,1:2];V11
V12=as.matrix(V[1:2,3]);V12   
V21=V[3,1:2]; V21     
V22=V[3,3];V22

V11.2=V11-(V12%*%solve(V22)%*%V21); V11.2
E=diag(diag(V11.2)) ;E
F=sqrt(solve(E));F
RP=F%*%V11.2%*%F; RP    


# EJERCICIO 5 

V=matrix(c(7,2,1,2,7,-1,1,-1,4),nrow=3,byrow=TRUE);V
A=chol(V) ; A
t(A)%*%A


# EJERCICIO 6 

V=matrix(c(3,1,1,3),nrow=2,byrow=TRUE);V
A=chol(V) ; A
solve(t(A))



