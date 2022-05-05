#nuevo ejemplito
remove(list = ls())

set.seed(30)
#creamos base de datos
data=expand.grid(A=c("a1","a2","a3"),B=c("b1","b2"),r=c(1,2))
#le agrega la respuesta
data$Y=rnorm(nrow(data),30,3)
data=data[,-3]

modelo1=lm(Y~A+B,data=data)
model.matrix(modelo1)

#agrego interacciones, aparecerian 2 columnitas mas que serain los productos entre las columnas
modelo2=lm(Y~A+B+A*B,data=data)
model.matrix(modelo2)

anova(modelo2)
# cuenta de donde salen los grados de libertad, los de la interaccion de la cantidad de columnas de cada una de las variables individuales.

#como calculamos las medias de los tratamientos
#medias de los niveles y luego de las medias marginales

model.matrix(modelo2)[12,]
#matriz de incidencia para una observacion particular 


#matriz para devovler las medias

M= matrix(c(1,0,0,0,0,0,
            1,0,0,1,0,0,
            1,1,0,0,0,0,
            1,1,0,1,1,0,
            1,0,1,0,0,0,
            1,0,1,1,0,1),byrow=T,ncol=length(modelo2$coefficients))
rownames(M)=c("a1b1","a1b2","a2b1","a2b2","a3b1","a3b2")


#devuelve las medias de cada tratamiento 
M%*%modelo2$coefficients

#matriz K que por la de medias me permite obtener las medias marginales
K=matrix(c(1/2,1/2,0,0,0,0, #primeras 2 marginal a1 
    0,0,1/2,1/2,0,0, #3 y 4 marginal a2
    0,0,0,0,1/2,1/2),ncol=6, byrow=T) #marginal a3

#obtengo las medias marginales de A
K%*%M%*%modelo2$coefficients

#matriz facil de combinacion de medias
Q=matrix(c(1,-1,0,
           1,0,-1),ncol=3, byrow=T)

b=modelo2$coefficients
#esta es la H del test de hipotesis para el tratamiento A, efecto del factor principal A
H=Q%*%K%*%M

#Obtengo el estadistico para calcular el p valor para esta prueba de hipotesis.
W=(t(H%*%b)%*%solve(H%*%vcov(modelo2)%*%t(H))%*%(H%*%b))/nrow(H)

1-pf(W,nrow(H),nrow(data)-length(b))                      
#es nrow(data)




#lo mismo para el efecto principal de B
#matriz K que por la de medias me permite obtener las medias marginales
K=matrix(c(1/3,0,1/3,0,1/3,0, #marginal b1 
           0,1/3,0,1/3,0,1/3 #marginal b2
           ),ncol=6, byrow=T)

#obtengo las medias marginales de A
K%*%M%*%modelo2$coefficients

#matriz facil de combinacion de medias
Q=matrix(c(1,-1),ncol=2, byrow=T)

b=modelo2$coefficients
#esta es la H del test de hipotesis para el tratamiento A, efecto del factor principal A
H=Q%*%K%*%M

#Obtengo el estadistico para calcular el p valor para esta prueba de hipotesis.
W=(t(H%*%b)%*%solve(H%*%vcov(modelo2)%*%t(H))%*%(H%*%b))/nrow(H)

1-pf(W,nrow(H),nrow(data)-length(b))                          
#es nrow(data)




#y para la interaccion=
M

Q=matrix(c(1,-1,-1,1,0,0,
    1,-1,0,0,-1,1),byrow=T,ncol=6)

H=Q%*%M
#interesante da que los dos coeficinetes de los terminos de interaccion tienen que ser 0

W=(t(H%*%b)%*%solve(H%*%vcov(modelo2)%*%t(H))%*%(H%*%b))/nrow(H)
1-pf(W,nrow(H),nrow(data)-length(b))


#errores estandares para las A

K=matrix(c(1/2,1/2,0,0,0,0, #primeras 2 marginal a1 
           0,0,1/2,1/2,0,0, #3 y 4 marginal a2
           0,0,0,0,1/2,1/2),ncol=6, byrow=T) #marginal a3

#obtengo las medias marginales de A
K%*%M%*%modelo2$coefficients

#errores estandares para las A
K%*%M%*%vcov(modelo2)%*%t(K%*%M)

sqrt(diag(K%*%M%*%vcov(modelo2)%*%t(K%*%M)))
# modelo homocedastico dan los mismos errores estandares.




#y para la interaccion= PRUEBO NUEVO
M%*%modelo2$coefficients

Q=matrix(c(1,-1,-1,1,0,0,
           1,-1,0,0,-1,+1),byrow=T,ncol=6)

H=Q%*%M
#interesante da que los dos coeficinetes de los terminos de interaccion tienen que ser 0

W=(t(H%*%b)%*%solve(H%*%vcov(modelo2)%*%t(H))%*%(H%*%b))/nrow(H)
1-pf(W,nrow(H),nrow(data)-length(b))

