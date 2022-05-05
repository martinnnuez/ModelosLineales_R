#Taller de software pero mas modelos lineales
remove(list = ls())
setwd("C:/Users/Jose Maria/Desktop/MODELOS LINEALES/SEMANA 4")
#   
data<-read.table("InteraccionyCov.txt",header = T,sep="\t",dec = ",")

modelo=lm(Y~A+B+A*B+x+B*x,data=data)

model.matrix(modelo)

medias=NULL
#corro para 2,4,6
x=6
b=modelo$coefficients

M=matrix(c(1,0,0,0,x,0,0,0,
           1,0,0,1,x,0,0,x,
           1,1,0,0,x,0,0,0,
           1,1,0,1,x,1,0,x,
           1,0,1,0,x,0,0,0,
           1,0,1,1,x,0,0,x),byrow=T,ncol=length(modelo$coefficients))

M%*%b

rownames(M)=c("a1_b1","a1_b2","a2_b1","a2_b2","a3_b1","a3_b2")
medias=rbind(medias,cbind(x,M%*%b))

#ahora lo quiero separar en b1 y a1, aplicando strsplit

MID=do.call("rbind",lapply(rownames(medias),function(x) strsplit(x,split="_")[[1]]))

mediascompletas=data.frame(MID,medias)

library(ggplot2)
ggplot(mediascompletas, aes(X1,V2,group=X2,color=x))+
  geom_point()


#como se si esos coeficientes son significativos
anova(modelo)
#conclucion covariable no interactua con b y no tiene efecto
b=modelo$coefficients

#H para la interaccion con la covariable 
H=matrix(c(0,0,0,0,0,0,0,1),ncol=8,byrow=T)
W=(t(H%*%b)%*%solve(H%*%vcov(modelo)%*%t(H))%*%(H%*%b))/nrow(H)
1-pf(W,nrow(H),nrow(data)-length(b))   


#Para el efecto de la covariable
H=matrix(c(0,0,0,0,1,0,0,0),ncol=8,byrow=T)
W=(t(H%*%b)%*%solve(H%*%vcov(modelo)%*%t(H))%*%(H%*%b))/nrow(H)
1-pf(W,nrow(H),nrow(data)-length(b))   
#0.78878
#no coinciden los test con el de anova(modelo)
#como esta afectando la covariable, esta como efecto principal y como interaccion!
#existe efecto de la covariable, es solo ese 0 o son ambos cero, interaccion y efecto principal,
#necesito dos hipotesis


#esos resultados me dan igual que tabla gls, comparar con esa.


anova(modelo) 
#estas tablas muestra otras cosas
anova(modelo,type="marginal")
#las tablas no son las mismas 

library(nlme)
modelo=gls(Y~A+B+A*B+x+B*x,data=data)
anova(modelo,type="marginal")

model.matrix(modelo,data=data)
#forma distinta de codificar el modelo en la matriz de incidencia
#cuando le pusimos covariable e interaccion dejo de funcionar la otra cosa.

#pero la forma de identificar los test de hipotesis da los de cuadrados de tipo 3 
#comparamos con tabla que no era la correcta






#########arrancamos 18-11-2020
x=3
M=matrix(c(1,0,0,0,x,0,0,0,
           1,0,0,1,x,0,0,x,
           1,1,0,0,x,0,0,0,
           1,1,0,1,x,1,0,x,
           1,0,1,0,x,0,0,0,
           1,0,1,1,x,0,0,x),byrow=T,ncol=length(modelo$coefficients))
M%*%b
rownames(M)=c("a1_b1","a1_b2","a2_b1","a2_b2","a3_b1","a3_b2")

#diferencia de medias y error estandar entre a1 y a2 dado que x=3


MA=matrix(c(1/2,1/2,0,0,0,0, #primeras 2 marginal a1 
           0,0,1/2,1/2,0,0), #3 y 4 marginal a2),
           ncol=6, byrow=T) 

#obtube las medias marginales de A
MA%*%M%*%b

#Matriz de comparacion
Q=matrix(c(1,-1),ncol=2, byrow=T)

#Esta es la diferencia de medias
difmed=Q%*%MA%*%M%*%b

#obtengo H
H=Q%*%MA%*%M

#Calculo del error estandar
EE=sqrt(H%*%vcov(modelo)%*%t(H))


#estadistico t equivalente a prueba F
(1-pt(difmed/EE,3))*2
#da no significativa





#diferencia de medias y error estandar entre a1 y a2 dado b=2 y x=3
x=3
M=matrix(c(1,0,0,0,x,0,0,0,
           1,0,0,1,x,0,0,x,
           1,1,0,0,x,0,0,0,
           1,1,0,1,x,1,0,x,
           1,0,1,0,x,0,0,0,
           1,0,1,1,x,0,0,x),byrow=T,ncol=length(modelo$coefficients))

rownames(M)=c("a1_b1","a1_b2","a2_b1","a2_b2","a3_b1","a3_b2")
M%*%b

#ahora no hay que sacar marginales
#Si no que buscar la matriz que me permite comparar esas dos medias directamente
Q=matrix(c(0,1,0,-1,0,0),ncol=6, byrow=T)

#Esta es la direferncia de estas medias
difmed=Q%*%M%*%b

#obtengo H
H=Q%*%M

#Calculo del error estandar
EE=sqrt(H%*%vcov(modelo)%*%t(H))

#estadistico t= difmed/EE
(1-pt(difmed/EE,3))*2
#da no significativa

(difmed/EE)**2

#calculamos la F
W=(t(H%*%b)%*%solve(H%*%vcov(modelo)%*%t(H))%*%(H%*%b))/nrow(H)

#calculamos el p valor
1-pf(W,nrow(H),nrow(data)-length(b)) 

#T al cuadrado da lo mismo que la F


#diferencia de medias y error estandar entre a1 y a2 y a3 dado b=2 y x=3 EFECTO DEL TRATAMIENTO A CON B2 Y X=3
x=3
M=matrix(c(1,0,0,0,x,0,0,0,
           1,0,0,1,x,0,0,x,
           1,1,0,0,x,0,0,0,
           1,1,0,1,x,1,0,x,
           1,0,1,0,x,0,0,0,
           1,0,1,1,x,0,0,x),byrow=T,ncol=length(modelo$coefficients))

rownames(M)=c("a1_b1","a1_b2","a2_b1","a2_b2","a3_b1","a3_b2")
M%*%b

#ahora no hay que sacar marginales
#Si no que buscar la matriz que me permite comparar esas dos medias directamente
Q=matrix(c(0,1,0,-1,0,0,
           0,1,0,0,0,-1),ncol=6, byrow=T)


#Esta es la direferncia de estas medias
Q%*%M%*%b

#obtengo H
H=Q%*%M
dim(H)

#calculamos la F
W=(t(H%*%b)%*%solve(H%*%vcov(modelo)%*%t(H))%*%(H%*%b))/nrow(H)

#calculamos el p valor
1-pf(W,nrow(H),nrow(data)-length(b)) 

























################script del profe
rm(list=ls())
library(nlme)
Data=read.table("clipboard",header = T,sep="\t")
Data$A=C(Data$A,contr.sum)
Data$B=C(Data$B,contr.sum)

modelo=gls(Y~A+B+A*B+x+B*x,data=Data)
modelo$coefficients
(model.matrix(modelo,data=Data))


#                      b1                  b2
# a1                  u+b1*x              u +   Ib2*b2+b1*x+b2*x
# 
# a2            u +  Ia2*a2+b1*x         u +  Ia2*a2 + Ib2*b2 + Ia2b2*a2b2+b1*x+b2*x
# 
# a3            u +  Ia3*a3+b1*x         u +  Ia3*a3 + Ib2*b2 + Ia3b2*a3b2+b1*x+b2*x
Medias=NULL
x=6
b=modelo$coefficients
M=matrix(c(1,0,0,0,x,0,0,0,
           1,0,0,1,x,0,0,x,
           1,1,0,0,x,0,0,0,
           1,1,0,1,x,1,0,x,
           1,0,1,0,x,0,0,0,
           1,0,1,1,x,0,1,x),byrow=T,ncol=length(modelo$coefficients))
rownames(M)=c("a1_b1","a1_b2","a2_b1","a2_b2","a3_b1","a3_b2")
Medias=rbind(Medias,cbind(x,M%*%b))

MID=do.call("rbind",lapply(rownames(Medias),function(x) strsplit(x,split="_")[[1]]))
MediaCompletas=data.frame(MID,Medias)

library(ggplot2)
MediaComple
ggplot(MediaCompletas,aes(X1,V2,group=X2,color=x))+
  geom_point(size=2)
anova(modelo,type="marginal")

H=matrix(c(0,0,0,0,0,0,0,1),ncol=8,byrow=T)

miF=(t(H%*%b)%*%solve(H%*%vcov(modelo)%*%t(H))%*%(H%*%b))/nrow(H)
1-pf(miF,nrow(H),3)

H=matrix(c(0,0,0,0,1,0,0,0),ncol=8,byrow=T)
miF=(t(H%*%b)%*%solve(H%*%vcov(modelo)%*%t(H))%*%(H%*%b))/nrow(H)
1-pf(miF,nrow(H),3)

