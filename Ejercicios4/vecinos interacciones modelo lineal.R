setwd("C:/Users/Jose Maria/Desktop/MODELOS LINEALES/SEMANA 4")
datos<-read.table("vecinos.txt",header = T,sep="\t",dec = ",")

names(datos)

modelo=lm(Increm.~Especie,data=datos)
model.matrix(modelo)
datos
anova(modelo)

#cuando tengo en cuenta vecinos la variabilidad explicada por esta covariable, es casi toda la SCE del modelo sin incluirlo
modelo=lm(Increm.~Especie+vecinos,data=datos)
model.matrix(modelo)
datos

anova(modelo)

mib=modelo$coefficients
#me falta la columna de la categorica que toma como referencia

#Mostro en infostat que el valor que toma la covariable es el valor medio, cunado predice media para categoria
x=mean(datos$vecinos)
mediaFlexuosa= mib[1]+mib[2]+mib[3]*0+mib[4]*x
#cuando uno calcula las medias esas medias dependen de la covariable y esas hay que especificarlas.

M=matrix(c(1,0,1,x,#nigra
           1,1,0,x,#flexuosa
           1,0,0,x),ncol=length(mib),byrow=T)#chilensis #media de cada especie

#Devuelve las medias de las especies
M%*%mib

#Errores estandares
sqrt(diag(M%*%vcov(modelo)%*%t(M)))

#cuando hacemos x=0
x=0
#las medias dan mas altas

x=5
#cambia mis medias y mis errores estandares 


#Las diferencais se manteinen iguales pq es un modelo aditivo sin interaccion
#ahora incluimos interaccion
modelo=lm(Increm.~Especie+vecinos+Especie*vecinos,data=datos)
model.matrix(modelo)

mib=modelo$coefficients
M=matrix(c(1,0,1,x,0,x,#nigra
           1,1,0,x,x,0,#flexuosa
           1,0,0,x,0,0),ncol=length(mib),byrow=T)

#Devuelve las medias de las especies
M%*%mib

#Errores estandares
sqrt(diag(M%*%vcov(modelo)%*%t(M)))

unique(datos$Especie)

plotdatos=expand.grid(Especie=c("flexuosa" , "chilensis" ,"nigra"),vecinos=seq(0,5,0.1))

plotdatos$predict=predict(modelo,newdata=plotdatos)

library(ggplot2)

ggplot(plotdatos,aes(vecinos,predict,group=Especie))+
  geom_line(aes(color=plotdatos$Especie))+
  scale_color_discrete(name="Especie")

#aca estan las 3 rectas que ajusto el modelo.
