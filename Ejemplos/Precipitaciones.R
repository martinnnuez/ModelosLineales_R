#Ejemplo PRECIPITACIONES
getwd()
setwd("C:/Users/Jose Maria/Desktop/MODELOS LINEALES/Ejemplos en R Di Rienzo")
datos<-read.table("precipitaciones.txt",header=T,dec = ",")
datos

modelo1=lm(PT~Años,data=datos,na.action=na.omit)

#agrego columnita predichos a datos
datos$Pred=predict(modelo1)

library(ggplot2)
ggplot(datos,aes(x=Años,y=Pred,group=Estacion,col=Estacion))+
  geom_line()+ #la linea representa a las predicciones del modelo y los puntos los originales
  geom_point(aes(x=Años,y=PT))
#regresion unico que hubo efecto año, estacion no teine influencia
#por eso sale una recta comun para todos los datos
#esta es la recta que yo genere

#modelo2 agrego esacion
modelo2=lm(PT~Años+Estacion,data=datos,na.action=na.omit)
datos$Pred2=predict(modelo2)
#ahora x tiene indicadora de cada regresion, R elimina la primera osea A , entonces aparece indicadora de condicion B

model.matrix(modelo2)
#la de estacion vale 0 para A y 1 donde comienzan datos de B
modelo1
#Coefficients:
# (Intercept)    Años  
#578.684        2.508  # al incio era 578 promedio y aumentaba 2.5 por año
anova(modelo1)
#p valor    0.071 no significativo

#modelo2 agregamos estacion es distinta a la del modelo 1.
ggplot(datos,aes(x=Años,y=Pred2,group=Estacion,col=Estacion))+
  geom_line()+
  geom_point(aes(x=Años,y=PT))
#obtenemos en prediccion 2 rectas paralelas

modelo2
#Coefficients:
#  (Intercept)  Años    EstacionB  
#505.608        2.076      153.419  


#agrego interaccion
modelo3=lm(PT~Años+Estacion+Años*Estacion,data=datos,na.action=na.omit)
datos$Pred3=predict(modelo3)
#ahora x tiene indicadora de cada regresion, R elimina la primera osea A , entonces aparece indicadora de condicion B

model.matrix(modelo3)
datos$Pred3=predict(modelo3)

ggplot(datos,aes(x=Años,y=Pred3,group=Estacion,col=Estacion))+
  geom_line()+
  geom_point(aes(x=Años,y=PT))

modelo3
#Coefficients:
 # (Intercept)     Años  
#573.949          -1.328  
#EstacionB  Años:EstacionB  
#19.852           6.369  

modelo3$coefficients[1]+modelo3$coefficients[3]
#ordenada al origen de B

modelo3$coefficients[2]+modelo3$coefficients[4]
#pendiente de B

#CUAL ES LA PRECIPITACION ESPERADA PARA 1990 EN B
EB=matrix(c(1,40,1,40),ncol=4)
VE=EB%*%modelo3$coefficients
#EB ES LA MATRIZ TAL QUE PREMULTIPLICANDO LOS COEFICIENTES ESTIMADOS ME DEVUELVE LA ESPERANZA PARA ESE CASO PARTICULAR

n=length(predict(modelo3))#con cuantos datos trabajo


#error estandar
EE=sqrt(EB%*%vcov(modelo3)%*%t(EB))

#INTERVALO DE CONFIANAZ
LI=VE-qt(0.975,modelo3$df.residual)*EE
LS=VE+qt(0.975,modelo3$df.residual)*EE
c(LI,LS)


#MEDIAS PARA CADA AÑO
n=length(predict(modelo3))
EB=cbind(integer(da),seq(0,40),1,seq(0,40))

#PARA PODER COMPARAR AL ULTIMO
n=length(predict(modelo3))
EB=cbind(1,datos$Años,as.numeric(datos$Estacion=="B"),datos$Años*as.numeric(datos$Estacion=="B"))
##EB=matrix(c(1,seq(0,40),1,seq(0,40)),ncol=4,byrow=T)
VE=EB%*%modelo3$coefficients
EE=sqrt(diag(EB%*%vcov(modelo3)%*%t(EB)))
LI=VE-qt(0.975,modelo3$df.residual)*EE
LS=VE+qt(0.975,modelo3$df.residual)*EE
Limites=cbind(VE,LI,LS)
LimitesLm=predict(modelo3,interval="confidence",level=0.95)
all(signif(Limites,5)==signif(LimitesLm,5))

#Ahora vamos a hacer las curvas de intervalo de confianza
#agrgamos las columnas LI Y LS para poder graficarlas
datos$LI=LI
datos$LS=LS
head(datos)


ggplot(datos,aes(x=Años,y=Pred3,group=Estacion,col=Estacion))+
  geom_line()+
  geom_point(aes(x=Años,y=PT),size=1)+
  geom_line(aes(x=Años,y=LI),lty=2)+ #agrega limites inferiores
  geom_line(aes(x=Años,y=LS),lty=2)+
  geom_line(aes(x=Años,y=LIP),lty=3)+ 
  geom_line(aes(x=Años,y=LSP),lty=3)+
  facet_wrap(~Estacion)+
  theme(legend.position="none")


#intervalo de prediccion

LIP=VE-qt(0.975,modelo3$df.residual)*sqrt(EE^2+sigma(modelo3)^2)
LSP=VE+qt(0.975,modelo3$df.residual)*sqrt(EE^2+sigma(modelo3)^2)

datos$LIP=LIP
datos$LSP=LSP
