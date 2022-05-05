# ejercicio 5 guia actividades N°2

Y=matrix(c(3,2,4,9,6,7,2,6,5,8),nrow=10,byrow=TRUE);Y   
X=matrix(c(1,1,3,1,1,4,1,3,7,1,7,9,1,8,7,1,7,6,1,4,5,1,6,8,1,6,5,1,9,7),nrow=10,byrow=TRUE);X

X=matrix(c(rep(1,10),1,1,3,7,8,7,4,6,6,9,3,4,7,9,7,6,5,8,5,7), nrow=10,byrow=F)

Xmodelo=matrix(c(1,1,3,7,8,7,4,6,6,9,3,4,7,9,7,6,5,8,5,7), nrow=10,byrow=F)

modelo=lm(Y~Xmodelo)

#C#Estime lo parametros del modelo
b=solve(t(X)%*%X)%*%t(X)%*%Y
modelo$coefficients

#D)#Varianza del estimador beta
s2=t(Y-(X%*%b))%*%(Y-(X%*%b))/(nrow(X)-ncol(X))
covb=1.602763*(solve(t(X)%*%X))

vcov(modelo)

diag(covb)

#E#Intervalo de confianza al 95 para prediccion
x=matrix(c(1,5,7),nrow=1)
VE=x%*%b

#Varianza de esta prediccion, por ser combinacion lineal del b que tiene distirbucion normal es
var=x%*%vcov(modelo)%*%t(x)

#el intervalo de confianza
n=nrow(X)
p=ncol(X)
LI=VE-qt(0.975,(n-p))*sqrt(var)
LS=VE+qt(0.975,(n-p))*sqrt(var)
c(LI,LS)

#F#Menor amplitud intervalo confianza en las medias
apply(X,2,mean)
xmin=matrix(c(1,5.2,6.1),nrow=1)
Ve=xmin%*%b

Var=xmin%*%vcov(modelo)%*%t(xmin)
Li=Ve-qt(0.975,(n-p))*sqrt(Var)
Ls=Ve+qt(0.975,(n-p))*sqrt(Var)
c(Li,Ls)

#G#prueba de hipotesis b1=0
H=matrix(c(0,1,0),nrow=1)
h=matrix(c(0),nrow=1)

W=t((H%*%b)-h)%*%solve(H%*%vcov(modelo)%*%t(H))%*%(H%*%b-h)/nrow(H)

1-pf(W,nrow(H),length(Y)-length(b))

#H#prueba de hipotesis b1=b2=0
H=matrix(c(0,1,0,0,0,1),nrow=2,ncol=3,byrow=T)
h=matrix(c(0,0),nrow=2)

W=t((H%*%b)-h)%*%solve(H%*%vcov(modelo)%*%t(H))%*%(H%*%b-h)/nrow(H)

1-pf(W,nrow(H),length(Y)-length(b))


#I#prueba de hipotesis b1=b2    b1-b2=0
H=matrix(c(0,1,-1),nrow=1)
h=matrix(c(0),nrow=1)

W=t((H%*%b)-h)%*%solve(H%*%vcov(modelo)%*%t(H))%*%(H%*%b-h)/nrow(H)

1-pf(W,nrow(H),length(Y)-length(b))

#K#intervalo confianza b1
x=matrix(c(0,1,0),nrow=1)
b1=x%*%b

#Varianza de esta prediccion, por ser combinacion lineal del b que tiene distirbucion normal es
var=x%*%vcov(modelo)%*%t(x)

#el intervalo de confianza
n=nrow(X)
p=ncol(X)
LI=b1-qt(0.975,(n-p))*sqrt(var)
LS=b1+qt(0.975,(n-p))*sqrt(var)
c(LI,LS)


#L#intervalo confianza b1-b2
x=matrix(c(0,1,-1),nrow=1)
b1mb2=x%*%b

#Varianza de esta prediccion, por ser combinacion lineal del b que tiene distirbucion normal es
var=x%*%vcov(modelo)%*%t(x)

#el intervalo de confianza
n=nrow(X)
p=ncol(X)
LI=b1mb2-qt(0.975,(n-p))*sqrt(var)
LS=b1mb2+qt(0.975,(n-p))*sqrt(var)
c(LI,LS)

#CALCULO DE LAS SUMAS DE CUADRADOS
n=nrow(X)
p=ncol(X)
I=diag(1,nrow(X))
P=X%*%(solve(t(X)%*%X)%*%t(X))

U=matrix(rep(1,nrow(X)),ncol=1) #como la corrijo
P1=U%*%solve(t(U)%*%U)%*%t(U)

SCT=t(Y)%*%Y
SCMnc=t(Y)%*%P%*%Y
SCR=t(Y)%*%(I-P)%*%Y
SCbo=t(Y)%*%P1%*%Y
SCMc=t(Y)%*%(P-P1)%*%Y

##############################Ejercicio 6
Y=matrix(c(25,29.5,34.8,28.4,33,38.9),ncol=1)
X=matrix(c(rep(1,6),0,1,2,0,1,2,0,0,0,1,1,1),ncol=3,nrow=6,byrow=F)

Xmodelo=matrix(c(0,1,2,0,1,2,0,0,0,1,1,1),ncol=2,nrow=6,byrow=F)

modelo=lm(Y~Xmodelo)

#B#Estime lo parametros del modelo
b=solve(t(X)%*%X)%*%t(X)%*%Y
modelo$coefficients

#C#Varianza del estimador beta
s2=t(Y-(X%*%b))%*%(Y-(X%*%b))/(nrow(X)-ncol(X))
covb=0.1702778*(solve(t(X)%*%X))

vcov(modelo)

diag(vcov(modelo))


#D#prueba de hipotesis
H=matrix(c(0,1,-2),nrow=1)
h=matrix(c(0),nrow=1)

W=t((H%*%b)-h)%*%solve(H%*%vcov(modelo)%*%t(H))%*%(H%*%b-h)/nrow(H)

1-pf(W,nrow(H),length(Y)-length(b))

SCH=t(H%*%b-h)%*%solve(H%*%solve(t(X)%*%X)%*%t(H))%*%(H%*%b-h)

#E#prueba de hipotesis
H=matrix(c(0,1,-1),nrow=1)
h=matrix(c(1.5),nrow=1)

W=t((H%*%b)-h)%*%solve(H%*%vcov(modelo)%*%t(H))%*%(H%*%b-h)/nrow(H)

1-pf(W,nrow(H),length(Y)-length(b))

SCH=t(H%*%b-h)%*%solve(H%*%solve(t(X)%*%X)%*%t(H))%*%(H%*%b-h)
