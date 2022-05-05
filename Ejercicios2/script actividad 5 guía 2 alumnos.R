# ejercicio 5 guia actividades N°2

Y=matrix(c(3,2,4,9,6,7,2,6,5,8),nrow=10,byrow=TRUE);Y   
X=matrix(c(1,1,3,1,1,4,1,3,7,1,7,9,1,8,7,1,7,6,1,4,5,1,6,8,1,6,5,1,9,7),nrow=10,byrow=TRUE);X
Xt=t(X)
Yt=t(Y)

# INCISO C)
b=(solve(Xt%*%X))%*%Xt%*%Y; b  # estimación de betas
P=X%*%(solve(Xt%*%X))%*%Xt; P   #  matriz de proyección P
Leverage=diag(P); Leverage
sum(Leverage)
V=((Yt%*%Y)-(Yt%*%P%*%Y))/7; V  #  estimación de varianza

# INCISO D)
Vb=c(V)*(solve(Xt%*%X)); Vb     # estimación de la varianza de los betas
of=c(V)*solve(Xt%*%X)

# INCISO E)
xit=matrix(c(1,5,7),nrow=1,byrow=TRUE); xit
xi=t(xit)
yihat=xit%*%b; yihat
hii=xit%*%(solve(Xt%*%X))%*%xi; hii
             
LI= yihat-qt(0.975,df=7)%*%(sqrt(V*hii)); LI            
LS= yihat+qt(0.975,df=7)%*%(sqrt(V*hii)); LS
LS-LI   # amplitud intervalo



# INCISO F)
xit=matrix(c(1,5.2,6.1),nrow=1,byrow=TRUE); xit
xi=t(xit)
yihat=xit%*%b; yihat
hii=xit%*%(solve(Xt%*%X))%*%xi; hii


LI= yihat-qt(0.975,df=7)%*%(sqrt(V*hii)); LI            
LS= yihat+qt(0.975,df=7)%*%(sqrt(V*hii)); LS
LS-LI   # amplitud intervalo
 


# INCISO G)
Ho=matrix(c(0,1,0),nrow=1,byrow=TRUE); Ho
SCHo=(t(Ho%*%b))*(solve(Ho%*%(solve(Xt%*%X))%*%(t(Ho))))*(Ho%*%b);  SCHo



# INCISO H)
Ho=matrix(c(0,1,0,0,0,1),nrow=2,byrow=TRUE); Ho
SCHo=(t(Ho%*%b))%*%(solve(Ho%*%(solve(Xt%*%X))%*%(t(Ho))))%*%(Ho%*%b);  SCHo


# INCISO I)
Ho=matrix(c(0,1,-1),nrow=1,byrow=TRUE); Ho
SCHo=(t(Ho%*%b))*(solve(Ho%*%(solve(Xt%*%X))%*%(t(Ho))))*(Ho%*%b);  SCHo



# INCISO J)

# CÁLCULO DE SUMAS DE CUADRADOS

#  suma de cuadrado total
SCt=Yt%*%Y; SCt  

#  suma de cuadrado explicada por el modelo
SCm=(Yt%*%P%*%Y); SCm

#  suma de cuadrado residual
SCr=(Yt%*%Y)-(Yt%*%P%*%Y); SCr

SCm + SCr


#  suma de cuadrado debida a Bo
Xo=as.matrix(X[1:10,1]);Xo
Xot=t(Xo); 
Po=Xo%*%(solve(Xot%*%Xo))%*%Xot; 
SCBo=Yt%*%Po%*%Y; SCBo

#  suma de cuadrado del modelo corregida por Bo
SCmc=Yt%*%(P-Po)%*%Y; SCmc

SCm -SCBo  

#  suma de cuadrado total corregida por Bo
SCtc=SCt - SCBo  ; SCtc


#inciso g)

Ho=matrix(c(0,1,0),nrow=1,byrow=TRUE); Ho
SCHo=(t(Ho%*%b))*(solve(Ho%*%(solve(Xt%*%X))%*%(t(Ho))))*(Ho%*%b);  SCHo
n=10
p=qr(X)$rank ; p
q=qr(Ho)$rank ; q

F=((n-p)/q)*SCHo/SCr ; F

pvalue= pf(F, q,n-p, lower.tail=FALSE) ; pvalue


# INCISO h)
Ho=matrix(c(0,1,0,0,0,1),nrow=2,byrow=TRUE); Ho
SCHo=(t(Ho%*%b))%*%(solve(Ho%*%(solve(Xt%*%X))%*%(t(Ho))))%*%(Ho%*%b);  SCHo

n=10
p=qr(X)$rank ; p
q=qr(Ho)$rank ; q

F=((n-p)/q)*SCHo/SCr ; F

pvalue= pf(F, q,n-p, lower.tail=FALSE) ; pvalue


# INCISO i)
Ho=matrix(c(0,1,-1),nrow=1,byrow=TRUE); Ho
SCHo=(t(Ho%*%b))*(solve(Ho%*%(solve(Xt%*%X))%*%(t(Ho))))*(Ho%*%b);  SCHo

n=10
p=qr(X)$rank ; p
q=qr(Ho)$rank ; q

F=((n-p)/q)*SCHo/SCr ; F

pvalue= pf(F, q,n-p, lower.tail=FALSE) ; pvalue




# INCISO K) intervalo de confianza para Beta1
l=matrix(c(0,1,0),nrow=1,byrow=TRUE) ; l
l%*%b

l%*%b - qt(0.975,df=7)%*%sqrt(l%*%Vb%*%t(l))
l%*%b + qt(0.975,df=7)%*%sqrt(l%*%Vb%*%t(l))


# INCISO L) intervalo de confianza para Beta1 - Beta2
l=matrix(c(0,1,-1),nrow=1,byrow=TRUE)  ; l
l%*%b

l%*%b - qt(0.975,df=7)%*%sqrt(l%*%Vb%*%t(l))
l%*%b + qt(0.975,df=7)%*%sqrt(l%*%Vb%*%t(l))






