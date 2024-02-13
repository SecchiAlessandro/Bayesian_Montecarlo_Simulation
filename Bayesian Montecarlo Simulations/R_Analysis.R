rm(list=ls())
library(fastDummies)
library(BAS)
library(rjags)
library(jagsUI)

#### LETTURA DATI####
ames <- read.csv("ameshouse.csv", sep=";")
house <- ames[which(complete.cases(ames)),] #pulizia missing values




#### PRIMO STUDIO COVARIATES####
#Formattazione dati per analisi
ames.numeric <- dplyr::select_if(house, is.numeric) #selezione variabili solo numeriche
x <- as.matrix(ames.numeric[,1:36])
x <- scale(x)
y <- ames.numeric$SalePrice
y <- log(y) 

#Istogramma e normale con media di y e sd di y
hist(y, prob=T, breaks=20)
curve(dnorm(x,mean=mean(y),sd=sd(y),log=FALSE),add=T, col='red', lwd=2)

#Prima lettura dei dati per possibile eliminazione di alcune covariates poco influenti
summary(x)
dataframe <- data.frame( x=as.matrix(x))
par(mfrow=c(3,1))
boxplot(dataframe[1:12])
boxplot(dataframe[13:25])
boxplot(dataframe[26:36])

#Adesso scriviamo il nuovo modello inserendo solo le variabili numeriche piu' rilevanti
x <- x[,-c(15, 10 ,18, 21, 22, 24, 26, 28, 30:36)]
summary(x)




#### ANALISI CORRELAZIONE TRA LE 21 VARIABILI PER ULTERIORE PULIZIA ####
cormat <- round(cor(x),2)
head(cormat)

library(reshape2)
melted_cormat <- melt(cormat)

# Get lower triangle of the correlation matrix
get_lower_tri<-function(cormat){
  cormat[upper.tri(cormat)] <- NA
  return(cormat)
}
# Get upper triangle of the correlation matrix
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}
reorder_cormat <- function(cormat){
  # Use correlation between variables as distance
  dd <- as.dist((1-cormat)/2)
  hc <- hclust(dd)
  cormat <-cormat[hc$order, hc$order]
}
# Reorder the correlation matrix
cormat <- reorder_cormat(cormat)
upper_tri <- get_upper_tri(cormat)
# Melt the correlation matrix
melted_cormat <- melt(upper_tri, na.rm = TRUE)
library(ggplot2)
# Create a ggheatmap
ggheatmap <- ggplot(melted_cormat, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Pearson\nCorrelation") +
  theme_minimal()+ # minimal theme
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))+
  coord_fixed()
# Print the heatmap
ggheatmap + 
  geom_text(aes(Var2, Var1, label = value), color = "black", size = 4) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    legend.justification = c(1, 0),
    legend.position = c(0.6, 0.7),
    legend.direction = "horizontal")+
  guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                               title.position = "top", title.hjust = 0.5))





#### DATA FRAME SCREMATO ####
x <- x[,-c(12, 18, 19)] #rimozione di ulteriori 3 covariates altamente correlate tra loro
numeric <- data.frame( x=as.matrix(x))

n <- nrow(numeric)
par(mfrow=c(3,6), mar=c(4.2,2,2,1.2))
for (j in 1:18){ plot(numeric[,j], y, xlab=names(numeric)[j],
                     pch=19, col="sienna4", xaxt="n", yaxt="n")
abline(lm(y~numeric[,j]), col="green")}


#Analisi ulteriori covariates categoriche
ames.factor <- dplyr::select_if(house, is.factor) #selezione delle sole variabili categoriche
summary(ames.factor) #chiarire  quali come sono distribuite le dummies
dummies <- fastDummies::dummy_cols(ames.factor,remove_first_dummy = TRUE, remove_selected_columns = TRUE) #trasformazione in dummy

dummies.pulito <- dummies[,-c(2, 12, 49, 51, 55, 56, 57, 73, 76, 77, 81, 82, 84, 98, 118, 122, 124, 129, 135, 146, 148, 151, 155, 156, 165, 180, 194, 196 )]

#Unione dei data frame numerico+dummies
completa <- data.frame(y=numeric , x=dummies.pulito)
#DATAFRAME PULITO DA MISSING VALUES E COVARIATES POCO INFLUENTI





#### EFFETTUARE CONFRONTO DIC CON I TRE METODI BIC, AIC E ZG PER VEDERE QUAL'E' IL PIU' EFFICACE####

ames.comparazione <- completa
#creare un nuovo modello con le covariates risultate dalla prima scrematura e poi fare i tre metodi per calcolare il DIC


# --BIC--
ames.BIC = bas.lm(y ~ ., data = ames.comparazione,
                      prior = "BIC", method = "MCMC", MCMC.iterations = 10000 ,  
                      modelprior = uniform())
summary(ames.BIC)
best.BIC = which.max(ames.BIC$logmarg)
bestmodelBIC = ames.BIC$which[[best.BIC]]
bestmodelBIC

X.BIC <- ames.comparazione[,c(bestmodelBIC)]

ames.BIC$logmarg[best.BIC] #Troviamo il logmarg del modello BIC 


# --AIC--
ames.AIC = bas.lm(y ~ ., data = ames.comparazione,
                        prior = "AIC", method = "MCMC", MCMC.iterations = 10000 ,  
                        modelprior = uniform())
summary(ames.AIC)
best.AIC = which.max(ames.AIC$logmarg)
bestmodelAIC = ames.AIC$which[[best.AIC]]
bestmodelAIC

X.AIC <- ames.comparazione[,c(bestmodelAIC)] #Troviamo il logmarg del modello AIC 


# --ZG-prior--
ames.ZG = bas.lm(y ~ ., data = ames.comparazione,
                       prior = "g-prior", method = "MCMC", MCMC.iterations = 10000 , 
                       modelprior = uniform())
summary(ames.ZG)
best.ZG = which.max(ames.ZG$logmarg)
bestmodelZG = ames.ZG$which[[best.ZG]]
bestmodelZG

X.ZG <- ames.comparazione[,c(bestmodelZG)] #Troviamo il logmarg del modello ZS 





#### USO DI JAGS PER DEFINIRE IL VALORE DIC TRAMITE METODI BIC, AIC, ZS####

model.string <- "model{
    # prior
    c2 <- n
    # prior means
    for (j in 1:P){ mu.beta[j] <- 0.0 }
   
    # calculation of xtx
    for (i in 1:P){ for (j in 1:P){
       inverse.V[i,j] <- inprod( x[,i] , x[,j] )
    }}
    for(i in 1:P){ for (j in 1:P){
      prior.T[i,j] <- inverse.V[i,j] * tau /c2
    }}
   
    # likelihood
   
        for (i in 1:n){
            y[i] ~ dnorm( mu[i], tau ) # stochastic componenent
            mu[i] <- inprod( beta[], x[i,] )
        }
       
    # prior distributions
    beta[1:P] ~ dmnorm( mu.beta[], prior.T[,] )
    tau    ~ dgamma( 0.01, 0.01 )
    s2    <- 1/tau
     s <-sqrt(s2)
}"


#DIC - AIC
library(jagsUI)

model <- textConnection(model.string)

x <- model.matrix(~., X.AIC)
n <- nrow(x)
p <- ncol(x)

jags.m <- jagsUI::jags(model,
                       inits=NULL,
                       data = list('n' = n,
                                   'P' = p,
                                   'x' = x, 
                                   'y' = y),
                       n.chains = 2,
                       n.iter = 5000,
                       n.burnin = 2000,
                       n.adapt = 1000,
                       parameters.to.save="beta")
jags.m 
summary(jags.m)
#VALORE DIC~1021


#DIC - ZG
library(jagsUI)

model <- textConnection(model.string)

x <- model.matrix(~., X.ZG)
n <- nrow(x)
p <- ncol(x)

jags.m <- jagsUI::jags(model,
                       inits=NULL,
                       data = list('n' = n,
                                   'P' = p,
                                   'x' = x, 
                                   'y' = y),
                       n.chains = 2,
                       n.iter = 5000,
                       n.burnin = 2000,
                       n.adapt = 1000,
                       parameters.to.save="beta")
jags.m 
summary(jags.m)
#VALORE DIC~961


#DIC - BIC
library(jagsUI)

model <- textConnection(model.string)

x <- model.matrix(~., X.BIC)
n <- nrow(x)
p <- ncol(x)

jags.m <- jagsUI::jags(model,
                       inits=NULL,
                       data = list('n' = n,
                                   'P' = p,
                                   'x' = x, 
                                   'y' = y),
                       n.chains = 2,
                       n.iter = 5000,
                       n.burnin = 2000,
                       n.adapt = 1000,
                       parameters.to.save="beta")
jags.m 
summary(jags.m)
#VALORE DIC~948

#SCELTA CON CRITERIO BIC
#Diagnostic check
par(mfrow=c(1,1))
plot(jags.m) # traceplots e densita'
#Finite comparazioni del DIC






#### RIDUZIONE NUMERO VARIABILI CON METODO BIC ####

library(BAS)

#Prima iterazione

ames.BIC1 = bas.lm(y ~ ., data = X.BIC,
                      prior = "BIC", method = "MCMC", MCMC.iterations = 100000 ,  # We only fit 1 model
                      modelprior = uniform())
summary(ames.BIC1)
best.BIC1 = which.max(ames.BIC1$logmarg)
bestmodelBIC1 = ames.BIC1$which[[best.BIC1]]
bestmodelBIC1

X.BIC1 <- X.BIC[,c(bestmodelBIC1)]

#Seconda iterazione

ames.BIC2 = bas.lm(y ~ ., data = X.BIC1,
                   prior = "BIC", method = "MCMC", MCMC.iterations = 100000 ,  # We only fit 1 model
                   modelprior = uniform())
summary(ames.BIC2)
best.BIC2 = which.max(ames.BIC2$logmarg)
bestmodelBIC2 = ames.BIC2$which[[best.BIC2]]
bestmodelBIC2

X.BIC2 <- X.BIC1[,c(bestmodelBIC2)]
ames.BIC2$logmarg[best.BIC2] #Troviamo il logmarg del modello BIC2 (~-673)

#riordine variabili in funzione della probabilita' marginale di inclusione
ordinato = order(ames.BIC2[["probne0"]])
#scelta delle 20 covariates piu' influenti nella costruzione del modello
k <- ncol(X.BIC2)
ames.primeventiBIC = ordinato[(k-19):k]

ames.finale <- X.BIC2[,c(ames.primeventiBIC-1)]


#BIC FINALE PER STIMA LOGMARG CON LE 20 COVARIATES MIGLIORI
ames.BICfin = bas.lm(y ~ ., data = ames.finale,
                   prior = "BIC", method = "MCMC", MCMC.iterations = 100000 ,  # We only fit 1 model
                   modelprior = uniform())
summary(ames.BICfin)
best.BICfin = which.max(ames.BICfin$logmarg)
ames.BICfin$logmarg[best.BICfin] #~-690
#MODELLO FINALE TROVATO, adesso possiamo procedere con il calcolo dei vari coefficenti 
x <- model.matrix(~., ames.finale)





#### DIC DEL MODELLO FINALE CON SOLE 20 VARIABILI ####
library(jagsUI)

model <- textConnection(model.string)

x <- model.matrix(~., ames.finale)
n <- nrow(x)
p <- ncol(x)

jags.m <- jagsUI::jags(model,
                       inits=NULL,
                       data = list('n' = n,
                                   'P' = p,
                                   'x' = x, 
                                   'y' = y),
                       n.chains = 2,
                       n.iter = 5000,
                       n.burnin = 2000,
                       n.adapt = 1000,
                       parameters.to.save="beta")
jags.m 
summary(jags.m)
#VALORE DIC~915





#### MODELLO LINEARE E COEFFICIENTI ####
#Modello lineare
x.lin <- as.matrix(ames.finale)
linear <- lm(y~x.lin)
summary(linear)
X <- x.lin

pllot(linear)

#COEFFICIENTI
ames.coef = coef(ames.BICfin)
#Ricavo intervalli del Credible Interval
out = confint(ames.coef)[, 1:2]

#Combinazione risultati e creazione di un sommario
coef.BIC = cbind(ames.coef$postmean, ames.coef$postsd, out)
names = c("post mean", "post sd", colnames(out))
colnames(coef.BIC) = names
coef.BIC

#Plot dei beta
par(mfrow = c(2, 2), col.lab = "darkgrey", col.axis = "darkgrey", col = "darkgrey") 
plot(ames.coef, subset = 2:5, ask = F)
par(mfrow = c(2, 2), col.lab = "darkgrey", col.axis = "darkgrey", col = "darkgrey") 
plot(ames.coef, subset = 18:21, ask = F)
confint(ames.coef, parm = 2:21)




#### JAGS NNIG #####
library(rjags)
#Modello con Normal-Normal InverseGamma
n <- nrow(X)
p <- ncol(X)
model_string.AMES.NNIG <- "model{

  # Likelihood
  for(i in 1:n){
    Y[i]   ~ dnorm(mu[i],inv.var)
    mu[i] <- alpha + inprod(X[i,],beta[])
  }

  # Prior for beta
  for(j in 1:p){
    beta[j] ~ dnorm(0,inv.var.b)
  }

  # Prior for the inverse variance
  inv.var   ~ dgamma(0.01, 0.01)
  inv.var.b ~ dgamma(0.01, 0.01)
  alpha     ~ dnorm(0, 0.01)

}"


model_AMES.NNIG <- jags.model(textConnection(model_string.AMES.NNIG), 
                       data = list(Y=y,n=n,p=p,X=X), n.chains=2)

update(model_AMES.NNIG, 10000)

AMES.NNIG <- coda.samples(model_AMES.NNIG, 
                      variable.names=c("beta","alpha"), 
                      n.iter=20000,)

summary(AMES.NNIG)
plot(AMES.NNIG)




#### TEST SENSITIVITA' ####
#Cambio (a,b) Gamma prior
model_string.AMES.NNIG2 <- "model{

  # Likelihood
  for(i in 1:n){
    Y[i]   ~ dnorm(mu[i],inv.var)
    mu[i] <- alpha + inprod(X[i,],beta[])
  }

  # Prior for beta
  for(j in 1:p){
    beta[j] ~ dnorm(0,inv.var.b)
  }

  # Prior for the inverse variance
  inv.var   ~ dgamma(10, 10)
  inv.var.b ~ dgamma(10, 10)
  alpha     ~ dnorm(0, 0.01)

}"


model_AMES.NNIG2 <- jags.model(textConnection(model_string.AMES.NNIG2), 
                              data = list(Y=y,n=n,p=p,X=X), n.chains=2)

update(model_AMES.NNIG2, 10000)

AMES.NNIG2 <- coda.samples(model_AMES.NNIG2, 
                      variable.names=c("beta","alpha"), 
                      n.iter=20000)

summary(AMES.NNIG2)
plot(AMES.NNIG2)

#Cambio numero iterazioni/burn in
model_string.AMES.NNIG3 <- "model{

  # Likelihood
  for(i in 1:n){
    Y[i]   ~ dnorm(mu[i],inv.var)
    mu[i] <- alpha + inprod(X[i,],beta[])
  }

  # Prior for beta
  for(j in 1:p){
    beta[j] ~ dnorm(0,inv.var.b)
  }

  # Prior for the inverse variance
  inv.var   ~ dgamma(0.01, 0.01)
  inv.var.b ~ dgamma(0.01, 0.01)
  alpha     ~ dnorm(0, 0.01)

}"


model_AMES.NNIG3 <- jags.model(textConnection(model_string.AMES.NNIG3), 
                               data = list(Y=y,n=n,p=p,X=X))

update(model_AMES.NNIG3, 100) #Burnin 100

AMES.NNIG3 <- coda.samples(model_AMES.NNIG3, 
                           variable.names=c("beta","alpha"), 
                           n.iter=10000)
summary(AMES.NNIG3)
plot(AMES.NNIG3)

#CONFRONTO SENSITIVITIES
#Prime 12 beta
s1 <- AMES.NNIG[[1]]
s2 <- AMES.NNIG2[[1]]
s3 <- AMES.NNIG3[[1]]
par(mfrow=c(3,4))
for(index in 1:12){
  d1 <- density(s1[,index])
  d2 <- density(s2[,index])
  d3 <- density(s3[,index])
  mx <- max(d1$y,d2$y,d3$y)
  plot(d1,ylim=c(0,mx),xlab="Beta",ylab="Posterior density",main=colnames(x)[index],col="darkgreen")
  lines(d2,col="blue")
  lines(d3,col="red")
}
#Ultime 8 beta
for(index in 13:21){
  d1 <- density(s1[,index])
  d2 <- density(s2[,index])
  d3 <- density(s3[,index])
  mx <- max(d1$y,d2$y,d3$y)
  plot(d1,ylim=c(0,mx),xlab="Beta",ylab="Posterior density",main=colnames(x)[index],col="darkgreen")
  lines(d2,col="blue")
  lines(d3,col="red")
}
#verde standard, blu diversa a,b della Gamma prior, rosso diverse iterazioni/burnin


#Confronto alpha traceplots 
plot(AMES.NNIG[,c('alpha')], main="alpha NNIG")
plot(AMES.NNIG2[,c('alpha')], main="alpha NNIG2")
plot(AMES.NNIG3[,c('alpha')], main="alpha NNIG3")

#Confronto beta 8 traceplots 
plot(AMES.NNIG[,c('beta[8]')], main="Beta[8]")
plot(AMES.NNIG2[,c('beta[8]')], main="Beta[8],2")
plot(AMES.NNIG3[,c('beta[8]')], main="Beta[8],3")


#### JAGS LASSO PRIOR ####
library(rjags)
model_string.AMES.LASSO <- "model{

  # Likelihood
  for(i in 1:n){
    Y[i]   ~ dnorm(mu[i],inv.var)
    mu[i] <- alpha + inprod(X[i,],beta[])
  }

  # Prior for beta
  for(j in 1:p){
    beta[j] ~ ddexp(0,inv.var.b)
  }

  # Prior for the inverse variance
  inv.var   ~ dgamma(0.01, 0.01)
  inv.var.b ~ dgamma(0.01, 0.01)
  alpha     ~ dnorm(0, 0.01)

}"

DATA = list(Y=y,n=n,p=p,X=X)
model_AMES.LASSO <- jags.model(textConnection(model_string.AMES.LASSO), data =DATA, n.chains=2 ) 

update(model_AMES.LASSO, 2000)
AMES.LASSO <- coda.samples(model_AMES.LASSO,variable.names=c("beta","alpha"),n.iter=10000)
summary(AMES.LASSO)
plot(AMES.LASSO)

#VISUALIZZAZIONE BETA
#Prime 12 beta
s1 <- AMES.NNIG[[1]]
s2 <- AMES.LASSO[[1]]
par(mfrow=c(3,4))
for(index in 1:12){
  d1 <- density(s1[,index])
  d2 <- density(s2[,index])
  mx <- max(d1$y,d2$y)
  plot(d1,ylim=c(0,mx),xlab="Beta",ylab="Posterior density",main=colnames(x)[index],col="blue")
  lines(d2,col="red")
}


#Ultime 8 beta
for(index in 13:21){
  d1 <- density(s1[,index])
  d2 <- density(s2[,index])
  mx <- max(d1$y,d2$y)
  plot(d1,ylim=c(0,mx),xlab="Beta",ylab="Posterior density",main=colnames(x)[index],col="blue")
  lines(d2,col="red")
}
#blu NNIG, rosso LASSO


#alpha
plot(AMES.NNIG[,c('alpha')], main="alpha NNIG")
plot(AMES.LASSO[,c('alpha')], main="alpha LASSO")

#beta
plot(AMES.NNIG[,c('beta[12]')], main="beta[12], NNIG")
plot(AMES.NNIG[,c('beta[13]')], main="beta[13], NNIG")
plot(AMES.NNIG[,c('beta[14]')], main="beta[14], NNIG")

plot(AMES.LASSO[,c('beta[12]')], main="beta[12], LASSO")
plot(AMES.LASSO[,c('beta[13]')], main="beta[13], LASSO")
plot(AMES.LASSO[,c('beta[14]')], main="beta[14], LASSO")



#### PREDICTION ####
set.seed(2020)
sample <- sample(1:nrow(X), 15)
xreg <- X[-sample,]
xreal <- X[sample,]

nreg <- nrow(xreg)
p <- ncol(xreg)

yreg <- y[-sample] #response delle prime 568 variabili
yreal <- y[sample] #response y reale dei dati
nreal <- nrow(xreal)

#Modello predictive
model_string.AMES.PREDICTIVE <- "model{

  # Likelihood
  for(i in 1:nreg){
    Y[i]   ~ dnorm(mu[i],inv.var)
    mu[i] <- alpha + inprod(xreg[i,],beta[])
  }

  # Prior for beta
  for(j in 1:p){
    beta[j] ~ dnorm(0,inv.var.b)
  }

  # Prior for the inverse variance
  inv.var   ~ dgamma(0.01, 0.01)
  inv.var.b ~ dgamma(0.01, 0.01)
  alpha     ~ dnorm(0, 0.01)

  # Predictive
  for(k in 1:nreal){
    Yp[k]   ~ dnorm(mup[k],inv.var)
    mup[k] <- alpha + inprod(xreal[k,],beta[])
  }

}"


model_AMES.PREDICTIVE <- jags.model(textConnection(model_string.AMES.PREDICTIVE), 
                              data = list(Y=yreg, nreg=nreg, nreal=nreal, p=p, xreal=xreal, xreg=xreg), n.chains=3)

update(model_AMES.PREDICTIVE, 5000)

AMES.PREDICTIVE <- coda.samples(model_AMES.PREDICTIVE, 
                          variable.names=c("Yp"), 
                          n.iter=10000)

su <- summary(AMES.PREDICTIVE) #sommario modello predictive
su
quantiles<-su$quantiles
quantile1<-su$quantiles[1:15, 1]
quantile2<-su$quantiles[1:15, 5]

#chiedere alice
(ypredicted <- AMES.PREDICTIVE[,c('Yp[1]')][[1]])
d <- density(ypredicted)
plot(d, main="Predictive1")
legend("topright",legend=c("Y Real", "Y Predicted", "quantiles"),col=c("red", "blue", "darkgreen"), lty=3, cex=0.8)
abline(v=yreal[1],lty=3,lwd=1, col="red") #Y REAL
abline(v=mean(ypredicted),lty=3,lwd=1, col="blue") #Y PREDICTED


abline(v=quantile1[1],lty=2,lwd=1, col="darkgreen") #INTERVALLO DI CONFIDENZA
abline(v=quantile2[1],lty=2,lwd=1, col="darkgreen") #INTERVALLO DI CONFIDENZA


par(mfrow=c(1,1))
lmpred = su$statistics[1:(nrow(x)-npredizioni),1]
ypred = as.matrix(lmpred) #ultime 15 response (y) del nostro modello predictive
#autocorrelation factor
acf(ypred,lag.max=50, main="ACF")


#Confronto Real vs Prediction
plot(yreal,col="red",ylim=c(11.5,13),main = "Real vs prediction",ylab="Saleprice",lwd=2)
lines(ypred,type="p", col="blue",lwd=2)
legend("topright",legend=c("Real", "Predicted"),col=c("red", "blue"), lty=3, cex=0.8)


#MSE
yest = ypred-yreal
n = length(yest)
MSE = 1/ (n - 2) * sum((yest ^ 2))
MSE


#PLOT CURVE REALI VS PRED
hist(ypred, prob=T, breaks=10)
curve(dnorm(x,mean=mean(yreal),sd=sd(yreal),log=FALSE),add=T, col='blue', lwd=1)
curve(dnorm(x,mean=mean(ypred),sd=sd(ypred),log=FALSE),add=T, col='red', lwd=1)
legend("topright",legend=c("Real", "Plug-in"),col=c("red", "blue"), lty=1, cex=0.8)


