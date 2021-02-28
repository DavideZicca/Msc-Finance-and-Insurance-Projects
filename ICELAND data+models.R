# Packages to install on first launch
#install.packages("demography")
#install.packages("StMoMo")
#install.packages("MTS")
#install.packages("devtools")
#install_github("robjhyndman/demography")

library(demography)
library("StMoMo")
library("MTS")
library(devtools)
library(fanplot)
library(xlsx) #it allows to read directly an xlsx file
library(writexl) #it allows to write directly an xlsx file

##### Data setting
a.min <- 20
a.max <- 99
omega <-121
A.fit <-c(a.min:a.max)
y.fit.min <- 1978
y.fit.max <- 2018
Y.fit <- c(y.fit.min:y.fit.max)
y.pred <- 100
Y.pred <- c((y.fit.max+1):(y.fit.max+y.pred))
n.sim <- 1000
n.boot <- 20  #you can also use 50 bootstrap but it'll require a lot of time
int_vect <- c(rep(0.02,omega-a.min)) # vector of interest rates (flat)
v_vect <- cumprod((1+int_vect)^-1) # vector of discount factors
#####


####################################
################### Lee-Carter Model
####################################
ISLdata <- hmd.mx(country = "ISL", username = paste0("#######"),
                  password = "#######", label = "Iceland") #create HMD profile
ISLmStMoMo <- StMoMoData(ISLdata, series = "male")
ISLfStMoMo <- StMoMoData(ISLdata, series = "female")
ISLmRates <- ISLmStMoMo$Dxt/ISLmStMoMo$Ext
ISLmRates <- ISLmRates[A.fit+1,tail(ISLmStMoMo$years+1,length(Y.fit))-ISLmStMoMo$years[1]]
ISLfRates <- ISLfStMoMo$Dxt/ISLfStMoMo$Ext
ISLfRates <- ISLfRates[A.fit+1,tail(ISLfStMoMo$years+1,length(Y.fit))-ISLfStMoMo$years[1]]

#ISLfStMoMo and ISLmStMoMo have a central Exposure by default
LCfit_ISLm <- fit(lc(), data = ISLmStMoMo, ages.fit = a.min:a.max,
                  years.fit=y.fit.min:y.fit.max)

LCfit_ISLf <- fit(lc(), data = ISLfStMoMo, ages.fit = a.min:a.max,
                  years.fit=y.fit.min:y.fit.max)
LCfit_ISLm
LCfit_ISLf
#using fit(lc()) without specifing anything, 
#you're using lc(link = c("log"), const = c("sum"))
#use the summary

summary(LCfit_ISLm)
summary(LCfit_ISLf)

BIC(LCfit_ISLm)
BIC(LCfit_ISLf)

plot(LCfit_ISLm)
plot(LCfit_ISLf)

par(mfrow=c(2, 3))
plot(LCfit_ISLm$ax, main="Alpha - Male population")
plot(LCfit_ISLm$bx, main="Beta - Male population")
plot(t(LCfit_ISLm$kt), main="Kappa - Male population")
plot(LCfit_ISLf$ax, main="Alpha - Female population")
plot(LCfit_ISLf$bx, main="Beta - Female population")
plot(t(LCfit_ISLf$kt), main="Kappa - Female population")

LCres_ISLm <- residuals(LCfit_ISLm)        
LCres_ISLf <- residuals(LCfit_ISLf)        
plot(LCres_ISLm)
plot(LCres_ISLf)
plot(LCres_ISLm , type = "signplot")
plot(LCres_ISLf , type = "signplot")
plot(LCres_ISLm , type = "colourmap")
plot(LCres_ISLf , type = "colourmap")

####################################
################### CBD Model with logs (by default cbd uses the logit)
####################################

#central2initial Transform StMoMoData from central to initial exposures. 
#Initial exposures are computed by adding one half of the deaths to the central exposures.
#Usage: central2initial(data)

CBDlogfit_ISLm <- fit(cbd(link = "log"), data = ISLmStMoMo, ages.fit = a.min:a.max,
                   years.fit=y.fit.min:y.fit.max)

CBDlogfit_ISLf <- fit(cbd(link = "log"), data = ISLfStMoMo, ages.fit = a.min:a.max,
                   years.fit=y.fit.min:y.fit.max)

CBDlogfit_ISLm
CBDlogfit_ISLf
summary(CBDlogDfit_ISLm)
summary(CBDlogfit_ISLf)
CBDlogfit_ISLf$loglik
CBDlogfit_ISLm$loglik
BIC(CBDlogfit_ISLm)
BIC(CBDlogfit_ISLf)

par(mar = c(2, 2, 2, 2))
plot(CBDlogfit_ISLm)
plot(CBDlogfit_ISLf)

#plot residuals of CBD model
CBDlogres_ISLm <- residuals(CBDlogfit_ISLm)        
CBDlogres_ISLf <- residuals(CBDlogfit_ISLf)        
plot(CBDlogres_ISLm)
plot(CBDlogres_ISLf)
plot(CBDlogres_ISLm , type = "signplot")
plot(CBDlogres_ISLf , type = "signplot")
plot(CBDlogres_ISLm , type = "colourmap")
plot(CBDlogres_ISLf , type = "colourmap")

####################################
################### Renshaw and Haberman (Lee-Carter with cohorts) mortality model using logs
####################################
RHfit_ISLm <- fit(rh(), data = ISLmStMoMo, ages.fit = A.fit, years.fit = Y.fit )

RHfit_ISLf <- fit(rh(), data = ISLfStMoMo, ages.fit = A.fit, years.fit = Y.fit)

RHfit_ISLm
RHfit_ISLf

summary(RHfit_ISLm)
summary(RHfit_ISLf)
RHfit_ISLf$loglik
RHfit_ISLm$loglik
BIC(RHfit_ISLm)
BIC(RHfit_ISLf)


par(mar = c(2, 2, 2, 2))
plot(RHfit_ISLm)
plot(RHfit_ISLf)


#plot residuals of RH model
RHres_ISLm <- residuals(RHfit_ISLm)        
RHres_ISLf <- residuals(RHfit_ISLf)        
plot(RHres_ISLm)
plot(RHres_ISLf)
plot(RHres_ISLm , type = "signplot")
plot(RHres_ISLf , type = "signplot")
plot(RHres_ISLm , type = "colourmap")
plot(RHres_ISLf , type = "colourmap")

####################################
################### Comparing models using the weights and as
################### link function the logit one
###################################
ISLmIniData <- central2initial(ISLmStMoMo)
ISLfIniData <- central2initial(ISLfStMoMo)
wxt= genWeightMat(ages = A.fit, years = Y.fit, clip = 3 )
ISLmIniData

#Male LC model with logit
LClogit <- lc(link = "logit")
LClogit$gnmFormula
LCfitlogit_ISLm<-fit(LClogit, data = ISLmIniData, ages.fit = A.fit, 
                     years.fit = Y.fit, wxt = wxt)
head(LCfitlogit_ISLm)
summary(LCfitlogit_ISLm)
BIC(LCfitlogit_ISLm)

par(mar = c(2, 2, 2, 2))
plot(LCfitlogit_ISLm)

par(mfrow=c(2, 3))
plot(LCfitlogit_ISLm$ax, main="Alpha - Male population")
plot(LCfitlogit_ISLm$bx, main="Beta - Male population")
plot(t(LCfitlogit_ISLm$kt), main="Kappa - Male population")


LCresLOGIT_ISLm <- residuals(LCfitlogit_ISLm)
plot(LCresLOGIT_ISLm)
plot(LCresLOGIT_ISLm , type = "signplot")
plot(LCresLOGIT_ISLm , type = "colourmap")

#Female LC model with logit

LCfitlogit_ISLf<-fit(LClogit, data = ISLfIniData, ages.fit = A.fit, 
                     years.fit = Y.fit, wxt = wxt)

summary(LCfitlogit_ISLf)
BIC(LCfitlogit_ISLf)

par(mar = c(2, 2, 2, 2))
plot(LCfitlogit_ISLf)

plot(LCfitlogit_ISLf$ax, main="Alpha - Female population")
plot(LCfitlogit_ISLf$bx, main="Beta - Female population")
plot(t(LCfitlogit_ISLf$kt), main="Kappa - Female population")

LCresLOGIT_ISLf <- residuals(LCfitlogit_ISLf) 
plot(LCresLOGIT_ISLf)
plot(LCresLOGIT_ISLf , type = "signplot")
plot(LCresLOGIT_ISLf , type = "colourmap")

#Male CBD model with logit
CBDlogit <- cbd(link = "logit")

CBDfitlogit_ISLm<-fit(CBDlogit, data = ISLmIniData, ages.fit = A.fit, 
                     years.fit = Y.fit, wxt = wxt)
CBDfitlogit_ISLm
summary(CBDfitlogit_ISLm)
BIC(CBDfitlogit_ISLm)

par(mar = c(2, 2, 2, 2))
plot(CBDfitlogit_ISLm)

par(mfrow=c(1, 2))
plot(CBDfitlogit_ISLm$bx, main="Beta - Male population")
plot(t(CBDfitlogit_ISLm$kt), main="Kappa - Male population")


CBDresLOGIT_ISLm <- residuals(CBDfitlogit_ISLm)
plot(CBDresLOGIT_ISLm)
par(mfrow=c(1, 1))

plot(CBDresLOGIT_ISLm , type = "signplot")
plot(CBDresLOGIT_ISLm , type = "colourmap")


#Female CBD model with logit

CBDfitlogit_ISLf<-fit(CBDlogit, data = ISLfIniData, ages.fit = A.fit, 
                      years.fit = Y.fit, wxt = wxt)

summary(CBDfitlogit_ISLf)
BIC(CBDfitlogit_ISLf)

par(mar = c(2, 2, 2, 2))
plot(CBDfitlogit_ISLf)

par(mfrow=c(1, 2))
plot(CBDfitlogit_ISLf$bx, main="Beta - Female population")
plot(t(CBDfitlogit_ISLf$kt), main="Kappa - Female population")

CBDresLOGIT_ISLf <- residuals(CBDfitlogit_ISLf)
plot(CBDresLOGIT_ISLf)
par(mfrow=c(1,1))
plot(CBDresLOGIT_ISLf , type = "signplot")
plot(CBDresLOGIT_ISLf , type = "colourmap")

#Male RH model with logit

RHlogit <- rh(link = "logit", cohortAgeFun = "1")

RHfitlogit_ISLm<-fit(RHlogit, data = ISLmIniData, ages.fit = A.fit, 
                      years.fit = Y.fit, wxt = wxt)

RHfitlogit_ISLm
summary(RHfitlogit_ISLm)
BIC(RHfitlogit_ISLm)

par(mar = c(2, 2, 2, 2))
plot(RHfitlogit_ISLm)

par(mfrow=c(2, 3))
plot(RHfitlogit_ISLm$ax, main="Alpha - Male population")
plot(RHfitlogit_ISLm$bx, main="Beta - Male population")
plot(t(RHfitlogit_ISLm$kt), main="Kappa - Male population")


RHresLOGIT_ISLm <- residuals(RHfitlogit_ISLm)
plot(RHresLOGIT_ISLm)
plot(RHresLOGIT_ISLm , type = "signplot")
plot(RHresLOGIT_ISLm , type = "colourmap")


#Female RH model with logit

RHfitlogit_ISLf<-fit(RHlogit, data = ISLfIniData, ages.fit = A.fit, 
                     years.fit = Y.fit, wxt = wxt)

RHfitlogit_ISLf
summary(RHfitlogit_ISLf)
BIC(RHfitlogit_ISLf)

par(mar = c(2, 2, 2, 2))
plot(RHfitlogit_ISLf)

par(mfrow=c(2, 3))
plot(RHfitlogit_ISLf$ax, main="Alpha - Female population")
plot(RHfitlogit_ISLf$bx, main="Beta - Female population")
plot(t(RHfitlogit_ISLf$kt), main="Kappa - Female population")


RHresLOGIT_ISLf <- residuals(RHfitlogit_ISLf)
plot(RHresLOGIT_ISLf)
plot(RHresLOGIT_ISLf , type = "signplot")
plot(RHresLOGIT_ISLf , type = "colourmap")


#Male M7 model with logit

M7logit <- m7(link = "logit")
M7fitlogitm <- fit(M7logit, data = ISLmIniData, ages.fit = A.fit, 
              years.fit = Y.fit, wxt = wxt)

M7fitlogitm
summary(M7fitlogitm)
BIC(M7fitlogitm)

par(mar = c(2, 2, 2, 2))
plot(M7fitlogitm)

par(mfrow=c(2, 2))
plot(M7fitlogitm$bx, main="Beta - Male population")
plot(t(M7fitlogitm$kt), main="Kappa - Male population")


M7reslogitm <- residuals(M7fitlogitm)
plot(M7reslogitm)
plot(M7reslogitm , type = "signplot")
plot(M7reslogitm , type = "colourmap")


#Female M7 model with logit

M7logit <- m7(link = "logit")
M7fitlogitf <- fit(M7logit, data = ISLfIniData, ages.fit = A.fit, 
                   years.fit = Y.fit, wxt = wxt)

M7fitlogitf
summary(M7fitlogitf)
BIC(M7fitlogitf)

par(mar = c(2, 2, 2, 2))
plot(M7fitlogitf)

par(mfrow=c(1, 2))
plot(M7fitlogitf$bx, main="Beta - Female population")
plot(t(M7fitlogitf$kt), main="Kappa - Female population")

M7reslogitf <- residuals(M7fitlogitf)
plot(M7reslogitf)
par(mfrow=c(1, 1))
plot(M7reslogitf , type = "signplot")
plot(M7reslogitf , type = "colourmap")

#Male APC model with logit

APClogit <- apc(link = "logit")
APCfitlogitm <- fit(APClogit, data = ISLmIniData, ages.fit = A.fit, 
                   years.fit = Y.fit, wxt = wxt)

APCfitlogitm
summary(APCfitlogitm)
BIC(APCfitlogitm)

par(mar = c(2, 2, 2, 2))
plot(APCfitlogitm)

par(mfrow=c(1, 3))
plot(APCfitlogitm$ax, main="Alpha - Male population")
plot(APCfitlogitm$bx, main="Beta - Male population")
plot(t(APCfitlogitm$kt), main="Kappa - Male population")


APCreslogitm <- residuals(APCfitlogitm)
plot(APCreslogitm)
par(mfrow=c(1, 1))
plot(APCreslogitm , type = "signplot")
plot(APCreslogitm , type = "colourmap")


#Female APC model with logit

APClogit <- apc(link = "logit")
APCfitlogitf <- fit(APClogit, data = ISLfIniData, ages.fit = A.fit, 
                    years.fit = Y.fit, wxt = wxt)

APCfitlogitf
summary(APCfitlogitf)
BIC(APCfitlogitf)

par(mar = c(2, 2, 2, 2))
plot(APCfitlogitf)

par(mfrow=c(1, 3))
plot(APCfitlogitf$ax, main="Alpha - Female population")
plot(APCfitlogitf$bx, main="Beta - Female population")
plot(t(APCfitlogitf$kt), main="Kappa - Female population")


APCreslogitf <- residuals(APCfitlogitf)
plot(APCreslogitf)
par(mfrow=c(1, 1))
plot(APCreslogitf , type = "signplot")
plot(APCreslogitf , type = "colourmap")

#####################################
##### PLAT model
#####################################
f2 <- function(x, ages) mean(ages) - x
constPlat <- function(ax, bx, kt, b0x, gc, wxt, ages){
  nYears <- dim(wxt)[2]
  x <- ages
  t <- 1:nYears
  c <- (1 - tail(ages, 1)):(nYears - ages[1])
  xbar <- mean(x)
  phiReg <- lm(gc ~ 1 + c + I(c ^ 2), na.action = na.omit)
  phi <- coef(phiReg)
  gc <- gc - phi[1] - phi[2] * c - phi[3] * c ^ 2
  kt[2, ] <- kt[2, ] + 2 * phi[3] * t
  kt[1, ] <- kt[1, ] + phi[2] * t + phi[3] * (t ^ 2 - 2 * xbar * t)
  ax <- ax + phi[1] - phi[2] * x + phi[3] * x ^ 2
  ci <- rowMeans(kt, na.rm = TRUE)
  ax <- ax + ci[1] + ci[2] * (xbar - x)
  kt[1, ] <- kt[1, ] - ci[1]
  kt[2, ] <- kt[2, ] - ci[2]
  list(ax = ax, bx = bx, kt = kt, b0x = b0x, gc = gc)
}
PLAT <- StMoMo(link = "logit", staticAgeFun = TRUE,
               periodAgeFun = c("1", f2), cohortAgeFun = "1", constFun = constPlat)

#Male PLAT model with logit
PLATfitlogit_ISLm <- fit(PLAT, data = ISLmIniData, ages.fit = A.fit,years.fit=Y.fit, wxt = wxt)
PLATfitlogit_ISLm

par(mar = c(2, 2, 2, 2))
plot(PLATfitlogit_ISLm, parametricbx = F)

par(mfrow=c(1, 3))
plot(PLATfitlogit_ISLm$ax, main="Alpha - Male population")
plot(PLATfitlogit_ISLm$bx, main="Beta - Male population")
plot(t(PLATfitlogit_ISLm$kt), main="Kappa - Male population")


PLATreslogitm <- residuals(PLATfitlogit_ISLm)
plot(PLATreslogitm)
par(mfrow=c(1, 1))
plot(PLATreslogitm , type = "signplot")
plot(PLATreslogitm , type = "colourmap")

#Female PLAT model with logit
PLATfitlogit_ISLf <- fit(PLAT, data = ISLfIniData, ages.fit = A.fit,years.fit=Y.fit, wxt = wxt)
PLATfitlogit_ISLf
par(mar = c(2, 2, 2, 2))
plot(PLATfitlogit_ISLf, parametricbx = F)

par(mfrow=c(1, 3))
plot(PLATfitlogit_ISLf$ax, main="Alpha - Female population")
plot(PLATfitlogit_ISLf$bx, main="Beta - Female population")
plot(t(PLATfitlogit_ISLf$kt), main="Kappa - Female population")


PLATreslogitf <- residuals(PLATfitlogit_ISLf)
plot(PLATreslogitf)
par(mfrow=c(1, 1))
plot(PLATreslogitf , type = "signplot")
plot(PLATreslogitf , type = "colourmap")


#####################################
######### Models comparison
#####################################

#quantitative comparative using BIC and log-likelihood

#comparing MALE results with logs
loglikemalelogs= c(RHfit_ISLm$loglik,CBDlogfit_ISLm$loglik, LCfit_ISLm$loglik)
BICmalelogs= c(BIC(RHfit_ISLm), BIC(CBDlogfit_ISLm), BIC(LCfit_ISLm))

rowSloglikemale = c("RH", "CBD", "LC")
loglikeResultsmale= data.frame(loglikemalelogs, row.names = rowSloglikemale )
loglikeResultsmale

loglikeResultsmale[which.max(loglikeResultsmale$loglikemalelogs),]

rowSBICmale= c("RH", "CBD", "LC")
BIClogsResultsmale= data.frame(BICmalelogs, row.names = rowSBICmale )
BIClogsResultsmale

BIClogsResultsmale[which.min(BIClogsResultsmale$BICmalelogs),]

#comparing FEMALE results with logs

loglikefemalelogs= c(RHfit_ISLf$loglik,CBDlogfit_ISLf$loglik, LCfit_ISLf$loglik)
BICfemalelogs= c(BIC(RHfit_ISLf), BIC(CBDlogfit_ISLf), BIC(LCfit_ISLf))

rowSloglikefemale = c("RH", "CBD", "LC")
loglikeResultsfemale= data.frame(loglikefemalelogs, row.names = rowSloglikefemale )
loglikeResultsfemale
loglikeResultsfemale[which.max(loglikeResultsfemale$loglikefemalelogs),]

rowSBICfemale= c("RH", "CBD", "LC")
BIClogsResultsfemale= data.frame(BICfemalelogs, row.names = rowSBICfemale )
BIClogsResultsfemale
BIClogsResultsfemale[which.min(BIClogsResultsfemale$BICfemalelogs),]


#comparing MALE results with LOGIT
loglikemaleLOGIT= c(RHfitlogit_ISLm$loglik,CBDfitlogit_ISLm$loglik, 
                   LCfitlogit_ISLm$loglik, APCfitlogitm$loglik, M7fitlogitm$loglik, PLATfitlogit_ISLm$loglik)
BICmaleLOGIT= c(BIC(RHfitlogit_ISLm), BIC(CBDfitlogit_ISLm), BIC(LCfitlogit_ISLm), 
               BIC(APCfitlogitm),BIC(M7fitlogitm), BIC(PLATfitlogit_ISLm))

rowSLOGITlikemale = c("RH", "CBD", "LC", "APC", "M7", "PLAT")
LOGITlikeResultsmale= data.frame(loglikemaleLOGIT, row.names = rowSLOGITlikemale )
LOGITlikeResultsmale

LOGITlikeResultsmale[which.max(LOGITlikeResultsmale$loglikemaleLOGIT),]

rowSBICLOGITmale= c("RH", "CBD", "LC", "APC", "M7", "PLAT")
BICLOGITResultsmale= data.frame(BICmaleLOGIT, row.names = rowSBICLOGITmale )
BICLOGITResultsmale

BICLOGITResultsmale[which.min(BICLOGITResultsmale$BICmaleLOGIT),]
#comparing FEMALE results with LOGIT

loglikefemaleLOGIT= c(RHfitlogit_ISLf$loglik,CBDfitlogit_ISLf$loglik, 
                    LCfitlogit_ISLf$loglik, APCfitlogitf$loglik, 
                    M7fitlogitf$loglik, PLATfitlogit_ISLf$loglik)
BICfemaleLOGIT= c(BIC(RHfitlogit_ISLf), BIC(CBDfitlogit_ISLf), BIC(LCfitlogit_ISLf), 
                BIC(APCfitlogitf),BIC(M7fitlogitf), BIC(PLATfitlogit_ISLf))

rowSLOGITlikefemale = c("RH", "CBD", "LC", "APC", "M7", "PLAT")
LOGITlikeResultsfemale= data.frame(loglikefemaleLOGIT, row.names = rowSLOGITlikefemale )
LOGITlikeResultsfemale

LOGITlikeResultsfemale[which.max(LOGITlikeResultsfemale$loglikefemaleLOGIT),]

rowSBICLOGITfemale= c("RH", "CBD", "LC", "APC", "M7", "PLAT")
BICLOGITResultsfemale= data.frame(BICfemaleLOGIT, row.names = rowSBICLOGITfemale )
BICLOGITResultsfemale
BICLOGITResultsfemale[which.min(BICLOGITResultsfemale$BICfemaleLOGIT),]


#number of parameters in the Logit models

RHfitlogit_ISLm$npar
CBDfitlogit_ISLm$npar
LCfitlogit_ISLm$npar
APCfitlogitm$npar
M7fitlogitm$npar
PLATfitlogit_ISLm$npar

n.parameters=c(RHfitlogit_ISLm$npar,CBDfitlogit_ISLm$npar,LCfitlogit_ISLm$npar,APCfitlogitm$npar,
               M7fitlogitm$npar,PLATfitlogit_ISLm$npar)
modelsparameters=data.frame(n.parameters,row.names = rowSBICLOGITfemale )

RHfitlogit_ISLf$npar
CBDfitlogit_ISLf$npar
LCfitlogit_ISLf$npar
APCfitlogitf$npar
M7fitlogitf$npar
PLATfitlogit_ISLf$npar
