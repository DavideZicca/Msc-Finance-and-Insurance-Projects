source("Functions.r") #uploading functions from an other script in the directory

####################################
################### Lee-Carter model with logit
####################################
################### Forecast
LCfor_ISLm <- forecast(LCfitlogit_ISLm, h=y.pred)
### Alternative istructions
#LCfor_ISLmArima <- forecast(LCfitlogit_ISLm, h=y.pred, kt.method = "iarima",  kt.order = c(1, 1, 2))  # ARIMA order specified
#LCfor_ISLmBArima <- forecast(LCfitlogit_ISLm, h=y.pred, kt.method = "iarima",  kt.order = NULL)  # auto.arima
LCfor_ISLf <- forecast(LCfitlogit_ISLf, h=y.pred)
plot(LCfor_ISLm, only.kt = TRUE)
plot(LCfor_ISLf, only.kt = TRUE)
#####
par(mfrow=c(1, 2), oma=c(2,0,0,0))
plot(c(a.min:a.max),log(LCfor_ISLm$rates[,y.pred]), bty="l", xlab="Ages", ylab="Male death rates (log scale)",
  ylim=c(min(log(LCfor_ISLm$rates[,y.pred])),max(log((LCfitlogit_ISLm$Dxt/LCfitlogit_ISLm$Ext)[,length(Y.fit)]))), type="l")
lines(c(a.min:a.max),log((LCfitlogit_ISLm$Dxt/LCfitlogit_ISLm$Ext)[,length(Y.fit)]), type="p")
plot(c(a.min:a.max),log(LCfor_ISLf$rates[,y.pred]), bty="l", xlab="Ages", ylab="Female death rates (log scale)",
  ylim=c(min(log(LCfor_ISLf$rates[,y.pred])),max(log((LCfitlogit_ISLf$Dxt/LCfitlogit_ISLf$Ext)[,length(Y.fit)]))), type="l")
lines(c(a.min:a.max),log((LCfitlogit_ISLf$Dxt/LCfitlogit_ISLf$Ext)[,length(Y.fit)]), type="p")
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend("bottom", c("last observed", "last projected"), pch=c(1,-1), lty=c(0,1), lwd=c(1,1), col=c(1:1), cex=0.9)
dev.off()

rates_LClogitfit_ISLm <- fitted(LCfitlogit_ISLm, type = "rates")
rates_LClogitfor_ISLm <- LCfor_ISLm$rates
#rates_LClogit_ISLm <- cbind(rates_LClogitfit_ISLm,rates_LClogitfor_ISLm)
q_LC_ISLm <- 1- exp(-rates_LClogit_ISLm)
q_LC_ISLm.ext  <- extrapolation.fit(q_LC_ISLm)

rates_LClogitfit_ISLf <- fitted(LCfitlogit_ISLf, type = "rates")
rates_LClogitfor_ISLf <- LCfor_ISLf$rates
#rates_LC_ISLf <- cbind(rates_LCfit_ISLf,rates_LCfor_ISLf)
q_LC_ISLf <- 1- exp(-rates_LC_ISLf)
q_LC_ISLf.ext  <- extrapolation.fit(q_LC_ISLf)

write.table(q_LC_ISLm.ext,file="q_LC.ISLm.txt",sep=",")
write.table(q_LC_ISLf.ext,file="q_LC.ISLf.txt",sep=",")
q_LC.ISLm<-as.data.frame(q_LC_ISLm.ext)
q_LC.ISLf<-as.data.frame(q_LC_ISLf.ext)

plot(extractCohort(q_LC_ISLm.ext, cohort = 1980), type = "p", log = "y", xlab = "age", ylab = "q(x)",
 main = "Mortality rates for a cohort (log scale)")
lines(extractCohort(q_LC_ISLm, cohort = 1980))
plot(extractCohort(log(q_LC_ISLm.ext/(1-q_LC_ISLm.ext)), cohort = 1980), type = "p", xlab = "age", ylab = "logit(q(x))",
 main = "Logit transform of mortality rates for a cohort")
lines(extractCohort(log(q_LC_ISLm/(1-q_LC_ISLm)), cohort = 1980))
############################

####################################
######## 1 year survival probability, n year survival probability and life expectancy
####################################
p_LC_ISLm.ext <- 1-q_LC_ISLm.ext # 1 year survival probability
p0n_LC_ISLm.ext <- apply(p_LC_ISLm.ext, 2, cumprod) # n year survival probability
ex_LC_ISLm.ext <- life.exp(q_LC_ISLm.ext)

p_LC_ISLf.ext <- 1-q_LC_ISLf.ext # 1 year survival probability
p0n_LC_ISLf.ext <- apply(p_LC_ISLf.ext, 2, cumprod) # n year survival probability
ex_LC_ISLf.ext <- life.exp(q_LC_ISLf.ext)

write.table(ex_LC_ISLm.ext,file="ex_LC.ISLm.txt",sep=",")
write.table(ex_LC_ISLf.ext,file="ex_LC.ISLf.txt",sep=",")
p_LC.ISLm<-as.data.frame(p_LC_ISLm.ext)
p0n_LC.ISLm<-as.data.frame(p0n_LC_ISLm.ext)
ex_LC.ISLm<-as.data.frame(ex_LC_ISLm.ext)
p_LC.ISLf<-as.data.frame(p_LC_ISLf.ext)
p0n_LC.ISLf<-as.data.frame(p0n_LC_ISLf.ext)
ex_LC.ISLf<-as.data.frame(ex_LC_ISLf.ext)
####################################

####################################
##### SIMULATION WITH RANDOM WALK WITH DRIFT
####################################
LCsim_ISLm.mrwd <- simulate(LCfitlogit_ISLm, nsim = n.sim, h=y.pred)
### Alternative istructions
#LCsim_ISLmArima <- forecast(LCfitlogit_ISLm, nsim = n.sim, h=y.pred, kt.method = "iarima",  kt.order = c(1, 1, 2))  # ARIMA order specified
#LCsim_ISLmBArima <- forecast(LCfitlogit_ISLm, nsim = n.sim, h=y.pred, kt.method = "iarima",  kt.order = NULL)  # auto.arima
rates_LC_ISLm.st <- LCsim_ISLm.mrwd$rates
q_LC_ISLm.st <- 1- exp(-rates_LC_ISLm.st)
q_LC_ISLm.st.ext <-  extrapolation.sim(q_LC_ISLm.st)

LCsim_ISLf.mrwd <- simulate(LCfitlogit_ISLf, nsim = n.sim, h=y.pred)
rates_LC_ISLf.st <- LCsim_ISLf.mrwd$rates
q_LC_ISLf.st <- 1- exp(-rates_LC_ISLf.st)
q_LC_ISLf.st.ext <-  extrapolation.sim(q_LC_ISLf.st)

par(mfrow=c(1, 2))
plot(LCfitlogit_ISLm$years, LCfitlogit_ISLm$kt[1, ], xlim = range(LCfitlogit_ISLm$years, LCsim_ISLm.mrwd$kt.s$years),
  ylim = range(LCfitlogit_ISLm$kt, LCsim_ISLm.mrwd$kt.s$sim), type = "l", xlab = "year", ylab = "kt",
  main = "Lee-Carter: Period index (mrwd)")
matlines(LCsim_ISLm.mrwd$kt.s$years, LCsim_ISLm.mrwd$kt.s$sim[1, , ], type = "l", lty = 1)
plot(LCfitlogit_ISLm$years, (LCfitlogit_ISLm$Dxt / LCfitlogit_ISLm$Ext)["65", ], xlim = range(LCfitlogit_ISLm$years, LCsim_ISLm.mrwd$years),
  ylim = range((LCfitlogit_ISLm$Dxt / LCfitlogit_ISLm$Ext)["65", ], LCsim_ISLm.mrwd$rates["65", , ]), type = "l", xlab = "year", ylab = "rate",
  main = "Lee-Carter: Simulated mortality rates at age 65")
matlines(LCsim_ISLm.mrwd$years, LCsim_ISLm.mrwd$rates["65", , ], type = "l", lty = 1)
dev.off()

library("fanplot")
probs = c(2.5, 10, 25, 50, 75, 90, 97.5)
par(mfrow=c(1, 1))
matplot(LCfitlogit_ISLm$years, t(q_LC_ISLm[c("65", "75", "85"),c(1:length(Y.fit)) ]),
 xlim = c(y.fit.min, (y.fit.min+y.pred)), ylim = c(0.001, 0.2), pch = 20, col = "black",
 log = "y", xlab = "year", ylab = "male mortality rate (log scale)")
fan(t(LCsim_ISLm.mrwd$rates["65", , ]), start = y.fit.max+1, probs = probs, n.fan = 4,
 fan.col = colorRampPalette(c("black", "white")), ln = NULL)
fan(t(LCsim_ISLm.mrwd$rates["75", , ]), start = y.fit.max+1, probs = probs, n.fan = 4,
 fan.col = colorRampPalette(c("red", "white")), ln = NULL)
fan(t(LCsim_ISLm.mrwd$rates["85", , ]), start = y.fit.max+1, probs = probs, n.fan = 4,
 fan.col = colorRampPalette(c("blue", "white")), ln = NULL)
text(y.fit.min+4, q_LC_ISLm[c("65", "75", "85"), 30],
 labels = c("x = 65", "x = 75", "x = 85"))
matplot(LCfitlogit_ISLf$years, t(q_LC_ISLf[c("65", "75", "85"),c(1:length(Y.fit)) ]),
 xlim = c(y.fit.min, (y.fit.min+y.pred)), ylim = c(0.001, 0.2), pch = 20, col = "black",
 log = "y", xlab = "year", ylab = "female mortality rate (log scale)")
fan(t(LCsim_ISLf.mrwd$rates["65", , ]), start = y.fit.max+1, probs = probs, n.fan = 4,
 fan.col = colorRampPalette(c("black", "white")), ln = NULL)
fan(t(LCsim_ISLf.mrwd$rates["75", , ]), start = y.fit.max+1, probs = probs, n.fan = 4,
 fan.col = colorRampPalette(c("red", "white")), ln = NULL)
fan(t(LCsim_ISLf.mrwd$rates["85", , ]), start = y.fit.max+1, probs = probs, n.fan = 4,
 fan.col = colorRampPalette(c("blue", "white")), ln = NULL)
text(y.fit.min+4, q_LC_ISLf[c("65", "75", "85"), 30],
 labels = c("x = 65", "x = 75", "x = 85"))
dev.off()
####################################

####################################
####### Bootstrap
####################################
LCboot_ISLm <- bootstrap(LCfitlogit_ISLm, nBoot = n.boot, type = "semiparametric")
plot(LCboot_ISLm, nCol = 3)
LCsim_ISLm.boot <- simulate(LCboot_ISLm, nsim = n.sim/n.boot, h = y.pred)
rates_LC_ISLm.boot.st <- LCsim_ISLm.boot$rates
q_LC_ISLm.boot.st <- 1- exp(-rates_LC_ISLm.boot.st)
q_LC_ISLm.boot.st.ext <-  extrapolation.sim(q_LC_ISLm.boot.st)

LCboot_ISLf <- bootstrap(LCfitlogit_ISLf, nBoot = n.boot, type = "semiparametric")
plot(LCboot_ISLf, nCol = 3)
LCsim_ISLf.boot <- simulate(LCboot_ISLf, nsim = n.sim/n.boot, h = y.pred)
rates_LC_ISLf.boot.st <- LCsim_ISLf.boot$rates
q_LC_ISLf.boot.st <- 1- exp(-rates_LC_ISLf.boot.st)
q_LC_ISLf.boot.st.ext <-  extrapolation.sim(q_LC_ISLf.boot.st)
dev.off()

####################################

####################################
#### Confidence intervals at 90%
####################################
##### 1 year death probability
conf.lev <- 0.9
q_LC_ISLm.q95 <- matrix(0, nrow=(omega-a.min),ncol=y.pred, dimnames=list(c(a.min:(omega-1)), c((tail(Y.fit,1)+1):(tail(Y.fit,1)+y.pred))))
q_LC_ISLm.q05 <- matrix(0, nrow=(omega-a.min),ncol=y.pred, dimnames=list(c(a.min:(omega-1)), c((tail(Y.fit,1)+1):(tail(Y.fit,1)+y.pred))))
q_LC_ISLm.boot.q95 <- matrix(0, nrow=(omega-a.min),ncol=y.pred, dimnames=list(c(a.min:(omega-1)), c((tail(Y.fit,1)+1):(tail(Y.fit,1)+y.pred))))
q_LC_ISLm.boot.q05 <- matrix(0, nrow=(omega-a.min),ncol=y.pred, dimnames=list(c(a.min:(omega-1)), c((tail(Y.fit,1)+1):(tail(Y.fit,1)+y.pred))))
for (j in 1:(omega-a.min)){
  for (k in 1:y.pred){
      q_LC_ISLm.q95[j,k] <- quantile(q_LC_ISLm.st.ext[j,k,], probs=0.5+conf.lev/2)
      q_LC_ISLm.q05[j,k] <- quantile(q_LC_ISLm.st.ext[j,k,], probs=(1-conf.lev)/2)
      q_LC_ISLm.boot.q95[j,k] <- quantile(q_LC_ISLm.boot.st.ext[j,k,], probs=0.5+conf.lev/2)
      q_LC_ISLm.boot.q05[j,k] <- quantile(q_LC_ISLm.boot.st.ext[j,k,], probs=(1-conf.lev)/2)
  }
}

q_LC_ISLf.q95 <- matrix(0, nrow=(omega-a.min),ncol=y.pred, dimnames=list(c(a.min:(omega-1)), c((tail(Y.fit,1)+1):(tail(Y.fit,1)+y.pred))))
q_LC_ISLf.q05 <- matrix(0, nrow=(omega-a.min),ncol=y.pred, dimnames=list(c(a.min:(omega-1)), c((tail(Y.fit,1)+1):(tail(Y.fit,1)+y.pred))))
q_LC_ISLf.boot.q95 <- matrix(0, nrow=(omega-a.min),ncol=y.pred, dimnames=list(c(a.min:(omega-1)), c((tail(Y.fit,1)+1):(tail(Y.fit,1)+y.pred))))
q_LC_ISLf.boot.q05 <- matrix(0, nrow=(omega-a.min),ncol=y.pred, dimnames=list(c(a.min:(omega-1)), c((tail(Y.fit,1)+1):(tail(Y.fit,1)+y.pred))))
for (j in 1:(omega-a.min)){
  for (k in 1:y.pred){
      q_LC_ISLf.q95[j,k] <- quantile(q_LC_ISLf.st.ext[j,k,], probs=0.5+conf.lev/2)
      q_LC_ISLf.q05[j,k] <- quantile(q_LC_ISLf.st.ext[j,k,], probs=(1-conf.lev)/2)
      q_LC_ISLf.boot.q95[j,k] <- quantile(q_LC_ISLf.boot.st.ext[j,k,], probs=0.5+conf.lev/2)
      q_LC_ISLf.boot.q05[j,k] <- quantile(q_LC_ISLf.boot.st.ext[j,k,], probs=(1-conf.lev)/2)
  }
}
#write.table(q_LC_ISLm.q95,file="q_LC.ISLm.q95.txt",sep=",")
#write.table(q_LC_ISLm.q05,file="q_LC.ISLm.q05.txt",sep=",")
#write.table(q_LC_ISLm.boot.q95,file="q_LC.ISLm.boot.q95.txt",sep=",")
#write.table(q_LC_ISLm.boot.q05,file="q_LC.ISLm.boot.q05.txt",sep=",")
#write.table(q_LC_ISLf.q95,file="q_LC.ISLf.q95.txt",sep=",")
#write.table(q_LC_ISLf.q05,file="q_LC.ISLf.q05.txt",sep=",")
#write.table(q_LC_ISLf.boot.q95,file="q_LC.ISLf.boot.q95.txt",sep=",")
#write.table(q_LC_ISLf.boot.q05,file="q_LC.ISLf.boot.q05.txt",sep=",")
q_LC.ISLm.q95<-as.data.frame(q_LC_ISLm.q95)
q_LC.ISLm.q05<-as.data.frame(q_LC_ISLm.q05)
q_LC.ISLm.boot.q95<-as.data.frame(q_LC_ISLm.boot.q95)
q_LC.ISLm.boot.q05<-as.data.frame(q_LC_ISLm.boot.q05)
q_LC.ISLf.q95<-as.data.frame(q_LC_ISLf.q95)
q_LC.ISLf.q05<-as.data.frame(q_LC_ISLf.q05)
q_LC.ISLf.boot.q95<-as.data.frame(q_LC_ISLf.boot.q95)
q_LC.ISLf.boot.q05<-as.data.frame(q_LC_ISLf.boot.q05)
#male
plot(q_LC_ISLm.ext["65",],x=c(Y.fit,Y.pred), type="l", lty=1) # evolution with t (same age)
lines(q_LC_ISLm.q95["65",], x=Y.pred, type="l", lty=2)
lines(q_LC_ISLm.q05["65",], x=Y.pred, type="l", lty=2)
lines(q_LC_ISLm.boot.q95["65",], x=Y.pred, type="l", lty=2, col=2)
lines(q_LC_ISLm.boot.q05["65",], x=Y.pred, type="l", lty=2, col=2)

plot(q_LC_ISLm.ext["20",],x=c(Y.fit,Y.pred), type="l", lty=1) # evolution with t (same age)
lines(q_LC_ISLm.q95["20",], x=Y.pred, type="l", lty=2)
lines(q_LC_ISLm.q05["20",], x=Y.pred, type="l", lty=2)
lines(q_LC_ISLm.boot.q95["20",], x=Y.pred, type="l", lty=2, col=2)
lines(q_LC_ISLm.boot.q05["20",], x=Y.pred, type="l", lty=2, col=2)
#female
plot(q_LC_ISLf.ext["65",],x=c(Y.fit,Y.pred), type="l", lty=1) # evolution with t (same age)
lines(q_LC_ISLf.q95["65",], x=Y.pred, type="l", lty=2)
lines(q_LC_ISLf.q05["65",], x=Y.pred, type="l", lty=2)
lines(q_LC_ISLf.boot.q95["65",], x=Y.pred, type="l", lty=2, col=2)
lines(q_LC_ISLf.boot.q05["65",], x=Y.pred, type="l", lty=2, col=2)

plot(q_LC_ISLf.ext["20",],x=c(Y.fit,Y.pred), type="l", lty=1) # evolution with t (same age)
lines(q_LC_ISLf.q95["20",], x=Y.pred, type="l", lty=2)
lines(q_LC_ISLf.q05["20",], x=Y.pred, type="l", lty=2)
lines(q_LC_ISLf.boot.q95["20",], x=Y.pred, type="l", lty=2, col=2)
lines(q_LC_ISLf.boot.q05["20",], x=Y.pred, type="l", lty=2, col=2)
dev.off()
####################################

##### 1 year survival probability, n year survival probability and life expectancy confidence interval
p_LC_ISLm.q95 <- 1-q_LC_ISLm.q95# 1 year survival probability
p0n_LC_ISLm.q95 <- apply(p_LC_ISLm.q95, 2, cumprod) # n year survival probability
ex_LC_ISLm.q95 <- life.exp(q_LC_ISLm.q95) # life expectancy

p_LC_ISLm.q05 <- 1-q_LC_ISLm.q05# 1 year survival probability
p0n_LC_ISLm.q05 <- apply(p_LC_ISLm.q05, 2, cumprod) # n year survival probability
ex_LC_ISLm.q05 <- life.exp(q_LC_ISLm.q05) # life expectancy

p_LC_ISLm.boot.q95 <- 1-q_LC_ISLm.boot.q95# 1 year survival probability
p0n_LC_ISLm.boot.q95 <- apply(p_LC_ISLm.boot.q95, 2, cumprod) # n year survival probability
ex_LC_ISLm.boot.q95 <- life.exp(q_LC_ISLm.boot.q95) # life expectancy

p_LC_ISLm.boot.q05 <- 1-q_LC_ISLm.boot.q05# 1 year survival probability
p0n_LC_ISLm.boot.q05 <- apply(p_LC_ISLm.boot.q05, 2, cumprod) # n year survival probability
ex_LC_ISLm.boot.q05 <- life.exp(q_LC_ISLm.boot.q05) # life expectancy

p_LC_ISLf.q95 <- 1-q_LC_ISLf.q95# 1 year survival probability
p0n_LC_ISLf.q95 <- apply(p_LC_ISLf.q95, 2, cumprod) # n year survival probability
ex_LC_ISLf.q95 <- life.exp(q_LC_ISLf.q95) # life expectancy

p_LC_ISLf.q05 <- 1-q_LC_ISLf.q05# 1 year survival probability
p0n_LC_ISLf.q05 <- apply(p_LC_ISLf.q05, 2, cumprod) # n year survival probability
ex_LC_ISLf.q05 <- life.exp(q_LC_ISLf.q05) # life expectancy

p_LC_ISLf.boot.q95 <- 1-q_LC_ISLf.boot.q95# 1 year survival probability
p0n_LC_ISLf.boot.q95 <- apply(p_LC_ISLf.boot.q95, 2, cumprod) # n year survival probability
ex_LC_ISLf.boot.q95 <- life.exp(q_LC_ISLf.boot.q95) # life expectancy

p_LC_ISLf.boot.q05 <- 1-q_LC_ISLf.boot.q05# 1 year survival probability
p0n_LC_ISLf.boot.q05 <- apply(p_LC_ISLf.boot.q05, 2, cumprod) # n year survival probability
ex_LC_ISLf.boot.q05 <- life.exp(q_LC_ISLf.boot.q05) # life expectancy

ex_LC.ISLm.q05<-as.data.frame(ex_LC_ISLm.q05)
ex_LC.ISLm.q95<-as.data.frame(ex_LC_ISLm.q95)
ex_LC.ISLm.boot.q05<-as.data.frame(ex_LC_ISLm.boot.q05)
ex_LC.ISLm.boot.q95<-as.data.frame(ex_LC_ISLm.boot.q95)
ex_LC.ISLf.q05<-as.data.frame(ex_LC_ISLf.q05)
ex_LC.ISLf.q95<-as.data.frame(ex_LC_ISLf.q95)
ex_LC.ISLf.boot.q05<-as.data.frame(ex_LC_ISLf.boot.q05)
ex_LC.ISLf.boot.q95<-as.data.frame(ex_LC_ISLf.boot.q95)
####################################

#################################
##### Present value of stocastic annuities
#################################
p_LC_ISLm.st.ext <- 1-q_LC_ISLm.st.ext # 1 year survival probability
ann_LC_ISLm.st <- annuity.st(q_LC_ISLm.st.ext,v_vect,0.9)
ann_LC_ISLm.mean <- ann_LC_ISLm.st[[1]]
ann_LC_ISLm.q95 <- ann_LC_ISLm.st[[2]]
ann_LC_ISLm.q05 <- ann_LC_ISLm.st[[3]]

########
par(mfrow=c(1, 2))
plot(ann_LC_ISLm.mean["65",], type="l", lty=1) # evolution with t (same age)
lines(ann_LC_ISLm.q95["65",], type="l", lty=2)
lines(ann_LC_ISLm.q05["65",], type="l", lty=2)
plot(ann_LC_ISLm.mean[,"2050"], type="l", lty=1) # evolution with x (same year)
lines(ann_LC_ISLm.q95[,"2050"], type="l", lty=2)
lines(ann_LC_ISLm.q05[,"2050"], type="l", lty=2)
dev.off()

########

p_LC_ISLf.st.ext <- 1-q_LC_ISLf.st.ext # 1 year survival probability
ann_LC_ISLf.st <- annuity.st(q_LC_ISLf.st.ext,v_vect,0.9)
ann_LC_ISLf.mean <- ann_LC_ISLf.st[[1]]
ann_LC_ISLf.q95 <- ann_LC_ISLf.st[[2]]
ann_LC_ISLf.q05 <- ann_LC_ISLf.st[[3]]

########
par(mfrow=c(1, 2))
plot(ann_LC_ISLf.mean["65",], type="l", lty=1) # evolution with t (same age)
lines(ann_LC_ISLf.q95["65",], type="l", lty=2)
lines(ann_LC_ISLf.q05["65",], type="l", lty=2)
plot(ann_LC_ISLf.mean[,"2050"], type="l", lty=1) # evolution with x (same year)
lines(ann_LC_ISLf.q95[,"2050"], type="l", lty=2)
lines(ann_LC_ISLf.q05[,"2050"], type="l", lty=2)
dev.off()

########

ann_LC.ISLm.m <- as.data.frame(ann_LC_ISLm.mean)
ann_LC.ISLm.q95 <- as.data.frame(ann_LC_ISLm.q95)
ann_LC.ISLm.q05 <- as.data.frame(ann_LC_ISLm.q05)
ann_LC.ISLf.m <- as.data.frame(ann_LC_ISLf.mean)
ann_LC.ISLf.q95 <-as.data.frame(ann_LC_ISLf.q95)
ann_LC.ISLf.q05<- as.data.frame(ann_LC_ISLf.q05)

##### Present value of stocastic annuities : II method
ann_LC.ISLm.c <- as.data.frame(annuity(q_LC_ISLm.ext, v_vect))
ann_LC.ISLm.q05q <- as.data.frame(annuity(q_LC_ISLm.q95, v_vect))
ann_LC.ISLm.q95q <- as.data.frame(annuity(q_LC_ISLm.q05, v_vect))
ann_LC.ISLf.c <- as.data.frame(annuity(q_LC_ISLf.ext, v_vect))
ann_LC.ISLf.q05q <- as.data.frame(annuity(q_LC_ISLf.q95, v_vect))
ann_LC.ISLf.q95q <- as.data.frame(annuity(q_LC_ISLf.q05, v_vect))
#################################

#################################
##### writing results on an excel file
#################################
write_xlsx(q_LC.ISLm, "q_LC.ISLm.xlsx")
write_xlsx(p_LC.ISLm, "p_LC.ISLm.xlsx")
write_xlsx(p0n_LC.ISLm, "p0n_LC.ISLm.xlsx")

write_xlsx(q_LC.ISLf, "q_LC.ISLf.xlsx")
write_xlsx(p_LC.ISLf, "p_LC.ISLf.xlsx")
write_xlsx(p0n_LC.ISLf, "p0n_LC.ISLf.xlsx")

write_xlsx(q_LC.ISLm.q95, "q_LC.ISLm.q95.xlsx")
write_xlsx(q_LC.ISLm.q05, "q_LC.ISLm.q05.xlsx")
write_xlsx(q_LC.ISLm.boot.q95, "q_LC.ISLm.boot.q95.xlsx")
write_xlsx(q_LC.ISLm.boot.q05, "q_LC.ISLm.boot.q05.xlsx")

write_xlsx(q_LC.ISLf.q95, "q_LC.ISLf.q95.xlsx")
write_xlsx(q_LC.ISLf.q05, "q_LC.ISLf.q05.xlsx")
write_xlsx(q_LC.ISLf.boot.q95, "q_LC.ISLf.boot.q95.xlsx")
write_xlsx(q_LC.ISLf.boot.q05, "q_LC.ISLf.boot.q05.xlsx")

write_xlsx(list(ex_LC.ISLm.boot.q05,ex_LC.ISLm.q05,ex_LC.ISLm,ex_LC.ISLm.q95,ex_LC.ISLm.boot.q95),"ex_LC.ISLm.xlsx")
write_xlsx(list(ex_LC.ISLf.boot.q05,ex_LC.ISLf.q05,ex_LC.ISLf,ex_LC.ISLf.q95,ex_LC.ISLf.boot.q95),"ex_LC.ISLf.xlsx")

write_xlsx(list(ann_LC.ISLm.q05,ann_LC.ISLm.m,ann_LC.ISLm.q95),"ann_LC.ISLm.xlsx")
write_xlsx(list(ann_LC.ISLf.q05,ann_LC.ISLf.m,ann_LC.ISLf.q95),"ann_LC.ISLf.xlsx")

write_xlsx(list(ann_LC.ISLm.q05q,ann_LC.ISLm.c,ann_LC.ISLm.q95q),"ann_LC.ISLmq.xlsx")
write_xlsx(list(ann_LC.ISLf.q05q,ann_LC.ISLf.c,ann_LC.ISLf.q95q),"ann_LC.ISLfq.xlsx")

