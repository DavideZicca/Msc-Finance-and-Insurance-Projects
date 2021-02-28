source("Functions.r") #uploading functions from an other script in the directory

####################################
################### PLAT model with logit
####################################
################### Forecast
auto.arima(PLATfitlogit_ISLm$gc)
auto.arima(PLATfitlogit_ISLf$gc)

PLATfor_ISLm <- forecast(PLATfitlogit_ISLm, h=y.pred, gc.order = c(1, 0, 1))
# default ARIMA order for gc. Alternative order can be choosen via auto.arima
PLATfor_ISLf <- forecast(PLATfitlogit_ISLf, h=y.pred, gc.order = c(2, 0, 0))
plot(PLATfor_ISLm, only.kt = TRUE)
plot(PLATfor_ISLf, only.kt = TRUE)
plot(PLATfor_ISLm, only.gc = TRUE)
plot(PLATfor_ISLf, only.gc = TRUE)
#####
par(mfrow=c(1, 2), oma=c(2,0,0,0))
plot(c(a.min:a.max),log(PLATfor_ISLm$rates[,y.pred]), bty="l", xlab="Ages", ylab="Male death rates (log scale)",
  ylim=c(min(log(PLATfor_ISLm$rates[,y.pred])),max(log((PLATfitlogit_ISLm$Dxt/PLATfitlogit_ISLm$Ext)[,length(Y.fit)]))), type="l")
lines(c(a.min:a.max),log((PLATfitlogit_ISLm$Dxt/PLATfitlogit_ISLm$Ext)[,length(Y.fit)]), type="p")
plot(c(a.min:a.max),log(PLATfor_ISLf$rates[,y.pred]), bty="l", xlab="Ages", ylab="Female death rates (log scale)",
  ylim=c(min(log(PLATfor_ISLf$rates[,y.pred])),max(log((PLATfitlogit_ISLf$Dxt/PLATfitlogit_ISLf$Ext)[,length(Y.fit)]))), type="l")
lines(c(a.min:a.max),log((PLATfitlogit_ISLf$Dxt/PLATfitlogit_ISLf$Ext)[,length(Y.fit)]), type="p")
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend("bottom", c("last observed", "last projected"), pch=c(1,-1), lty=c(0,1), lwd=c(1,1), col=c(1:1), cex=0.9)
dev.off()
########
rates_PLATfit_ISLm <- fitted(PLATfitlogit_ISLm, type = "rates")
ifelse(is.na(rates_PLATfit_ISLm),ISLmRates,rates_PLATfit_ISLm)
rates_PLATfor_ISLm <- PLATfor_ISLm$rates
rates_PLAT_ISLm <- cbind(rates_PLATfit_ISLm,rates_PLATfor_ISLm)
#rates_PLAT_ISLm <- cbind(ISLmRates,rates_PLATfor_ISLm)
q_PLAT_ISLm <- 1- exp(-rates_PLAT_ISLm)
q_PLAT_ISLm.ext  <- extrapolation.fit(q_PLAT_ISLm)

rates_PLATfit_ISLf <- fitted(PLATfitlogit_ISLf, type = "rates")
ifelse(is.na(rates_PLATfit_ISLf),ISLfRates,rates_PLATfit_ISLf)
rates_PLATfor_ISLf <- PLATfor_ISLf$rates
rates_PLAT_ISLf <- cbind(rates_PLATfit_ISLf,rates_PLATfor_ISLf)
#rates_PLAT_ISLf <- cbind(ISLfRates,rates_PLATfor_ISLf)
q_PLAT_ISLf <- 1- exp(-rates_PLAT_ISLf)
q_PLAT_ISLf.ext  <- extrapolation.fit(q_PLAT_ISLf)

#write.table(q_PLAT_ISLm.ext,file="q_PLAT.ISLm.txt",sep=",")
#write.table(q_PLAT_ISLf.ext,file="q_PLAT.ISLf.txt",sep=",")
q_PLAT.ISLm<-as.data.frame(q_PLAT_ISLm.ext)
q_PLAT.ISLf<-as.data.frame(q_PLAT_ISLf.ext)

plot(extractCohort(q_PLAT_ISLm.ext, cohort = 1980), type = "p", log = "y", xlab = "age", ylab = "q(x)",
     main = "Mortality rates for a cohort (log scale)")
lines(extractCohort(q_PLAT_ISLm, cohort = 1980))

plot(extractCohort(log(q_PLAT_ISLm.ext/(1-q_PLAT_ISLm.ext)), cohort = 1980), type = "p", xlab = "age", ylab = "logit(q(x))",
     main = "Logit transform of mortality rates for a cohort")
lines(extractCohort(log(q_PLAT_ISLm/(1-q_PLAT_ISLm)), cohort = 1980))
dev.off()
############################

####################################
######## 1 year survival probability, n year survival probability and life expectancy
####################################
p_PLAT_ISLm.ext <- 1-q_PLAT_ISLm.ext # 1 year survival probability
p0n_PLAT_ISLm.ext <- apply(p_PLAT_ISLm.ext, 2, cumprod) # n year survival probability
ex_PLAT_ISLm.ext <- life.exp(q_PLAT_ISLm.ext)

p_PLAT_ISLf.ext <- 1-q_PLAT_ISLf.ext # 1 year survival probability
p0n_PLAT_ISLf.ext <- apply(p_PLAT_ISLf.ext, 2, cumprod) # n year survival probability
ex_PLAT_ISLf.ext <- life.exp(q_PLAT_ISLf.ext)

write.table(ex_PLAT_ISLm.ext,file="ex_PLAT.ISLm.txt",sep=",")
write.table(ex_PLAT_ISLf.ext,file="ex_PLAT.ISLf.txt",sep=",")
p_PLAT.ISLm<-as.data.frame(p_PLAT_ISLm.ext)
p0n_PLAT.ISLm<-as.data.frame(p0n_PLAT_ISLm.ext)
ex_PLAT.ISLm<-as.data.frame(ex_PLAT_ISLm.ext)
p_PLAT.ISLf<-as.data.frame(p_PLAT_ISLf.ext)
p0n_PLAT.ISLf<-as.data.frame(p0n_PLAT_ISLf.ext)
ex_PLAT.ISLf<-as.data.frame(ex_PLAT_ISLf.ext)
####################################

####################################
##### SIMULATION WITH RANDOM WALK WITH DRIFT
####################################
PLATsim_ISLm.mrwd <- simulate(PLATfitlogit_ISLm, nsim = n.sim, h=y.pred, gc.order = c(1, 0, 1))
### Alternative istructions
#PLATsim_ISLmArima <- forecast(PLATfitlogit_ISLm, nsim = n.sim, h=y.pred, kt.method = "iarima",  kt.order = c(1, 1, 2))  # ARIMA order specified
#PLATsim_ISLmBArima <- forecast(PLATfitlogit_ISLm, nsim = n.sim, h=y.pred, kt.method = "iarima",  kt.order = NULL)  # auto.arima
rates_PLAT_ISLm.st <- PLATsim_ISLm.mrwd$rates
q_PLAT_ISLm.st <- 1- exp(-rates_PLAT_ISLm.st)
q_PLAT_ISLm.st.ext <-  extrapolation.sim(q_PLAT_ISLm.st)

PLATsim_ISLf.mrwd <- simulate(PLATfitlogit_ISLf, nsim = n.sim, h=y.pred, gc.order = c(2, 0, 0))
rates_PLAT_ISLf.st <- PLATsim_ISLf.mrwd$rates
q_PLAT_ISLf.st <- 1- exp(-rates_PLAT_ISLf.st)
q_PLAT_ISLf.st.ext <-  extrapolation.sim(q_PLAT_ISLf.st)

par(mfrow=c(1, 2))
plot(PLATfitlogit_ISLm$years, PLATfitlogit_ISLm$kt[1, ], xlim = range(PLATfitlogit_ISLm$years, PLATsim_ISLm.mrwd$kt.s$years),
     ylim = range(PLATfitlogit_ISLm$kt, PLATsim_ISLm.mrwd$kt.s$sim), type = "l", xlab = "year", ylab = "kt",
     main = "PLAT: Period index (mrwd)")
matlines(PLATsim_ISLm.mrwd$kt.s$years, PLATsim_ISLm.mrwd$kt.s$sim[1, , ], type = "l", lty = 1)
plot(PLATfitlogit_ISLm$years, (PLATfitlogit_ISLm$Dxt / PLATfitlogit_ISLm$Ext)["65", ], xlim = range(PLATfitlogit_ISLm$years, PLATsim_ISLm.mrwd$years),
     ylim = range((PLATfitlogit_ISLm$Dxt / PLATfitlogit_ISLm$Ext)["65", ], PLATsim_ISLm.mrwd$rates["65", , ]), type = "l", xlab = "year", ylab = "rate",
     main = "PLAT: Simulated mortality rates at age 65")
matlines(PLATsim_ISLm.mrwd$years, PLATsim_ISLm.mrwd$rates["65", , ], type = "l", lty = 1)
dev.off()
########

library("fanplot")
probs = c(2.5, 10, 25, 50, 75, 90, 97.5)
par(mfrow=c(1, 2))
matplot(PLATfitlogit_ISLm$years, t(q_PLAT_ISLm[c("65", "75", "85"),c(1:length(Y.fit)) ]),
        xlim = c(y.fit.min, (y.fit.min+y.pred)), ylim = c(0.001, 0.2), pch = 20, col = "black",
        log = "y", xlab = "year", ylab = "male mortality rate (log scale)")
fan(t(PLATsim_ISLm.mrwd$rates["65", , ]), start = y.fit.max+1, probs = probs, n.fan = 4,
    fan.col = colorRampPalette(c("black", "white")), ln = NULL)
fan(t(PLATsim_ISLm.mrwd$rates["75", , ]), start = y.fit.max+1, probs = probs, n.fan = 4,
    fan.col = colorRampPalette(c("red", "white")), ln = NULL)
fan(t(PLATsim_ISLm.mrwd$rates["85", , ]), start = y.fit.max+1, probs = probs, n.fan = 4,
    fan.col = colorRampPalette(c("blue", "white")), ln = NULL)
text(y.fit.min+4, q_PLAT_ISLm[c("65", "75", "85"), 30],
     labels = c("x = 65", "x = 75", "x = 85"))
matplot(PLATfitlogit_ISLf$years, t(q_PLAT_ISLf[c("65", "75", "85"),c(1:length(Y.fit)) ]),
        xlim = c(y.fit.min, (y.fit.min+y.pred)), ylim = c(0.001, 0.2), pch = 20, col = "black",
        log = "y", xlab = "year", ylab = "female mortality rate (log scale)")
fan(t(PLATsim_ISLf.mrwd$rates["65", , ]), start = y.fit.max+1, probs = probs, n.fan = 4,
    fan.col = colorRampPalette(c("black", "white")), ln = NULL)
fan(t(PLATsim_ISLf.mrwd$rates["75", , ]), start = y.fit.max+1, probs = probs, n.fan = 4,
    fan.col = colorRampPalette(c("red", "white")), ln = NULL)
fan(t(PLATsim_ISLf.mrwd$rates["85", , ]), start = y.fit.max+1, probs = probs, n.fan = 4,
    fan.col = colorRampPalette(c("blue", "white")), ln = NULL)
text(y.fit.min+4, q_PLAT_ISLf[c("65", "75", "85"), 30],
     labels = c("x = 65", "x = 75", "x = 85"))
dev.off()
####################################

####################################
####### Bootstrap
####################################
PLATboot_ISLm <- bootstrap(PLATfitlogit_ISLm, nBoot = n.boot, type = "semiparametric")
par(mar = c(2, 2, 2, 2))
plot(PLATboot_ISLm, nCol = 3)
PLATsim_ISLm.boot <- simulate(PLATboot_ISLm, nsim = n.sim/n.boot, h = y.pred, gc.order = c(1, 0,1))
rates_PLAT_ISLm.boot.st <- PLATsim_ISLm.boot$rates
q_PLAT_ISLm.boot.st <- 1- exp(-rates_PLAT_ISLm.boot.st)
q_PLAT_ISLm.boot.st.ext <-  extrapolation.sim(q_PLAT_ISLm.boot.st)

PLATboot_ISLf <- bootstrap(PLATfitlogit_ISLf, nBoot = n.boot, type = "semiparametric")
plot(PLATboot_ISLf, nCol = 3)
PLATsim_ISLf.boot <- simulate(PLATboot_ISLf, nsim = n.sim/n.boot, h = y.pred, gc.order = c(2, 0, 0))
rates_PLAT_ISLf.boot.st <- PLATsim_ISLf.boot$rates
q_PLAT_ISLf.boot.st <- 1- exp(-rates_PLAT_ISLf.boot.st)
q_PLAT_ISLf.boot.st.ext <-  extrapolation.sim(q_PLAT_ISLf.boot.st)
dev.off()
####################################

####################################
#### Confidence intervals at 90%
####################################
##### 1 year death probability
conf.lev <- 0.9
q_PLAT_ISLm.q95 <- matrix(0, nrow=(omega-a.min),ncol=y.pred, dimnames=list(c(a.min:(omega-1)), c((tail(Y.fit,1)+1):(tail(Y.fit,1)+y.pred))))
q_PLAT_ISLm.q05 <- matrix(0, nrow=(omega-a.min),ncol=y.pred, dimnames=list(c(a.min:(omega-1)), c((tail(Y.fit,1)+1):(tail(Y.fit,1)+y.pred))))
q_PLAT_ISLm.boot.q95 <- matrix(0, nrow=(omega-a.min),ncol=y.pred, dimnames=list(c(a.min:(omega-1)), c((tail(Y.fit,1)+1):(tail(Y.fit,1)+y.pred))))
q_PLAT_ISLm.boot.q05 <- matrix(0, nrow=(omega-a.min),ncol=y.pred, dimnames=list(c(a.min:(omega-1)), c((tail(Y.fit,1)+1):(tail(Y.fit,1)+y.pred))))
for (j in 1:(omega-a.min)){
  for (k in 1:y.pred){
    q_PLAT_ISLm.q95[j,k] <- quantile(q_PLAT_ISLm.st.ext[j,k,], probs=0.5+conf.lev/2)
    q_PLAT_ISLm.q05[j,k] <- quantile(q_PLAT_ISLm.st.ext[j,k,], probs=(1-conf.lev)/2)
    q_PLAT_ISLm.boot.q95[j,k] <- quantile(q_PLAT_ISLm.boot.st.ext[j,k,], probs=0.5+conf.lev/2)
    q_PLAT_ISLm.boot.q05[j,k] <- quantile(q_PLAT_ISLm.boot.st.ext[j,k,], probs=(1-conf.lev)/2)
  }
}

q_PLAT_ISLf.q95 <- matrix(0, nrow=(omega-a.min),ncol=y.pred, dimnames=list(c(a.min:(omega-1)), c((tail(Y.fit,1)+1):(tail(Y.fit,1)+y.pred))))
q_PLAT_ISLf.q05 <- matrix(0, nrow=(omega-a.min),ncol=y.pred, dimnames=list(c(a.min:(omega-1)), c((tail(Y.fit,1)+1):(tail(Y.fit,1)+y.pred))))
q_PLAT_ISLf.boot.q95 <- matrix(0, nrow=(omega-a.min),ncol=y.pred, dimnames=list(c(a.min:(omega-1)), c((tail(Y.fit,1)+1):(tail(Y.fit,1)+y.pred))))
q_PLAT_ISLf.boot.q05 <- matrix(0, nrow=(omega-a.min),ncol=y.pred, dimnames=list(c(a.min:(omega-1)), c((tail(Y.fit,1)+1):(tail(Y.fit,1)+y.pred))))
for (j in 1:(omega-a.min)){
  for (k in 1:y.pred){
    q_PLAT_ISLf.q95[j,k] <- quantile(q_PLAT_ISLf.st.ext[j,k,], probs=0.5+conf.lev/2)
    q_PLAT_ISLf.q05[j,k] <- quantile(q_PLAT_ISLf.st.ext[j,k,], probs=(1-conf.lev)/2)
    q_PLAT_ISLf.boot.q95[j,k] <- quantile(q_PLAT_ISLf.boot.st.ext[j,k,], probs=0.5+conf.lev/2)
    q_PLAT_ISLf.boot.q05[j,k] <- quantile(q_PLAT_ISLf.boot.st.ext[j,k,], probs=(1-conf.lev)/2)
  }
}
#write.table(q_PLAT_ISLm.q95,file="q_PLAT.ISLm.q95.txt",sep=",")
#write.table(q_PLAT_ISLm.q05,file="q_PLAT.ISLm.q05.txt",sep=",")
#write.table(q_PLAT_ISLm.boot.q95,file="q_PLAT.ISLm.boot.q95.txt",sep=",")
#write.table(q_PLAT_ISLm.boot.q05,file="q_PLAT.ISLm.boot.q05.txt",sep=",")
#write.table(q_PLAT_ISLf.q95,file="q_PLAT.ISLf.q95.txt",sep=",")
#write.table(q_PLAT_ISLf.q05,file="q_PLAT.ISLf.q05.txt",sep=",")
#write.table(q_PLAT_ISLf.boot.q95,file="q_PLAT.ISLf.boot.q95.txt",sep=",")
#write.table(q_PLAT_ISLf.boot.q05,file="q_PLAT.ISLf.boot.q05.txt",sep=",")
q_PLAT.ISLm.q95<-as.data.frame(q_PLAT_ISLm.q95)
q_PLAT.ISLm.q05<-as.data.frame(q_PLAT_ISLm.q05)
q_PLAT.ISLm.boot.q95<-as.data.frame(q_PLAT_ISLm.boot.q95)
q_PLAT.ISLm.boot.q05<-as.data.frame(q_PLAT_ISLm.boot.q05)
q_PLAT.ISLf.q95<-as.data.frame(q_PLAT_ISLf.q95)
q_PLAT.ISLf.q05<-as.data.frame(q_PLAT_ISLf.q05)
q_PLAT.ISLf.boot.q95<-as.data.frame(q_PLAT_ISLf.boot.q95)
q_PLAT.ISLf.boot.q05<-as.data.frame(q_PLAT_ISLf.boot.q05)

plot(q_PLAT_ISLm.ext["65",],x=c(Y.fit,Y.pred), type="l", lty=1) # evolution with t (same age)
lines(q_PLAT_ISLm.q95["65",], x=Y.pred, type="l", lty=2)
lines(q_PLAT_ISLm.q05["65",], x=Y.pred, type="l", lty=2)
lines(q_PLAT_ISLm.boot.q95["65",], x=Y.pred, type="l", lty=2, col=2)
lines(q_PLAT_ISLm.boot.q05["65",], x=Y.pred, type="l", lty=2, col=2)

plot(q_PLAT_ISLm.ext["20",],x=c(Y.fit,Y.pred), type="l", lty=1) # evolution with t (same age)
lines(q_PLAT_ISLm.q95["20",], x=Y.pred, type="l", lty=2)
lines(q_PLAT_ISLm.q05["20",], x=Y.pred, type="l", lty=2)
lines(q_PLAT_ISLm.boot.q95["20",], x=Y.pred, type="l", lty=2, col=2)
lines(q_PLAT_ISLm.boot.q05["20",], x=Y.pred, type="l", lty=2, col=2)
dev.off()
####################################

##### 1 year survival probability, n year survival probability and life expectancy confidence interval
p_PLAT_ISLm.q95 <- 1-q_PLAT_ISLm.q95# 1 year survival probability
p0n_PLAT_ISLm.q95 <- apply(p_PLAT_ISLm.q95, 2, cumprod) # n year survival probability
ex_PLAT_ISLm.q95 <- life.exp(q_PLAT_ISLm.q95) # life expectancy

p_PLAT_ISLm.q05 <- 1-q_PLAT_ISLm.q05# 1 year survival probability
p0n_PLAT_ISLm.q05 <- apply(p_PLAT_ISLm.q05, 2, cumprod) # n year survival probability
ex_PLAT_ISLm.q05 <- life.exp(q_PLAT_ISLm.q05) # life expectancy

p_PLAT_ISLm.boot.q95 <- 1-q_PLAT_ISLm.boot.q95# 1 year survival probability
p0n_PLAT_ISLm.boot.q95 <- apply(p_PLAT_ISLm.boot.q95, 2, cumprod) # n year survival probability
ex_PLAT_ISLm.boot.q95 <- life.exp(q_PLAT_ISLm.boot.q95) # life expectancy

p_PLAT_ISLm.boot.q05 <- 1-q_PLAT_ISLm.boot.q05# 1 year survival probability
p0n_PLAT_ISLm.boot.q05 <- apply(p_PLAT_ISLm.boot.q05, 2, cumprod) # n year survival probability
ex_PLAT_ISLm.boot.q05 <- life.exp(q_PLAT_ISLm.boot.q05) # life expectancy

p_PLAT_ISLf.q95 <- 1-q_PLAT_ISLf.q95# 1 year survival probability
p0n_PLAT_ISLf.q95 <- apply(p_PLAT_ISLf.q95, 2, cumprod) # n year survival probability
ex_PLAT_ISLf.q95 <- life.exp(q_PLAT_ISLf.q95) # life expectancy

p_PLAT_ISLf.q05 <- 1-q_PLAT_ISLf.q05# 1 year survival probability
p0n_PLAT_ISLf.q05 <- apply(p_PLAT_ISLf.q05, 2, cumprod) # n year survival probability
ex_PLAT_ISLf.q05 <- life.exp(q_PLAT_ISLf.q05) # life expectancy

p_PLAT_ISLf.boot.q95 <- 1-q_PLAT_ISLf.boot.q95# 1 year survival probability
p0n_PLAT_ISLf.boot.q95 <- apply(p_PLAT_ISLf.boot.q95, 2, cumprod) # n year survival probability
ex_PLAT_ISLf.boot.q95 <- life.exp(q_PLAT_ISLf.boot.q95) # life expectancy

p_PLAT_ISLf.boot.q05 <- 1-q_PLAT_ISLf.boot.q05# 1 year survival probability
p0n_PLAT_ISLf.boot.q05 <- apply(p_PLAT_ISLf.boot.q05, 2, cumprod) # n year survival probability
ex_PLAT_ISLf.boot.q05 <- life.exp(q_PLAT_ISLf.boot.q05) # life expectancy

ex_PLAT.ISLm.q05<-as.data.frame(ex_PLAT_ISLm.q05)
ex_PLAT.ISLm.q95<-as.data.frame(ex_PLAT_ISLm.q95)
ex_PLAT.ISLm.boot.q05<-as.data.frame(ex_PLAT_ISLm.boot.q05)
ex_PLAT.ISLm.boot.q95<-as.data.frame(ex_PLAT_ISLm.boot.q95)
ex_PLAT.ISLf.q05<-as.data.frame(ex_PLAT_ISLf.q05)
ex_PLAT.ISLf.q95<-as.data.frame(ex_PLAT_ISLf.q95)
ex_PLAT.ISLf.boot.q05<-as.data.frame(ex_PLAT_ISLf.boot.q05)
ex_PLAT.ISLf.boot.q95<-as.data.frame(ex_PLAT_ISLf.boot.q95)
####################################

#################################
##### Present value of stocastic annuities
#################################
p_PLAT_ISLm.st.ext <- 1-q_PLAT_ISLm.st.ext # 1 year survival probability
ann_PLAT_ISLm.st <- annuity.st(q_PLAT_ISLm.st.ext,v_vect,0.9)
ann_PLAT_ISLm.mean <- ann_PLAT_ISLm.st[[1]]
ann_PLAT_ISLm.q95 <- ann_PLAT_ISLm.st[[2]]
ann_PLAT_ISLm.q05 <- ann_PLAT_ISLm.st[[3]]

########
par(mfrow=c(1, 2))
plot(ann_PLAT_ISLm.mean["65",], type="l", lty=1) # evolution with t (same age)
lines(ann_PLAT_ISLm.q95["65",], type="l", lty=2)
lines(ann_PLAT_ISLm.q05["65",], type="l", lty=2)
plot(ann_PLAT_ISLm.mean[,"2050"], type="l", lty=1) # evolution with x (same year)
lines(ann_PLAT_ISLm.q95[,"2050"], type="l", lty=2)
lines(ann_PLAT_ISLm.q05[,"2050"], type="l", lty=2)
dev.off()
########

p_PLAT_ISLf.st.ext <- 1-q_PLAT_ISLf.st.ext # 1 year survival probability
ann_PLAT_ISLf.st <- annuity.st(q_PLAT_ISLf.st.ext,v_vect,0.9)
ann_PLAT_ISLf.mean <- ann_PLAT_ISLf.st[[1]]
ann_PLAT_ISLf.q95 <- ann_PLAT_ISLf.st[[2]]
ann_PLAT_ISLf.q05 <- ann_PLAT_ISLf.st[[3]]

########
par(mfrow=c(1, 2))
plot(ann_PLAT_ISLf.mean["65",], type="l", lty=1) # evolution with t (same age)
lines(ann_PLAT_ISLf.q95["65",], type="l", lty=2)
lines(ann_PLAT_ISLf.q05["65",], type="l", lty=2)
plot(ann_PLAT_ISLf.mean[,"2050"], type="l", lty=1) # evolution with x (same year)
lines(ann_PLAT_ISLf.q95[,"2050"], type="l", lty=2)
lines(ann_PLAT_ISLf.q05[,"2050"], type="l", lty=2)
dev.off()
########

ann_PLAT.ISLm.m <- as.data.frame(ann_PLAT_ISLm.mean)
ann_PLAT.ISLm.q95 <- as.data.frame(ann_PLAT_ISLm.q95)
ann_PLAT.ISLm.q05 <- as.data.frame(ann_PLAT_ISLm.q05)
ann_PLAT.ISLf.m <- as.data.frame(ann_PLAT_ISLf.mean)
ann_PLAT.ISLf.q95 <-as.data.frame(ann_PLAT_ISLf.q95)
ann_PLAT.ISLf.q05<- as.data.frame(ann_PLAT_ISLf.q05)

##### Present value of stocastic annuities : II method
ann_PLAT.ISLm.c <- as.data.frame(annuity(q_PLAT_ISLm.ext, v_vect))
ann_PLAT.ISLm.q05q <- as.data.frame(annuity(q_PLAT_ISLm.q95, v_vect))
ann_PLAT.ISLm.q95q <- as.data.frame(annuity(q_PLAT_ISLm.q05, v_vect))
ann_PLAT.ISLf.c <- as.data.frame(annuity(q_PLAT_ISLf.ext, v_vect))
ann_PLAT.ISLf.q05q <- as.data.frame(annuity(q_PLAT_ISLf.q95, v_vect))
ann_PLAT.ISLf.q95q <- as.data.frame(annuity(q_PLAT_ISLf.q05, v_vect))
#################################

#################################
##### writing results on an excel file
#################################
write_xlsx(q_PLAT.ISLm, "q_PLAT.ISLm.xlsx")
write_xlsx(p_PLAT.ISLm, "p_PLAT.ISLm.xlsx")
write_xlsx(p0n_PLAT.ISLm, "p0n_PLAT.ISLm.xlsx")

write_xlsx(q_PLAT.ISLf, "q_PLAT.ISLf.xlsx")
write_xlsx(p_PLAT.ISLf, "p_PLAT.ISLf.xlsx")
write_xlsx(p0n_PLAT.ISLf, "p0n_PLAT.ISLf.xlsx")

write_xlsx(q_PLAT.ISLm.q95, "q_PLAT.ISLm.q95.xlsx")
write_xlsx(q_PLAT.ISLm.q05, "q_PLAT.ISLm.q05.xlsx")
write_xlsx(q_PLAT.ISLm.boot.q95, "q_PLAT.ISLm.boot.q95.xlsx")
write_xlsx(q_PLAT.ISLm.boot.q05, "q_PLAT.ISLm.boot.q05.xlsx")

write_xlsx(q_PLAT.ISLf.q95, "q_PLAT.ISLf.q95.xlsx")
write_xlsx(q_PLAT.ISLf.q05, "q_PLAT.ISLf.q05.xlsx")
write_xlsx(q_PLAT.ISLf.boot.q95, "q_PLAT.ISLf.boot.q95.xlsx")
write_xlsx(q_PLAT.ISLf.boot.q05, "q_PLAT.ISLf.boot.q05.xlsx")

write_xlsx(list(ex_PLAT.ISLm.boot.q05,ex_PLAT.ISLm.q05,ex_PLAT.ISLm,ex_PLAT.ISLm.q95,ex_PLAT.ISLm.boot.q95),"ex_PLAT.ISLm.xlsx")
write_xlsx(list(ex_PLAT.ISLf.boot.q05,ex_PLAT.ISLf.q05,ex_PLAT.ISLf,ex_PLAT.ISLf.q95,ex_PLAT.ISLf.boot.q95),"ex_PLAT.ISLf.xlsx")

write_xlsx(list(ann_PLAT.ISLm.q05,ann_PLAT.ISLm.m,ann_PLAT.ISLm.q95),"ann_PLAT.ISLm.xlsx")
write_xlsx(list(ann_PLAT.ISLf.q05,ann_PLAT.ISLf.m,ann_PLAT.ISLf.q95),"ann_PLAT.ISLf.xlsx")

write_xlsx(list(ann_PLAT.ISLm.q05q,ann_PLAT.ISLm.c,ann_PLAT.ISLm.q95q),"ann_PLAT.ISLmq.xlsx")
write_xlsx(list(ann_PLAT.ISLf.q05q,ann_PLAT.ISLf.c,ann_PLAT.ISLf.q95q),"ann_PLAT.ISLfq.xlsx")
