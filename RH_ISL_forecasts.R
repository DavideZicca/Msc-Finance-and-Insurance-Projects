source("Functions.r") #uploading functions from an other script in the directory

####################################
################### Renshaw-Haberman model with logit
####################################
################### Forecast
auto.arima(RHfitlogit_ISLm$gc)
auto.arima(RHfitlogit_ISLf$gc)

RHfor_ISLm <- forecast(RHfitlogit_ISLm, h=y.pred, gc.order = c(2, 2, 2))
# default ARIMA order for gc. Alternative order can be choosen via auto.arima
RHfor_ISLf <- forecast(RHfitlogit_ISLf, h=y.pred, gc.order = c(0, 1, 1))
plot(RHfor_ISLm, only.kt = TRUE)
plot(RHfor_ISLf, only.kt = TRUE)
plot(RHfor_ISLm, only.gc = TRUE)
plot(RHfor_ISLf, only.gc = TRUE)
#####
par(mfrow=c(1, 2), oma=c(2,0,0,0))
plot(c(a.min:a.max),log(RHfor_ISLm$rates[,y.pred]), bty="l", xlab="Ages", ylab="Male death rates (log scale)",
  ylim=c(min(log(RHfor_ISLm$rates[,y.pred])),max(log((RHfitlogit_ISLm$Dxt/RHfitlogit_ISLm$Ext)[,length(Y.fit)]))), type="l")
lines(c(a.min:a.max),log((RHfitlogit_ISLm$Dxt/RHfitlogit_ISLm$Ext)[,length(Y.fit)]), type="p")
plot(c(a.min:a.max),log(RHfor_ISLf$rates[,y.pred]), bty="l", xlab="Ages", ylab="Female death rates (log scale)",
  ylim=c(min(log(RHfor_ISLf$rates[,y.pred])),max(log((RHfitlogit_ISLf$Dxt/RHfitlogit_ISLf$Ext)[,length(Y.fit)]))), type="l")
lines(c(a.min:a.max),log((RHfitlogit_ISLf$Dxt/RHfitlogit_ISLf$Ext)[,length(Y.fit)]), type="p")
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend("bottom", c("last observed", "last projected"), pch=c(1,-1), lty=c(0,1), lwd=c(1,1), col=c(1:1), cex=0.9)

dev.off()
########

rates_RHfit_ISLm <- fitted(RHfitlogit_ISLm, type = "rates")
ifelse(is.na(rates_RHfit_ISLm),ISLmRates,rates_RHfit_ISLm)
rates_RHfor_ISLm <- RHfor_ISLm$rates
rates_RH_ISLm <- cbind(rates_RHfit_ISLm,rates_RHfor_ISLm)
#rates_RH_ISLm <- cbind(ISLmRates,rates_RHfor_ISLm)
q_RH_ISLm <- 1- exp(-rates_RH_ISLm)
q_RH_ISLm.ext  <- extrapolation.fit(q_RH_ISLm)

rates_RHfit_ISLf <- fitted(RHfitlogit_ISLf, type = "rates")
ifelse(is.na(rates_RHfit_ISLf),ISLfRates,rates_RHfit_ISLf)
rates_RHfor_ISLf <- RHfor_ISLf$rates
rates_RH_ISLf <- cbind(rates_RHfit_ISLf,rates_RHfor_ISLf)
#rates_RH_ISLf <- cbind(ISLfRates,rates_RHfor_ISLf)
q_RH_ISLf <- 1- exp(-rates_RH_ISLf)
q_RH_ISLf.ext  <- extrapolation.fit(q_RH_ISLf)

#write.table(q_RH_ISLm.ext,file="q_RH.ISLm.txt",sep=",")
#write.table(q_RH_ISLf.ext,file="q_RH.ISLf.txt",sep=",")
q_RH.ISLm<-as.data.frame(q_RH_ISLm.ext)
q_RH.ISLf<-as.data.frame(q_RH_ISLf.ext)

plot(extractCohort(q_RH_ISLm.ext, cohort = 1980), type = "p", log = "y", xlab = "age", ylab = "q(x)",
     main = "Mortality rates for a cohort (log scale)")
lines(extractCohort(q_RH_ISLm, cohort = 1980))

plot(extractCohort(log(q_RH_ISLm.ext/(1-q_RH_ISLm.ext)), cohort = 1980), type = "p", xlab = "age", ylab = "logit(q(x))",
     main = "Logit transform of mortality rates for a cohort")
lines(extractCohort(log(q_RH_ISLm/(1-q_RH_ISLm)), cohort = 1980))
dev.off()
########
############################

####################################
######## 1 year survival probability, n year survival probability and life expectancy
####################################
p_RH_ISLm.ext <- 1-q_RH_ISLm.ext # 1 year survival probability
p0n_RH_ISLm.ext <- apply(p_RH_ISLm.ext, 2, cumprod) # n year survival probability
ex_RH_ISLm.ext <- life.exp(q_RH_ISLm.ext)

p_RH_ISLf.ext <- 1-q_RH_ISLf.ext # 1 year survival probability
p0n_RH_ISLf.ext <- apply(p_RH_ISLf.ext, 2, cumprod) # n year survival probability
ex_RH_ISLf.ext <- life.exp(q_RH_ISLf.ext)

write.table(ex_RH_ISLm.ext,file="ex_RH.ISLm.txt",sep=",")
write.table(ex_RH_ISLf.ext,file="ex_RH.ISLf.txt",sep=",")
p_RH.ISLm<-as.data.frame(p_RH_ISLm.ext)
p0n_RH.ISLm<-as.data.frame(p0n_RH_ISLm.ext)
ex_RH.ISLm<-as.data.frame(ex_RH_ISLm.ext)
p_RH.ISLf<-as.data.frame(p_RH_ISLf.ext)
p0n_RH.ISLf<-as.data.frame(p0n_RH_ISLf.ext)
ex_RH.ISLf<-as.data.frame(ex_RH_ISLf.ext)
####################################

####################################
##### SIMULATION WITH RANDOM WALK WITH DRIFT
####################################
RHsim_ISLm.mrwd <- simulate(RHfitlogit_ISLm, nsim = n.sim, h=y.pred, gc.order = c(2, 2, 2))
### Alternative istructions
#RHsim_ISLmArima <- forecast(RHfitlogit_ISLm, nsim = n.sim, h=y.pred, kt.method = "iarima",  kt.order = c(1, 1, 2))  # ARIMA order specified
#RHsim_ISLmBArima <- forecast(RHfitlogit_ISLm, nsim = n.sim, h=y.pred, kt.method = "iarima",  kt.order = NULL)  # auto.arima
rates_RH_ISLm.st <- RHsim_ISLm.mrwd$rates
q_RH_ISLm.st <- 1- exp(-rates_RH_ISLm.st)
q_RH_ISLm.st.ext <-  extrapolation.sim(q_RH_ISLm.st)

RHsim_ISLf.mrwd <- simulate(RHfitlogit_ISLf, nsim = n.sim, h=y.pred, gc.order = c(0, 1, 1))
rates_RH_ISLf.st <- RHsim_ISLf.mrwd$rates
q_RH_ISLf.st <- 1- exp(-rates_RH_ISLf.st)
q_RH_ISLf.st.ext <-  extrapolation.sim(q_RH_ISLf.st)

par(mfrow=c(1, 2))
plot(RHfitlogit_ISLm$years, RHfitlogit_ISLm$kt[1, ], xlim = range(RHfitlogit_ISLm$years, RHsim_ISLm.mrwd$kt.s$years),
     ylim = range(RHfitlogit_ISLm$kt, RHsim_ISLm.mrwd$kt.s$sim), type = "l", xlab = "year", ylab = "kt",
     main = "Renshaw and Haberman : Period index (mrwd)")
matlines(RHsim_ISLm.mrwd$kt.s$years, RHsim_ISLm.mrwd$kt.s$sim[1, , ], type = "l", lty = 1)
plot(RHfitlogit_ISLm$years, (RHfitlogit_ISLm$Dxt / RHfitlogit_ISLm$Ext)["65", ], xlim = range(RHfitlogit_ISLm$years, RHsim_ISLm.mrwd$years),
     ylim = range((RHfitlogit_ISLm$Dxt / RHfitlogit_ISLm$Ext)["65", ], RHsim_ISLm.mrwd$rates["65", , ]), type = "l", xlab = "year", ylab = "rate",
     main = "Renshaw and Haberman : Simulated mortality rates at age 65")
matlines(RHsim_ISLm.mrwd$years, RHsim_ISLm.mrwd$rates["65", , ], type = "l", lty = 1)
dev.off()
########

library("fanplot")
probs = c(2.5, 10, 25, 50, 75, 90, 97.5)
par(mfrow=c(1, 2))
matplot(RHfitlogit_ISLm$years, t(q_RH_ISLm[c("65", "75", "85"),c(1:length(Y.fit)) ]),
        xlim = c(y.fit.min, (y.fit.min+y.pred)), ylim = c(0.001, 0.2), pch = 20, col = "black",
        log = "y", xlab = "year", ylab = "male mortality rate (log scale)")
fan(t(RHsim_ISLm.mrwd$rates["65", , ]), start = y.fit.max+1, probs = probs, n.fan = 4,
    fan.col = colorRampPalette(c("black", "white")), ln = NULL)
fan(t(RHsim_ISLm.mrwd$rates["75", , ]), start = y.fit.max+1, probs = probs, n.fan = 4,
    fan.col = colorRampPalette(c("red", "white")), ln = NULL)
fan(t(RHsim_ISLm.mrwd$rates["85", , ]), start = y.fit.max+1, probs = probs, n.fan = 4,
    fan.col = colorRampPalette(c("blue", "white")), ln = NULL)
text(y.fit.min+4, q_RH_ISLm[c("65", "75", "85"), 30],
     labels = c("x = 65", "x = 75", "x = 85"))
matplot(RHfitlogit_ISLf$years, t(q_RH_ISLf[c("65", "75", "85"),c(1:length(Y.fit)) ]),
        xlim = c(y.fit.min, (y.fit.min+y.pred)), ylim = c(0.001, 0.2), pch = 20, col = "black",
        log = "y", xlab = "year", ylab = "female mortality rate (log scale)")
fan(t(RHsim_ISLf.mrwd$rates["65", , ]), start = y.fit.max+1, probs = probs, n.fan = 4,
    fan.col = colorRampPalette(c("black", "white")), ln = NULL)
fan(t(RHsim_ISLf.mrwd$rates["75", , ]), start = y.fit.max+1, probs = probs, n.fan = 4,
    fan.col = colorRampPalette(c("red", "white")), ln = NULL)
fan(t(RHsim_ISLf.mrwd$rates["85", , ]), start = y.fit.max+1, probs = probs, n.fan = 4,
    fan.col = colorRampPalette(c("blue", "white")), ln = NULL)
text(y.fit.min+4, q_RH_ISLf[c("65", "75", "85"), 30],
     labels = c("x = 65", "x = 75", "x = 85"))
dev.off()
####################################

####################################
####### Bootstrap
####################################
RHboot_ISLm <- bootstrap(RHfitlogit_ISLm, nBoot = n.boot, type = "residual")
plot(RHboot_ISLm, nCol = 3)
RHsim_ISLm.boot <- simulate(RHboot_ISLm, nsim = n.sim/n.boot, h = y.pred, gc.order = c(2, 2, 2))
rates_RH_ISLm.boot.st <- RHsim_ISLm.boot$rates
q_RH_ISLm.boot.st <- 1- exp(-rates_RH_ISLm.boot.st)
q_RH_ISLm.boot.st.ext <-  extrapolation.sim(q_RH_ISLm.boot.st)

RHboot_ISLf <- bootstrap(RHfitlogit_ISLf, nBoot = n.boot, type = "residual")
plot(RHboot_ISLf, nCol = 3)
RHsim_ISLf.boot <- simulate(RHboot_ISLf, nsim = n.sim/n.boot, h = y.pred, gc.order = c(0, 1, 1))
rates_RH_ISLf.boot.st <- RHsim_ISLf.boot$rates
q_RH_ISLf.boot.st <- 1- exp(-rates_RH_ISLf.boot.st)
q_RH_ISLf.boot.st.ext <-  extrapolation.sim(q_RH_ISLf.boot.st)
dev.off()

####################################

####################################
#### Confidence intervals at 90%
####################################
##### 1 year death probability
conf.lev <- 0.9
q_RH_ISLm.q95 <- matrix(0, nrow=(omega-a.min),ncol=y.pred, dimnames=list(c(a.min:(omega-1)), c((tail(Y.fit,1)+1):(tail(Y.fit,1)+y.pred))))
q_RH_ISLm.q05 <- matrix(0, nrow=(omega-a.min),ncol=y.pred, dimnames=list(c(a.min:(omega-1)), c((tail(Y.fit,1)+1):(tail(Y.fit,1)+y.pred))))
q_RH_ISLm.boot.q95 <- matrix(0, nrow=(omega-a.min),ncol=y.pred, dimnames=list(c(a.min:(omega-1)), c((tail(Y.fit,1)+1):(tail(Y.fit,1)+y.pred))))
q_RH_ISLm.boot.q05 <- matrix(0, nrow=(omega-a.min),ncol=y.pred, dimnames=list(c(a.min:(omega-1)), c((tail(Y.fit,1)+1):(tail(Y.fit,1)+y.pred))))
for (j in 1:(omega-a.min)){
  for (k in 1:y.pred){
    q_RH_ISLm.q95[j,k] <- quantile(q_RH_ISLm.st.ext[j,k,], probs=0.5+conf.lev/2)
    q_RH_ISLm.q05[j,k] <- quantile(q_RH_ISLm.st.ext[j,k,], probs=(1-conf.lev)/2)
    q_RH_ISLm.boot.q95[j,k] <- quantile(q_RH_ISLm.boot.st.ext[j,k,], probs=0.5+conf.lev/2)
    q_RH_ISLm.boot.q05[j,k] <- quantile(q_RH_ISLm.boot.st.ext[j,k,], probs=(1-conf.lev)/2)
  }
}

q_RH_ISLf.q95 <- matrix(0, nrow=(omega-a.min),ncol=y.pred, dimnames=list(c(a.min:(omega-1)), c((tail(Y.fit,1)+1):(tail(Y.fit,1)+y.pred))))
q_RH_ISLf.q05 <- matrix(0, nrow=(omega-a.min),ncol=y.pred, dimnames=list(c(a.min:(omega-1)), c((tail(Y.fit,1)+1):(tail(Y.fit,1)+y.pred))))
q_RH_ISLf.boot.q95 <- matrix(0, nrow=(omega-a.min),ncol=y.pred, dimnames=list(c(a.min:(omega-1)), c((tail(Y.fit,1)+1):(tail(Y.fit,1)+y.pred))))
q_RH_ISLf.boot.q05 <- matrix(0, nrow=(omega-a.min),ncol=y.pred, dimnames=list(c(a.min:(omega-1)), c((tail(Y.fit,1)+1):(tail(Y.fit,1)+y.pred))))
for (j in 1:(omega-a.min)){
  for (k in 1:y.pred){
    q_RH_ISLf.q95[j,k] <- quantile(q_RH_ISLf.st.ext[j,k,], probs=0.5+conf.lev/2)
    q_RH_ISLf.q05[j,k] <- quantile(q_RH_ISLf.st.ext[j,k,], probs=(1-conf.lev)/2)
    q_RH_ISLf.boot.q95[j,k] <- quantile(q_RH_ISLf.boot.st.ext[j,k,], probs=0.5+conf.lev/2)
    q_RH_ISLf.boot.q05[j,k] <- quantile(q_RH_ISLf.boot.st.ext[j,k,], probs=(1-conf.lev)/2)
  }
}
#write.table(q_RH_ISLm.q95,file="q_RH.ISLm.q95.txt",sep=",")
#write.table(q_RH_ISLm.q05,file="q_RH.ISLm.q05.txt",sep=",")
#write.table(q_RH_ISLm.boot.q95,file="q_RH.ISLm.boot.q95.txt",sep=",")
#write.table(q_RH_ISLm.boot.q05,file="q_RH.ISLm.boot.q05.txt",sep=",")
#write.table(q_RH_ISLf.q95,file="q_RH.ISLf.q95.txt",sep=",")
#write.table(q_RH_ISLf.q05,file="q_RH.ISLf.q05.txt",sep=",")
#write.table(q_RH_ISLf.boot.q95,file="q_RH.ISLf.boot.q95.txt",sep=",")
#write.table(q_RH_ISLf.boot.q05,file="q_RH.ISLf.boot.q05.txt",sep=",")
q_RH.ISLm.q95<-as.data.frame(q_RH_ISLm.q95)
q_RH.ISLm.q05<-as.data.frame(q_RH_ISLm.q05)
q_RH.ISLm.boot.q95<-as.data.frame(q_RH_ISLm.boot.q95)
q_RH.ISLm.boot.q05<-as.data.frame(q_RH_ISLm.boot.q05)
q_RH.ISLf.q95<-as.data.frame(q_RH_ISLf.q95)
q_RH.ISLf.q05<-as.data.frame(q_RH_ISLf.q05)
q_RH.ISLf.boot.q95<-as.data.frame(q_RH_ISLf.boot.q95)
q_RH.ISLf.boot.q05<-as.data.frame(q_RH_ISLf.boot.q05)

plot(q_RH_ISLm.ext["65",],x=c(Y.fit,Y.pred), type="l", lty=1) # evolution with t (same age)
lines(q_RH_ISLm.q95["65",], x=Y.pred, type="l", lty=2)
lines(q_RH_ISLm.q05["65",], x=Y.pred, type="l", lty=2)
lines(q_RH_ISLm.boot.q95["65",], x=Y.pred, type="l", lty=2, col=2)
lines(q_RH_ISLm.boot.q05["65",], x=Y.pred, type="l", lty=2, col=2)

plot(q_RH_ISLm.ext["20",],x=c(Y.fit,Y.pred), type="l", lty=1) # evolution with t (same age)
lines(q_RH_ISLm.q95["20",], x=Y.pred, type="l", lty=2)
lines(q_RH_ISLm.q05["20",], x=Y.pred, type="l", lty=2)
lines(q_RH_ISLm.boot.q95["20",], x=Y.pred, type="l", lty=2, col=2)
lines(q_RH_ISLm.boot.q05["20",], x=Y.pred, type="l", lty=2, col=2)
dev.off()
####################################

##### 1 year survival probability, n year survival probability and life expectancy confidence interval
p_RH_ISLm.q95 <- 1-q_RH_ISLm.q95# 1 year survival probability
p0n_RH_ISLm.q95 <- apply(p_RH_ISLm.q95, 2, cumprod) # n year survival probability
ex_RH_ISLm.q95 <- life.exp(q_RH_ISLm.q95) # life expectancy

p_RH_ISLm.q05 <- 1-q_RH_ISLm.q05# 1 year survival probability
p0n_RH_ISLm.q05 <- apply(p_RH_ISLm.q05, 2, cumprod) # n year survival probability
ex_RH_ISLm.q05 <- life.exp(q_RH_ISLm.q05) # life expectancy

p_RH_ISLm.boot.q95 <- 1-q_RH_ISLm.boot.q95# 1 year survival probability
p0n_RH_ISLm.boot.q95 <- apply(p_RH_ISLm.boot.q95, 2, cumprod) # n year survival probability
ex_RH_ISLm.boot.q95 <- life.exp(q_RH_ISLm.boot.q95) # life expectancy

p_RH_ISLm.boot.q05 <- 1-q_RH_ISLm.boot.q05# 1 year survival probability
p0n_RH_ISLm.boot.q05 <- apply(p_RH_ISLm.boot.q05, 2, cumprod) # n year survival probability
ex_RH_ISLm.boot.q05 <- life.exp(q_RH_ISLm.boot.q05) # life expectancy

p_RH_ISLf.q95 <- 1-q_RH_ISLf.q95# 1 year survival probability
p0n_RH_ISLf.q95 <- apply(p_RH_ISLf.q95, 2, cumprod) # n year survival probability
ex_RH_ISLf.q95 <- life.exp(q_RH_ISLf.q95) # life expectancy

p_RH_ISLf.q05 <- 1-q_RH_ISLf.q05# 1 year survival probability
p0n_RH_ISLf.q05 <- apply(p_RH_ISLf.q05, 2, cumprod) # n year survival probability
ex_RH_ISLf.q05 <- life.exp(q_RH_ISLf.q05) # life expectancy

p_RH_ISLf.boot.q95 <- 1-q_RH_ISLf.boot.q95# 1 year survival probability
p0n_RH_ISLf.boot.q95 <- apply(p_RH_ISLf.boot.q95, 2, cumprod) # n year survival probability
ex_RH_ISLf.boot.q95 <- life.exp(q_RH_ISLf.boot.q95) # life expectancy

p_RH_ISLf.boot.q05 <- 1-q_RH_ISLf.boot.q05# 1 year survival probability
p0n_RH_ISLf.boot.q05 <- apply(p_RH_ISLf.boot.q05, 2, cumprod) # n year survival probability
ex_RH_ISLf.boot.q05 <- life.exp(q_RH_ISLf.boot.q05) # life expectancy

ex_RH.ISLm.q05<-as.data.frame(ex_RH_ISLm.q05)
ex_RH.ISLm.q95<-as.data.frame(ex_RH_ISLm.q95)
ex_RH.ISLm.boot.q05<-as.data.frame(ex_RH_ISLm.boot.q05)
ex_RH.ISLm.boot.q95<-as.data.frame(ex_RH_ISLm.boot.q95)
ex_RH.ISLf.q05<-as.data.frame(ex_RH_ISLf.q05)
ex_RH.ISLf.q95<-as.data.frame(ex_RH_ISLf.q95)
ex_RH.ISLf.boot.q05<-as.data.frame(ex_RH_ISLf.boot.q05)
ex_RH.ISLf.boot.q95<-as.data.frame(ex_RH_ISLf.boot.q95)
####################################

#################################
##### Present value of stocastic annuities
#################################
p_RH_ISLm.st.ext <- 1-q_RH_ISLm.st.ext # 1 year survival probability
ann_RH_ISLm.st <- annuity.st(q_RH_ISLm.st.ext,v_vect,0.9)
ann_RH_ISLm.mean <- ann_RH_ISLm.st[[1]]
ann_RH_ISLm.q95 <- ann_RH_ISLm.st[[2]]
ann_RH_ISLm.q05 <- ann_RH_ISLm.st[[3]]

########
par(mfrow=c(1, 2))
plot(ann_RH_ISLm.mean["65",], type="l", lty=1) # evolution with t (same age)
lines(ann_RH_ISLm.q95["65",], type="l", lty=2)
lines(ann_RH_ISLm.q05["65",], type="l", lty=2)
plot(ann_RH_ISLm.mean[,"2050"], type="l", lty=1) # evolution with x (same year)
lines(ann_RH_ISLm.q95[,"2050"], type="l", lty=2)
lines(ann_RH_ISLm.q05[,"2050"], type="l", lty=2)
dev.off()
########

p_RH_ISLf.st.ext <- 1-q_RH_ISLf.st.ext # 1 year survival probability
ann_RH_ISLf.st <- annuity.st(q_RH_ISLf.st.ext,v_vect,0.9)
ann_RH_ISLf.mean <- ann_RH_ISLf.st[[1]]
ann_RH_ISLf.q95 <- ann_RH_ISLf.st[[2]]
ann_RH_ISLf.q05 <- ann_RH_ISLf.st[[3]]

########
par(mfrow=c(1, 2))
plot(ann_RH_ISLf.mean["65",], type="l", lty=1) # evolution with t (same age)
lines(ann_RH_ISLf.q95["65",], type="l", lty=2)
lines(ann_RH_ISLf.q05["65",], type="l", lty=2)
plot(ann_RH_ISLf.mean[,"2050"], type="l", lty=1) # evolution with x (same year)
lines(ann_RH_ISLf.q95[,"2050"], type="l", lty=2)
lines(ann_RH_ISLf.q05[,"2050"], type="l", lty=2)
dev.off()
########

ann_RH.ISLm.m <- as.data.frame(ann_RH_ISLm.mean)
ann_RH.ISLm.q95 <- as.data.frame(ann_RH_ISLm.q95)
ann_RH.ISLm.q05 <- as.data.frame(ann_RH_ISLm.q05)
ann_RH.ISLf.m <- as.data.frame(ann_RH_ISLf.mean)
ann_RH.ISLf.q95 <-as.data.frame(ann_RH_ISLf.q95)
ann_RH.ISLf.q05<- as.data.frame(ann_RH_ISLf.q05)

##### Present value of stocastic annuities : II method
ann_RH.ISLm.c <- as.data.frame(annuity(q_RH_ISLm.ext, v_vect))
ann_RH.ISLm.q05q <- as.data.frame(annuity(q_RH_ISLm.q95, v_vect))
ann_RH.ISLm.q95q <- as.data.frame(annuity(q_RH_ISLm.q05, v_vect))
ann_RH.ISLf.c <- as.data.frame(annuity(q_RH_ISLf.ext, v_vect))
ann_RH.ISLf.q05q <- as.data.frame(annuity(q_RH_ISLf.q95, v_vect))
ann_RH.ISLf.q95q <- as.data.frame(annuity(q_RH_ISLf.q05, v_vect))
#################################

#################################
##### writing results on an excel file
#################################
write_xlsx(q_RH.ISLm, "q_RH.ISLm.xlsx")
write_xlsx(p_RH.ISLm, "p_RH.ISLm.xlsx")
write_xlsx(p0n_RH.ISLm, "p0n_RH.ISLm.xlsx")

write_xlsx(q_RH.ISLf, "q_RH.ISLf.xlsx")
write_xlsx(p_RH.ISLf, "p_RH.ISLf.xlsx")
write_xlsx(p0n_RH.ISLf, "p0n_RH.ISLf.xlsx")

write_xlsx(q_RH.ISLm.q95, "q_RH.ISLm.q95.xlsx")
write_xlsx(q_RH.ISLm.q05, "q_RH.ISLm.q05.xlsx")
write_xlsx(q_RH.ISLm.boot.q95, "q_RH.ISLm.boot.q95.xlsx")
write_xlsx(q_RH.ISLm.boot.q05, "q_RH.ISLm.boot.q05.xlsx")

write_xlsx(q_RH.ISLf.q95, "q_RH.ISLf.q95.xlsx")
write_xlsx(q_RH.ISLf.q05, "q_RH.ISLf.q05.xlsx")
write_xlsx(q_RH.ISLf.boot.q95, "q_RH.ISLf.boot.q95.xlsx")
write_xlsx(q_RH.ISLf.boot.q05, "q_RH.ISLf.boot.q05.xlsx")

write_xlsx(list(ex_RH.ISLm.boot.q05,ex_RH.ISLm.q05,ex_RH.ISLm,ex_RH.ISLm.q95,ex_RH.ISLm.boot.q95),"ex_RH.ISLm.xlsx")
write_xlsx(list(ex_RH.ISLf.boot.q05,ex_RH.ISLf.q05,ex_RH.ISLf,ex_RH.ISLf.q95,ex_RH.ISLf.boot.q95),"ex_RH.ISLf.xlsx")

write_xlsx(list(ann_RH.ISLm.q05,ann_RH.ISLm.m,ann_RH.ISLm.q95),"ann_RH.ISLm.xlsx")
write_xlsx(list(ann_RH.ISLf.q05,ann_RH.ISLf.m,ann_RH.ISLf.q95),"ann_RH.ISLf.xlsx")

write_xlsx(list(ann_RH.ISLm.q05q,ann_RH.ISLm.c,ann_RH.ISLm.q95q),"ann_RH.ISLmq.xlsx")
write_xlsx(list(ann_RH.ISLf.q05q,ann_RH.ISLf.c,ann_RH.ISLf.q95q),"ann_RH.ISLfq.xlsx")

