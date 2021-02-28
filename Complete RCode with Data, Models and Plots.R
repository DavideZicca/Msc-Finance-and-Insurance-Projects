#DATA
library(tidyquant)

# Downloading SP500 using the library tidyquant

SP500= getSymbols("^GSPC", from = '1987-01-01',
                  to = "1998-01-01", auto.assign = FALSE)
head(SP500)

logret = diff(log(SP500$GSPC.Adjusted))
log_ret <- logret[-1]

plot(log_ret, length(log_ret))
#Stationarity test
library(tseries)
adf.test(log_ret, alternative = "stationary")

#ARCH EFFECT

library(fDMA)
archtest(as.vector(log_ret))

#DUMMY VARIABLES

n=length(log_ret)
#10/22/87

e1=which(row.names(as.matrix(log_ret))=="1987-10-20")
f1=which(row.names(as.matrix(log_ret))=="1987-10-21")
g1=which(row.names(as.matrix(log_ret))=="1987-10-22")
h1=which(row.names(as.matrix(log_ret))=="1987-10-23")
j1=which(row.names(as.matrix(log_ret))=="1987-10-26")

d1=rep(0,n)

d1[e1]=1
d1[f1]=1
d1[g1]=1
d1[h1]=1
d1[j1]=1

#1/14/88

e2=which(row.names(as.matrix(log_ret))=="1988-01-12")
f2=which(row.names(as.matrix(log_ret))=="1988-01-13")
g2=which(row.names(as.matrix(log_ret))=="1988-01-14")
h2=which(row.names(as.matrix(log_ret))=="1988-01-15")
j2=which(row.names(as.matrix(log_ret))=="1988-01-18")

d2=rep(0,n)

d2[e2]=1
d2[f2]=1
d2[g2]=1
d2[h2]=1
d2[j2]=1

#2/4/88
e3=which(row.names(as.matrix(log_ret))=="1988-02-02")
f3=which(row.names(as.matrix(log_ret))=="1988-02-03")
g3=which(row.names(as.matrix(log_ret))=="1988-02-04")
h3=which(row.names(as.matrix(log_ret))=="1988-02-05")
j3=which(row.names(as.matrix(log_ret))=="1988-02-08")

d3=rep(0,n)

d3[e3]=1
d3[f3]=1
d3[g3]=1
d3[h3]=1
d3[j3]=1

#10/20/88

e4=which(row.names(as.matrix(log_ret))=="1988-10-18")
f4=which(row.names(as.matrix(log_ret))=="1988-10-19")
g4=which(row.names(as.matrix(log_ret))=="1988-10-20")
h4=which(row.names(as.matrix(log_ret))=="1988-10-21")
j4=which(row.names(as.matrix(log_ret))=="1988-10-24")

d4=rep(0,n)

d4[e4]=1
d4[f4]=1
d4[g4]=1
d4[h4]=1
d4[j4]=1
#8/1/90

e5=which(row.names(as.matrix(log_ret))=="1990-07-30")
f5=which(row.names(as.matrix(log_ret))=="1990-07-31")
g5=which(row.names(as.matrix(log_ret))=="1990-08-01")
h5=which(row.names(as.matrix(log_ret))=="1990-08-02")
j5=which(row.names(as.matrix(log_ret))=="1990-08-03")

d5=rep(0,n)

d5[e5]=1
d5[f5]=1
d5[g5]=1
d5[h5]=1
d5[j5]=1
#7/22/96

e6=which(row.names(as.matrix(log_ret))=="1996-07-18")
f6=which(row.names(as.matrix(log_ret))=="1996-07-19")
g6=which(row.names(as.matrix(log_ret))=="1996-07-22")
h6=which(row.names(as.matrix(log_ret))=="1996-07-23")
j6=which(row.names(as.matrix(log_ret))=="1996-07-24")

d6=rep(0,n)

d6[e6]=1
d6[f6]=1
d6[g6]=1
d6[h6]=1
d6[j6]=1
#3/3/97
e7=which(row.names(as.matrix(log_ret))=="1997-02-27")
f7=which(row.names(as.matrix(log_ret))=="1997-02-28")
g7=which(row.names(as.matrix(log_ret))=="1997-03-03")
h7=which(row.names(as.matrix(log_ret))=="1997-03-04")
j7=which(row.names(as.matrix(log_ret))=="1997-03-05")

d7=rep(0,n)

d7[e7]=1
d7[f7]=1
d7[g7]=1
d7[h7]=1
d7[j7]=1

df_dummies=matrix(c(d1,d2,d3,d4,d5,d6,d7), ncol = 7)

# GARCH MODELS assuming a Normal distribution
#sGARCH
library(rugarch)

s_garchMod <- ugarchspec(mean.model = list(armaOrder = c(1, 0), include.mean = TRUE
), 
variance.model = list(model = 'sGARCH', 
                      garchOrder = c(1, 1),
                      external.regressors = df_dummies
),
distribution.model = "norm")

s_garchFit <- ugarchfit(spec=s_garchMod, data=log_ret)
s_garchFit


## Results review 1
coef(s_garchFit)

s_rhat <- s_garchFit@fit$fitted.values
plot.ts(s_rhat)
s_hhat <- ts(s_garchFit@fit$sigma^2)
plot.ts(s_hhat)
dev.off()
## Results review 2
fit.val     <- coef(s_garchFit)
fit.sd      <- diag(vcov(s_garchFit))
true.val = s_garchFit@fit$tval

fit.conf.lb <- fit.val + qnorm(0.025) * fit.sd
fit.conf.ub <- fit.val + qnorm(0.975) * fit.sd
print(fit.val)

print(fit.sd)

print(true.val)
plot(true.val, pch = 1, col = "red",
     ylim = range(c(fit.conf.lb, fit.conf.ub, true.val)),
     xlab = "", ylab = "", axes = FALSE)
box(); axis(1, at = 1:length(fit.val), labels = names(fit.val)); axis(2)
points(coef(s_garchFit), col = "blue", pch = 7)
for (i in 1:length(fit.val)) {
  lines(c(i,i), c(fit.conf.lb[i], fit.conf.ub[i]))
}
legend( "topleft", legend = c("true value", "estimate", "confidence interval"),
        col = c("red", "blue", 1), pch = c(1, 7, NA), lty = c(NA, NA, 1), inset = 0.01)

dev.off()
par(mfrow=c(2, 3))
par(mar = c(2, 2, 2, 2))
plot(s_garchFit,which="all")
dev.off()

#gjrGARCH
gjr_garchMod <- ugarchspec(
  variance.model=list(model="gjrGARCH",
                      garchOrder=c(1,1), external.regressors = df_dummies ),
  mean.model=list(armaOrder=c(1,0),
                  include.mean=TRUE
  ), 
  distribution.model="norm"
)
gjr_garchFit <- ugarchfit(spec=gjr_garchMod, data=log_ret)
coef(gjr_garchFit)
gjr_garchFit
gjr_rhat <- gjr_garchFit@fit$fitted.values
plot.ts(gjr_rhat)
gjr_hhat <- ts(gjr_garchFit@fit$sigma^2)
plot.ts(gjr_hhat)
dev.off()
par(mfrow=c(2, 3))
par(mar = c(2, 2, 2, 2))
plot(gjr_garchFit,which="all")
dev.off()

#eGARCH
e_garchMod <- ugarchspec(
  variance.model=list(model="eGARCH",
                      garchOrder=c(1,1), external.regressors = df_dummies ),
  mean.model=list(armaOrder=c(1,0),
                  include.mean=TRUE
  ), 
  distribution.model="norm"
)
e_garchFit <- ugarchfit(spec=e_garchMod, data=log_ret)
coef(e_garchFit)
e_garchFit
e_rhat <- e_garchFit@fit$fitted.values
plot.ts(e_rhat)
e_hhat <- ts(e_garchFit@fit$sigma^2)
plot.ts(e_hhat)
dev.off()
par(mfrow=c(2, 3))
par(mar = c(2, 2, 2, 2))
plot(e_garchFit,which="all")
dev.off()
# apARCH
apARCHMod <- ugarchspec(variance.model=list(model="apARCH",
                                             garchOrder=c(1,1),
                                            external.regressors = df_dummies ),
                         mean.model=list(armaOrder=c(1,0)), 
                         distribution.model="norm")
apARCHFit <- ugarchfit(spec=apARCHMod, data= log_ret)
coef(apARCHFit)
apARCHFit
apARCH_rhat <- apARCHFit@fit$fitted.values
plot.ts(apARCH_rhat)
apARCH_hhat <- ts(apARCHFit@fit$sigma^2)
plot.ts(apARCH_hhat)
dev.off()
par(mfrow=c(2, 3))
par(mar = c(2, 2, 2, 2))
plot(apARCHFit,which="all")
dev.off()

# IGARCH 
i_garchMod <- ugarchspec(variance.model=list(model="iGARCH",
                                              garchOrder=c(1,1),
                                              external.regressors = df_dummies ),
                          mean.model=list(armaOrder=c(1,0)), 
                          distribution.model="norm")
i_garchFit <- ugarchfit(spec=i_garchMod, data= log_ret)
coef(i_garchFit)
i_garchFit
i_rhat <- i_garchFit@fit$fitted.values
plot.ts(i_rhat)
i_hhat <- ts(i_garchFit@fit$sigma^2)
plot.ts(i_hhat)
dev.off()
par(mfrow=c(2, 3))
par(mar = c(2, 2, 2, 2))
plot(i_garchFit,which="all")
dev.off()

# GARCH MODELS assuming a Student's t-distribution
# Student's t-sGARCH --->  tsGARCH
library(rugarch)

ts_garchMod <- ugarchspec(mean.model = list(armaOrder = c(1, 0), include.mean = TRUE
), 
variance.model = list(model = 'sGARCH', 
                      garchOrder = c(1, 1),
                      external.regressors = df_dummies
),
distribution.model = "std")

ts_garchFit <- ugarchfit(spec=ts_garchMod, data=log_ret)
ts_garchFit
## Results review 1
coef(ts_garchFit)

ts_rhat <- ts_garchFit@fit$fitted.values
plot.ts(ts_rhat)
ts_hhat <- ts(ts_garchFit@fit$sigma^2)
plot.ts(ts_hhat)
dev.off()
## Results review 2
tfit.val     <- coef(ts_garchFit)
tfit.sd      <- diag(vcov(ts_garchFit))
ttrue.val = ts_garchFit@fit$tval

tfit.conf.lb <- tfit.val + qnorm(0.025) * tfit.sd
tfit.conf.ub <- tfit.val + qnorm(0.975) * tfit.sd
print(tfit.val)

print(tfit.sd)

print(ttrue.val)
plot(ttrue.val, pch = 1, col = "red",
     ylim = range(c(tfit.conf.lb, tfit.conf.ub, ttrue.val)),
     xlab = "", ylab = "", axes = TRUE)
box(); axis(1, at = 1:length(tfit.val), labels = names(tfit.val)); axis(2)
points(coef(ts_garchFit), col = "blue", pch = 7)
for (i in 1:length(tfit.val)) {
  lines(c(i,i), c(tfit.conf.lb[i], tfit.conf.ub[i]))
}
legend( "topleft", legend = c("true value", "estimate", "confidence interval"),
        col = c("red", "blue", 1), pch = c(1, 7, NA), lty = c(NA, NA, 1), inset = 0.01)
dev.off()
par(mfrow=c(2, 3))
par(mar = c(2, 2, 2, 2))
plot(ts_garchFit,which="all")
dev.off()

#Student's t-gjrGARCH --->tgjrGARCH
tgjr_garchMod <- ugarchspec(
  variance.model=list(model="gjrGARCH",
                      garchOrder=c(1,1), external.regressors = df_dummies ),
  mean.model=list(armaOrder=c(1,0),
                  include.mean=TRUE
  ), 
  distribution.model="std"
)
tgjr_garchFit <- ugarchfit(spec=tgjr_garchMod, data=log_ret)
coef(tgjr_garchFit)
tgjr_garchFit
tgjr_rhat <- tgjr_garchFit@fit$fitted.values
plot.ts(tgjr_rhat)
tgjr_hhat <- ts(tgjr_garchFit@fit$sigma^2)
plot.ts(tgjr_hhat)
dev.off()
par(mfrow=c(2, 3))
par(mar = c(2, 2, 2, 2))
plot(tgjr_garchFit,which="all")
dev.off()

#Student's t-eGARCH --->teGARCH
te_garchMod <- ugarchspec(
  variance.model=list(model="eGARCH",
                      garchOrder=c(1,1), external.regressors = df_dummies ),
  mean.model=list(armaOrder=c(1,0),
                  include.mean=TRUE
  ), 
  distribution.model="std"
)
te_garchFit <- ugarchfit(spec=te_garchMod, data=log_ret)
coef(te_garchFit)
te_garchFit
te_rhat <- te_garchFit@fit$fitted.values
plot.ts(te_rhat)
te_hhat <- ts(te_garchFit@fit$sigma^2)
plot.ts(te_hhat)
dev.off()
par(mfrow=c(2, 3))
par(mar = c(2, 2, 2, 2))
plot(te_garchFit,which="all")
dev.off()

#Student's t-apARCH ---> tapARCH 
tapARCHMod <- ugarchspec(variance.model=list(model="apARCH",
                                            garchOrder=c(1,1),
                                            external.regressors = df_dummies ),
                        mean.model=list(armaOrder=c(1,0)), 
                        distribution.model="std")
tapARCHFit <- ugarchfit(spec=tapARCHMod, data= log_ret)
coef(tapARCHFit)
tapARCHFit
tapARCH_rhat <- tapARCHFit@fit$fitted.values
plot.ts(tapARCH_rhat)
tapARCH_hhat <- ts(tapARCHFit@fit$sigma^2)
plot.ts(tapARCH_hhat)
dev.off()
par(mfrow=c(2, 3))
par(mar = c(2, 2, 2, 2))
plot(tapARCHFit,which="all")
dev.off()

#Student's t-tiGARCH ---> tiGARCH 
ti_garchMod <- ugarchspec(variance.model=list(model="iGARCH",
                                             garchOrder=c(1,1),
                                             external.regressors = df_dummies ),
                         mean.model=list(armaOrder=c(1,0), include.mean=TRUE),
                         distribution.model="std")
ti_garchFit <- ugarchfit(spec=ti_garchMod, data= log_ret)
coef(ti_garchFit)
ti_garchFit
ti_rhat <- ti_garchFit@fit$fitted.values
plot.ts(ti_rhat)
ti_hhat <- ts(ti_garchFit@fit$sigma^2)
plot.ts(ti_hhat)
dev.off()
par(mfrow=c(2, 3))
par(mar = c(2, 2, 2, 2))
plot(ti_garchFit,which="all")
dev.off()

#Choose best GARCH model
#Normal
AIC= c(infocriteria(s_garchFit)[1], infocriteria(gjr_garchFit)[1], 
       infocriteria(e_garchFit)[1], infocriteria(apARCHFit)[1],  infocriteria(i_garchFit)[1] )
rowSAIC= c("s_garch", "gjr_garchFit", "e_garchFit", "apARCHFit", "i_garchFit")


AIC_Results= data.frame(AIC, row.names = rowSAIC)
AIC_Results
BestAICNorm=AIC_Results[which.min(AIC_Results$AIC),]
BestAICNorm

BIC= c(infocriteria(s_garchFit)[2], infocriteria(gjr_garchFit)[2], 
       infocriteria(e_garchFit)[2], infocriteria(apARCHFit)[2],  infocriteria(i_garchFit)[2] )
rowSBIC= c("s_garch", "gjr_garchFit", "e_garchFit", "apARCHFit", "i_garchFit")
BIC_Results= data.frame(BIC, row.names = rowSBIC)
BIC_Results
BestBICNorm=BIC_Results[which.min(BIC_Results$BIC),]
BestBICNorm

Best_Norm=c(BestAICNorm,BestBICNorm)
rowSBest_Norm= c("AICbest","BICbest")
Norm_Results= data.frame(Best_Norm, row.names =rowSBest_Norm)
Norm_Results

#Student-t
tAIC= c(infocriteria(ts_garchFit)[1], infocriteria(tgjr_garchFit)[1], 
       infocriteria(te_garchFit)[1], infocriteria(tapARCHFit)[1],  infocriteria(ti_garchFit)[1] )
trowSAIC= c("ts_garch", "tgjr_garchFit", "te_garchFit", "tapARCHFit", "ti_garchFit")


tAIC_Results= data.frame(tAIC, row.names = trowSAIC)
tAIC_Results
tAIC_Results[which.min(tAIC_Results$tAIC),]


tBIC= c(infocriteria(ts_garchFit)[2], infocriteria(tgjr_garchFit)[2], 
       infocriteria(te_garchFit)[2], infocriteria(tapARCHFit)[2],  infocriteria(ti_garchFit)[2]  )
trowSBIC= c("ts_garch", "tgjr_garchFit", "te_garchFit", "tapARCHFit", "ti_garchFit")
tBIC_Results= data.frame(tBIC, row.names = trowSBIC)
tBIC_Results
tBIC_Results[which.min(tBIC_Results$tBIC),]

#Goodness of Fit
Distribution=c("Normal", "Normal", "Normal", "Normal", "Normal",
               "Student's-t", "Student's-t", "Student's-t", "Student's-t", "Student's-t" )
rowsGARCH=c("s_garchFit", "gjr_garchFit","e_garchFit","apARCHFit", "i_garchFit",
            "ts_garchFit", "tgjr_garchFit", "te_garchFit", "tapARCHFit", "ti_garchFit")
Best=c("o", "o", "*", "o", "o",
       "o", "o", "*", "o","o")
final_AIC=c(AIC,tAIC)
final_BIC=c(BIC,tBIC)
Results <- data.frame("AIC" = final_AIC, "BIC" = final_BIC , "Distribution" = Distribution, 
                      "Best"= Best ,row.names = rowsGARCH)
Results


#Validation of the model
residual_bestModel= residuals(te_garchFit, standardize = TRUE)
Box.test(abs(residual_bestModel), 10, type = "Ljung-Box")
acf(abs(residual_bestModel), 10)
