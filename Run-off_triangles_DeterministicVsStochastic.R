#BATTLE CREEK MUT INS CO
############################################################################################
# 					Development Triangles							 #
############################################################################################

n<-10
Claims <- data.frame(originf = factor(rep(1988:1997, n:1)),
                     dev=sequence(n:1),
                     inc.paid=
                       c(1087,1103,1133,1136,1140,1150,1142,1142,1142,1142,
                         1136,1200,1229,1214,1273,1264,1270,1270,1270,
                         1411,1464,1524,1532,1577,1557,1558,1546,
                         1164,1146,1175,1161,1161,1162,1165,
                         920,867,854,838,833,843,
                         745,734,731,729,732,
                         666,718,739,747,
                         737,808,769,
                         530,699,
                         766))

(inc.triangle <- with(Claims, {
  M <- matrix(nrow=n, ncol=n,
              dimnames=list(origin=levels(originf), dev=1:n))
  M[cbind(originf, dev)] <- inc.paid
  M
}))

(cum.triangle <- t(apply(inc.triangle, 1, cumsum)))

(latest.paid <- cum.triangle[row(cum.triangle) == n - col(cum.triangle) + 1])

Claims$cum.paid <- cum.triangle[with(Claims, cbind(originf, dev))]

op <- par(fig=c(0,0.5,0,1), cex=0.8, oma=c(0,0,0,0))

with(Claims, {
  interaction.plot(x.factor=dev, trace.factor=originf, response=inc.paid,
                   fun=sum, type="b", bty="n", legend=FALSE); axis(1, at=1:n)
  par(fig=c(0.45,1,0,1), new=TRUE, cex=0.8, oma=c(0,0,0,0))
  interaction.plot(x.factor=dev, trace.factor=originf, response=cum.paid,
                   fun=sum, type="b", bty="n"); axis(1,at=1:n)
})

mtext("Incremental and cumulative claims development",
      side=3, outer=TRUE, line=-3, cex = 1.1, font=2)

par(op)

library(lattice)

xyplot(cum.paid ~ dev | originf, data=Claims, t="b", layout=c(6,2),
       as.table=TRUE, main="Cumulative claims development")

############################################################################################
# 					Deterministic Reserving Methods					 #
# 					    Chain-Ladder Algorithm						 #
############################################################################################

f <- sapply((n-1):1, function(i) {
  sum( cum.triangle[1:i, n-i+1] ) / sum( cum.triangle[1:i, n-i] )
})

tail <- 1

(f <- c(f, tail))

full.triangle <- cum.triangle

for(k in 1:(n-1)){
  full.triangle[(n-k+1):n, k+1] <- full.triangle[(n-k+1):n,k]*f[k]
}

full.triangle

(ultimate.paid <- full.triangle[,n])

(ldf <- rev(cumprod(rev(f))))

(dev.pattern <- 1/ldf)

(reserve <- sum (latest.paid * (ldf - 1)))
#  otherwise
#sum(ultimate.paid - latest.paid)

a <- ultimate.paid

(b <- c(dev.pattern[1], diff(dev.pattern)))

(X.hat <- a %*% t(b))
#Suppose the expected loss cost for the 1997 origin year is 20,000
(BF1997 <- ultimate.paid[n] * dev.pattern[1] + 20000 * (1 - dev.pattern[1]))

## Tail Factors ##

dat <- data.frame(lf1=log(f[-c(1,n)]-1), dev=2:(n-1))

(m <- lm(lf1 ~ dev , data=dat))

plot(lf1 ~ dev, main="log(f - 1) ~ dev", data=dat, bty="n")

abline(m)

sigma <- summary(m)$sigma

extrapolation <- predict(m, data.frame(dev=n:100))

(tail <- prod(exp(extrapolation + 0.5*sigma^2) + 1))

library(ChainLadder)

ata(cum.triangle)

############################################################################################
#						Stochastic Reserving Models					 #
#				Chain-Ladder in the Context of Linear Regression			 #
############################################################################################

names(Claims)[3:4] <- c("inc.paid.k", "cum.paid.k")

ids <- with(Claims, cbind(originf, dev))

Claims <- within(Claims,{
  cum.paid.kp1 <- cbind(cum.triangle[,-1], NA)[ids]
  inc.paid.kp1 <- cbind(inc.triangle[,-1], NA)[ids]
  devf <- factor(dev)
}
)

delta <- 0:2

ATA <- sapply(delta, function(d)
  coef(lm(cum.paid.kp1 ~ 0 + cum.paid.k : devf,
          weights=1/cum.paid.k^d, data=Claims))
)

dimnames(ATA)[[2]] <- paste("Delta = ", delta)

ATA

xyplot(cum.paid.kp1 ~ cum.paid.k | devf,
       data=subset(Claims, dev < (n-1)),
       main="Age-to-age developments", as.table=TRUE,
       scales=list(relation="free"),
       key=list(columns=2, lines=list(lty=1:4, type="l"),
                text=list(lab=c("lm(y ~ x)",
                                "lm(y ~ 0 + x)",
                                "lm(y ~ 0 + x, w=1/x)",
                                "lm(y ~ 0 + x, w=1/x^2)"))),
       panel=function(x,y,...){
         panel.xyplot(x,y,...)
         if(length(x)>1){
           panel.abline(lm(y ~ x), lty=1)
           panel.abline(lm(y ~ 0 + x), lty=2)
           panel.abline(lm(y ~ 0 + x, weights=1/x), lty=3)
           panel.abline(lm(y ~ 0 + x, weights=1/x^2), lty=4)
         }
       }
)

############################################################################################
# 							Mack Model							 #
############################################################################################

library(ChainLadder)

(mack <- MackChainLadder(cum.triangle, weights=1, alpha=1,est.sigma="Mack"))

plot(mack, lattice=TRUE, layout=c(6,2))

plot(mack)
#GLM
preg <- glm(inc.paid.k ~ originf + devf,
            data=Claims, family=poisson(link = "log"))

summary(preg)

allClaims <- data.frame(origin = sort(rep(1988:1997, n)),
                        dev = rep(1:n,n))

allClaims <- within(allClaims, {
  devf <- factor(dev)
  cal <- origin + dev - 1
  originf <- factor(origin)
})

(pred.inc.tri <- t(matrix(predict(preg,type="response",
                                  newdata=allClaims), n, n)))

sum(predict(preg,type="response", newdata=subset(allClaims, cal > 1997)))

df <- c(0, coef(preg)[(n+1):(2*n-1)])

sapply(2:7, function(i) sum(exp(df[1:i]))/sum(exp(df[1:(i-1)])))

library(AER)

dispersiontest(preg)

## Quantifying Uncertainty in GLMs  ##

summary(odpreg <- glm(inc.paid.k ~ originf + devf, data=Claims,
                      family=quasipoisson))

mu.hat <- predict(odpreg, newdata=allClaims, type="response")*(allClaims$cal>1997)
phi <- summary(odpreg)$dispersion
Sigma <- vcov(odpreg)
model.formula <- as.formula(paste("~", formula(odpreg)[3]))
# Future design matrix
X <- model.matrix(model.formula, data=allClaims)
Cov.eta <- X%*% Sigma %*%t(X)
sqrt(phi * sum(mu.hat) + t(mu.hat) %*% Cov.eta %*% mu.hat)
op <- par(mfrow=c(2,2), oma = c(0, 0, 3, 0))
par(mar=c(2,2,2,2))
plot(preg)
par(op)
(odp <- glmReserve(as.triangle(inc.triangle), var.power=1, cum=FALSE)) 

############################################################################################
# 						Bootstrap Chain-Ladder						 #
############################################################################################

set.seed(1)

(B <- BootChainLadder(cum.triangle, R=1000, process.distr="od.pois"))

plot(B)

quantile(B, c(0.75,0.95,0.99, 0.995))

library(fitdistrplus)

(fit <- fitdist(B$IBNR.Totals[B$IBNR.Totals>0], "lnorm"))

plot(fit)

qlnorm(0.995, fit$estimate["meanlog"], fit$estimate["sdlog"])

ny <- (col(inc.triangle) == (nrow(inc.triangle) - row(inc.triangle) + 2))
paid.ny <- apply(B$IBNR.Triangles, 3,
                 function(x){
                   next.year.paid <- x[col(x) == (nrow(x) - row(x) + 2)]
                   sum(next.year.paid)
                 })
paid.ny.995 <- B$IBNR.Triangles[,,order(paid.ny)[round(B$R*0.995)]]
inc.triangle.ny <- inc.triangle
(inc.triangle.ny[ny] <- paid.ny.995[ny])

##  Reserving Based on Log-Incremental Payments  ##

Claims <- within(Claims, {
  log.inc <- log(inc.paid.k)
  cal <- as.numeric(levels(originf))[originf] + dev - 1
})

Claims <- within(Claims, {
  d1 <- ifelse(dev < 2, 1, 0)
  d27 <- ifelse(dev < 2, 0, dev - 1)
})

summary(fit1 <- lm(log.inc ~ originf + d1 + d27, data=Claims))

Claims <- within(Claims, {
  a6 <- ifelse(originf == 1996, 1, 0)
  a7 <- ifelse(originf == 1997, 1, 0)
})
summary(fit2 <- lm(log.inc ~ a6 + a7 + d1 + d27, data=Claims))

op <- par(mfrow=c(2,2), oma = c(0, 0, 3, 0))
plot(fit2)
par(op)

shapiro.test(fit2$residuals)

resPlot <- function(model, data){
  xvals <- list(
    fitted = model[["fitted.values"]],
    origin = as.numeric(levels(data$originf))[data$originf],
    cal=data$cal, dev=data$dev
  )
  op <- par(mfrow=c(2,2), oma = c(0, 0, 3, 0))
  for(i in 1:4){
    plot.default(rstandard(model) ~ xvals[[i]] ,
                 main=paste("Residuals vs", names(xvals)[i] ),
                 xlab=names(xvals)[i], ylab="Standardized residuals")
    panel.smooth(y=rstandard(model), x=xvals[[i]])
    abline(h=0, lty=2)
  }
  mtext(as.character(model$call)[2], outer = TRUE, cex = 1.2)
  par(op)
}

resPlot(fit2, Claims)

Claims <- within(Claims, {
  p34 <- ifelse(cal < 1995 & cal > 1989, cal-1989, 0)
})

summary(fit3 <- update(fit2, ~ . + p34, data=Claims))

resPlot(fit3, Claims)

log.incr.predict <- function(model, newdata){
  Pred <- predict(model, newdata=newdata, se.fit=TRUE)
  Y <- Pred$fit
  VarY <- Pred$se.fit^2 + Pred$residual.scale^2
  P <- exp(Y + VarY/2)
  VarP <- P^2*(exp(VarY)-1)
  seP <- sqrt(VarP)
  model.formula <- as.formula(paste("~", formula(model)[3]))
  mframe <- model.frame(model.formula, data=newdata)
  X <- model.matrix(model.formula, data=newdata)
  varcovar <- X %*% vcov(model) %*% t(X)
  CoVar <- sweep(sweep((exp(varcovar)-1), 1, P, "*"), 2, P, "*")
  CoVar[col(CoVar)==row(CoVar)] <- 0
  Total.SE <- sqrt(sum(CoVar) + sum(VarP))
  Total.Reserve <- sum(P)
  Incr=data.frame(newdata, Y, VarY, P, seP, CV=seP/P)
  out <- list(Forecast=Incr,
              Totals=data.frame(Total.Reserve,
                                Total.SE=Total.SE,
                                CV=Total.SE/Total.Reserve))
  return(out)
}

tail.years <-9

fdat <- data.frame(
  origin=rep(1988:1997, n+tail.years),
  dev=rep(1:(n+tail.years), each=n)
)

fdat <- within(fdat, {
  cal <- origin + dev - 1
  a7 <- ifelse(origin == 1997, 1, 0)
  a6 <- ifelse(origin == 1996, 1, 0)
  originf <- factor(origin)
  p34 <- ifelse(cal < 1995 & cal > 1989, cal-1989, 0)
  d1 <- ifelse(dev < 2, 1, 0)
  d27 <- ifelse(dev < 2, 0, dev - 1)
})

reserve2 <- log.incr.predict(fit2, subset(fdat, cal>1997))
reserve2$Totals

reserve3 <- log.incr.predict(fit3, subset(fdat, cal>1997))
reserve3$Totals

round(xtabs(P ~ origin + dev, reserve3$Forecast))

round(summary(MackChainLadder(cum.triangle, est.sigma="Mack",
                              tail=1.05, tail.se=0.02))$Totals,2)


Incurred= cum.triangle
n<-10
Claimspaid <- data.frame(originf = factor(rep(1988:1997, n:1)),
                         dev=sequence(n:1),
                         cum.paid=
                           c(705,978,1066,1098,1116,1147,1142,1142,1142,1142,
                             695,1046,1122,1179,1242,1249,1270,1270,1270,
                             849,1166,1320,1385,1518,1528,1574,1546,
                             679,903,1044,1090,1148,1150,1165,
                             553,741,790,810,831,843,
                             443,597,645,696,722,
                             397,556,655,699,
                             439,626,702,
                             388,504,
                             480))

(cumpaid.triangle <- with(Claimspaid, {
  M <- matrix(nrow=n, ncol=n,
              dimnames=list(origin=levels(originf), dev=1:n))
  M[cbind(originf, dev)] <- cum.paid
  M
}))

(cum.paid.triangle <- t(apply(cumpaid.triangle, 1, cumsum)))


Paid=cum.paid.triangle

MCL=MunichChainLadder(Paid=Paid ,Incurred= Incurred, 
                      est.sigmaP = "log-linear", est.sigmaI = "log-linear", 
                      tailP=FALSE, tailI=FALSE)
MCL
plot(MCL)