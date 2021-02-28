#####################################
##### PLAT model
#####################################
f2 <- function(x, ages) mean(ages) - x
constPlat <- function(ax, bx, kt, b0x, gc, wxt, ages){
 nYears <- dim(wxt)[2]
 x <- ages
 t <- 1:nYears
 c <- (1 - tail(ages, 1)):(nYears - ages[1])
 xbar <- mean(x)a
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
#####################################

####################################
#### EXTRAPOLATION qxt UP TO EXTREME AGE VIA EXTRAPOLATION
####################################
extrapolation.fit <- function(q){
logit_q <- log(q/(1-q))
#y = logit_q = ax+b
a1_log <- min((a.max-9),90) # starting age linear model
a2_log <- a.max # ending age linear model
l_vet1 <- c((a1_log-a.min+1):(a2_log-a.min+1))
lm(logit_q[l_vet1,]~A.fit[l_vet1])
B <- lm(logit_q[l_vet1,]~A.fit[l_vet1])$coeff[1,]
A <- lm(logit_q[l_vet1,]~A.fit[l_vet1])$coeff[2,]

a1_ext <- min((a.max-4),100) #starting age extrapolation
a2_ext <- omega-1 # ending age extrapolation
Age1 <- c(a1_ext:a2_ext)
Age2 <- c(a.min:a2_ext)
q_tail=array(0,dim=c(length(Age1),(length(Y.fit)+y.pred)), dimnames=list(Age1,c(Y.fit,Y.pred)))
for (j in Age1){
    q_tail[(j-Age1[1]+1),] <- exp(A*j+B)/(1+exp(A*j+B))
}
l_vet3 <- c(1:(a1_ext-a.min))
l_vet4 <- c(1:(a2_ext-a.min+1))
q.ext=array(0,dim=c(length(c(a.min:a2_ext)),(length(Y.fit)+y.pred)),
   dimnames=list(c(a.min:a2_ext),c(Y.fit,Y.pred)))
for (j in 1:(length(Y.fit)+y.pred)){
  q.ext[,j] <- c(q[l_vet3,j],q_tail[,j])
}
q.ext[(omega-a.min),]=1
return(q.ext)
}
#####################################

#####################################
#### EXTRAPOLATION SIMULATED qxt UP TO EXTREME AGE VIA EXTRAPOLATION
#####################################
extrapolation.sim <- function(q.st){
logit_q.st <- log(q.st/(1-q.st))
#y = logit_q.st = ax+b
B <- matrix(0, nrow=n.sim, ncol=y.pred)
A <- matrix(0, nrow=n.sim, ncol=y.pred)
a1_log <- min((a.max-9),90) # starting age graduation
a2_log <- a.max # ending age graduation
l_vet1 <- c((a1_log-a.min+1):(a2_log-a.min+1))
for (k in 1:n.sim){
  lm(logit_q.st[l_vet1,,k]~A.fit[l_vet1])
  B[k,] <- lm(logit_q.st[l_vet1,,k]~A.fit[l_vet1])$coeff[1,1:y.pred]
  A[k,] <- lm(logit_q.st[l_vet1,,k]~A.fit[l_vet1])$coeff[2,1:y.pred]
}

a1_ext <- min((a.max-4),100) #starting age extrapolation
a2_ext <- omega-1 # ending age extrapolation
Age1 <- c(a1_ext:a2_ext)
Age2 <- c(a.min:a2_ext)
q_tail.st=array(0,dim=c(length(Age1),y.pred,n.sim), dimnames=list(Age1,Y.pred,
  c(1:n.sim)))
for (j in Age1){
  for (k in 1:n.sim){
    q_tail.st[(j-Age1[1]+1),,k] <- exp(A[k,]*j+B[k,])/(1+exp(A[k,]*j+B[k,]))
  }
}
l_vet3 <- c(1:(a1_ext-a.min))
l_vet4 <- c(1:(a2_ext-a.min+1))
q.st.ext=array(0,dim=c(length(c(a.min:a2_ext)),y.pred,n.sim),
   dimnames=list(c(a.min:a2_ext),Y.pred,c(1:n.sim)))
for (j in 1:y.pred){
  for (k in 1:n.sim){
  q.st.ext[,j,k] <- c(q.st[l_vet3,j,k],q_tail.st[,j,k])
  }
}
q.st.ext[(omega-a.min),,]=1
return(q.st.ext)
}
#####################################

####################################
######## Life expectancy
####################################
life.exp <- function(q){
p <- 1-q # 1 year survival probability
p0n <- apply(p, 2, cumprod) # n year survival probability
ex <- p0n
ex[1,] <- apply(p0n, 2, sum)+0.5
ex[2,] <- (apply(p0n, 2, sum)-p0n[1,])/p0n[1,]+0.5
for (j in 3:(omega-a.min)){
      ex[j,]<-(apply(p0n, 2, sum)-apply(p0n[1:(j-1),],2,sum))/p0n[(j-1),]+0.5
}
return(ex)
}

####################################
##### Present value of annuities from a matrix of death probabilities q(x,t)
####################################
annuity <- function(q,v_vect){
#q <- q_ITAm.ext
p <- 1-q
p0n <- apply(p,2,cumprod)
E0n <-  p0n*v_vect
ann0 <-  apply(E0n, 2, cumsum) # matrix of present values of annuity starting from age 0
annx <- array(0,dim=c(nrow(q),ncol(q)), dimnames=list(rownames(q),colnames(q)))
annx[1,] <- ann0[omega-a.min,]
for (j in 2:(omega-a.min)){
      annx[j,] <- (ann0[omega-a.min,]-ann0[j-1,])/E0n[j-1,] # matrix of present values of annuity starting from age x
}
return(annx)
}
####################################


####################################
##### Present value of annuities from N. sim matrices of death probabilities q(x,t)
####################################
annuity.st <- function(q.st,v_vect, conf.lev){
p.st <- 1-q.st # 1 year survival probability
p0n.st <- apply(p.st, c(2,3), cumprod) # n year survival probability from age 0
E0n.st <-  p0n.st*v_vect
ann0.st <-  apply(E0n.st, c(2,3), cumsum) # matrix of present values of annuity starting from age 0
annx.st <- array(0,dim=c((omega-a.min),y.pred,n.sim), dimnames=list(c(a.min:(omega-1)),Y.pred,c(1:n.sim)))
annx.st[1,,] <- ann0.st[omega-a.min,,]
for (j in 2:(omega-a.min)){
      annx.st[j,,] <- (ann0.st[omega-a.min,,]-ann0.st[j-1,,])/E0n.st[j-1,,] # matrix of present values of annuity starting from age x
}
ann.mean <- apply(annx.st, c(1,2),mean) # expected present value of annuity starting form age x
ann.q95 <- matrix(0, nrow=(omega-a.min),ncol=y.pred, dimnames=list(c(a.min:(omega-1)), c((tail(Y.fit,1)+1):(tail(Y.fit,1)+y.pred))))
ann.q05 <- matrix(0, nrow=(omega-a.min),ncol=y.pred, dimnames=list(c(a.min:(omega-1)), c((tail(Y.fit,1)+1):(tail(Y.fit,1)+y.pred))))
for (j in 1:(omega-a.min)){
  for (k in 1:y.pred){
      ann.q95[j,k] <- quantile(annx.st[j,k,], probs=0.5+conf.lev/2)
      ann.q05[j,k] <- quantile(annx.st[j,k,], probs=(1-conf.lev)/2)
  }
}
return(list(ann.mean,ann.q95,ann.q05))
}
####################################

