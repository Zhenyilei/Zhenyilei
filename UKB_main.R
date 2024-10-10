# Set working directory
rm(list = ls())
#setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(quantreg)
library(VGAM)
data1 <- read.csv("C:/Users/leizh/OneDrive/Desktop/UKB_HCPA_FA/UKB_FA_dat.csv")


hist(data1$Age_Instance2)
bn<-5.5
x0<-50
plot(data1$Age_Instance2,data1$ACR.L)
tt<-rq(ACR.L~Age_Instance2,data=data1,tau = 0.9)
mt<-summary(tt)
mt$coefficients[2,1]

test_data <- data1[(-bn <= data1$Age - x0 & data1$Age - x0 <= bn), ]
test_data$Age=(test_data$Age_Instance2 - x0) / bn
tt<-rq(ACR.L~Age_Instance2,data=test_data,tau = 0.9)
tt1<-rq(ACR.L~Age_Instance2,data=test_data,tau = 0.9)
summary(tt1)
mt<-summary(tt)

plot(test_data$Age_Instance2,test_data$ACR.L)

new_data <- data.frame(Age_Instance2 = 50)
predicted_value <- predict(tt, newdata = new_data)
abline(a=mt$coefficients[1,1],b=bn*mt$coefficients[2,1],col="red",lwd=2)
##### Kernel #####

Epanechnikov <- function(u) {
  return(3 / 4 * (1 - u ^ 2) * (-1 <= u & u <= 1))
}

Normal.Kernel <- function(u, k = 2) {
  return(dnorm(u, mean = 0, sd = 1 / k) * (pnorm(k) - pnorm(-k)))
}



##### Local linear quantile regression #####
LLQR <- function(bn, x0, x, y,tau=0.9,Kernel=Epanechnikov) {
  data <- data.frame(Age = x, Test = y)
  ## draw the data within the bandwidth
  test_data <- data[(abs(data$Age - x0) <= bn), ]
  ## calculate weights 
  weights <- Kernel((test_data$Age - x0) / bn)
  ## standradlize the data
  test_data$Standrad_Age<-(test_data$Age - x0) / bn
  ## linear quantile regression
  rq_temp <- rq(Test ~ Standrad_Age, data = test_data, weight = weights, tau = tau)
  beta_bn_beta_hat <- as.vector(rq_temp$coefficients)
  
  res<-rq_temp$residuals
  local_n<-nrow(test_data)
  ## get the density of residuals
  dens <- density(res)
  ## get the index of 90% quantile of residuals 
  idx <- which.min(abs(dens$x - quantile(res,tau)))
  ## get the density of value of 90% quantile
  u_t<-dens$y[idx]
  return(list(beta_0 = beta_bn_beta_hat[1],
              bn_beta_1 = beta_bn_beta_hat[2],
              u_t=u_t,local_n=local_n))
}


##### Get the parameter for given bn #####
Get_beta<-function(bn,x,y,N_set,tau=0.9){
  N<- length(N_set)
  beta_df <- data.frame(
    Age = N_set,beta_hat_0 = numeric(N),bn_beta_hat_1 = numeric(N),
    beta_check_0 = numeric(N),bn_beta_check_1 = numeric(N),
    u_t=numeric(N),local_n=numeric(N)
  )
  
  for (i in seq_along(N_set)) {
    beta_hat<- LLQR(bn, N_set[i] , x, y,tau=tau)
    beta_check<-LLQR(sqrt(2)*bn,N_set[i],x,y,tau=tau)
    
    beta_df$beta_hat_0[i] <- beta_hat$beta_0
    beta_df$bn_beta_hat_1[i] <- beta_hat$bn_beta_1
    beta_df$beta_check_0[i] <- beta_check$beta_0
    beta_df$bn_beta_check_1[i] <- beta_check$bn_beta_1
    beta_df$u_t[i]<-beta_hat$u_t
    beta_df$local_n[i]<-beta_hat$local_n
  }
  return(beta_df)
}

##### calculate the Hill estimator #####
getHnhat <- function(X, alphan = 0.1) {
  n <- length(X)
  
  kn <- floor(n * alphan)
  
  X_sorted <- sort(X, decreasing = TRUE)
  
  H_hat <- sum(log(X_sorted[1:kn]) - log(X_sorted[kn + 1])) / kn
  
  return(H_hat)
}



N_set<-seq(45,82,by=0.025)
N<-length(N_set)
temp_beta<-Get_beta(bn=5.5,data1$Age_Instance2,data1$ACR.L,N_set,tau=0.9)


res_list <- list()
##Get_local_residuals(5.5, 50,temp_beta$beta_hat_0[201],temp_beta$bn_beta_hat_1[201], data1$Age_Instance2, data1$ACR.L,tau=0.5)$local_res
for (i in 1:38) {
  idd = which(data1$Age_Instance2 == (i+44))
  idx<-which(temp_beta$Age==i+44)
  id_data<-which(data1$Age_Instance2==i+44)
  
  res_list[[i]] <- data1$ACR.L[idd]-temp_beta$beta_hat_0[idx]
}
res <- unname(unlist(res_list))


Hhat=getHnhat(abs(res))


#### Simulation for T_n ####
Simulation1<-function(num,Hhat,tau=0.9,bn=5.5){
  result<-matrix(,nrow = N,ncol = num)
  true_beta<-qfrechet(p = tau, shape = 1/Hhat)
  for (i in 1:num) {
    set.seed(i)
    print(i)
    Sim<-rfrechet(40414, shape = 1/Hhat)
    temp_beta<-Get_beta(bn,data1$Age_Instance2,Sim,N_set,tau=tau)
    beta_dif<-abs(2*temp_beta$beta_hat_0-temp_beta$beta_check_0-true_beta)
    result[,i]<-(temp_beta$u_t*sqrt(temp_beta$local_n * bn * (tau)))*beta_dif
  }
  thres<-apply(result,1, quantile,probs=0.95)
  return(thres)  
}


thres<- Simulation1(150,Hhat,tau=0.9,bn=5.5)

write.csv(thres, "thres90.csv", row.names = FALSE)


#### 95
N_set<-seq(45,82,by=0.025)
N<-length(N_set)
temp_beta<-Get_beta(bn=5.5,data1$Age_Instance2,data1$ACR.L,N_set,tau=0.95)


res_list <- list()
for (i in 1:38) {
  idd = which(data1$Age_Instance2 == (i+44))
  idx<-which(temp_beta$Age==i+44)
  id_data<-which(data1$Age_Instance2==i+44)
  
  res_list[[i]] <- data1$ACR.L[idd]-temp_beta$beta_hat_0[idx]
}
res <- unname(unlist(res_list))

Hhat=getHnhat(abs(res))


thres_95<- Simulation1(150,Hhat,tau=0.95,bn=5.5)
write.csv(thres_95, "thres95.csv", row.names = FALSE)



#### 99
N_set<-seq(45,82,by=0.025)
N<-length(N_set)
temp_beta<-Get_beta(bn=5.5,data1$Age_Instance2,data1$ACR.L,N_set,tau=0.99)


res_list <- list()
for (i in 1:38) {
  idd = which(data1$Age_Instance2 == (i+44))
  idx<-which(temp_beta$Age==i+44)
  id_data<-which(data1$Age_Instance2==i+44)
  
  res_list[[i]] <- data1$ACR.L[idd]-temp_beta$beta_hat_0[idx]
}
res <- unname(unlist(res_list))

Hhat=getHnhat(abs(res))


thres_99<- Simulation1(150,Hhat,tau=0.99,bn=5.5)
write.csv(thres_99, "thres99.csv", row.names = FALSE)

#### 5
N_set<-seq(45,82,by=0.025)
N<-length(N_set)
temp_beta<-Get_beta(bn=5.5,data1$Age_Instance2,data1$ACR.L,N_set,tau=0.05)


res_list <- list()
for (i in 1:38) {
  idd = which(data1$Age_Instance2 == (i+44))
  idx<-which(temp_beta$Age==i+44)
  id_data<-which(data1$Age_Instance2==i+44)
  
  res_list[[i]] <- data1$ACR.L[idd]-temp_beta$beta_hat_0[idx]
}
res <- unname(unlist(res_list))

Hhat=getHnhat(abs(res))


thres_5<- Simulation1(150,Hhat,tau=0.05,bn=5.5)
write.csv(thres_05, "thres05.csv", row.names = FALSE)

#### 10
N_set<-seq(45,82,by=0.025)
N<-length(N_set)
temp_beta<-Get_beta(bn=5.5,data1$Age_Instance2,data1$ACR.L,N_set,tau=0.10)


res_list <- list()
for (i in 1:38) {
  idd = which(data1$Age_Instance2 == (i+44))
  idx<-which(temp_beta$Age==i+44)
  id_data<-which(data1$Age_Instance2==i+44)
  
  res_list[[i]] <- data1$ACR.L[idd]-temp_beta$beta_hat_0[idx]
}
res <- unname(unlist(res_list))

Hhat=getHnhat(abs(res))


thres_10<- Simulation1(150,Hhat,tau=0.10,bn=5.5)
write.csv(thres_10, "thres10.csv", row.names = FALSE)



tau=0.9

gap_t<-abs(thres_t/ (temp_beta$u_t*sqrt(temp_beta$local_n* bn * (tau))))
gap<-abs(thres/ (temp_beta$u_t*sqrt(temp_beta$local_n* bn * (tau))))

plot(data1$Age_Instance2,data1$ACR.L)
tt<-rq(ACR.L~Age_Instance2,data=data1,tau = 0.9)
summary(tt)


abline(tt,col="red",lwd=1)

polygon(c(temp_beta$Age, rev(temp_beta$Age)), c(2*temp_beta$beta_hat_0-temp_beta$beta_check_0-gap,
                                                rev(2*temp_beta$beta_hat_0-temp_beta$beta_check_0+gap)), col = "gray")

unique(ci_high)

polygon(c(45:82, rev(45:82)), c(unique(ci_high),
                                                rev(unique(ci_low))), col =  rgb(0, 0, 1, 0.2))

lines(N_set,2*temp_beta$beta_hat_0-temp_beta$beta_check_0,col='green',lwd=1)

# View the results
print(results_df)

plot(data1$Age_Instance2,data1$ACR.L)
tt<-rq(ACR.L~Age_Instance2,data=data1,tau = 0.9)
summary(tt)

abline(tt,col="red",lwd=1)

polygon(c(results_df$Age, rev(results_df$Age)), c(results_df$upper, rev(results_df$lower)), col = "gray")

lines(N_set,results_df$true_pred,col='green',lwd=1)

temp_beta_95<-Get_beta(bn=5.5,data1$Age_Instance2,data1$ACR.L,N_set,tau=0.95)


rfrechet(1000,shape = 5)
hist(rfrechet(1000,shape = 1/Hhat,scale = 1))

## Simulation

Simula<-function(N=1000,lower,upper, repe=50,tau=0.9){
  gap=upper-lower+1
  result<-matrix(0,nrow=N,ncol=gap)
  for (i in 1:N ) {
    set.seed(i)
    y=rep(sin(2*(1:gap)*pi/gap),repe)
    x=rep((lower:upper),repe)
    y1=rfrechet(gap*repe,shape = 10)
    #y1=rnorm(gap*repe)
    y2=y+y1
    temp_beta<-Get_beta(bn=5.5,x,y2,seq(45,82))
    result[i,]<-2*temp_beta$beta_hat_0-temp_beta$beta_check_0
  }
  true_q <- apply(result, 2, function(col) {
    quantile(col, probs = c(0.025, 0.975))  
  })
  return(true_q)
}



sim<-Simula(lower=45,upper=82)

set.seed(1245)
y=rep(sin(2*pi*(1:38)/38),50)
x=rep(c(45:82),50)
plot(x,y)
y1=rfrechet(1900,shape = 10)
#y1=rnorm(1900)
y2=y+y1
plot(x,y2)
temp_beta<-Get_beta(bn=5.5,x,y2,seq(45,82))
res_list <- list()
##Get_local_residuals(5.5, 50,temp_beta$beta_hat_0[201],temp_beta$bn_beta_hat_1[201], data1$Age_Instance2, data1$ACR.L,tau=0.5)$local_res
for (i in 1:38) {
  idd = which(x == (i+44))
  res_list[[i]] <- y[idd]-temp_beta$beta_hat_0[i]
}
res <- unname(unlist(res_list))


Hhat=getHnhat(abs(res))



#### Simulation for T_n ####
Simulation1<-function(num,Hhat){
  tau=0.9 
  bn=5.5
  result<-matrix(,nrow = 38,ncol = num)
  true_beta<-qfrechet(tau, shape = 1/Hhat)
  for (i in 1:num) {
    set.seed(i)
    print(i)
    Sim<-rfrechet(100, shape = 1/Hhat)
    temp_beta<-Get_beta(bn,x,Sim,seq(45,82))
    beta_dif<-(2*temp_beta$beta_hat_0-temp_beta$beta_check_0-true_beta)
    result[,i]<-(temp_beta$u_t*sqrt(temp_beta$local_n * bn * (tau)))*beta_dif
  }
  thres<-apply(result,1, quantile,probs=0.95)
  return(thres)  
}

tau=0.9
thres<- Simulation1(150,Hhat)
gap<-abs(thres/ (temp_beta$u_t*sqrt(temp_beta$local_n* bn * (tau))))

plot(x,y2)

polygon(c(45:82, rev(45:82)), c(sim[1,], rev(sim[2,])), col = "gray")
                                  rev(2*temp_beta$beta_hat_0-temp_beta$beta_check_0-gap)), col = rgb(0, 0, 1, 0.2))

lines(N_set,results_df$true_pred,col='green',lwd=1)
