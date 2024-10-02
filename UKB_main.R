# Set working directory
rm(list = ls())
#setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(quantreg)
library(VGAM)
data1 <- read.csv("~/Desktop/Ting zhang Project/UKB_HCPA_FA/UKB_FA_dat.csv")

sort(unique(data2$age))

hist(data1$Age_Instance2)
bn<-5.5
x0<-50
plot(data1$Age_Instance2,data1$ACR.L)
tt<-rq(ACR.L~Age_Instance2,data=data1,tau = 0.9)
mt<-summary(tt)
mt$coefficients[2,1]

test_data <- data1[(-bn <= data1$Age - x0 & data1$Age - x0 <= bn), ]
test_data$Age=(test_data$Age_Instance2 - x0) / bn
tt<-rq(ACR.L~Age,data=test_data,tau = 0.9)
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
  
  return(list(beta_0 = beta_bn_beta_hat[1],
              bn_beta_1 = beta_bn_beta_hat[2]))
}


##### Get the parameter for given bn #####
Get_beta<-function(bn,x,y,N_set,tau=0.9){
  N<- length(N_set)
  beta_df <- data.frame(
    Age = N_set,beta_hat_0 = numeric(N),bn_beta_hat_1 = numeric(N),
    beta_check_0 = numeric(N),bn_beta_check_1 = numeric(N)
  )

  for (i in seq_along(N_set)) {
    beta_hat<- LLQR(bn, N_set[i] , x, y,tau=tau)
    beta_check<-LLQR(sqrt(2)*bn,N_set[i],x,y,tau=tau)
    
    beta_df$beta_hat_0[i] <- beta_hat$beta_0
    beta_df$bn_beta_hat_1[i] <- beta_hat$bn_beta_1
    beta_df$beta_check_0[i] <- beta_check$beta_0
    beta_df$bn_beta_check_1[i] <- beta_check$bn_beta_1
  }
  return(beta_df)
}
N_set<-seq(45,82,by=0.025)
N<-length(N_set)
temp_beta<-Get_beta(bn=5.5,data1$Age_Instance2,data1$ACR.L,N_set)

##### Get  local residuals and density of residuals #####
Get_local_residuals <- function(bn,x0,beta_0,beta_1,x,y,tau=0.9) {
  ## get the index of residual at age x0
  idx0<- x==x0
  ## get index of residual within bandwidth
  id<-(abs(x - x0) <= bn)
  ## draw the data within the bandwidth
  test_x<-x[id]
  test_y<-y[id]
  nrow<-length(test_x)
  ## standradlize the data
  test_stand_x<-(test_x - x0) / bn
  ## calculate residuals 
  res<-test_y-beta_0-beta_1*test_stand_x
  ## get the density of residuals
  dens <- density(res)
  ## get the index of 90% quantile of residuals 
  idx <- which.min(abs(dens$x - quantile(res,tau)))
  ## get the density of value of 90% quantile
  u_t<-dens$y[idx]
  return(list(u_t=u_t,nrow=nrow))
}

for (i in 1:N) {
  temp_beta$u_t[i]<-Get_local_residuals(5.5,N_set[i],temp_beta$beta_hat_0[i],temp_beta$bn_beta_hat_1[i],
                                        data1$Age_Instance2, data1$ACR.L,tau=0.9)$u_t
  temp_beta$nrow[i]<-Get_local_residuals(5.5,N_set[i],temp_beta$beta_hat_0[i],temp_beta$bn_beta_hat_1[i],
                                        data1$Age_Instance2, data1$ACR.L,tau=0.9)$nrow
}


res_list <- list()
##Get_local_residuals(5.5, 50,temp_beta$beta_hat_0[201],temp_beta$bn_beta_hat_1[201], data1$Age_Instance2, data1$ACR.L,tau=0.5)$local_res
for (i in 1:38) {
  idx<-which(temp_beta$Age==i+44)
  id_data<-which(data1$Age_Instance2==i+44)
  res_list[[i]] <- data1$ACR.L-temp_beta$beta_hat_0[idx]
}
res <- unname(unlist(res_list))


##### calculate the Hill estimator #####
getHnhat <- function(X, alphan = 0.1) {
  n <- length(X)
  
  kn <- floor(n * alphan)

  X_sorted <- sort(X, decreasing = TRUE)

  H_hat <- sum(log(X_sorted[1:kn]) - log(X_sorted[kn + 1])) / kn
  
  return(H_hat)
}

Hhat=getHnhat(abs(res))



#### Simulation for T_n ####
Simulation1<-function(num,Hhat){
  tau=0.9
  bn=5.5
  result<-matrix(,nrow = N,ncol = num)
  true_beta<-qfrechet(tau, shape = 1/Hhat)
  for (i in 1:num) {
    set.seed(i)
    print(i)
    Sim<-rfrechet(40414, shape = 1/Hhat)
    temp_beta<-Get_beta(bn,data1$Age_Instance2,Sim,N_set)
    for (j in 1:N) {
      para<-Get_local_residuals(bn,N_set[j],temp_beta$beta_hat_0[j],temp_beta$bn_beta_hat_1[j],
                                data1$Age_Instance2, data1$ACR.L,tau=0.9)
      temp_beta$u_t[j]<-para$u_t
      temp_beta$nrow[j]<-para$nrow
    }
    beta_dif<-(2*temp_beta$beta_hat_0-temp_beta$beta_check_0-true_beta)
    result[,i]<-(temp_beta$u_t*sqrt(temp_beta$nrow * bn * (1-tau)))*beta_dif
  }
  thres<-apply(result,1, quantile,probs=0.95)
  return(thres)  
}

thres<- Simulation1(150,Hhat)


 # View the results
print(results_df)

plot(data1$Age_Instance2,data1$ACR.L)
tt<-rq(ACR.L~Age_Instance2,data=data1,tau = 0.9)
summary(tt)

abline(tt,col="red",lwd=1)

polygon(c(results_df$Age, rev(results_df$Age)), c(results_df$upper, rev(results_df$lower)), col = "gray")

lines(N_set,results_df$true_pred,col='green',lwd=1)


gt<-summary(test_md)
gt$residuals
quantile(gt$residuals,0.0)
quantile(test_md$residuals,0.9)*sqrt(3185)



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
    result[i, ] <- sapply(lower:upper, function(j) {
      test_function_2(5,j,x,y2)$true_pred
    })
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

res_list <- list()

for (i in 1:38) {
  res_list[[i]] <- Get_residuals(5, i + 44, x,y2)
}
res <- unname(unlist(res_list))

Hhat=getHnhat(abs(res))
thres <- quantile(rfrechet(100000, shape = 1/Hhat),0.95)
results_df <- data.frame(Age = integer(), true_pred = numeric(), lower_p = numeric(), upper_p = numeric())
N_set<-seq(1:1000)*37/1000+45
for (i in N_set) {
  result <- test_function_2(5, i , x, y2)
  results_df <- rbind(results_df, data.frame(Age =  i , true_pred = result$true_pred, lower_p = result$lower_p, upper_p = result$upper_p))
}

plot(x,y2)

polygon(c(45:82, rev(45:82)), c(sim[1,], rev(sim[2,])), col = "gray")

polygon(c(results_df$Age, rev(results_df$Age)), c(results_df$upper, rev(results_df$lower)), col = "red")

lines(N_set,results_df$true_pred,col='green',lwd=1)


