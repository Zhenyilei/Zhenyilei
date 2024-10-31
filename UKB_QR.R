
# Set working directory
rm(list = ls())
## Windows path
setwd("C:/Users/leizh/OneDrive/Desktop/work/UKB_HCPA_FA")

library(quantreg)
library(extraDistr)
##A function that computes the high quantile regression estimate and its asymptotic covariance matrix estimate
HQRbetahat = function(Y,X,tau)
{
  n = length(Y)
  p = ncol(X)
  m = rq(Y ~ X-1,tau = tau)
  betahat = as.numeric(m$coefficients)
  Uhat = m$residuals
  ##Estimating the tail density
  dUhat = density(Uhat)
  fn0 = approx(dUhat$x,dUhat$y,xout = 0)$y
  taun = sqrt(n/tau/(1-tau))*fn0
  ##Asymptotic covariance matrix using the central limit theorem
  Sigman = t(X)%*%X/n
  aCovbetahatiid = solve(Sigman)/taun^2
  list(betahat = betahat,aCovbetahatiid = aCovbetahatiid)
}

##A function that tests the homogeneity of high quantile regression estimates
HQHTest = function(Y,X,tau1,tau2,A,cA)
{
  n = nrow(X)
  p = ncol(X)
  m1 = rq(Y ~ X-1,tau = tau1)
  Uhat1 = m1$residuals
  m2 = rq(Y ~ X-1,tau = tau2)
  Uhat2 = m2$residuals
  ##Estimating the cross quantile correlation
  IndUhat1 = as.numeric(Uhat1 > 0)
  IndUhat2 = as.numeric(Uhat2 > 0)
  lrvrho12 = as.numeric(ccf(Uhat1,Uhat2,type = "correlation",lag.max = 0,plot = FALSE)$acf)
  ##Estimating the tail density
  dUhat1 = density(Uhat1)
  fn01 = approx(dUhat1$x,dUhat1$y,xout = 0)$y
  dUhat2 = density(Uhat2)
  fn02 = approx(dUhat2$x,dUhat2$y,xout = 0)$y
  ##Compute the tau_n
  ##Make adjustment to accommodate for fixed quantiles
  taun1 = sqrt(n/tau1/(1-tau1))*fn01
  taun2 = sqrt(n/tau2/(1-tau2))*fn02
  ##Compute the test statistic
  Sigman = t(X)%*%X/n
  Phiinv = solve(A%*%solve(Sigman)%*%t(A))/(1/taun1^2 + 1/taun2^2 - 2*lrvrho12/(taun1*taun2))
  thetanA0 = A%*%(as.numeric(m1$coefficients)-as.numeric(m2$coefficients))-cA
  ##Computes the tail index estimate
  Uhat1ordered = Uhat1[order(Uhat1,decreasing = TRUE)] - min(Uhat1) + 1
  k1 = as.integer(n*(1-tau1)-1)
  Hhat1 = sum(log(Uhat1ordered[1:k1]/Uhat1ordered[k1]))/k1
  Uhat2ordered = Uhat2[order(Uhat2,decreasing = TRUE)] - min(Uhat2) + 1
  k2 = as.integer(n*(1-tau2)-1)
  Hhat2 = sum(log(Uhat2ordered[1:k2]/Uhat2ordered[k2]))/k2
  ##Outputs
  list(thetanA0 = thetanA0,Phiinv = Phiinv,PTestStatistic = t(thetanA0)%*%Phiinv%*%thetanA0,Hhat1 = Hhat1,Hhat2 = Hhat2)
}

HQHTPolynomialSA = function(Y,X,tau1,tau2,Hhat,A,cA,nsim)
{
  n = length(Y)
  ##H = 1/\lambda for Frechet distribution
  teststatcirc = rep(0,nsim)
  for(k in 1:nsim)
  {
    Ycirc = rfrechet(n,lambda = 1/Hhat,mu = 0,sigma = 1)
    teststatcirc[k] = HQHTest(Ycirc,X,tau1,tau2,A,cA)$PTestStatistic
  }
  teststatcirc
}



##Data Analysis
Data = read.csv("data.csv",header = TRUE)
X = cbind(1,Data$Age)
Y = Data$EXL

##High quantile regression estimate and the associated standard error
tau = 0.9
m = HQRbetahat(Y,X,tau)
m$betahat
sqrt(diag(m$aCovbetahatiid))
##The associated 95% confidence interval
m$betahat - qnorm(0.975)*sqrt(diag(m$aCovbetahatiid))
m$betahat + qnorm(0.975)*sqrt(diag(m$aCovbetahatiid))

tau = 0.5
m = HQRbetahat(Y,X,tau)
m$betahat
sqrt(diag(m$aCovbetahatiid))
##The associated 95% confidence interval
m$betahat - qnorm(0.975)*sqrt(diag(m$aCovbetahatiid))
m$betahat + qnorm(0.975)*sqrt(diag(m$aCovbetahatiid))

tau = 0.1
m = HQRbetahat(Y,X,tau)
m$betahat
sqrt(diag(m$aCovbetahatiid))
##The associated 95% confidence interval
m$betahat - qnorm(0.975)*sqrt(diag(m$aCovbetahatiid))
m$betahat + qnorm(0.975)*sqrt(diag(m$aCovbetahatiid))



##Testing for slope homogeneity
tau1 = 0.9
tau2 = 0.5
A = matrix(c(0,1),nrow = 1,ncol = 2)
cA = 0
testoutput = HQHTest(Y,X,tau1,tau2,A,cA)
pchisq(as.numeric(testoutput$PTestStatistic),df = nrow(A),lower.tail = FALSE)

Hhat = mean(c(testoutput$Hhat1,testoutput$Hhat2))
nsim = 100
PTcirc = HQHTPolynomialSA(Y,X,tau1,tau2,Hhat,A,cA,nsim)
##namePTcirc = paste("Results/SAYX","kappa",round(1/Hhat),"tau1",tau1*100,"tau2",tau2*100,".txt",sep = "")
##write.table(as.matrix(PTcirc),namePTcirc,quote = FALSE,row.names = FALSE,col.names = FALSE,append = TRUE)
##PTcyclonecirc = read.table(nameBRHQcyclone,header = FALSE)$V1
sum(PTcirc > as.numeric(testoutput$PTestStatistic))/length(PTcirc)




data<-read.csv("UKB_FA_dat.csv")


var_list <- c("EX.L","EX.R", "SFO.R","SFO.L", "PTR.R", "PTR.L","TAP.R", "TAP.L",
                   "ALIC.R","ALIC.L", "PLIC.R","PLIC.L", "SLF.R", "SCC")

for (var in var_list) {
  print(head(data[[var]]))
}

Y <- data[["EX.L"]] 

X = cbind(1,data$Age_Instance2)

tau_values <- seq(0.1, 0.9, by = 0.1)


compute_trend_results <- function(var, tau) {
  Y <- data[[var]]  # Set Y to the current variable
  
  m <- HQRbetahat(Y, X, tau)
  
  # Calculate betahat and CI
  betahat <- m$betahat[2]
  se <- sqrt(diag(m$aCovbetahatiid))[2]
  ci_lower <- betahat - qnorm(0.975) * se
  ci_upper <- betahat + qnorm(0.975) * se
  
  # Return results as a list
  return(list(betahat = betahat, ci_lower = ci_lower, ci_upper = ci_upper))
}

# Use lapply to iterate over each variable
results_list <- lapply(var_list, function(var) {
  # Use sapply to iterate over each tau for the given variable
  sapply(tau_values, function(tau) {
    compute_results(var, tau)
  }, simplify = FALSE)
})

# Name the results list by variable names
names(results_list) <- var_list
results_df <- data.frame()
var <- "EX.L"
results_list[[var]]
# Loop through each variable and tau value in results_list to populate results_df
for (var in names(results_list)) {
  for (i in seq_along(results_list[[var]])) {
    res <- results_list[[var]][[i]]
    
    # Create a data frame row for the current variable and tau
    row <- data.frame(
      variable = var,
      tau = i / 10,                     # Assuming tau = 0.1, 0.2, ..., 0.9
      betahat = res$betahat[2],         # Extracting the second parameter
      ci_lower = res$ci_lower[2],       # Second parameter lower CI
      ci_upper = res$ci_upper[2]        # Second parameter upper CI
    )
    
    # Append the row to results_df
    results_df <- rbind(results_df, row)
  }
}

write.csv(results_df,"trend_and_ci.csv",row.names = FALSE)
library(ggplot2)

# Plot betahat and confidence intervals by tau, faceted by variable
ggplot(results_df, aes(x = tau, y = betahat, group = variable)) +
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), fill = "blue", alpha = 0.2) + # Shaded CI region
  geom_line(color = "black") +                                # Line connecting betahat points
  geom_point(color = "black", size = 2) +                     # Solid points for betahat
  facet_wrap(~ variable, scales = "free_y",axes = "all_x") +                # Separate plot for each variable
  labs(x = "Tau", y = "Trend and 95% CI",                  # Labels for axes
       title = "Trend and Confidence Intervals by Tau for Each Variable") +
  scale_x_continuous(breaks = tau_values, limits = c(0.1, 0.9)) + # Consistent x-axis across facets
  theme_minimal()+                          # Add more space between facets
  theme(axis.text.x = element_text(angle = 45, hjust = 1),     # Rotate x-axis labels for better visibility
        strip.text = element_text(size = 10),                  # Adjust facet label text size
        plot.title = element_text(size = 14, face = "bold"))   # Customize title appearance


boxplot(EX.L ~ Age_Instance2, data = data, outline = TRUE, # outline = FALSE hides outliers in boxplot
        xlab = "Age", ylab = "EX.L", main = "Boxplot with Quantile Regression Lines",
        col = "lightgrey", border = "black")
quan_test <- c(0.1, 0.5, 0.9)
colors <- c("blue", "red", "green")  

for (i in 1:3) {
  tau <- quan_test[i]
  color <- colors[i]
  
  
  model <- rq(EX.L ~ Age_Instance2, tau = tau, data = data)
  
  # Add the quantile regression lines to the plot by matching factor levels on the x-axis
  age_levels <- 1:length(unique(data$Age_Instance2))
  fitted_values <- coef(model)[1] + coef(model)[2] * age_levels
  points(x = age_levels,
         y = fitted_values, type = "l", col = color, lwd = 2)
}

# Add a legend to label each quantile line
legend("topright", legend = c("10th Percentile", "Median (50th)", "90th Percentile"),
       col = colors, lwd = 2, cex = 0.8)
