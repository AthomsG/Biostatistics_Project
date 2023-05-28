# Model 1
log_normal_model <- survreg(Surv(time=data$TIME, event = data$CEN, type='right')
                            ~ factor(DRUG) + AGE + factor(HEP) + factor(ASI) + SERUM, 
                            data = data, 
                            dist = "lognormal")

step(log_normal_model)

summary(log_normal_model)

# Model 2
log_normal_model <- survreg(Surv(time=data$TIME, event = data$CEN, type='right')
                            ~ AGE + factor(HEP) + factor(ASI) + SERUM, 
                            data = data, 
                            dist = "lognormal")

summary(log_normal_model)




# PATIENTS WHO DIED

# Model 1
log_normal_model <- survreg(Surv(time=complete_data$TIME, event = complete_data$CEN)
                            ~ factor(DRUG) + AGE + factor(HEP) + factor(ASI) + SERUM, 
                            data = complete_data, 
                            dist = "lognormal")
summary(log_normal_model)
AIC(log_normal_model)
BIC(log_normal_model)

# Model 2
log_normal_model <- survreg(Surv(time=complete_data$TIME, event = complete_data$CEN)
                            ~ factor(DRUG) + AGE + factor(ASI) + SERUM, 
                            data = complete_data, 
                            dist = "lognormal")
summary(log_normal_model)

# Model 3
log_normal_model <- survreg(Surv(time=complete_data$TIME, event = complete_data$CEN)
                            ~ AGE + factor(ASI) + SERUM, 
                            data = complete_data, 
                            dist = "lognormal")
summary(log_normal_model)




# Calculate the residuals - residuals_c_data.pdf - 6x10

residuals_log_normal <- residuals(log_normal_model, type = "deviance")
{
  layout(matrix(c(1, 2), nrow = 2, byrow = TRUE), heights = c(2, 1))
  
  # Create a Q-Q plot of the residuals
  qqnorm(residuals_log_normal, main = "", col = "#cc5500", cex.main = 1.2) #cc5500 #2c7bb6
  abline(a = 0, b = 1, col = "red", lty = 'dashed')
  title("Q-Q Plot of Residuals", cex.main = 1.2)
  grid()
  
  # Plot the residuals
  plot(residuals_log_normal, type = "p", ylab = "Deviance Residuals", col = "#cc5500", cex.lab = 1)
  abline(h = 0, col = "red", lty = 2)
  grid()
  #hist(residuals_log_normal)
}

{
  par(mfrow = c(1, 1))
  # Create a Q-Q plot of the residuals
  qqnorm(residuals_log_normal, main="")
  abline(a = 0, b = 1, col = "red", lty='dashed')
  title("Q-Q Plot of Residuals")
  grid()
}

# Perform Shapiro-Wilk test -> Do the same before, to check the assumptions
print(shapiro.test(residuals_log_normal))
