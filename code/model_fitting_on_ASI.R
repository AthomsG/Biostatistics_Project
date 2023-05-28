################################
# LOG-NORMAL REGRESSION MODELS # WALD TEST ON SUMMARY!
################################

# COMPLETE DATA

# Model 1
log_normal_model <- survreg(Surv(time=complete_data$TIME, event = complete_data$CEN)
                            ~ factor(DRUG) + AGE + factor(HEP) + factor(ASI) + SERUM, 
                            data = complete_data, 
                            dist = "lognormal")
summary(log_normal_model)

step(log_normal_model)

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

# Model 4
log_normal_model <- survreg(Surv(time=complete_data$TIME, event = complete_data$CEN)
                            ~ factor(ASI) + SERUM, 
                            data = complete_data, 
                            dist = "lognormal")
summary(log_normal_model)

# Calculate the residuals

residuals_log_normal <- residuals(log_normal_model, type = "deviance")
{
  par(mfrow = c(2, 1))
  # Plot the residuals
  plot(residuals_log_normal, type = "p", ylab = "Deviance Residuals", main = "Deviance Residuals Plot")
  abline(h = 0, col = "red", lty = 2)  # Add a horizontal line at zero
  
  hist(residuals_log_normal)
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




























# ALL PATIENT DATA

# Model 1
log_normal_model <- survreg(Surv(time=data$TIME, event = data$CEN)
                            ~ factor(DRUG) + AGE + factor(HEP) + factor(ASI) + SERUM, 
                            data = data, 
                            dist = "lognormal")
summary(log_normal_model)

# Model 2
log_normal_model <- survreg(Surv(time=data$TIME, event = data$CEN)
                            ~ AGE + factor(HEP) + factor(ASI) + SERUM, 
                            data = data, 
                            dist = "lognormal")
summary(log_normal_model)

# Calculate the residuals

residuals_log_normal <- residuals(log_normal_model, type = "deviance")
{
  par(mfrow = c(2, 1))
  # Plot the residuals
  plot(residuals_log_normal, type = "p", ylab = "Deviance Residuals", main = "Deviance Residuals Plot")
  abline(h = 0, col = "red", lty = 2)  # Add a horizontal line at zero
  
  hist(residuals_log_normal)
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






























# ONLY ASI

yASI <- data[data$ASI=="Yes",]
nASI <- data[data$ASI=="No",]

# Model 1
log_normal_model <- survreg(Surv(time=yASI$TIME, event = yASI$CEN)
                            ~ factor(DRUG) + AGE + factor(HEP) + SERUM, 
                            data = yASI, 
                            dist = "lognormal")
summary(log_normal_model)

# Model 2
log_normal_model <- survreg(Surv(time=yASI$TIME, event = yASI$CEN)
                            ~ factor(DRUG) + AGE + SERUM, 
                            data = yASI, 
                            dist = "lognormal")
summary(log_normal_model)

# Model 3
log_normal_model <- survreg(Surv(time=yASI$TIME, event = yASI$CEN)
                            ~ AGE + SERUM, 
                            data = yASI, 
                            dist = "lognormal")
summary(log_normal_model)

# Model 4
log_normal_model <- survreg(Surv(time=yASI$TIME, event = yASI$CEN)
                            ~ SERUM, 
                            data = yASI, 
                            dist = "lognormal")
summary(log_normal_model)