# The log-normal regression model does not seem appropriate to the model as is. The assumption of log(Y) normality does not hold.
# in fact, this variable follows a bell shaped distribution with a longer left tail.
# Let's test if we can divide the dataset into two datasets based on ASI

#####################
# DATA ORGANIZATION #
#####################

data=read.table(file="DatasetM.txt", header=TRUE)
# Non-informative variable
data <- subset(data, select = -ID)

# Convert Age in days to years (Easier to interpret)
data$AGE <- data$AGE / 365
data$TIME <- data$TIME / 365

# Transform numerals to Strings (Easier to interpret)
data$ASI  <- ifelse(data$ASI == 0, 'No', ifelse(data$ASI == 1, 'Yes', data$ASI))
data$HEP  <- ifelse(data$HEP == 0, 'No', ifelse(data$HEP == 1, 'Yes', data$HEP))
#data$CEN  <- ifelse(data$CEN == 0, 'Censored', ifelse(data$CEN == 1, 'Death', data$CEN))
data$DRUG <- ifelse(data$DRUG == 1, 'D-Penicillamine', ifelse(data$DRUG == 2, 'Placebo', data$DRUG))

yASI <- data[data$ASI=="Yes",]
nASI <- data[data$ASI=="No",]

{
  par(mfrow=c(1, 2))
  hist(log(yASI$TIME), main='ASI = Yes', xlab='log(TIME)')
  hist(log(nASI$TIME), main='ASI = No' , xlab='log(TIME)')
  par(mfrow=c(1, 1))
}

################################
# LOG-NORMAL REGRESSION MODELS #
################################

# Model 1 -> ASI
log_normal_model <- survreg(Surv(time=yASI$TIME, event = yASI$CEN)
                            ~ factor(DRUG) + AGE + factor(HEP) + SERUM, 
                            data = yASI, 
                            dist = "lognormal")

AIC(log_normal_model)
BIC(log_normal_model)
summary(log_normal_model)

# Model 2 -> ASI
log_normal_model <- survreg(Surv(time=yASI$TIME, event = yASI$CEN)
                            ~ SERUM, 
                            data = yASI, 
                            dist = "lognormal")

AIC(log_normal_model)
BIC(log_normal_model)
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

# ---> Residuals from the fit to the ASI population, according to the Shapiro-Wilk test follow a normal distribution

# Model 1 -> NO ASI
log_normal_model <- survreg(Surv(time=nASI$TIME, event = nASI$CEN, type='right')
                            ~ AGE + factor(HEP) + SERUM, 
                            data = nASI, 
                            dist = "lognormal")

AIC(log_normal_model)
BIC(log_normal_model)
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



# WILL NOW SEE IF I CAN SEPARATE THE POPULATION AGAIN, BASED ON HEP

yHEP <- nASI[nASI$HEP=="Yes",]
nHEP <- nASI[nASI$HEP=="No",]

{
  par(mfrow=c(1, 2))
  hist(log(yHEP$TIME))
  hist(log(nHEP$TIME))
  par(mfrow=c(1, 1))
}

# Model 1 -> HEP
log_normal_model <- survreg(Surv(time=yHEP$TIME, event = yHEP$CEN, type = 'right')
                            ~ AGE + SERUM, 
                            data = yHEP, 
                            dist = "lognormal")

AIC(log_normal_model)
BIC(log_normal_model)
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

# Model 1 -> nHEP
log_normal_model <- survreg(Surv(time=nHEP$TIME, event = nHEP$CEN, type='right')
                            ~ AGE + SERUM, 
                            data = nHEP, 
                            dist = "lognormal")

AIC(log_normal_model)
BIC(log_normal_model)
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
