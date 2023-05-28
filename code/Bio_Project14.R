library(ggplot2)
library(survival)
library(ggfortify)
library(gridExtra)

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

############################
# EXPORATORY DATA ANALYSIS #
############################

# Initial Exploratory Data Analysis on Patients who died (Disregard Censored data)

# Filter out censored data
complete_data <- data[data$CEN == 1, ] # Patients that have died

# Exploratory Data Analysis on Patients all data (Include Censored data)

# Create boxplot for Presence of Asictes
plot1 <- ggplot(complete_data, aes(x = ASI, y = TIME)) +
  geom_boxplot(fill = "#2c7bb6", color = "black", alpha = 0.8) +
  labs(x = "Presence of Asictes", y = "Survival Time (years)") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text = element_text(size = 13),
    legend.title = element_blank(),
    legend.text = element_text(size = 12)
  )

# Create boxplot for Presence of hepatomegaly
plot2 <- ggplot(complete_data, aes(x = HEP, y = TIME)) +
  geom_boxplot(fill = "#2c7bb6", color = "black", alpha = 0.8) +
  labs(x = "Presence of hepatomegaly") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_blank(),
    axis.text = element_text(size = 13),
    legend.title = element_blank(),
    legend.text = element_text(size = 12)
  )

# Create boxplot for Drug
plot3 <- ggplot(complete_data, aes(x = DRUG, y = TIME)) +
  geom_boxplot(fill = "#2c7bb6", color = "black", alpha = 0.8) +
  labs(x = "Drug") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_blank(),
    axis.text = element_text(size = 13),
    legend.title = element_blank(),
    legend.text = element_text(size = 12)
  )

print(plot1)

# Arrange plots side by side
grid.arrange(plot1, plot2, plot3, ncol = 3)

# 8 x 6
{
  # Set the layout to a 3x2 grid
  par(mfrow = c(2, 3))
  
  # Bar plot for Drug
  barplot(table(complete_data$DRUG),
          main = "Drug",
          ylab = "Frequency",
          col = "#2c7bb6",
          border = "white",
          cex.main = 1.8,
          cex.names= 1.3,
          cex.axis = 1.6,
          cex.lab = 1.5)
  
  # Bar plot for Presence of Asictes
  barplot(table(complete_data$ASI),
          main = "Presence of\n Asictes",
          ylab = "Frequency",
          col = "#2c7bb6",
          border = "white",
          cex.main = 1.8,
          cex.names= 1.6,
          cex.axis = 1.6,
          cex.lab = 1.5)
  
  # Bar plot for Presence of Hepatomegaly
  barplot(table(complete_data$HEP),
          main = "Presence of\n Hepatomegaly",
          ylab = "Frequency",
          col = "#2c7bb6",
          border = "white",
          cex.main = 1.8,
          cex.names= 1.6,
          cex.axis = 1.6,
          cex.lab = 1.5)
  
  # Histogram for Serum
  hist(complete_data$SERUM,
       breaks = 15,
       main = "Serum",
       xlab = "Serum cholesterol in mg/dl",
       col = "#2c7bb6",
       border = "white",
       cex.main = 1.8,
       cex.axis = 1.6,
       cex.lab = 1.5)
  
  # Histogram for Age
  hist(complete_data$AGE,
       breaks = 15,
       main = "Age",
       xlab = "Age in years",
       col = "#2c7bb6",
       border = "white",
       cex.main = 1.8,
       cex.axis = 1.6,
       cex.lab = 1.5)
  
  # Histogram for Survival Time
  hist(complete_data$TIME,
       breaks = 15,
       main = "Survival Time",
       xlab = "Survival Time in years",
       col = "#2c7bb6",
       border = "white",
       cex.main = 1.8,
       cex.axis = 1.6,
       cex.lab = 1.5)
  
  # Reset the layout to a 1x1 grid
  par(mfrow = c(1, 1))
}

# Survival Probability Plot

# Calculate survival probabilities using Kaplan-Meier estimator
surv_object <- Surv(time = complete_data$TIME, event = complete_data$CEN)
km_fit <- survfit(surv_object ~ DRUG, complete_data)
# Plot Kaplan-Meier survival curves for each drug
autoplot(km_fit, xlab = "Time (in years)", ylab = "Survival Probability", main="Survival Function of patients that died") +
  theme_minimal(base_size=16) +
  theme(panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "gray", linetype = "dashed"),
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.position = "top") +
  xlim(0, 12.2) + 
  ylim(0, 1)

km_fit <- survfit(surv_object ~ HEP, complete_data)
# Plot Kaplan-Meier survival curves for each drug
autoplot(km_fit, xlab = "Time (in years)", ylab = "Survival Probability", main="Survival Function of patients that died\nregarding having hepatomegaly") +
  theme_minimal(base_size=16) +
  theme(panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "gray", linetype = "dashed"),
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.position = "top") +
  xlim(0, 12.2) + 
  ylim(0, 1)

# Calculate survival probabilities using Kaplan-Meier estimator
km_fit <- survfit(surv_object ~ ASI, complete_data)
# Plot Kaplan-Meier survival curves for each drug
autoplot(km_fit, xlab = "Time (in years)", ylab = "Survival Probability", main="Survival Function of patients that died\nregarding having ascites") +
  theme_minimal(base_size=16) +
  theme(panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "gray", linetype = "dashed"),
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.position = "top") +
  xlim(0, 12.2) + 
  ylim(0, 1)

# Dataset Exploratory Data Analysis (Include Censored Data)

# Create boxplot for Presence of Asictes
plot1 <- ggplot(data, aes(x = ASI, y = TIME)) +
  geom_boxplot(fill = "#cc5500", color = "black", alpha = 0.8) +
  labs(x = "Presence of Asictes", y = "Survival Time (years)") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text = element_text(size = 13),
    legend.title = element_blank(),
    legend.text = element_text(size = 12)
  )

# Create boxplot for Presence of hepatomegaly
plot2 <- ggplot(data, aes(x = HEP, y = TIME)) +
  geom_boxplot(fill = "#cc5500", color = "black", alpha = 0.8) +
  labs(x = "Presence of hepatomegaly") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_blank(),
    axis.text = element_text(size = 13),
    legend.title = element_blank(),
    legend.text = element_text(size = 12)
  )

# Create boxplot for Drug
plot3 <- ggplot(data, aes(x = DRUG, y = TIME)) +
  geom_boxplot(fill = "#cc5500", color = "black", alpha = 0.8) +
  labs(x = "Drug") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_blank(),
    axis.text = element_text(size = 13),
    legend.title = element_blank(),
    legend.text = element_text(size = 12)
  )

# Arrange plots side by side
grid.arrange(plot1, plot2, plot3, ncol = 3)

# 8 x 6
{
# Set the layout to a 3x2 grid
par(mfrow = c(2, 3))

# Bar plot for Drug
barplot(table(data$DRUG),
        main = "Drug",
        ylab = "Frequency",
        col = "#cc5500",
        border = "white",
        cex.main = 1.8,
        cex.names= 1.3,
        cex.axis = 1.6,
        cex.lab = 1.5)

# Bar plot for Presence of Asictes
barplot(table(data$ASI),
        main = "Presence of\n Asictes",
        ylab = "Frequency",
        col = "#cc5500",
        border = "white",
        cex.main = 1.8,
        cex.names= 1.6,
        cex.axis = 1.6,
        cex.lab = 1.5)

# Bar plot for Presence of Hepatomegaly
barplot(table(data$HEP),
        main = "Presence of\n Hepatomegaly",
        ylab = "Frequency",
        col = "#cc5500",
        border = "white",
        cex.main = 1.8,
        cex.names= 1.6,
        cex.axis = 1.6,
        cex.lab = 1.5)

# Histogram for Serum
hist(data$SERUM,
     breaks = 15,
     main = "Serum",
     xlab = "Serum cholesterol in mg/dl",
     col = "#cc5500",
     border = "white",
     cex.main = 1.8,
     cex.axis = 1.6,
     cex.lab = 1.5)

# Histogram for Age
hist(data$AGE,
     breaks = 15,
     main = "Age",
     xlab = "Age in years",
     col = "#cc5500",
     border = "white",
     cex.main = 1.8,
     cex.axis = 1.6,
     cex.lab = 1.5)

# Histogram for Survival Time
hist(data$TIME,
     breaks = 15,
     main = "Survival Time",
     xlab = "Survival Time in years",
     col = "#cc5500",
     border = "white",
     cex.main = 1.8,
     cex.axis = 1.6,
     cex.lab = 1.5)

# Reset the layout to a 1x1 grid
par(mfrow = c(1, 1))
}

# Calculate the Pearson correlation coefficient
correlation <- cor(data$TIME, data$SERUM, method = "pearson")

# Create a scatter plot
plot(data$TIME, data$SERUM, pch = 16, xlab = "TIME", ylab = "SERUM",
     main = paste("Pearson Correlation:", round(correlation, 2)))

# Add a trendline (optional)
abline(lm(data$SERUM ~ data$TIME), col = "red")


# Survival Probability Plot

# Calculate survival probabilities using Kaplan-Meier estimator
surv_object <- Surv(time = data$TIME, event = data$CEN)
km_fit <- survfit(surv_object ~ DRUG, data)
# Plot Kaplan-Meier survival curves for each drug
autoplot(km_fit, xlab = "Time (in years)", ylab = "Survival Probability", main="Survival Function of all patients") +
  theme_minimal(base_size=16) +
  theme(panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "gray", linetype = "dashed"),
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.position = "top") +
  xlim(0, 12.2) + 
  ylim(0, 1)


# Calculate survival probabilities using Kaplan-Meier estimator
km_fit <- survfit(surv_object ~ HEP, data)
# Plot Kaplan-Meier survival curves for each drug
autoplot(km_fit, xlab = "Time (in years)", ylab = "Survival Probability", main="Survival Function of all patients\nregarding having hepatomegaly") +
  theme_minimal(base_size=16) +
  theme(panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "gray", linetype = "dashed"),
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.position = "top") +
  xlim(0, 12.2) + 
  ylim(0, 1)



# Calculate survival probabilities using Kaplan-Meier estimator
km_fit <- survfit(surv_object ~ ASI, data)
# Plot Kaplan-Meier survival curves for each drug
autoplot(km_fit, xlab = "Time (in years)", ylab = "Survival Probability", main="Survival Function of all patients\nregarding having ascites") +
  theme_minimal(base_size=16) +
  theme(panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "gray", linetype = "dashed"),
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.position = "top") +
  xlim(0, 12.2) + 
  ylim(0, 1)


################################
# LOG-NORMAL REGRESSION MODELS # WALD TEST ON SUMMARY!
################################

# LOG-NORMAL FOR THE COMPLETE DATA ################################################################

# Model 1
log_normal_model <- survreg(Surv(time=complete_data$TIME, event = complete_data$CEN)
                            ~ factor(DRUG) + AGE + factor(HEP) + factor(ASI) + SERUM, 
                            data = complete_data, 
                            dist = "lognormal")

AIC(log_normal_model)
BIC(log_normal_model)
summary(log_normal_model)

# Model 2
log_normal_model <- survreg(Surv(time=complete_data$TIME, event = complete_data$CEN)
                            ~ factor(HEP) + factor(ASI) + SERUM, 
                            data = complete_data, 
                            dist = "lognormal")

AIC(log_normal_model)
BIC(log_normal_model)
summary(log_normal_model)

# Model 3
log_normal_model <- survreg(Surv(time=complete_data$TIME, event = complete_data$CEN)
                            ~ factor(ASI) + SERUM, 
                            data = complete_data, 
                            dist = "lognormal")

AIC(log_normal_model)
BIC(log_normal_model)
summary(log_normal_model)

# LOG-NORMAL FOR THE CENSORED DATA (including complete data) ######################################

# Model 1
log_normal_model <- survreg(Surv(time=data$TIME, event = data$CEN, type='right')
                            ~ factor(DRUG) + AGE + factor(HEP) + factor(ASI) + SERUM, 
                            data = data, 
                            dist = "lognormal")

AIC(log_normal_model)
BIC(log_normal_model)

summary(log_normal_model)

step(log_normal_model)

# Model 2
log_normal_model <- survreg(Surv(time=data$TIME, event = data$CEN, type='right')
                            ~ factor(HEP) + factor(ASI) + SERUM, 
                            data = data, 
                            dist = "lognormal")
AIC(log_normal_model)
BIC(log_normal_model)
summary(log_normal_model)

# Model 3
log_normal_model <- survreg(Surv(time=data$TIME, event = data$CEN, type='right')
                            ~ factor(ASI) + SERUM, 
                            data = data, 
                            dist = "lognormal")
AIC(log_normal_model)
BIC(log_normal_model)
summary(log_normal_model)

# Additional analysis or visualization of residuals can be performed here

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