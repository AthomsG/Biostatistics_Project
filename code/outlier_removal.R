# REMOVE OUTLIERS AND SEE HOW IT GOES

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
data$DRUG <- ifelse(data$DRUG == 1, 'D-Penicillamine', ifelse(data$DRUG == 2, 'Placebo', data$DRUG))

# Subset the data for ASI == "YES"
subset_data <- data[data$ASI == 'Yes', ]

# Extract the TIME values for ASI == "YES"
time_values <- subset_data$TIME

# Calculate the lower and upper quartiles
q1 <- quantile(time_values, 0.25)
q3 <- quantile(time_values, 0.75)

# Calculate the interquartile range (IQR)
iqr <- q3 - q1

# Define the lower and upper boundaries for outliers
lower_bound <- q1 - 1.5 * iqr
upper_bound <- q3 + 1.5 * iqr

# Detect outliers
outliers <- subset_data$TIME[subset_data$TIME < lower_bound | subset_data$TIME > upper_bound]

# Remove outliers from data
data <- subset(data, !(data$ASI == 'Yes' & (data$TIME < lower_bound | data$TIME > upper_bound)))


# TROUBLE SHOOTING
par(mfrow = c(3, 1))
hist(log(data$TIME), main='log(Survival Time)', xlab='log(time)')

# Subset the data based on ASI value
data_ASI0 <- data[data$ASI == 'No', ]
data_ASI1 <- data[data$ASI == 'Yes', ]

# Histogram for ASI == 'No'
hist(log(data_ASI0$TIME), main = 'log(Survival Time) - ASI = No', xlab = 'log(time)')

# Histogram for ASI == 'Yes'
hist(log(data_ASI1$TIME), main = 'log(Survival Time) - ASI = Yes', xlab = 'log(time)')

# Reset the plot layout to default
par(mfrow = c(1, 1))

# Create a new dataframe with log-transformed survival time and ASI variable
df <- data.frame(log_time = log(data$TIME), ASI = data$ASI)

# Plot the histograms with transparency
ggplot(df, aes(x = log_time, fill = ASI)) +
  geom_histogram(alpha = 0.5, position = "identity", bins=15) +
  labs(x = "log(time)", y = "Frequency", title = "Histogram of log(Survival Time)") +
  scale_fill_manual(values = c("No" = "blue", "Yes" = "red"), labels = c("No", "Yes")) +
  theme_minimal()

# Create a new column named "Group" in the data dataframe
data$Group <- ifelse(data$ASI == "No" & data$HEP == "No", "No", "Yes")

# Create a new dataframe with log-transformed survival time and ASI variable
df <- data.frame(log_time = log(data$TIME), Group = data$Group)

# Plot the histograms with transparency
ggplot(df, aes(x = log_time, fill = Group)) +
  geom_histogram(alpha = 0.5, position = "identity", bins=15) +
  labs(x = "log(time)", y = "Frequency", title = "Histogram of log(Survival Time)") +
  scale_fill_manual(values = c("No" = "blue", "Yes" = "red"), labels = c("No", "Yes")) +
  theme_minimal()
