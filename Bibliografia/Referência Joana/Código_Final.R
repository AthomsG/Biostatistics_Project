library(pca3d)
library(readr)
library(psych)
library(GGally)
library(ggplot2)
library(pdp)
library(RVAideMemoire)
library(dplyr)
library(PerformanceAnalytics)
library(rospca)
library(carData)
library(car)
library(tidyverse)
library(splitstackshape)
library(vctrs)
library(ggplot2)
library(lattice)
library(caret)
library(randomForest)
library(rpart)
library(rpart.plot)
library(ggcorrplot)
require(foreign)
require(ggplot2)
require(MASS)
require(Hmisc)
require(reshape2)
require(car)
require(brant)

dat<-read.csv("C:/Users/Fernando/Desktop/Bio/Maternal Health Risk Data Set.csv")

names(dat) <- c('Age', 'SystolicBP', 'DiastolicBP', 'BS', 'BodyTemp', 'HeartRate', 'RiskLevel')

#Summary and descriptive statistics
summary(dat)
describe(dat)

dat$RiskLevel<-as.factor(dat$RiskLevel)
dat$Age <- as.numeric(dat$Age)
dat$SystolicBP <- as.numeric(dat$SystolicBP)
dat$DiastolicBP <- as.numeric(dat$DiastolicBP)
dat$BS <- as.numeric(dat$BS)
dat$BodyTemp <- as.numeric(dat$BodyTemp)
dat$HeartRate <- as.numeric(dat$HeartRate)


#Missing Values
sum(is.na(dat))   #no missing values

#RPCA for Outlier Detection
resR <- robpca(
  dat[,-length(dat)],
  k=0,
  alpha = 0.95
  )

dat$Colour[dat$RiskLevel=="low risk"]="green"
dat$Colour[dat$RiskLevel=="mid risk"]="blue"
dat$Colour[dat$RiskLevel=="high risk"]="red"

plot(x=resR$sd,
     y=resR$od,
     xlim = c(0,10),
     ylim = c(0,100),
     pch=16,
     col= dat$Colour, 
     xlab = "Score Distance", 
     ylab = "Orthogonal Distance",
     )

legend(x="topright",
       legend=c("Low Risk","Mid Risk","High Risk")
      )

abline(h = resR$cutoff.od,v = resR$cutoff.sd)
dat$Colour<-NULL

################
dat$sd<-resR$sd#
dat$od<-resR$od#
dat$sd<-NULL   #    
dat$od<-NULL   #
################

dat<-dat[-c(500, 909), ]
rownames(dat) <- NULL

#PCA 3d for Visualization
pca <- prcomp(dat[,-length(dat)], scale.=TRUE)
gr <- factor(dat$RiskLevel)

round((pca$rotation),3)
summary(pca)
summary(gr)

pca3d(pca, 
      group=gr, 
      legend="topleft",
      shape = "sphere",
      palette = c("red","green","blue")
      )

snapshotPCA3d(file="first_plot.png")
snapshotPCA3d(file="second_plot.png")
snapshotPCA3d(file="third_plot.png")

#Histograms
for (i in 1:(ncol(dat)-1)) {
    hist(dat[,c(i)], 
         xlab = capture.output(cat((colnames(dat)[i]),"values")),
         main=capture.output(cat("Histogram of",colnames(dat)[i], "per Risk Level"))
         )
}

#Distributions
for (i in 1:(ncol(dat)-1)) {
  plot <- ggplot(data=dat, aes(x = dat[,c(i)], fill = RiskLevel)) +
  geom_density(alpha = 0.5) +
  ggtitle(capture.output(cat("Distribution of",colnames(dat)[i], "per Risk Level"))) + 
  labs(x = capture.output(cat(colnames(dat)[i],"values")))
  
  print(plot)
  rm(plot)
}

#BoxPlots
for (i in 1:(ncol(dat)-1)) {
  bp <- ggplot(data = dat, aes(x = RiskLevel,y = dat[,c(i)], group = RiskLevel, fill = RiskLevel)) +
    geom_boxplot(show.legend = FALSE) +
    ggtitle(capture.output(cat("Boxplot of",colnames(dat)[i], "Variable per Risk Level"))) +
    ylab(capture.output(cat(colnames(dat)[i],"values")))
  
  print(bp)
  rm(bp)
}

#Percentage of low, mid and high risk individuals in the dataset
round(table(dat$RiskLevel)/length(dat$RiskLevel)*100,1)


#Normality Test (Shapiro-Wilk test for factor levels)
swAge<-byf.shapiro(Age~RiskLevel,data=dat)
swSystolicBP<-byf.shapiro(SystolicBP~RiskLevel,data=dat)
swDiastolicBP<-byf.shapiro(DiastolicBP~RiskLevel,data=dat)
swBS<-byf.shapiro(BS~RiskLevel,data=dat)
swBodyTemp<-byf.shapiro(BodyTemp~RiskLevel,data=dat)
swHeartRate<-byf.shapiro(HeartRate~RiskLevel,data=dat)

swAge
swSystolicBP
swDiastolicBP
swBS
swBodyTemp
swHeartRate

#Levene Test of Homogeneity of Variances
#The Levene test is less sensitive than the Bartlett test to departures from normality
tlevene1<-leveneTest(dat$Age~RiskLevel, data=dat)
tlevene2<-leveneTest(SystolicBP ~ RiskLevel, data=dat)
tlevene3<-leveneTest(DiastolicBP ~ RiskLevel, data=dat)
tlevene4<-leveneTest(BS ~ RiskLevel, data=dat)
tlevene5<-leveneTest(BodyTemp ~ RiskLevel, data=dat)
tlevene6<-leveneTest(HeartRate ~ RiskLevel, data=dat)

tlevene1$`Pr(>F)`
tlevene2$`Pr(>F)`
tlevene3$`Pr(>F)`
tlevene4$`Pr(>F)`
tlevene5$`Pr(>F)`
tlevene6$`Pr(>F)`

# Figner-Killeen Test of Homogeneity of Variances
# non-parametric test of Homogeneity of Variances
tfligner1<-fligner.test(Age ~ RiskLevel, data=dat)
tfligner2<-fligner.test(SystolicBP ~ RiskLevel, data=dat)
tfligner3<-fligner.test(DiastolicBP ~ RiskLevel, data=dat)
tfligner4<-fligner.test(BS ~ RiskLevel, data=dat)
tfligner5<-fligner.test(BodyTemp ~ RiskLevel, data=dat)
tfligner6<-fligner.test(HeartRate ~ RiskLevel, data=dat)

tfligner1$p.value
tfligner2$p.value
tfligner3$p.value
tfligner4$p.value
tfligner5$p.value
tfligner6$p.value


#Kruskal-Wallis
kwAgw <- kruskal.test(Age ~ RiskLevel, data = dat)
kwSystBP<-kruskal.test(SystolicBP ~ RiskLevel,data=dat)
kwDiasBP<-kruskal.test(DiastolicBP ~ RiskLevel,data=dat)
kwBS<-kruskal.test(BS~RiskLevel,data=dat)
kwBodyTemp<-kruskal.test(BodyTemp~RiskLevel,data=dat)
kwHR<-kruskal.test(HeartRate~RiskLevel,data=dat)

kwAgw$p.value
kwSystBP$p.value
kwDiasBP$p.value
kwBS$p.value
kwBodyTemp$p.value
kwHR$p.value

#Pairwise t-tests
ptAge <- with(dat,pairwise.t.test(dat$Age, RiskLevel,pool.sd=FALSE,p.adjust.method = "bonferroni"))
ptSystolicBP <- with(dat,pairwise.t.test(dat$SystolicBP, RiskLevel, pool.sd=FALSE, p.adjust.method = "bonferroni"))
ptDiastolicBP <- with(dat,pairwise.t.test(dat$DiastolicBP, RiskLevel, pool.sd=FALSE, p.adjust.method = "bonferroni"))
ptBS <- with(dat,pairwise.t.test(dat$BS, RiskLevel, pool.sd=FALSE, p.adjust.method = "bonferroni"))
ptBodyTemp <- with(dat,pairwise.t.test(dat$BodyTemp, RiskLevel, pool.sd=FALSE, p.adjust.method = "bonferroni"))
ptHeartRate <- with(dat,pairwise.t.test(dat$HeartRate, RiskLevel, pool.sd=FALSE, p.adjust.method = "bonferroni"))

ptAge
ptDiastolicBP
ptHeartRate
ptBodyTemp
ptSystolicBP
ptBS

ggcorrplot(cor(as.matrix(dat[,1:6])),hc.order = TRUE, 
           type = "lower",
           lab = TRUE,)

corr
dat$sd<-NULL
dat$od<-NULL

## split dat into train and test using stratified sampling
d <- rownames_to_column(dat, var = "id") %>% mutate_at(vars(id), as.integer)
training <- d %>% stratified(., group = "RiskLevel", size = 0.90)
dim(training)

## proportion check
prop.table(table(training$Species)) 

testing <- d[-training$id, ]
dim(testing)
prop.table(table(testing$RiskLevel)) 


## Modelling

##Data splitting and Decision tree was made in Python notebook
write.csv(dat,"dat_python.csv")

#In that notebook we created X_train, X_test, y_train, y_test and now we will import them:

X_train <- read.csv("C:/Users/Fernando/Desktop/Bio/X_train.csv")
X_test <- read.csv("C:/Users/Fernando/Desktop/Bio/X_test.csv")
y_train <- read.csv("C:/Users/Fernando/Desktop/Bio/y_train.csv")
y_test <- read.csv("C:/Users/Fernando/Desktop/Bio/y_test.csv")

#OLR and OPR

X_train$X <- NULL
X_test$X <- NULL
y_train$X <- NULL
y_test$X <- NULL

train <- cbind(X_train,y_train)
test <- cbind(X_test,y_test)

train$Age<-as.numeric(train$Age)
train$SystolicBP<-as.numeric(train$SystolicBP)
train$BS<-as.numeric(train$BS)
train$BodyTemp<-as.numeric(train$BodyTemp)
train$HeartRate<-as.numeric(train$HeartRate)

test$Age<-as.numeric(test$Age)
test$SystolicBP<-as.numeric(test$SystolicBP)
test$BS<-as.numeric(test$BS)
test$BodyTemp<-as.numeric(test$BodyTemp)
test$HeartRate<-as.numeric(test$HeartRate)


train$RiskLevel <- factor(train$RiskLevel, levels = c("low risk","mid risk", "high risk") ,ordered = TRUE)
test$RiskLevel <- factor(test$RiskLevel, levels = c("low risk","mid risk", "high risk") ,ordered = TRUE)

train_lm<-train
train_lm$RiskLevel<-as.numeric(train_lm$RiskLevel)

lm_regression <- lm(RiskLevel ~ Age + SystolicBP + DiastolicBP + BS + BodyTemp + HeartRate, 
                    data = train_lm)

vif <- vif(lm_regression)
vif
#conclude that we have to take off Diastolic

train <- train[,-3]

ordered_logistic <- polr(RiskLevel ~ Age + SystolicBP + BS + BodyTemp + HeartRate, data = train, Hess=TRUE)
summary(ordered_logistic)

ctable <- coef(summary(ordered_logistic))

## calculate and store p values
p <- pnorm(abs(ctable[, "t value"]), lower.tail = FALSE) * 2

## combined table
ctable <- cbind(ctable, "p value" = p)

profile_ci_olr <- confint(ordered_logistic, level = 0.95) # default method gives profiled CIs #####1
round(profile_ci_olr,4)

normal_ci_olr <- confint.default(ordered_logistic) # CIs assuming normality  ######2
round(normal_ci_olr,4)

# odds ratios
exp(coef(ordered_logistic))

## OR and CI
orci<- exp(cbind(OR = coef(ordered_logistic), normal_ci_olr))  #####3
round(orci,4)

typeof(train$RiskLevel)

myprobit <- polr(RiskLevel ~ Age + SystolicBP  + BS + BodyTemp + HeartRate,
                 data = train,method = "probit", Hess = TRUE)

## model summary
summary(myprobit)

ctable_probit <- coef(summary(myprobit))

#with(myprobit, pchisq(null.deviance - deviance, df.null - df.residual, lower.tail = FALSE))

## calculate and store p values
p_probit <- pnorm(abs(ctable_probit[, "t value"]), lower.tail = FALSE) * 2
ctable_probit <- cbind(ctable_probit, "p value" = p_probit)
ctable_probit

profile_ci_probit <- confint(myprobit, level = 0.95) # default method gives profiled CIs
normal_ci_probit <- confint.default(myprobit) # CIs assuming normality
round(profile_ci_probit,4)
round(normal_ci_probit,4)

brant(ordered_logistic)

#Predictions 
y_test$RiskLevel_numeric <- y_test$RiskLevel

y_test$RiskLevel_numeric[y_test$RiskLevel_numeric == "low risk"] = 1
y_test$RiskLevel_numeric[y_test$RiskLevel_numeric == "mid risk"] = 2
y_test$RiskLevel_numeric[y_test$RiskLevel_numeric == "high risk"] = 3


ord_log_pred <-  predict(ordered_logistic, X_test)
table(ord_log_pred, y_test$RiskLevel_numeric)

myprobit_pred <-  predict(myprobit, X_test)
table(myprobit_pred, y_test$RiskLevel_numeric)