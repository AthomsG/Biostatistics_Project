library("psych")
library("survival")
library("ggplot2")
library("ggfortify")
library("ranger")
library("dplyr")
library("brms")
library("flexsurv")
library("rstan")
library("tidyverse")
library("modelr")
library("rjags")
library("functional")
library("coda")
library("sm")
library("loo")
library("bayesplot")
library('corrplot')
library("PerformanceAnalytics")
library("ggcorrplot")
library("lares")
library("psych")
patient_info <- read.table("/Users/davidsousa/Downloads/ProjectCS/Report3_Dataset.txt")

names(patient_info) <- c("treat", "celltype","time","status","karno","age")

vet <- mutate(patient_info,
              treat = factor(treat,labels=c("standard","test")),
              celltype = factor(celltype, labels =c("squamous", "small cell", "adeno", "large")))

#PARTE CLÁSSICA
corrplot(cor(patient_info[c(3,5,6)]))
corrplot.mixed(cor(patient_info[c(3,5,6)]),
               upper = "square",
               lower = "number",
               addgrid.col = "black",
               tl.col = "black")
chart.Correlation(cor(patient_info[c(3,5,6)]), histogram = TRUE, method = "pearson")
ggcorrplot(cor(patient_info[c(3,5,6)]), 
           method = "circle",
           type = "lower",
           outline.color = "black",
           lab_size = 6)
corr_cross(patient_info[c(3,5,6)], rm.na = T, max_pvalue = 0.43, top = 3, grid = T)

pairs.panels(corr(patient_info[c(3,5,6)]),
             smooth = TRUE,      
             scale = FALSE,      
             density = TRUE,     
             ellipses = TRUE,   
             method = "pearson", 
             pch = 21,          
             lm = FALSE,         
             cor = TRUE,         
             jiggle = FALSE,     
             factor = 2,         
             hist.col = 4,       
             stars = TRUE,       
             ci = TRUE)  

boxplot(time~treat,
        data=vet,
        main="Different boxplots for each type of treatment",
        xlab="Treatment",
        ylab="Survival Time",
        col="blue",
        border="brown"
)
mtext(paste("Outliers: ", paste(out, collapse = ", ")))

km <- with(patient_info,Surv(time,status))

km_fit <- survfit(Surv(time, status) ~ treat, data = vet)

autoplot(km_fit)

km_trt_fit <- survfit(Surv(time, status) ~ celltype, data=vet)

autoplot(km_trt_fit)

km_ST_fit <- survfit(Surv(time,status) ~ ST, data = vet)

autoplot(km_ST_fit) 

cox <- coxph(Surv(time, status) ~ treat + celltype + karno + age  , data = vet)

summary(cox)

cox_fit <- survfit(cox)

autoplot(cox_fit)

aa_fit <-aareg(Surv(time, status) ~ treat + celltype +
                 karno  + age , 
               data = vet)
summary(aa_fit)

autoplot(aa_fit)

wei_fit_reg <- survreg(Surv(time, status) ~ treat + celltype + karno + age  , data = vet, dist = "weibull")

autoplot(wei_fit_reg)

surv <- seq(.99, .01, by = -.01)
t <- predict(wei_fit_reg, type = "quantile", p = 1 - surv)
surv_wb <- data.frame(time = t, surv = surv, 
                      upper = NA, lower = NA, std.err = NA)

wei_fit_2 <- flexsurvreg(Surv(time, status) ~ treat + celltype + karno + age  , data = vet, dist = "weibull")

plot(wei_fit_2)

#PARTE BAYESIANA

#CENSURA DOS TEMPOS COM INDICADOR DE MORTE
cens <- matrix(c(vet$time, rep(NA, length(vet$time))),
               nrow = length(vet$time), ncol = 2)

vet$time[vet$status == 0 ] <- NA
is.censored <-as.numeric(is.na(vet$time))


#MATRIX COM OS DADOS EM QUE OS FATORES FORAM TRANSFORMADOS EM BINÁRIO
X <- model.matrix(~ treat + celltype + karno + age, data = vet)

#CRIAÇÃO DO MODELO JAGS
d.jags <- list(n = nrow(vet), time = vet$time, cens = cens, X = X, is.censored = is.censored, Nbetas = ncol(X))
i.jags <- function() {list(beta = rnorm(ncol(X)), alpha = runif(1))}
p.jags <- c("beta","alpha","loglh")

m1 <- jags.model(data = d.jags, file = "/Users/davidsousa/Downloads/ProjectCS/AFT.txt", inits = i.jags, n.chains = 2)

update(m1,1000)

#CRIAÇÃO DOS DADOS PARA FAZER OS TRACES ASSOCIADOS A QUESTÕES DE CONVERGÊNCIA
res <- coda.samples(m1,variable.names = p.jags, n.iter = 10000,thin = 2)
results_in_matrix <- as.matrix(res)
error <-jags.samples(m1, c("pD", "deviance","loglh"), 10000, type="mean")
summary(res)
gelman.diag(res)
heidel.diag(res)
raftery.diag(res)
par(mfrow = c(2,4))
densplot(res[,c(1:8)])
traceplot(res[,c(1:8)])
dic.samples(m1,1000,type="pD")

#Geração das posterior samples
result <- as.mcmc(do.call(rbind,res))

alpha <- result[,1] ; beta1 <- result[,2] ; beta2 <- result[,3] ; beta3 <- result[,4]
beta4 <- result[,5]; beta5 <- result[,6] ; beta6 <- result[,7]; beta7 <- result[,8]; loglh <- result[,9]

#Cálculo apróximado de AIC e BIC e cálculo preciso de DIC

AIC_1 <- 2*8 - 2*max(loglh)
BIC_1 <- 8*log(137) - 2*max(loglh)
DIC_1 <- -4*mean(loglh) + 2*likeli_on_mean_1

likeli_on_mean_1 <- 0
for (i in 1:137) {
  if (is.numeric(patient_info[i,][4]) == 1){
    likeli_on_mean_1 <- likeli_on_mean_1 + log(DensityF(as.numeric(patient_info[i,][3]),mean(alpha),mean(beta1),mean(beta2),mean(beta3),mean(beta4),mean(beta5),mean(beta6),mean(beta7),
                                                        as.numeric(X[i,][2]),as.numeric(X[i,][3]),as.numeric(X[i,][4]),as.numeric(X[i,][5]),as.numeric(X[i,][6]),as.numeric(X[i,][7])))
  } else {
    likeli_on_mean_1 <- likeli_on_mean_1 + log(SurvivalF(as.numeric(patient_info[i,][3]),mean(alpha),mean(beta1),mean(beta2),mean(beta3),mean(beta4),mean(beta5),mean(beta6),mean(beta7),
                                                         as.numeric(X[i,][2]),as.numeric(X[i,][3]),as.numeric(X[i,][4]),as.numeric(X[i,][5]),as.numeric(X[i,][6]),as.numeric(X[i,][7])))
  }
}



#Vamos considerar um modelo mais simples em que removemos a idade

X_age<- model.matrix(~ treat + celltype + karno, data = vet)

d.jags_age <- list(n = nrow(vet), time = vet$time, cens = cens, X = X_age, is.censored = is.censored, Nbetas = ncol(X_age))
i.jags_age <- function() {list(beta = rnorm(ncol(X_age)), alpha = runif(1))}
p.jags_age <- c("beta","alpha","loglh")

m3 <- jags.model(data = d.jags_age, file = "/Users/davidsousa/Downloads/ProjectCS/AFT.txt", inits = i.jags_age, n.chains = 2)

update(m3,1000)

res_age <- coda.samples(m3,variable.names = p.jags_age, n.iter = 10000,thin = 2)
summary(res_agel)
gelman.diag(res_age)
heidel.diag(res_age)
raftery.diag(res_age)
par(mfrow = c(2,4))
densplot(res_age)
traceplot(res_age[,c(1:7)])

result_2 <- as.mcmc(do.call(rbind,res_age))

# Resultado das distribuições posteriori
alpha_2 <- result_2[,1] ; beta1_2 <- result_2[,2] ; beta2_2 <- result_2[,3] ; beta3_2 <- result_2[,4]
beta4_2 <- result_2[,5]; beta5_2 <- result_2[,6] ; beta6_2 <- result_2[,7]; loglh_2 <- result_2[,8]; 

#Cálculo aproximado de AIC e BIC 

AIC_2 <- 2*7 - 2*max(loglh_2)
BIC_2 <- 7*log(137) - 2*max(loglh_2)
DIC_2 <- -4*mean(loglh_2) + 2*likeli_on_mean_2
  

likeli_on_mean_2 <- 0
for (i in 1:137) {
  if (is.numeric(patient_info[i,][4]) == 1){
    likeli_on_mean_2 <- likeli_on_mean_2 + log(DensityF(as.numeric(patient_info[i,][3]),mean(alpha_2),mean(beta1_2),mean(beta2_2),mean(beta3_2),mean(beta4_2),mean(beta5_2),mean(beta6_2),0,
                                                        as.numeric(X[i,][2]),as.numeric(X[i,][3]),as.numeric(X[i,][4]),as.numeric(X[i,][5]),as.numeric(X[i,][6]),as.numeric(X[i,][7])))
  } else {
    likeli_on_mean_2 <- likeli_on_mean_2 + log(SurvivalF(as.numeric(patient_info[i,][3]),mean(alpha_2),mean(beta1_2),mean(beta2_2),mean(beta3_2),mean(beta4_2),mean(beta5_2),mean(beta6_2),0,
                                                         as.numeric(X[i,][2]),as.numeric(X[i,][3]),as.numeric(X[i,][4]),as.numeric(X[i,][5]),as.numeric(X[i,][6]),as.numeric(X[i,][7])))
  }
}

#Vamos considerar um modelo ainda mais simples, onde se remove a idade e reagrupamos as células

vet_re <- patient_info
for (i in 1:137) {
  if(vet_re$celltype[i] == 1 | vet_re$celltype[i]==4){
    vet_re$celltype[i] <- 'C_1'
  }
  else if( vet_re$celltype[i] == 2 | vet_re$celltype[i]==3){
    vet_re$celltype[i] <- 'C_2'
  }
}

vet_re_1 <- mutate(vet_re,
                   treat = factor(treat,labels=c("standard","test")),
                   celltype = factor(celltype, labels =c("C_1", "C_2")))
X_age_cellgroup <- model.matrix(~treat + celltype + karno, data = vet_re_1)


d.jags_age_cellgroup <- list(n = nrow(vet_re_1), time = vet$time, cens = cens, X = X_age_cellgroup, is.censored = is.censored, Nbetas = ncol(X_age_cellgroup))
i.jags_age_cellgroup <- function() {list(beta = rnorm(ncol(X_age_cellgroup)), alpha = runif(1))}
p.jags_age_cellgroup <- c("beta","alpha","loglh","time")

m4 <- jags.model(data = d.jags_age_cellgroup, file = "/Users/davidsousa/Downloads/ProjectCS/AFT.txt", inits = i.jags_age_cellgroup, n.chains = 2)

update(m4,1000)

res_age_cellgroup <- coda.samples(m4,variable.names = p.jags_age_cellgroup, n.iter = 10000,thin = 2)
summary(res_age_cellgroup)
gelman.diag(res_age_cellgroup)
heidel.diag(res_age_cellgroup)
raftery.diag(res_age_cellgroup)
par(mfrow = c(2,4))
densplot(res_age_cellgroup[,1:5])
traceplot(res_age_cellgroup[,1:5])

result_3 <- as.mcmc(do.call(rbind,res_age_cellgroup))


#posteriors
alpha_3 <- result_3[,1] ; beta1_3 <- result_3[,2] ; beta2_3 <- result_3[,3] ; beta3_3 <- result_3[,4]
beta4_3 <- result_3[,5]; loglh_3 <- result_3[,6] ; 

#Cálculo aproximado de AIC, BIC e cálculo preciso de DIC
AIC_3 <- 2*5 - 2*max(loglh_3)
BIC_3 <- 5*log(137) - 2*max(loglh_3)
DIC_3 <- -4*mean(loglh_3) + 2*likeli_on_mean_3

likeli_on_mean_3 <- 0
for (i in 1:137) {
  if (is.numeric(patient_info[i,][4]) == 1){
    likeli_on_mean_3 <- likeli_on_mean_3 + log(DensityF(as.numeric(patient_info[i,][3]),mean(alpha_3),mean(beta1_3),mean(beta2_3),mean(beta3_3),0,0,mean(beta4_3),0,
                                                        as.numeric(X[i,][2]),as.numeric(X[i,][3])+as.numeric(X[i,][4]),0,0,as.numeric(X[i,][6]),0))
  } else {
    likeli_on_mean_3 <- likeli_on_mean_3 + log(SurvivalF(as.numeric(patient_info[i,][3]),mean(alpha_3),mean(beta1_3),mean(beta2_3),mean(beta3_3),0,0,mean(beta4_3),0,
                                                         as.numeric(X[i,][2]),as.numeric(X[i,][3])+as.numeric(X[i,][4]),0,0,as.numeric(X[i,][6]),0))
  }
}

#Predictions

#ind iguais com células C_2 contra células C_1

Rm_C2_vs_C1 <- exp(beta3_3)
summary(Rm_C2_vs_C1)
plot(Rm_C2_vs_C1)

#ind iguais , os primeiros com tratamento especial e os segundos com o usual

RM_special_vs_usual <- exp(beta2_3)
summary(RM_special_vs_usual)


#ind iguais do grupo C1 com tratamento especial vs C2 com tratamento usual 

Rm_small_vs_large <- exp(beta2_3 -beta3_3)
summary(Rm_small_vs_large)



SurvivalF <- function(t,alpha, beta1,beta2,beta3,beta4,X1,X2,X3) {
  out <- exp(- t^(alpha)*exp(-1*alpha*(beta1 + beta2*X1 + beta3*X2 + beta4*X3 )))
  return(out)
}

DensityF <- function(t,alpha, beta1,beta2,beta3,beta4,X1,X2,X3) {
  out <- alpha*exp(-1*alpha*(beta1 + beta2*X1 + beta3*X2 + beta4*X3 ))*t^(alpha)*exp(- t^(alpha)*exp(-1*alpha*(beta1 + beta2*X1 + beta3*X2 + beta4*X3)))
  return(out)
}

MedianF <- function(alpha,beta1,beta2,beta3,beta4,X1,X2,X3) {
  out <- log(2)^(1/alpha)/exp(-1*(beta1 + beta2*X1 + beta3*X2 + beta4*X3 ))
  return(out)
}



# Vamos dar plot com a média dos betas e do alpha


Time <-  seq(0,1500, by = 0.5)

#Cálculo de vários casos, onde variamos 

caso_1 <- SurvivalF(Time,mean(alpha_3),mean(beta1_3),mean(beta2_3),mean(beta3_3),mean(beta4_3),1,0,40)
caso_2 <-SurvivalF(Time,mean(alpha_3),mean(beta1_3),mean(beta2_3),mean(beta3_3),mean(beta4_3),0,0,40)

caso_3 <- SurvivalF(Time,mean(alpha_3),mean(beta1_3),mean(beta2_3),mean(beta3_3),mean(beta4_3),1,0,60)
caso_4 <-SurvivalF(Time,mean(alpha_3),mean(beta1_3),mean(beta2_3),mean(beta3_3),mean(beta4_3),0,0,60)

caso_5 <- SurvivalF(Time,mean(alpha_3),mean(beta1_3),mean(beta2_3),mean(beta3_3),mean(beta4_3),1,0,80)
caso_6 <-SurvivalF(Time,mean(alpha_3),mean(beta1_3),mean(beta2_3),mean(beta3_3),mean(beta4_3),0,0,80)

plot(Time,caso_1 ,
     main="",
     ylab="Survival probability",
     type="l",
     col=1)
lines(Time,caso_2, col=2)
lines(Time,caso_3, col=3)
lines(Time,caso_4, col=4)
lines(Time,caso_5, col=5)
lines(Time,caso_6, col=6)
legend("topright",
       c("Experimental(40)","Usual(40)","Experimental(60)","Usual(60)","Experimental(80)","Usual(80)"),
       fill=1:6
)

Median_1 <- MedianF(mean(alpha_3),mean(beta1_3),mean(beta2_3),mean(beta3_3),mean(beta4_3),1,0,40)
Median_2 <- MedianF(mean(alpha_3),mean(beta1_3),mean(beta2_3),mean(beta3_3),mean(beta4_3),0,0,40)
Median_3 <- MedianF(mean(alpha_3),mean(beta1_3),mean(beta2_3),mean(beta3_3),mean(beta4_3),1,0,60)
Median_4 <- MedianF(mean(alpha_3),mean(beta1_3),mean(beta2_3),mean(beta3_3),mean(beta4_3),0,0,60)
Median_5 <- MedianF(mean(alpha_3),mean(beta1_3),mean(beta2_3),mean(beta3_3),mean(beta4_3),1,0,80)
Median_6 <- MedianF(mean(alpha_3),mean(beta1_3),mean(beta2_3),mean(beta3_3),mean(beta4_3),0,0,80)

#40/60/80 karno , células grupo C2 com e sem tratamento especial

caso_7 <- SurvivalF(Time,mean(alpha_3),mean(beta1_3),mean(beta2_3),mean(beta3_3),mean(beta4_3),1,1,40)
caso_8 <-SurvivalF(Time,mean(alpha_3),mean(beta1_3),mean(beta2_3),mean(beta3_3),mean(beta4_3),0,1,40)

caso_9 <- SurvivalF(Time,mean(alpha_3),mean(beta1_3),mean(beta2_3),mean(beta3_3),mean(beta4_3),1,1,60)
caso_10 <-SurvivalF(Time,mean(alpha_3),mean(beta1_3),mean(beta2_3),mean(beta3_3),mean(beta4_3),0,1,60)

caso_11 <- SurvivalF(Time,mean(alpha_3),mean(beta1_3),mean(beta2_3),mean(beta3_3),mean(beta4_3),1,1,80)
caso_12 <-SurvivalF(Time,mean(alpha_3),mean(beta1_3),mean(beta2_3),mean(beta3_3),mean(beta4_3),0,1,80)

plot(Time,caso_7 ,
     main="",
     ylab="Survival probability",
     type="l",
     col=1)
lines(Time,caso_8, col=2)
lines(Time,caso_9, col=3)
lines(Time,caso_10, col=4)
lines(Time,caso_11, col=5)
lines(Time,caso_12, col=6)
legend("topright",
       c("Experimental(40)","Usual(40)","Experimental(60)","Usual(60)","Experimental(80)","Usual(80)"),
       fill=1:6
)

Median_7 <- MedianF(mean(alpha_3),mean(beta1_3),mean(beta2_3),mean(beta3_3),mean(beta4_3),1,1,40)
Median_8 <- MedianF(mean(alpha_3),mean(beta1_3),mean(beta2_3),mean(beta3_3),mean(beta4_3),0,1,40)
Median_9 <- MedianF(mean(alpha_3),mean(beta1_3),mean(beta2_3),mean(beta3_3),mean(beta4_3),1,1,60)
Median_10 <- MedianF(mean(alpha_3),mean(beta1_3),mean(beta2_3),mean(beta3_3),mean(beta4_3),0,1,60)
Median_11 <- MedianF(mean(alpha_3),mean(beta1_3),mean(beta2_3),mean(beta3_3),mean(beta4_3),1,1,80)
Median_12 <- MedianF(mean(alpha_3),mean(beta1_3),mean(beta2_3),mean(beta3_3),mean(beta4_3),0,1,80)


