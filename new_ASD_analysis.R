library(ggplot2)
library(dplyr)
library(rms)
library(twang)
library(MatchIt)
library(CBPS)
library(lubridate)
library(lmtest) #coeftest
library(sandwich) #vcovCL
library(cobalt)
library(dbarts)
library(xtable)

### Helper functions
compute_OR_matching <- function(m.out){
  m.data <- match.data(m.out) #matched data
  fit <- glm(Autism~ treatment, data = m.data, family = binomial(link = "logit"))
  ct <- coeftest(fit, vcov. = vcovCL, cluster = ~subclass)
  a <- matrix(exp(cbind(coef(fit), confint(ct)))[2,],1,3)
  return(cbind(a,ct[2,4]))
}

draw_ps_hist <- function(ps,treatment){
  datap <- data.frame(
    type = as.factor(treatment),
    value = ps)
  
  p <- datap %>%
    ggplot( aes(x=value, fill=type)) +
    geom_histogram( color="#e9ecef", alpha=0.6, position = 'identity',bins = 10) +
    scale_fill_manual(values=c("#69b3a2", "#404080")) +
    labs(fill="")+
    ylab("Frequency") +
    xlab("Propensity")
  return(p)
}


###Load Data
## -----------------------------------------------------------
data <- read.csv("data dryad.csv")
data <- subset(data, select= -c(Mom.Fake.ID,Dad.Fake.ID))
summary(data)
#Omit fathers younger than 15
ind <- which(data$Father.Age.At.Birth<15,arr.ind = TRUE)
data[ind,] <- NA
#Compute month of conception and age of conception of parents
data$ï..Month.of.Birth <- stringr::str_to_title(data$ï..Month.of.Birth)
months <- data$ï..Month.of.Birth
data$numeric_MOB <- match(months, month.name)
data$numberic_date <- paste0("2004-",as.character(data$numeric_MOB))
data$numberic_date <- as.Date(paste(data$numberic_date,"-15",sep=""))
data$numeric_DOC <- data$numberic_date - as.difftime(data$Gestational.Age, unit="weeks")
data$numeric_MOC2 <- month(as.POSIXlt(data$numeric_DOC, format="%Y-%m-%d"))
data$treatment <- ifelse(data$numeric_MOC2<=2 | data$numeric_MOC2>=11,1,0)
data$Mom.Age.At.Birth <- data$Mom.Age.At.Birth+0.5-data$Gestational.Age/(12*4.286)
data$Father.Age.At.Birth <- data$Father.Age.At.Birth+0.5-data$Gestational.Age/(12*4.286)

#Categorical and Ordinal variables
## -----------------------------------------------------------
new_data <- na.omit(subset(data, select=-c(ï..Month.of.Birth,numberic_date,numeric_DOC,numeric_MOB,numeric_MOC2,Birth.Weight.10.gr,Gestational.Age))) #keep only pre-treatment covariates
new_data$SES.Points <- factor(new_data$SES.Points, ordered = TRUE)
new_data$District.Number <- factor(new_data$District.Number)
colnames(new_data)[c(1,4)] <- c("Mom.Age.At.Conception","Father.Age.At.Conception")
summary(new_data)


#Add dummy variables
## -----------------------------------------------------------
design_matrix <- model.matrix(Autism~., new_data)[,-1]
d <- ncol(design_matrix)
new_data <- data.frame(cbind(design_matrix,new_data$Autism))
colnames(new_data)[d+1] <- "Autism"

#Pre-treatment balance
## -----------------------------------------------------------
summary(new_data)
# comparison between treatment and control groups - before matching or weigthing
m.out0 <- matchit(treatment~  Mom.Age.At.Conception + Father.Age.At.Conception + Number.of.birhts + SES.Points.L + SES.Points.Q + SES.Points.C
                  + SES.Points.4+ SES.Points.5+ SES.Points.6+ SES.Points.7+ SES.Points.8+ SES.Points.9+ SES.Points.10
                  + District.Number6 + District.Number7 + District.Number24 + District.Number30, data = new_data,
                  method = NULL, distance = "glm")
summary(m.out0)
plot(summary(m.out0))
love.plot(m.out0, stars="raw",stats = "mean.diffs", thresholds = c(m = .1),title="Unadjusted Covariate Balance",drop.distance = TRUE)
# ggsave("unadjusted_balance.pdf")
# ggsave("unadjusted_balance.eps")


### IPW methods

## 1. Simple Logistic Regression
## -----------------------------------------------------------
log_reg_model <- glm(treatment ~. -Autism , family=binomial(link='logit'),data=new_data )
summary(log_reg_model)
## probabilities ##
fitted.prob <- predict(log_reg_model ,newdata=new_data,type='response')
## Check overlap ##
p <- draw_ps_hist(ps=fitted.prob, treatment = new_data$treatment) +
    ggtitle("Propensity Estimates for Treated (1) vs. Control (0)")
print(p)

## Calibration ## 
val.prob(fitted.prob,new_data$treatment)

## 2.Quadratic Logistic Regression (including also interaction terms)
## -----------------------------------------------------------
log_reg_model2 <- glm(treatment ~.^2 , family=binomial(link='logit'),data=subset(new_data,select=-Autism))
summary(log_reg_model2)
## probabilities ##
fitted.prob2 <- predict(log_reg_model2 ,newdata=subset(new_data,select=-Autism),type='response')
## Check overlap ##
p <- draw_ps_hist(ps=fitted.prob2, treatment = new_data$treatment) +
  ggtitle("Quadratic LR Propensity Estimates for Treated (1) vs. Control (0)")
print(p)

## Calibration ## 
val.prob(fitted.prob2,new_data$treatment)

## 3.GBM
## -----------------------------------------------------------
propensity <- ps(formula= treatment ~.-Autism, data=new_data, estimand="ATT")
propensity3 <- propensity$ps$ks.mean.ATT
## Check overlap ##
p <- draw_ps_hist(ps=propensity3, treatment = new_data$treatment) +
  ggtitle("GBM1 Propensity Estimates for Treated (1) vs. Control (0)")
print(p)
## Calibration ## 
val.prob(propensity3,new_data$treatment)

propensity3b <- propensity$ps$es.mean.ATT
## Check overlap ##
p <- draw_ps_hist(ps=propensity3b, treatment = new_data$treatment) +
  ggtitle("GBM2 Propensity Estimates for Treated (1) vs. Control (0)")
print(p)
## Calibration ## 
val.prob(propensity3b,new_data$treatment)

## 4.CBPS
## -----------------------------------------------------------
cbps <-CBPS(treatment ~.-Autism, data=new_data,standardize=FALSE,method="exact",ATT=1)
propensity4 <- cbps$fitted.values
## Check overlap ##
p <- draw_ps_hist(ps=propensity4, treatment = new_data$treatment) +
  ggtitle("CBPS Propensity Estimates for Treated (1) vs. Control (0)")
print(p)
## Calibration ## 
val.prob(propensity4,new_data$treatment)


## -----------------------------------------------------------
## Compute ATT weights for each propensity model ##
weights1 <- new_data$treatment+fitted.prob*(1-new_data$treatment)/(1-fitted.prob)
weights2 <- new_data$treatment+fitted.prob2*(1-new_data$treatment)/(1-fitted.prob2)
weights3 <- new_data$treatment+propensity3*(1-new_data$treatment)/(1-propensity3)
weights3b <- new_data$treatment+propensity3b*(1-new_data$treatment)/(1-propensity3b)
weights4 <- new_data$treatment+propensity4*(1-new_data$treatment)/(1-propensity4)

## Check for balance ##
covs <- subset(new_data, select = -c(treatment, Autism))
# 1. Logistic regression
baltab <- bal.tab(covs, treat = new_data$treatment, weights = weights1, disp = c("means", "sds"), un = TRUE, 
                  stats = c("mean.diffs", "variance.ratios"))
love.plot(baltab, stats = "mean.diffs", thresholds = c(m = .1), stars = "raw", title="Covariate Balance after LR IPW",abs = TRUE)
# ggsave("LR_IPW_balance.pdf")
# ggsave("LR_IPW_balance.eps")

# 2. Quadratic logistic regression
baltab2 <- bal.tab(covs, treat = new_data$treatment, weights = weights2, disp = c("means", "sds"), un = TRUE, 
                   stats = c("mean.diffs", "variance.ratios"))
love.plot(baltab2, stats = "mean.diffs", thresholds = c(m = .1), stars = "raw",title="Covariate Balance after Quadratic LR IPW",abs = TRUE)
# ggsave("QLR_IPW_balance.pdf")
# ggsave("QLR_IPW_balance.eps")

# 3a. GBM1
baltab3a <- bal.tab(covs, treat = new_data$treatment, weights = weights3, disp = c("means", "sds"), un = TRUE, 
                    stats = c("mean.diffs", "variance.ratios"))
love.plot(baltab3a, stats = "mean.diffs", thresholds = c(m = .1), stars = "raw",title="Covariate Balance after GBM1 IPW",abs = TRUE)
# ggsave("GBM1_IPW_balance.pdf")
# ggsave("GBM1_IPW_balance.eps")

# 3b. GBM2
baltab3b <- bal.tab(covs, treat = new_data$treatment, weights = weights3b, disp = c("means", "sds"), un = TRUE, 
                    stats = c("mean.diffs", "variance.ratios"))
love.plot(baltab3b, stats = "mean.diffs", thresholds = c(m = .1), stars = "raw",title="Covariate Balance after GBM2 IPW",abs = TRUE)
# ggsave("GBM2_IPW_balance.pdf")
# ggsave("GBM2_IPW_balance.eps")

# 4. Covariate balancing propensity score
baltab4 <- bal.tab(covs, treat = new_data$treatment, weights = weights4, disp = c("means", "sds"), un = TRUE, 
                   stats = c("mean.diffs", "variance.ratios"))
love.plot(baltab4, stats = "mean.diffs", thresholds = c(m = .1), stars = "raw",title="Covariate Balance after CBPS IPW",abs = TRUE)
# ggsave("CBPS_IPW_balance.pdf")
# ggsave("CBPS_IPW_balance.eps")

## -----------------------------------------------------------
#trimmimg for the chosen propensity
new_data <- cbind(new_data,propensity4)
new_data_T <- subset(new_data, treatment==1)
new_data_C <- subset(new_data, treatment==0)
min_prop_T <- min(new_data_T$propensity4)
max_prop_T <- max(new_data_T$propensity4)
min_prop_C <- min(new_data_C$propensity4)
max_prop_C <- max(new_data_C$propensity4)
min_common_support <- max(min_prop_C,min_prop_T)
max_common_support <- min(max_prop_C,max_prop_T)
ind_trimmed <- which(new_data$propensity4>=min_common_support & new_data$propensity4<=max_common_support)
trimmed_data <- new_data[ind_trimmed,]

#plot propensity after trimming and chcek for better overlap
p <- draw_ps_hist(ps=trimmed_data$propensity4, treatment = trimmed_data$treatment) +
  ggtitle("Propensity Estimates after trimming for Treated (1) vs. Control (0)")
print(p)
summary(trimmed_data$propensity4)

## -----------------------------------------------------------
## Final ATT weights after trimming
propensity4 <- trimmed_data$propensity4
W=trimmed_data$treatment+propensity4*(1-trimmed_data$treatment)/(1-propensity4) #for ATT
fit4 <- glm(Autism~ treatment, data = trimmed_data, family = binomial(link = "logit"),weights = W)
summary(fit4)
ct4 <- coeftest(fit4, vcov. = vcovCL)
#OR
exp(coef(fit4))
#CI
exp(confint(ct4))
#both
a <- matrix(exp(cbind(coef(fit4), confint(ct4)))[2,],1,3)
IPW_cbps <- cbind(a,ct4[2,4]) 
IPW_df <- data.frame(matrix(NA,1,5))
IPW_df[1,2:5] <- data.frame(IPW_cbps)
IPW_df[1,1] <- "IPW on CBPS"
colnames(IPW_df) <- c("Method", "OR", "2.5% CI", "97.5% CI","p-value")
#print(xtable(IPW_df, type = "latex"), file = "IPW_table.tex")

### Matching
matching_df <- data.frame(matrix(NA, 10,5))
colnames(matching_df) <- c("Method", "OR", "2.5% CI", "97.5% CI","p-value")

## -----------------------------------------------------------
#nn matching on logit propensity score
m.out1 <- matchit(treatment~  Mom.Age.At.Conception + Father.Age.At.Conception + Number.of.birhts + SES.Points.L + SES.Points.Q + SES.Points.C
                  + SES.Points.4+ SES.Points.5+ SES.Points.6+ SES.Points.7+ SES.Points.8+ SES.Points.9+ SES.Points.10
                  + District.Number6 + District.Number7 + District.Number24 + District.Number30, data = new_data,
                  method = "nearest", distance = "glm",discard = "both")
##assess balance
summary(m.out1,standardize = TRUE)
plot(summary(m.out1,standardize = TRUE))
love.plot(m.out1, binary = "std", thresholds = c(m = .1), drop.distance = TRUE,title = "Covariate Balance after NN Matching on Logit Propensity")
# ggsave("m1.eps")
# ggsave("m1.pdf")
#Estimate OR
matching_df[1,2:5] <- compute_OR_matching(m.out1)
matching_df[1,1] <- "NN on logit propensity"

#1-nn matching on probit propensity score
m.out3 <- matchit(treatment~  Mom.Age.At.Conception + Father.Age.At.Conception + Number.of.birhts + SES.Points.L + SES.Points.Q + SES.Points.C
                  + SES.Points.4+ SES.Points.5+ SES.Points.6+ SES.Points.7+ SES.Points.8+ SES.Points.9+ SES.Points.10
                  + District.Number6 + District.Number7 + District.Number24 + District.Number30, data = new_data,
                  method = "nearest", distance = "glm",link="probit",discard = "both")
##assess balance
summary(m.out3,standardize = TRUE)
plot(summary(m.out3,standardize = TRUE))
love.plot(m.out3,binary = "std", thresholds = c(m = .1), drop.distance = TRUE,title = "Covariate Balance after NN Matching on Probit Propensity")
# ggsave("m3.eps")
# ggsave("m3.pdf")
#Estimate OR
matching_df[2,2:5] <- compute_OR_matching(m.out3)
matching_df[2,1] <- "NN on probit propensity"


#1-nn matching on gam propensity score
m.out4 <- matchit(treatment~  Mom.Age.At.Conception + Father.Age.At.Conception + Number.of.birhts + SES.Points.L + SES.Points.Q + SES.Points.C
                  + SES.Points.4+ SES.Points.5+ SES.Points.6+ SES.Points.7+ SES.Points.8+ SES.Points.9+ SES.Points.10
                  + District.Number6 + District.Number7 + District.Number24 + District.Number30, data = new_data,
                  method = "nearest", distance = "gam",discard = "both")
##assess balance
summary(m.out4,standardize = TRUE)
plot(summary(m.out4,standardize = TRUE))
love.plot(m.out4,binary = "std", thresholds = c(m = .1), drop.distance = TRUE,title = "Covariate Balance after NN Matching on GAM Propensity")
# ggsave("m4.eps")
# ggsave("m4.pdf")
#Estimate OR
matching_df[3,2:5] <- compute_OR_matching(m.out4)
matching_df[3,1] <- "NN on GAM propensity"


#1-nn matching on random forest propensity score
m.out5 <- matchit(treatment~  Mom.Age.At.Conception + Father.Age.At.Conception + Number.of.birhts + SES.Points.L + SES.Points.Q + SES.Points.C
                  + SES.Points.4+ SES.Points.5+ SES.Points.6+ SES.Points.7+ SES.Points.8+ SES.Points.9+ SES.Points.10
                  + District.Number6 + District.Number7 + District.Number24 + District.Number30, data = new_data,
                  method = "nearest", distance = "randomforest",discard = "both")
##assess balance
summary(m.out5,standardize = TRUE)
plot(summary(m.out5,standardize = TRUE))
love.plot(m.out5,binary = "std", thresholds = c(m = .1), drop.distance = TRUE,title = "Covariate Balance after NN Matching on random forest Propensity")
# ggsave("m5.eps")
# ggsave("m5.pdf")
#Estimate OR
matching_df[10,2:5] <- compute_OR_matching(m.out5)
matching_df[10,1] <- "NN on random forest propensity"


#1-nn matching on GBM propensity score
m.out2 <- matchit(treatment~  Mom.Age.At.Conception + Father.Age.At.Conception + Number.of.birhts + SES.Points.L + SES.Points.Q + SES.Points.C
                  + SES.Points.4+ SES.Points.5+ SES.Points.6+ SES.Points.7+ SES.Points.8+ SES.Points.9+ SES.Points.10
                  + District.Number6 + District.Number7 + District.Number24 + District.Number30, data = new_data,
                  method = "nearest", distance = propensity3,discard = "both")
##assess balance
summary(m.out2)
plot(summary(m.out2))
love.plot(m.out2,binary = "std", thresholds = c(m = .1), drop.distance = TRUE,title = "Covariate Balance after NN Matching on GBM Propensity")
# ggsave("m2.pdf")
# ggsave("m2.eps")
#Estimate OR
matching_df[9,2:5] <- compute_OR_matching(m.out2)
matching_df[9,1] <- "NN on GBM propensity"


#1-NN matching on mahalanobis distance
m.out7 <- matchit(treatment ~  Mom.Age.At.Conception + Father.Age.At.Conception + Number.of.birhts + SES.Points.L + SES.Points.Q + SES.Points.C
                  + SES.Points.4+ SES.Points.5+ SES.Points.6+ SES.Points.7+ SES.Points.8+ SES.Points.9+ SES.Points.10
                  + District.Number6 + District.Number7 + District.Number24 + District.Number30, data = new_data,
                  method = "nearest", distance = "mahalanobis")
##assess balance
summary(m.out7)
plot(summary(m.out7))
love.plot(m.out7,binary = "std", thresholds = c(m = .1), drop.distance = TRUE, title = "Covariate Balance after NN Matching on Mahalanobis Distance")
# ggsave("m7.pdf")
# ggsave("m7.eps")
#Estimate OR
matching_df[5,2:5] <- compute_OR_matching(m.out7)
matching_df[5,1] <- "NN on mahalanobis"


#1-NN matching on CBPS propensity
m.out8 <- matchit(treatment ~  Mom.Age.At.Conception + Father.Age.At.Conception + Number.of.birhts + SES.Points.L + SES.Points.Q + SES.Points.C
                  + SES.Points.4+ SES.Points.5+ SES.Points.6+ SES.Points.7+ SES.Points.8+ SES.Points.9+ SES.Points.10
                  + District.Number6 + District.Number7 + District.Number24 + District.Number30, data = new_data,
                  method = "nearest", distance = "cbps")
##assess balance
summary(m.out8)
plot(summary(m.out8))
love.plot(m.out8,binary = "std", thresholds = c(m = .1), drop.distance = TRUE,title = "Covariate Balance after NN Matching on CBPS propensity")
# ggsave("m8.pdf")
# ggsave("m8.eps")
#Estimate OR
matching_df[4,2:5] <- compute_OR_matching(m.out8)
matching_df[4,1] <- "NN on CBPS propensity"

matching_df_new <- matching_df[1:5,]
#print(xtable(matching_df_new, type = "latex"), file = "matching_table.tex")