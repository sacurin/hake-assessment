# Title: Empirical weight at age pacific hake MSE
# Purpose: Incorporate variability in year and cohort in the projection period of hake
# Author: Sandra Curin-Osorio, Kristin Marshall, Aaron Berger & Kelli Johnson
# method: Fit a non spatial temporal generalized linear mixed effects model (GLMM) with sdmTMB
# Project: Hake MSE
# Date: 5/3/2024

library(gfdata)
library(gfplot)
library(ggmagnify)
library(hakedataUSA)
library(rnaturalearthhires)
library(remotes)
library(devtools)
library(glmm)
library(sdmTMB)
library(latex2exp)
library(tidyverse)
library(ggeffects)
library(visreg)
library(mgcv)
library(tidymv)
library(fs)
library(ggpubr)
library(ggeffects)
library(broom)
#devtools::install_github("pacific-hake/hake-assessment") Didn't work
#devtools::install("C:/GitHub/hake-assessment", dependencies = TRUE)


# The first part of the code checks if the GAM we will use generates different results when run twice.
# The second part of the code incorporate the effects or variability of year, cohort, and their combination into the empirical weight-at-age during the projection periods.
# The third part of the code produce 100 simulation of year+cohort effect on the weight of pacific hake.
# The fourth part of the code involves plotting the predicted weights for the three EWAA matrices.


rm(list = c("plot_map"))
# Load the WAA matrix (only use fishery data)/where your hake-assessment/vignettes/setup-bridge-models.Rmd (Kelli code)line 458
setwd("C:\\GitHub\\hake-assessment\\doc")
devtools::load_all()

#######################################################################################################################
# First part
# 1- Verify whether running the same model give different results and diagnostics upon subsequent executions
# The code used a sdmTMB to replicates the empirical weight at age used on the hake assessment. In the stock assessment
# the empirical weight at age is modeled using GAM. The code predict the weight given the age
# M1: w~s(age)+(1|fcohort)+(1|fyear)+sex, we used the same GAM from the stock assessment to be consistent (Kelli Johnson)
########################################################################################################################
m1 <- sdmTMB::sdmTMB(
  data = weight_at_age_df  |>
    dplyr::filter(sex != "U") |>
    dplyr::mutate(
      age = ifelse(age > 15, 15, age),
      cohort = year - age,
      fyear = as.factor(year),
      fcohort = as.factor(cohort)
    ),
  formula = weight ~ 1 + s(age) + (1|fcohort) + (1|fyear) + sex,
  family = lognormal(link = "log"),
  spatial = "off",
  spatiotemporal = "off",
  time = "year",
  control = sdmTMB::sdmTMBcontrol(newton_loops = 1)
)

#Print the m1 model fit
print(m1)

# Function to calculate the variance-covariance matrix of m1
vcov.sdmTMB <- function(object, complete = FALSE, ...) {
  sdr <- object$sd_report
  v <- sdr$cov.fixed
  fe <- tidy(object)$term
  nm <- colnames(v)
  i <- grepl("^b_j$", nm)
  if (sum(i)) {
    if (sum(i) == length(fe)) { # should always be true
      nm[i] <- fe
    }
  }
  colnames(v) <- nm
  rownames(v) <- nm
  if (isTRUE(complete)) {
    return(v)
  } else {
    return(v[i,i,drop=FALSE])
  }
}
# Extracting covariance of m1
random_vcov <- vcov.sdmTMB(m1)


## Extracting fixed effect parameters of m1 as a data frame:
fixed_eff_m1<-tidy(m1,conf.int=TRUE)


#Dispersion parameters for random effect of m1
random_par_m1 <- tidy(m1,"ran_pars",conf.int=TRUE)

# term    estimate std.error conf.low conf.high
# <chr>      <dbl>     <dbl>    <dbl>     <dbl>
#   1 phi       0.261   0.000375   0.260      0.262 # A dispersion parameter for a distribution
# 2 sigma_G   0.199   0.0186     0.166      0.239   # IID random intercept variance (fcohort)
# 3 sigma_G   0.0986  0.0103     0.0803     0.121   # IID random intercept variance (fyear)

#random effect values by variable
random_eff_val_m1<-tidy(m1,"ran_vals",conf.int=TRUE)

# Extracting dispersion parameter of m1
dispersion_param_m1 <- 0.26
print(dispersion_param_m1)

# Extracting ML criterion at convergence of m1
ML_criterion_m1 <- -154961.173  #this value is coming from the output
print(ML_criterion_m1)


#Diagnostic of m1 using sanity function of sdmTMB.sanity checks include:
# i-  Verifying that the model converged properly,
# ii- Checking if the estimated parameters are within expected ranges
sanity(m1)
predictions_m1 <- predict(m1)
head(predictions_m1)
predictions_m1$resids <- residuals(m1) # randomized quantile residuals
resid<-ggplot(predictions_m1, aes(age, weight, col = resids)) + scale_colour_gradient2() +
  geom_point() + facet_wrap(~year)
hist(predictions_m1$resids)


## Residual of m1
residuals_m1 <- residuals(m1)
residuals_m1 <- residuals_m1[is.finite(residuals_m1)]
#r <- residuals(m1, "mle-mvn", mcmc_warmup = 100, mcmc_iter = 101)

png(fs::path("C:/GitHub/hake-assessment/Tables_figures_Sandra/weight-at-age-qq_m1.png"))
qqnorm(residuals_m1); qqline(residuals_m1)
dev.off()
saveRDS(
  m1,
  file = fs::path("C:/GitHub/hake-assessment/Tables_figures_Sandra/weight-at-age-model_m1.rds")
)


# create prediction grid of m1
pred_grid1 <- expand.grid(
  year = unique(m1[["data"]][["year"]]),
  sex = unique(m1[["data"]][["sex"]]),
  age = 0:15
) |>
  dplyr::mutate(cohort = year - age) |>
  dplyr::filter(cohort %in% unique(m1[["data"]][["cohort"]])) |>
  dplyr::mutate(
    fyear = as.factor(year),
    fcohort = as.factor(cohort)
  )

# Get estimates of m1
preds1 <- predict(m1, newdata = pred_grid1)
preds1 <- preds1 |>
  dplyr::mutate(est_weight = exp(est))
saveRDS(preds,
        file = fs::path("C:/GitHub/hake-assessment/Tables_figures_Sandra/weight-at-age-preds_m1.rds")
)


# Plot of predicted weight of m1
pr_m1<-ggplot(preds, aes(age, est_weight, colour = sex)) + geom_point() +
  facet_wrap(~sex)+labs(y="Estimated weight")

# make predicts EWAA
ewaa_m1 <- preds1 |>
  dplyr::group_by(year, age) |>
  dplyr::summarise(
    n = dplyr::n(),
    pred_weight = mean(exp(est)),
    upper=quantile(exp(est),probs = 0.975),
    lower = quantile(exp(est), probs = 0.025),
    run="1"
  ) |>
  as.data.frame() |>
  dplyr::select(year, age, pred_weight,upper,lower,run)


## Visualize predict weight at age by year
ewaa_m1<-ggplot(ewaa_m1, aes(age, pred_weight)) + geom_point()+labs(x = "Age", y = "Predict weight")

## Visualize predict weight at age by year
weight_OP1<-ggplot(ewaa_m1, aes(x = factor(year), y = pred_weight, color = factor(age),
                                group = factor(age))) +
  geom_text(aes(label=round(age,2)), size = 4.5) +
  geom_line(alpha = 0.85) +ylim(0, 2.5)+
  theme_bw() +
  labs(x = "Year", y = "Weight") + geom_ribbon(data=ewaa_m1,aes(ymin=upper, ymax=lower,fill=factor(age)),alpha=0.1, linetype = "dotted") +
  theme(axis.title = element_text(size = 17),
        axis.text = element_text(size = 15),
        legend.title = element_text(size = 17),
        legend.text = element_text(size = 15),
        legend.position = "none",
        axis.title.x = element_blank(),
        axis.text.x=element_blank(),
        axis.ticks=element_blank())+ annotate("text", x = 32, y = 2.5, parse = FALSE,label = "weight ~ 1 + s(age) + (1|fcohort) +(1|fyear) + sex",size=5)

# M1 Akaike information criterion summary
summary(AIC(m1))
AIC(m1)


#############################################################################################################################################
#Checking if running m1 a second time produces different or the same outputs.
# m11=m1 (w~s(age)+(1|fcohort)+(1|fyear)+sex)
############################################################################################################################################

#Status quo model (Kelli)
m11 <- sdmTMB::sdmTMB(
  data = weight_at_age_df  |>
    dplyr::filter(sex != "U") |>
    dplyr::mutate(
      age = ifelse(age > 15, 15, age),
      cohort = year - age,
      fyear = as.factor(year),
      fcohort = as.factor(cohort)
    ),
  formula = weight ~ 1 + s(age) + (1|fcohort) + (1|fyear) + sex,
  family = lognormal(link = "log"),
  spatial = "off",
  spatiotemporal = "off",
  time = "year",
  control = sdmTMB::sdmTMBcontrol(newton_loops = 1)
)

print(m11)

# Extracting fixed effect parameters of m11 as a data frame
fixed_eff_m11<-tidy(m11,conf.int=TRUE)


#Extracting Dispersion parameters for random effect of m11 as a data frame
random_par_m11 <- tidy(m11,"ran_pars",conf.int=TRUE)
print(random_par_m11)
# term    estimate std.error conf.low conf.high
# <chr>      <dbl>     <dbl>    <dbl>     <dbl>
# 1 phi       0.261   0.000375   0.260      0.262 # A dispersion parameter for a distribution
# 2 sigma_G   0.199   0.0186     0.166      0.239   # IID random intercept variance (fcohort)
# 3 sigma_G   0.0986  0.0103     0.0803     0.121   # IID random intercept variance (fyear)

#Extracting random effect values by variable of m11
random_eff_val_m11<-tidy(m11,"ran_vals",conf.int=TRUE)
print(random_eff_val_m11,n=110)


# Extracting dispersion parameter of the m1
dispersion_param_m11 <- 0.26
print(dispersion_param_m11)

# Extracting ML criterion at convergence of m11
ML_criterion_m11 <- -154961.173  #this value is coming from the output
print(ML_criterion_m11)


# In summary, we are extracting the random effects from two mixed-effects models (m1 and m11), and merges the results from both models into a single data frame
# for further analysis.
tidy(m11,"ran_vals")
print(tidy(m11,"ran_vals"),n=110)
ran_eff_par_m1<-tidy(m1,"ran_vals")
ran_eff_par_m11<-tidy(m11,"ran_vals")
ran_eff_par_all<-merge(ran_eff_par_m1,ran_eff_par_m11)

#Perform several sanity checks
sanity(m11)

#MODEL DIAGNOSTIC
residuals_m11 <- residuals(m11)
residuals_m11 <- residuals_m11[is.finite(residuals_m11)]
qqnorm(residuals_m11); qqline(residuals_m11)
residuals_plot<-ggplot(predictions_m11, aes(age, weight, col = resids)) + scale_colour_gradient2() +
  geom_point() + facet_wrap(~year)

# create prediction grid
pred_grid11 <- expand.grid(
  year = unique(m11[["data"]][["year"]]),
  sex = unique(m11[["data"]][["sex"]]),
  age = 0:15
) |>
  dplyr::mutate(cohort = year - age) |>
  dplyr::filter(cohort %in% unique(m11[["data"]][["cohort"]])) |>
  dplyr::mutate(
    fyear = as.factor(year),
    fcohort = as.factor(cohort)
  )

# Get estimates
preds11 <- predict(m11, newdata = pred_grid11)
preds11 <- preds11 |>
  dplyr::mutate(est_weight = exp(est))
head(preds11)
preds11_Cohot_plot<-ggplot(preds11, aes(x = fcohort, y = est_weight)) +
  geom_line() +  labs(title = "Cohort effect",
       x = "Cohort",
       y = "Predicted weight") +
  facet_wrap(~ age, scales = "free_x")+theme(axis.text.x = element_text(angle = 90,hjust = 1,size=5))


# make predict EWAA m11
ewaa_m11 <- preds11 |>
  dplyr::group_by(year, age) |>
  dplyr::summarise(
    n = dplyr::n(),
    pred_weight = mean(exp(est)),
    upper=quantile(exp(est),probs = 0.975),
    lower = quantile(exp(est), probs = 0.025),
    run="2"
  ) |>
  as.data.frame() |>
  dplyr::select(year, age, pred_weight,upper,lower,run)
head(ewaa_m11)

# Ggplot visualization of the predicted weight at age over the years, using m11: weight ~ 1 + s(age)+ (1|fyear) + cohort
weight_OP11<-ggplot(ewaa_m11, aes(x = factor(year), y = pred_weight, color = factor(age),
                                  group = factor(age))) +
  geom_text(aes(label=round(age,2)), size = 4.5) +
  geom_line(alpha = 0.85) +ylim(0, 2.5)+
  theme_bw() +
  labs(x = "Year", y = "Weight") + geom_ribbon(data=ewaa_m1,aes(ymin=upper, ymax=lower,fill=factor(age)),alpha=0.1, linetype = "dotted") +
  theme(axis.title = element_text(size = 17),
        axis.text = element_text(size = 15),
        legend.title = element_text(size = 17),
        legend.text = element_text(size = 15),
        axis.text.x = element_text(angle = 90,hjust = 1,size=8),
        legend.position = "none",
        axis.title.x = element_blank(),
        axis.ticks=element_blank())+ annotate("text", x = 32, y = 2.5, parse = FALSE,label = "weight ~ 1 + s(age) + (1|fcohort) +(1|fyear) + sex",size=5)

# summary of the Akaike Information Criterion (AIC) for the fitted model m11
summary(AIC(m11))
AIC(m11)

# Empirical weight-at-age comparison between m1 and m11
ewaa_comp1<-data.frame(ewaa_m1)
ewaa_comp11<-data.frame(ewaa_m11)

# Merge the two data frames based on 'year' and 'age'
merged_data <- rbind(ewaa_comp1, ewaa_comp11)

# Print the merged data with the 'model' column
print(merged_data)
head(merged_data)


# Plotting m1 and m11 to verify if the estimates from GAM-m1 are identical to those from GAM-m11.
# Plot of Age versus predicted weight from 1975 to 2000
plot_1975_2000 <- ggplot(merged_data[merged_data$year <= 2000, ], aes(x = age, y = pred_weight, color = run)) +
  geom_line() +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = run), alpha = 0.3) +
  labs(title = "Running comparison (1975-2000)",
       x = "Age",
       y = "Predicted weight",
       color = "Run",
       fill = "Run") +
  scale_color_manual(values = c("1" = "blue", "2" = "red")) +
  scale_fill_manual(values = c("1" = "blue", "2" = "red")) +
  facet_wrap(~ year, scales = "free_x")


# Filter data for the years 2001-2023 and create plot age versus predicted weight
plot_2001_2023 <- ggplot(merged_data[merged_data$year >= 2001, ], aes(x = age, y = pred_weight, color = run)) +
  geom_line() +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = run), alpha = 0.3) +
  labs(title = "Running comparison (2001-2023)",
       x = "Age",
       y = "Predicted weight",
       color = "Run",
       fill = "Run") +
  scale_color_manual(values = c("1" = "blue", "2" = "red")) +
  scale_fill_manual(values = c("1" = "blue", "2" = "red")) +
  facet_wrap(~ year, scales = "free_x")


Comp_re_full_dat<-ggplot(merged_data, aes(x = age, y = pred_weight, color = run)) +
  geom_line() +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = run), alpha = 0.3) +
  labs(title = "Running comparison by Year",
       x = "Age",
       y = "Predicted weight",
       color = "Run",
       fill = "Run") +
  scale_color_manual(values = c("1" = "blue", "2" = "red")) +
  scale_fill_manual(values = c("1" = "blue", "2" = "red")) +
  facet_wrap(~ year, scales = "free_x")

# Print plots
print(plot_1975_2000)
print(plot_2001_2023)
print(Comp_re_full_dat)


#############################################################################################################################################

# 2 part: Status quo + variability in the projection period of the MSE hake

#####################################################################
# We use this GAM: w~s(age)+(1|fcohort)+(1|fyear)+sex to predict weight-at-age, incorporating random effects for cohort and year.
# Then we extend or add years to the projection period (1975 to 2050)
# Historical period: 1982:2018 Historical years
# Projection period: 2019:2020

#################################################################################################################################################
# Incorporation of variability on the projection period

# Fit the model with a GAM using sdmTMB
m3 <- sdmTMB::sdmTMB(
  data = weight_at_age_df %>%
    dplyr::filter(sex != "U") %>%
    dplyr::mutate(
      age = ifelse(age > 15, 15, age),
      cohort = year - age,
      fyear = as.factor(year),
      fcohort = as.factor(cohort)
    ),
  formula = weight ~ 1 + s(age) + (1|fcohort) + (1|fyear) + sex,
  family = lognormal(link = "log"),
  spatial = "off",
  spatiotemporal = "off",
  time = "year",
  control = sdmTMB::sdmTMBcontrol(newton_loops = 1)
)

# Perform model diagnostics
residuals_m3 <- residuals(m3)
residuals_m3 <- residuals_m3[is.finite(residuals_m3)]
qqnorm(residuals_m3); qqline(residuals_m3)

#Extract coefficients of GAM
tidy(m3)

#Confidence intervals on the fixed effects
tidy(m3,conf.int=TRUE)

#Confidence intervals on the random effects
tidy(m3,"ran_pars",conf.int=TRUE)

tidy(m3,"ran_vals")
print(tidy(m3,"ran_vals"),n=114)
r_eff_val_m3<-data.frame(tidy(m3,"ran_vals")) #random effect values (reffval)

#Perform several sanity checks
sanity(m3)

# Expand the grid for projection period
pred_historic_m3 <- expand.grid(
  year = unique(m3[["data"]][["year"]]),
  sex = unique(m3[["data"]][["sex"]]),
  age = 0:15
) %>%
  dplyr::mutate(cohort = year - age) %>%
  dplyr::filter(cohort %in% unique(m3[["data"]][["cohort"]])) %>%
  dplyr::mutate(
    fyear = as.factor(year),
    fcohort = as.factor(cohort)
  ) %>%
  select(-cohort)  # Remove the cohort column
head(pred_historic_m3)

pred_historic_m3$real_year <- pred_historic_m3$year
unique(pred_historic_m3$real_year)
#expand the grid from 2024 to 2040, this is because the data is until 2023
pred_future_m3 <- expand.grid(real_year = as.character(2024:2040), sex = unique(m3[["data"]][["sex"]]),age=0:15,fyear = unique(m3[["data"]][["fyear"]]),fcohort = unique(m3[["data"]][["fcohort"]]))

#Projection period
year_map_m3 <- data.frame(real_year = as.character(2024:2040), year = sample(unique(m3[["data"]][["year"]]), size=length(2024:2040), replace=T))

#left join year and pred_future_m3 data
pred_future_m3 <- dplyr::left_join(pred_future_m3, year_map_m3)
head.matrix(pred_future_m3)

#Combining the historical and future data frame
pred_all_m3 <- rbind(pred_historic_m3, pred_future_m3)

# Get estimates for all the grid (1982--2040)
preds_m3 <- predict(m3, newdata=pred_all_m3)
preds_m3 <- preds_m3 %>%
  dplyr::mutate(est_weight = est)

head(preds_m3)

# Calculate the predicted weight and statistical quantities on the original logarithmic scale, such as the mean estimate, for the entire grid.
ewaa_m3 <- preds_m3 %>%
  dplyr::group_by(real_year, age) %>%
  dplyr::summarise(
    n = dplyr::n(),
    pred_weight = mean(est), #Calculate predicted weight on original scale (log)
    upper=quantile(est,probs = 0.975),
    lower = quantile(est, probs = 0.025),
    sd=sd(est)
  ) %>%
  as.data.frame() %>%
  dplyr::select(real_year, age, pred_weight,upper,lower,sd)
head(ewaa_m3,n=10)

########################################################################################################
# Creating three weight-at-age matrices (year, cohort, and year+cohort effect)
# Each simulation is set with different seed for reproducibility
# Matrix 1: Add variability into predicted weight matrix based on fyear effects, using standard deviation of year effect (0.0986)
# Matrix 2: Add variability into predicted weight matrix based on fcohort effects, using standard deviation of year effect (0.199)
# Matrix 3. Add variability into predicted weight matrix based on fyear+fcohort effects,
# Rename the predicted weight-at-age matrix 3 to SS3 format as "ewaa.ss"
######################################################################################################
#Sim 0
set.seed(123)

r_eff_val_m3<-data.frame(tidy(m3,"ran_vals")) #  random effects

#Extracting cohort random effects values from the GAM, those values are estimated on the original scale (log).
fyear_est<-r_eff_val_m3[r_eff_val_m3$term %in% paste0("fyear_", 1975:2023), ]
fcohort_est<-r_eff_val_m3[r_eff_val_m3$term %in% paste0("fcohort_", 1961:2022), ]

# Add fyear variability during the projection period by applying the GAM standard deviation of the year effect (with a standard deviation of 0.0986 and a mean of 0).
mean_fyear<-mean(fyear_est$estimate)
sd_fyear<-sd(fyear_est$estimate)

# Step by step building the fcohort matrix using the standard deviations of the past (GAM fcohort random effect values) and the deviations of the future with N(mean=0, sd=0.199)
# Extract the random val (estimate) of cohort from 2008 to 2022
mean_fcohort<-mean(fcohort_est$estimate)
sd_fcohort<-sd(fcohort_est$estimate)

# fyear: extract only fyear and cohort values
fyear_est$estimate <- fyear_est$estimate
fcohort_est$estimate <- fcohort_est$estimate

# Extract only the 'estimate' from 2008 to 2022 of fcohort values of random effect to build diagonal matrix from 2023-2040
fcohort_estimate_values <- fcohort_est$estimate[48:62] # fcohort random values from 2008 to 2022

# Print the extracted estimate values
print(fcohort_estimate_values)

start_year <- 2024
end_year <- 2040
nT <- end_year - start_year + 1 # Number of years in the projection period (2024-2040=17 yrs)
devs_fyear_past <- rep(0, 782)               #  Past standar deviation of fyear, past dim 782 *5 years 1975-2023

devs_fyear_d<-rnorm(nT,mean=0, sd=sd_fyear)  # future standard deviation of fyear
hist(devs_fyear_d)
devs_fyear_fut<-rep(devs_fyear_d, each =16)

#Add fcohort varibility
start_year_23 <- 2023
end_year_23 <- 2040
nT_23 <- end_year_23 - start_year_23 + 1 # Number of years in the projection period (2023-2040=18 yrs)
devs_fcohort_d<- rnorm(nT_23,mean=0, sd=sd_fcohort) # Future devs fcohort
hist(devs_fcohort_d)

# Combine both vectors into a single vector starting with devs_fcohort_d noraml distribuite with mean=0 and sd=0.199
devs_fcohort_past_future <- c(fcohort_estimate_values,devs_fcohort_d )

# Define the years vector
years_cohort <- c(2008:2040)

# Building a diagonal matrix, step by step
# Extract the random effect values from random effect GAM and the sd for 2023 to 2040 normal distribution mean=- and sd=-.199 (this is coming from GAM)
# Create a data frame with years in the and corresponding estimate values
estimate_fcohort_all <- data.frame(year = years_cohort, estimate = devs_fcohort_past_future)

# Number of columns in the matrix (equal to the length of years)
num_cols <- length(estimate_fcohort_all$year)

# Create an empty matrix filled with zeros
result_matrix <- matrix(0, nrow = 17, ncol = num_cols)

# Assign the first row to be the years (converted to integers with reduced decimal places)
result_matrix[1, ] <- format(as.integer(estimate_fcohort_all$year), nsmall = 0)

# Assign the second row to be the estimates
result_matrix[2, ] <- estimate_fcohort_all$estimate

# Print the resulting matrix with integrate sd from the past and future second line
print(result_matrix)

# building the matrix with the fcohort estimates to be repeated diagonally

#i<-result_matrix[2,1] # first value of second row and first column   "-0.201972185011327"
#i<-result_matrix[3,2] # first value third row value and second column "0", I want the -0.201972185011327 in this position
#i<-result_matrix[4,3]  # first value of the fourth row and third column

# Define the dimensions of the result_matrix
num_rows <- nrow(result_matrix)
num_cols <- ncol(result_matrix)

# Define the value to repeat
value_to_repeat <- result_matrix[2,]  # Assuming you want to repeat the string "2" and "1"


# Determine the starting positions for diagonal repetition # 3 # 2
start_row <- 3
start_col <- 2
end_row <- num_rows
end_col <- num_cols

# Fill the specified cells with the repeating pattern [2,1]
for (i in start_row:end_row) {
  for (j in start_col:end_col) {
    if (j >= (i - start_row + start_col)) {
      # Calculate the index within the repeating pattern
      pattern_index <- (j - (i - start_row + start_col)) %% 32 + 1
      result_matrix[i, j] <- value_to_repeat[pattern_index]
    }
  }
}

# Print the updated result_matrix
print(result_matrix)


# creating new columns with fcohort_devs
# Extract 2023 to 2040 columns (numeric values only)
extracted_values <- as.numeric(result_matrix[2:17, 16:33])  # Convert selected part of matrix to numeric

# Remove the first element from the extracted_values vector
extracted_values <- extracted_values[-1]
length(extracted_values)

# Show the resulting vector
print(extracted_values)

#Create a zeros in the past
devs_fcohort_past <- rep(0, 767)    # Past devs fcohort (1975 to 2022 = 48 yrs), dim 767 * 5
# Combine both cohort vectors into a single column
devs_fcohort_all <- c(devs_fcohort_past, extracted_values)

# Combine both fyear vectors into a single column
devs_fyear_all <- c(devs_fyear_past, devs_fyear_fut)

# Add a new column with correct standard deviation of fyear and fcohort
ewaa_m3$fyear_devs <- devs_fyear_all
ewaa_m3$fcohort_devs <- devs_fcohort_all
head(ewaa_m3,n=10)

# Add new column with predicted weight with adjustment to cohort effects, applying an exponential transformation
ewaa_m3$ewaa_fyear <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fyear_devs)
head(ewaa_m3,n=10)
hist(ewaa_m3$ewaa_fyear)

 #Add new column with predicted weight with adjustment to cohort effects, applying an exponential transformation
ewaa_m3$ewaa_fcohort <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fcohort_devs)
head(ewaa_m3,n=10)
print(ewaa_m3)

 #Add new column with predicted weight adjustments or corrections based on year and cohort effects, applying an exponential transformation
ewaa_m3$ewaa_fyear_fcoh <- exp(ewaa_m3$pred_weight)*exp(ewaa_m3$fyear_devs + ewaa_m3$fcohort_devs)
head(ewaa_m3,n=10)
print(ewaa_m3)

# selected the predict weight with fyear and cohort effect only.
ewaa_m3_sim0 <- ewaa_m3[, c("real_year", "age", "ewaa_fyear_fcoh")]
head(ewaa_m3_sim0,n=10)

#rename the columns as SS3 format ewaa.ss
ewaa_m3_sim0 <- ewaa_m3_sim0 %>%
  rename(year = real_year, pred_weight = ewaa_fyear_fcoh)
colnames(ewaa_m3_sim0) <- c("year", "age", "pred_weight")
head(ewaa_m3_sim0,n=10)
unique(ewaa_m3_sim0$year)
str(ewaa_m3_sim0)
# Convert the year column to numeric if necessary
ewaa_m3_sim1$year <- as.numeric(ewaa_m3_sim1$year)

# Run the pad_weight_at_age function and inspect the output
padded_data <- pad_weight_at_age(ewaa_m3_sim1)
str(padded_data)

hakedataUSA:::write_wtatage_file(file="ewaa_sim0.ss",data=pad_weight_at_age(ewaa_m3_sim1),maturity = maturity_at_age)

#########################################################################################
 #Sim 1
set.seed(456)
r_eff_val_m3<-data.frame(tidy(m3,"ran_vals")) #  random effects

#Random effects from sdmTMB are estimated on the original scale.
fyear_est<-r_eff_val_m3[r_eff_val_m3$term %in% paste0("fyear_", 1975:2023), ]
fcohort_est<-r_eff_val_m3[r_eff_val_m3$term %in% paste0("fcohort_", 1961:2022), ]
mean_fyear<-mean(fyear_est$estimate)
sd_fyear<-sd(fyear_est$estimate)

mean_fcohort<-mean(fcohort_est$estimate)
sd_fcohort<-sd(fcohort_est$estimate)

# fyear: extract only fyear and cohort values
fyear_est$estimate <- fyear_est$estimate
fcohort_est$estimate <- fcohort_est$estimate

# Extract only the 'estimate' from 2008 to 2022 of fcohort values of random effect to build diagonal matrix from 2023-2040
fcohort_estimate_values <- fcohort_est$estimate[48:62] # fcohort random values from 2008 to 2022

start_year <- 2024
end_year <- 2040
nT <- end_year - start_year + 1              # Number of years in the projection period (2024-2040=17 yrs)
devs_fyear_past <- rep(0, 782)               #  Past standar deviation of fyear, past dim 782 *5 years 1975-2023
devs_fyear1<-rnorm(nT,mean=0, sd=sd_fyear)   # future standard deviation of fyear
devs_fyear_fut1<-rep(devs_fyear1, each =16)

#Add fcohort varibility
start_year_23 <- 2023
end_year_23 <- 2040
nT_23 <- end_year_23 - start_year_23 + 1 # Number of years in the projection period (2023-2040=18 yrs)
devs_fcohort1<- rnorm(nT_23,mean=0, sd=sd_fcohort) # Future devs fcohort

# Define the vectors
# Combine both vectors into a single vector starting with devs_fcohort_d noraml distribuite with mean=0 and sd=0.199
devs_fcohort_past_future1 <- c(fcohort_estimate_values,devs_fcohort1 )

# Define the years vector
years_cohort <- c(2008:2040)


## Building a diagonal matrix
# Extract the random effect values from random effect GAM and the sd for 2023 to 2040 normal distribution mean=- and sd=-.199 (this is coming from GAM)
# Create a data frame with years in the and corresponding estimate values
estimate_fcohort_all1 <- data.frame(year = years_cohort, estimate = devs_fcohort_past_future1)

# Number of columns in the matrix (equal to the length of years)
num_cols <- length(estimate_fcohort_all1$year)

# Create an empty matrix filled with zeros
result_matrix1 <- matrix(0, nrow = 17, ncol = num_cols)

# Assign the first row to be the years (converted to integers with reduced decimal places)
result_matrix1[1, ] <- format(as.integer(estimate_fcohort_all1$year), nsmall = 0)

# Assign the second row to be the estimates
result_matrix1[2, ] <- estimate_fcohort_all1$estimate

# Print the resulting matrix with integrate sd from the past and future second line
print(result_matrix1)

# Define the dimensions of the result_matrix
num_rows1 <- nrow(result_matrix1)
num_cols1 <- ncol(result_matrix1)

# Define the value to repeat
value_to_repeat1 <- result_matrix1[2,]  # Assuming you want to repeat the string "2" and "1"

# Determine the starting positions for diagonal repetition # 3 # 2
start_row1 <- 3
start_col1 <- 2
end_row1 <- num_rows1
end_col1 <- num_cols1

# Fill the specified cells with the repeating pattern [2,1]
for (i in start_row1:end_row1) {
  for (j in start_col1:end_col1) {
    if (j >= (i - start_row1 + start_col1)) {
      # Calculate the index within the repeating pattern
      pattern_index1 <- (j - (i - start_row1 + start_col1)) %% 32 + 1
      result_matrix1[i, j] <- value_to_repeat1[pattern_index1]
    }
  }
}

# Print the updated result_matrix
print(result_matrix1)


# creating new columns with fcohort_devs
# Extract 2023 to 2040 columns (numeric values only)
extracted_values1 <- as.numeric(result_matrix1[2:17, 16:33])  # Convert selected part of matrix to numeric
# Remove the first element from the extracted_values vector
extracted_values1 <- extracted_values1[-1]

#Create a zeros in the past
devs_fcohort_past <- rep(0, 767)    # Past devs fcohort (1975 to 2022 = 48 yrs) dim 767 * 5


# Combine both cohort vectors into a single column
devs_fcohort_all1 <- c(devs_fcohort_past, extracted_values1)

# Combine both fyear vectors into a single column
devs_fyear_all1 <- c(devs_fyear_past, devs_fyear_fut1)

# Add a new column with correct standard deviation of fyear and fcohort
ewaa_m3$fyear_devs1 <- devs_fyear_all1
ewaa_m3$fcohort_devs1 <- devs_fcohort_all1

# # Add new column with predicted weight with adjustment to fyear effects, applying an exponential transformation
ewaa_m3$ewaa_fyear1 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fyear_devs1)

# Add new column with predicted weight with adjustment to cohort effects, applying an exponential transformation
ewaa_m3$ewaa_fcohort1 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fcohort_devs1)

# Add new column with predicted weight with adjustment to year + cohort effects, applying an exponential transformation
ewaa_m3$ewaa_fyear_fcoh1 <- exp(ewaa_m3$pred_weight)*exp(ewaa_m3$fyear_devs1 + ewaa_m3$fcohort_devs1)
head(ewaa_m3)

# selected the predict weight with fyear and cohort effect only.
ewaa_m3_sim1 <- ewaa_m3[, c("real_year", "age", "ewaa_fyear_fcoh1")]
head(ewaa_m3_sim1,n=10)

#rename the columns as SS3 format ewaa.ss
ewaa_m3_sim1 <- ewaa_m3_sim1 %>%
  rename(year = real_year, pred_weight = ewaa_fyear_fcoh1)
colnames(ewaa_m3_sim1) <- c("year", "age", "pred_weight")
head(ewaa_m3_sim1,n=10)

# Convert the year column to numeric if necessary
ewaa_m3_sim1$year <- as.numeric(ewaa_m3_sim1$year)

# Run the pad_weight_at_age function and inspect the output
padded_data1 <- pad_weight_at_age(ewaa_m3_sim1)
str(padded_data1)

hakedataUSA:::write_wtatage_file(file="ewaa_sim1.ss",data=pad_weight_at_age(ewaa_m3_sim1),maturity = maturity_at_age)


##########################################################################################
 #Sim 2
set.seed(398)
devs_fyear2<-rnorm(nT,mean=0, sd=sd_fyear)  # future standard deviation of fyear
devs_fyear_fut2<-rep(devs_fyear2, each =16)

#Add fcohort varibility
devs_fcohort2<- rnorm(nT_23,mean=0, sd=sd_fcohort) # Future devs fcohort

# Combine both vectors into a single vector starting with devs_fcohort_d noraml distribuite with mean=0 and sd=0.199
devs_fcohort_past_future2 <- c(fcohort_estimate_values,devs_fcohort2)

##building a diagonal matrix
estimate_fcohort_all2 <- data.frame(year = years_cohort, estimate = devs_fcohort_past_future2)

# Create an empty matrix filled with zeros
result_matrix2 <- matrix(0, nrow = 17, ncol = num_cols)

# Assign the first row to be the years (converted to integers with reduced decimal places)
result_matrix2[1, ] <- format(as.integer(estimate_fcohort_all2$year), nsmall = 0)

# Assign the second row to be the estimates
result_matrix2[2, ] <- estimate_fcohort_all2$estimate

# Print the resulting matrix with integrate sd from the past and future second line
print(result_matrix2)

# Define the dimensions of the result_matrix
num_rows2 <- nrow(result_matrix2)
num_cols2 <- ncol(result_matrix2)

# Define the value to repeat
value_to_repeat2 <- result_matrix2[2,]  # Assuming you want to repeat the string "2" and "1"

# Determine the starting positions for diagonal repetition # 3 # 2
start_row2 <- 3
start_col2 <- 2
end_row2 <- num_rows2
end_col2 <- num_cols2

# Fill the specified cells with the repeating pattern [2,1]
for (i in start_row2:end_row2) {
  for (j in start_col2:end_col2) {
    if (j >= (i - start_row2 + start_col2)) {
      # Calculate the index within the repeating pattern
      pattern_index2 <- (j - (i - start_row2 + start_col2)) %% 32 + 1
      result_matrix2[i, j] <- value_to_repeat2[pattern_index2]
    }
  }
}

# Print the updated result_matrix
print(result_matrix2)


# creating new columns with fcohort_devs
# Extract 2023 to 2040 columns (numeric values only)
extracted_values2 <- as.numeric(result_matrix2[2:17, 16:33])  # Convert selected part of matrix to numeric
# Remove the first element from the extracted_values vector
extracted_values2 <- extracted_values2[-1]

# Combine both cohort vectors into a single column
devs_fcohort_all2 <- c(devs_fcohort_past, extracted_values2)

# Combine both fyear vectors into a single column
devs_fyear_all2 <- c(devs_fyear_past, devs_fyear_fut2)

# Add a new column with correct standard deviation of fyear and fcohort
ewaa_m3$fyear_devs2 <- devs_fyear_all2
ewaa_m3$fcohort_devs2 <- devs_fcohort_all2

# Add new column with predicted weight with adjustment to year effects, applying an exponential transformation
ewaa_m3$ewaa_fyear2 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fyear_devs2)

# Add new column with predicted weight with adjustment to cohort effects, applying an exponential transformation
ewaa_m3$ewaa_fcohort2 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fcohort_devs2)

# Add new column with predicted weight with adjustment to year + cohort effects, applying an exponential transformation
ewaa_m3$ewaa_fyear_fcoh2 <- exp(ewaa_m3$pred_weight)*exp(ewaa_m3$fyear_devs2 + ewaa_m3$fcohort_devs2)
head(ewaa_m3,n=800)

# selected the predict weight with fyear and cohort effect only.
ewaa_m3_sim2 <- ewaa_m3[, c("real_year", "age", "ewaa_fyear_fcoh2")]
head(ewaa_m3_sim2,n=800)

#rename the columns as SS3 format ewaa.ss
ewaa_m3_sim2 <- ewaa_m3_sim2 %>%
  rename(year = real_year, pred_weight = ewaa_fyear_fcoh2)
colnames(ewaa_m3_sim2) <- c("year", "age", "pred_weight")
head(ewaa_m3_sim2,n=800)

# Convert the year column to numeric if necessary
ewaa_m3_sim2$year <- as.numeric(ewaa_m3_sim2$year)

# Run the pad_weight_at_age function and inspect the output
padded_data2 <- pad_weight_at_age(ewaa_m3_sim2)
str(padded_data2)

hakedataUSA:::write_wtatage_file(file="ewaa_sim2.ss",data=pad_weight_at_age(ewaa_m3_sim2),maturity = maturity_at_age)

############################################################################################
 #Sim 3

set.seed(698)
devs_fyear3<-rnorm(nT,mean=0, sd=sd_fyear)  # future standard deviation of fyear
devs_fyear_fut3<-rep(devs_fyear3, each =16)

#Add fcohort varibility
devs_fcohort3<- rnorm(nT_23,mean=0, sd=sd_fcohort) # Future devs fcohort

# Combine both vectors into a single vector starting with devs_fcohort_d noraml distribuite with mean=0 and sd=0.199
devs_fcohort_past_future3 <- c(fcohort_estimate_values,devs_fcohort3)

##building a diagonal matrix
estimate_fcohort_all3 <- data.frame(year = years_cohort, estimate = devs_fcohort_past_future3)

# Create an empty matrix filled with zeros
result_matrix3 <- matrix(0, nrow = 17, ncol = num_cols)

# Assign the first row to be the years (converted to integers with reduced decimal places)
result_matrix3[1, ] <- format(as.integer(estimate_fcohort_all3$year), nsmall = 0)

# Assign the second row to be the estimates
result_matrix3[2, ] <- estimate_fcohort_all3$estimate

# Print the resulting matrix with integrate sd from the past and future second line
print(result_matrix3)

# Define the dimensions of the result_matrix
num_rows3 <- nrow(result_matrix3)
num_cols3 <- ncol(result_matrix3)

# Define the value to repeat
value_to_repeat3 <- result_matrix3[2,]  # Assuming you want to repeat the string "2" and "1"

# Determine the starting positions for diagonal repetition # 3 # 2
start_row3 <- 3
start_col3 <- 2
end_row3 <- num_rows3
end_col3 <- num_cols3

# Fill the specified cells with the repeating pattern [2,1]
for (i in start_row3:end_row3) {
  for (j in start_col3:end_col3) {
    if (j >= (i - start_row3 + start_col3)) {
      # Calculate the index within the repeating pattern
      pattern_index3 <- (j - (i - start_row3 + start_col3)) %% 32 + 1
      result_matrix3[i, j] <- value_to_repeat3[pattern_index3]
    }
  }
}

# Print the updated result_matrix
print(result_matrix3)


# creating new columns with fcohort_devs
# Extract 2023 to 2040 columns (numeric values only)
extracted_values3 <- as.numeric(result_matrix3[2:17, 16:33])  # Convert selected part of matrix to numeric
# Remove the first element from the extracted_values vector
extracted_values3 <- extracted_values3[-1]

# Combine both cohort vectors into a single column
devs_fcohort_all3 <- c(devs_fcohort_past, extracted_values3)

# Combine both fyear vectors into a single column
devs_fyear_all3 <- c(devs_fyear_past, devs_fyear_fut3)

# Add a new column with correct standard deviation of fyear and fcohort
ewaa_m3$fyear_devs3 <- devs_fyear_all3
ewaa_m3$fcohort_devs3 <- devs_fcohort_all3

# Add new column with predicted weight with adjustment to year effects, applying an exponential transformation
ewaa_m3$ewaa_fyear3 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fyear_devs3)

# Add new column with predicted weight with adjustment to cohort effects, applying an exponential transformation
ewaa_m3$ewaa_fcohort3 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fcohort_devs3)

# Add new column with predicted weight with adjustment to year + cohort effects, applying an exponential transformation
ewaa_m3$ewaa_fyear_fcoh3 <- exp(ewaa_m3$pred_weight)*exp(ewaa_m3$fyear_devs3 + ewaa_m3$fcohort_devs3)
head(ewaa_m3,n=800)

# selected the predict weight with fyear and cohort effect only.
ewaa_m3_sim3 <- ewaa_m3[, c("real_year", "age", "ewaa_fyear_fcoh3")]
head(ewaa_m3_sim3,n=800)

#rename the columns as SS3 format ewaa.ss
ewaa_m3_sim3 <- ewaa_m3_sim3 %>%
  rename(year = real_year, pred_weight = ewaa_fyear_fcoh3)
colnames(ewaa_m3_sim3) <- c("year", "age", "pred_weight")
head(ewaa_m3_sim3,n=800)

# Convert the year column to numeric if necessary
ewaa_m3_sim3$year <- as.numeric(ewaa_m3_sim3$year)

# Run the pad_weight_at_age function and inspect the output
padded_data3 <- pad_weight_at_age(ewaa_m3_sim3)
str(padded_data3)

hakedataUSA:::write_wtatage_file(file="ewaa_sim3.ss",data=pad_weight_at_age(ewaa_m3_sim3),maturity = maturity_at_age)

############################################################################################
 #sim 4
 set.seed(798)
devs_fyear4<-rnorm(nT,mean=0, sd=sd_fyear)  # future standard deviation of fyear
devs_fyear_fut4<-rep(devs_fyear4, each =16)

#Add fcohort varibility
devs_fcohort4<- rnorm(nT_23,mean=0, sd=sd_fcohort) # Future devs fcohort

# Combine both vectors into a single vector starting with devs_fcohort_d noraml distribuite with mean=0 and sd=0.199
devs_fcohort_past_future4 <- c(fcohort_estimate_values,devs_fcohort4)

##building a diagonal matrix
estimate_fcohort_all4 <- data.frame(year = years_cohort, estimate = devs_fcohort_past_future4)

# Create an empty matrix filled with zeros
result_matrix4 <- matrix(0, nrow = 17, ncol = num_cols)

# Assign the first row to be the years (converted to integers with reduced decimal places)
result_matrix4[1, ] <- format(as.integer(estimate_fcohort_all4$year), nsmall = 0)

# Assign the second row to be the estimates
result_matrix4[2, ] <- estimate_fcohort_all4$estimate

# Print the resulting matrix with integrate sd from the past and future second line
print(result_matrix4)

# Define the dimensions of the result_matrix
num_rows4 <- nrow(result_matrix4)
num_cols4 <- ncol(result_matrix4)

# Define the value to repeat
value_to_repeat4 <- result_matrix4[2,]  # Assuming you want to repeat the string "2" and "1"

# Determine the starting positions for diagonal repetition # 3 # 2
start_row4 <- 3
start_col4 <- 2
end_row4 <- num_rows4
end_col4 <- num_cols4

# Fill the specified cells with the repeating pattern [2,1]
for (i in start_row4:end_row4) {
  for (j in start_col4:end_col4) {
    if (j >= (i - start_row4 + start_col4)) {
      # Calculate the index within the repeating pattern
      pattern_index4 <- (j - (i - start_row4 + start_col4)) %% 32 + 1
      result_matrix4[i, j] <- value_to_repeat4[pattern_index4]
    }
  }
}

# Print the updated result_matrix
print(result_matrix4)


# creating new columns with fcohort_devs
# Extract 2023 to 2040 columns (numeric values only)
extracted_values4 <- as.numeric(result_matrix4[2:17, 16:33])  # Convert selected part of matrix to numeric
# Remove the first element from the extracted_values vector
extracted_values4 <- extracted_values4[-1]

# Combine both cohort vectors into a single column
devs_fcohort_all4 <- c(devs_fcohort_past, extracted_values4)

# Combine both fyear vectors into a single column
devs_fyear_all4 <- c(devs_fyear_past, devs_fyear_fut4)

# Add a new column with correct standard deviation of fyear and fcohort
ewaa_m3$fyear_devs4 <- devs_fyear_all4
ewaa_m3$fcohort_devs4 <- devs_fcohort_all4

# Add new column with predicted weight with adjustment to year effects, applying an exponential transformation
ewaa_m3$ewaa_fyear4 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fyear_devs4)

# Add new column with predicted weight with adjustment to cohort effects, applying an exponential transformation
ewaa_m3$ewaa_fcohort4 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fcohort_devs4)

# Add new column with predicted weight with adjustment to year + cohort effects, applying an exponential transformation
ewaa_m3$ewaa_fyear_fcoh4 <- exp(ewaa_m3$pred_weight)*exp(ewaa_m3$fyear_devs4 + ewaa_m3$fcohort_devs4)
head(ewaa_m3,n=800)

# selected the predict weight with fyear and cohort effect only.
ewaa_m3_sim4 <- ewaa_m3[, c("real_year", "age", "ewaa_fyear_fcoh4")]
head(ewaa_m3_sim4,n=800)

#rename the columns as SS3 format ewaa.ss
ewaa_m3_sim4 <- ewaa_m3_sim4 %>%
  rename(year = real_year, pred_weight = ewaa_fyear_fcoh4)
colnames(ewaa_m3_sim4) <- c("year", "age", "pred_weight")
head(ewaa_m3_sim4,n=800)

# Convert the year column to numeric if necessary
ewaa_m3_sim4$year <- as.numeric(ewaa_m3_sim4$year)

# Run the pad_weight_at_age function and inspect the output
padded_data4 <- pad_weight_at_age(ewaa_m3_sim4)
str(padded_data4)

hakedataUSA:::write_wtatage_file(file="ewaa_sim4.ss",data=pad_weight_at_age(ewaa_m3_sim4),maturity = maturity_at_age)

####################################################################################
 #Sim 5
set.seed(898)
devs_fyear5<-rnorm(nT,mean=0, sd=sd_fyear)  # future standard deviation of fyear
devs_fyear_fut5<-rep(devs_fyear5, each =16)

#Add fcohort varibility
devs_fcohort5<- rnorm(nT_23,mean=0, sd=sd_fcohort) # Future devs fcohort

# Combine both vectors into a single vector starting with devs_fcohort_d noraml distribuite with mean=0 and sd=0.199
devs_fcohort_past_future5 <- c(fcohort_estimate_values,devs_fcohort5)

##building a diagonal matrix
estimate_fcohort_all5 <- data.frame(year = years_cohort, estimate = devs_fcohort_past_future5)

# Create an empty matrix filled with zeros
result_matrix5 <- matrix(0, nrow = 17, ncol = num_cols)

# Assign the first row to be the years (converted to integers with reduced decimal places)
result_matrix5[1, ] <- format(as.integer(estimate_fcohort_all5$year), nsmall = 0)

# Assign the second row to be the estimates
result_matrix5[2, ] <- estimate_fcohort_all5$estimate

# Print the resulting matrix with integrate sd from the past and future second line
print(result_matrix5)

# Define the dimensions of the result_matrix
num_rows5 <- nrow(result_matrix5)
num_cols5 <- ncol(result_matrix5)

# Define the value to repeat
value_to_repeat5 <- result_matrix5[2,]  # Assuming you want to repeat the string "2" and "1"

# Determine the starting positions for diagonal repetition # 3 # 2
start_row5 <- 3
start_col5 <- 2
end_row5 <- num_rows5
end_col5 <- num_cols5

# Fill the specified cells with the repeating pattern [2,1]
for (i in start_row5:end_row5) {
  for (j in start_col5:end_col5) {
    if (j >= (i - start_row5 + start_col5)) {
      # Calculate the index within the repeating pattern
      pattern_index5 <- (j - (i - start_row5 + start_col5)) %% 32 + 1
      result_matrix5[i, j] <- value_to_repeat5[pattern_index5]
    }
  }
}

# Print the updated result_matrix
print(result_matrix5)


# creating new columns with fcohort_devs
# Extract 2023 to 2040 columns (numeric values only)
extracted_values5 <- as.numeric(result_matrix5[2:17, 16:33])  # Convert selected part of matrix to numeric
# Remove the first element from the extracted_values vector
extracted_values5 <- extracted_values5[-1]

# Combine both cohort vectors into a single column
devs_fcohort_all5 <- c(devs_fcohort_past, extracted_values5)

# Combine both fyear vectors into a single column
devs_fyear_all5 <- c(devs_fyear_past, devs_fyear_fut5)

# Add a new column with correct standard deviation of fyear and fcohort
ewaa_m3$fyear_devs5 <- devs_fyear_all5
ewaa_m3$fcohort_devs5 <- devs_fcohort_all5

# Add new column with predicted weight with adjustment to year effects, applying an exponential transformation
ewaa_m3$ewaa_fyear5 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fyear_devs5)

# Add new column with predicted weight with adjustment to cohort effects, applying an exponential transformation
ewaa_m3$ewaa_fcohort5 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fcohort_devs5)

# Add new column with predicted weight with adjustment to year + cohort effects, applying an exponential transformation
ewaa_m3$ewaa_fyear_fcoh5 <- exp(ewaa_m3$pred_weight)*exp(ewaa_m3$fyear_devs5 + ewaa_m3$fcohort_devs5)
head(ewaa_m3,n=800)

# selected the predict weight with fyear and cohort effect only.
ewaa_m3_sim5 <- ewaa_m3[, c("real_year", "age", "ewaa_fyear_fcoh5")]
head(ewaa_m3_sim5,n=800)

#rename the columns as SS3 format ewaa.ss
ewaa_m3_sim5 <- ewaa_m3_sim5 %>%
  rename(year = real_year, pred_weight = ewaa_fyear_fcoh5)
colnames(ewaa_m3_sim5) <- c("year", "age", "pred_weight")
head(ewaa_m3_sim5,n=800)

# Convert the year column to numeric if necessary
ewaa_m3_sim5$year <- as.numeric(ewaa_m3_sim5$year)

# Run the pad_weight_at_age function and inspect the output
padded_data5 <- pad_weight_at_age(ewaa_m3_sim5)
str(padded_data5)

hakedataUSA:::write_wtatage_file(file="ewaa_sim5.ss",data=pad_weight_at_age(ewaa_m3_sim5),maturity = maturity_at_age)

 #############################################################################################
 #Sim 6
set.seed(998)
devs_fyear6<-rnorm(nT,mean=0, sd=sd_fyear)  # future standard deviation of fyear
devs_fyear_fut6<-rep(devs_fyear6, each =16)

#Add fcohort varibility
devs_fcohort6<- rnorm(nT_23,mean=0, sd=sd_fcohort) # Future devs fcohort

# Combine both vectors into a single vector starting with devs_fcohort_d noraml distribuite with mean=0 and sd=0.199
devs_fcohort_past_future6 <- c(fcohort_estimate_values,devs_fcohort6)

##building a diagonal matrix
estimate_fcohort_all6 <- data.frame(year = years_cohort, estimate = devs_fcohort_past_future6)

# Create an empty matrix filled with zeros
result_matrix6 <- matrix(0, nrow = 17, ncol = num_cols)

# Assign the first row to be the years (converted to integers with reduced decimal places)
result_matrix6[1, ] <- format(as.integer(estimate_fcohort_all6$year), nsmall = 0)

# Assign the second row to be the estimates
result_matrix6[2, ] <- estimate_fcohort_all6$estimate

# Print the resulting matrix with integrate sd from the past and future second line
print(result_matrix6)

# Define the dimensions of the result_matrix
num_rows6 <- nrow(result_matrix6)
num_cols6 <- ncol(result_matrix6)

# Define the value to repeat
value_to_repeat6 <- result_matrix6[2,]  # Assuming you want to repeat the string "2" and "1"

# Determine the starting positions for diagonal repetition # 3 # 2
start_row6 <- 3
start_col6 <- 2
end_row6 <- num_rows6
end_col6 <- num_cols6

# Fill the specified cells with the repeating pattern [2,1]
for (i in start_row6:end_row6) {
  for (j in start_col6:end_col6) {
    if (j >= (i - start_row6 + start_col6)) {
      # Calculate the index within the repeating pattern
      pattern_index6 <- (j - (i - start_row6 + start_col6)) %% 32 + 1
      result_matrix6[i, j] <- value_to_repeat6[pattern_index6]
    }
  }
}

# Print the updated result_matrix
print(result_matrix6)

# creating new columns with fcohort_devs
# Extract 2023 to 2040 columns (numeric values only)
extracted_values6 <- as.numeric(result_matrix6[2:17, 16:33])  # Convert selected part of matrix to numeric
# Remove the first element from the extracted_values vector
extracted_values6 <- extracted_values6[-1]

# Combine both cohort vectors into a single column
devs_fcohort_all6 <- c(devs_fcohort_past, extracted_values6)

# Combine both fyear vectors into a single column
devs_fyear_all6 <- c(devs_fyear_past, devs_fyear_fut6)

# Add a new column with correct standard deviation of fyear and fcohort
ewaa_m3$fyear_devs6 <- devs_fyear_all6
ewaa_m3$fcohort_devs6 <- devs_fcohort_all6

# Add new column with predicted weight with adjustment to year effects, applying an exponential transformation
ewaa_m3$ewaa_fyear6 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fyear_devs6)

# Add new column with predicted weight with adjustment to cohort effects, applying an exponential transformation
ewaa_m3$ewaa_fcohort6 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fcohort_devs6)

# Add new column with predicted weight with adjustment to year + cohort effects, applying an exponential transformation
ewaa_m3$ewaa_fyear_fcoh6 <- exp(ewaa_m3$pred_weight)*exp(ewaa_m3$fyear_devs6 + ewaa_m3$fcohort_devs6)
head(ewaa_m3,n=800)

# selected the predict weight with fyear and cohort effect only.
ewaa_m3_sim6 <- ewaa_m3[, c("real_year", "age", "ewaa_fyear_fcoh6")]
head(ewaa_m3_sim6,n=800)

#rename the columns as SS3 format ewaa.ss
ewaa_m3_sim6 <- ewaa_m3_sim6 %>%
  rename(year = real_year, pred_weight = ewaa_fyear_fcoh6)
colnames(ewaa_m3_sim6) <- c("year", "age", "pred_weight")
head(ewaa_m3_sim6,n=800)

# Convert the year column to numeric if necessary
ewaa_m3_sim6$year <- as.numeric(ewaa_m3_sim6$year)

# Run the pad_weight_at_age function and inspect the output
padded_data6 <- pad_weight_at_age(ewaa_m3_sim6)
str(padded_data6)

hakedataUSA:::write_wtatage_file(file="ewaa_sim6.ss",data=pad_weight_at_age(ewaa_m3_sim6),maturity = maturity_at_age)


############################################################################################
 #Sim 7
set.seed(466)
devs_fyear7<-rnorm(nT,mean=0, sd=sd_fyear)  # future standard deviation of fyear
devs_fyear_fut7<-rep(devs_fyear7, each =16)

#Add fcohort varibility
devs_fcohort7<- rnorm(nT_23,mean=0, sd=sd_fcohort) # Future devs fcohort

# Combine both vectors into a single vector starting with devs_fcohort_d noraml distribuite with mean=0 and sd=0.199
devs_fcohort_past_future7 <- c(fcohort_estimate_values,devs_fcohort7)

##building a diagonal matrix
estimate_fcohort_all7 <- data.frame(year = years_cohort, estimate = devs_fcohort_past_future7)

# Create an empty matrix filled with zeros
result_matrix7 <- matrix(0, nrow = 17, ncol = num_cols)

# Assign the first row to be the years (converted to integers with reduced decimal places)
result_matrix7[1, ] <- format(as.integer(estimate_fcohort_all7$year), nsmall = 0)

# Assign the second row to be the estimates
result_matrix7[2, ] <- estimate_fcohort_all7$estimate

# Print the resulting matrix with integrate sd from the past and future second line
print(result_matrix7)

# Define the dimensions of the result_matrix
num_rows7 <- nrow(result_matrix7)
num_cols7 <- ncol(result_matrix7)

# Define the value to repeat
value_to_repeat7 <- result_matrix7[2,]  # Assuming you want to repeat the string "2" and "1"

# Determine the starting positions for diagonal repetition # 3 # 2
start_row7 <- 3
start_col7 <- 2
end_row7 <- num_rows7
end_col7 <- num_cols7

# Fill the specified cells with the repeating pattern [2,1]
for (i in start_row7:end_row7) {
  for (j in start_col7:end_col7) {
    if (j >= (i - start_row7 + start_col7)) {
      # Calculate the index within the repeating pattern
      pattern_index7 <- (j - (i - start_row7 + start_col7)) %% 32 + 1
      result_matrix7[i, j] <- value_to_repeat7[pattern_index7]
    }
  }
}

# Print the updated result_matrix
print(result_matrix7)

# creating new columns with fcohort_devs
# Extract 2023 to 2040 columns (numeric values only)
extracted_values7 <- as.numeric(result_matrix7[2:17, 16:33])  # Convert selected part of matrix to numeric
# Remove the first element from the extracted_values vector
extracted_values7 <- extracted_values7[-1]

# Combine both cohort vectors into a single column
devs_fcohort_all7 <- c(devs_fcohort_past, extracted_values7)

# Combine both fyear vectors into a single column
devs_fyear_all7 <- c(devs_fyear_past, devs_fyear_fut7)

# Add a new column with correct standard deviation of fyear and fcohort
ewaa_m3$fyear_devs7 <- devs_fyear_all7
ewaa_m3$fcohort_devs7 <- devs_fcohort_all7

# Add new column with predicted weight with adjustment to year effects, applying an exponential transformation
ewaa_m3$ewaa_fyear7 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fyear_devs7)

# Add new column with predicted weight with adjustment to cohort effects, applying an exponential transformation
ewaa_m3$ewaa_fcohort7 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fcohort_devs7)

# Add new column with predicted weight with adjustment to year + cohort effects, applying an exponential transformation
ewaa_m3$ewaa_fyear_fcoh7 <- exp(ewaa_m3$pred_weight)*exp(ewaa_m3$fyear_devs7 + ewaa_m3$fcohort_devs7)
head(ewaa_m3,n=800)

# selected the predict weight with fyear and cohort effect only.
ewaa_m3_sim7 <- ewaa_m3[, c("real_year", "age", "ewaa_fyear_fcoh7")]
head(ewaa_m3_sim7,n=800)

#rename the columns as SS3 format ewaa.ss
ewaa_m3_sim7 <- ewaa_m3_sim7 %>%
  rename(year = real_year, pred_weight = ewaa_fyear_fcoh7)
colnames(ewaa_m3_sim7) <- c("year", "age", "pred_weight")
head(ewaa_m3_sim7,n=800)

# Convert the year column to numeric if necessary
ewaa_m3_sim7$year <- as.numeric(ewaa_m3_sim7$year)

# Run the pad_weight_at_age function and inspect the output
padded_data7 <- pad_weight_at_age(ewaa_m3_sim7)
str(padded_data7)

hakedataUSA:::write_wtatage_file(file="ewaa_sim7.ss",data=pad_weight_at_age(ewaa_m3_sim7),maturity = maturity_at_age)

#############################################################################################
# Sim 8
set.seed(355)
devs_fyear8<-rnorm(nT,mean=0, sd=sd_fyear)  # future standard deviation of fyear
devs_fyear_fut8<-rep(devs_fyear8, each =16)

#Add fcohort varibility
devs_fcohort8<- rnorm(nT_23,mean=0, sd=sd_fcohort) # Future devs fcohort

# Combine both vectors into a single vector starting with devs_fcohort_d noraml distribuite with mean=0 and sd=0.199
devs_fcohort_past_future8 <- c(fcohort_estimate_values,devs_fcohort8)

##building a diagonal matrix
estimate_fcohort_all8 <- data.frame(year = years_cohort, estimate = devs_fcohort_past_future8)

# Create an empty matrix filled with zeros
result_matrix8 <- matrix(0, nrow = 17, ncol = num_cols)

# Assign the first row to be the years (converted to integers with reduced decimal places)
result_matrix8[1, ] <- format(as.integer(estimate_fcohort_all8$year), nsmall = 0)

# Assign the second row to be the estimates
result_matrix8[2, ] <- estimate_fcohort_all8$estimate

# Print the resulting matrix with integrate sd from the past and future second line
print(result_matrix8)

# Define the dimensions of the result_matrix
num_rows8 <- nrow(result_matrix8)
num_cols8 <- ncol(result_matrix8)

# Define the value to repeat
value_to_repeat8 <- result_matrix8[2,]  # Assuming you want to repeat the string "2" and "1"

# Determine the starting positions for diagonal repetition # 3 # 2
start_row8 <- 3
start_col8 <- 2
end_row8 <- num_rows8
end_col8 <- num_cols8

# Fill the specified cells with the repeating pattern [2,1]
for (i in start_row8:end_row8) {
  for (j in start_col8:end_col8) {
    if (j >= (i - start_row8 + start_col8)) {
      # Calculate the index within the repeating pattern
      pattern_index8 <- (j - (i - start_row8 + start_col8)) %% 32 + 1
      result_matrix8[i, j] <- value_to_repeat8[pattern_index8]
    }
  }
}

# Print the updated result_matrix
print(result_matrix8)

# creating new columns with fcohort_devs
# Extract 2023 to 2040 columns (numeric values only)
extracted_values8 <- as.numeric(result_matrix8[2:17, 16:33])  # Convert selected part of matrix to numeric
# Remove the first element from the extracted_values vector
extracted_values8 <- extracted_values8[-1]

# Combine both cohort vectors into a single column
devs_fcohort_all8 <- c(devs_fcohort_past, extracted_values8)

# Combine both fyear vectors into a single column
devs_fyear_all8 <- c(devs_fyear_past, devs_fyear_fut8)

# Add a new column with correct standard deviation of fyear and fcohort
ewaa_m3$fyear_devs8 <- devs_fyear_all8
ewaa_m3$fcohort_devs8 <- devs_fcohort_all8

# Add new column with predicted weight with adjustment to year effects, applying an exponential transformation
ewaa_m3$ewaa_fyear8 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fyear_devs8)

# Add new column with predicted weight with adjustment to cohort effects, applying an exponential transformation
ewaa_m3$ewaa_fcohort8 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fcohort_devs8)

# Add new column with predicted weight with adjustment to year + cohort effects, applying an exponential transformation
ewaa_m3$ewaa_fyear_fcoh8 <- exp(ewaa_m3$pred_weight)*exp(ewaa_m3$fyear_devs8 + ewaa_m3$fcohort_devs8)
head(ewaa_m3,n=800)

# selected the predict weight with fyear and cohort effect only.
ewaa_m3_sim8 <- ewaa_m3[, c("real_year", "age", "ewaa_fyear_fcoh8")]
head(ewaa_m3_sim8,n=800)

#rename the columns as SS3 format ewaa.ss
ewaa_m3_sim8 <- ewaa_m3_sim8 %>%
  rename(year = real_year, pred_weight = ewaa_fyear_fcoh8)
colnames(ewaa_m3_sim8) <- c("year", "age", "pred_weight")
head(ewaa_m3_sim8,n=800)

# Convert the year column to numeric if necessary
ewaa_m3_sim8$year <- as.numeric(ewaa_m3_sim8$year)

# Run the pad_weight_at_age function and inspect the output
padded_data8 <- pad_weight_at_age(ewaa_m3_sim8)
str(padded_data8)

hakedataUSA:::write_wtatage_file(file="ewaa_sim8.ss",data=pad_weight_at_age(ewaa_m3_sim8),maturity = maturity_at_age)

##############################################################################################
 #Sim 9
set.seed(811)
devs_fyear9<-rnorm(nT,mean=0, sd=sd_fyear)  # future standard deviation of fyear
devs_fyear_fut9<-rep(devs_fyear9, each =16)

#Add fcohort varibility
devs_fcohort9<- rnorm(nT_23,mean=0, sd=sd_fcohort) # Future devs fcohort

# Combine both vectors into a single vector starting with devs_fcohort_d noraml distribuite with mean=0 and sd=0.199
devs_fcohort_past_future9 <- c(fcohort_estimate_values,devs_fcohort9)

##building a diagonal matrix
estimate_fcohort_all9 <- data.frame(year = years_cohort, estimate = devs_fcohort_past_future9)

# Create an empty matrix filled with zeros
result_matrix9 <- matrix(0, nrow = 17, ncol = num_cols)

# Assign the first row to be the years (converted to integers with reduced decimal places)
result_matrix9[1, ] <- format(as.integer(estimate_fcohort_all9$year), nsmall = 0)

# Assign the second row to be the estimates
result_matrix9[2, ] <- estimate_fcohort_all9$estimate

# Print the resulting matrix with integrate sd from the past and future second line
print(result_matrix9)

# Define the dimensions of the result_matrix
num_rows9 <- nrow(result_matrix9)
num_cols9 <- ncol(result_matrix9)

# Define the value to repeat
value_to_repeat9 <- result_matrix9[2,]  # Assuming you want to repeat the string "2" and "1"

# Determine the starting positions for diagonal repetition # 3 # 2
start_row9 <- 3
start_col9 <- 2
end_row9 <- num_rows9
end_col9 <- num_cols9

# Fill the specified cells with the repeating pattern [2,1]
for (i in start_row9:end_row9) {
  for (j in start_col9:end_col9) {
    if (j >= (i - start_row9 + start_col9)) {
      # Calculate the index within the repeating pattern
      pattern_index9 <- (j - (i - start_row9 + start_col9)) %% 32 + 1
      result_matrix9[i, j] <- value_to_repeat9[pattern_index9]
    }
  }
}

# Print the updated result_matrix
print(result_matrix9)

# creating new columns with fcohort_devs
# Extract 2023 to 2040 columns (numeric values only)
extracted_values9 <- as.numeric(result_matrix9[2:17, 16:33])  # Convert selected part of matrix to numeric
# Remove the first element from the extracted_values vector
extracted_values9 <- extracted_values9[-1]

# Combine both cohort vectors into a single column
devs_fcohort_all9 <- c(devs_fcohort_past, extracted_values9)

# Combine both fyear vectors into a single column
devs_fyear_all9 <- c(devs_fyear_past, devs_fyear_fut9)

# Add a new column with correct standard deviation of fyear and fcohort
ewaa_m3$fyear_devs9 <- devs_fyear_all9
ewaa_m3$fcohort_devs9 <- devs_fcohort_all9

# Add new column with predicted weight with adjustment to year effects, applying an exponential transformation
ewaa_m3$ewaa_fyear9 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fyear_devs9)

# Add new column with predicted weight with adjustment to cohort effects, applying an exponential transformation
ewaa_m3$ewaa_fcohort9 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fcohort_devs9)

# Add new column with predicted weight with adjustment to year + cohort effects, applying an exponential transformation
ewaa_m3$ewaa_fyear_fcoh9 <- exp(ewaa_m3$pred_weight)*exp(ewaa_m3$fyear_devs9 + ewaa_m3$fcohort_devs9)
head(ewaa_m3,n=800)

# selected the predict weight with fyear and cohort effect only.
ewaa_m3_sim9 <- ewaa_m3[, c("real_year", "age", "ewaa_fyear_fcoh9")]
head(ewaa_m3_sim9,n=800)

#rename the columns as SS3 format ewaa.ss
ewaa_m3_sim9 <- ewaa_m3_sim9 %>%
  rename(year = real_year, pred_weight = ewaa_fyear_fcoh9)
colnames(ewaa_m3_sim9) <- c("year", "age", "pred_weight")
head(ewaa_m3_sim9,n=800)

# Convert the year column to numeric if necessary
ewaa_m3_sim9$year <- as.numeric(ewaa_m3_sim9$year)

# Run the pad_weight_at_age function and inspect the output
padded_data9 <- pad_weight_at_age(ewaa_m3_sim9)
str(padded_data9)

hakedataUSA:::write_wtatage_file(file="ewaa_sim9.ss",data=pad_weight_at_age(ewaa_m3_sim9),maturity = maturity_at_age)

##########################################################################################
#Sim 10
set.seed(441)
devs_fyear10<-rnorm(nT,mean=0, sd=sd_fyear)  # future standard deviation of fyear
devs_fyear_fut10<-rep(devs_fyear10, each =16)

#Add fcohort varibility
devs_fcohort10<- rnorm(nT_23,mean=0, sd=sd_fcohort) # Future devs fcohort

# Combine both vectors into a single vector starting with devs_fcohort_d noraml distribuite with mean=0 and sd=0.199
devs_fcohort_past_future10 <- c(fcohort_estimate_values,devs_fcohort10)

##building a diagonal matrix
estimate_fcohort_all10 <- data.frame(year = years_cohort, estimate = devs_fcohort_past_future10)

# Create an empty matrix filled with zeros
result_matrix10 <- matrix(0, nrow = 17, ncol = num_cols)

# Assign the first row to be the years (converted to integers with reduced decimal places)
result_matrix10[1, ] <- format(as.integer(estimate_fcohort_all10$year), nsmall = 0)

# Assign the second row to be the estimates
result_matrix10[2, ] <- estimate_fcohort_all10$estimate

# Print the resulting matrix with integrate sd from the past and future second line
print(result_matrix10)

# Define the dimensions of the result_matrix
num_rows10 <- nrow(result_matrix10)
num_cols10 <- ncol(result_matrix10)

# Define the value to repeat
value_to_repeat10 <- result_matrix10[2,]  # Assuming you want to repeat the string "2" and "1"

# Determine the starting positions for diagonal repetition # 3 # 2
start_row10 <- 3
start_col10 <- 2
end_row10 <- num_rows10
end_col10 <- num_cols10

# Fill the specified cells with the repeating pattern [2,1]
for (i in start_row10:end_row10) {
  for (j in start_col10:end_col10) {
    if (j >= (i - start_row10 + start_col10)) {
      # Calculate the index within the repeating pattern
      pattern_index10 <- (j - (i - start_row10 + start_col10)) %% 32 + 1
      result_matrix10[i, j] <- value_to_repeat10[pattern_index10]
    }
  }
}

# Print the updated result_matrix
print(result_matrix10)

# creating new columns with fcohort_devs
# Extract 2023 to 2040 columns (numeric values only)
extracted_values10 <- as.numeric(result_matrix10[2:17, 16:33])  # Convert selected part of matrix to numeric
# Remove the first element from the extracted_values vector
extracted_values10 <- extracted_values10[-1]

# Combine both cohort vectors into a single column
devs_fcohort_all10 <- c(devs_fcohort_past, extracted_values10)

# Combine both fyear vectors into a single column
devs_fyear_all10 <- c(devs_fyear_past, devs_fyear_fut10)

# Add a new column with correct standard deviation of fyear and fcohort
ewaa_m3$fyear_devs10 <- devs_fyear_all10
ewaa_m3$fcohort_devs10 <- devs_fcohort_all10

# Add new column with predicted weight with adjustment to year effects, applying an exponential transformation
ewaa_m3$ewaa_fyear10 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fyear_devs10)

# Add new column with predicted weight with adjustment to cohort effects, applying an exponential transformation
ewaa_m3$ewaa_fcohort10 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fcohort_devs10)

# Add new column with predicted weight with adjustment to year cohort effects, applying an exponential transformation
ewaa_m3$ewaa_fyear_fcoh10 <- exp(ewaa_m3$pred_weight)*exp(ewaa_m3$fyear_devs10 + ewaa_m3$fcohort_devs10)
head(ewaa_m3,n=800)

# selected the predict weight with fyear and cohort effect only.
ewaa_m3_sim10 <- ewaa_m3[, c("real_year", "age", "ewaa_fyear_fcoh10")]
head(ewaa_m3_sim10,n=800)

#rename the columns as SS3 format ewaa.ss
ewaa_m3_sim10 <- ewaa_m3_sim10 %>%
  rename(year = real_year, pred_weight = ewaa_fyear_fcoh10)
colnames(ewaa_m3_sim10) <- c("year", "age", "pred_weight")
head(ewaa_m3_sim10,n=800)

# Convert the year column to numeric if necessary
ewaa_m3_sim10$year <- as.numeric(ewaa_m3_sim10$year)

# Run the pad_weight_at_age function and inspect the output
padded_data10 <- pad_weight_at_age(ewaa_m3_sim10)
str(padded_data10)

hakedataUSA:::write_wtatage_file(file="ewaa_sim10.ss",data=pad_weight_at_age(ewaa_m3_sim10),maturity = maturity_at_age)


########################################################################################
# Sim 11
set.seed(355)
devs_fyear11 <- rnorm(nT, mean=0, sd=sd_fyear)
devs_fyear_fut11 <- rep(devs_fyear11, each=16)
devs_fcohort11 <- rnorm(nT_23, mean=0, sd=sd_fcohort)
devs_fcohort_past_future11 <- c(fcohort_estimate_values, devs_fcohort11)
estimate_fcohort_all11 <- data.frame(year = years_cohort, estimate = devs_fcohort_past_future11)
result_matrix11 <- matrix(0, nrow = 17, ncol = num_cols)
result_matrix11[1, ] <- format(as.integer(estimate_fcohort_all11$year), nsmall = 0)
result_matrix11[2, ] <- estimate_fcohort_all11$estimate
print(result_matrix11)
num_rows11 <- nrow(result_matrix11)
num_cols11 <- ncol(result_matrix11)
value_to_repeat11 <- result_matrix11[2,]
start_row11 <- 3
start_col11 <- 2
end_row11 <- num_rows11
end_col11 <- num_cols11
for (i in start_row11:end_row11) {
  for (j in start_col11:end_col11) {
    if (j >= (i - start_row11 + start_col11)) {
      pattern_index11 <- (j - (i - start_row11 + start_col11)) %% 32 + 1
      result_matrix11[i, j] <- value_to_repeat11[pattern_index11]
    }
  }
}
print(result_matrix11)
extracted_values11 <- as.numeric(result_matrix11[2:17, 16:33])
extracted_values11 <- extracted_values11[-1]
devs_fcohort_all11 <- c(devs_fcohort_past, extracted_values11)
devs_fyear_all11 <- c(devs_fyear_past, devs_fyear_fut11)
ewaa_m3$fyear_devs11 <- devs_fyear_all11
ewaa_m3$fcohort_devs11 <- devs_fcohort_all11
# Add new column with predicted weight with adjustment to year effects, applying an exponential transformation
ewaa_m3$ewaa_fyear11 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fyear_devs11)
# Add new column with predicted weight with adjustment to cohort effects, applying an exponential transformation
ewaa_m3$ewaa_fcohort11 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fcohort_devs11)
# Add new column with predicted weight with adjustment to year+ cohort effects, applying an exponential transformation
ewaa_m3$ewaa_fyear_fcoh11 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fyear_devs11 + ewaa_m3$fcohort_devs11)
head(ewaa_m3, n=800)

# selected the predict weight with fyear and cohort effect only.
ewaa_m3_sim11 <- ewaa_m3[, c("real_year", "age", "ewaa_fyear_fcoh11")]
head(ewaa_m3_sim11,n=800)

#rename the columns as SS3 format ewaa.ss
ewaa_m3_sim11 <- ewaa_m3_sim11 %>%
  rename(year = real_year, pred_weight = ewaa_fyear_fcoh11)
colnames(ewaa_m3_sim11) <- c("year", "age", "pred_weight")
head(ewaa_m3_sim11,n=800)

# Convert the year column to numeric if necessary
ewaa_m3_sim11$year <- as.numeric(ewaa_m3_sim11$year)

# Run the pad_weight_at_age function and inspect the output
padded_data11 <- pad_weight_at_age(ewaa_m3_sim11)
str(padded_data11)

hakedataUSA:::write_wtatage_file(file="ewaa_sim11.ss",data=pad_weight_at_age(ewaa_m3_sim11),maturity = maturity_at_age)


####################################################################################
# Sim 12
set.seed(661)
devs_fyear12 <- rnorm(nT, mean=0, sd=sd_fyear)
devs_fyear_fut12 <- rep(devs_fyear12, each=16)
devs_fcohort12 <- rnorm(nT_23, mean=0, sd=sd_fcohort)
devs_fcohort_past_future12 <- c(fcohort_estimate_values, devs_fcohort12)
estimate_fcohort_all12 <- data.frame(year = years_cohort, estimate = devs_fcohort_past_future12)
result_matrix12 <- matrix(0, nrow = 17, ncol = num_cols)
result_matrix12[1, ] <- format(as.integer(estimate_fcohort_all12$year), nsmall = 0)
result_matrix12[2, ] <- estimate_fcohort_all12$estimate
print(result_matrix12)
num_rows12 <- nrow(result_matrix12)
num_cols12 <- ncol(result_matrix12)
value_to_repeat12 <- result_matrix12[2,]
start_row12 <- 3
start_col12 <- 2
end_row12 <- num_rows12
end_col12 <- num_cols12
for (i in start_row12:end_row12) {
  for (j in start_col12:end_col12) {
    if (j >= (i - start_row12 + start_col12)) {
      pattern_index12 <- (j - (i - start_row12 + start_col12)) %% 32 + 1
      result_matrix12[i, j] <- value_to_repeat12[pattern_index12]
    }
  }
}
print(result_matrix12)
extracted_values12 <- as.numeric(result_matrix12[2:17, 16:33])
extracted_values12 <- extracted_values12[-1]
devs_fcohort_all12 <- c(devs_fcohort_past, extracted_values12)
devs_fyear_all12 <- c(devs_fyear_past, devs_fyear_fut12)
ewaa_m3$fyear_devs12 <- devs_fyear_all12
ewaa_m3$fcohort_devs12 <- devs_fcohort_all12
ewaa_m3$ewaa_fyear12 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fyear_devs12)
ewaa_m3$ewaa_fcohort12 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fcohort_devs12)
ewaa_m3$ewaa_fyear_fcoh12 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fyear_devs12 + ewaa_m3$fcohort_devs12)
head(ewaa_m3, n=800)

# selected the predict weight with fyear and cohort effect only.
ewaa_m3_sim12 <- ewaa_m3[, c("real_year", "age", "ewaa_fyear_fcoh12")]
head(ewaa_m3_sim12,n=800)

#rename the columns as SS3 format ewaa.ss
ewaa_m3_sim12 <- ewaa_m3_sim12 %>%
  rename(year = real_year, pred_weight = ewaa_fyear_fcoh12)
colnames(ewaa_m3_sim12) <- c("year", "age", "pred_weight")
head(ewaa_m3_sim12,n=800)

# Convert the year column to numeric if necessary
ewaa_m3_sim12$year <- as.numeric(ewaa_m3_sim12$year)

# Run the pad_weight_at_age function and inspect the output
padded_data12 <- pad_weight_at_age(ewaa_m3_sim12)
str(padded_data12)

hakedataUSA:::write_wtatage_file(file="ewaa_sim12.ss",data=pad_weight_at_age(ewaa_m3_sim12),maturity = maturity_at_age)

##############################################################################
# Sim 13
set.seed(871)
devs_fyear13 <- rnorm(nT, mean=0, sd=sd_fyear)
devs_fyear_fut13 <- rep(devs_fyear13, each=16)
devs_fcohort13 <- rnorm(nT_23, mean=0, sd=sd_fcohort)
devs_fcohort_past_future13 <- c(fcohort_estimate_values, devs_fcohort13)
estimate_fcohort_all13 <- data.frame(year = years_cohort, estimate = devs_fcohort_past_future13)
result_matrix13 <- matrix(0, nrow = 17, ncol = num_cols)
result_matrix13[1, ] <- format(as.integer(estimate_fcohort_all13$year), nsmall = 0)
result_matrix13[2, ] <- estimate_fcohort_all13$estimate
print(result_matrix13)
num_rows13 <- nrow(result_matrix13)
num_cols13 <- ncol(result_matrix13)
value_to_repeat13 <- result_matrix13[2,]
start_row13 <- 3
start_col13 <- 2
end_row13 <- num_rows13
end_col13 <- num_cols13
for (i in start_row13:end_row13) {
  for (j in start_col13:end_col13) {
    if (j >= (i - start_row13 + start_col13)) {
      pattern_index13 <- (j - (i - start_row13 + start_col13)) %% 32 + 1
      result_matrix13[i, j] <- value_to_repeat13[pattern_index13]
    }
  }
}
print(result_matrix13)
extracted_values13 <- as.numeric(result_matrix13[2:17, 16:33])
extracted_values13 <- extracted_values13[-1]
devs_fcohort_all13 <- c(devs_fcohort_past, extracted_values13)
devs_fyear_all13 <- c(devs_fyear_past, devs_fyear_fut13)
ewaa_m3$fyear_devs13 <- devs_fyear_all13
ewaa_m3$fcohort_devs13 <- devs_fcohort_all13
ewaa_m3$ewaa_fyear13 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fyear_devs13)
ewaa_m3$ewaa_fcohort13 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fcohort_devs13)
ewaa_m3$ewaa_fyear_fcoh13 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fyear_devs13 + ewaa_m3$fcohort_devs13)
head(ewaa_m3, n=800)

# selected the predict weight with fyear and cohort effect only.
ewaa_m3_sim13 <- ewaa_m3[, c("real_year", "age", "ewaa_fyear_fcoh13")]
head(ewaa_m3_sim13,n=800)

#rename the columns as SS3 format ewaa.ss
ewaa_m3_sim13 <- ewaa_m3_sim13 %>%
  rename(year = real_year, pred_weight = ewaa_fyear_fcoh13)
colnames(ewaa_m3_sim13) <- c("year", "age", "pred_weight")
head(ewaa_m3_sim13,n=800)

# Convert the year column to numeric if necessary
ewaa_m3_sim13$year <- as.numeric(ewaa_m3_sim13$year)

# Run the pad_weight_at_age function and inspect the output
padded_data13 <- pad_weight_at_age(ewaa_m3_sim13)
str(padded_data13)

hakedataUSA:::write_wtatage_file(file="ewaa_sim13.ss",data=pad_weight_at_age(ewaa_m3_sim13),maturity = maturity_at_age)

###################################################################################
# Sim 14
set.seed(981)
devs_fyear14 <- rnorm(nT, mean=0, sd=sd_fyear)
devs_fyear_fut14 <- rep(devs_fyear14, each=16)
devs_fcohort14 <- rnorm(nT_23, mean=0, sd=sd_fcohort)
devs_fcohort_past_future14 <- c(fcohort_estimate_values, devs_fcohort14)
estimate_fcohort_all14 <- data.frame(year = years_cohort, estimate = devs_fcohort_past_future14)
result_matrix14 <- matrix(0, nrow = 17, ncol = num_cols)
result_matrix14[1, ] <- format(as.integer(estimate_fcohort_all14$year), nsmall = 0)
result_matrix14[2, ] <- estimate_fcohort_all14$estimate
print(result_matrix14)
num_rows14 <- nrow(result_matrix14)
num_cols14 <- ncol(result_matrix14)
value_to_repeat14 <- result_matrix14[2,]
start_row14 <- 3
start_col14 <- 2
end_row14 <- num_rows14
end_col14 <- num_cols14
for (i in start_row14:end_row14) {
  for (j in start_col14:end_col14) {
    if (j >= (i - start_row14 + start_col14)) {
      pattern_index14 <- (j - (i - start_row14 + start_col14)) %% 32 + 1
      result_matrix14[i, j] <- value_to_repeat14[pattern_index14]
    }
  }
}
print(result_matrix14)
extracted_values14 <- as.numeric(result_matrix14[2:17, 16:33])
extracted_values14 <- extracted_values14[-1]
devs_fcohort_all14 <- c(devs_fcohort_past, extracted_values14)
devs_fyear_all14 <- c(devs_fyear_past, devs_fyear_fut14)
ewaa_m3$fyear_devs14 <- devs_fyear_all14
ewaa_m3$fcohort_devs14 <- devs_fcohort_all14
ewaa_m3$ewaa_fyear14 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fyear_devs14)
ewaa_m3$ewaa_fcohort14 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fcohort_devs14)
ewaa_m3$ewaa_fyear_fcoh14 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fyear_devs14 + ewaa_m3$fcohort_devs14)
head(ewaa_m3, n=800)

# selected the predict weight with fyear and cohort effect only.
ewaa_m3_sim14 <- ewaa_m3[, c("real_year", "age", "ewaa_fyear_fcoh14")]
head(ewaa_m3_sim14,n=800)

#rename the columns as SS3 format ewaa.ss
ewaa_m3_sim14 <- ewaa_m3_sim14 %>%
  rename(year = real_year, pred_weight = ewaa_fyear_fcoh14)
colnames(ewaa_m3_sim14) <- c("year", "age", "pred_weight")
head(ewaa_m3_sim14,n=800)

# Convert the year column to numeric if necessary
ewaa_m3_sim14$year <- as.numeric(ewaa_m3_sim14$year)

# Run the pad_weight_at_age function and inspect the output
padded_data14 <- pad_weight_at_age(ewaa_m3_sim14)
str(padded_data14)

hakedataUSA:::write_wtatage_file(file="ewaa_sim14.ss",data=pad_weight_at_age(ewaa_m3_sim14),maturity = maturity_at_age)

##############################################################################
# Sim 15
set.seed(1091)
devs_fyear15 <- rnorm(nT, mean=0, sd=sd_fyear)
devs_fyear_fut15 <- rep(devs_fyear15, each=16)
devs_fcohort15 <- rnorm(nT_23, mean=0, sd=sd_fcohort)
devs_fcohort_past_future15 <- c(fcohort_estimate_values, devs_fcohort15)
estimate_fcohort_all15 <- data.frame(year = years_cohort, estimate = devs_fcohort_past_future15)
result_matrix15 <- matrix(0, nrow = 17, ncol = num_cols)
result_matrix15[1, ] <- format(as.integer(estimate_fcohort_all15$year), nsmall = 0)
result_matrix15[2, ] <- estimate_fcohort_all15$estimate
print(result_matrix15)
num_rows15 <- nrow(result_matrix15)
num_cols15 <- ncol(result_matrix15)
value_to_repeat15 <- result_matrix15[2,]
start_row15 <- 3
start_col15 <- 2
end_row15 <- num_rows15
end_col15 <- num_cols15
for (i in start_row15:end_row15) {
  for (j in start_col15:end_col15) {
    if (j >= (i - start_row15 + start_col15)) {
      pattern_index15 <- (j - (i - start_row15 + start_col15)) %% 32 + 1
      result_matrix15[i, j] <- value_to_repeat15[pattern_index15]
    }
  }
}
print(result_matrix15)
extracted_values15 <- as.numeric(result_matrix15[2:17, 16:33])
extracted_values15 <- extracted_values15[-1]
devs_fcohort_all15 <- c(devs_fcohort_past, extracted_values15)
devs_fyear_all15 <- c(devs_fyear_past, devs_fyear_fut15)
ewaa_m3$fyear_devs15 <- devs_fyear_all15
ewaa_m3$fcohort_devs15 <- devs_fcohort_all15
ewaa_m3$ewaa_fyear15 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fyear_devs15)
ewaa_m3$ewaa_fcohort15 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fcohort_devs15)
ewaa_m3$ewaa_fyear_fcoh15 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fyear_devs15 + ewaa_m3$fcohort_devs15)
head(ewaa_m3, n=800)

# selected the predict weight with fyear and cohort effect only.
ewaa_m3_sim15 <- ewaa_m3[, c("real_year", "age", "ewaa_fyear_fcoh15")]
head(ewaa_m3_sim15,n=800)

#rename the columns as SS3 format ewaa.ss
ewaa_m3_sim15 <- ewaa_m3_sim15 %>%
  rename(year = real_year, pred_weight = ewaa_fyear_fcoh15)
colnames(ewaa_m3_sim15) <- c("year", "age", "pred_weight")
head(ewaa_m3_sim15,n=800)

# Convert the year column to numeric if necessary
ewaa_m3_sim15$year <- as.numeric(ewaa_m3_sim15$year)

# Run the pad_weight_at_age function and inspect the output
padded_data15 <- pad_weight_at_age(ewaa_m3_sim15)
str(padded_data15)

hakedataUSA:::write_wtatage_file(file="ewaa_sim15.ss",data=pad_weight_at_age(ewaa_m3_sim15),maturity = maturity_at_age)

###################################################################################
# Sim 16
set.seed(1201)
devs_fyear16 <- rnorm(nT, mean=0, sd=sd_fyear)
devs_fyear_fut16 <- rep(devs_fyear16, each=16)
devs_fcohort16 <- rnorm(nT_23, mean=0, sd=sd_fcohort)
devs_fcohort_past_future16 <- c(fcohort_estimate_values, devs_fcohort16)
estimate_fcohort_all16 <- data.frame(year = years_cohort, estimate = devs_fcohort_past_future16)
result_matrix16 <- matrix(0, nrow = 17, ncol = num_cols)
result_matrix16[1, ] <- format(as.integer(estimate_fcohort_all16$year), nsmall = 0)
result_matrix16[2, ] <- estimate_fcohort_all16$estimate
print(result_matrix16)
num_rows16 <- nrow(result_matrix16)
num_cols16 <- ncol(result_matrix16)
value_to_repeat16 <- result_matrix16[2,]
start_row16 <- 3
start_col16 <- 2
end_row16 <- num_rows16
end_col16 <- num_cols16
for (i in start_row16:end_row16) {
  for (j in start_col16:end_col16) {
    if (j >= (i - start_row16 + start_col16)) {
      pattern_index16 <- (j - (i - start_row16 + start_col16)) %% 32 + 1
      result_matrix16[i, j] <- value_to_repeat16[pattern_index16]
    }
  }
}
print(result_matrix16)
extracted_values16 <- as.numeric(result_matrix16[2:17, 16:33])
extracted_values16 <- extracted_values16[-1]
devs_fcohort_all16 <- c(devs_fcohort_past, extracted_values16)
devs_fyear_all16 <- c(devs_fyear_past, devs_fyear_fut16)
ewaa_m3$fyear_devs16 <- devs_fyear_all16
ewaa_m3$fcohort_devs16 <- devs_fcohort_all16
ewaa_m3$ewaa_fyear16 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fyear_devs16)
ewaa_m3$ewaa_fcohort16 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fcohort_devs16)
ewaa_m3$ewaa_fyear_fcoh16 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fyear_devs16 + ewaa_m3$fcohort_devs16)
head(ewaa_m3, n=800)

# selected the predict weight with fyear and cohort effect only.
ewaa_m3_sim16 <- ewaa_m3[, c("real_year", "age", "ewaa_fyear_fcoh16")]
head(ewaa_m3_sim16,n=800)

#rename the columns as SS3 format ewaa.ss
ewaa_m3_sim16 <- ewaa_m3_sim16 %>%
  rename(year = real_year, pred_weight = ewaa_fyear_fcoh16)
colnames(ewaa_m3_sim16) <- c("year", "age", "pred_weight")
head(ewaa_m3_sim16,n=800)

# Convert the year column to numeric if necessary
ewaa_m3_sim16$year <- as.numeric(ewaa_m3_sim16$year)

# Run the pad_weight_at_age function and inspect the output
padded_data16 <- pad_weight_at_age(ewaa_m3_sim16)
str(padded_data16)

hakedataUSA:::write_wtatage_file(file="ewaa_sim16.ss",data=pad_weight_at_age(ewaa_m3_sim16),maturity = maturity_at_age)

######################################################################################
# Sim 17
set.seed(1311)
devs_fyear17 <- rnorm(nT, mean=0, sd=sd_fyear)
devs_fyear_fut17 <- rep(devs_fyear17, each=16)
devs_fcohort17 <- rnorm(nT_23, mean=0, sd=sd_fcohort)
devs_fcohort_past_future17 <- c(fcohort_estimate_values, devs_fcohort17)
estimate_fcohort_all17 <- data.frame(year = years_cohort, estimate = devs_fcohort_past_future17)
result_matrix17 <- matrix(0, nrow = 17, ncol = num_cols)
result_matrix17[1, ] <- format(as.integer(estimate_fcohort_all17$year), nsmall = 0)
result_matrix17[2, ] <- estimate_fcohort_all17$estimate
print(result_matrix17)
num_rows17 <- nrow(result_matrix17)
num_cols17 <- ncol(result_matrix17)
value_to_repeat17 <- result_matrix17[2,]
start_row17 <- 3
start_col17 <- 2
end_row17 <- num_rows17
end_col17 <- num_cols17
for (i in start_row17:end_row17) {
  for (j in start_col17:end_col17) {
    if (j >= (i - start_row17 + start_col17)) {
      pattern_index17 <- (j - (i - start_row17 + start_col17)) %% 32 + 1
      result_matrix17[i, j] <- value_to_repeat17[pattern_index17]
    }
  }
}
print(result_matrix17)
extracted_values17 <- as.numeric(result_matrix17[2:17, 16:33])
extracted_values17 <- extracted_values17[-1]
devs_fcohort_all17 <- c(devs_fcohort_past, extracted_values17)
devs_fyear_all17 <- c(devs_fyear_past, devs_fyear_fut17)
ewaa_m3$fyear_devs17 <- devs_fyear_all17
ewaa_m3$fcohort_devs17 <- devs_fcohort_all17
ewaa_m3$ewaa_fyear17 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fyear_devs17)
ewaa_m3$ewaa_fcohort17 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fcohort_devs17)
ewaa_m3$ewaa_fyear_fcoh17 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fyear_devs17 + ewaa_m3$fcohort_devs17)
head(ewaa_m3, n=800)

# selected the predict weight with fyear and cohort effect only.
ewaa_m3_sim17 <- ewaa_m3[, c("real_year", "age", "ewaa_fyear_fcoh17")]
head(ewaa_m3_sim17,n=800)

#rename the columns as SS3 format ewaa.ss
ewaa_m3_sim17 <- ewaa_m3_sim17 %>%
  rename(year = real_year, pred_weight = ewaa_fyear_fcoh17)
colnames(ewaa_m3_sim17) <- c("year", "age", "pred_weight")
head(ewaa_m3_sim17,n=800)

# Convert the year column to numeric if necessary
ewaa_m3_sim17$year <- as.numeric(ewaa_m3_sim17$year)

# Run the pad_weight_at_age function and inspect the output
padded_data17 <- pad_weight_at_age(ewaa_m3_sim17)
str(padded_data17)

hakedataUSA:::write_wtatage_file(file="ewaa_sim17.ss",data=pad_weight_at_age(ewaa_m3_sim17),maturity = maturity_at_age)

####################################################################################
# Sim 18
set.seed(1421)
devs_fyear18 <- rnorm(nT, mean=0, sd=sd_fyear)
devs_fyear_fut18 <- rep(devs_fyear18, each=16)
devs_fcohort18 <- rnorm(nT_23, mean=0, sd=sd_fcohort)
devs_fcohort_past_future18 <- c(fcohort_estimate_values, devs_fcohort18)
estimate_fcohort_all18 <- data.frame(year = years_cohort, estimate = devs_fcohort_past_future18)
result_matrix18 <- matrix(0, nrow = 17, ncol = num_cols)
result_matrix18[1, ] <- format(as.integer(estimate_fcohort_all18$year), nsmall = 0)
result_matrix18[2, ] <- estimate_fcohort_all18$estimate
print(result_matrix18)
num_rows18 <- nrow(result_matrix18)
num_cols18 <- ncol(result_matrix18)
value_to_repeat18 <- result_matrix18[2,]
start_row18 <- 3
start_col18 <- 2
end_row18 <- num_rows18
end_col18 <- num_cols18
for (i in start_row18:end_row18) {
  for (j in start_col18:end_col18) {
    if (j >= (i - start_row18 + start_col18)) {
      pattern_index18 <- (j - (i - start_row18 + start_col18)) %% 32 + 1
      result_matrix18[i, j] <- value_to_repeat18[pattern_index18]
    }
  }
}
print(result_matrix18)
extracted_values18 <- as.numeric(result_matrix18[2:17, 16:33])
extracted_values18 <- extracted_values18[-1]
devs_fcohort_all18 <- c(devs_fcohort_past, extracted_values18)
devs_fyear_all18 <- c(devs_fyear_past, devs_fyear_fut18)
ewaa_m3$fyear_devs18 <- devs_fyear_all18
ewaa_m3$fcohort_devs18 <- devs_fcohort_all18
ewaa_m3$ewaa_fyear18 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fyear_devs18)
ewaa_m3$ewaa_fcohort18 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fcohort_devs18)
ewaa_m3$ewaa_fyear_fcoh18 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fyear_devs18 + ewaa_m3$fcohort_devs18)
head(ewaa_m3, n=800)

# selected the predict weight with fyear and cohort effect only.
ewaa_m3_sim18 <- ewaa_m3[, c("real_year", "age", "ewaa_fyear_fcoh18")]
head(ewaa_m3_sim18,n=800)

#rename the columns as SS3 format ewaa.ss
ewaa_m3_sim18 <- ewaa_m3_sim18 %>%
  rename(year = real_year, pred_weight = ewaa_fyear_fcoh18)
colnames(ewaa_m3_sim18) <- c("year", "age", "pred_weight")
head(ewaa_m3_sim18,n=800)

# Convert the year column to numeric if necessary
ewaa_m3_sim18$year <- as.numeric(ewaa_m3_sim18$year)

# Run the pad_weight_at_age function and inspect the output
padded_data18 <- pad_weight_at_age(ewaa_m3_sim18)
str(padded_data18)

hakedataUSA:::write_wtatage_file(file="ewaa_sim18.ss",data=pad_weight_at_age(ewaa_m3_sim18),maturity = maturity_at_age)

###################################################################################
# Sim 19
set.seed(1531)
devs_fyear19 <- rnorm(nT, mean=0, sd=sd_fyear)
devs_fyear_fut19 <- rep(devs_fyear19, each=16)
devs_fcohort19 <- rnorm(nT_23, mean=0, sd=sd_fcohort)
devs_fcohort_past_future19 <- c(fcohort_estimate_values, devs_fcohort19)
estimate_fcohort_all19 <- data.frame(year = years_cohort, estimate = devs_fcohort_past_future19)
result_matrix19 <- matrix(0, nrow = 17, ncol = num_cols)
result_matrix19[1, ] <- format(as.integer(estimate_fcohort_all19$year), nsmall = 0)
result_matrix19[2, ] <- estimate_fcohort_all19$estimate
print(result_matrix19)
num_rows19 <- nrow(result_matrix19)
num_cols19 <- ncol(result_matrix19)
value_to_repeat19 <- result_matrix19[2,]
start_row19 <- 3
start_col19 <- 2
end_row19 <- num_rows19
end_col19 <- num_cols19
for (i in start_row19:end_row19) {
  for (j in start_col19:end_col19) {
    if (j >= (i - start_row19 + start_col19)) {
      pattern_index19 <- (j - (i - start_row19 + start_col19)) %% 32 + 1
      result_matrix19[i, j] <- value_to_repeat19[pattern_index19]
    }
  }
}
print(result_matrix19)
extracted_values19 <- as.numeric(result_matrix19[2:17, 16:33])
extracted_values19 <- extracted_values19[-1]
devs_fcohort_all19 <- c(devs_fcohort_past, extracted_values19)
devs_fyear_all19 <- c(devs_fyear_past, devs_fyear_fut19)
ewaa_m3$fyear_devs19 <- devs_fyear_all19
ewaa_m3$fcohort_devs19 <- devs_fcohort_all19
ewaa_m3$ewaa_fyear19 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fyear_devs19)
ewaa_m3$ewaa_fcohort19 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fcohort_devs19)
ewaa_m3$ewaa_fyear_fcoh19 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fyear_devs19 + ewaa_m3$fcohort_devs19)
head(ewaa_m3, n=800)

# selected the predict weight with fyear and cohort effect only.
ewaa_m3_sim19 <- ewaa_m3[, c("real_year", "age", "ewaa_fyear_fcoh19")]
head(ewaa_m3_sim19,n=800)

#rename the columns as SS3 format ewaa.ss
ewaa_m3_sim19 <- ewaa_m3_sim19 %>%
  rename(year = real_year, pred_weight = ewaa_fyear_fcoh19)
colnames(ewaa_m3_sim19) <- c("year", "age", "pred_weight")
head(ewaa_m3_sim19,n=800)

# Convert the year column to numeric if necessary
ewaa_m3_sim19$year <- as.numeric(ewaa_m3_sim19$year)

# Run the pad_weight_at_age function and inspect the output
padded_data19 <- pad_weight_at_age(ewaa_m3_sim19)
str(padded_data19)

hakedataUSA:::write_wtatage_file(file="ewaa_sim19.ss",data=pad_weight_at_age(ewaa_m3_sim19),maturity = maturity_at_age)

##############################################################################
# Sim 20
set.seed(1641)
devs_fyear20 <- rnorm(nT, mean=0, sd=sd_fyear)
devs_fyear_fut20 <- rep(devs_fyear20, each=16)
devs_fcohort20 <- rnorm(nT_23, mean=0, sd=sd_fcohort)
devs_fcohort_past_future20 <- c(fcohort_estimate_values, devs_fcohort20)
estimate_fcohort_all20 <- data.frame(year = years_cohort, estimate = devs_fcohort_past_future20)
result_matrix20 <- matrix(0, nrow = 17, ncol = num_cols)
result_matrix20[1, ] <- format(as.integer(estimate_fcohort_all20$year), nsmall = 0)
result_matrix20[2, ] <- estimate_fcohort_all20$estimate
print(result_matrix20)
num_rows20 <- nrow(result_matrix20)
num_cols20 <- ncol(result_matrix20)
value_to_repeat20 <- result_matrix20[2,]
start_row20 <- 3
start_col20 <- 2
end_row20 <- num_rows20
end_col20 <- num_cols20
for (i in start_row20:end_row20) {
  for (j in start_col20:end_col20) {
    if (j >= (i - start_row20 + start_col20)) {
      pattern_index20 <- (j - (i - start_row20 + start_col20)) %% 32 + 1
      result_matrix20[i, j] <- value_to_repeat20[pattern_index20]
    }
  }
}
print(result_matrix20)
extracted_values20 <- as.numeric(result_matrix20[2:17, 16:33])
extracted_values20 <- extracted_values20[-1]
devs_fcohort_all20 <- c(devs_fcohort_past, extracted_values20)
devs_fyear_all20 <- c(devs_fyear_past, devs_fyear_fut20)
ewaa_m3$fyear_devs20 <- devs_fyear_all20
ewaa_m3$fcohort_devs20 <- devs_fcohort_all20
ewaa_m3$ewaa_fyear20 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fyear_devs20)
ewaa_m3$ewaa_fcohort20 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fcohort_devs20)
ewaa_m3$ewaa_fyear_fcoh20 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fyear_devs20 + ewaa_m3$fcohort_devs20)
head(ewaa_m3, n=800)

# selected the predict weight with fyear and cohort effect only.
ewaa_m3_sim20 <- ewaa_m3[, c("real_year", "age", "ewaa_fyear_fcoh20")]
head(ewaa_m3_sim20,n=800)

#rename the columns as SS3 format ewaa.ss
ewaa_m3_sim20 <- ewaa_m3_sim20 %>%
  rename(year = real_year, pred_weight = ewaa_fyear_fcoh20)
colnames(ewaa_m3_sim20) <- c("year", "age", "pred_weight")
head(ewaa_m3_sim20,n=800)

# Convert the year column to numeric if necessary
ewaa_m3_sim20$year <- as.numeric(ewaa_m3_sim20$year)

# Run the pad_weight_at_age function and inspect the output
padded_data20 <- pad_weight_at_age(ewaa_m3_sim20)
str(padded_data20)

hakedataUSA:::write_wtatage_file(file="ewaa_sim20.ss",data=pad_weight_at_age(ewaa_m3_sim20),maturity = maturity_at_age)

#############################################################################
# Sim 21
set.seed(1751)
devs_fyear21 <- rnorm(nT, mean=0, sd=sd_fyear)
devs_fyear_fut21 <- rep(devs_fyear21, each=16)
devs_fcohort21 <- rnorm(nT_23, mean=0, sd=sd_fcohort)
devs_fcohort_past_future21 <- c(fcohort_estimate_values, devs_fcohort21)
estimate_fcohort_all21 <- data.frame(year = years_cohort, estimate = devs_fcohort_past_future21)
result_matrix21 <- matrix(0, nrow = 17, ncol = num_cols)
result_matrix21[1, ] <- format(as.integer(estimate_fcohort_all21$year), nsmall = 0)
result_matrix21[2, ] <- estimate_fcohort_all21$estimate
print(result_matrix21)
num_rows21 <- nrow(result_matrix21)
num_cols21 <- ncol(result_matrix21)
value_to_repeat21 <- result_matrix21[2,]
start_row21 <- 3
start_col21 <- 2
end_row21 <- num_rows21
end_col21 <- num_cols21
for (i in start_row21:end_row21) {
  for (j in start_col21:end_col21) {
    if (j >= (i - start_row21 + start_col21)) {
      pattern_index21 <- (j - (i - start_row21 + start_col21)) %% 32 + 1
      result_matrix21[i, j] <- value_to_repeat21[pattern_index21]
    }
  }
}
print(result_matrix21)
extracted_values21 <- as.numeric(result_matrix21[2:17, 16:33])
extracted_values21 <- extracted_values21[-1]
devs_fcohort_all21 <- c(devs_fcohort_past, extracted_values21)
devs_fyear_all21 <- c(devs_fyear_past, devs_fyear_fut21)
ewaa_m3$fyear_devs21 <- devs_fyear_all21
ewaa_m3$fcohort_devs21 <- devs_fcohort_all21
ewaa_m3$ewaa_fyear21 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fyear_devs21)
ewaa_m3$ewaa_fcohort21 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fcohort_devs21)
ewaa_m3$ewaa_fyear_fcoh21 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fyear_devs21 + ewaa_m3$fcohort_devs21)
head(ewaa_m3, n=800)

# selected the predict weight with fyear and cohort effect only.
ewaa_m3_sim21 <- ewaa_m3[, c("real_year", "age", "ewaa_fyear_fcoh21")]
head(ewaa_m3_sim21,n=800)

#rename the columns as SS3 format ewaa.ss
ewaa_m3_sim21 <- ewaa_m3_sim21 %>%
  rename(year = real_year, pred_weight = ewaa_fyear_fcoh21)
colnames(ewaa_m3_sim21) <- c("year", "age", "pred_weight")
head(ewaa_m3_sim21,n=800)

# Convert the year column to numeric if necessary
ewaa_m3_sim21$year <- as.numeric(ewaa_m3_sim21$year)

# Run the pad_weight_at_age function and inspect the output
padded_data21 <- pad_weight_at_age(ewaa_m3_sim21)
str(padded_data21)

hakedataUSA:::write_wtatage_file(file="ewaa_sim21.ss",data=pad_weight_at_age(ewaa_m3_sim21),maturity = maturity_at_age)

#############################################################################
# Sim 22
set.seed(1861)
devs_fyear22 <- rnorm(nT, mean=0, sd=sd_fyear)
devs_fyear_fut22 <- rep(devs_fyear22, each=16)
devs_fcohort22 <- rnorm(nT_23, mean=0, sd=sd_fcohort)
devs_fcohort_past_future22 <- c(fcohort_estimate_values, devs_fcohort22)
estimate_fcohort_all22 <- data.frame(year = years_cohort, estimate = devs_fcohort_past_future22)
result_matrix22 <- matrix(0, nrow = 17, ncol = num_cols)
result_matrix22[1, ] <- format(as.integer(estimate_fcohort_all22$year), nsmall = 0)
result_matrix22[2, ] <- estimate_fcohort_all22$estimate
print(result_matrix22)
num_rows22 <- nrow(result_matrix22)
num_cols22 <- ncol(result_matrix22)
value_to_repeat22 <- result_matrix22[2,]
start_row22 <- 3
start_col22 <- 2
end_row22 <- num_rows22
end_col22 <- num_cols22
for (i in start_row22:end_row22) {
  for (j in start_col22:end_col22) {
    if (j >= (i - start_row22 + start_col22)) {
      pattern_index22 <- (j - (i - start_row22 + start_col22)) %% 32 + 1
      result_matrix22[i, j] <- value_to_repeat22[pattern_index22]
    }
  }
}
print(result_matrix22)
extracted_values22 <- as.numeric(result_matrix22[2:17, 16:33])
extracted_values22 <- extracted_values22[-1]
devs_fcohort_all22 <- c(devs_fcohort_past, extracted_values22)
devs_fyear_all22 <- c(devs_fyear_past, devs_fyear_fut22)
ewaa_m3$fyear_devs22 <- devs_fyear_all22
ewaa_m3$fcohort_devs22 <- devs_fcohort_all22
ewaa_m3$ewaa_fyear22 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fyear_devs22)
ewaa_m3$ewaa_fcohort22 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fcohort_devs22)
ewaa_m3$ewaa_fyear_fcoh22 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fyear_devs22 + ewaa_m3$fcohort_devs22)
head(ewaa_m3, n=800)

# selected the predict weight with fyear and cohort effect only.
ewaa_m3_sim22 <- ewaa_m3[, c("real_year", "age", "ewaa_fyear_fcoh22")]
head(ewaa_m3_sim22,n=800)

#rename the columns as SS3 format ewaa.ss
ewaa_m3_sim22 <- ewaa_m3_sim22 %>%
  rename(year = real_year, pred_weight = ewaa_fyear_fcoh22)
colnames(ewaa_m3_sim22) <- c("year", "age", "pred_weight")
head(ewaa_m3_sim22,n=800)

# Convert the year column to numeric if necessary
ewaa_m3_sim22$year <- as.numeric(ewaa_m3_sim22$year)

# Run the pad_weight_at_age function and inspect the output
padded_data22 <- pad_weight_at_age(ewaa_m3_sim22)
str(padded_data22)

hakedataUSA:::write_wtatage_file(file="ewaa_sim22.ss",data=pad_weight_at_age(ewaa_m3_sim22),maturity = maturity_at_age)

###################################################################################
# Sim 23
set.seed(1971)
devs_fyear23 <- rnorm(nT, mean=0, sd=sd_fyear)
devs_fyear_fut23 <- rep(devs_fyear23, each=16)
devs_fcohort23 <- rnorm(nT_23, mean=0, sd=sd_fcohort)
devs_fcohort_past_future23 <- c(fcohort_estimate_values, devs_fcohort23)
estimate_fcohort_all23 <- data.frame(year = years_cohort, estimate = devs_fcohort_past_future23)
result_matrix23 <- matrix(0, nrow = 17, ncol = num_cols)
result_matrix23[1, ] <- format(as.integer(estimate_fcohort_all23$year), nsmall = 0)
result_matrix23[2, ] <- estimate_fcohort_all23$estimate
print(result_matrix23)
num_rows23 <- nrow(result_matrix23)
num_cols23 <- ncol(result_matrix23)
value_to_repeat23 <- result_matrix23[2,]
start_row23 <- 3
start_col23 <- 2
end_row23 <- num_rows23
end_col23 <- num_cols23
for (i in start_row23:end_row23) {
  for (j in start_col23:end_col23) {
    if (j >= (i - start_row23 + start_col23)) {
      pattern_index23 <- (j - (i - start_row23 + start_col23)) %% 32 + 1
      result_matrix23[i, j] <- value_to_repeat23[pattern_index23]
    }
  }
}
print(result_matrix23)
extracted_values23 <- as.numeric(result_matrix23[2:17, 16:33])
extracted_values23 <- extracted_values23[-1]
devs_fcohort_all23 <- c(devs_fcohort_past, extracted_values23)
devs_fyear_all23 <- c(devs_fyear_past, devs_fyear_fut23)
ewaa_m3$fyear_devs23 <- devs_fyear_all23
ewaa_m3$fcohort_devs23 <- devs_fcohort_all23
ewaa_m3$ewaa_fyear23 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fyear_devs23)
ewaa_m3$ewaa_fcohort23 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fcohort_devs23)
ewaa_m3$ewaa_fyear_fcoh23 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fyear_devs23 + ewaa_m3$fcohort_devs23)
head(ewaa_m3, n=800)

# selected the predict weight with fyear and cohort effect only.
ewaa_m3_sim23 <- ewaa_m3[, c("real_year", "age", "ewaa_fyear_fcoh23")]
head(ewaa_m3_sim23,n=800)

#rename the columns as SS3 format ewaa.ss
ewaa_m3_sim23 <- ewaa_m3_sim23 %>%
  rename(year = real_year, pred_weight = ewaa_fyear_fcoh23)
colnames(ewaa_m3_sim23) <- c("year", "age", "pred_weight")
head(ewaa_m3_sim23,n=800)

# Convert the year column to numeric if necessary
ewaa_m3_sim23$year <- as.numeric(ewaa_m3_sim23$year)

# Run the pad_weight_at_age function and inspect the output
padded_data23 <- pad_weight_at_age(ewaa_m3_sim23)
str(padded_data23)

hakedataUSA:::write_wtatage_file(file="ewaa_sim23.ss",data=pad_weight_at_age(ewaa_m3_sim23),maturity = maturity_at_age)

##############################################################################
# Sim 24
set.seed(2081)
devs_fyear24 <- rnorm(nT, mean=0, sd=sd_fyear)
devs_fyear_fut24 <- rep(devs_fyear24, each=16)
devs_fcohort24 <- rnorm(nT_23, mean=0, sd=sd_fcohort)
devs_fcohort_past_future24 <- c(fcohort_estimate_values, devs_fcohort24)
estimate_fcohort_all24 <- data.frame(year = years_cohort, estimate = devs_fcohort_past_future24)
result_matrix24 <- matrix(0, nrow = 17, ncol = num_cols)
result_matrix24[1, ] <- format(as.integer(estimate_fcohort_all24$year), nsmall = 0)
result_matrix24[2, ] <- estimate_fcohort_all24$estimate
print(result_matrix24)
num_rows24 <- nrow(result_matrix24)
num_cols24 <- ncol(result_matrix24)
value_to_repeat24 <- result_matrix24[2,]
start_row24 <- 3
start_col24 <- 2
end_row24 <- num_rows24
end_col24 <- num_cols24
for (i in start_row24:end_row24) {
  for (j in start_col24:end_col24) {
    if (j >= (i - start_row24 + start_col24)) {
      pattern_index24 <- (j - (i - start_row24 + start_col24)) %% 32 + 1
      result_matrix24[i, j] <- value_to_repeat24[pattern_index24]
    }
  }
}
print(result_matrix24)
extracted_values24 <- as.numeric(result_matrix24[2:17, 16:33])
extracted_values24 <- extracted_values24[-1]
devs_fcohort_all24 <- c(devs_fcohort_past, extracted_values24)
devs_fyear_all24 <- c(devs_fyear_past, devs_fyear_fut24)
ewaa_m3$fyear_devs24 <- devs_fyear_all24
ewaa_m3$fcohort_devs24 <- devs_fcohort_all24
ewaa_m3$ewaa_fyear24 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fyear_devs24)
ewaa_m3$ewaa_fcohort24 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fcohort_devs24)
ewaa_m3$ewaa_fyear_fcoh24 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fyear_devs24 + ewaa_m3$fcohort_devs24)
head(ewaa_m3, n=800)

# selected the predict weight with fyear and cohort effect only.
ewaa_m3_sim24 <- ewaa_m3[, c("real_year", "age", "ewaa_fyear_fcoh24")]
head(ewaa_m3_sim24,n=800)

#rename the columns as SS3 format ewaa.ss
ewaa_m3_sim24 <- ewaa_m3_sim24 %>%
  rename(year = real_year, pred_weight = ewaa_fyear_fcoh24)
colnames(ewaa_m3_sim24) <- c("year", "age", "pred_weight")
head(ewaa_m3_sim24,n=800)

# Convert the year column to numeric if necessary
ewaa_m3_sim24$year <- as.numeric(ewaa_m3_sim24$year)

# Run the pad_weight_at_age function and inspect the output
padded_data24 <- pad_weight_at_age(ewaa_m3_sim24)
str(padded_data24)

hakedataUSA:::write_wtatage_file(file="ewaa_sim24.ss",data=pad_weight_at_age(ewaa_m3_sim24),maturity = maturity_at_age)

#################################################################################
# Sim 25
set.seed(2191)
devs_fyear25 <- rnorm(nT, mean=0, sd=sd_fyear)
devs_fyear_fut25 <- rep(devs_fyear25, each=16)
devs_fcohort25 <- rnorm(nT_23, mean=0, sd=sd_fcohort)
devs_fcohort_past_future25 <- c(fcohort_estimate_values, devs_fcohort25)
estimate_fcohort_all25 <- data.frame(year = years_cohort, estimate = devs_fcohort_past_future25)
result_matrix25 <- matrix(0, nrow = 17, ncol = num_cols)
result_matrix25[1, ] <- format(as.integer(estimate_fcohort_all25$year), nsmall = 0)
result_matrix25[2, ] <- estimate_fcohort_all25$estimate
print(result_matrix25)
num_rows25 <- nrow(result_matrix25)
num_cols25 <- ncol(result_matrix25)
value_to_repeat25 <- result_matrix25[2,]
start_row25 <- 3
start_col25 <- 2
end_row25 <- num_rows25
end_col25 <- num_cols25
for (i in start_row25:end_row25) {
  for (j in start_col25:end_col25) {
    if (j >= (i - start_row25 + start_col25)) {
      pattern_index25 <- (j - (i - start_row25 + start_col25)) %% 32 + 1
      result_matrix25[i, j] <- value_to_repeat25[pattern_index25]
    }
  }
}
print(result_matrix25)
extracted_values25 <- as.numeric(result_matrix25[2:17, 16:33])
extracted_values25 <- extracted_values25[-1]
devs_fcohort_all25 <- c(devs_fcohort_past, extracted_values25)
devs_fyear_all25 <- c(devs_fyear_past, devs_fyear_fut25)
ewaa_m3$fyear_devs25 <- devs_fyear_all25
ewaa_m3$fcohort_devs25 <- devs_fcohort_all25
ewaa_m3$ewaa_fyear25 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fyear_devs25)
ewaa_m3$ewaa_fcohort25 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fcohort_devs25)
ewaa_m3$ewaa_fyear_fcoh25 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fyear_devs25 + ewaa_m3$fcohort_devs25)
head(ewaa_m3, n=800)

# selected the predict weight with fyear and cohort effect only.
ewaa_m3_sim25 <- ewaa_m3[, c("real_year", "age", "ewaa_fyear_fcoh25")]
head(ewaa_m3_sim25,n=800)

#rename the columns as SS3 format ewaa.ss
ewaa_m3_sim25 <- ewaa_m3_sim25 %>%
  rename(year = real_year, pred_weight = ewaa_fyear_fcoh25)
colnames(ewaa_m3_sim25) <- c("year", "age", "pred_weight")
head(ewaa_m3_sim25,n=800)

# Convert the year column to numeric if necessary
ewaa_m3_sim25$year <- as.numeric(ewaa_m3_sim25$year)

# Run the pad_weight_at_age function and inspect the output
padded_data25 <- pad_weight_at_age(ewaa_m3_sim25)
str(padded_data25)

hakedataUSA:::write_wtatage_file(file="ewaa_sim25.ss",data=pad_weight_at_age(ewaa_m3_sim25),maturity = maturity_at_age)

##############################################################################
# Sim 26
set.seed(2251)
devs_fyear26 <- rnorm(nT, mean=0, sd=sd_fyear)
devs_fyear_fut26 <- rep(devs_fyear26, each=16)
devs_fcohort26 <- rnorm(nT_23, mean=0, sd=sd_fcohort)
devs_fcohort_past_future26 <- c(fcohort_estimate_values, devs_fcohort26)
estimate_fcohort_all26 <- data.frame(year = years_cohort, estimate = devs_fcohort_past_future26)
result_matrix26 <- matrix(0, nrow = 17, ncol = num_cols)
result_matrix26[1, ] <- format(as.integer(estimate_fcohort_all26$year), nsmall = 0)
result_matrix26[2, ] <- estimate_fcohort_all26$estimate
print(result_matrix26)
num_rows26 <- nrow(result_matrix26)
num_cols26 <- ncol(result_matrix26)
value_to_repeat26 <- result_matrix26[2,]
start_row26 <- 3
start_col26 <- 2
end_row26 <- num_rows26
end_col26 <- num_cols26
for (i in start_row26:end_row26) {
  for (j in start_col26:end_col26) {
    if (j >= (i - start_row26 + start_col26)) {
      pattern_index26 <- (j - (i - start_row26 + start_col26)) %% 32 + 1
      result_matrix26[i, j] <- value_to_repeat26[pattern_index26]
    }
  }
}
print(result_matrix26)
extracted_values26 <- as.numeric(result_matrix26[2:17, 16:33])
extracted_values26 <- extracted_values26[-1]
devs_fcohort_all26 <- c(devs_fcohort_past, extracted_values26)
devs_fyear_all26 <- c(devs_fyear_past, devs_fyear_fut26)
ewaa_m3$fyear_devs26 <- devs_fyear_all26
ewaa_m3$fcohort_devs26 <- devs_fcohort_all26
ewaa_m3$ewaa_fyear26 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fyear_devs26)
ewaa_m3$ewaa_fcohort26 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fcohort_devs26)
ewaa_m3$ewaa_fyear_fcoh26 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fyear_devs26 + ewaa_m3$fcohort_devs26)
head(ewaa_m3, n=800)

# selected the predict weight with fyear and cohort effect only.
ewaa_m3_sim26 <- ewaa_m3[, c("real_year", "age", "ewaa_fyear_fcoh26")]
head(ewaa_m3_sim26,n=800)

#rename the columns as SS3 format ewaa.ss
ewaa_m3_sim26 <- ewaa_m3_sim26 %>%
  rename(year = real_year, pred_weight = ewaa_fyear_fcoh26)
colnames(ewaa_m3_sim26) <- c("year", "age", "pred_weight")
head(ewaa_m3_sim26,n=800)

# Convert the year column to numeric if necessary
ewaa_m3_sim26$year <- as.numeric(ewaa_m3_sim26$year)

# Run the pad_weight_at_age function and inspect the output
padded_data26 <- pad_weight_at_age(ewaa_m3_sim26)
str(padded_data26)

hakedataUSA:::write_wtatage_file(file="ewaa_sim26.ss",data=pad_weight_at_age(ewaa_m3_sim26),maturity = maturity_at_age)

################################################################################
# Sim 27
set.seed(2407)
devs_fyear27 <- rnorm(nT, mean = 0, sd = sd_fyear)
devs_fyear_fut27 <- rep(devs_fyear27, each = 16)
devs_fcohort27 <- rnorm(nT_23, mean = 0, sd = sd_fcohort)
devs_fcohort_past_future27 <- c(fcohort_estimate_values, devs_fcohort27)
estimate_fcohort_all27 <- data.frame(year = years_cohort, estimate = devs_fcohort_past_future27)
result_matrix27 <- matrix(0, nrow = 17, ncol = num_cols)
result_matrix27[1, ] <- format(as.integer(estimate_fcohort_all27$year), nsmall = 0)
result_matrix27[2, ] <- estimate_fcohort_all27$estimate
print(result_matrix27)
num_rows27 <- nrow(result_matrix27)
num_cols27 <- ncol(result_matrix27)
value_to_repeat27 <- result_matrix27[2, ]
start_row27 <- 3
start_col27 <- 2
end_row27 <- num_rows27
end_col27 <- num_cols27
for (i in start_row27:end_row27) {
  for (j in start_col27:end_col27) {
    if (j >= (i - start_row27 + start_col27)) {
      pattern_index27 <- (j - (i - start_row27 + start_col27)) %% 32 + 1
      result_matrix27[i, j] <- value_to_repeat27[pattern_index27]
    }
  }
}
print(result_matrix27)
extracted_values27 <- as.numeric(result_matrix27[2:17, 16:33])
extracted_values27 <- extracted_values27[-1]
devs_fcohort_all27 <- c(devs_fcohort_past, extracted_values27)
devs_fyear_all27 <- c(devs_fyear_past, devs_fyear_fut27)
ewaa_m3$fyear_devs27 <- devs_fyear_all27
ewaa_m3$fcohort_devs27 <- devs_fcohort_all27
ewaa_m3$ewaa_fyear27 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fyear_devs27)
ewaa_m3$ewaa_fcohort27 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fcohort_devs27)
ewaa_m3$ewaa_fyear_fcoh27 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fyear_devs27 + ewaa_m3$fcohort_devs27)
head(ewaa_m3, n = 800)

# selected the predict weight with fyear and cohort effect only.
ewaa_m3_sim27 <- ewaa_m3[, c("real_year", "age", "ewaa_fyear_fcoh27")]
head(ewaa_m3_sim27,n=800)

#rename the columns as SS3 format ewaa.ss
ewaa_m3_sim27 <- ewaa_m3_sim27 %>%
  rename(year = real_year, pred_weight = ewaa_fyear_fcoh27)
colnames(ewaa_m3_sim27) <- c("year", "age", "pred_weight")
head(ewaa_m3_sim27,n=800)

# Convert the year column to numeric if necessary
ewaa_m3_sim27$year <- as.numeric(ewaa_m3_sim27$year)

# Run the pad_weight_at_age function and inspect the output
padded_data27 <- pad_weight_at_age(ewaa_m3_sim27)
str(padded_data27)

hakedataUSA:::write_wtatage_file(file="ewaa_sim27.ss",data=pad_weight_at_age(ewaa_m3_sim27),maturity = maturity_at_age)

#######################################################################################
# Sim 28
set.seed(2521)
devs_fyear28 <- rnorm(nT, mean = 0, sd = sd_fyear)
devs_fyear_fut28 <- rep(devs_fyear28, each = 16)
devs_fcohort28 <- rnorm(nT_23, mean = 0, sd = sd_fcohort)
devs_fcohort_past_future28 <- c(fcohort_estimate_values, devs_fcohort28)
estimate_fcohort_all28 <- data.frame(year = years_cohort, estimate = devs_fcohort_past_future28)
result_matrix28 <- matrix(0, nrow = 17, ncol = num_cols)
result_matrix28[1, ] <- format(as.integer(estimate_fcohort_all28$year), nsmall = 0)
result_matrix28[2, ] <- estimate_fcohort_all28$estimate
print(result_matrix28)
num_rows28 <- nrow(result_matrix28)
num_cols28 <- ncol(result_matrix28)
value_to_repeat28 <- result_matrix28[2, ]
start_row28 <- 3
start_col28 <- 2
end_row28 <- num_rows28
end_col28 <- num_cols28
for (i in start_row28:end_row28) {
  for (j in start_col28:end_col28) {
    if (j >= (i - start_row28 + start_col28)) {
      pattern_index28 <- (j - (i - start_row28 + start_col28)) %% 32 + 1
      result_matrix28[i, j] <- value_to_repeat28[pattern_index28]
    }
  }
}
print(result_matrix28)
extracted_values28 <- as.numeric(result_matrix28[2:17, 16:33])
extracted_values28 <- extracted_values28[-1]
devs_fcohort_all28 <- c(devs_fcohort_past, extracted_values28)
devs_fyear_all28 <- c(devs_fyear_past, devs_fyear_fut28)
ewaa_m3$fyear_devs28 <- devs_fyear_all28
ewaa_m3$fcohort_devs28 <- devs_fcohort_all28
ewaa_m3$ewaa_fyear28 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fyear_devs28)
ewaa_m3$ewaa_fcohort28 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fcohort_devs28)
ewaa_m3$ewaa_fyear_fcoh28 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fyear_devs28 + ewaa_m3$fcohort_devs28)
head(ewaa_m3, n = 800)

# selected the predict weight with fyear and cohort effect only.
ewaa_m3_sim28 <- ewaa_m3[, c("real_year", "age", "ewaa_fyear_fcoh28")]
head(ewaa_m3_sim28,n=800)

#rename the columns as SS3 format ewaa.ss
ewaa_m3_sim28 <- ewaa_m3_sim28 %>%
  rename(year = real_year, pred_weight = ewaa_fyear_fcoh28)
colnames(ewaa_m3_sim28) <- c("year", "age", "pred_weight")
head(ewaa_m3_sim28,n=800)

# Convert the year column to numeric if necessary
ewaa_m3_sim28$year <- as.numeric(ewaa_m3_sim28$year)

# Run the pad_weight_at_age function and inspect the output
padded_data28 <- pad_weight_at_age(ewaa_m3_sim28)
str(padded_data28)

hakedataUSA:::write_wtatage_file(file="ewaa_sim28.ss",data=pad_weight_at_age(ewaa_m3_sim28),maturity = maturity_at_age)

############################################################################
# Sim 29
set.seed(2683)
devs_fyear29 <- rnorm(nT, mean = 0, sd = sd_fyear)
devs_fyear_fut29 <- rep(devs_fyear29, each = 16)
devs_fcohort29 <- rnorm(nT_23, mean = 0, sd = sd_fcohort)
devs_fcohort_past_future29 <- c(fcohort_estimate_values, devs_fcohort29)
estimate_fcohort_all29 <- data.frame(year = years_cohort, estimate = devs_fcohort_past_future29)
result_matrix29 <- matrix(0, nrow = 17, ncol = num_cols)
result_matrix29[1, ] <- format(as.integer(estimate_fcohort_all29$year), nsmall = 0)
result_matrix29[2, ] <- estimate_fcohort_all29$estimate
print(result_matrix29)
num_rows29 <- nrow(result_matrix29)
num_cols29 <- ncol(result_matrix29)
value_to_repeat29 <- result_matrix29[2, ]
start_row29 <- 3
start_col29 <- 2
end_row29 <- num_rows29
end_col29 <- num_cols29
for (i in start_row29:end_row29) {
  for (j in start_col29:end_col29) {
    if (j >= (i - start_row29 + start_col29)) {
      pattern_index29 <- (j - (i - start_row29 + start_col29)) %% 32 + 1
      result_matrix29[i, j] <- value_to_repeat29[pattern_index29]
    }
  }
}
print(result_matrix29)
extracted_values29 <- as.numeric(result_matrix29[2:17, 16:33])
extracted_values29 <- extracted_values29[-1]
devs_fcohort_all29 <- c(devs_fcohort_past, extracted_values29)
devs_fyear_all29 <- c(devs_fyear_past, devs_fyear_fut29)
ewaa_m3$fyear_devs29 <- devs_fyear_all29
ewaa_m3$fcohort_devs29 <- devs_fcohort_all29
ewaa_m3$ewaa_fyear29 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fyear_devs29)
ewaa_m3$ewaa_fcohort29 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fcohort_devs29)
ewaa_m3$ewaa_fyear_fcoh29 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fyear_devs29 + ewaa_m3$fcohort_devs29)
head(ewaa_m3, n = 800)

# selected the predict weight with fyear and cohort effect only.
ewaa_m3_sim29 <- ewaa_m3[, c("real_year", "age", "ewaa_fyear_fcoh29")]
head(ewaa_m3_sim29,n=800)

#rename the columns as SS3 format ewaa.ss
ewaa_m3_sim29 <- ewaa_m3_sim29 %>%
  rename(year = real_year, pred_weight = ewaa_fyear_fcoh29)
colnames(ewaa_m3_sim29) <- c("year", "age", "pred_weight")
head(ewaa_m3_sim29,n=800)

# Convert the year column to numeric if necessary
ewaa_m3_sim29$year <- as.numeric(ewaa_m3_sim29$year)

# Run the pad_weight_at_age function and inspect the output
padded_data29 <- pad_weight_at_age(ewaa_m3_sim29)
str(padded_data29)

hakedataUSA:::write_wtatage_file(file="ewaa_sim29.ss",data=pad_weight_at_age(ewaa_m3_sim29),maturity = maturity_at_age)

###########################################################################
# Sim 30
set.seed(2769)
devs_fyear30 <- rnorm(nT, mean = 0, sd = sd_fyear)
devs_fyear_fut30 <- rep(devs_fyear30, each = 16)
devs_fcohort30 <- rnorm(nT_23, mean = 0, sd = sd_fcohort)
devs_fcohort_past_future30 <- c(fcohort_estimate_values, devs_fcohort30)
estimate_fcohort_all30 <- data.frame(year = years_cohort, estimate = devs_fcohort_past_future30)
result_matrix30 <- matrix(0, nrow = 17, ncol = num_cols)
result_matrix30[1, ] <- format(as.integer(estimate_fcohort_all30$year), nsmall = 0)
result_matrix30[2, ] <- estimate_fcohort_all30$estimate
print(result_matrix30)
num_rows30 <- nrow(result_matrix30)
num_cols30 <- ncol(result_matrix30)
value_to_repeat30 <- result_matrix30[2, ]
start_row30 <- 3
start_col30 <- 2
end_row30 <- num_rows30
end_col30 <- num_cols30
for (i in start_row30:end_row30) {
  for (j in start_col30:end_col30) {
    if (j >= (i - start_row30 + start_col30)) {
      pattern_index30 <- (j - (i - start_row30 + start_col30)) %% 32 + 1
      result_matrix30[i, j] <- value_to_repeat30[pattern_index30]
    }
  }
}
print(result_matrix30)
extracted_values30 <- as.numeric(result_matrix30[2:17, 16:33])
extracted_values30 <- extracted_values30[-1]
devs_fcohort_all30 <- c(devs_fcohort_past, extracted_values30)
devs_fyear_all30 <- c(devs_fyear_past, devs_fyear_fut30)
ewaa_m3$fyear_devs30 <- devs_fyear_all30
ewaa_m3$fcohort_devs30 <- devs_fcohort_all30
ewaa_m3$ewaa_fyear30 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fyear_devs30)
ewaa_m3$ewaa_fcohort30 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fcohort_devs30)
ewaa_m3$ewaa_fyear_fcoh30 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fyear_devs30 + ewaa_m3$fcohort_devs30)
head(ewaa_m3, n = 800)

# selected the predict weight with fyear and cohort effect only.
ewaa_m3_sim30 <- ewaa_m3[, c("real_year", "age", "ewaa_fyear_fcoh30")]
head(ewaa_m3_sim30,n=800)

#rename the columns as SS3 format ewaa.ss
ewaa_m3_sim30 <- ewaa_m3_sim30 %>%
  rename(year = real_year, pred_weight = ewaa_fyear_fcoh30)
colnames(ewaa_m3_sim30) <- c("year", "age", "pred_weight")
head(ewaa_m3_sim30,n=800)

# Convert the year column to numeric if necessary
ewaa_m3_sim30$year <- as.numeric(ewaa_m3_sim30$year)

# Run the pad_weight_at_age function and inspect the output
padded_data30 <- pad_weight_at_age(ewaa_m3_sim30)
str(padded_data30)

hakedataUSA:::write_wtatage_file(file="ewaa_sim30.ss",data=pad_weight_at_age(ewaa_m3_sim30),maturity = maturity_at_age)

#######################################################################################
# Sim 31
set.seed(2891)
devs_fyear31 <- rnorm(nT, mean = 0, sd = sd_fyear)
devs_fyear_fut31 <- rep(devs_fyear31, each = 16)
devs_fcohort31 <- rnorm(nT_23, mean = 0, sd = sd_fcohort)
devs_fcohort_past_future31 <- c(fcohort_estimate_values, devs_fcohort31)
estimate_fcohort_all31 <- data.frame(year = years_cohort, estimate = devs_fcohort_past_future31)
result_matrix31 <- matrix(0, nrow = 17, ncol = num_cols)
result_matrix31[1, ] <- format(as.integer(estimate_fcohort_all31$year), nsmall = 0)
result_matrix31[2, ] <- estimate_fcohort_all31$estimate
print(result_matrix31)
num_rows31 <- nrow(result_matrix31)
num_cols31 <- ncol(result_matrix31)
value_to_repeat31 <- result_matrix31[2, ]
start_row31 <- 3
start_col31 <- 2
end_row31 <- num_rows31
end_col31 <- num_cols31
for (i in start_row31:end_row31) {
  for (j in start_col31:end_col31) {
    if (j >= (i - start_row31 + start_col31)) {
      pattern_index31 <- (j - (i - start_row31 + start_col31)) %% 32 + 1
      result_matrix31[i, j] <- value_to_repeat31[pattern_index31]
    }
  }
}
print(result_matrix31)
extracted_values31 <- as.numeric(result_matrix31[2:17, 16:33])
extracted_values31 <- extracted_values31[-1]
devs_fcohort_all31 <- c(devs_fcohort_past, extracted_values31)
devs_fyear_all31 <- c(devs_fyear_past, devs_fyear_fut31)
ewaa_m3$fyear_devs31 <- devs_fyear_all31
ewaa_m3$fcohort_devs31 <- devs_fcohort_all31
ewaa_m3$ewaa_fyear31 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fyear_devs31)
ewaa_m3$ewaa_fcohort31 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fcohort_devs31)
ewaa_m3$ewaa_fyear_fcoh31 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fyear_devs31 + ewaa_m3$fcohort_devs31)
head(ewaa_m3, n = 800)

# selected the predict weight with fyear and cohort effect only.
ewaa_m3_sim31 <- ewaa_m3[, c("real_year", "age", "ewaa_fyear_fcoh31")]
head(ewaa_m3_sim31,n=800)

#rename the columns as SS3 format ewaa.ss
ewaa_m3_sim31 <- ewaa_m3_sim31 %>%
  rename(year = real_year, pred_weight = ewaa_fyear_fcoh31)
colnames(ewaa_m3_sim31) <- c("year", "age", "pred_weight")
head(ewaa_m3_sim31,n=800)

# Convert the year column to numeric if necessary
ewaa_m3_sim31$year <- as.numeric(ewaa_m3_sim31$year)

# Run the pad_weight_at_age function and inspect the output
padded_data31 <- pad_weight_at_age(ewaa_m3_sim31)
str(padded_data31)

hakedataUSA:::write_wtatage_file(file="ewaa_sim31.ss",data=pad_weight_at_age(ewaa_m3_sim31),maturity = maturity_at_age)

#############################################################################
# Sim 32
set.seed(3001)
devs_fyear32 <- rnorm(nT, mean = 0, sd = sd_fyear)
devs_fyear_fut32 <- rep(devs_fyear32, each = 16)
devs_fcohort32 <- rnorm(nT_23, mean = 0, sd = sd_fcohort)
devs_fcohort_past_future32 <- c(fcohort_estimate_values, devs_fcohort32)
estimate_fcohort_all32 <- data.frame(year = years_cohort, estimate = devs_fcohort_past_future32)
result_matrix32 <- matrix(0, nrow = 17, ncol = num_cols)
result_matrix32[1, ] <- format(as.integer(estimate_fcohort_all32$year), nsmall = 0)
result_matrix32[2, ] <- estimate_fcohort_all32$estimate
print(result_matrix32)
num_rows32 <- nrow(result_matrix32)
num_cols32 <- ncol(result_matrix32)
value_to_repeat32 <- result_matrix32[2, ]
start_row32 <- 3
start_col32 <- 2
end_row32 <- num_rows32
end_col32 <- num_cols32
for (i in start_row32:end_row32) {
  for (j in start_col32:end_col32) {
    if (j >= (i - start_row32 + start_col32)) {
      pattern_index32 <- (j - (i - start_row32 + start_col32)) %% 32 + 1
      result_matrix32[i, j] <- value_to_repeat32[pattern_index32]
    }
  }
}
print(result_matrix32)
extracted_values32 <- as.numeric(result_matrix32[2:17, 16:33])
extracted_values32 <- extracted_values32[-1]
devs_fcohort_all32 <- c(devs_fcohort_past, extracted_values32)
devs_fyear_all32 <- c(devs_fyear_past, devs_fyear_fut32)
ewaa_m3$fyear_devs32 <- devs_fyear_all32
ewaa_m3$fcohort_devs32 <- devs_fcohort_all32
ewaa_m3$ewaa_fyear32 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fyear_devs32)
ewaa_m3$ewaa_fcohort32 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fcohort_devs32)
ewaa_m3$ewaa_fyear_fcoh32 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fyear_devs32 + ewaa_m3$fcohort_devs32)
head(ewaa_m3, n = 800)

# selected the predict weight with fyear and cohort effect only.
ewaa_m3_sim32 <- ewaa_m3[, c("real_year", "age", "ewaa_fyear_fcoh32")]
head(ewaa_m3_sim32,n=800)

#rename the columns as SS3 format ewaa.ss
ewaa_m3_sim32 <- ewaa_m3_sim32 %>%
  rename(year = real_year, pred_weight = ewaa_fyear_fcoh32)
colnames(ewaa_m3_sim32) <- c("year", "age", "pred_weight")
head(ewaa_m3_sim32,n=800)

# Convert the year column to numeric if necessary
ewaa_m3_sim32$year <- as.numeric(ewaa_m3_sim32$year)

# Run the pad_weight_at_age function and inspect the output
padded_data32 <- pad_weight_at_age(ewaa_m3_sim32)
str(padded_data32)

hakedataUSA:::write_wtatage_file(file="ewaa_sim32.ss",data=pad_weight_at_age(ewaa_m3_sim32),maturity = maturity_at_age)

#################################################################################
# Sim 33
set.seed(3141)
devs_fyear33 <- rnorm(nT, mean = 0, sd = sd_fyear)
devs_fyear_fut33 <- rep(devs_fyear33, each = 16)
devs_fcohort33 <- rnorm(nT_23, mean = 0, sd = sd_fcohort)
devs_fcohort_past_future33 <- c(fcohort_estimate_values, devs_fcohort33)
estimate_fcohort_all33 <- data.frame(year = years_cohort, estimate = devs_fcohort_past_future33)
result_matrix33 <- matrix(0, nrow = 17, ncol = num_cols)
result_matrix33[1, ] <- format(as.integer(estimate_fcohort_all33$year), nsmall = 0)
result_matrix33[2, ] <- estimate_fcohort_all33$estimate
print(result_matrix33)
num_rows33 <- nrow(result_matrix33)
num_cols33 <- ncol(result_matrix33)
value_to_repeat33 <- result_matrix33[2, ]
start_row33 <- 3
start_col33 <- 2
end_row33 <- num_rows33
end_col33 <- num_cols33
for (i in start_row33:end_row33) {
  for (j in start_col33:end_col33) {
    if (j >= (i - start_row33 + start_col33)) {
      pattern_index33 <- (j - (i - start_row33 + start_col33)) %% 32 + 1
      result_matrix33[i, j] <- value_to_repeat33[pattern_index33]
    }
  }
}
print(result_matrix33)
extracted_values33 <- as.numeric(result_matrix33[2:17, 16:33])
extracted_values33 <- extracted_values33[-1]
devs_fcohort_all33 <- c(devs_fcohort_past, extracted_values33)
devs_fyear_all33 <- c(devs_fyear_past, devs_fyear_fut33)
ewaa_m3$fyear_devs33 <- devs_fyear_all33
ewaa_m3$fcohort_devs33 <- devs_fcohort_all33
ewaa_m3$ewaa_fyear33 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fyear_devs33)
ewaa_m3$ewaa_fcohort33 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fcohort_devs33)
ewaa_m3$ewaa_fyear_fcoh33 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fyear_devs33 + ewaa_m3$fcohort_devs33)
head(ewaa_m3, n = 800)

# selected the predict weight with fyear and cohort effect only.
ewaa_m3_sim33 <- ewaa_m3[, c("real_year", "age", "ewaa_fyear_fcoh33")]
head(ewaa_m3_sim33,n=800)

#rename the columns as SS3 format ewaa.ss
ewaa_m3_sim33 <- ewaa_m3_sim33 %>%
  rename(year = real_year, pred_weight = ewaa_fyear_fcoh33)
colnames(ewaa_m3_sim33) <- c("year", "age", "pred_weight")
head(ewaa_m3_sim33,n=800)

# Convert the year column to numeric if necessary
ewaa_m3_sim33$year <- as.numeric(ewaa_m3_sim33$year)

# Run the pad_weight_at_age function and inspect the output
padded_data33 <- pad_weight_at_age(ewaa_m3_sim33)
str(padded_data33)

hakedataUSA:::write_wtatage_file(file="ewaa_sim33.ss",data=pad_weight_at_age(ewaa_m3_sim33),maturity = maturity_at_age)

#######################################################################
# Sim 34
set.seed(3313)
devs_fyear34 <- rnorm(nT, mean = 0, sd = sd_fyear)
devs_fyear_fut34 <- rep(devs_fyear34, each = 16)
devs_fcohort34 <- rnorm(nT_23, mean = 0, sd = sd_fcohort)
devs_fcohort_past_future34 <- c(fcohort_estimate_values, devs_fcohort34)
estimate_fcohort_all34 <- data.frame(year = years_cohort, estimate = devs_fcohort_past_future34)
result_matrix34 <- matrix(0, nrow = 17, ncol = num_cols)
result_matrix34[1, ] <- format(as.integer(estimate_fcohort_all34$year), nsmall = 0)
result_matrix34[2, ] <- estimate_fcohort_all34$estimate
print(result_matrix34)
num_rows34 <- nrow(result_matrix34)
num_cols34 <- ncol(result_matrix34)
value_to_repeat34 <- result_matrix34[2, ]
start_row34 <- 3
start_col34 <- 2
end_row34 <- num_rows34
end_col34 <- num_cols34
for (i in start_row34:end_row34) {
  for (j in start_col34:end_col34) {
    if (j >= (i - start_row34 + start_col34)) {
      pattern_index34 <- (j - (i - start_row34 + start_col34)) %% 32 + 1
      result_matrix34[i, j] <- value_to_repeat34[pattern_index34]
    }
  }
}
print(result_matrix34)
extracted_values34 <- as.numeric(result_matrix34[2:17, 16:33])
extracted_values34 <- extracted_values34[-1]
devs_fcohort_all34 <- c(devs_fcohort_past, extracted_values34)
devs_fyear_all34 <- c(devs_fyear_past, devs_fyear_fut34)
ewaa_m3$fyear_devs34 <- devs_fyear_all34
ewaa_m3$fcohort_devs34 <- devs_fcohort_all34
ewaa_m3$ewaa_fyear34 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fyear_devs34)
ewaa_m3$ewaa_fcohort34 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fcohort_devs34)
ewaa_m3$ewaa_fyear_fcoh34 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fyear_devs34 + ewaa_m3$fcohort_devs34)
head(ewaa_m3, n = 800)

# selected the predict weight with fyear and cohort effect only.
ewaa_m3_sim34 <- ewaa_m3[, c("real_year", "age", "ewaa_fyear_fcoh34")]
head(ewaa_m3_sim34,n=800)

#rename the columns as SS3 format ewaa.ss
ewaa_m3_sim34 <- ewaa_m3_sim34 %>%
  rename(year = real_year, pred_weight = ewaa_fyear_fcoh34)
colnames(ewaa_m3_sim34) <- c("year", "age", "pred_weight")
head(ewaa_m3_sim34,n=800)

# Convert the year column to numeric if necessary
ewaa_m3_sim34$year <- as.numeric(ewaa_m3_sim34$year)

# Run the pad_weight_at_age function and inspect the output
padded_data34 <- pad_weight_at_age(ewaa_m3_sim34)
str(padded_data34)

hakedataUSA:::write_wtatage_file(file="ewaa_sim34.ss",data=pad_weight_at_age(ewaa_m3_sim34),maturity = maturity_at_age)

###########################################################################################
# Sim 35
set.seed(3547)
devs_fyear35 <- rnorm(nT, mean = 0, sd = sd_fyear)
devs_fyear_fut35 <- rep(devs_fyear35, each = 16)
devs_fcohort35 <- rnorm(nT_23, mean = 0, sd = sd_fcohort)
devs_fcohort_past_future35 <- c(fcohort_estimate_values, devs_fcohort35)
estimate_fcohort_all35 <- data.frame(year = years_cohort, estimate = devs_fcohort_past_future35)
result_matrix35 <- matrix(0, nrow = 17, ncol = num_cols)
result_matrix35[1, ] <- format(as.integer(estimate_fcohort_all35$year), nsmall = 0)
result_matrix35[2, ] <- estimate_fcohort_all35$estimate
print(result_matrix35)
num_rows35 <- nrow(result_matrix35)
num_cols35 <- ncol(result_matrix35)
value_to_repeat35 <- result_matrix35[2, ]
start_row35 <- 3
start_col35 <- 2
end_row35 <- num_rows35
end_col35 <- num_cols35
for (i in start_row35:end_row35) {
  for (j in start_col35:end_col35) {
    if (j >= (i - start_row35 + start_col35)) {
      pattern_index35 <- (j - (i - start_row35 + start_col35)) %% 32 + 1
      result_matrix35[i, j] <- value_to_repeat35[pattern_index35]
    }
  }
}
print(result_matrix35)
extracted_values35 <- as.numeric(result_matrix35[2:17, 16:33])
extracted_values35 <- extracted_values35[-1]
devs_fcohort_all35 <- c(devs_fcohort_past, extracted_values35)
devs_fyear_all35 <- c(devs_fyear_past, devs_fyear_fut35)
ewaa_m3$fyear_devs35 <- devs_fyear_all35
ewaa_m3$fcohort_devs35 <- devs_fcohort_all35
ewaa_m3$ewaa_fyear35 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fyear_devs35)
ewaa_m3$ewaa_fcohort35 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fcohort_devs35)
ewaa_m3$ewaa_fyear_fcoh35 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fyear_devs35 + ewaa_m3$fcohort_devs35)
head(ewaa_m3, n = 800)

# selected the predict weight with fyear and cohort effect only.
ewaa_m3_sim35 <- ewaa_m3[, c("real_year", "age", "ewaa_fyear_fcoh35")]
head(ewaa_m3_sim35,n=800)

#rename the columns as SS3 format ewaa.ss
ewaa_m3_sim35 <- ewaa_m3_sim35 %>%
  rename(year = real_year, pred_weight = ewaa_fyear_fcoh35)
colnames(ewaa_m3_sim35) <- c("year", "age", "pred_weight")
head(ewaa_m3_sim35,n=800)

# Convert the year column to numeric if necessary
ewaa_m3_sim35$year <- as.numeric(ewaa_m3_sim35$year)

# Run the pad_weight_at_age function and inspect the output
padded_data35 <- pad_weight_at_age(ewaa_m3_sim35)
str(padded_data35)

hakedataUSA:::write_wtatage_file(file="ewaa_sim35.ss",data=pad_weight_at_age(ewaa_m3_sim35),maturity = maturity_at_age)

##############################################################################
# Sim 36
set.seed(3691)
devs_fyear36 <- rnorm(nT, mean = 0, sd = sd_fyear)
devs_fyear_fut36 <- rep(devs_fyear36, each = 16)
devs_fcohort36 <- rnorm(nT_23, mean = 0, sd = sd_fcohort)
devs_fcohort_past_future36 <- c(fcohort_estimate_values, devs_fcohort36)
estimate_fcohort_all36 <- data.frame(year = years_cohort, estimate = devs_fcohort_past_future36)
result_matrix36 <- matrix(0, nrow = 17, ncol = num_cols)
result_matrix36[1, ] <- format(as.integer(estimate_fcohort_all36$year), nsmall = 0)
result_matrix36[2, ] <- estimate_fcohort_all36$estimate
print(result_matrix36)
num_rows36 <- nrow(result_matrix36)
num_cols36 <- ncol(result_matrix36)
value_to_repeat36 <- result_matrix36[2, ]
start_row36 <- 3
start_col36 <- 2
end_row36 <- num_rows36
end_col36 <- num_cols36
for (i in start_row36:end_row36) {
  for (j in start_col36:end_col36) {
    if (j >= (i - start_row36 + start_col36)) {
      pattern_index36 <- (j - (i - start_row36 + start_col36)) %% 32 + 1
      result_matrix36[i, j] <- value_to_repeat36[pattern_index36]
    }
  }
}
print(result_matrix36)
extracted_values36 <- as.numeric(result_matrix36[2:17, 16:33])
extracted_values36 <- extracted_values36[-1]
devs_fcohort_all36 <- c(devs_fcohort_past, extracted_values36)
devs_fyear_all36 <- c(devs_fyear_past, devs_fyear_fut36)
ewaa_m3$fyear_devs36 <- devs_fyear_all36
ewaa_m3$fcohort_devs36 <- devs_fcohort_all36
ewaa_m3$ewaa_fyear36 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fyear_devs36)
ewaa_m3$ewaa_fcohort36 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fcohort_devs36)
ewaa_m3$ewaa_fyear_fcoh36 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fyear_devs36 + ewaa_m3$fcohort_devs36)
head(ewaa_m3, n = 800)

# selected the predict weight with fyear and cohort effect only.
ewaa_m3_sim36 <- ewaa_m3[, c("real_year", "age", "ewaa_fyear_fcoh36")]
head(ewaa_m3_sim36,n=800)

#rename the columns as SS3 format ewaa.ss
ewaa_m3_sim36 <- ewaa_m3_sim36 %>%
  rename(year = real_year, pred_weight = ewaa_fyear_fcoh36)
colnames(ewaa_m3_sim36) <- c("year", "age", "pred_weight")
head(ewaa_m3_sim36,n=800)

# Convert the year column to numeric if necessary
ewaa_m3_sim36$year <- as.numeric(ewaa_m3_sim36$year)

# Run the pad_weight_at_age function and inspect the output
padded_data36 <- pad_weight_at_age(ewaa_m3_sim36)
str(padded_data36)

hakedataUSA:::write_wtatage_file(file="ewaa_sim36.ss",data=pad_weight_at_age(ewaa_m3_sim36),maturity = maturity_at_age)

##############################################################################

# Sim 37
set.seed(3733)
devs_fyear37 <- rnorm(nT, mean = 0, sd = sd_fyear)
devs_fyear_fut37 <- rep(devs_fyear37, each = 16)
devs_fcohort37 <- rnorm(nT_23, mean = 0, sd = sd_fcohort)
devs_fcohort_past_future37 <- c(fcohort_estimate_values, devs_fcohort37)
estimate_fcohort_all37 <- data.frame(year = years_cohort, estimate = devs_fcohort_past_future37)
result_matrix37 <- matrix(0, nrow = 17, ncol = num_cols)
result_matrix37[1, ] <- format(as.integer(estimate_fcohort_all37$year), nsmall = 0)
result_matrix37[2, ] <- estimate_fcohort_all37$estimate
print(result_matrix37)
num_rows37 <- nrow(result_matrix37)
num_cols37 <- ncol(result_matrix37)
value_to_repeat37 <- result_matrix37[2, ]
start_row37 <- 3
start_col37 <- 2
end_row37 <- num_rows37
end_col37 <- num_cols37
for (i in start_row37:end_row37) {
  for (j in start_col37:end_col37) {
    if (j >= (i - start_row37 + start_col37)) {
      pattern_index37 <- (j - (i - start_row37 + start_col37)) %% 32 + 1
      result_matrix37[i, j] <- value_to_repeat37[pattern_index37]
    }
  }
}
print(result_matrix37)
extracted_values37 <- as.numeric(result_matrix37[2:17, 16:33])
extracted_values37 <- extracted_values37[-1]
devs_fcohort_all37 <- c(devs_fcohort_past, extracted_values37)
devs_fyear_all37 <- c(devs_fyear_past, devs_fyear_fut37)
ewaa_m3$fyear_devs37 <- devs_fyear_all37
ewaa_m3$fcohort_devs37 <- devs_fcohort_all37
ewaa_m3$ewaa_fyear37 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fyear_devs37)
ewaa_m3$ewaa_fcohort37 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fcohort_devs37)
ewaa_m3$ewaa_fyear_fcoh37 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fyear_devs37 + ewaa_m3$fcohort_devs37)
head(ewaa_m3, n = 800)

# selected the predict weight with fyear and cohort effect only.
ewaa_m3_sim37 <- ewaa_m3[, c("real_year", "age", "ewaa_fyear_fcoh37")]
head(ewaa_m3_sim37,n=800)

#rename the columns as SS3 format ewaa.ss
ewaa_m3_sim37 <- ewaa_m3_sim37 %>%
  rename(year = real_year, pred_weight = ewaa_fyear_fcoh37)
colnames(ewaa_m3_sim37) <- c("year", "age", "pred_weight")
head(ewaa_m3_sim37,n=800)

# Convert the year column to numeric if necessary
ewaa_m3_sim37$year <- as.numeric(ewaa_m3_sim37$year)

# Run the pad_weight_at_age function and inspect the output
padded_data37 <- pad_weight_at_age(ewaa_m3_sim37)
str(padded_data37)

hakedataUSA:::write_wtatage_file(file="ewaa_sim37.ss",data=pad_weight_at_age(ewaa_m3_sim37),maturity = maturity_at_age)


##############################################################################


# Sim 38
set.seed(3833)
devs_fyear38 <- rnorm(nT, mean = 0, sd = sd_fyear)
devs_fyear_fut38 <- rep(devs_fyear38, each = 16)
devs_fcohort38 <- rnorm(nT_23, mean = 0, sd = sd_fcohort)
devs_fcohort_past_future38 <- c(fcohort_estimate_values, devs_fcohort38)
estimate_fcohort_all38 <- data.frame(year = years_cohort, estimate = devs_fcohort_past_future38)
result_matrix38 <- matrix(0, nrow = 17, ncol = num_cols)
result_matrix38[1, ] <- format(as.integer(estimate_fcohort_all38$year), nsmall = 0)
result_matrix38[2, ] <- estimate_fcohort_all38$estimate
print(result_matrix38)
num_rows38 <- nrow(result_matrix38)
num_cols38 <- ncol(result_matrix38)
value_to_repeat38 <- result_matrix38[2, ]
start_row38 <- 3
start_col38 <- 2
end_row38 <- num_rows38
end_col38 <- num_cols38
for (i in start_row38:end_row38) {
  for (j in start_col38:end_col38) {
    if (j >= (i - start_row38 + start_col38)) {
      pattern_index38 <- (j - (i - start_row38 + start_col38)) %% 32 + 1
      result_matrix38[i, j] <- value_to_repeat38[pattern_index38]
    }
  }
}
print(result_matrix38)
extracted_values38 <- as.numeric(result_matrix38[2:17, 16:33])
extracted_values38 <- extracted_values38[-1]
devs_fcohort_all38 <- c(devs_fcohort_past, extracted_values38)
devs_fyear_all38 <- c(devs_fyear_past, devs_fyear_fut38)
ewaa_m3$fyear_devs38 <- devs_fyear_all38
ewaa_m3$fcohort_devs38 <- devs_fcohort_all38
ewaa_m3$ewaa_fyear38 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fyear_devs38)
ewaa_m3$ewaa_fcohort38 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fcohort_devs38)
ewaa_m3$ewaa_fyear_fcoh38 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fyear_devs38 + ewaa_m3$fcohort_devs38)
head(ewaa_m3, n = 800)

# selected the predict weight with fyear and cohort effect only.
ewaa_m3_sim38 <- ewaa_m3[, c("real_year", "age", "ewaa_fyear_fcoh38")]
head(ewaa_m3_sim38,n=800)

#rename the columns as SS3 format ewaa.ss
ewaa_m3_sim38 <- ewaa_m3_sim38 %>%
  rename(year = real_year, pred_weight = ewaa_fyear_fcoh38)
colnames(ewaa_m3_sim38) <- c("year", "age", "pred_weight")
head(ewaa_m3_sim38,n=800)

# Convert the year column to numeric if necessary
ewaa_m3_sim38$year <- as.numeric(ewaa_m3_sim38$year)

# Run the pad_weight_at_age function and inspect the output
padded_data38 <- pad_weight_at_age(ewaa_m3_sim38)
str(padded_data38)

hakedataUSA:::write_wtatage_file(file="ewaa_sim38.ss",data=pad_weight_at_age(ewaa_m3_sim38),maturity = maturity_at_age)
#######################################################################################################################

# Sim 39
set.seed(3967)
devs_fyear39 <- rnorm(nT, mean = 0, sd = sd_fyear)
devs_fyear_fut39 <- rep(devs_fyear39, each = 16)
devs_fcohort39 <- rnorm(nT_23, mean = 0, sd = sd_fcohort)
devs_fcohort_past_future39 <- c(fcohort_estimate_values, devs_fcohort39)
estimate_fcohort_all39 <- data.frame(year = years_cohort, estimate = devs_fcohort_past_future39)
result_matrix39 <- matrix(0, nrow = 17, ncol = num_cols)
result_matrix39[1, ] <- format(as.integer(estimate_fcohort_all39$year), nsmall = 0)
result_matrix39[2, ] <- estimate_fcohort_all39$estimate
print(result_matrix39)
num_rows39 <- nrow(result_matrix39)
num_cols39 <- ncol(result_matrix39)
value_to_repeat39 <- result_matrix39[2, ]
start_row39 <- 3
start_col39 <- 2
end_row39 <- num_rows39
end_col39 <- num_cols39
for (i in start_row39:end_row39) {
  for (j in start_col39:end_col39) {
    if (j >= (i - start_row39 + start_col39)) {
      pattern_index39 <- (j - (i - start_row39 + start_col39)) %% 32 + 1
      result_matrix39[i, j] <- value_to_repeat39[pattern_index39]
    }
  }
}
print(result_matrix39)
extracted_values39 <- as.numeric(result_matrix39[2:17, 16:33])
extracted_values39 <- extracted_values39[-1]
devs_fcohort_all39 <- c(devs_fcohort_past, extracted_values39)
devs_fyear_all39 <- c(devs_fyear_past, devs_fyear_fut39)
ewaa_m3$fyear_devs39 <- devs_fyear_all39
ewaa_m3$fcohort_devs39 <- devs_fcohort_all39
ewaa_m3$ewaa_fyear39 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fyear_devs39)
ewaa_m3$ewaa_fcohort39 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fcohort_devs39)
ewaa_m3$ewaa_fyear_fcoh39 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fyear_devs39 + ewaa_m3$fcohort_devs39)
head(ewaa_m3, n = 800)

# selected the predict weight with fyear and cohort effect only.
ewaa_m3_sim39 <- ewaa_m3[, c("real_year", "age", "ewaa_fyear_fcoh39")]
head(ewaa_m3_sim39,n=800)

#rename the columns as SS3 format ewaa.ss
ewaa_m3_sim39 <- ewaa_m3_sim39 %>%
  rename(year = real_year, pred_weight = ewaa_fyear_fcoh39)
colnames(ewaa_m3_sim39) <- c("year", "age", "pred_weight")
head(ewaa_m3_sim39,n=800)

# Convert the year column to numeric if necessary
ewaa_m3_sim39$year <- as.numeric(ewaa_m3_sim39$year)

# Run the pad_weight_at_age function and inspect the output
padded_data39 <- pad_weight_at_age(ewaa_m3_sim39)
str(padded_data39)

hakedataUSA:::write_wtatage_file(file="ewaa_sim39.ss",data=pad_weight_at_age(ewaa_m3_sim39),maturity = maturity_at_age)
#######################################################################################################################

# Sim 40
set.seed(4032)
devs_fyear40 <- rnorm(nT, mean = 0, sd = sd_fyear)
devs_fyear_fut40 <- rep(devs_fyear40, each = 16)
devs_fcohort40 <- rnorm(nT_23, mean = 0, sd = sd_fcohort)
devs_fcohort_past_future40 <- c(fcohort_estimate_values, devs_fcohort40)
estimate_fcohort_all40 <- data.frame(year = years_cohort, estimate = devs_fcohort_past_future40)
result_matrix40 <- matrix(0, nrow = 17, ncol = num_cols)
result_matrix40[1, ] <- format(as.integer(estimate_fcohort_all40$year), nsmall = 0)
result_matrix40[2, ] <- estimate_fcohort_all40$estimate
print(result_matrix40)
num_rows40 <- nrow(result_matrix40)
num_cols40 <- ncol(result_matrix40)
value_to_repeat40 <- result_matrix40[2, ]
start_row40 <- 3
start_col40 <- 2
end_row40 <- num_rows40
end_col40 <- num_cols40
for (i in start_row40:end_row40) {
  for (j in start_col40:end_col40) {
    if (j >= (i - start_row40 + start_col40)) {
      pattern_index40 <- (j - (i - start_row40 + start_col40)) %% 32 + 1
      result_matrix40[i, j] <- value_to_repeat40[pattern_index40]
    }
  }
}
print(result_matrix40)
extracted_values40 <- as.numeric(result_matrix40[2:17, 16:33])
extracted_values40 <- extracted_values40[-1]
devs_fcohort_all40 <- c(devs_fcohort_past, extracted_values40)
devs_fyear_all40 <- c(devs_fyear_past, devs_fyear_fut40)
ewaa_m3$fyear_devs40 <- devs_fyear_all40
ewaa_m3$fcohort_devs40 <- devs_fcohort_all40
ewaa_m3$ewaa_fyear40 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fyear_devs40)
ewaa_m3$ewaa_fcohort40 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fcohort_devs40)
ewaa_m3$ewaa_fyear_fcoh40 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fyear_devs40 + ewaa_m3$fcohort_devs40)
head(ewaa_m3, n = 10)

# selected the predict weight with fyear and cohort effect only.
ewaa_m3_sim40 <- ewaa_m3[, c("real_year", "age", "ewaa_fyear_fcoh40")]
head(ewaa_m3_sim40,n=800)

#rename the columns as SS3 format ewaa.ss
ewaa_m3_sim40 <- ewaa_m3_sim40 %>%
  rename(year = real_year, pred_weight = ewaa_fyear_fcoh40)
colnames(ewaa_m3_sim40) <- c("year", "age", "pred_weight")
head(ewaa_m3_sim40,n=800)

# Convert the year column to numeric if necessary
ewaa_m3_sim40$year <- as.numeric(ewaa_m3_sim40$year)

# Run the pad_weight_at_age function and inspect the output
padded_data40 <- pad_weight_at_age(ewaa_m3_sim40)
str(padded_data40)

hakedataUSA:::write_wtatage_file(file="ewaa_sim40.ss",data=pad_weight_at_age(ewaa_m3_sim40),maturity = maturity_at_age)
#######################################################################################################################

# Sim 41
set.seed(4113)
devs_fyear41 <- rnorm(nT, mean = 0, sd = sd_fyear)
devs_fyear_fut41 <- rep(devs_fyear41, each = 16)
devs_fcohort41 <- rnorm(nT_23, mean = 0, sd = sd_fcohort)
devs_fcohort_past_future41 <- c(fcohort_estimate_values, devs_fcohort41)
estimate_fcohort_all41 <- data.frame(year = years_cohort, estimate = devs_fcohort_past_future41)
result_matrix41 <- matrix(0, nrow = 17, ncol = num_cols)
result_matrix41[1, ] <- format(as.integer(estimate_fcohort_all41$year), nsmall = 0)
result_matrix41[2, ] <- estimate_fcohort_all41$estimate
print(result_matrix41)
num_rows41 <- nrow(result_matrix41)
num_cols41 <- ncol(result_matrix41)
value_to_repeat41 <- result_matrix41[2, ]
start_row41 <- 3
start_col41 <- 2
end_row41 <- num_rows41
end_col41 <- num_cols41
for (i in start_row41:end_row41) {
  for (j in start_col41:end_col41) {
    if (j >= (i - start_row41 + start_col41)) {
      pattern_index41 <- (j - (i - start_row41 + start_col41)) %% 32 + 1
      result_matrix41[i, j] <- value_to_repeat41[pattern_index41]
    }
  }
}
print(result_matrix41)
extracted_values41 <- as.numeric(result_matrix41[2:17, 16:33])
extracted_values41 <- extracted_values41[-1]
devs_fcohort_all41 <- c(devs_fcohort_past, extracted_values41)
devs_fyear_all41 <- c(devs_fyear_past, devs_fyear_fut41)
ewaa_m3$fyear_devs41 <- devs_fyear_all41
ewaa_m3$fcohort_devs41 <- devs_fcohort_all41
ewaa_m3$ewaa_fyear41 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fyear_devs41)
ewaa_m3$ewaa_fcohort41 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fcohort_devs41)
ewaa_m3$ewaa_fyear_fcoh41 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fyear_devs41 + ewaa_m3$fcohort_devs41)
head(ewaa_m3, n = 800)
# selected the predict weight with fyear and cohort effect only.
ewaa_m3_sim41 <- ewaa_m3[, c("real_year", "age", "ewaa_fyear_fcoh41")]
head(ewaa_m3_sim41,n=800)

#rename the columns as SS3 format ewaa.ss
ewaa_m3_sim41 <- ewaa_m3_sim41 %>%
  rename(year = real_year, pred_weight = ewaa_fyear_fcoh41)
colnames(ewaa_m3_sim41) <- c("year", "age", "pred_weight")
head(ewaa_m3_sim41,n=800)

# Convert the year column to numeric if necessary
ewaa_m3_sim41$year <- as.numeric(ewaa_m3_sim41$year)

# Run the pad_weight_at_age function and inspect the output
padded_data41 <- pad_weight_at_age(ewaa_m3_sim41)
str(padded_data41)

hakedataUSA:::write_wtatage_file(file="ewaa_sim41.ss",data=pad_weight_at_age(ewaa_m3_sim41),maturity = maturity_at_age)
#######################################################################################################################

# Sim 42
set.seed(4211)
devs_fyear42 <- rnorm(nT, mean = 0, sd = sd_fyear)
devs_fyear_fut42 <- rep(devs_fyear42, each = 16)
devs_fcohort42 <- rnorm(nT_23, mean = 0, sd = sd_fcohort)
devs_fcohort_past_future42 <- c(fcohort_estimate_values, devs_fcohort42)
estimate_fcohort_all42 <- data.frame(year = years_cohort, estimate = devs_fcohort_past_future42)
result_matrix42 <- matrix(0, nrow = 17, ncol = num_cols)
result_matrix42[1, ] <- format(as.integer(estimate_fcohort_all42$year), nsmall = 0)
result_matrix42[2, ] <- estimate_fcohort_all42$estimate
print(result_matrix42)
num_rows42 <- nrow(result_matrix42)
num_cols42 <- ncol(result_matrix42)
value_to_repeat42 <- result_matrix42[2, ]
start_row42 <- 3
start_col42 <- 2
end_row42 <- num_rows42
end_col42 <- num_cols42
for (i in start_row42:end_row42) {
  for (j in start_col42:end_col42) {
    if (j >= (i - start_row42 + start_col42)) {
      pattern_index42 <- (j - (i - start_row42 + start_col42)) %% 32 + 1
      result_matrix42[i, j] <- value_to_repeat42[pattern_index42]
    }
  }
}
print(result_matrix42)
extracted_values42 <- as.numeric(result_matrix42[2:17, 16:33])
extracted_values42 <- extracted_values42[-1]
devs_fcohort_all42 <- c(devs_fcohort_past, extracted_values42)
devs_fyear_all42 <- c(devs_fyear_past, devs_fyear_fut42)
ewaa_m3$fyear_devs42 <- devs_fyear_all42
ewaa_m3$fcohort_devs42 <- devs_fcohort_all42
ewaa_m3$ewaa_fyear42 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fyear_devs42)
ewaa_m3$ewaa_fcohort42 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fcohort_devs42)
ewaa_m3$ewaa_fyear_fcoh42 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fyear_devs42 + ewaa_m3$fcohort_devs42)
head(ewaa_m3, n = 800)

# selected the predict weight with fyear and cohort effect only.
ewaa_m3_sim42 <- ewaa_m3[, c("real_year", "age", "ewaa_fyear_fcoh42")]
head(ewaa_m3_sim42,n=800)

#rename the columns as SS3 format ewaa.ss
ewaa_m3_sim42 <- ewaa_m3_sim42 %>%
  rename(year = real_year, pred_weight = ewaa_fyear_fcoh42)
colnames(ewaa_m3_sim42) <- c("year", "age", "pred_weight")
head(ewaa_m3_sim42,n=800)

# Convert the year column to numeric if necessary
ewaa_m3_sim42$year <- as.numeric(ewaa_m3_sim42$year)

# Run the pad_weight_at_age function and inspect the output
padded_data42 <- pad_weight_at_age(ewaa_m3_sim42)
str(padded_data42)

hakedataUSA:::write_wtatage_file(file="ewaa_sim42.ss",data=pad_weight_at_age(ewaa_m3_sim42),maturity = maturity_at_age)
#######################################################################################################################

# Sim 43
set.seed(4307)
devs_fyear43 <- rnorm(nT, mean = 0, sd = sd_fyear)
devs_fyear_fut43 <- rep(devs_fyear43, each = 16)
devs_fcohort43 <- rnorm(nT_23, mean = 0, sd = sd_fcohort)
devs_fcohort_past_future43 <- c(fcohort_estimate_values, devs_fcohort43)
estimate_fcohort_all43 <- data.frame(year = years_cohort, estimate = devs_fcohort_past_future43)
result_matrix43 <- matrix(0, nrow = 17, ncol = num_cols)
result_matrix43[1, ] <- format(as.integer(estimate_fcohort_all43$year), nsmall = 0)
result_matrix43[2, ] <- estimate_fcohort_all43$estimate
print(result_matrix43)
num_rows43 <- nrow(result_matrix43)
num_cols43 <- ncol(result_matrix43)
value_to_repeat43 <- result_matrix43[2, ]
start_row43 <- 3
start_col43 <- 2
end_row43 <- num_rows43
end_col43 <- num_cols43
for (i in start_row43:end_row43) {
  for (j in start_col43:end_col43) {
    if (j >= (i - start_row43 + start_col43)) {
      pattern_index43 <- (j - (i - start_row43 + start_col43)) %% 32 + 1
      result_matrix43[i, j] <- value_to_repeat43[pattern_index43]
    }
  }
}
print(result_matrix43)
extracted_values43 <- as.numeric(result_matrix43[2:17, 16:33])
extracted_values43 <- extracted_values43[-1]
devs_fcohort_all43 <- c(devs_fcohort_past, extracted_values43)
devs_fyear_all43 <- c(devs_fyear_past, devs_fyear_fut43)
ewaa_m3$fyear_devs43 <- devs_fyear_all43
ewaa_m3$fcohort_devs43 <- devs_fcohort_all43
ewaa_m3$ewaa_fyear43 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fyear_devs43)
ewaa_m3$ewaa_fcohort43 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fcohort_devs43)
ewaa_m3$ewaa_fyear_fcoh43 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fyear_devs43 + ewaa_m3$fcohort_devs43)
head(ewaa_m3, n = 800)

# selected the predict weight with fyear and cohort effect only.
ewaa_m3_sim43 <- ewaa_m3[, c("real_year", "age", "ewaa_fyear_fcoh43")]
head(ewaa_m3_sim43, n=800)

# rename the columns as SS3 format ewaa.ss
ewaa_m3_sim43 <- ewaa_m3_sim43 %>%
  rename(year = real_year, pred_weight = ewaa_fyear_fcoh43)
colnames(ewaa_m3_sim43) <- c("year", "age", "pred_weight")
head(ewaa_m3_sim43, n=800)

# Convert the year column to numeric if necessary
ewaa_m3_sim43$year <- as.numeric(ewaa_m3_sim43$year)

# Run the pad_weight_at_age function and inspect the output
padded_data43 <- pad_weight_at_age(ewaa_m3_sim43)
str(padded_data43)

hakedataUSA:::write_wtatage_file(file="ewaa_sim43.ss", data=padded_data43, maturity = maturity_at_age)


#######################################################
# Sim 44
set.seed(4423)
devs_fyear44 <- rnorm(nT, mean = 0, sd = sd_fyear)
devs_fyear_fut44 <- rep(devs_fyear44, each = 16)
devs_fcohort44 <- rnorm(nT_23, mean = 0, sd = sd_fcohort)
devs_fcohort_past_future44 <- c(fcohort_estimate_values, devs_fcohort44)
estimate_fcohort_all44 <- data.frame(year = years_cohort, estimate = devs_fcohort_past_future44)
result_matrix44 <- matrix(0, nrow = 17, ncol = num_cols)
result_matrix44[1, ] <- format(as.integer(estimate_fcohort_all44$year), nsmall = 0)
result_matrix44[2, ] <- estimate_fcohort_all44$estimate
print(result_matrix44)
num_rows44 <- nrow(result_matrix44)
num_cols44 <- ncol(result_matrix44)
value_to_repeat44 <- result_matrix44[2, ]
start_row44 <- 3
start_col44 <- 2
end_row44 <- num_rows44
end_col44 <- num_cols44
for (i in start_row44:end_row44) {
  for (j in start_col44:end_col44) {
    if (j >= (i - start_row44 + start_col44)) {
      pattern_index44 <- (j - (i - start_row44 + start_col44)) %% 32 + 1
      result_matrix44[i, j] <- value_to_repeat44[pattern_index44]
    }
  }
}
print(result_matrix44)
extracted_values44 <- as.numeric(result_matrix44[2:17, 16:33])
extracted_values44 <- extracted_values44[-1]
devs_fcohort_all44 <- c(devs_fcohort_past, extracted_values44)
devs_fyear_all44 <- c(devs_fyear_past, devs_fyear_fut44)
ewaa_m3$fyear_devs44 <- devs_fyear_all44
ewaa_m3$fcohort_devs44 <- devs_fcohort_all44
ewaa_m3$ewaa_fyear44 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fyear_devs44)
ewaa_m3$ewaa_fcohort44 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fcohort_devs44)
ewaa_m3$ewaa_fyear_fcoh44 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fyear_devs44 + ewaa_m3$fcohort_devs44)
head(ewaa_m3, n = 10)


# selected the predict weight with fyear and cohort effect only.
ewaa_m3_sim44 <- ewaa_m3[, c("real_year", "age", "ewaa_fyear_fcoh44")]
head(ewaa_m3_sim44, n=800)

# rename the columns as SS3 format ewaa.ss
ewaa_m3_sim44 <- ewaa_m3_sim44 %>%
  rename(year = real_year, pred_weight = ewaa_fyear_fcoh44)
colnames(ewaa_m3_sim44) <- c("year", "age", "pred_weight")
head(ewaa_m3_sim44, n=800)

# Convert the year column to numeric if necessary
ewaa_m3_sim44$year <- as.numeric(ewaa_m3_sim44$year)

# Run the pad_weight_at_age function and inspect the output
padded_data44 <- pad_weight_at_age(ewaa_m3_sim44)
str(padded_data44)

hakedataUSA:::write_wtatage_file(file="ewaa_sim44.ss", data=padded_data44, maturity = maturity_at_age)

#################################################################
# Sim 45
set.seed(4531)
devs_fyear45 <- rnorm(nT, mean = 0, sd = sd_fyear)
devs_fyear_fut45 <- rep(devs_fyear45, each = 16)
devs_fcohort45 <- rnorm(nT_23, mean = 0, sd = sd_fcohort)
devs_fcohort_past_future45 <- c(fcohort_estimate_values, devs_fcohort45)
estimate_fcohort_all45 <- data.frame(year = years_cohort, estimate = devs_fcohort_past_future45)
result_matrix45 <- matrix(0, nrow = 17, ncol = num_cols)
result_matrix45[1, ] <- format(as.integer(estimate_fcohort_all45$year), nsmall = 0)
result_matrix45[2, ] <- estimate_fcohort_all45$estimate
print(result_matrix45)
num_rows45 <- nrow(result_matrix45)
num_cols45 <- ncol(result_matrix45)
value_to_repeat45 <- result_matrix45[2, ]
start_row45 <- 3
start_col45 <- 2
end_row45 <- num_rows45
end_col45 <- num_cols45
for (i in start_row45:end_row45) {
  for (j in start_col45:end_col45) {
    if (j >= (i - start_row45 + start_col45)) {
      pattern_index45 <- (j - (i - start_row45 + start_col45)) %% 32 + 1
      result_matrix45[i, j] <- value_to_repeat45[pattern_index45]
    }
  }
}
print(result_matrix45)
extracted_values45 <- as.numeric(result_matrix45[2:17, 16:33])
extracted_values45 <- extracted_values45[-1]
devs_fcohort_all45 <- c(devs_fcohort_past, extracted_values45)
devs_fyear_all45 <- c(devs_fyear_past, devs_fyear_fut45)
ewaa_m3$fyear_devs45 <- devs_fyear_all45
ewaa_m3$fcohort_devs45 <- devs_fcohort_all45
ewaa_m3$ewaa_fyear45 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fyear_devs45)
ewaa_m3$ewaa_fcohort45 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fcohort_devs45)
ewaa_m3$ewaa_fyear_fcoh45 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fyear_devs45 + ewaa_m3$fcohort_devs45)
head(ewaa_m3, n = 10)

# selected the predict weight with fyear and cohort effect only.
ewaa_m3_sim45 <- ewaa_m3[, c("real_year", "age", "ewaa_fyear_fcoh45")]
head(ewaa_m3_sim45, n=800)

# rename the columns as SS3 format ewaa.ss
ewaa_m3_sim45 <- ewaa_m3_sim45 %>%
  rename(year = real_year, pred_weight = ewaa_fyear_fcoh45)
colnames(ewaa_m3_sim45) <- c("year", "age", "pred_weight")
head(ewaa_m3_sim45, n=800)

# Convert the year column to numeric if necessary
ewaa_m3_sim45$year <- as.numeric(ewaa_m3_sim45$year)

# Run the pad_weight_at_age function and inspect the output
padded_data45 <- pad_weight_at_age(ewaa_m3_sim45)
str(padded_data45)

hakedataUSA:::write_wtatage_file(file="ewaa_sim45.ss", data=padded_data45, maturity = maturity_at_age)

##################################################################################################
# Sim 46
set.seed(4612)
devs_fyear46 <- rnorm(nT, mean = 0, sd = sd_fyear)
devs_fyear_fut46 <- rep(devs_fyear46, each = 16)
devs_fcohort46 <- rnorm(nT_23, mean = 0, sd = sd_fcohort)
devs_fcohort_past_future46 <- c(fcohort_estimate_values, devs_fcohort46)
estimate_fcohort_all46 <- data.frame(year = years_cohort, estimate = devs_fcohort_past_future46)
result_matrix46 <- matrix(0, nrow = 17, ncol = num_cols)
result_matrix46[1, ] <- format(as.integer(estimate_fcohort_all46$year), nsmall = 0)
result_matrix46[2, ] <- estimate_fcohort_all46$estimate
print(result_matrix46)
num_rows46 <- nrow(result_matrix46)
num_cols46 <- ncol(result_matrix46)
value_to_repeat46 <- result_matrix46[2, ]
start_row46 <- 3
start_col46 <- 2
end_row46 <- num_rows46
end_col46 <- num_cols46
for (i in start_row46:end_row46) {
  for (j in start_col46:end_col46) {
    if (j >= (i - start_row46 + start_col46)) {
      pattern_index46 <- (j - (i - start_row46 + start_col46)) %% 32 + 1
      result_matrix46[i, j] <- value_to_repeat46[pattern_index46]
    }
  }
}
print(result_matrix46)
extracted_values46 <- as.numeric(result_matrix46[2:17, 16:33])
extracted_values46 <- extracted_values46[-1]
devs_fcohort_all46 <- c(devs_fcohort_past, extracted_values46)
devs_fyear_all46 <- c(devs_fyear_past, devs_fyear_fut46)
ewaa_m3$fyear_devs46 <- devs_fyear_all46
ewaa_m3$fcohort_devs46 <- devs_fcohort_all46
ewaa_m3$ewaa_fyear46 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fyear_devs46)
ewaa_m3$ewaa_fcohort46 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fcohort_devs46)
ewaa_m3$ewaa_fyear_fcoh46 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fyear_devs46 + ewaa_m3$fcohort_devs46)
head(ewaa_m3, n = 800)

# selected the predict weight with fyear and cohort effect only.
ewaa_m3_sim46 <- ewaa_m3[, c("real_year", "age", "ewaa_fyear_fcoh46")]
head(ewaa_m3_sim46, n=800)

# rename the columns as SS3 format ewaa.ss
ewaa_m3_sim46 <- ewaa_m3_sim46 %>%
  rename(year = real_year, pred_weight = ewaa_fyear_fcoh46)
colnames(ewaa_m3_sim46) <- c("year", "age", "pred_weight")
head(ewaa_m3_sim46, n=800)

# Convert the year column to numeric if necessary
ewaa_m3_sim46$year <- as.numeric(ewaa_m3_sim46$year)

# Run the pad_weight_at_age function and inspect the output
padded_data46 <- pad_weight_at_age(ewaa_m3_sim46)
str(padded_data46)

hakedataUSA:::write_wtatage_file(file="ewaa_sim46.ss", data=padded_data46, maturity = maturity_at_age)

########################################################################################################
# Sim 47
set.seed(4723)
devs_fyear47 <- rnorm(nT, mean = 0, sd = sd_fyear)
devs_fyear_fut47 <- rep(devs_fyear47, each = 16)
devs_fcohort47 <- rnorm(nT_23, mean = 0, sd = sd_fcohort)
devs_fcohort_past_future47 <- c(fcohort_estimate_values, devs_fcohort47)
estimate_fcohort_all47 <- data.frame(year = years_cohort, estimate = devs_fcohort_past_future47)
result_matrix47 <- matrix(0, nrow = 17, ncol = num_cols)
result_matrix47[1, ] <- format(as.integer(estimate_fcohort_all47$year), nsmall = 0)
result_matrix47[2, ] <- estimate_fcohort_all47$estimate
print(result_matrix47)
num_rows47 <- nrow(result_matrix47)
num_cols47 <- ncol(result_matrix47)
value_to_repeat47 <- result_matrix47[2, ]
start_row47 <- 3
start_col47 <- 2
end_row47 <- num_rows47
end_col47 <- num_cols47
for (i in start_row47:end_row47) {
  for (j in start_col47:end_col47) {
    if (j >= (i - start_row47 + start_col47)) {
      pattern_index47 <- (j - (i - start_row47 + start_col47)) %% 32 + 1
      result_matrix47[i, j] <- value_to_repeat47[pattern_index47]
    }
  }
}
print(result_matrix47)
extracted_values47 <- as.numeric(result_matrix47[2:17, 16:33])
extracted_values47 <- extracted_values47[-1]
devs_fcohort_all47 <- c(devs_fcohort_past, extracted_values47)
devs_fyear_all47 <- c(devs_fyear_past, devs_fyear_fut47)
ewaa_m3$fyear_devs47 <- devs_fyear_all47
ewaa_m3$fcohort_devs47 <- devs_fcohort_all47
ewaa_m3$ewaa_fyear47 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fyear_devs47)
ewaa_m3$ewaa_fcohort47 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fcohort_devs47)
ewaa_m3$ewaa_fyear_fcoh47 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fyear_devs47 + ewaa_m3$fcohort_devs47)
head(ewaa_m3, n = 10)

# selected the predict weight with fyear and cohort effect only.
ewaa_m3_sim47 <- ewaa_m3[, c("real_year", "age", "ewaa_fyear_fcoh47")]
head(ewaa_m3_sim47, n=800)

# rename the columns as SS3 format ewaa.ss
ewaa_m3_sim47 <- ewaa_m3_sim47 %>%
  rename(year = real_year, pred_weight = ewaa_fyear_fcoh47)
colnames(ewaa_m3_sim47) <- c("year", "age", "pred_weight")
head(ewaa_m3_sim47, n=800)

# Convert the year column to numeric if necessary
ewaa_m3_sim47$year <- as.numeric(ewaa_m3_sim47$year)

# Run the pad_weight_at_age function and inspect the output
padded_data47 <- pad_weight_at_age(ewaa_m3_sim47)
str(padded_data47)

hakedataUSA:::write_wtatage_file(file="ewaa_sim47.ss", data=padded_data47, maturity = maturity_at_age)
#####################################################################################################

# Sim 48
set.seed(4819)
devs_fyear48 <- rnorm(nT, mean = 0, sd = sd_fyear)
devs_fyear_fut48 <- rep(devs_fyear48, each = 16)
devs_fcohort48 <- rnorm(nT_23, mean = 0, sd = sd_fcohort)
devs_fcohort_past_future48 <- c(fcohort_estimate_values, devs_fcohort48)
estimate_fcohort_all48 <- data.frame(year = years_cohort, estimate = devs_fcohort_past_future48)
result_matrix48 <- matrix(0, nrow = 17, ncol = num_cols)
result_matrix48[1, ] <- format(as.integer(estimate_fcohort_all48$year), nsmall = 0)
result_matrix48[2, ] <- estimate_fcohort_all48$estimate
print(result_matrix48)
num_rows48 <- nrow(result_matrix48)
num_cols48 <- ncol(result_matrix48)
value_to_repeat48 <- result_matrix48[2, ]
start_row48 <- 3
start_col48 <- 2
end_row48 <- num_rows48
end_col48 <- num_cols48
for (i in start_row48:end_row48) {
  for (j in start_col48:end_col48) {
    if (j >= (i - start_row48 + start_col48)) {
      pattern_index48 <- (j - (i - start_row48 + start_col48)) %% 32 + 1
      result_matrix48[i, j] <- value_to_repeat48[pattern_index48]
    }
  }
}
print(result_matrix48)
extracted_values48 <- as.numeric(result_matrix48[2:17, 16:33])
extracted_values48 <- extracted_values48[-1]
devs_fcohort_all48 <- c(devs_fcohort_past, extracted_values48)
devs_fyear_all48 <- c(devs_fyear_past, devs_fyear_fut48)
ewaa_m3$fyear_devs48 <- devs_fyear_all48
ewaa_m3$fcohort_devs48 <- devs_fcohort_all48
ewaa_m3$ewaa_fyear48 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fyear_devs48)
ewaa_m3$ewaa_fcohort48 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fcohort_devs48)
ewaa_m3$ewaa_fyear_fcoh48 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fyear_devs48 + ewaa_m3$fcohort_devs48)
head(ewaa_m3, n = 10)

# selected the predict weight with fyear and cohort effect only.
ewaa_m3_sim48 <- ewaa_m3[, c("real_year", "age", "ewaa_fyear_fcoh48")]
head(ewaa_m3_sim48, n=800)

# rename the columns as SS3 format ewaa.ss
ewaa_m3_sim48 <- ewaa_m3_sim48 %>%
  rename(year = real_year, pred_weight = ewaa_fyear_fcoh48)
colnames(ewaa_m3_sim48) <- c("year", "age", "pred_weight")
head(ewaa_m3_sim48, n=800)

# Convert the year column to numeric if necessary
ewaa_m3_sim48$year <- as.numeric(ewaa_m3_sim48$year)

# Run the pad_weight_at_age function and inspect the output
padded_data48 <- pad_weight_at_age(ewaa_m3_sim48)
str(padded_data48)

hakedataUSA:::write_wtatage_file(file="ewaa_sim48.ss", data=padded_data48, maturity = maturity_at_age)

#########################################################################################
# Sim 49
set.seed(4925)
devs_fyear49 <- rnorm(nT, mean = 0, sd = sd_fyear)
devs_fyear_fut49 <- rep(devs_fyear49, each = 16)
devs_fcohort49 <- rnorm(nT_23, mean = 0, sd = sd_fcohort)
devs_fcohort_past_future49 <- c(fcohort_estimate_values, devs_fcohort49)
estimate_fcohort_all49 <- data.frame(year = years_cohort, estimate = devs_fcohort_past_future49)
result_matrix49 <- matrix(0, nrow = 17, ncol = num_cols)
result_matrix49[1, ] <- format(as.integer(estimate_fcohort_all49$year), nsmall = 0)
result_matrix49[2, ] <- estimate_fcohort_all49$estimate
print(result_matrix49)
num_rows49 <- nrow(result_matrix49)
num_cols49 <- ncol(result_matrix49)
value_to_repeat49 <- result_matrix49[2, ]
start_row49 <- 3
start_col49 <- 2
end_row49 <- num_rows49
end_col49 <- num_cols49
for (i in start_row49:end_row49) {
  for (j in start_col49:end_col49) {
    if (j >= (i - start_row49 + start_col49)) {
      pattern_index49 <- (j - (i - start_row49 + start_col49)) %% 32 + 1
      result_matrix49[i, j] <- value_to_repeat49[pattern_index49]
    }
  }
}
print(result_matrix49)
extracted_values49 <- as.numeric(result_matrix49[2:17, 16:33])
extracted_values49 <- extracted_values49[-1]
devs_fcohort_all49 <- c(devs_fcohort_past, extracted_values49)
devs_fyear_all49 <- c(devs_fyear_past, devs_fyear_fut49)
ewaa_m3$fyear_devs49 <- devs_fyear_all49
ewaa_m3$fcohort_devs49 <- devs_fcohort_all49
ewaa_m3$ewaa_fyear49 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fyear_devs49)
ewaa_m3$ewaa_fcohort49 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fcohort_devs49)
ewaa_m3$ewaa_fyear_fcoh49 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fyear_devs49 + ewaa_m3$fcohort_devs49)
head(ewaa_m3, n = 10)

# selected the predict weight with fyear and cohort effect only.
ewaa_m3_sim49 <- ewaa_m3[, c("real_year", "age", "ewaa_fyear_fcoh49")]
head(ewaa_m3_sim49, n=800)

# rename the columns as SS3 format ewaa.ss
ewaa_m3_sim49 <- ewaa_m3_sim49 %>%
  rename(year = real_year, pred_weight = ewaa_fyear_fcoh49)
colnames(ewaa_m3_sim49) <- c("year", "age", "pred_weight")
head(ewaa_m3_sim49, n=800)

# Convert the year column to numeric if necessary
ewaa_m3_sim49$year <- as.numeric(ewaa_m3_sim49$year)

# Run the pad_weight_at_age function and inspect the output
padded_data49 <- pad_weight_at_age(ewaa_m3_sim49)
str(padded_data49)

hakedataUSA:::write_wtatage_file(file="ewaa_sim49.ss", data=padded_data49, maturity = maturity_at_age)

#################################################################################
# Sim 50
set.seed(5049)
devs_fyear50 <- rnorm(nT, mean = 0, sd = sd_fyear)
devs_fyear_fut50 <- rep(devs_fyear50, each = 16)
devs_fcohort50 <- rnorm(nT_23, mean = 0, sd = sd_fcohort)
devs_fcohort_past_future50 <- c(fcohort_estimate_values, devs_fcohort50)
estimate_fcohort_all50 <- data.frame(year = years_cohort, estimate = devs_fcohort_past_future50)
result_matrix50 <- matrix(0, nrow = 17, ncol = num_cols)
result_matrix50[1, ] <- format(as.integer(estimate_fcohort_all50$year), nsmall = 0)
result_matrix50[2, ] <- estimate_fcohort_all50$estimate
print(result_matrix50)
num_rows50 <- nrow(result_matrix50)
num_cols50 <- ncol(result_matrix50)
value_to_repeat50 <- result_matrix50[2, ]
start_row50 <- 3
start_col50 <- 2
end_row50 <- num_rows50
end_col50 <- num_cols50
for (i in start_row50:end_row50) {
  for (j in start_col50:end_col50) {
    if (j >= (i - start_row50 + start_col50)) {
      pattern_index50 <- (j - (i - start_row50 + start_col50)) %% 32 + 1
      result_matrix50[i, j] <- value_to_repeat50[pattern_index50]
    }
  }
}
print(result_matrix50)
extracted_values50 <- as.numeric(result_matrix50[2:17, 16:33])
extracted_values50 <- extracted_values50[-1]
devs_fcohort_all50 <- c(devs_fcohort_past, extracted_values50)
devs_fyear_all50 <- c(devs_fyear_past, devs_fyear_fut50)
ewaa_m3$fyear_devs50 <- devs_fyear_all50
ewaa_m3$fcohort_devs50 <- devs_fcohort_all50
ewaa_m3$ewaa_fyear50 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fyear_devs50)
ewaa_m3$ewaa_fcohort50 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fcohort_devs50)
ewaa_m3$ewaa_fyear_fcoh50 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fyear_devs50 + ewaa_m3$fcohort_devs50)
head(ewaa_m3, n = 10)

# selected the predict weight with fyear and cohort effect only.
ewaa_m3_sim50 <- ewaa_m3[, c("real_year", "age", "ewaa_fyear_fcoh50")]
head(ewaa_m3_sim50, n=800)

# rename the columns as SS3 format ewaa.ss
ewaa_m3_sim50 <- ewaa_m3_sim50 %>%
  rename(year = real_year, pred_weight = ewaa_fyear_fcoh50)
colnames(ewaa_m3_sim50) <- c("year", "age", "pred_weight")
head(ewaa_m3_sim50, n=800)

# Convert the year column to numeric if necessary
ewaa_m3_sim50$year <- as.numeric(ewaa_m3_sim50$year)

# Run the pad_weight_at_age function and inspect the output
padded_data50 <- pad_weight_at_age(ewaa_m3_sim50)
str(padded_data50)

hakedataUSA:::write_wtatage_file(file="ewaa_sim50.ss", data=padded_data50, maturity = maturity_at_age)

########################################################################
# Sim 51
set.seed(5127)
devs_fyear51 <- rnorm(nT, mean = 0, sd = sd_fyear)
devs_fyear_fut51 <- rep(devs_fyear51, each = 16)
devs_fcohort51 <- rnorm(nT_23, mean = 0, sd = sd_fcohort)
devs_fcohort_past_future51 <- c(fcohort_estimate_values, devs_fcohort51)
estimate_fcohort_all51 <- data.frame(year = years_cohort, estimate = devs_fcohort_past_future51)
result_matrix51 <- matrix(0, nrow = 17, ncol = num_cols)
result_matrix51[1, ] <- format(as.integer(estimate_fcohort_all51$year), nsmall = 0)
result_matrix51[2, ] <- estimate_fcohort_all51$estimate
print(result_matrix51)
num_rows51 <- nrow(result_matrix51)
num_cols51 <- ncol(result_matrix51)
value_to_repeat51 <- result_matrix51[2, ]
start_row51 <- 3
start_col51 <- 2
end_row51 <- num_rows51
end_col51 <- num_cols51
for (i in start_row51:end_row51) {
  for (j in start_col51:end_col51) {
    if (j >= (i - start_row51 + start_col51)) {
      pattern_index51 <- (j - (i - start_row51 + start_col51)) %% 32 + 1
      result_matrix51[i, j] <- value_to_repeat51[pattern_index51]
    }
  }
}
print(result_matrix51)
extracted_values51 <- as.numeric(result_matrix51[2:17, 16:33])
extracted_values51 <- extracted_values51[-1]
devs_fcohort_all51 <- c(devs_fcohort_past, extracted_values51)
devs_fyear_all51 <- c(devs_fyear_past, devs_fyear_fut51)
ewaa_m3$fyear_devs51 <- devs_fyear_all51
ewaa_m3$fcohort_devs51 <- devs_fcohort_all51
ewaa_m3$ewaa_fyear51 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fyear_devs51)
ewaa_m3$ewaa_fcohort51 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fcohort_devs51)
ewaa_m3$ewaa_fyear_fcoh51 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fyear_devs51 + ewaa_m3$fcohort_devs51)
head(ewaa_m3, n = 10)

# selected the predict weight with fyear and cohort effect only.
ewaa_m3_sim51 <- ewaa_m3[, c("real_year", "age", "ewaa_fyear_fcoh51")]
head(ewaa_m3_sim51, n=800)

# rename the columns as SS3 format ewaa.ss
ewaa_m3_sim51 <- ewaa_m3_sim51 %>%
  rename(year = real_year, pred_weight = ewaa_fyear_fcoh51)
colnames(ewaa_m3_sim51) <- c("year", "age", "pred_weight")
head(ewaa_m3_sim51, n=800)

# Convert the year column to numeric if necessary
ewaa_m3_sim51$year <- as.numeric(ewaa_m3_sim51$year)

# Run the pad_weight_at_age function and inspect the output
padded_data51 <- pad_weight_at_age(ewaa_m3_sim51)
str(padded_data51)

hakedataUSA:::write_wtatage_file(file="ewaa_sim51.ss", data=padded_data51, maturity = maturity_at_age)


##########################################################################
# Sim 52
set.seed(528)
devs_fyear52 <- rnorm(nT, mean = 0, sd = sd_fyear)
devs_fyear_fut52 <- rep(devs_fyear52, each = 16)
devs_fcohort52 <- rnorm(nT_23, mean = 0, sd = sd_fcohort)
devs_fcohort_past_future52 <- c(fcohort_estimate_values, devs_fcohort52)
estimate_fcohort_all52 <- data.frame(year = years_cohort, estimate = devs_fcohort_past_future52)
result_matrix52 <- matrix(0, nrow = 17, ncol = num_cols)
result_matrix52[1, ] <- format(as.integer(estimate_fcohort_all52$year), nsmall = 0)
result_matrix52[2, ] <- estimate_fcohort_all52$estimate
print(result_matrix52)
num_rows52 <- nrow(result_matrix52)
num_cols52 <- ncol(result_matrix52)
value_to_repeat52 <- result_matrix52[2, ]
start_row52 <- 3
start_col52 <- 2
end_row52 <- num_rows52
end_col52 <- num_cols52
for (i in start_row52:end_row52) {
  for (j in start_col52:end_col52) {
    if (j >= (i - start_row52 + start_col52)) {
      pattern_index52 <- (j - (i - start_row52 + start_col52)) %% 32 + 1
      result_matrix52[i, j] <- value_to_repeat52[pattern_index52]
    }
  }
}
print(result_matrix52)
extracted_values52 <- as.numeric(result_matrix52[2:17, 16:33])
extracted_values52 <- extracted_values52[-1]
devs_fcohort_all52 <- c(devs_fcohort_past, extracted_values52)
devs_fyear_all52 <- c(devs_fyear_past, devs_fyear_fut52)
ewaa_m3$fyear_devs52 <- devs_fyear_all52
ewaa_m3$fcohort_devs52 <- devs_fcohort_all52
ewaa_m3$ewaa_fyear52 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fyear_devs52)
ewaa_m3$ewaa_fcohort52 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fcohort_devs52)
ewaa_m3$ewaa_fyear_fcoh52 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fyear_devs52 + ewaa_m3$fcohort_devs52)
head(ewaa_m3, n = 10)

# Selecciona el peso predicho con el efecto de año y cohorte solamente.
ewaa_m3_sim52 <- ewaa_m3[, c("real_year", "age", "ewaa_fyear_fcoh52")]
head(ewaa_m3_sim52, n = 800)

# Renombra las columnas al formato SS3 ewaa.ss
ewaa_m3_sim52 <- ewaa_m3_sim52 %>%
  rename(year = real_year, pred_weight = ewaa_fyear_fcoh52)
colnames(ewaa_m3_sim52) <- c("year", "age", "pred_weight")
head(ewaa_m3_sim52, n = 800)

# Convierte la columna de año a numérica si es necesario
ewaa_m3_sim52$year <- as.numeric(ewaa_m3_sim52$year)

# Ejecuta la función pad_weight_at_age e inspecciona la salida
padded_data52 <- pad_weight_at_age(ewaa_m3_sim52)
str(padded_data52)

# Escribe el archivo con el formato SS3
hakedataUSA:::write_wtatage_file(file = "ewaa_sim52.ss", data = padded_data52, maturity = maturity_at_age)

#############################################################################
# Sim 53
set.seed(531)
devs_fyear53 <- rnorm(nT, mean = 0, sd = sd_fyear)
devs_fyear_fut53 <- rep(devs_fyear53, each = 16)
devs_fcohort53 <- rnorm(nT_23, mean = 0, sd = sd_fcohort)
devs_fcohort_past_future53 <- c(fcohort_estimate_values, devs_fcohort53)
estimate_fcohort_all53 <- data.frame(year = years_cohort, estimate = devs_fcohort_past_future53)
result_matrix53 <- matrix(0, nrow = 17, ncol = num_cols)
result_matrix53[1, ] <- format(as.integer(estimate_fcohort_all53$year), nsmall = 0)
result_matrix53[2, ] <- estimate_fcohort_all53$estimate
print(result_matrix53)
num_rows53 <- nrow(result_matrix53)
num_cols53 <- ncol(result_matrix53)
value_to_repeat53 <- result_matrix53[2, ]
start_row53 <- 3
start_col53 <- 2
end_row53 <- num_rows53
end_col53 <- num_cols53
for (i in start_row53:end_row53) {
  for (j in start_col53:end_col53) {
    if (j >= (i - start_row53 + start_col53)) {
      pattern_index53 <- (j - (i - start_row53 + start_col53)) %% 32 + 1
      result_matrix53[i, j] <- value_to_repeat53[pattern_index53]
    }
  }
}
print(result_matrix53)
extracted_values53 <- as.numeric(result_matrix53[2:17, 16:33])
extracted_values53 <- extracted_values53[-1]
devs_fcohort_all53 <- c(devs_fcohort_past, extracted_values53)
devs_fyear_all53 <- c(devs_fyear_past, devs_fyear_fut53)
ewaa_m3$fyear_devs53 <- devs_fyear_all53
ewaa_m3$fcohort_devs53 <- devs_fcohort_all53
ewaa_m3$ewaa_fyear53 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fyear_devs53)
ewaa_m3$ewaa_fcohort53 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fcohort_devs53)
ewaa_m3$ewaa_fyear_fcoh53 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fyear_devs53 + ewaa_m3$fcohort_devs53)
head(ewaa_m3, n = 10)

# Selecciona el peso predicho con el efecto de año y cohorte solamente.
ewaa_m3_sim53 <- ewaa_m3[, c("real_year", "age", "ewaa_fyear_fcoh53")]
head(ewaa_m3_sim53, n = 800)

# Renombra las columnas al formato SS3 ewaa.ss
ewaa_m3_sim53 <- ewaa_m3_sim53 %>%
  rename(year = real_year, pred_weight = ewaa_fyear_fcoh53)
colnames(ewaa_m3_sim53) <- c("year", "age", "pred_weight")
head(ewaa_m3_sim53, n = 800)

# Convierte la columna de año a numérica si es necesario
ewaa_m3_sim53$year <- as.numeric(ewaa_m3_sim53$year)

# Ejecuta la función pad_weight_at_age e inspecciona la salida
padded_data53 <- pad_weight_at_age(ewaa_m3_sim53)
str(padded_data53)

# Escribe el archivo con el formato SS3
hakedataUSA:::write_wtatage_file(file = "ewaa_sim53.ss", data = padded_data53, maturity = maturity_at_age)

################################################################################
# Sim 54
set.seed(549)
devs_fyear54 <- rnorm(nT, mean = 0, sd = sd_fyear)
devs_fyear_fut54 <- rep(devs_fyear54, each = 16)
devs_fcohort54 <- rnorm(nT_23, mean = 0, sd = sd_fcohort)
devs_fcohort_past_future54 <- c(fcohort_estimate_values, devs_fcohort54)
estimate_fcohort_all54 <- data.frame(year = years_cohort, estimate = devs_fcohort_past_future54)
result_matrix54 <- matrix(0, nrow = 17, ncol = num_cols)
result_matrix54[1, ] <- format(as.integer(estimate_fcohort_all54$year), nsmall = 0)
result_matrix54[2, ] <- estimate_fcohort_all54$estimate
print(result_matrix54)
num_rows54 <- nrow(result_matrix54)
num_cols54 <- ncol(result_matrix54)
value_to_repeat54 <- result_matrix54[2, ]
start_row54 <- 3
start_col54 <- 2
end_row54 <- num_rows54
end_col54 <- num_cols54
for (i in start_row54:end_row54) {
  for (j in start_col54:end_col54) {
    if (j >= (i - start_row54 + start_col54)) {
      pattern_index54 <- (j - (i - start_row54 + start_col54)) %% 32 + 1
      result_matrix54[i, j] <- value_to_repeat54[pattern_index54]
    }
  }
}
print(result_matrix54)
extracted_values54 <- as.numeric(result_matrix54[2:17, 16:33])
extracted_values54 <- extracted_values54[-1]
devs_fcohort_all54 <- c(devs_fcohort_past, extracted_values54)
devs_fyear_all54 <- c(devs_fyear_past, devs_fyear_fut54)
ewaa_m3$fyear_devs54 <- devs_fyear_all54
ewaa_m3$fcohort_devs54 <- devs_fcohort_all54
ewaa_m3$ewaa_fyear54 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fyear_devs54)
ewaa_m3$ewaa_fcohort54 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fcohort_devs54)
ewaa_m3$ewaa_fyear_fcoh54 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fyear_devs54 + ewaa_m3$fcohort_devs54)
head(ewaa_m3, n = 10)

# Selecciona el peso predicho con el efecto de año y cohorte solamente.
ewaa_m3_sim54 <- ewaa_m3[, c("real_year", "age", "ewaa_fyear_fcoh54")]
head(ewaa_m3_sim54, n = 800)

# Renombra las columnas al formato SS3 ewaa.ss
ewaa_m3_sim54 <- ewaa_m3_sim54 %>%
  rename(year = real_year, pred_weight = ewaa_fyear_fcoh54)
colnames(ewaa_m3_sim54) <- c("year", "age", "pred_weight")
head(ewaa_m3_sim54, n = 800)

# Convierte la columna de año a numérica si es necesario
ewaa_m3_sim54$year <- as.numeric(ewaa_m3_sim54$year)

# Ejecuta la función pad_weight_at_age e inspecciona la salida
padded_data54 <- pad_weight_at_age(ewaa_m3_sim54)
str(padded_data54)

# Escribe el archivo con el formato SS3
hakedataUSA:::write_wtatage_file(file = "ewaa_sim54.ss", data = padded_data54, maturity = maturity_at_age)

############################################################################
# Sim 55
set.seed(556)
devs_fyear55 <- rnorm(nT, mean = 0, sd = sd_fyear)
devs_fyear_fut55 <- rep(devs_fyear55, each = 16)
devs_fcohort55 <- rnorm(nT_23, mean = 0, sd = sd_fcohort)
devs_fcohort_past_future55 <- c(fcohort_estimate_values, devs_fcohort55)
estimate_fcohort_all55 <- data.frame(year = years_cohort, estimate = devs_fcohort_past_future55)
result_matrix55 <- matrix(0, nrow = 17, ncol = num_cols)
result_matrix55[1, ] <- format(as.integer(estimate_fcohort_all55$year), nsmall = 0)
result_matrix55[2, ] <- estimate_fcohort_all55$estimate
print(result_matrix55)
num_rows55 <- nrow(result_matrix55)
num_cols55 <- ncol(result_matrix55)
value_to_repeat55 <- result_matrix55[2, ]
start_row55 <- 3
start_col55 <- 2
end_row55 <- num_rows55
end_col55 <- num_cols55
for (i in start_row55:end_row55) {
  for (j in start_col55:end_col55) {
    if (j >= (i - start_row55 + start_col55)) {
      pattern_index55 <- (j - (i - start_row55 + start_col55)) %% 32 + 1
      result_matrix55[i, j] <- value_to_repeat55[pattern_index55]
    }
  }
}
print(result_matrix55)
extracted_values55 <- as.numeric(result_matrix55[2:17, 16:33])
extracted_values55 <- extracted_values55[-1]
devs_fcohort_all55 <- c(devs_fcohort_past, extracted_values55)
devs_fyear_all55 <- c(devs_fyear_past, devs_fyear_fut55)
ewaa_m3$fyear_devs55 <- devs_fyear_all55
ewaa_m3$fcohort_devs55 <- devs_fcohort_all55
ewaa_m3$ewaa_fyear55 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fyear_devs55)
ewaa_m3$ewaa_fcohort55 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fcohort_devs55)
ewaa_m3$ewaa_fyear_fcoh55 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fyear_devs55 + ewaa_m3$fcohort_devs55)
head(ewaa_m3, n = 10)

# Selecciona el peso predicho con el efecto de año y cohorte solamente.
ewaa_m3_sim55 <- ewaa_m3[, c("real_year", "age", "ewaa_fyear_fcoh55")]
head(ewaa_m3_sim55, n = 800)

# Renombra las columnas al formato SS3 ewaa.ss
ewaa_m3_sim55 <- ewaa_m3_sim55 %>%
  rename(year = real_year, pred_weight = ewaa_fyear_fcoh55)
colnames(ewaa_m3_sim55) <- c("year", "age", "pred_weight")
head(ewaa_m3_sim55, n = 800)

# Convierte la columna de año a numérica si es necesario
ewaa_m3_sim55$year <- as.numeric(ewaa_m3_sim55$year)

# Ejecuta la función pad_weight_at_age e inspecciona la salida
padded_data55 <- pad_weight_at_age(ewaa_m3_sim55)
str(padded_data55)

# Escribe el archivo con el formato SS3
hakedataUSA:::write_wtatage_file(file = "ewaa_sim55.ss", data = padded_data55, maturity = maturity_at_age)

# Sim 56
set.seed(563)
devs_fyear56 <- rnorm(nT, mean = 0, sd = sd_fyear)
devs_fyear_fut56 <- rep(devs_fyear56, each = 16)
devs_fcohort56 <- rnorm(nT_23, mean = 0, sd = sd_fcohort)
devs_fcohort_past_future56 <- c(fcohort_estimate_values, devs_fcohort56)
estimate_fcohort_all56 <- data.frame(year = years_cohort, estimate = devs_fcohort_past_future56)
result_matrix56 <- matrix(0, nrow = 17, ncol = num_cols)
result_matrix56[1, ] <- format(as.integer(estimate_fcohort_all56$year), nsmall = 0)
result_matrix56[2, ] <- estimate_fcohort_all56$estimate
print(result_matrix56)
num_rows56 <- nrow(result_matrix56)
num_cols56 <- ncol(result_matrix56)
value_to_repeat56 <- result_matrix56[2, ]
start_row56 <- 3
start_col56 <- 2
end_row56 <- num_rows56
end_col56 <- num_cols56
for (i in start_row56:end_row56) {
  for (j in start_col56:end_col56) {
    if (j >= (i - start_row56 + start_col56)) {
      pattern_index56 <- (j - (i - start_row56 + start_col56)) %% 32 + 1
      result_matrix56[i, j] <- value_to_repeat56[pattern_index56]
    }
  }
}
print(result_matrix56)
extracted_values56 <- as.numeric(result_matrix56[2:17, 16:33])
extracted_values56 <- extracted_values56[-1]
devs_fcohort_all56 <- c(devs_fcohort_past, extracted_values56)
devs_fyear_all56 <- c(devs_fyear_past, devs_fyear_fut56)
ewaa_m3$fyear_devs56 <- devs_fyear_all56
ewaa_m3$fcohort_devs56 <- devs_fcohort_all56
ewaa_m3$ewaa_fyear56 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fyear_devs56)
ewaa_m3$ewaa_fcohort56 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fcohort_devs56)
ewaa_m3$ewaa_fyear_fcoh56 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fyear_devs56 + ewaa_m3$fcohort_devs56)
head(ewaa_m3, n = 10)

# Selecciona el peso predicho con el efecto de año y cohorte solamente.
ewaa_m3_sim56 <- ewaa_m3[, c("real_year", "age", "ewaa_fyear_fcoh56")]
head(ewaa_m3_sim56, n = 800)

# Renombra las columnas al formato SS3 ewaa.ss
ewaa_m3_sim56 <- ewaa_m3_sim56 %>%
  rename(year = real_year, pred_weight = ewaa_fyear_fcoh56)
colnames(ewaa_m3_sim56) <- c("year", "age", "pred_weight")
head(ewaa_m3_sim56, n = 800)

# Convierte la columna de año a numérica si es necesario
ewaa_m3_sim56$year <- as.numeric(ewaa_m3_sim56$year)

# Ejecuta la función pad_weight_at_age e inspecciona la salida
padded_data56 <- pad_weight_at_age(ewaa_m3_sim56)
str(padded_data56)

# Escribe el archivo con el formato SS3
hakedataUSA:::write_wtatage_file(file = "ewaa_sim56.ss", data = padded_data56, maturity = maturity_at_age)

####################################################################
# Sim 57
set.seed(577)
devs_fyear57 <- rnorm(nT, mean = 0, sd = sd_fyear)
devs_fyear_fut57 <- rep(devs_fyear57, each = 16)
devs_fcohort57 <- rnorm(nT_23, mean = 0, sd = sd_fcohort)
devs_fcohort_past_future57 <- c(fcohort_estimate_values, devs_fcohort57)
estimate_fcohort_all57 <- data.frame(year = years_cohort, estimate = devs_fcohort_past_future57)
result_matrix57 <- matrix(0, nrow = 17, ncol = num_cols)
result_matrix57[1, ] <- format(as.integer(estimate_fcohort_all57$year), nsmall = 0)
result_matrix57[2, ] <- estimate_fcohort_all57$estimate
print(result_matrix57)
num_rows57 <- nrow(result_matrix57)
num_cols57 <- ncol(result_matrix57)
value_to_repeat57 <- result_matrix57[2, ]
start_row57 <- 3
start_col57 <- 2
end_row57 <- num_rows57
end_col57 <- num_cols57
for (i in start_row57:end_row57) {
  for (j in start_col57:end_col57) {
    if (j >= (i - start_row57 + start_col57)) {
      pattern_index57 <- (j - (i - start_row57 + start_col57)) %% 32 + 1
      result_matrix57[i, j] <- value_to_repeat57[pattern_index57]
    }
  }
}
print(result_matrix57)
extracted_values57 <- as.numeric(result_matrix57[2:17, 16:33])
extracted_values57 <- extracted_values57[-1]
devs_fcohort_all57 <- c(devs_fcohort_past, extracted_values57)
devs_fyear_all57 <- c(devs_fyear_past, devs_fyear_fut57)
ewaa_m3$fyear_devs57 <- devs_fyear_all57
ewaa_m3$fcohort_devs57 <- devs_fcohort_all57
ewaa_m3$ewaa_fyear57 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fyear_devs57)
ewaa_m3$ewaa_fcohort57 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fcohort_devs57)
ewaa_m3$ewaa_fyear_fcoh57 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fyear_devs57 + ewaa_m3$fcohort_devs57)
head(ewaa_m3, n = 10)

# Selecciona el peso predicho con el efecto de año y cohorte solamente.
ewaa_m3_sim57 <- ewaa_m3[, c("real_year", "age", "ewaa_fyear_fcoh57")]
head(ewaa_m3_sim57, n = 800)

# Renombra las columnas al formato SS3 ewaa.ss
ewaa_m3_sim57 <- ewaa_m3_sim57 %>%
  rename(year = real_year, pred_weight = ewaa_fyear_fcoh57)
colnames(ewaa_m3_sim57) <- c("year", "age", "pred_weight")
head(ewaa_m3_sim57, n = 800)

# Convierte la columna de año a numérica si es necesario
ewaa_m3_sim57$year <- as.numeric(ewaa_m3_sim57$year)

# Ejecuta la función pad_weight_at_age e inspecciona la salida
padded_data57 <- pad_weight_at_age(ewaa_m3_sim57)
str(padded_data57)

# Escribe el archivo con el formato SS3
hakedataUSA:::write_wtatage_file(file = "ewaa_sim57.ss", data = padded_data57, maturity = maturity_at_age)

###############################################################################
# Sim 58
set.seed(583)
devs_fyear58 <- rnorm(nT, mean = 0, sd = sd_fyear)
devs_fyear_fut58 <- rep(devs_fyear58, each = 16)
devs_fcohort58 <- rnorm(nT_23, mean = 0, sd = sd_fcohort)
devs_fcohort_past_future58 <- c(fcohort_estimate_values, devs_fcohort58)
estimate_fcohort_all58 <- data.frame(year = years_cohort, estimate = devs_fcohort_past_future58)
result_matrix58 <- matrix(0, nrow = 17, ncol = num_cols)
result_matrix58[1, ] <- format(as.integer(estimate_fcohort_all58$year), nsmall = 0)
result_matrix58[2, ] <- estimate_fcohort_all58$estimate
print(result_matrix58)
num_rows58 <- nrow(result_matrix58)
num_cols58 <- ncol(result_matrix58)
value_to_repeat58 <- result_matrix58[2, ]
start_row58 <- 3
start_col58 <- 2
end_row58 <- num_rows58
end_col58 <- num_cols58
for (i in start_row58:end_row58) {
  for (j in start_col58:end_col58) {
    if (j >= (i - start_row58 + start_col58)) {
      pattern_index58 <- (j - (i - start_row58 + start_col58)) %% 32 + 1
      result_matrix58[i, j] <- value_to_repeat58[pattern_index58]
    }
  }
}
print(result_matrix58)
extracted_values58 <- as.numeric(result_matrix58[2:17, 16:33])
extracted_values58 <- extracted_values58[-1]
devs_fcohort_all58 <- c(devs_fcohort_past, extracted_values58)
devs_fyear_all58 <- c(devs_fyear_past, devs_fyear_fut58)
ewaa_m3$fyear_devs58 <- devs_fyear_all58
ewaa_m3$fcohort_devs58 <- devs_fcohort_all58
ewaa_m3$ewaa_fyear58 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fyear_devs58)
ewaa_m3$ewaa_fcohort58 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fcohort_devs58)
ewaa_m3$ewaa_fyear_fcoh58 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fyear_devs58 + ewaa_m3$fcohort_devs58)
head(ewaa_m3, n = 10)

# Selecciona el peso predicho con el efecto de año y cohorte solamente.
ewaa_m3_sim58 <- ewaa_m3[, c("real_year", "age", "ewaa_fyear_fcoh58")]
head(ewaa_m3_sim58, n = 800)

# Renombra las columnas al formato SS3 ewaa.ss
ewaa_m3_sim58 <- ewaa_m3_sim58 %>%
  rename(year = real_year, pred_weight = ewaa_fyear_fcoh58)
colnames(ewaa_m3_sim58) <- c("year", "age", "pred_weight")
head(ewaa_m3_sim58, n = 800)

# Convierte la columna de año a numérica si es necesario
ewaa_m3_sim58$year <- as.numeric(ewaa_m3_sim58$year)

# Ejecuta la función pad_weight_at_age e inspecciona la salida
padded_data58 <- pad_weight_at_age(ewaa_m3_sim58)
str(padded_data58)

# Escribe el archivo con el formato SS3
hakedataUSA:::write_wtatage_file(file = "ewaa_sim58.ss", data = padded_data58, maturity = maturity_at_age)

##########################################################################
# Sim 59
set.seed(599)
devs_fyear59 <- rnorm(nT, mean = 0, sd = sd_fyear)
devs_fyear_fut59 <- rep(devs_fyear59, each = 16)
devs_fcohort59 <- rnorm(nT_23, mean = 0, sd = sd_fcohort)
devs_fcohort_past_future59 <- c(fcohort_estimate_values, devs_fcohort59)
estimate_fcohort_all59 <- data.frame(year = years_cohort, estimate = devs_fcohort_past_future59)
result_matrix59 <- matrix(0, nrow = 17, ncol = num_cols)
result_matrix59[1, ] <- format(as.integer(estimate_fcohort_all59$year), nsmall = 0)
result_matrix59[2, ] <- estimate_fcohort_all59$estimate
print(result_matrix59)
num_rows59 <- nrow(result_matrix59)
num_cols59 <- ncol(result_matrix59)
value_to_repeat59 <- result_matrix59[2, ]
start_row59 <- 3
start_col59 <- 2
end_row59 <- num_rows59
end_col59 <- num_cols59
for (i in start_row59:end_row59) {
  for (j in start_col59:end_col59) {
    if (j >= (i - start_row59 + start_col59)) {
      pattern_index59 <- (j - (i - start_row59 + start_col59)) %% 32 + 1
      result_matrix59[i, j] <- value_to_repeat59[pattern_index59]
    }
  }
}
print(result_matrix59)
extracted_values59 <- as.numeric(result_matrix59[2:17, 16:33])
extracted_values59 <- extracted_values59[-1]
devs_fcohort_all59 <- c(devs_fcohort_past, extracted_values59)
devs_fyear_all59 <- c(devs_fyear_past, devs_fyear_fut59)
ewaa_m3$fyear_devs59 <- devs_fyear_all59
ewaa_m3$fcohort_devs59 <- devs_fcohort_all59
ewaa_m3$ewaa_fyear59 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fyear_devs59)
ewaa_m3$ewaa_fcohort59 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fcohort_devs59)
ewaa_m3$ewaa_fyear_fcoh59 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fyear_devs59 + ewaa_m3$fcohort_devs59)
head(ewaa_m3, n = 10)

# Selecciona el peso predicho con el efecto de año y cohorte solamente.
ewaa_m3_sim59 <- ewaa_m3[, c("real_year", "age", "ewaa_fyear_fcoh59")]
head(ewaa_m3_sim59, n = 800)

# Renombra las columnas al formato SS3 ewaa.ss
ewaa_m3_sim59 <- ewaa_m3_sim59 %>%
  rename(year = real_year, pred_weight = ewaa_fyear_fcoh59)
colnames(ewaa_m3_sim59) <- c("year", "age", "pred_weight")
head(ewaa_m3_sim59, n = 800)

# Convierte la columna de año a numérica si es necesario
ewaa_m3_sim59$year <- as.numeric(ewaa_m3_sim59$year)

# Ejecuta la función pad_weight_at_age e inspecciona la salida
padded_data59 <- pad_weight_at_age(ewaa_m3_sim59)
str(padded_data59)

# Escribe el archivo con el formato SS3
hakedataUSA:::write_wtatage_file(file = "ewaa_sim59.ss", data = padded_data59, maturity = maturity_at_age)

#############################################################################
# Sim 60
set.seed(601)
devs_fyear60 <- rnorm(nT, mean = 0, sd = sd_fyear)
devs_fyear_fut60 <- rep(devs_fyear60, each = 16)
devs_fcohort60 <- rnorm(nT_23, mean = 0, sd = sd_fcohort)
devs_fcohort_past_future60 <- c(fcohort_estimate_values, devs_fcohort60)
estimate_fcohort_all60 <- data.frame(year = years_cohort, estimate = devs_fcohort_past_future60)
result_matrix60 <- matrix(0, nrow = 17, ncol = num_cols)
result_matrix60[1, ] <- format(as.integer(estimate_fcohort_all60$year), nsmall = 0)
result_matrix60[2, ] <- estimate_fcohort_all60$estimate
print(result_matrix60)
num_rows60 <- nrow(result_matrix60)
num_cols60 <- ncol(result_matrix60)
value_to_repeat60 <- result_matrix60[2, ]
start_row60 <- 3
start_col60 <- 2
end_row60 <- num_rows60
end_col60 <- num_cols60
for (i in start_row60:end_row60) {
  for (j in start_col60:end_col60) {
    if (j >= (i - start_row60 + start_col60)) {
      pattern_index60 <- (j - (i - start_row60 + start_col60)) %% 32 + 1
      result_matrix60[i, j] <- value_to_repeat60[pattern_index60]
    }
  }
}
print(result_matrix60)
extracted_values60 <- as.numeric(result_matrix60[2:17, 16:33])
extracted_values60 <- extracted_values60[-1]
devs_fcohort_all60 <- c(devs_fcohort_past, extracted_values60)
devs_fyear_all60 <- c(devs_fyear_past, devs_fyear_fut60)
ewaa_m3$fyear_devs60 <- devs_fyear_all60
ewaa_m3$fcohort_devs60 <- devs_fcohort_all60
ewaa_m3$ewaa_fyear60 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fyear_devs60)
ewaa_m3$ewaa_fcohort60 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fcohort_devs60)
ewaa_m3$ewaa_fyear_fcoh60 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fyear_devs60 + ewaa_m3$fcohort_devs60)
head(ewaa_m3, n = 10)

# Selecciona el peso predicho con el efecto de año y cohorte solamente.
ewaa_m3_sim60 <- ewaa_m3[, c("real_year", "age", "ewaa_fyear_fcoh60")]
head(ewaa_m3_sim60, n = 800)

# Renombra las columnas al formato SS3 ewaa.ss
ewaa_m3_sim60 <- ewaa_m3_sim60 %>%
  rename(year = real_year, pred_weight = ewaa_fyear_fcoh60)
colnames(ewaa_m3_sim60) <- c("year", "age", "pred_weight")
head(ewaa_m3_sim60, n = 800)

# Convierte la columna de año a numérica si es necesario
ewaa_m3_sim60$year <- as.numeric(ewaa_m3_sim60$year)

# Ejecuta la función pad_weight_at_age e inspecciona la salida
padded_data60 <- pad_weight_at_age(ewaa_m3_sim60)
str(padded_data60)

# Escribe el archivo con el formato SS3
hakedataUSA:::write_wtatage_file(file = "ewaa_sim60.ss", data = padded_data60, maturity = maturity_at_age)

#############################################################################
# Sim 61
set.seed(617)
devs_fyear61 <- rnorm(nT, mean = 0, sd = sd_fyear)
devs_fyear_fut61 <- rep(devs_fyear61, each = 16)
devs_fcohort61 <- rnorm(nT_23, mean = 0, sd = sd_fcohort)
devs_fcohort_past_future61 <- c(fcohort_estimate_values, devs_fcohort61)
estimate_fcohort_all61 <- data.frame(year = years_cohort, estimate = devs_fcohort_past_future61)
result_matrix61 <- matrix(0, nrow = 17, ncol = num_cols)
result_matrix61[1, ] <- format(as.integer(estimate_fcohort_all61$year), nsmall = 0)
result_matrix61[2, ] <- estimate_fcohort_all61$estimate
print(result_matrix61)
num_rows61 <- nrow(result_matrix61)
num_cols61 <- ncol(result_matrix61)
value_to_repeat61 <- result_matrix61[2, ]
start_row61 <- 3
start_col61 <- 2
end_row61 <- num_rows61
end_col61 <- num_cols61
for (i in start_row61:end_row61) {
  for (j in start_col61:end_col61) {
    if (j >= (i - start_row61 + start_col61)) {
      pattern_index61 <- (j - (i - start_row61 + start_col61)) %% 32 + 1
      result_matrix61[i, j] <- value_to_repeat61[pattern_index61]
    }
  }
}
print(result_matrix61)
extracted_values61 <- as.numeric(result_matrix61[2:17, 16:33])
extracted_values61 <- extracted_values61[-1]
devs_fcohort_all61 <- c(devs_fcohort_past, extracted_values61)
devs_fyear_all61 <- c(devs_fyear_past, devs_fyear_fut61)
ewaa_m3$fyear_devs61 <- devs_fyear_all61
ewaa_m3$fcohort_devs61 <- devs_fcohort_all61
ewaa_m3$ewaa_fyear61 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fyear_devs61)
ewaa_m3$ewaa_fcohort61 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fcohort_devs61)
ewaa_m3$ewaa_fyear_fcoh61 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fyear_devs61 + ewaa_m3$fcohort_devs61)
head(ewaa_m3, n = 10)

# Selecciona el peso predicho con el efecto de año y cohorte solamente.
ewaa_m3_sim61 <- ewaa_m3[, c("real_year", "age", "ewaa_fyear_fcoh61")]
head(ewaa_m3_sim61, n = 800)

# Renombra las columnas al formato SS3 ewaa.ss
ewaa_m3_sim61 <- ewaa_m3_sim61 %>%
  rename(year = real_year, pred_weight = ewaa_fyear_fcoh61)
colnames(ewaa_m3_sim61) <- c("year", "age", "pred_weight")
head(ewaa_m3_sim61, n = 800)

# Convierte la columna de año a numérica si es necesario
ewaa_m3_sim61$year <- as.numeric(ewaa_m3_sim61$year)

# Ejecuta la función pad_weight_at_age e inspecciona la salida
padded_data61 <- pad_weight_at_age(ewaa_m3_sim61)
str(padded_data61)

# Escribe el archivo con el formato SS3
hakedataUSA:::write_wtatage_file(file = "ewaa_sim61.ss", data = padded_data61, maturity = maturity_at_age)

###############################################
# Sim 62
set.seed(627)
devs_fyear62 <- rnorm(nT, mean = 0, sd = sd_fyear)
devs_fyear_fut62 <- rep(devs_fyear62, each = 16)
devs_fcohort62 <- rnorm(nT_23, mean = 0, sd = sd_fcohort)
devs_fcohort_past_future62 <- c(fcohort_estimate_values, devs_fcohort62)
estimate_fcohort_all62 <- data.frame(year = years_cohort, estimate = devs_fcohort_past_future62)
result_matrix62 <- matrix(0, nrow = 17, ncol = num_cols)
result_matrix62[1, ] <- format(as.integer(estimate_fcohort_all62$year), nsmall = 0)
result_matrix62[2, ] <- estimate_fcohort_all62$estimate
print(result_matrix62)
num_rows62 <- nrow(result_matrix62)
num_cols62 <- ncol(result_matrix62)
value_to_repeat62 <- result_matrix62[2, ]
start_row62 <- 3
start_col62 <- 2
end_row62 <- num_rows62
end_col62 <- num_cols62
for (i in start_row62:end_row62) {
  for (j in start_col62:end_col62) {
    if (j >= (i - start_row62 + start_col62)) {
      pattern_index62 <- (j - (i - start_row62 + start_col62)) %% 32 + 1
      result_matrix62[i, j] <- value_to_repeat62[pattern_index62]
    }
  }
}
print(result_matrix62)
extracted_values62 <- as.numeric(result_matrix62[2:17, 16:33])
extracted_values62 <- extracted_values62[-1]
devs_fcohort_all62 <- c(devs_fcohort_past, extracted_values62)
devs_fyear_all62 <- c(devs_fyear_past, devs_fyear_fut62)
ewaa_m3$fyear_devs62 <- devs_fyear_all62
ewaa_m3$fcohort_devs62 <- devs_fcohort_all62
ewaa_m3$ewaa_fyear62 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fyear_devs62)
ewaa_m3$ewaa_fcohort62 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fcohort_devs62)
ewaa_m3$ewaa_fyear_fcoh62 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fyear_devs62 + ewaa_m3$fcohort_devs62)
head(ewaa_m3, n = 10)

# Selecciona el peso predicho con el efecto de año y cohorte solamente.
ewaa_m3_sim62 <- ewaa_m3[, c("real_year", "age", "ewaa_fyear_fcoh62")]
head(ewaa_m3_sim62, n = 800)

# Renombra las columnas al formato SS3 ewaa.ss
ewaa_m3_sim62 <- ewaa_m3_sim62 %>%
  rename(year = real_year, pred_weight = ewaa_fyear_fcoh62)
colnames(ewaa_m3_sim62) <- c("year", "age", "pred_weight")
head(ewaa_m3_sim62, n = 800)

# Convierte la columna de año a numérica si es necesario
ewaa_m3_sim62$year <- as.numeric(ewaa_m3_sim62$year)

# Ejecuta la función pad_weight_at_age e inspecciona la salida
padded_data62 <- pad_weight_at_age(ewaa_m3_sim62)
str(padded_data62)

# Escribe el archivo con el formato SS3
hakedataUSA:::write_wtatage_file(file = "ewaa_sim62.ss", data = padded_data62, maturity = maturity_at_age)

########################################################################
# Sim 63
set.seed(631)
devs_fyear63 <- rnorm(nT, mean = 0, sd = sd_fyear)
devs_fyear_fut63 <- rep(devs_fyear63, each = 16)
devs_fcohort63 <- rnorm(nT_23, mean = 0, sd = sd_fcohort)
devs_fcohort_past_future63 <- c(fcohort_estimate_values, devs_fcohort63)
estimate_fcohort_all63 <- data.frame(year = years_cohort, estimate = devs_fcohort_past_future63)
result_matrix63 <- matrix(0, nrow = 17, ncol = num_cols)
result_matrix63[1, ] <- format(as.integer(estimate_fcohort_all63$year), nsmall = 0)
result_matrix63[2, ] <- estimate_fcohort_all63$estimate
print(result_matrix63)
num_rows63 <- nrow(result_matrix63)
num_cols63 <- ncol(result_matrix63)
value_to_repeat63 <- result_matrix63[2, ]
start_row63 <- 3
start_col63 <- 2
end_row63 <- num_rows63
end_col63 <- num_cols63
for (i in start_row63:end_row63) {
  for (j in start_col63:end_col63) {
    if (j >= (i - start_row63 + start_col63)) {
      pattern_index63 <- (j - (i - start_row63 + start_col63)) %% 32 + 1
      result_matrix63[i, j] <- value_to_repeat63[pattern_index63]
    }
  }
}
print(result_matrix63)
extracted_values63 <- as.numeric(result_matrix63[2:17, 16:33])
extracted_values63 <- extracted_values63[-1]
devs_fcohort_all63 <- c(devs_fcohort_past, extracted_values63)
devs_fyear_all63 <- c(devs_fyear_past, devs_fyear_fut63)
ewaa_m3$fyear_devs63 <- devs_fyear_all63
ewaa_m3$fcohort_devs63 <- devs_fcohort_all63
ewaa_m3$ewaa_fyear63 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fyear_devs63)
ewaa_m3$ewaa_fcohort63 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fcohort_devs63)
ewaa_m3$ewaa_fyear_fcoh63 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fyear_devs63 + ewaa_m3$fcohort_devs63)
head(ewaa_m3, n = 10)

# Selecciona el peso predicho con el efecto de año y cohorte solamente.
ewaa_m3_sim63 <- ewaa_m3[, c("real_year", "age", "ewaa_fyear_fcoh63")]
head(ewaa_m3_sim63, n = 800)

# Renombra las columnas al formato SS3 ewaa.ss
ewaa_m3_sim63 <- ewaa_m3_sim63 %>%
  rename(year = real_year, pred_weight = ewaa_fyear_fcoh63)
colnames(ewaa_m3_sim63) <- c("year", "age", "pred_weight")
head(ewaa_m3_sim63, n = 800)

# Convierte la columna de año a numérica si es necesario
ewaa_m3_sim63$year <- as.numeric(ewaa_m3_sim63$year)

# Ejecuta la función pad_weight_at_age e inspecciona la salida
padded_data63 <- pad_weight_at_age(ewaa_m3_sim63)
str(padded_data63)

# Escribe el archivo con el formato SS3
hakedataUSA:::write_wtatage_file(file = "ewaa_sim63.ss", data = padded_data63, maturity = maturity_at_age)

#########################################################################
# Sim 64
set.seed(643)
devs_fyear64 <- rnorm(nT, mean = 0, sd = sd_fyear)
devs_fyear_fut64 <- rep(devs_fyear64, each = 16)
devs_fcohort64 <- rnorm(nT_23, mean = 0, sd = sd_fcohort)
devs_fcohort_past_future64 <- c(fcohort_estimate_values, devs_fcohort64)
estimate_fcohort_all64 <- data.frame(year = years_cohort, estimate = devs_fcohort_past_future64)
result_matrix64 <- matrix(0, nrow = 17, ncol = num_cols)
result_matrix64[1, ] <- format(as.integer(estimate_fcohort_all64$year), nsmall = 0)
result_matrix64[2, ] <- estimate_fcohort_all64$estimate
print(result_matrix64)
num_rows64 <- nrow(result_matrix64)
num_cols64 <- ncol(result_matrix64)
value_to_repeat64 <- result_matrix64[2, ]
start_row64 <- 3
start_col64 <- 2
end_row64 <- num_rows64
end_col64 <- num_cols64
for (i in start_row64:end_row64) {
  for (j in start_col64:end_col64) {
    if (j >= (i - start_row64 + start_col64)) {
      pattern_index64 <- (j - (i - start_row64 + start_col64)) %% 32 + 1
      result_matrix64[i, j] <- value_to_repeat64[pattern_index64]
    }
  }
}
print(result_matrix64)
extracted_values64 <- as.numeric(result_matrix64[2:17, 16:33])
extracted_values64 <- extracted_values64[-1]
devs_fcohort_all64 <- c(devs_fcohort_past, extracted_values64)
devs_fyear_all64 <- c(devs_fyear_past, devs_fyear_fut64)
ewaa_m3$fyear_devs64 <- devs_fyear_all64
ewaa_m3$fcohort_devs64 <- devs_fcohort_all64
ewaa_m3$ewaa_fyear64 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fyear_devs64)
ewaa_m3$ewaa_fcohort64 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fcohort_devs64)
ewaa_m3$ewaa_fyear_fcoh64 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fyear_devs64 + ewaa_m3$fcohort_devs64)
head(ewaa_m3, n = 10)

# Selecciona el peso predicho con el efecto de año y cohorte solamente.
ewaa_m3_sim64 <- ewaa_m3[, c("real_year", "age", "ewaa_fyear_fcoh64")]
head(ewaa_m3_sim64, n = 800)

# Renombra las columnas al formato SS3 ewaa.ss
ewaa_m3_sim64 <- ewaa_m3_sim64 %>%
  rename(year = real_year, pred_weight = ewaa_fyear_fcoh64)
colnames(ewaa_m3_sim64) <- c("year", "age", "pred_weight")
head(ewaa_m3_sim64, n = 800)

# Convierte la columna de año a numérica si es necesario
ewaa_m3_sim64$year <- as.numeric(ewaa_m3_sim64$year)

# Ejecuta la función pad_weight_at_age e inspecciona la salida
padded_data64 <- pad_weight_at_age(ewaa_m3_sim64)
str(padded_data64)

# Escribe el archivo con el formato SS3
hakedataUSA:::write_wtatage_file(file = "ewaa_sim64.ss", data = padded_data64, maturity = maturity_at_age)

####################################################
# Sim 65
set.seed(659)
devs_fyear65 <- rnorm(nT, mean = 0, sd = sd_fyear)
devs_fyear_fut65 <- rep(devs_fyear65, each = 16)
devs_fcohort65 <- rnorm(nT_23, mean = 0, sd = sd_fcohort)
devs_fcohort_past_future65 <- c(fcohort_estimate_values, devs_fcohort65)
estimate_fcohort_all65 <- data.frame(year = years_cohort, estimate = devs_fcohort_past_future65)
result_matrix65 <- matrix(0, nrow = 17, ncol = num_cols)
result_matrix65[1, ] <- format(as.integer(estimate_fcohort_all65$year), nsmall = 0)
result_matrix65[2, ] <- estimate_fcohort_all65$estimate
print(result_matrix65)
num_rows65 <- nrow(result_matrix65)
num_cols65 <- ncol(result_matrix65)
value_to_repeat65 <- result_matrix65[2, ]
start_row65 <- 3
start_col65 <- 2
end_row65 <- num_rows65
end_col65 <- num_cols65
for (i in start_row65:end_row65) {
  for (j in start_col65:end_col65) {
    if (j >= (i - start_row65 + start_col65)) {
      pattern_index65 <- (j - (i - start_row65 + start_col65)) %% 32 + 1
      result_matrix65[i, j] <- value_to_repeat65[pattern_index65]
    }
  }
}
print(result_matrix65)
extracted_values65 <- as.numeric(result_matrix65[2:17, 16:33])
extracted_values65 <- extracted_values65[-1]
devs_fcohort_all65 <- c(devs_fcohort_past, extracted_values65)
devs_fyear_all65 <- c(devs_fyear_past, devs_fyear_fut65)
ewaa_m3$fyear_devs65 <- devs_fyear_all65
ewaa_m3$fcohort_devs65 <- devs_fcohort_all65
ewaa_m3$ewaa_fyear65 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fyear_devs65)
ewaa_m3$ewaa_fcohort65 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fcohort_devs65)
ewaa_m3$ewaa_fyear_fcoh65 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fyear_devs65 + ewaa_m3$fcohort_devs65)
head(ewaa_m3, n = 10)

# Selecciona el peso predicho con el efecto de año y cohorte solamente.
ewaa_m3_sim65 <- ewaa_m3[, c("real_year", "age", "ewaa_fyear_fcoh65")]
head(ewaa_m3_sim65, n = 800)

# Renombra las columnas al formato SS3 ewaa.ss
ewaa_m3_sim65 <- ewaa_m3_sim65 %>%
  rename(year = real_year, pred_weight = ewaa_fyear_fcoh65)
colnames(ewaa_m3_sim65) <- c("year", "age", "pred_weight")
head(ewaa_m3_sim65, n = 800)

# Convierte la columna de año a numérica si es necesario
ewaa_m3_sim65$year <- as.numeric(ewaa_m3_sim65$year)

# Ejecuta la función pad_weight_at_age e inspecciona la salida
padded_data65 <- pad_weight_at_age(ewaa_m3_sim65)
str(padded_data65)

# Escribe el archivo con el formato SS3
hakedataUSA:::write_wtatage_file(file = "ewaa_sim65.ss", data = padded_data65, maturity = maturity_at_age)

#####################################################################
# Sim 66
set.seed(661)
devs_fyear66 <- rnorm(nT, mean = 0, sd = sd_fyear)
devs_fyear_fut66 <- rep(devs_fyear66, each = 16)
devs_fcohort66 <- rnorm(nT_23, mean = 0, sd = sd_fcohort)
devs_fcohort_past_future66 <- c(fcohort_estimate_values, devs_fcohort66)
estimate_fcohort_all66 <- data.frame(year = years_cohort, estimate = devs_fcohort_past_future66)
result_matrix66 <- matrix(0, nrow = 17, ncol = num_cols)
result_matrix66[1, ] <- format(as.integer(estimate_fcohort_all66$year), nsmall = 0)
result_matrix66[2, ] <- estimate_fcohort_all66$estimate
print(result_matrix66)
num_rows66 <- nrow(result_matrix66)
num_cols66 <- ncol(result_matrix66)
value_to_repeat66 <- result_matrix66[2, ]
start_row66 <- 3
start_col66 <- 2
end_row66 <- num_rows66
end_col66 <- num_cols66
for (i in start_row66:end_row66) {
  for (j in start_col66:end_col66) {
    if (j >= (i - start_row66 + start_col66)) {
      pattern_index66 <- (j - (i - start_row66 + start_col66)) %% 32 + 1
      result_matrix66[i, j] <- value_to_repeat66[pattern_index66]
    }
  }
}
print(result_matrix66)
extracted_values66 <- as.numeric(result_matrix66[2:17, 16:33])
extracted_values66 <- extracted_values66[-1]
devs_fcohort_all66 <- c(devs_fcohort_past, extracted_values66)
devs_fyear_all66 <- c(devs_fyear_past, devs_fyear_fut66)
ewaa_m3$fyear_devs66 <- devs_fyear_all66
ewaa_m3$fcohort_devs66 <- devs_fcohort_all66
ewaa_m3$ewaa_fyear66 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fyear_devs66)
ewaa_m3$ewaa_fcohort66 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fcohort_devs66)
ewaa_m3$ewaa_fyear_fcoh66 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fyear_devs66 + ewaa_m3$fcohort_devs66)
head(ewaa_m3, n = 10)

# Selecciona el peso predicho con el efecto de año y cohorte solamente.
ewaa_m3_sim66 <- ewaa_m3[, c("real_year", "age", "ewaa_fyear_fcoh66")]
head(ewaa_m3_sim66, n = 800)

# Renombra las columnas al formato SS3 ewaa.ss
ewaa_m3_sim66 <- ewaa_m3_sim66 %>%
  rename(year = real_year, pred_weight = ewaa_fyear_fcoh66)
colnames(ewaa_m3_sim66) <- c("year", "age", "pred_weight")
head(ewaa_m3_sim66, n = 800)

# Convierte la columna de año a numérica si es necesario
ewaa_m3_sim66$year <- as.numeric(ewaa_m3_sim66$year)

# Ejecuta la función pad_weight_at_age e inspecciona la salida
padded_data66 <- pad_weight_at_age(ewaa_m3_sim66)
str(padded_data66)

# Escribe el archivo con el formato SS3
hakedataUSA:::write_wtatage_file(file = "ewaa_sim66.ss", data = padded_data66, maturity = maturity_at_age)

#####################################################################
# Sim 67
set.seed(673)
devs_fyear67 <- rnorm(nT, mean = 0, sd = sd_fyear)
devs_fyear_fut67 <- rep(devs_fyear67, each = 16)
devs_fcohort67 <- rnorm(nT_23, mean = 0, sd = sd_fcohort)
devs_fcohort_past_future67 <- c(fcohort_estimate_values, devs_fcohort67)
estimate_fcohort_all67 <- data.frame(year = years_cohort, estimate = devs_fcohort_past_future67)
result_matrix67 <- matrix(0, nrow = 17, ncol = num_cols)
result_matrix67[1, ] <- format(as.integer(estimate_fcohort_all67$year), nsmall = 0)
result_matrix67[2, ] <- estimate_fcohort_all67$estimate
print(result_matrix67)
num_rows67 <- nrow(result_matrix67)
num_cols67 <- ncol(result_matrix67)
value_to_repeat67 <- result_matrix67[2, ]
start_row67 <- 3
start_col67 <- 2
end_row67 <- num_rows67
end_col67 <- num_cols67
for (i in start_row67:end_row67) {
  for (j in start_col67:end_col67) {
    if (j >= (i - start_row67 + start_col67)) {
      pattern_index67 <- (j - (i - start_row67 + start_col67)) %% 32 + 1
      result_matrix67[i, j] <- value_to_repeat67[pattern_index67]
    }
  }
}
print(result_matrix67)
extracted_values67 <- as.numeric(result_matrix67[2:17, 16:33])
extracted_values67 <- extracted_values67[-1]
devs_fcohort_all67 <- c(devs_fcohort_past, extracted_values67)
devs_fyear_all67 <- c(devs_fyear_past, devs_fyear_fut67)
ewaa_m3$fyear_devs67 <- devs_fyear_all67
ewaa_m3$fcohort_devs67 <- devs_fcohort_all67
ewaa_m3$ewaa_fyear67 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fyear_devs67)
ewaa_m3$ewaa_fcohort67 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fcohort_devs67)
ewaa_m3$ewaa_fyear_fcoh67 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fyear_devs67 + ewaa_m3$fcohort_devs67)
head(ewaa_m3, n = 10)

# Selecciona el peso predicho con el efecto de año y cohorte solamente.
ewaa_m3_sim67 <- ewaa_m3[, c("real_year", "age", "ewaa_fyear_fcoh67")]
head(ewaa_m3_sim67, n = 800)

# Renombra las columnas al formato SS3 ewaa.ss
ewaa_m3_sim67 <- ewaa_m3_sim67 %>%
  rename(year = real_year, pred_weight = ewaa_fyear_fcoh67)
colnames(ewaa_m3_sim67) <- c("year", "age", "pred_weight")
head(ewaa_m3_sim67, n = 800)

# Convierte la columna de año a numérica si es necesario
ewaa_m3_sim67$year <- as.numeric(ewaa_m3_sim67$year)

# Ejecuta la función pad_weight_at_age e inspecciona la salida
padded_data67 <- pad_weight_at_age(ewaa_m3_sim67)
str(padded_data67)

# Escribe el archivo con el formato SS3
hakedataUSA:::write_wtatage_file(file = "ewaa_sim67.ss", data = padded_data67, maturity = maturity_at_age)

###########################################################################
# Sim 68
set.seed(683)
devs_fyear68 <- rnorm(nT, mean = 0, sd = sd_fyear)
devs_fyear_fut68 <- rep(devs_fyear68, each = 16)
devs_fcohort68 <- rnorm(nT_23, mean = 0, sd = sd_fcohort)
devs_fcohort_past_future68 <- c(fcohort_estimate_values, devs_fcohort68)
estimate_fcohort_all68 <- data.frame(year = years_cohort, estimate = devs_fcohort_past_future68)
result_matrix68 <- matrix(0, nrow = 17, ncol = num_cols)
result_matrix68[1, ] <- format(as.integer(estimate_fcohort_all68$year), nsmall = 0)
result_matrix68[2, ] <- estimate_fcohort_all68$estimate
print(result_matrix68)
num_rows68 <- nrow(result_matrix68)
num_cols68 <- ncol(result_matrix68)
value_to_repeat68 <- result_matrix68[2, ]
start_row68 <- 3
start_col68 <- 2
end_row68 <- num_rows68
end_col68 <- num_cols68
for (i in start_row68:end_row68) {
  for (j in start_col68:end_col68) {
    if (j >= (i - start_row68 + start_col68)) {
      pattern_index68 <- (j - (i - start_row68 + start_col68)) %% 32 + 1
      result_matrix68[i, j] <- value_to_repeat68[pattern_index68]
    }
  }
}
print(result_matrix68)
extracted_values68 <- as.numeric(result_matrix68[2:17, 16:33])
extracted_values68 <- extracted_values68[-1]
devs_fcohort_all68 <- c(devs_fcohort_past, extracted_values68)
devs_fyear_all68 <- c(devs_fyear_past, devs_fyear_fut68)
ewaa_m3$fyear_devs68 <- devs_fyear_all68
ewaa_m3$fcohort_devs68 <- devs_fcohort_all68
ewaa_m3$ewaa_fyear68 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fyear_devs68)
ewaa_m3$ewaa_fcohort68 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fcohort_devs68)
ewaa_m3$ewaa_fyear_fcoh68 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fyear_devs68 + ewaa_m3$fcohort_devs68)
head(ewaa_m3, n = 10)

# Selecciona el peso predicho con el efecto de año y cohorte solamente.
ewaa_m3_sim68 <- ewaa_m3[, c("real_year", "age", "ewaa_fyear_fcoh68")]
head(ewaa_m3_sim68, n = 800)

# Renombra las columnas al formato SS3 ewaa.ss
ewaa_m3_sim68 <- ewaa_m3_sim68 %>%
  rename(year = real_year, pred_weight = ewaa_fyear_fcoh68)
colnames(ewaa_m3_sim68) <- c("year", "age", "pred_weight")
head(ewaa_m3_sim68, n = 800)

# Convierte la columna de año a numérica si es necesario
ewaa_m3_sim68$year <- as.numeric(ewaa_m3_sim68$year)

# Ejecuta la función pad_weight_at_age e inspecciona la salida
padded_data68 <- pad_weight_at_age(ewaa_m3_sim68)
str(padded_data68)

# Escribe el archivo con el formato SS3
hakedataUSA:::write_wtatage_file(file = "ewaa_sim68.ss", data = padded_data68, maturity = maturity_at_age)

##########################################################################
# Sim 69
set.seed(691)
devs_fyear69 <- rnorm(nT, mean = 0, sd = sd_fyear)
devs_fyear_fut69 <- rep(devs_fyear69, each = 16)
devs_fcohort69 <- rnorm(nT_23, mean = 0, sd = sd_fcohort)
devs_fcohort_past_future69 <- c(fcohort_estimate_values, devs_fcohort69)
estimate_fcohort_all69 <- data.frame(year = years_cohort, estimate = devs_fcohort_past_future69)
result_matrix69 <- matrix(0, nrow = 17, ncol = num_cols)
result_matrix69[1, ] <- format(as.integer(estimate_fcohort_all69$year), nsmall = 0)
result_matrix69[2, ] <- estimate_fcohort_all69$estimate
print(result_matrix69)
num_rows69 <- nrow(result_matrix69)
num_cols69 <- ncol(result_matrix69)
value_to_repeat69 <- result_matrix69[2, ]
start_row69 <- 3
start_col69 <- 2
end_row69 <- num_rows69
end_col69 <- num_cols69
for (i in start_row69:end_row69) {
  for (j in start_col69:end_col69) {
    if (j >= (i - start_row69 + start_col69)) {
      pattern_index69 <- (j - (i - start_row69 + start_col69)) %% 32 + 1
      result_matrix69[i, j] <- value_to_repeat69[pattern_index69]
    }
  }
}
print(result_matrix69)
extracted_values69 <- as.numeric(result_matrix69[2:17, 16:33])
extracted_values69 <- extracted_values69[-1]
devs_fcohort_all69 <- c(devs_fcohort_past, extracted_values69)
devs_fyear_all69 <- c(devs_fyear_past, devs_fyear_fut69)
ewaa_m3$fyear_devs69 <- devs_fyear_all69
ewaa_m3$fcohort_devs69 <- devs_fcohort_all69
ewaa_m3$ewaa_fyear69 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fyear_devs69)
ewaa_m3$ewaa_fcohort69 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fcohort_devs69)
ewaa_m3$ewaa_fyear_fcoh69 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fyear_devs69 + ewaa_m3$fcohort_devs69)
head(ewaa_m3, n = 10)

# Selecciona el peso predicho con el efecto de año y cohorte solamente.
ewaa_m3_sim69 <- ewaa_m3[, c("real_year", "age", "ewaa_fyear_fcoh69")]
head(ewaa_m3_sim69, n = 800)

# Renombra las columnas al formato SS3 ewaa.ss
ewaa_m3_sim69 <- ewaa_m3_sim69 %>%
  rename(year = real_year, pred_weight = ewaa_fyear_fcoh69)
colnames(ewaa_m3_sim69) <- c("year", "age", "pred_weight")
head(ewaa_m3_sim69, n = 800)

# Convierte la columna de año a numérica si es necesario
ewaa_m3_sim69$year <- as.numeric(ewaa_m3_sim69$year)

# Ejecuta la función pad_weight_at_age e inspecciona la salida
padded_data69 <- pad_weight_at_age(ewaa_m3_sim69)
str(padded_data69)

# Escribe el archivo con el formato SS3
hakedataUSA:::write_wtatage_file(file = "ewaa_sim69.ss", data = padded_data69, maturity = maturity_at_age)

####################################################################
# Sim 70
set.seed(701)
devs_fyear70 <- rnorm(nT, mean = 0, sd = sd_fyear)
devs_fyear_fut70 <- rep(devs_fyear70, each = 16)
devs_fcohort70 <- rnorm(nT_23, mean = 0, sd = sd_fcohort)
devs_fcohort_past_future70 <- c(fcohort_estimate_values, devs_fcohort70)
estimate_fcohort_all70 <- data.frame(year = years_cohort, estimate = devs_fcohort_past_future70)
result_matrix70 <- matrix(0, nrow = 17, ncol = num_cols)
result_matrix70[1, ] <- format(as.integer(estimate_fcohort_all70$year), nsmall = 0)
result_matrix70[2, ] <- estimate_fcohort_all70$estimate
print(result_matrix70)
num_rows70 <- nrow(result_matrix70)
num_cols70 <- ncol(result_matrix70)
value_to_repeat70 <- result_matrix70[2, ]
start_row70 <- 3
start_col70 <- 2
end_row70 <- num_rows70
end_col70 <- num_cols70
for (i in start_row70:end_row70) {
  for (j in start_col70:end_col70) {
    if (j >= (i - start_row70 + start_col70)) {
      pattern_index70 <- (j - (i - start_row70 + start_col70)) %% 32 + 1
      result_matrix70[i, j] <- value_to_repeat70[pattern_index70]
    }
  }
}
print(result_matrix70)
extracted_values70 <- as.numeric(result_matrix70[2:17, 16:33])
extracted_values70 <- extracted_values70[-1]
devs_fcohort_all70 <- c(devs_fcohort_past, extracted_values70)
devs_fyear_all70 <- c(devs_fyear_past, devs_fyear_fut70)
ewaa_m3$fyear_devs70 <- devs_fyear_all70
ewaa_m3$fcohort_devs70 <- devs_fcohort_all70
ewaa_m3$ewaa_fyear70 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fyear_devs70)
ewaa_m3$ewaa_fcohort70 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fcohort_devs70)
ewaa_m3$ewaa_fyear_fcoh70 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fyear_devs70 + ewaa_m3$fcohort_devs70)
head(ewaa_m3, n = 10)

# Selecciona el peso predicho con el efecto de año y cohorte solamente.
ewaa_m3_sim70 <- ewaa_m3[, c("real_year", "age", "ewaa_fyear_fcoh70")]
head(ewaa_m3_sim70, n = 800)

# Renombra las columnas al formato SS3 ewaa.ss
ewaa_m3_sim70 <- ewaa_m3_sim70 %>%
  rename(year = real_year, pred_weight = ewaa_fyear_fcoh70)
colnames(ewaa_m3_sim70) <- c("year", "age", "pred_weight")
head(ewaa_m3_sim70, n = 800)

# Convierte la columna de año a numérica si es necesario
ewaa_m3_sim70$year <- as.numeric(ewaa_m3_sim70$year)

# Ejecuta la función pad_weight_at_age e inspecciona la salida
padded_data70 <- pad_weight_at_age(ewaa_m3_sim70)
str(padded_data70)

# Escribe el archivo con el formato SS3
hakedataUSA:::write_wtatage_file(file = "ewaa_sim70.ss", data = padded_data70, maturity = maturity_at_age)


######################################################
# Sim 71
set.seed(719)
devs_fyear71 <- rnorm(nT, mean = 0, sd = sd_fyear)
devs_fyear_fut71 <- rep(devs_fyear71, each = 16)
devs_fcohort71 <- rnorm(nT_23, mean = 0, sd = sd_fcohort)
devs_fcohort_past_future71 <- c(fcohort_estimate_values, devs_fcohort71)
estimate_fcohort_all71 <- data.frame(year = years_cohort, estimate = devs_fcohort_past_future71)
result_matrix71 <- matrix(0, nrow = 17, ncol = num_cols)
result_matrix71[1, ] <- format(as.integer(estimate_fcohort_all71$year), nsmall = 0)
result_matrix71[2, ] <- estimate_fcohort_all71$estimate
print(result_matrix71)
num_rows71 <- nrow(result_matrix71)
num_cols71 <- ncol(result_matrix71)
value_to_repeat71 <- result_matrix71[2, ]
start_row71 <- 3
start_col71 <- 2
end_row71 <- num_rows71
end_col71 <- num_cols71
for (i in start_row71:end_row71) {
  for (j in start_col71:end_col71) {
    if (j >= (i - start_row71 + start_col71)) {
      pattern_index71 <- (j - (i - start_row71 + start_col71)) %% 32 + 1
      result_matrix71[i, j] <- value_to_repeat71[pattern_index71]
    }
  }
}
print(result_matrix71)
extracted_values71 <- as.numeric(result_matrix71[2:17, 16:33])
extracted_values71 <- extracted_values71[-1]
devs_fcohort_all71 <- c(devs_fcohort_past, extracted_values71)
devs_fyear_all71 <- c(devs_fyear_past, devs_fyear_fut71)
ewaa_m3$fyear_devs71 <- devs_fyear_all71
ewaa_m3$fcohort_devs71 <- devs_fcohort_all71
ewaa_m3$ewaa_fyear71 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fyear_devs71)
ewaa_m3$ewaa_fcohort71 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fcohort_devs71)
ewaa_m3$ewaa_fyear_fcoh71 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fyear_devs71 + ewaa_m3$fcohort_devs71)
head(ewaa_m3, n = 10)

# Selecciona el peso predicho con el efecto de año y cohorte solamente.
ewaa_m3_sim71 <- ewaa_m3[, c("real_year", "age", "ewaa_fyear_fcoh71")]
head(ewaa_m3_sim71, n = 800)

# Renombra las columnas al formato SS3 ewaa.ss
ewaa_m3_sim71 <- ewaa_m3_sim71 %>%
  rename(year = real_year, pred_weight = ewaa_fyear_fcoh71)
colnames(ewaa_m3_sim71) <- c("year", "age", "pred_weight")
head(ewaa_m3_sim71, n = 800)

# Convierte la columna de año a numérica si es necesario
ewaa_m3_sim71$year <- as.numeric(ewaa_m3_sim71$year)

# Ejecuta la función pad_weight_at_age e inspecciona la salida
padded_data71 <- pad_weight_at_age(ewaa_m3_sim71)
str(padded_data71)

# Escribe el archivo con el formato SS3
hakedataUSA:::write_wtatage_file(file = "ewaa_sim71.ss", data = padded_data71, maturity = maturity_at_age)

######################################################################
# Sim 72
set.seed(727)
devs_fyear72 <- rnorm(nT, mean = 0, sd = sd_fyear)
devs_fyear_fut72 <- rep(devs_fyear72, each = 16)
devs_fcohort72 <- rnorm(nT_23, mean = 0, sd = sd_fcohort)
devs_fcohort_past_future72 <- c(fcohort_estimate_values, devs_fcohort72)
estimate_fcohort_all72 <- data.frame(year = years_cohort, estimate = devs_fcohort_past_future72)
result_matrix72 <- matrix(0, nrow = 17, ncol = num_cols)
result_matrix72[1, ] <- format(as.integer(estimate_fcohort_all72$year), nsmall = 0)
result_matrix72[2, ] <- estimate_fcohort_all72$estimate
print(result_matrix72)
num_rows72 <- nrow(result_matrix72)
num_cols72 <- ncol(result_matrix72)
value_to_repeat72 <- result_matrix72[2, ]
start_row72 <- 3
start_col72 <- 2
end_row72 <- num_rows72
end_col72 <- num_cols72
for (i in start_row72:end_row72) {
  for (j in start_col72:end_col72) {
    if (j >= (i - start_row72 + start_col72)) {
      pattern_index72 <- (j - (i - start_row72 + start_col72)) %% 32 + 1
      result_matrix72[i, j] <- value_to_repeat72[pattern_index72]
    }
  }
}
print(result_matrix72)
extracted_values72 <- as.numeric(result_matrix72[2:17, 16:33])
extracted_values72 <- extracted_values72[-1]
devs_fcohort_all72 <- c(devs_fcohort_past, extracted_values72)
devs_fyear_all72 <- c(devs_fyear_past, devs_fyear_fut72)
ewaa_m3$fyear_devs72 <- devs_fyear_all72
ewaa_m3$fcohort_devs72 <- devs_fcohort_all72
ewaa_m3$ewaa_fyear72 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fyear_devs72)
ewaa_m3$ewaa_fcohort72 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fcohort_devs72)
ewaa_m3$ewaa_fyear_fcoh72 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fyear_devs72 + ewaa_m3$fcohort_devs72)
head(ewaa_m3, n = 10)

# Selecciona el peso predicho con el efecto de año y cohorte solamente.
ewaa_m3_sim72 <- ewaa_m3[, c("real_year", "age", "ewaa_fyear_fcoh72")]
head(ewaa_m3_sim72, n = 800)

# Renombra las columnas al formato SS3 ewaa.ss
ewaa_m3_sim72 <- ewaa_m3_sim72 %>%
  rename(year = real_year, pred_weight = ewaa_fyear_fcoh72)
colnames(ewaa_m3_sim72) <- c("year", "age", "pred_weight")
head(ewaa_m3_sim72, n = 800)

# Convierte la columna de año a numérica si es necesario
ewaa_m3_sim72$year <- as.numeric(ewaa_m3_sim72$year)

# Ejecuta la función pad_weight_at_age e inspecciona la salida
padded_data72 <- pad_weight_at_age(ewaa_m3_sim72)
str(padded_data72)

# Escribe el archivo con el formato SS3
hakedataUSA:::write_wtatage_file(file = "ewaa_sim72.ss", data = padded_data72, maturity = maturity_at_age)
######################################################################
# Sim 73
set.seed(733)
devs_fyear73 <- rnorm(nT, mean = 0, sd = sd_fyear)
devs_fyear_fut73 <- rep(devs_fyear73, each = 16)
devs_fcohort73 <- rnorm(nT_23, mean = 0, sd = sd_fcohort)
devs_fcohort_past_future73 <- c(fcohort_estimate_values, devs_fcohort73)
estimate_fcohort_all73 <- data.frame(year = years_cohort, estimate = devs_fcohort_past_future73)
result_matrix73 <- matrix(0, nrow = 17, ncol = num_cols)
result_matrix73[1, ] <- format(as.integer(estimate_fcohort_all73$year), nsmall = 0)
result_matrix73[2, ] <- estimate_fcohort_all73$estimate
print(result_matrix73)
num_rows73 <- nrow(result_matrix73)
num_cols73 <- ncol(result_matrix73)
value_to_repeat73 <- result_matrix73[2, ]
start_row73 <- 3
start_col73 <- 2
end_row73 <- num_rows73
end_col73 <- num_cols73
for (i in start_row73:end_row73) {
  for (j in start_col73:end_col73) {
    if (j >= (i - start_row73 + start_col73)) {
      pattern_index73 <- (j - (i - start_row73 + start_col73)) %% 32 + 1
      result_matrix73[i, j] <- value_to_repeat73[pattern_index73]
    }
  }
}
print(result_matrix73)
extracted_values73 <- as.numeric(result_matrix73[2:17, 16:33])
extracted_values73 <- extracted_values73[-1]
devs_fcohort_all73 <- c(devs_fcohort_past, extracted_values73)
devs_fyear_all73 <- c(devs_fyear_past, devs_fyear_fut73)
ewaa_m3$fyear_devs73 <- devs_fyear_all73
ewaa_m3$fcohort_devs73 <- devs_fcohort_all73
ewaa_m3$ewaa_fyear73 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fyear_devs73)
ewaa_m3$ewaa_fcohort73 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fcohort_devs73)
ewaa_m3$ewaa_fyear_fcoh73 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fyear_devs73 + ewaa_m3$fcohort_devs73)
head(ewaa_m3, n = 10)

# Selecciona el peso predicho con el efecto de año y cohorte solamente.
ewaa_m3_sim73 <- ewaa_m3[, c("real_year", "age", "ewaa_fyear_fcoh73")]
head(ewaa_m3_sim73, n = 800)

# Renombra las columnas al formato SS3 ewaa.ss
ewaa_m3_sim73 <- ewaa_m3_sim73 %>%
  rename(year = real_year, pred_weight = ewaa_fyear_fcoh73)
colnames(ewaa_m3_sim73) <- c("year", "age", "pred_weight")
head(ewaa_m3_sim73, n = 800)

# Convierte la columna de año a numérica si es necesario
ewaa_m3_sim73$year <- as.numeric(ewaa_m3_sim73$year)

# Ejecuta la función pad_weight_at_age e inspecciona la salida
padded_data73 <- pad_weight_at_age(ewaa_m3_sim73)
str(padded_data73)

# Escribe el archivo con el formato SS3
hakedataUSA:::write_wtatage_file(file = "ewaa_sim73.ss", data = padded_data73, maturity = maturity_at_age)

###########################################################################
# Sim 74
set.seed(743)
devs_fyear74 <- rnorm(nT, mean = 0, sd = sd_fyear)
devs_fyear_fut74 <- rep(devs_fyear74, each = 16)
devs_fcohort74 <- rnorm(nT_23, mean = 0, sd = sd_fcohort)
devs_fcohort_past_future74 <- c(fcohort_estimate_values, devs_fcohort74)
estimate_fcohort_all74 <- data.frame(year = years_cohort, estimate = devs_fcohort_past_future74)
result_matrix74 <- matrix(0, nrow = 17, ncol = num_cols)
result_matrix74[1, ] <- format(as.integer(estimate_fcohort_all74$year), nsmall = 0)
result_matrix74[2, ] <- estimate_fcohort_all74$estimate
print(result_matrix74)
num_rows74 <- nrow(result_matrix74)
num_cols74 <- ncol(result_matrix74)
value_to_repeat74 <- result_matrix74[2, ]
start_row74 <- 3
start_col74 <- 2
end_row74 <- num_rows74
end_col74 <- num_cols74
for (i in start_row74:end_row74) {
  for (j in start_col74:end_col74) {
    if (j >= (i - start_row74 + start_col74)) {
      pattern_index74 <- (j - (i - start_row74 + start_col74)) %% 32 + 1
      result_matrix74[i, j] <- value_to_repeat74[pattern_index74]
    }
  }
}
print(result_matrix74)
extracted_values74 <- as.numeric(result_matrix74[2:17, 16:33])
extracted_values74 <- extracted_values74[-1]
devs_fcohort_all74 <- c(devs_fcohort_past, extracted_values74)
devs_fyear_all74 <- c(devs_fyear_past, devs_fyear_fut74)
ewaa_m3$fyear_devs74 <- devs_fyear_all74
ewaa_m3$fcohort_devs74 <- devs_fcohort_all74
ewaa_m3$ewaa_fyear74 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fyear_devs74)
ewaa_m3$ewaa_fcohort74 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fcohort_devs74)
ewaa_m3$ewaa_fyear_fcoh74 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fyear_devs74 + ewaa_m3$fcohort_devs74)
head(ewaa_m3, n = 10)

# Selecciona el peso predicho con el efecto de año y cohorte solamente.
ewaa_m3_sim74 <- ewaa_m3[, c("real_year", "age", "ewaa_fyear_fcoh74")]
head(ewaa_m3_sim74, n = 800)

# Renombra las columnas al formato SS3 ewaa.ss
ewaa_m3_sim74 <- ewaa_m3_sim74 %>%
  rename(year = real_year, pred_weight = ewaa_fyear_fcoh74)
colnames(ewaa_m3_sim74) <- c("year", "age", "pred_weight")
head(ewaa_m3_sim74, n = 800)

# Convierte la columna de año a numérica si es necesario
ewaa_m3_sim74$year <- as.numeric(ewaa_m3_sim74$year)

# Ejecuta la función pad_weight_at_age e inspecciona la salida
padded_data74 <- pad_weight_at_age(ewaa_m3_sim74)
str(padded_data74)

# Escribe el archivo con el formato SS3
hakedataUSA:::write_wtatage_file(file = "ewaa_sim74.ss", data = padded_data74, maturity = maturity_at_age)

#######################################################################################################
# Sim 75
set.seed(751)
devs_fyear75 <- rnorm(nT, mean = 0, sd = sd_fyear)
devs_fyear_fut75 <- rep(devs_fyear75, each = 16)
devs_fcohort75 <- rnorm(nT_23, mean = 0, sd = sd_fcohort)
devs_fcohort_past_future75 <- c(fcohort_estimate_values, devs_fcohort75)
estimate_fcohort_all75 <- data.frame(year = years_cohort, estimate = devs_fcohort_past_future75)
result_matrix75 <- matrix(0, nrow = 17, ncol = num_cols)
result_matrix75[1, ] <- format(as.integer(estimate_fcohort_all75$year), nsmall = 0)
result_matrix75[2, ] <- estimate_fcohort_all75$estimate
print(result_matrix75)
num_rows75 <- nrow(result_matrix75)
num_cols75 <- ncol(result_matrix75)
value_to_repeat75 <- result_matrix75[2, ]
start_row75 <- 3
start_col75 <- 2
end_row75 <- num_rows75
end_col75 <- num_cols75
for (i in start_row75:end_row75) {
  for (j in start_col75:end_col75) {
    if (j >= (i - start_row75 + start_col75)) {
      pattern_index75 <- (j - (i - start_row75 + start_col75)) %% 32 + 1
      result_matrix75[i, j] <- value_to_repeat75[pattern_index75]
    }
  }
}
print(result_matrix75)
extracted_values75 <- as.numeric(result_matrix75[2:17, 16:33])
extracted_values75 <- extracted_values75[-1]
devs_fcohort_all75 <- c(devs_fcohort_past, extracted_values75)
devs_fyear_all75 <- c(devs_fyear_past, devs_fyear_fut75)
ewaa_m3$fyear_devs75 <- devs_fyear_all75
ewaa_m3$fcohort_devs75 <- devs_fcohort_all75
ewaa_m3$ewaa_fyear75 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fyear_devs75)
ewaa_m3$ewaa_fcohort75 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fcohort_devs75)
ewaa_m3$ewaa_fyear_fcoh75 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fyear_devs75 + ewaa_m3$fcohort_devs75)
head(ewaa_m3, n = 10)

# Selecciona el peso predicho con el efecto de año y cohorte solamente.
ewaa_m3_sim75 <- ewaa_m3[, c("real_year", "age", "ewaa_fyear_fcoh75")]
head(ewaa_m3_sim75, n = 800)

# Renombra las columnas al formato SS3 ewaa.ss
ewaa_m3_sim75 <- ewaa_m3_sim75 %>%
  rename(year = real_year, pred_weight = ewaa_fyear_fcoh75)
colnames(ewaa_m3_sim75) <- c("year", "age", "pred_weight")
head(ewaa_m3_sim75, n = 800)

# Convierte la columna de año a numérica si es necesario
ewaa_m3_sim75$year <- as.numeric(ewaa_m3_sim75$year)

# Ejecuta la función pad_weight_at_age e inspecciona la salida
padded_data75 <- pad_weight_at_age(ewaa_m3_sim75)
str(padded_data75)

# Escribe el archivo con el formato SS3
hakedataUSA:::write_wtatage_file(file = "ewaa_sim75.ss", data = padded_data75, maturity = maturity_at_age)

##############################################################################
# Sim 76
set.seed(761)
devs_fyear76 <- rnorm(nT, mean = 0, sd = sd_fyear)
devs_fyear_fut76 <- rep(devs_fyear76, each = 16)
devs_fcohort76 <- rnorm(nT_23, mean = 0, sd = sd_fcohort)
devs_fcohort_past_future76 <- c(fcohort_estimate_values, devs_fcohort76)
estimate_fcohort_all76 <- data.frame(year = years_cohort, estimate = devs_fcohort_past_future76)
result_matrix76 <- matrix(0, nrow = 17, ncol = num_cols)
result_matrix76[1, ] <- format(as.integer(estimate_fcohort_all76$year), nsmall = 0)
result_matrix76[2, ] <- estimate_fcohort_all76$estimate
print(result_matrix76)
num_rows76 <- nrow(result_matrix76)
num_cols76 <- ncol(result_matrix76)
value_to_repeat76 <- result_matrix76[2, ]
start_row76 <- 3
start_col76 <- 2
end_row76 <- num_rows76
end_col76 <- num_cols76
for (i in start_row76:end_row76) {
  for (j in start_col76:end_col76) {
    if (j >= (i - start_row76 + start_col76)) {
      pattern_index76 <- (j - (i - start_row76 + start_col76)) %% 32 + 1
      result_matrix76[i, j] <- value_to_repeat76[pattern_index76]
    }
  }
}
print(result_matrix76)
extracted_values76 <- as.numeric(result_matrix76[2:17, 16:33])
extracted_values76 <- extracted_values76[-1]
devs_fcohort_all76 <- c(devs_fcohort_past, extracted_values76)
devs_fyear_all76 <- c(devs_fyear_past, devs_fyear_fut76)
ewaa_m3$fyear_devs76 <- devs_fyear_all76
ewaa_m3$fcohort_devs76 <- devs_fcohort_all76
ewaa_m3$ewaa_fyear76 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fyear_devs76)
ewaa_m3$ewaa_fcohort76 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fcohort_devs76)
ewaa_m3$ewaa_fyear_fcoh76 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fyear_devs76 + ewaa_m3$fcohort_devs76)
head(ewaa_m3, n = 10)

# Selecciona el peso predicho con el efecto de año y cohorte solamente.
ewaa_m3_sim76 <- ewaa_m3[, c("real_year", "age", "ewaa_fyear_fcoh76")]
head(ewaa_m3_sim76, n = 800)

# Renombra las columnas al formato SS3 ewaa.ss
ewaa_m3_sim76 <- ewaa_m3_sim76 %>%
  rename(year = real_year, pred_weight = ewaa_fyear_fcoh76)
colnames(ewaa_m3_sim76) <- c("year", "age", "pred_weight")
head(ewaa_m3_sim76, n = 800)

# Convierte la columna de año a numérica si es necesario
ewaa_m3_sim76$year <- as.numeric(ewaa_m3_sim76$year)

# Ejecuta la función pad_weight_at_age e inspecciona la salida
padded_data76 <- pad_weight_at_age(ewaa_m3_sim76)
str(padded_data76)

# Escribe el archivo con el formato SS3
hakedataUSA:::write_wtatage_file(file = "ewaa_sim76.ss", data = padded_data76, maturity = maturity_at_age)

#####################################################################
# Sim 77
set.seed(773)
devs_fyear77 <- rnorm(nT, mean = 0, sd = sd_fyear)
devs_fyear_fut77 <- rep(devs_fyear77, each = 16)
devs_fcohort77 <- rnorm(nT_23, mean = 0, sd = sd_fcohort)
devs_fcohort_past_future77 <- c(fcohort_estimate_values, devs_fcohort77)
estimate_fcohort_all77 <- data.frame(year = years_cohort, estimate = devs_fcohort_past_future77)
result_matrix77 <- matrix(0, nrow = 17, ncol = num_cols)
result_matrix77[1, ] <- format(as.integer(estimate_fcohort_all77$year), nsmall = 0)
result_matrix77[2, ] <- estimate_fcohort_all77$estimate
print(result_matrix77)
num_rows77 <- nrow(result_matrix77)
num_cols77 <- ncol(result_matrix77)
value_to_repeat77 <- result_matrix77[2, ]
start_row77 <- 3
start_col77 <- 2
end_row77 <- num_rows77
end_col77 <- num_cols77
for (i in start_row77:end_row77) {
  for (j in start_col77:end_col77) {
    if (j >= (i - start_row77 + start_col77)) {
      pattern_index77 <- (j - (i - start_row77 + start_col77)) %% 32 + 1
      result_matrix77[i, j] <- value_to_repeat77[pattern_index77]
    }
  }
}
print(result_matrix77)
extracted_values77 <- as.numeric(result_matrix77[2:17, 16:33])
extracted_values77 <- extracted_values77[-1]
devs_fcohort_all77 <- c(devs_fcohort_past, extracted_values77)
devs_fyear_all77 <- c(devs_fyear_past, devs_fyear_fut77)
ewaa_m3$fyear_devs77 <- devs_fyear_all77
ewaa_m3$fcohort_devs77 <- devs_fcohort_all77
ewaa_m3$ewaa_fyear77 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fyear_devs77)
ewaa_m3$ewaa_fcohort77 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fcohort_devs77)
ewaa_m3$ewaa_fyear_fcoh77 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fyear_devs77 + ewaa_m3$fcohort_devs77)
head(ewaa_m3, n = 10)

# Selecciona el peso predicho con el efecto de año y cohorte solamente.
ewaa_m3_sim77 <- ewaa_m3[, c("real_year", "age", "ewaa_fyear_fcoh77")]
head(ewaa_m3_sim77, n = 800)

# Renombra las columnas al formato SS3 ewaa.ss
ewaa_m3_sim77 <- ewaa_m3_sim77 %>%
  rename(year = real_year, pred_weight = ewaa_fyear_fcoh77)
colnames(ewaa_m3_sim77) <- c("year", "age", "pred_weight")
head(ewaa_m3_sim77, n = 800)

# Convierte la columna de año a numérica si es necesario
ewaa_m3_sim77$year <- as.numeric(ewaa_m3_sim77$year)

# Ejecuta la función pad_weight_at_age e inspecciona la salida
padded_data77 <- pad_weight_at_age(ewaa_m3_sim77)
str(padded_data77)

# Escribe el archivo con el formato SS3
hakedataUSA:::write_wtatage_file(file = "ewaa_sim77.ss", data = padded_data77, maturity = maturity_at_age)

############################################################
# Sim 78
set.seed(787)
devs_fyear78 <- rnorm(nT, mean = 0, sd = sd_fyear)
devs_fyear_fut78 <- rep(devs_fyear78, each = 16)
devs_fcohort78 <- rnorm(nT_23, mean = 0, sd = sd_fcohort)
devs_fcohort_past_future78 <- c(fcohort_estimate_values, devs_fcohort78)
estimate_fcohort_all78 <- data.frame(year = years_cohort, estimate = devs_fcohort_past_future78)
result_matrix78 <- matrix(0, nrow = 17, ncol = num_cols)
result_matrix78[1, ] <- format(as.integer(estimate_fcohort_all78$year), nsmall = 0)
result_matrix78[2, ] <- estimate_fcohort_all78$estimate
print(result_matrix78)
num_rows78 <- nrow(result_matrix78)
num_cols78 <- ncol(result_matrix78)
value_to_repeat78 <- result_matrix78[2, ]
start_row78 <- 3
start_col78 <- 2
end_row78 <- num_rows78
end_col78 <- num_cols78
for (i in start_row78:end_row78) {
  for (j in start_col78:end_col78) {
    if (j >= (i - start_row78 + start_col78)) {
      pattern_index78 <- (j - (i - start_row78 + start_col78)) %% 32 + 1
      result_matrix78[i, j] <- value_to_repeat78[pattern_index78]
    }
  }
}
print(result_matrix78)
extracted_values78 <- as.numeric(result_matrix78[2:17, 16:33])
extracted_values78 <- extracted_values78[-1]
devs_fcohort_all78 <- c(devs_fcohort_past, extracted_values78)
devs_fyear_all78 <- c(devs_fyear_past, devs_fyear_fut78)
ewaa_m3$fyear_devs78 <- devs_fyear_all78
ewaa_m3$fcohort_devs78 <- devs_fcohort_all78
ewaa_m3$ewaa_fyear78 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fyear_devs78)
ewaa_m3$ewaa_fcohort78 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fcohort_devs78)
ewaa_m3$ewaa_fyear_fcoh78 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fyear_devs78 + ewaa_m3$fcohort_devs78)
head(ewaa_m3, n = 10)

# Selecciona el peso predicho con el efecto de año y cohorte solamente.
ewaa_m3_sim78 <- ewaa_m3[, c("real_year", "age", "ewaa_fyear_fcoh78")]
head(ewaa_m3_sim78, n = 800)

# Renombra las columnas al formato SS3 ewaa.ss
ewaa_m3_sim78 <- ewaa_m3_sim78 %>%
  rename(year = real_year, pred_weight = ewaa_fyear_fcoh78)
colnames(ewaa_m3_sim78) <- c("year", "age", "pred_weight")
head(ewaa_m3_sim78, n = 800)

# Convierte la columna de año a numérica si es necesario
ewaa_m3_sim78$year <- as.numeric(ewaa_m3_sim78$year)

# Ejecuta la función pad_weight_at_age e inspecciona la salida
padded_data78 <- pad_weight_at_age(ewaa_m3_sim78)
str(padded_data78)

# Escribe el archivo con el formato SS3
hakedataUSA:::write_wtatage_file(file = "ewaa_sim78.ss", data = padded_data78, maturity = maturity_at_age)

# Sim 79
set.seed(797)
devs_fyear79 <- rnorm(nT, mean = 0, sd = sd_fyear)
devs_fyear_fut79 <- rep(devs_fyear79, each = 16)
devs_fcohort79 <- rnorm(nT_23, mean = 0, sd = sd_fcohort)
devs_fcohort_past_future79 <- c(fcohort_estimate_values, devs_fcohort79)
estimate_fcohort_all79 <- data.frame(year = years_cohort, estimate = devs_fcohort_past_future79)
result_matrix79 <- matrix(0, nrow = 17, ncol = num_cols)
result_matrix79[1, ] <- format(as.integer(estimate_fcohort_all79$year), nsmall = 0)
result_matrix79[2, ] <- estimate_fcohort_all79$estimate
print(result_matrix79)
num_rows79 <- nrow(result_matrix79)
num_cols79 <- ncol(result_matrix79)
value_to_repeat79 <- result_matrix79[2, ]
start_row79 <- 3
start_col79 <- 2
end_row79 <- num_rows79
end_col79 <- num_cols79
for (i in start_row79:end_row79) {
  for (j in start_col79:end_col79) {
    if (j >= (i - start_row79 + start_col79)) {
      pattern_index79 <- (j - (i - start_row79 + start_col79)) %% 32 + 1
      result_matrix79[i, j] <- value_to_repeat79[pattern_index79]
    }
  }
}
print(result_matrix79)
extracted_values79 <- as.numeric(result_matrix79[2:17, 16:33])
extracted_values79 <- extracted_values79[-1]
devs_fcohort_all79 <- c(devs_fcohort_past, extracted_values79)
devs_fyear_all79 <- c(devs_fyear_past, devs_fyear_fut79)
ewaa_m3$fyear_devs79 <- devs_fyear_all79
ewaa_m3$fcohort_devs79 <- devs_fcohort_all79
ewaa_m3$ewaa_fyear79 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fyear_devs79)
ewaa_m3$ewaa_fcohort79 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fcohort_devs79)
ewaa_m3$ewaa_fyear_fcoh79 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fyear_devs79 + ewaa_m3$fcohort_devs79)
head(ewaa_m3, n = 10)

# Selecciona el peso predicho con el efecto de año y cohorte solamente.
ewaa_m3_sim79 <- ewaa_m3[, c("real_year", "age", "ewaa_fyear_fcoh79")]
head(ewaa_m3_sim79, n = 800)

# Renombra las columnas al formato SS3 ewaa.ss
ewaa_m3_sim79 <- ewaa_m3_sim79 %>%
  rename(year = real_year, pred_weight = ewaa_fyear_fcoh79)
colnames(ewaa_m3_sim79) <- c("year", "age", "pred_weight")
head(ewaa_m3_sim79, n = 800)

# Convierte la columna de año a numérica si es necesario
ewaa_m3_sim79$year <- as.numeric(ewaa_m3_sim79$year)

# Ejecuta la función pad_weight_at_age e inspecciona la salida
padded_data79 <- pad_weight_at_age(ewaa_m3_sim79)
str(padded_data79)

# Escribe el archivo con el formato SS3
hakedataUSA:::write_wtatage_file(file = "ewaa_sim79.ss", data = padded_data79, maturity = maturity_at_age)
###############################################################
# Sim 80
set.seed(809)
devs_fyear80 <- rnorm(nT, mean = 0, sd = sd_fyear)
devs_fyear_fut80 <- rep(devs_fyear80, each = 16)
devs_fcohort80 <- rnorm(nT_23, mean = 0, sd = sd_fcohort)
devs_fcohort_past_future80 <- c(fcohort_estimate_values, devs_fcohort80)
estimate_fcohort_all80 <- data.frame(year = years_cohort, estimate = devs_fcohort_past_future80)
result_matrix80 <- matrix(0, nrow = 17, ncol = num_cols)
result_matrix80[1, ] <- format(as.integer(estimate_fcohort_all80$year), nsmall = 0)
result_matrix80[2, ] <- estimate_fcohort_all80$estimate
print(result_matrix80)
num_rows80 <- nrow(result_matrix80)
num_cols80 <- ncol(result_matrix80)
value_to_repeat80 <- result_matrix80[2, ]
start_row80 <- 3
start_col80 <- 2
end_row80 <- num_rows80
end_col80 <- num_cols80
for (i in start_row80:end_row80) {
  for (j in start_col80:end_col80) {
    if (j >= (i - start_row80 + start_col80)) {
      pattern_index80 <- (j - (i - start_row80 + start_col80)) %% 32 + 1
      result_matrix80[i, j] <- value_to_repeat80[pattern_index80]
    }
  }
}
print(result_matrix80)
extracted_values80 <- as.numeric(result_matrix80[2:17, 16:33])
extracted_values80 <- extracted_values80[-1]
devs_fcohort_all80 <- c(devs_fcohort_past, extracted_values80)
devs_fyear_all80 <- c(devs_fyear_past, devs_fyear_fut80)
ewaa_m3$fyear_devs80 <- devs_fyear_all80
ewaa_m3$fcohort_devs80 <- devs_fcohort_all80
ewaa_m3$ewaa_fyear80 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fyear_devs80)
ewaa_m3$ewaa_fcohort80 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fcohort_devs80)
ewaa_m3$ewaa_fyear_fcoh80 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fyear_devs80 + ewaa_m3$fcohort_devs80)
head(ewaa_m3, n = 10)

# Selecciona el peso predicho con el efecto de año y cohorte solamente.
ewaa_m3_sim80 <- ewaa_m3[, c("real_year", "age", "ewaa_fyear_fcoh80")]
head(ewaa_m3_sim80, n = 800)

# Renombra las columnas al formato SS3 ewaa.ss
ewaa_m3_sim80 <- ewaa_m3_sim80 %>%
  rename(year = real_year, pred_weight = ewaa_fyear_fcoh80)
colnames(ewaa_m3_sim80) <- c("year", "age", "pred_weight")
head(ewaa_m3_sim80, n = 800)

# Convierte la columna de año a numérica si es necesario
ewaa_m3_sim80$year <- as.numeric(ewaa_m3_sim80$year)

# Ejecuta la función pad_weight_at_age e inspecciona la salida
padded_data80 <- pad_weight_at_age(ewaa_m3_sim80)
str(padded_data80)

# Escribe el archivo con el formato SS3
hakedataUSA:::write_wtatage_file(file = "ewaa_sim80.ss", data = padded_data80, maturity = maturity_at_age)

#################################################################
# Sim 81
set.seed(811)
devs_fyear81 <- rnorm(nT, mean = 0, sd = sd_fyear)
devs_fyear_fut81 <- rep(devs_fyear81, each = 16)
devs_fcohort81 <- rnorm(nT_23, mean = 0, sd = sd_fcohort)
devs_fcohort_past_future81 <- c(fcohort_estimate_values, devs_fcohort81)
estimate_fcohort_all81 <- data.frame(year = years_cohort, estimate = devs_fcohort_past_future81)
result_matrix81 <- matrix(0, nrow = 17, ncol = num_cols)
result_matrix81[1, ] <- format(as.integer(estimate_fcohort_all81$year), nsmall = 0)
result_matrix81[2, ] <- estimate_fcohort_all81$estimate
print(result_matrix81)
num_rows81 <- nrow(result_matrix81)
num_cols81 <- ncol(result_matrix81)
value_to_repeat81 <- result_matrix81[2, ]
start_row81 <- 3
start_col81 <- 2
end_row81 <- num_rows81
end_col81 <- num_cols81
for (i in start_row81:end_row81) {
  for (j in start_col81:end_col81) {
    if (j >= (i - start_row81 + start_col81)) {
      pattern_index81 <- (j - (i - start_row81 + start_col81)) %% 32 + 1
      result_matrix81[i, j] <- value_to_repeat81[pattern_index81]
    }
  }
}
print(result_matrix81)
extracted_values81 <- as.numeric(result_matrix81[2:17, 16:33])
extracted_values81 <- extracted_values81[-1]
devs_fcohort_all81 <- c(devs_fcohort_past, extracted_values81)
devs_fyear_all81 <- c(devs_fyear_past, devs_fyear_fut81)
ewaa_m3$fyear_devs81 <- devs_fyear_all81
ewaa_m3$fcohort_devs81 <- devs_fcohort_all81
ewaa_m3$ewaa_fyear81 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fyear_devs81)
ewaa_m3$ewaa_fcohort81 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fcohort_devs81)
ewaa_m3$ewaa_fyear_fcoh81 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fyear_devs81 + ewaa_m3$fcohort_devs81)
head(ewaa_m3, n = 10)

# Selecciona el peso predicho con el efecto de año y cohorte solamente.
ewaa_m3_sim81 <- ewaa_m3[, c("real_year", "age", "ewaa_fyear_fcoh81")]
head(ewaa_m3_sim81, n = 800)

# Renombra las columnas al formato SS3 ewaa.ss
ewaa_m3_sim81 <- ewaa_m3_sim81 %>%
  rename(year = real_year, pred_weight = ewaa_fyear_fcoh81)
colnames(ewaa_m3_sim81) <- c("year", "age", "pred_weight")
head(ewaa_m3_sim81, n = 800)

# Convierte la columna de año a numérica si es necesario
ewaa_m3_sim81$year <- as.numeric(ewaa_m3_sim81$year)

# Ejecuta la función pad_weight_at_age e inspecciona la salida
padded_data81 <- pad_weight_at_age(ewaa_m3_sim81)
str(padded_data81)

# Escribe el archivo con el formato SS3
hakedataUSA:::write_wtatage_file(file = "ewaa_sim81.ss", data = padded_data81, maturity = maturity_at_age)

##################################################################
# Sim 82
set.seed(822)
devs_fyear82 <- rnorm(nT, mean = 0, sd = sd_fyear)
devs_fyear_fut82 <- rep(devs_fyear82, each = 16)
devs_fcohort82 <- rnorm(nT_23, mean = 0, sd = sd_fcohort)
devs_fcohort_past_future82 <- c(fcohort_estimate_values, devs_fcohort82)
estimate_fcohort_all82 <- data.frame(year = years_cohort, estimate = devs_fcohort_past_future82)
result_matrix82 <- matrix(0, nrow = 17, ncol = num_cols)
result_matrix82[1, ] <- format(as.integer(estimate_fcohort_all82$year), nsmall = 0)
result_matrix82[2, ] <- estimate_fcohort_all82$estimate
print(result_matrix82)
num_rows82 <- nrow(result_matrix82)
num_cols82 <- ncol(result_matrix82)
value_to_repeat82 <- result_matrix82[2, ]
start_row82 <- 3
start_col82 <- 2
end_row82 <- num_rows82
end_col82 <- num_cols82
for (i in start_row82:end_row82) {
  for (j in start_col82:end_col82) {
    if (j >= (i - start_row82 + start_col82)) {
      pattern_index82 <- (j - (i - start_row82 + start_col82)) %% 32 + 1
      result_matrix82[i, j] <- value_to_repeat82[pattern_index82]
    }
  }
}
print(result_matrix82)
extracted_values82 <- as.numeric(result_matrix82[2:17, 16:33])
extracted_values82 <- extracted_values82[-1]
devs_fcohort_all82 <- c(devs_fcohort_past, extracted_values82)
devs_fyear_all82 <- c(devs_fyear_past, devs_fyear_fut82)
ewaa_m3$fyear_devs82 <- devs_fyear_all82
ewaa_m3$fcohort_devs82 <- devs_fcohort_all82
ewaa_m3$ewaa_fyear82 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fyear_devs82)
ewaa_m3$ewaa_fcohort82 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fcohort_devs82)
ewaa_m3$ewaa_fyear_fcoh82 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fyear_devs82 + ewaa_m3$fcohort_devs82)
head(ewaa_m3, n = 10)

# Selecciona el peso predicho con el efecto de año y cohorte solamente.
ewaa_m3_sim82 <- ewaa_m3[, c("real_year", "age", "ewaa_fyear_fcoh82")]
head(ewaa_m3_sim82, n = 800)

# Renombra las columnas al formato SS3 ewaa.ss
ewaa_m3_sim82 <- ewaa_m3_sim82 %>%
  rename(year = real_year, pred_weight = ewaa_fyear_fcoh82)
colnames(ewaa_m3_sim82) <- c("year", "age", "pred_weight")
head(ewaa_m3_sim82, n = 800)

# Convierte la columna de año a numérica si es necesario
ewaa_m3_sim82$year <- as.numeric(ewaa_m3_sim82$year)

# Ejecuta la función pad_weight_at_age e inspecciona la salida
padded_data82 <- pad_weight_at_age(ewaa_m3_sim82)
str(padded_data82)

# Escribe el archivo con el formato SS3
hakedataUSA:::write_wtatage_file(file = "ewaa_sim82.ss", data = padded_data82, maturity = maturity_at_age)

##############################################################
# Sim 83
set.seed(833)
devs_fyear83 <- rnorm(nT, mean = 0, sd = sd_fyear)
devs_fyear_fut83 <- rep(devs_fyear83, each = 16)
devs_fcohort83 <- rnorm(nT_23, mean = 0, sd = sd_fcohort)
devs_fcohort_past_future83 <- c(fcohort_estimate_values, devs_fcohort83)
estimate_fcohort_all83 <- data.frame(year = years_cohort, estimate = devs_fcohort_past_future83)
result_matrix83 <- matrix(0, nrow = 17, ncol = num_cols)
result_matrix83[1, ] <- format(as.integer(estimate_fcohort_all83$year), nsmall = 0)
result_matrix83[2, ] <- estimate_fcohort_all83$estimate
print(result_matrix83)
num_rows83 <- nrow(result_matrix83)
num_cols83 <- ncol(result_matrix83)
value_to_repeat83 <- result_matrix83[2, ]
start_row83 <- 3
start_col83 <- 2
end_row83 <- num_rows83
end_col83 <- num_cols83
for (i in start_row83:end_row83) {
  for (j in start_col83:end_col83) {
    if (j >= (i - start_row83 + start_col83)) {
      pattern_index83 <- (j - (i - start_row83 + start_col83)) %% 32 + 1
      result_matrix83[i, j] <- value_to_repeat83[pattern_index83]
    }
  }
}
print(result_matrix83)
extracted_values83 <- as.numeric(result_matrix83[2:17, 16:33])
extracted_values83 <- extracted_values83[-1]
devs_fcohort_all83 <- c(devs_fcohort_past, extracted_values83)
devs_fyear_all83 <- c(devs_fyear_past, devs_fyear_fut83)
ewaa_m3$fyear_devs83 <- devs_fyear_all83
ewaa_m3$fcohort_devs83 <- devs_fcohort_all83
ewaa_m3$ewaa_fyear83 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fyear_devs83)
ewaa_m3$ewaa_fcohort83 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fcohort_devs83)
ewaa_m3$ewaa_fyear_fcoh83 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fyear_devs83 + ewaa_m3$fcohort_devs83)
head(ewaa_m3, n = 10)

# Selecciona el peso predicho con el efecto de año y cohorte solamente.
ewaa_m3_sim83 <- ewaa_m3[, c("real_year", "age", "ewaa_fyear_fcoh83")]
head(ewaa_m3_sim83, n = 800)

# Renombra las columnas al formato SS3 ewaa.ss
ewaa_m3_sim83 <- ewaa_m3_sim83 %>%
  rename(year = real_year, pred_weight = ewaa_fyear_fcoh83)
colnames(ewaa_m3_sim83) <- c("year", "age", "pred_weight")
head(ewaa_m3_sim83, n = 800)

# Convierte la columna de año a numérica si es necesario
ewaa_m3_sim83$year <- as.numeric(ewaa_m3_sim83$year)

# Ejecuta la función pad_weight_at_age e inspecciona la salida
padded_data83 <- pad_weight_at_age(ewaa_m3_sim83)
str(padded_data83)

# Escribe el archivo con el formato SS3
hakedataUSA:::write_wtatage_file(file = "ewaa_sim83.ss", data = padded_data83, maturity = maturity_at_age)

########################################################################
# Sim 84
set.seed(844)
devs_fyear84 <- rnorm(nT, mean = 0, sd = sd_fyear)
devs_fyear_fut84 <- rep(devs_fyear84, each = 16)
devs_fcohort84 <- rnorm(nT_23, mean = 0, sd = sd_fcohort)
devs_fcohort_past_future84 <- c(fcohort_estimate_values, devs_fcohort84)
estimate_fcohort_all84 <- data.frame(year = years_cohort, estimate = devs_fcohort_past_future84)
result_matrix84 <- matrix(0, nrow = 17, ncol = num_cols)
result_matrix84[1, ] <- format(as.integer(estimate_fcohort_all84$year), nsmall = 0)
result_matrix84[2, ] <- estimate_fcohort_all84$estimate
print(result_matrix84)
num_rows84 <- nrow(result_matrix84)
num_cols84 <- ncol(result_matrix84)
value_to_repeat84 <- result_matrix84[2, ]
start_row84 <- 3
start_col84 <- 2
end_row84 <- num_rows84
end_col84 <- num_cols84
for (i in start_row84:end_row84) {
  for (j in start_col84:end_col84) {
    if (j >= (i - start_row84 + start_col84)) {
      pattern_index84 <- (j - (i - start_row84 + start_col84)) %% 32 + 1
      result_matrix84[i, j] <- value_to_repeat84[pattern_index84]
    }
  }
}
print(result_matrix84)
extracted_values84 <- as.numeric(result_matrix84[2:17, 16:33])
extracted_values84 <- extracted_values84[-1]
devs_fcohort_all84 <- c(devs_fcohort_past, extracted_values84)
devs_fyear_all84 <- c(devs_fyear_past, devs_fyear_fut84)
ewaa_m3$fyear_devs84 <- devs_fyear_all84
ewaa_m3$fcohort_devs84 <- devs_fcohort_all84
ewaa_m3$ewaa_fyear84 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fyear_devs84)
ewaa_m3$ewaa_fcohort84 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fcohort_devs84)
ewaa_m3$ewaa_fyear_fcoh84 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fyear_devs84 + ewaa_m3$fcohort_devs84)
head(ewaa_m3, n = 10)

# Selecciona el peso predicho con el efecto de año y cohorte solamente.
ewaa_m3_sim84 <- ewaa_m3[, c("real_year", "age", "ewaa_fyear_fcoh84")]
head(ewaa_m3_sim84, n = 800)

# Renombra las columnas al formato SS3 ewaa.ss
ewaa_m3_sim84 <- ewaa_m3_sim84 %>%
  rename(year = real_year, pred_weight = ewaa_fyear_fcoh84)
colnames(ewaa_m3_sim84) <- c("year", "age", "pred_weight")
head(ewaa_m3_sim84, n = 800)

# Convierte la columna de año a numérica si es necesario
ewaa_m3_sim84$year <- as.numeric(ewaa_m3_sim84$year)

# Ejecuta la función pad_weight_at_age e inspecciona la salida
padded_data84 <- pad_weight_at_age(ewaa_m3_sim84)
str(padded_data84)

# Escribe el archivo con el formato SS3
hakedataUSA:::write_wtatage_file(file = "ewaa_sim84.ss", data = padded_data84, maturity = maturity_at_age)

###############################################################################
# Sim 85
set.seed(855)
devs_fyear85 <- rnorm(nT, mean = 0, sd = sd_fyear)
devs_fyear_fut85 <- rep(devs_fyear85, each = 16)
devs_fcohort85 <- rnorm(nT_23, mean = 0, sd = sd_fcohort)
devs_fcohort_past_future85 <- c(fcohort_estimate_values, devs_fcohort85)
estimate_fcohort_all85 <- data.frame(year = years_cohort, estimate = devs_fcohort_past_future85)
result_matrix85 <- matrix(0, nrow = 17, ncol = num_cols)
result_matrix85[1, ] <- format(as.integer(estimate_fcohort_all85$year), nsmall = 0)
result_matrix85[2, ] <- estimate_fcohort_all85$estimate
print(result_matrix85)
num_rows85 <- nrow(result_matrix85)
num_cols85 <- ncol(result_matrix85)
value_to_repeat85 <- result_matrix85[2, ]
start_row85 <- 3
start_col85 <- 2
end_row85 <- num_rows85
end_col85 <- num_cols85
for (i in start_row85:end_row85) {
  for (j in start_col85:end_col85) {
    if (j >= (i - start_row85 + start_col85)) {
      pattern_index85 <- (j - (i - start_row85 + start_col85)) %% 32 + 1
      result_matrix85[i, j] <- value_to_repeat85[pattern_index85]
    }
  }
}
print(result_matrix85)
extracted_values85 <- as.numeric(result_matrix85[2:17, 16:33])
extracted_values85 <- extracted_values85[-1]
devs_fcohort_all85 <- c(devs_fcohort_past, extracted_values85)
devs_fyear_all85 <- c(devs_fyear_past, devs_fyear_fut85)
ewaa_m3$fyear_devs85 <- devs_fyear_all85
ewaa_m3$fcohort_devs85 <- devs_fcohort_all85
ewaa_m3$ewaa_fyear85 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fyear_devs85)
ewaa_m3$ewaa_fcohort85 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fcohort_devs85)
ewaa_m3$ewaa_fyear_fcoh85 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fyear_devs85 + ewaa_m3$fcohort_devs85)
head(ewaa_m3, n = 10)

# Selecciona el peso predicho con el efecto de año y cohorte solamente.
ewaa_m3_sim85 <- ewaa_m3[, c("real_year", "age", "ewaa_fyear_fcoh85")]
head(ewaa_m3_sim85, n = 800)

# Renombra las columnas al formato SS3 ewaa.ss
ewaa_m3_sim85 <- ewaa_m3_sim85 %>%
  rename(year = real_year, pred_weight = ewaa_fyear_fcoh85)
colnames(ewaa_m3_sim85) <- c("year", "age", "pred_weight")
head(ewaa_m3_sim85, n = 800)

# Convierte la columna de año a numérica si es necesario
ewaa_m3_sim85$year <- as.numeric(ewaa_m3_sim85$year)

# Ejecuta la función pad_weight_at_age e inspecciona la salida
padded_data85 <- pad_weight_at_age(ewaa_m3_sim85)
str(padded_data85)

# Escribe el archivo con el formato SS3
hakedataUSA:::write_wtatage_file(file = "ewaa_sim85.ss", data = padded_data85, maturity = maturity_at_age)

#######################################################
# Sim 86
set.seed(866)
devs_fyear86 <- rnorm(nT, mean = 0, sd = sd_fyear)
devs_fyear_fut86 <- rep(devs_fyear86, each = 16)
devs_fcohort86 <- rnorm(nT_23, mean = 0, sd = sd_fcohort)
devs_fcohort_past_future86 <- c(fcohort_estimate_values, devs_fcohort86)
estimate_fcohort_all86 <- data.frame(year = years_cohort, estimate = devs_fcohort_past_future86)
result_matrix86 <- matrix(0, nrow = 17, ncol = num_cols)
result_matrix86[1, ] <- format(as.integer(estimate_fcohort_all86$year), nsmall = 0)
result_matrix86[2, ] <- estimate_fcohort_all86$estimate
print(result_matrix86)
num_rows86 <- nrow(result_matrix86)
num_cols86 <- ncol(result_matrix86)
value_to_repeat86 <- result_matrix86[2, ]
start_row86 <- 3
start_col86 <- 2
end_row86 <- num_rows86
end_col86 <- num_cols86
for (i in start_row86:end_row86) {
  for (j in start_col86:end_col86) {
    if (j >= (i - start_row86 + start_col86)) {
      pattern_index86 <- (j - (i - start_row86 + start_col86)) %% 32 + 1
      result_matrix86[i, j] <- value_to_repeat86[pattern_index86]
    }
  }
}
print(result_matrix86)
extracted_values86 <- as.numeric(result_matrix86[2:17, 16:33])
extracted_values86 <- extracted_values86[-1]
devs_fcohort_all86 <- c(devs_fcohort_past, extracted_values86)
devs_fyear_all86 <- c(devs_fyear_past, devs_fyear_fut86)
ewaa_m3$fyear_devs86 <- devs_fyear_all86
ewaa_m3$fcohort_devs86 <- devs_fcohort_all86
ewaa_m3$ewaa_fyear86 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fyear_devs86)
ewaa_m3$ewaa_fcohort86 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fcohort_devs86)
ewaa_m3$ewaa_fyear_fcoh86 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fyear_devs86 + ewaa_m3$fcohort_devs86)
head(ewaa_m3, n = 10)

# Selecciona el peso predicho con el efecto de año y cohorte solamente.
ewaa_m3_sim86 <- ewaa_m3[, c("real_year", "age", "ewaa_fyear_fcoh86")]
head(ewaa_m3_sim86, n = 800)

# Renombra las columnas al formato SS3 ewaa.ss
ewaa_m3_sim86 <- ewaa_m3_sim86 %>%
  rename(year = real_year, pred_weight = ewaa_fyear_fcoh86)
colnames(ewaa_m3_sim86) <- c("year", "age", "pred_weight")
head(ewaa_m3_sim86, n = 800)

# Convierte la columna de año a numérica si es necesario
ewaa_m3_sim86$year <- as.numeric(ewaa_m3_sim86$year)

# Ejecuta la función pad_weight_at_age e inspecciona la salida
padded_data86 <- pad_weight_at_age(ewaa_m3_sim86)
str(padded_data86)

# Escribe el archivo con el formato SS3
hakedataUSA:::write_wtatage_file(file = "ewaa_sim86.ss", data = padded_data86, maturity = maturity_at_age)
#######################################################################################
# Sim 87
set.seed(877)
devs_fyear87 <- rnorm(nT, mean = 0, sd = sd_fyear)
devs_fyear_fut87 <- rep(devs_fyear87, each = 16)
devs_fcohort87 <- rnorm(nT_23, mean = 0, sd = sd_fcohort)
devs_fcohort_past_future87 <- c(fcohort_estimate_values, devs_fcohort87)
estimate_fcohort_all87 <- data.frame(year = years_cohort, estimate = devs_fcohort_past_future87)
result_matrix87 <- matrix(0, nrow = 17, ncol = num_cols)
result_matrix87[1, ] <- format(as.integer(estimate_fcohort_all87$year), nsmall = 0)
result_matrix87[2, ] <- estimate_fcohort_all87$estimate
print(result_matrix87)
num_rows87 <- nrow(result_matrix87)
num_cols87 <- ncol(result_matrix87)
value_to_repeat87 <- result_matrix87[2, ]
start_row87 <- 3
start_col87 <- 2
end_row87 <- num_rows87
end_col87 <- num_cols87
for (i in start_row87:end_row87) {
  for (j in start_col87:end_col87) {
    if (j >= (i - start_row87 + start_col87)) {
      pattern_index87 <- (j - (i - start_row87 + start_col87)) %% 32 + 1
      result_matrix87[i, j] <- value_to_repeat87[pattern_index87]
    }
  }
}
print(result_matrix87)
extracted_values87 <- as.numeric(result_matrix87[2:17, 16:33])
extracted_values87 <- extracted_values87[-1]
devs_fcohort_all87 <- c(devs_fcohort_past, extracted_values87)
devs_fyear_all87 <- c(devs_fyear_past, devs_fyear_fut87)
ewaa_m3$fyear_devs87 <- devs_fyear_all87
ewaa_m3$fcohort_devs87 <- devs_fcohort_all87
ewaa_m3$ewaa_fyear87 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fyear_devs87)
ewaa_m3$ewaa_fcohort87 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fcohort_devs87)
ewaa_m3$ewaa_fyear_fcoh87 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fyear_devs87 + ewaa_m3$fcohort_devs87)
head(ewaa_m3, n = 10)

# Selecciona el peso predicho con el efecto de año y cohorte solamente.
ewaa_m3_sim87 <- ewaa_m3[, c("real_year", "age", "ewaa_fyear_fcoh87")]
head(ewaa_m3_sim87, n = 800)

# Renombra las columnas al formato SS3 ewaa.ss
ewaa_m3_sim87 <- ewaa_m3_sim87 %>%
  rename(year = real_year, pred_weight = ewaa_fyear_fcoh87)
colnames(ewaa_m3_sim87) <- c("year", "age", "pred_weight")
head(ewaa_m3_sim87, n = 800)

# Convierte la columna de año a numérica si es necesario
ewaa_m3_sim87$year <- as.numeric(ewaa_m3_sim87$year)

# Ejecuta la función pad_weight_at_age e inspecciona la salida
padded_data87 <- pad_weight_at_age(ewaa_m3_sim87)
str(padded_data87)

# Escribe el archivo con el formato SS3
hakedataUSA:::write_wtatage_file(file = "ewaa_sim87.ss", data = padded_data87, maturity = maturity_at_age)

#######################################################
# Sim 88
set.seed(888)
devs_fyear88 <- rnorm(nT, mean = 0, sd = sd_fyear)
devs_fyear_fut88 <- rep(devs_fyear88, each = 16)
devs_fcohort88 <- rnorm(nT_23, mean = 0, sd = sd_fcohort)
devs_fcohort_past_future88 <- c(fcohort_estimate_values, devs_fcohort88)
estimate_fcohort_all88 <- data.frame(year = years_cohort, estimate = devs_fcohort_past_future88)
result_matrix88 <- matrix(0, nrow = 17, ncol = num_cols)
result_matrix88[1, ] <- format(as.integer(estimate_fcohort_all88$year), nsmall = 0)
result_matrix88[2, ] <- estimate_fcohort_all88$estimate
print(result_matrix88)
num_rows88 <- nrow(result_matrix88)
num_cols88 <- ncol(result_matrix88)
value_to_repeat88 <- result_matrix88[2, ]
start_row88 <- 3
start_col88 <- 2
end_row88 <- num_rows88
end_col88 <- num_cols88
for (i in start_row88:end_row88) {
  for (j in start_col88:end_col88) {
    if (j >= (i - start_row88 + start_col88)) {
      pattern_index88 <- (j - (i - start_row88 + start_col88)) %% 32 + 1
      result_matrix88[i, j] <- value_to_repeat88[pattern_index88]
    }
  }
}
print(result_matrix88)
extracted_values88 <- as.numeric(result_matrix88[2:17, 16:33])
extracted_values88 <- extracted_values88[-1]
devs_fcohort_all88 <- c(devs_fcohort_past, extracted_values88)
devs_fyear_all88 <- c(devs_fyear_past, devs_fyear_fut88)
ewaa_m3$fyear_devs88 <- devs_fyear_all88
ewaa_m3$fcohort_devs88 <- devs_fcohort_all88
ewaa_m3$ewaa_fyear88 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fyear_devs88)
ewaa_m3$ewaa_fcohort88 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fcohort_devs88)
ewaa_m3$ewaa_fyear_fcoh88 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fyear_devs88 + ewaa_m3$fcohort_devs88)
head(ewaa_m3, n = 10)

# Selecciona el peso predicho con el efecto de año y cohorte solamente.
ewaa_m3_sim88 <- ewaa_m3[, c("real_year", "age", "ewaa_fyear_fcoh88")]
head(ewaa_m3_sim88, n = 800)

# Renombra las columnas al formato SS3 ewaa.ss
ewaa_m3_sim88 <- ewaa_m3_sim88 %>%
  rename(year = real_year, pred_weight = ewaa_fyear_fcoh88)
colnames(ewaa_m3_sim88) <- c("year", "age", "pred_weight")
head(ewaa_m3_sim88, n = 800)

# Convierte la columna de año a numérica si es necesario
ewaa_m3_sim88$year <- as.numeric(ewaa_m3_sim88$year)

# Ejecuta la función pad_weight_at_age e inspecciona la salida
padded_data88 <- pad_weight_at_age(ewaa_m3_sim88)
str(padded_data88)

# Escribe el archivo con el formato SS3
hakedataUSA:::write_wtatage_file(file = "ewaa_sim88.ss", data = padded_data88, maturity = maturity_at_age)

#################################################################
# Sim 89
set.seed(899)
devs_fyear89 <- rnorm(nT, mean = 0, sd = sd_fyear)
devs_fyear_fut89 <- rep(devs_fyear89, each = 16)
devs_fcohort89 <- rnorm(nT_23, mean = 0, sd = sd_fcohort)
devs_fcohort_past_future89 <- c(fcohort_estimate_values, devs_fcohort89)
estimate_fcohort_all89 <- data.frame(year = years_cohort, estimate = devs_fcohort_past_future89)
result_matrix89 <- matrix(0, nrow = 17, ncol = num_cols)
result_matrix89[1, ] <- format(as.integer(estimate_fcohort_all89$year), nsmall = 0)
result_matrix89[2, ] <- estimate_fcohort_all89$estimate
print(result_matrix89)
num_rows89 <- nrow(result_matrix89)
num_cols89 <- ncol(result_matrix89)
value_to_repeat89 <- result_matrix89[2, ]
start_row89 <- 3
start_col89 <- 2
end_row89 <- num_rows89
end_col89 <- num_cols89
for (i in start_row89:end_row89) {
  for (j in start_col89:end_col89) {
    if (j >= (i - start_row89 + start_col89)) {
      pattern_index89 <- (j - (i - start_row89 + start_col89)) %% 32 + 1
      result_matrix89[i, j] <- value_to_repeat89[pattern_index89]
    }
  }
}
print(result_matrix89)
extracted_values89 <- as.numeric(result_matrix89[2:17, 16:33])
extracted_values89 <- extracted_values89[-1]
devs_fcohort_all89 <- c(devs_fcohort_past, extracted_values89)
devs_fyear_all89 <- c(devs_fyear_past, devs_fyear_fut89)
ewaa_m3$fyear_devs89 <- devs_fyear_all89
ewaa_m3$fcohort_devs89 <- devs_fcohort_all89
ewaa_m3$ewaa_fyear89 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fyear_devs89)
ewaa_m3$ewaa_fcohort89 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fcohort_devs89)
ewaa_m3$ewaa_fyear_fcoh89 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fyear_devs89 + ewaa_m3$fcohort_devs89)
head(ewaa_m3, n = 10)

# Selecciona el peso predicho con el efecto de año y cohorte solamente.
ewaa_m3_sim89 <- ewaa_m3[, c("real_year", "age", "ewaa_fyear_fcoh89")]
head(ewaa_m3_sim89, n = 800)

# Renombra las columnas al formato SS3 ewaa.ss
ewaa_m3_sim89 <- ewaa_m3_sim89 %>%
  rename(year = real_year, pred_weight = ewaa_fyear_fcoh89)
colnames(ewaa_m3_sim89) <- c("year", "age", "pred_weight")
head(ewaa_m3_sim89, n = 800)

# Convierte la columna de año a numérica si es necesario
ewaa_m3_sim89$year <- as.numeric(ewaa_m3_sim89$year)

# Ejecuta la función pad_weight_at_age e inspecciona la salida
padded_data89 <- pad_weight_at_age(ewaa_m3_sim89)
str(padded_data89)

# Escribe el archivo con el formato SS3
hakedataUSA:::write_wtatage_file(file = "ewaa_sim89.ss", data = padded_data89, maturity = maturity_at_age)
###################################################################
# Sim 90
set.seed(900)
devs_fyear90 <- rnorm(nT, mean = 0, sd = sd_fyear)
devs_fyear_fut90 <- rep(devs_fyear90, each = 16)
devs_fcohort90 <- rnorm(nT_23, mean = 0, sd = sd_fcohort)
devs_fcohort_past_future90 <- c(fcohort_estimate_values, devs_fcohort90)
estimate_fcohort_all90 <- data.frame(year = years_cohort, estimate = devs_fcohort_past_future90)
result_matrix90 <- matrix(0, nrow = 17, ncol = num_cols)
result_matrix90[1, ] <- format(as.integer(estimate_fcohort_all90$year), nsmall = 0)
result_matrix90[2, ] <- estimate_fcohort_all90$estimate
print(result_matrix90)
num_rows90 <- nrow(result_matrix90)
num_cols90 <- ncol(result_matrix90)
value_to_repeat90 <- result_matrix90[2, ]
start_row90 <- 3
start_col90 <- 2
end_row90 <- num_rows90
end_col90 <- num_cols90
for (i in start_row90:end_row90) {
  for (j in start_col90:end_col90) {
    if (j >= (i - start_row90 + start_col90)) {
      pattern_index90 <- (j - (i - start_row90 + start_col90)) %% 32 + 1
      result_matrix90[i, j] <- value_to_repeat90[pattern_index90]
    }
  }
}
print(result_matrix90)
extracted_values90 <- as.numeric(result_matrix90[2:17, 16:33])
extracted_values90 <- extracted_values90[-1]
devs_fcohort_all90 <- c(devs_fcohort_past, extracted_values90)
devs_fyear_all90 <- c(devs_fyear_past, devs_fyear_fut90)
ewaa_m3$fyear_devs90 <- devs_fyear_all90
ewaa_m3$fcohort_devs90 <- devs_fcohort_all90
ewaa_m3$ewaa_fyear90 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fyear_devs90)
ewaa_m3$ewaa_fcohort90 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fcohort_devs90)
ewaa_m3$ewaa_fyear_fcoh90 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fyear_devs90 + ewaa_m3$fcohort_devs90)
head(ewaa_m3, n = 10)

# Selecciona el peso predicho con el efecto de año y cohorte solamente.
ewaa_m3_sim90 <- ewaa_m3[, c("real_year", "age", "ewaa_fyear_fcoh90")]
head(ewaa_m3_sim90, n = 800)

# Renombra las columnas al formato SS3 ewaa.ss
ewaa_m3_sim90 <- ewaa_m3_sim90 %>%
  rename(year = real_year, pred_weight = ewaa_fyear_fcoh90)
colnames(ewaa_m3_sim90) <- c("year", "age", "pred_weight")
head(ewaa_m3_sim90, n = 800)

# Convierte la columna de año a numérica si es necesario
ewaa_m3_sim90$year <- as.numeric(ewaa_m3_sim90$year)

# Ejecuta la función pad_weight_at_age e inspecciona la salida
padded_data90 <- pad_weight_at_age(ewaa_m3_sim90)
str(padded_data90)

# Escribe el archivo con el formato SS3
hakedataUSA:::write_wtatage_file(file = "ewaa_sim90.ss", data = padded_data90, maturity = maturity_at_age)

###############################################################
# Sim 91
set.seed(911)
devs_fyear91 <- rnorm(nT, mean = 0, sd = sd_fyear)
devs_fyear_fut91 <- rep(devs_fyear91, each = 16)
devs_fcohort91 <- rnorm(nT_23, mean = 0, sd = sd_fcohort)
devs_fcohort_past_future91 <- c(fcohort_estimate_values, devs_fcohort91)
estimate_fcohort_all91 <- data.frame(year = years_cohort, estimate = devs_fcohort_past_future91)
result_matrix91 <- matrix(0, nrow = 17, ncol = num_cols)
result_matrix91[1, ] <- format(as.integer(estimate_fcohort_all91$year), nsmall = 0)
result_matrix91[2, ] <- estimate_fcohort_all91$estimate
print(result_matrix91)
num_rows91 <- nrow(result_matrix91)
num_cols91 <- ncol(result_matrix91)
value_to_repeat91 <- result_matrix91[2, ]
start_row91 <- 3
start_col91 <- 2
end_row91 <- num_rows91
end_col91 <- num_cols91
for (i in start_row91:end_row91) {
  for (j in start_col91:end_col91) {
    if (j >= (i - start_row91 + start_col91)) {
      pattern_index91 <- (j - (i - start_row91 + start_col91)) %% 32 + 1
      result_matrix91[i, j] <- value_to_repeat91[pattern_index91]
    }
  }
}
print(result_matrix91)
extracted_values91 <- as.numeric(result_matrix91[2:17, 16:33])
extracted_values91 <- extracted_values91[-1]
devs_fcohort_all91 <- c(devs_fcohort_past, extracted_values91)
devs_fyear_all91 <- c(devs_fyear_past, devs_fyear_fut91)
ewaa_m3$fyear_devs91 <- devs_fyear_all91
ewaa_m3$fcohort_devs91 <- devs_fcohort_all91
ewaa_m3$ewaa_fyear91 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fyear_devs91)
ewaa_m3$ewaa_fcohort91 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fcohort_devs91)
ewaa_m3$ewaa_fyear_fcoh91 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fyear_devs91 + ewaa_m3$fcohort_devs91)
head(ewaa_m3, n = 10)

# Selecciona el peso predicho con el efecto de año y cohorte solamente.
ewaa_m3_sim91 <- ewaa_m3[, c("real_year", "age", "ewaa_fyear_fcoh91")]
head(ewaa_m3_sim91, n = 800)

# Renombra las columnas al formato SS3 ewaa.ss
ewaa_m3_sim91 <- ewaa_m3_sim91 %>%
  rename(year = real_year, pred_weight = ewaa_fyear_fcoh91)
colnames(ewaa_m3_sim91) <- c("year", "age", "pred_weight")
head(ewaa_m3_sim91, n = 800)

# Convierte la columna de año a numérica si es necesario
ewaa_m3_sim91$year <- as.numeric(ewaa_m3_sim91$year)

# Ejecuta la función pad_weight_at_age e inspecciona la salida
padded_data91 <- pad_weight_at_age(ewaa_m3_sim91)
str(padded_data91)

# Escribe el archivo con el formato SS3
hakedataUSA:::write_wtatage_file(file = "ewaa_sim91.ss", data = padded_data91, maturity = maturity_at_age)
#######################################################

# Sim 92
set.seed(920)
devs_fyear92 <- rnorm(nT, mean = 0, sd = sd_fyear)
devs_fyear_fut92 <- rep(devs_fyear92, each = 16)
devs_fcohort92 <- rnorm(nT_23, mean = 0, sd = sd_fcohort)
devs_fcohort_past_future92 <- c(fcohort_estimate_values, devs_fcohort92)
estimate_fcohort_all92 <- data.frame(year = years_cohort, estimate = devs_fcohort_past_future92)
result_matrix92 <- matrix(0, nrow = 17, ncol = num_cols)
result_matrix92[1, ] <- format(as.integer(estimate_fcohort_all92$year), nsmall = 0)
result_matrix92[2, ] <- estimate_fcohort_all92$estimate
print(result_matrix92)
num_rows92 <- nrow(result_matrix92)
num_cols92 <- ncol(result_matrix92)
value_to_repeat92 <- result_matrix92[2, ]
start_row92 <- 3
start_col92 <- 2
end_row92 <- num_rows92
end_col92 <- num_cols92
for (i in start_row92:end_row92) {
  for (j in start_col92:end_col92) {
    if (j >= (i - start_row92 + start_col92)) {
      pattern_index92 <- (j - (i - start_row92 + start_col92)) %% 32 + 1
      result_matrix92[i, j] <- value_to_repeat92[pattern_index92]
    }
  }
}
print(result_matrix92)
extracted_values92 <- as.numeric(result_matrix92[2:17, 16:33])
extracted_values92 <- extracted_values92[-1]
devs_fcohort_all92 <- c(devs_fcohort_past, extracted_values92)
devs_fyear_all92 <- c(devs_fyear_past, devs_fyear_fut92)
ewaa_m3$fyear_devs92 <- devs_fyear_all92
ewaa_m3$fcohort_devs92 <- devs_fcohort_all92
ewaa_m3$ewaa_fyear92 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fyear_devs92)
ewaa_m3$ewaa_fcohort92 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fcohort_devs92)
ewaa_m3$ewaa_fyear_fcoh92 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fyear_devs92 + ewaa_m3$fcohort_devs92)
head(ewaa_m3, n = 10)

# Selecciona el peso predicho con el efecto de año y cohorte solamente.
ewaa_m3_sim92 <- ewaa_m3[, c("real_year", "age", "ewaa_fyear_fcoh92")]
head(ewaa_m3_sim92, n = 800)

# Renombra las columnas al formato SS3 ewaa.ss
ewaa_m3_sim92 <- ewaa_m3_sim92 %>%
  rename(year = real_year, pred_weight = ewaa_fyear_fcoh92)
colnames(ewaa_m3_sim92) <- c("year", "age", "pred_weight")
head(ewaa_m3_sim92, n = 800)

# Convierte la columna de año a numérica si es necesario
ewaa_m3_sim92$year <- as.numeric(ewaa_m3_sim92$year)

# Ejecuta la función pad_weight_at_age e inspecciona la salida
padded_data92 <- pad_weight_at_age(ewaa_m3_sim92)
str(padded_data92)

# Escribe el archivo con el formato SS3
hakedataUSA:::write_wtatage_file(file = "ewaa_sim92.ss", data = padded_data92, maturity = maturity_at_age)

########################################
# Sim 93
set.seed(930)
devs_fyear93 <- rnorm(nT, mean = 0, sd = sd_fyear)
devs_fyear_fut93 <- rep(devs_fyear93, each = 16)
devs_fcohort93 <- rnorm(nT_23, mean = 0, sd = sd_fcohort)
devs_fcohort_past_future93 <- c(fcohort_estimate_values, devs_fcohort93)
estimate_fcohort_all93 <- data.frame(year = years_cohort, estimate = devs_fcohort_past_future93)
result_matrix93 <- matrix(0, nrow = 17, ncol = num_cols)
result_matrix93[1, ] <- format(as.integer(estimate_fcohort_all93$year), nsmall = 0)
result_matrix93[2, ] <- estimate_fcohort_all93$estimate
print(result_matrix93)
num_rows93 <- nrow(result_matrix93)
num_cols93 <- ncol(result_matrix93)
value_to_repeat93 <- result_matrix93[2, ]
start_row93 <- 3
start_col93 <- 2
end_row93 <- num_rows93
end_col93 <- num_cols93
for (i in start_row93:end_row93) {
  for (j in start_col93:end_col93) {
    if (j >= (i - start_row93 + start_col93)) {
      pattern_index93 <- (j - (i - start_row93 + start_col93)) %% 32 + 1
      result_matrix93[i, j] <- value_to_repeat93[pattern_index93]
    }
  }
}
print(result_matrix93)
extracted_values93 <- as.numeric(result_matrix93[2:17, 16:33])
extracted_values93 <- extracted_values93[-1]
devs_fcohort_all93 <- c(devs_fcohort_past, extracted_values93)
devs_fyear_all93 <- c(devs_fyear_past, devs_fyear_fut93)
ewaa_m3$fyear_devs93 <- devs_fyear_all93
ewaa_m3$fcohort_devs93 <- devs_fcohort_all93
ewaa_m3$ewaa_fyear93 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fyear_devs93)
ewaa_m3$ewaa_fcohort93 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fcohort_devs93)
ewaa_m3$ewaa_fyear_fcoh93 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fyear_devs93 + ewaa_m3$fcohort_devs93)
head(ewaa_m3, n = 10)

# Selecciona el peso predicho con el efecto de año y cohorte solamente.
ewaa_m3_sim93 <- ewaa_m3[, c("real_year", "age", "ewaa_fyear_fcoh93")]
head(ewaa_m3_sim93, n = 800)

# Renombra las columnas al formato SS3 ewaa.ss
ewaa_m3_sim93 <- ewaa_m3_sim93 %>%
  rename(year = real_year, pred_weight = ewaa_fyear_fcoh93)
colnames(ewaa_m3_sim93) <- c("year", "age", "pred_weight")
head(ewaa_m3_sim93, n = 800)

# Convierte la columna de año a numérica si es necesario
ewaa_m3_sim93$year <- as.numeric(ewaa_m3_sim93$year)

# Ejecuta la función pad_weight_at_age e inspecciona la salida
padded_data93 <- pad_weight_at_age(ewaa_m3_sim93)
str(padded_data93)

# Escribe el archivo con el formato SS3
hakedataUSA:::write_wtatage_file(file = "ewaa_sim93.ss", data = padded_data93, maturity = maturity_at_age)

############################################################################
# Sim 94
set.seed(94)
devs_fyear94 <- rnorm(nT, mean = 0, sd = sd_fyear)
devs_fyear_fut94 <- rep(devs_fyear94, each = 16)
devs_fcohort94 <- rnorm(nT_23, mean = 0, sd = sd_fcohort)
devs_fcohort_past_future94 <- c(fcohort_estimate_values, devs_fcohort94)
estimate_fcohort_all94 <- data.frame(year = years_cohort, estimate = devs_fcohort_past_future94)
result_matrix94 <- matrix(0, nrow = 17, ncol = num_cols)
result_matrix94[1, ] <- format(as.integer(estimate_fcohort_all94$year), nsmall = 0)
result_matrix94[2, ] <- estimate_fcohort_all94$estimate
print(result_matrix94)
num_rows94 <- nrow(result_matrix94)
num_cols94 <- ncol(result_matrix94)
value_to_repeat94 <- result_matrix94[2, ]
start_row94 <- 3
start_col94 <- 2
end_row94 <- num_rows94
end_col94 <- num_cols94
for (i in start_row94:end_row94) {
  for (j in start_col94:end_col94) {
    if (j >= (i - start_row94 + start_col94)) {
      pattern_index94 <- (j - (i - start_row94 + start_col94)) %% 32 + 1
      result_matrix94[i, j] <- value_to_repeat94[pattern_index94]
    }
  }
}
print(result_matrix94)
extracted_values94 <- as.numeric(result_matrix94[2:17, 16:33])
extracted_values94 <- extracted_values94[-1]
devs_fcohort_all94 <- c(devs_fcohort_past, extracted_values94)
devs_fyear_all94 <- c(devs_fyear_past, devs_fyear_fut94)
ewaa_m3$fyear_devs94 <- devs_fyear_all94
ewaa_m3$fcohort_devs94 <- devs_fcohort_all94
ewaa_m3$ewaa_fyear94 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fyear_devs94)
ewaa_m3$ewaa_fcohort94 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fcohort_devs94)
ewaa_m3$ewaa_fyear_fcoh94 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fyear_devs94 + ewaa_m3$fcohort_devs94)
head(ewaa_m3, n = 10)

# Selecciona el peso predicho con el efecto de año y cohorte solamente.
ewaa_m3_sim94 <- ewaa_m3[, c("real_year", "age", "ewaa_fyear_fcoh94")]
head(ewaa_m3_sim94, n = 800)

# Renombra las columnas al formato SS3 ewaa.ss
ewaa_m3_sim94 <- ewaa_m3_sim94 %>%
  rename(year = real_year, pred_weight = ewaa_fyear_fcoh94)
colnames(ewaa_m3_sim94) <- c("year", "age", "pred_weight")
head(ewaa_m3_sim94, n = 800)

# Convierte la columna de año a numérica si es necesario
ewaa_m3_sim94$year <- as.numeric(ewaa_m3_sim94$year)

# Ejecuta la función pad_weight_at_age e inspecciona la salida
padded_data94 <- pad_weight_at_age(ewaa_m3_sim94)
str(padded_data94)

# Escribe el archivo con el formato SS3
hakedataUSA:::write_wtatage_file(file = "ewaa_sim94.ss", data = padded_data94, maturity = maturity_at_age)
###########################################################

# Sim 95
set.seed(950)
devs_fyear95 <- rnorm(nT, mean = 0, sd = sd_fyear)
devs_fyear_fut95 <- rep(devs_fyear95, each = 16)
devs_fcohort95 <- rnorm(nT_23, mean = 0, sd = sd_fcohort)
devs_fcohort_past_future95 <- c(fcohort_estimate_values, devs_fcohort95)
estimate_fcohort_all95 <- data.frame(year = years_cohort, estimate = devs_fcohort_past_future95)
result_matrix95 <- matrix(0, nrow = 17, ncol = num_cols)
result_matrix95[1, ] <- format(as.integer(estimate_fcohort_all95$year), nsmall = 0)
result_matrix95[2, ] <- estimate_fcohort_all95$estimate
print(result_matrix95)
num_rows95 <- nrow(result_matrix95)
num_cols95 <- ncol(result_matrix95)
value_to_repeat95 <- result_matrix95[2, ]
start_row95 <- 3
start_col95 <- 2
end_row95 <- num_rows95
end_col95 <- num_cols95
for (i in start_row95:end_row95) {
  for (j in start_col95:end_col95) {
    if (j >= (i - start_row95 + start_col95)) {
      pattern_index95 <- (j - (i - start_row95 + start_col95)) %% 32 + 1
      result_matrix95[i, j] <- value_to_repeat95[pattern_index95]
    }
  }
}
print(result_matrix95)
extracted_values95 <- as.numeric(result_matrix95[2:17, 16:33])
extracted_values95 <- extracted_values95[-1]
devs_fcohort_all95 <- c(devs_fcohort_past, extracted_values95)
devs_fyear_all95 <- c(devs_fyear_past, devs_fyear_fut95)
ewaa_m3$fyear_devs95 <- devs_fyear_all95
ewaa_m3$fcohort_devs95 <- devs_fcohort_all95
ewaa_m3$ewaa_fyear95 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fyear_devs95)
ewaa_m3$ewaa_fcohort95 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fcohort_devs95)
ewaa_m3$ewaa_fyear_fcoh95 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fyear_devs95 + ewaa_m3$fcohort_devs95)
head(ewaa_m3, n = 10)

# Selecciona el peso predicho con el efecto de año y cohorte solamente.
ewaa_m3_sim95 <- ewaa_m3[, c("real_year", "age", "ewaa_fyear_fcoh95")]
head(ewaa_m3_sim95, n = 800)

# Renombra las columnas al formato SS3 ewaa.ss
ewaa_m3_sim95 <- ewaa_m3_sim95 %>%
  rename(year = real_year, pred_weight = ewaa_fyear_fcoh95)
colnames(ewaa_m3_sim95) <- c("year", "age", "pred_weight")
head(ewaa_m3_sim95, n = 800)

# Convierte la columna de año a numérica si es necesario
ewaa_m3_sim95$year <- as.numeric(ewaa_m3_sim95$year)

# Ejecuta la función pad_weight_at_age e inspecciona la salida
padded_data95 <- pad_weight_at_age(ewaa_m3_sim95)
str(padded_data95)

# Escribe el archivo con el formato SS3
hakedataUSA:::write_wtatage_file(file = "ewaa_sim95.ss", data = padded_data95, maturity = maturity_at_age)

#########################################################################
# Sim 96
set.seed(960)
devs_fyear96 <- rnorm(nT, mean = 0, sd = sd_fyear)
devs_fyear_fut96 <- rep(devs_fyear96, each = 16)
devs_fcohort96 <- rnorm(nT_23, mean = 0, sd = sd_fcohort)
devs_fcohort_past_future96 <- c(fcohort_estimate_values, devs_fcohort96)
estimate_fcohort_all96 <- data.frame(year = years_cohort, estimate = devs_fcohort_past_future96)
result_matrix96 <- matrix(0, nrow = 17, ncol = num_cols)
result_matrix96[1, ] <- format(as.integer(estimate_fcohort_all96$year), nsmall = 0)
result_matrix96[2, ] <- estimate_fcohort_all96$estimate
print(result_matrix96)
num_rows96 <- nrow(result_matrix96)
num_cols96 <- ncol(result_matrix96)
value_to_repeat96 <- result_matrix96[2, ]
start_row96 <- 3
start_col96 <- 2
end_row96 <- num_rows96
end_col96 <- num_cols96
for (i in start_row96:end_row96) {
  for (j in start_col96:end_col96) {
    if (j >= (i - start_row96 + start_col96)) {
      pattern_index96 <- (j - (i - start_row96 + start_col96)) %% 32 + 1
      result_matrix96[i, j] <- value_to_repeat96[pattern_index96]
    }
  }
}
print(result_matrix96)
extracted_values96 <- as.numeric(result_matrix96[2:17, 16:33])
extracted_values96 <- extracted_values96[-1]
devs_fcohort_all96 <- c(devs_fcohort_past, extracted_values96)
devs_fyear_all96 <- c(devs_fyear_past, devs_fyear_fut96)
ewaa_m3$fyear_devs96 <- devs_fyear_all96
ewaa_m3$fcohort_devs96 <- devs_fcohort_all96
ewaa_m3$ewaa_fyear96 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fyear_devs96)
ewaa_m3$ewaa_fcohort96 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fcohort_devs96)
ewaa_m3$ewaa_fyear_fcoh96 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fyear_devs96 + ewaa_m3$fcohort_devs96)
head(ewaa_m3, n = 10)

# Selecciona el peso predicho con el efecto de año y cohorte solamente.
ewaa_m3_sim96 <- ewaa_m3[, c("real_year", "age", "ewaa_fyear_fcoh96")]
head(ewaa_m3_sim96, n = 800)

# Renombra las columnas al formato SS3 ewaa.ss
ewaa_m3_sim96 <- ewaa_m3_sim96 %>%
  rename(year = real_year, pred_weight = ewaa_fyear_fcoh96)
colnames(ewaa_m3_sim96) <- c("year", "age", "pred_weight")
head(ewaa_m3_sim96, n = 800)

# Convierte la columna de año a numérica si es necesario
ewaa_m3_sim96$year <- as.numeric(ewaa_m3_sim96$year)

# Ejecuta la función pad_weight_at_age e inspecciona la salida
padded_data96 <- pad_weight_at_age(ewaa_m3_sim96)
str(padded_data96)

# Escribe el archivo con el formato SS3
hakedataUSA:::write_wtatage_file(file = "ewaa_sim96.ss", data = padded_data96, maturity = maturity_at_age)

########################################################################
# Sim 97
set.seed(97)
devs_fyear97 <- rnorm(nT, mean = 0, sd = sd_fyear)
devs_fyear_fut97 <- rep(devs_fyear97, each = 16)
devs_fcohort97 <- rnorm(nT_23, mean = 0, sd = sd_fcohort)
devs_fcohort_past_future97 <- c(fcohort_estimate_values, devs_fcohort97)
estimate_fcohort_all97 <- data.frame(year = years_cohort, estimate = devs_fcohort_past_future97)
result_matrix97 <- matrix(0, nrow = 17, ncol = num_cols)
result_matrix97[1, ] <- format(as.integer(estimate_fcohort_all97$year), nsmall = 0)
result_matrix97[2, ] <- estimate_fcohort_all97$estimate
print(result_matrix97)
num_rows97 <- nrow(result_matrix97)
num_cols97 <- ncol(result_matrix97)
value_to_repeat97 <- result_matrix97[2, ]
start_row97 <- 3
start_col97 <- 2
end_row97 <- num_rows97
end_col97 <- num_cols97
for (i in start_row97:end_row97) {
  for (j in start_col97:end_col97) {
    if (j >= (i - start_row97 + start_col97)) {
      pattern_index97 <- (j - (i - start_row97 + start_col97)) %% 32 + 1
      result_matrix97[i, j] <- value_to_repeat97[pattern_index97]
    }
  }
}
print(result_matrix97)
extracted_values97 <- as.numeric(result_matrix97[2:17, 16:33])
extracted_values97 <- extracted_values97[-1]
devs_fcohort_all97 <- c(devs_fcohort_past, extracted_values97)
devs_fyear_all97 <- c(devs_fyear_past, devs_fyear_fut97)
ewaa_m3$fyear_devs97 <- devs_fyear_all97
ewaa_m3$fcohort_devs97 <- devs_fcohort_all97
ewaa_m3$ewaa_fyear97 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fyear_devs97)
ewaa_m3$ewaa_fcohort97 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fcohort_devs97)
ewaa_m3$ewaa_fyear_fcoh97 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fyear_devs97 + ewaa_m3$fcohort_devs97)
head(ewaa_m3, n = 10)

# Selecciona el peso predicho con el efecto de año y cohorte solamente.
ewaa_m3_sim97 <- ewaa_m3[, c("real_year", "age", "ewaa_fyear_fcoh97")]
head(ewaa_m3_sim97, n = 800)

# Renombra las columnas al formato SS3 ewaa.ss
ewaa_m3_sim97 <- ewaa_m3_sim97 %>%
  rename(year = real_year, pred_weight = ewaa_fyear_fcoh97)
colnames(ewaa_m3_sim97) <- c("year", "age", "pred_weight")
head(ewaa_m3_sim97, n = 800)

# Convierte la columna de año a numérica si es necesario
ewaa_m3_sim97$year <- as.numeric(ewaa_m3_sim97$year)

# Ejecuta la función pad_weight_at_age e inspecciona la salida
padded_data97 <- pad_weight_at_age(ewaa_m3_sim97)
str(padded_data97)

# Escribe el archivo con el formato SS3
hakedataUSA:::write_wtatage_file(file = "ewaa_sim97.ss", data = padded_data97, maturity = maturity_at_age)

# Calcular la media y el intervalo de confianza de pred_weight para cada valor de age
mean_weight_by_age97 <- ewaa_m3_sim97 %>%
  group_by(age) %>%
  summarize(mean_pred_weight = mean(pred_weight, na.rm = TRUE),
            lower_ci = mean_pred_weight - qnorm(0.975) * sd(pred_weight) / sqrt(n()),
            upper_ci = mean_pred_weight + qnorm(0.975) * sd(pred_weight) / sqrt(n()))

# Crear el gráfico con intervalos de confianza
ewaa97<-ggplot(mean_weight_by_age97, aes(x = age, y = mean_pred_weight)) +
  geom_line() +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci), alpha = 0.2) +
  labs(x = "Age (years)", y = "Mean weight (kg)")

####################################################################
# Sim 98
set.seed(98)
devs_fyear98 <- rnorm(nT, mean = 0, sd = sd_fyear)
devs_fyear_fut98 <- rep(devs_fyear98, each = 16)
devs_fcohort98 <- rnorm(nT_23, mean = 0, sd = sd_fcohort)
devs_fcohort_past_future98 <- c(fcohort_estimate_values, devs_fcohort98)
estimate_fcohort_all98 <- data.frame(year = years_cohort, estimate = devs_fcohort_past_future98)
result_matrix98 <- matrix(0, nrow = 17, ncol = num_cols)
result_matrix98[1, ] <- format(as.integer(estimate_fcohort_all98$year), nsmall = 0)
result_matrix98[2, ] <- estimate_fcohort_all98$estimate
print(result_matrix98)
num_rows98 <- nrow(result_matrix98)
num_cols98 <- ncol(result_matrix98)
value_to_repeat98 <- result_matrix98[2, ]
start_row98 <- 3
start_col98 <- 2
end_row98 <- num_rows98
end_col98 <- num_cols98
for (i in start_row98:end_row98) {
  for (j in start_col98:end_col98) {
    if (j >= (i - start_row98 + start_col98)) {
      pattern_index98 <- (j - (i - start_row98 + start_col98)) %% 32 + 1
      result_matrix98[i, j] <- value_to_repeat98[pattern_index98]
    }
  }
}
print(result_matrix98)
extracted_values98 <- as.numeric(result_matrix98[2:17, 16:33])
extracted_values98 <- extracted_values98[-1]
devs_fcohort_all98 <- c(devs_fcohort_past, extracted_values98)
devs_fyear_all98 <- c(devs_fyear_past, devs_fyear_fut98)
ewaa_m3$fyear_devs98 <- devs_fyear_all98
ewaa_m3$fcohort_devs98 <- devs_fcohort_all98
ewaa_m3$ewaa_fyear98 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fyear_devs98)
ewaa_m3$ewaa_fcohort98 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fcohort_devs98)
ewaa_m3$ewaa_fyear_fcoh98 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fyear_devs98 + ewaa_m3$fcohort_devs98)
head(ewaa_m3, n = 10)

# Selecciona el peso predicho con el efecto de año y cohorte solamente.
ewaa_m3_sim98 <- ewaa_m3[, c("real_year", "age", "ewaa_fyear_fcoh98")]
head(ewaa_m3_sim98, n = 800)

# Renombra las columnas al formato SS3 ewaa.ss
ewaa_m3_sim98 <- ewaa_m3_sim98 %>%
  rename(year = real_year, pred_weight = ewaa_fyear_fcoh98)
colnames(ewaa_m3_sim98) <- c("year", "age", "pred_weight")
head(ewaa_m3_sim98, n = 800)

# Convierte la columna de año a numérica si es necesario
ewaa_m3_sim98$year <- as.numeric(ewaa_m3_sim98$year)

# Ejecuta la función pad_weight_at_age e inspecciona la salida
padded_data98 <- pad_weight_at_age(ewaa_m3_sim98)
str(padded_data98)

# Escribe el archivo con el formato SS3
hakedataUSA:::write_wtatage_file(file = "ewaa_sim98.ss", data = padded_data98, maturity = maturity_at_age)

# Calcular la media y el intervalo de confianza de pred_weight para cada valor de age
mean_weight_by_age98 <- ewaa_m3_sim98 %>%
  group_by(age) %>%
  summarize(mean_pred_weight = mean(pred_weight, na.rm = TRUE),
            lower_ci = mean_pred_weight - qnorm(0.975) * sd(pred_weight) / sqrt(n()),
            upper_ci = mean_pred_weight + qnorm(0.975) * sd(pred_weight) / sqrt(n()))

# Crear el gráfico con intervalos de confianza
ewaa98<-ggplot(mean_weight_by_age98, aes(x = age, y = mean_pred_weight)) +
  geom_line() +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci), alpha = 0.2) +
  labs(x = "Age (years)", y = "Mean Weight (kg)")
############################################################################
# Sim 99
set.seed(99)
devs_fyear99 <- rnorm(nT, mean = 0, sd = sd_fyear)
devs_fyear_fut99 <- rep(devs_fyear99, each = 16)
devs_fcohort99 <- rnorm(nT_23, mean = 0, sd = sd_fcohort)
devs_fcohort_past_future99 <- c(fcohort_estimate_values, devs_fcohort99)
estimate_fcohort_all99 <- data.frame(year = years_cohort, estimate = devs_fcohort_past_future99)
result_matrix99 <- matrix(0, nrow = 17, ncol = num_cols)
result_matrix99[1, ] <- format(as.integer(estimate_fcohort_all99$year), nsmall = 0)
result_matrix99[2, ] <- estimate_fcohort_all99$estimate
print(result_matrix99)
num_rows99 <- nrow(result_matrix99)
num_cols99 <- ncol(result_matrix99)
value_to_repeat99 <- result_matrix99[2, ]
start_row99 <- 3
start_col99 <- 2
end_row99 <- num_rows99
end_col99 <- num_cols99
for (i in start_row99:end_row99) {
  for (j in start_col99:end_col99) {
    if (j >= (i - start_row99 + start_col99)) {
      pattern_index99 <- (j - (i - start_row99 + start_col99)) %% 32 + 1
      result_matrix99[i, j] <- value_to_repeat99[pattern_index99]
    }
  }
}
print(result_matrix99)
extracted_values99 <- as.numeric(result_matrix99[2:17, 16:33])
extracted_values99 <- extracted_values99[-1]
devs_fcohort_all99 <- c(devs_fcohort_past, extracted_values99)
devs_fyear_all99 <- c(devs_fyear_past, devs_fyear_fut99)
ewaa_m3$fyear_devs99 <- devs_fyear_all99
ewaa_m3$fcohort_devs99 <- devs_fcohort_all99
ewaa_m3$ewaa_fyear99 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fyear_devs99)
ewaa_m3$ewaa_fcohort99 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fcohort_devs99)
ewaa_m3$ewaa_fyear_fcoh99 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fyear_devs99 + ewaa_m3$fcohort_devs99)
head(ewaa_m3, n = 10)

# Selecciona el peso predicho con el efecto de año y cohorte solamente.
ewaa_m3_sim99 <- ewaa_m3[, c("real_year", "age", "ewaa_fyear_fcoh99")]
head(ewaa_m3_sim99, n = 800)

# Renombra las columnas al formato SS3 ewaa.ss
ewaa_m3_sim99 <- ewaa_m3_sim99 %>%
  rename(year = real_year, pred_weight = ewaa_fyear_fcoh99)
colnames(ewaa_m3_sim99) <- c("year", "age", "pred_weight")
head(ewaa_m3_sim99, n = 800)

# Convierte la columna de año a numérica si es necesario
ewaa_m3_sim99$year <- as.numeric(ewaa_m3_sim99$year)

# Ejecuta la función pad_weight_at_age e inspecciona la salida
padded_data99 <- pad_weight_at_age(ewaa_m3_sim99)
str(padded_data99)

# Escribe el archivo con el formato SS3
hakedataUSA:::write_wtatage_file(file = "ewaa_sim99.ss", data = padded_data99, maturity = maturity_at_age)

# Calcular la media y el intervalo de confianza de pred_weight para cada valor de age
mean_weight_by_age99 <- ewaa_m3_sim99 %>%
  group_by(age) %>%
  summarize(mean_pred_weight = mean(pred_weight, na.rm = TRUE),
            lower_ci = mean_pred_weight - qnorm(0.975) * sd(pred_weight) / sqrt(n()),
            upper_ci = mean_pred_weight + qnorm(0.975) * sd(pred_weight) / sqrt(n()))

# Crear el gráfico con intervalos de confianza
ewaa99<-ggplot(mean_weight_by_age99, aes(x = age, y = mean_pred_weight)) +
  geom_line() +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci), alpha = 0.2) +
  labs(x = "Age (years)", y = "Mean Weight (kg)")
###########################################################################
# Sim 100
set.seed(100)
devs_fyear100 <- rnorm(nT, mean = 0, sd = sd_fyear)
devs_fyear_fut100 <- rep(devs_fyear100, each = 16)
devs_fcohort100 <- rnorm(nT_23, mean = 0, sd = sd_fcohort)
devs_fcohort_past_future100 <- c(fcohort_estimate_values, devs_fcohort100)
estimate_fcohort_all100 <- data.frame(year = years_cohort, estimate = devs_fcohort_past_future100)
result_matrix100 <- matrix(0, nrow = 17, ncol = num_cols)
result_matrix100[1, ] <- format(as.integer(estimate_fcohort_all100$year), nsmall = 0)
result_matrix100[2, ] <- estimate_fcohort_all100$estimate
print(result_matrix100)
num_rows100 <- nrow(result_matrix100)
num_cols100 <- ncol(result_matrix100)
value_to_repeat100 <- result_matrix100[2, ]
start_row100 <- 3
start_col100 <- 2
end_row100 <- num_rows100
end_col100 <- num_cols100
for (i in start_row100:end_row100) {
  for (j in start_col100:end_col100) {
    if (j >= (i - start_row100 + start_col100)) {
      pattern_index100 <- (j - (i - start_row100 + start_col100)) %% 32 + 1
      result_matrix100[i, j] <- value_to_repeat100[pattern_index100]
    }
  }
}
print(result_matrix100)
extracted_values100 <- as.numeric(result_matrix100[2:17, 16:33])
extracted_values100 <- extracted_values100[-1]
devs_fcohort_all100 <- c(devs_fcohort_past, extracted_values100)
devs_fyear_all100 <- c(devs_fyear_past, devs_fyear_fut100)
ewaa_m3$fyear_devs100 <- devs_fyear_all100
ewaa_m3$fcohort_devs100 <- devs_fcohort_all100
ewaa_m3$ewaa_fyear100 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fyear_devs100)
ewaa_m3$ewaa_fcohort100 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fcohort_devs100)
ewaa_m3$ewaa_fyear_fcoh100 <- exp(ewaa_m3$pred_weight) * exp(ewaa_m3$fyear_devs100 + ewaa_m3$fcohort_devs100)
head(ewaa_m3, n = 10)

# Selecciona el peso predicho con el efecto de año y cohorte solamente.
ewaa_m3_sim100 <- ewaa_m3[, c("real_year", "age", "ewaa_fyear_fcoh100")]
head(ewaa_m3_sim100, n = 800)

# Renombra las columnas al formato SS3 ewaa.ss
ewaa_m3_sim100 <- ewaa_m3_sim100 %>%
  rename(year = real_year, pred_weight = ewaa_fyear_fcoh100)
colnames(ewaa_m3_sim100) <- c("year", "age", "pred_weight")
head(ewaa_m3_sim100, n = 10)

# Ejecuta la función pad_weight_at_age e inspecciona la salida
padded_data100 <- pad_weight_at_age(ewaa_m3_sim100)
str(padded_data100)

# Escribe el archivo con el formato SS3
hakedataUSA:::write_wtatage_file(file = "ewaa_sim100.ss", data = padded_data100, maturity = maturity_at_age)

# Calcular la media y el intervalo de confianza de pred_weight para cada valor de age
mean_weight_by_age100 <- ewaa_m3_sim100 %>%
  group_by(age) %>%
  summarize(mean_pred_weight = mean(pred_weight, na.rm = TRUE),
            lower_ci = mean_pred_weight - qnorm(0.975) * sd(pred_weight) / sqrt(n()),
            upper_ci = mean_pred_weight + qnorm(0.975) * sd(pred_weight) / sqrt(n()))

# Crear el gráfico con intervalos de confianza
ggplot(mean_weight_by_age100, aes(x = age, y = mean_pred_weight)) +
  geom_line() +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci), alpha = 0.2) +
  labs(x = "Age (years)", y = "Mean Weight (kg)")

#############################################################################################
#ewaa Data processing (100sims)
############################################################################################

# ewaa year effect.
ewaa_fyear_pred_sel <- ewaa_m3[, c("real_year", "age", grep("^ewaa_fyear\\d+$", names(ewaa_m3), value = TRUE))]
ewaa_fyear_pred <- ewaa_fyear_pred_sel[, !grepl("ewaa_fyear_fcoh|fcohort_devs|fyear_devs|pred_weight|upper|lower|sd", names(ewaa_fyear_pred_sel))]
head(ewaa_fyear_pred,n=10)


# ewaa fcohort effect.
ewaa_fcohort_pred_sel <- ewaa_m3[, c("real_year", "age", grep("^ewaa_fcohort\\d+$", names(ewaa_m3), value = TRUE))]
ewaa_fcohort_pred <- ewaa_fcohort_pred_sel[, !grepl("ewaa_fyear_fcoh|ewaa_fyear|fcohort_devs|fyear_devs|pred_weight|upper|lower|sd", names(ewaa_fcohort_pred_sel))]
head(ewaa_fcohort_pred,n=2000)
ewaa_fcohort_pred1<-gather(ewaa_fcohort_pred,key="Sim", value="ewaa_fcohort", ewaa_fcohort1:ewaa_fcohort100)
head(ewaa_fcohort_pred1,n=50)

# ewaa fyear and fcohort effect.
ewaa_fyear_fcoh_pred <- ewaa_m3[, c("real_year", "age", grep("^ewaa_fyear_fcoh\\d+$", names(ewaa_m3), value = TRUE))]
head(ewaa_fyear_fcoh_pred_sel,n=10)


ewaa_fyear_fcoh_pred1<-gather(ewaa_fyear_fcoh_pred,key="Sim", value="ewaa_fyear_fcoh", ewaa_fyear_fcoh1:ewaa_fyear_fcoh100)
head(ewaa_fyear_fcoh_pred1,n=40)



#cohort effect on predict weight at age
ewaa_fyear<-ggplot(ewaa_fcohort_pred1, aes(x = age, y = c(ewaa_fcohort),color=Sim)) +
  geom_line() +
  facet_wrap(~real_year, scales = "free") +
  labs(x = "Age", y = "Weight") +
  ggtitle("Weight-at-age + SD fcohort")+theme(axis.title = element_text(size = 17),
                                                        axis.text = element_text(size = 15),
                                                        legend.title = element_text(size = 17),
                                                        legend.text = element_text(size = 15),
                                                        axis.text.x = element_text(size=8),
                                                        axis.text.y = element_text(size=8),
                                                        legend.position = "none",
                                                        axis.title.x = element_blank(),
                                                        axis.ticks=element_blank())



#fyear+fcohort effect on predict weight
ewaa_ewaa_fyear_fco<-ggplot(ewaa_fyear_fcoh_pred1, aes(x = age, y = c(ewaa_fyear_fcoh),color=Sim)) +
  geom_line() +
  facet_wrap(~real_year, scales = "free") +
  labs(x = "Age", y = "Weight") +
  ggtitle("Weight-at-age + SD fyear +SD fcohort")+theme(axis.title = element_text(size = 17),
                                            axis.text = element_text(size = 15),
                                            legend.title = element_text(size = 17),
                                            legend.text = element_text(size = 15),
                                            axis.text.x = element_text(size=8),
                                            axis.text.y = element_text(size=8),
                                            legend.position = "none",
                                            axis.title.x = element_blank(),
                                            axis.ticks=element_blank())

#####################################################################################################
 # Assuming 'random_eff_val_m1_final' is your dataframe that you want to export
 # Replace 'path_to_file.csv' with the desired path and filename for your CSV file

 # Export the dataframe to CSV
 write.csv(ewaa_m3, "C:/GitHub/hake-assessment//ewaa_m3.csv", row.names = FALSE)

 ############################################################################################
 # Specify the path to your CSV file
 dir <- "C:/GitHub/hake-assessment/ewaa_m3.csv"

 # Read the CSV file into a dataframe
 my_data <- read.csv(dir,sep = ",")

 # View the structure of the dataframe
 my_data1<-data.frame(my_data)
 head(my_data1,n=10)
color<-c()
for (i in 2:100) {
  ewaa_fyear <- ewaa_fyear + geom_line(aes_string(y = paste0("ewaa_fyear", i)), color = rainbow(100)[i])
}

 ewaa_fyear <- ggplot(my_data, aes(x = age)) +
   geom_line(aes(y = ewaa_fyear1), color = "blue") +
   geom_line(aes(y = ewaa_fyear2), color = "red") +
   geom_line(aes(y = ewaa_fyear3), color = "green") +
   geom_line(aes(y = ewaa_fyear4), color = "purple") +
   geom_line(aes(y = ewaa_fyear5), color = "orange") +
   geom_line(aes(y = ewaa_fyear6), color = "cyan") +
   geom_line(aes(y = ewaa_fyear7), color = "magenta") +
   geom_line(aes(y = ewaa_fyear8), color = "black") +
   # Continúa agregando más líneas para ewaa_fyear9 hasta ewaa_fyear100
   facet_wrap(~real_year, scales = "free") +
   labs(x = "Age", y = "Ewaa_fyear") +
   ggtitle("Ewaa_fyear") +
   theme(axis.title = element_text(size = 17),
         axis.text = element_text(size = 15),
         legend.title = element_text(size = 17),
         legend.text = element_text(size = 15),
         axis.text.x = element_text(size = 8),
         axis.text.y = element_text(size = 8),
         legend.position = "none",
         axis.title.x = element_blank(),
         axis.ticks = element_blank())

 print(ewaa_fyear)




 output_directory <- "C:/GitHub/hake-assessment/Tables_figures_Sandra/"
 ewaa_fyear1 <- paste0(output_directory, "ewaa_fyear_with_sim.png")
 ewaa_fyear<-ggplot(my_data1, aes(x = age, y = c(ewaa_fyear2))) +
   geom_point() +
   facet_wrap(~real_year, scales = "free") +
   labs(x = "Age", y = "Weight") +
   ggtitle("Weight-at-age + SD fyear")+theme(axis.title = element_text(size = 17),
                                            axis.text = element_text(size = 15),
                                            legend.title = element_text(size = 17),
                                            legend.text = element_text(size = 15),
                                            axis.text.x = element_text(size=8),
                                            axis.text.y = element_text(size=8),
                                            legend.position = "none",
                                            axis.title.x = element_blank(),
                                            axis.ticks=element_blank())
 # Save the plot to the specified directory
 ggsave(ewaa_fyear1, plot = ewaa_fyear, width = 12, height = 8, units = "in", dpi = 300)

 #Weight at age predicted with simulations
 ewaa_fcohort1 <- paste0(output_directory, "ewaa_fcohort_with_sim.png")
 ewaa_fcohort<-ggplot(my_data1, aes(x = age, y = ewaa_fcohort)) +
   geom_point() +
   facet_wrap(~real_year, scales = "free") +
   labs(x = "Age", y = "Weight") +
   ggtitle("Weight-at-age + SD fcohort")+theme(axis.title = element_text(size = 17),
                                               axis.text = element_text(size = 15),
                                               legend.title = element_text(size = 17),
                                               legend.text = element_text(size = 15),
                                               axis.text.x = element_text(size=8),
                                               axis.text.y = element_text(size=8),
                                               legend.position = "none",
                                               axis.title.x = element_blank(),
                                               axis.ticks=element_blank())
 ggsave(ewaa_fcohort1, plot = ewaa_fcohort, width = 12, height = 8, units = "in", dpi = 300)

 ewaa__fyear_fcoh1 <- paste0(output_directory, "ewaa_fyear_fcoh_with_sim.png")
  ewaa_fyear_fcoh<-ggplot(my_data1, aes(x = age, y = ewaa_fyear_fcoh)) +
   geom_point() +
   facet_wrap(~real_year, scales = "free") +
   labs(x = "Age", y = "Weight") +
   ggtitle("Weight-at-age + SD fyear +SD fcohort")+theme(axis.title = element_text(size = 17),
                                                         axis.text = element_text(size = 15),
                                                         legend.title = element_text(size = 17),
                                                         legend.text = element_text(size = 15),
                                                         axis.text.x = element_text(size=8),
                                                         axis.text.y = element_text(size=8),
                                                         legend.position = "none",
                                                         axis.title.x = element_blank(),
                                                         axis.ticks=element_blank())
  ggsave(ewaa__fyear_fcoh1, plot = ewaa_fyear_fcoh, width = 12, height = 8, units = "in", dpi = 300)



############################################################################################
# 4 part: Plotting
# Plot weight predict vs age, with variability on fyear
 ewaa_fyear_eff_age<-ggplot(my_data1, aes(age, ewaa_fcohort)) +
   geom_point() + facet_wrap(~real_year)+labs(x = "Age", y = "Weight")
###########################################################################################
 #predict weight at age fyear effect by simulation plot
 # Define the directory where you want to save the plot
 output_directory <- "C:/GitHub/hake-assessment/Tables_figures_Sandra/ewaa_fyear_sim/"
 #Sim 2
 output_file1 <- paste0(output_directory, "ewaa_fyear_sim1.png")
 ewaa_fyear0<-ggplot(ewaa_m3, aes(x = factor(real_year), y = ewaa_fyear, color = factor(age),
                                   group = factor(age))) +
   geom_text(aes(label=round(age,2)), size = 4.5) +
   geom_line(alpha = 0.85) +ylim(0, 2.5)+
   theme_bw() +
   labs(x = "Year", y = "Weight") + annotate("text", x = 62, y = 2.5, parse = FALSE,label = "Sim 1",size=5)+
   theme(axis.title = element_text(size = 17),
         axis.text = element_text(size = 15),
         legend.title = element_text(size = 17),
         legend.text = element_text(size = 15),
         axis.text.x = element_text(angle = 90,hjust = 1,size=8),
         legend.position = "none",
         axis.title.x = element_blank(),
         axis.ticks=element_blank())

 # Save the plot to the specified directory
 ggsave(output_file1, plot = ewaa_fyear0, width = 12, height = 8, units = "in", dpi = 300)

#Sim 2
 output_file2 <- paste0(output_directory, "ewaa_fyear_sim2.png")
 ewaa_fyear1<-ggplot(ewaa_m3, aes(x = factor(real_year), y = ewaa_fyear1, color = factor(age),
                                  group = factor(age))) +
   geom_text(aes(label=round(age,2)), size = 4.5) +
   geom_line(alpha = 0.85) +ylim(0, 2.5)+
   theme_bw() +
   labs(x = "Year", y = "Weight") + annotate("text", x = 62, y = 2.5, parse = FALSE,label = "Sim 2",size=5)+
   theme(axis.title = element_text(size = 17),
         axis.text = element_text(size = 15),
         legend.title = element_text(size = 17),
         legend.text = element_text(size = 15),
         axis.text.x = element_text(angle = 90,hjust = 1,size=8),
         legend.position = "none",
         axis.title.x = element_blank(),
         axis.ticks=element_blank())

 # Save the plot to the specified directory
 ggsave(output_file2, plot = ewaa_fyear1, width = 12, height = 8, units = "in", dpi = 300)


 #Sim 3
 output_file3 <- paste0(output_directory, "ewaa_fyear_sim3.png")
 ewaa_fyear2<-ggplot(ewaa_m3, aes(x = factor(real_year), y = ewaa_fyear2, color = factor(age),
                                  group = factor(age))) +
   geom_text(aes(label=round(age,2)), size = 4.5) +
   geom_line(alpha = 0.85) +ylim(0, 2.5)+
   theme_bw() +
   labs(x = "Year", y = "Weight") + annotate("text", x = 62, y = 2.5, parse = FALSE,label = "Sim 3",size=5)+
   theme(axis.title = element_text(size = 17),
         axis.text = element_text(size = 15),
         legend.title = element_text(size = 17),
         legend.text = element_text(size = 15),
         axis.text.x = element_text(angle = 90,hjust = 1,size=8),
         legend.position = "none",
         axis.title.x = element_blank(),
         axis.ticks=element_blank())
 # Save the plot to the specified directory
 ggsave(output_file3, plot = ewaa_fyear2, width = 12, height = 8, units = "in", dpi = 300)


 #sim 4
 output_file4 <- paste0(output_directory, "ewaa_fyear_sim4.png")
 ewaa_fyear3<-ggplot(ewaa_m3, aes(x = factor(real_year), y = ewaa_fyear3, color = factor(age),
                                  group = factor(age))) +
   geom_text(aes(label=round(age,2)), size = 4.5) +
   geom_line(alpha = 0.85) +ylim(0, 2.5)+
   theme_bw() +
   labs(x = "Year", y = "Weight") + annotate("text", x = 62, y = 2.5, parse = FALSE,label = "Sim 4",size=5)+
   theme(axis.title = element_text(size = 17),
         axis.text = element_text(size = 15),
         legend.title = element_text(size = 17),
         legend.text = element_text(size = 15),
         axis.text.x = element_text(angle = 90,hjust = 1,size=8),
         legend.position = "none",
         axis.title.x = element_blank(),
         axis.ticks=element_blank())
 # Save the plot to the specified directory
 ggsave(output_file4, plot = ewaa_fyear3, width = 12, height = 8, units = "in", dpi = 300)

 #Sim 5
 output_file5 <- paste0(output_directory, "ewaa_fyear_sim5.png")
 ewaa_fyear4<-ggplot(ewaa_m3, aes(x = factor(real_year), y = ewaa_fyear4, color = factor(age),
                                  group = factor(age))) +
   geom_text(aes(label=round(age,2)), size = 4.5) +
   geom_line(alpha = 0.85) +ylim(0, 2.5)+
   theme_bw() +
   labs(x = "Year", y = "Weight") + annotate("text", x = 62, y = 2.5, parse = FALSE,label = "Sim 5",size=5)+
   theme(axis.title = element_text(size = 17),
         axis.text = element_text(size = 15),
         legend.title = element_text(size = 17),
         legend.text = element_text(size = 15),
         axis.text.x = element_text(angle = 90,hjust = 1,size=8),
         legend.position = "none",
         axis.title.x = element_blank(),
         axis.ticks=element_blank())
 # Save the plot to the specified directory
 ggsave(output_file5, plot = ewaa_fyear4, width = 12, height = 8, units = "in", dpi = 300)

 #Sim 6
 output_file6 <- paste0(output_directory, "ewaa_fyear_sim6.png")
 ewaa_fyear5<-ggplot(ewaa_m3, aes(x = factor(real_year), y = ewaa_fyear5, color = factor(age),
                                  group = factor(age))) +
   geom_text(aes(label=round(age,2)), size = 4.5) +
   geom_line(alpha = 0.85) +ylim(0, 2.5)+
   theme_bw() +
   labs(x = "Year", y = "Weight") + annotate("text", x = 62, y = 2.5, parse = FALSE,label = "Sim 6",size=5)+
   theme(axis.title = element_text(size = 17),
         axis.text = element_text(size = 15),
         legend.title = element_text(size = 17),
         legend.text = element_text(size = 15),
         axis.text.x = element_text(angle = 90,hjust = 1,size=8),
         legend.position = "none",
         axis.title.x = element_blank(),
         axis.ticks=element_blank())
 # Save the plot to the specified directory
 ggsave(output_file6, plot = ewaa_fyear5, width = 12, height = 8, units = "in", dpi = 300)

 #Sim 7
 output_file7 <- paste0(output_directory, "ewaa_fyear_sim7.png")
 ewaa_fyear6<-ggplot(ewaa_m3, aes(x = factor(real_year), y = ewaa_fyear6, color = factor(age),
                                  group = factor(age))) +
   geom_text(aes(label=round(age,2)), size = 4.5) +
   geom_line(alpha = 0.85) +ylim(0, 2.5)+
   theme_bw() +
   labs(x = "Year", y = "Weight") + annotate("text", x = 62, y = 2.5, parse = FALSE,label = "Sim 7",size=5)+
   theme(axis.title = element_text(size = 17),
         axis.text = element_text(size = 15),
         legend.title = element_text(size = 17),
         legend.text = element_text(size = 15),
         axis.text.x = element_text(angle = 90,hjust = 1,size=8),
         legend.position = "none",
         axis.title.x = element_blank(),
         axis.ticks=element_blank())
 # Save the plot to the specified directory
 ggsave(output_file7, plot = ewaa_fyear6, width = 12, height = 8, units = "in", dpi = 300)

 #Sim 8
 output_file8 <- paste0(output_directory, "ewaa_fyear_sim8.png")
 ewaa_fyear7<-ggplot(ewaa_m3, aes(x = factor(real_year), y = ewaa_fyear7, color = factor(age),
                                  group = factor(age))) +
   geom_text(aes(label=round(age,2)), size = 4.5) +
   geom_line(alpha = 0.85) +ylim(0, 2.5)+
   theme_bw() +
   labs(x = "Year", y = "Weight") + annotate("text", x = 62, y = 2.5, parse = FALSE,label = "Sim 8",size=5)+
   theme(axis.title = element_text(size = 17),
         axis.text = element_text(size = 15),
         legend.title = element_text(size = 17),
         legend.text = element_text(size = 15),
         axis.text.x = element_text(angle = 90,hjust = 1,size=8),
         legend.position = "none",
         axis.title.x = element_blank(),
         axis.ticks=element_blank())

 # Save the plot to the specified directory
 ggsave(output_file8, plot = ewaa_fyear7, width = 12, height = 8, units = "in", dpi = 300)

 # Sim 9
 output_file9 <- paste0(output_directory, "ewaa_fyear_sim9.png")
 ewaa_fyear8<-ggplot(ewaa_m3, aes(x = factor(real_year), y = ewaa_fyear8, color = factor(age),
                                  group = factor(age))) +
   geom_text(aes(label=round(age,2)), size = 4.5) +
   geom_line(alpha = 0.85) +ylim(0, 2.5)+
   theme_bw() +
   labs(x = "Year", y = "Weight") + annotate("text", x = 62, y = 2.5, parse = FALSE,label = "Sim 9",size=5)+
   theme(axis.title = element_text(size = 17),
         axis.text = element_text(size = 15),
         legend.title = element_text(size = 17),
         legend.text = element_text(size = 15),
         axis.text.x = element_text(angle = 90,hjust = 1,size=8),
         legend.position = "none",
         axis.title.x = element_blank(),
         axis.ticks=element_blank())
 # Save the plot to the specified directory
 ggsave(output_file9, plot = ewaa_fyear8, width = 12, height = 8, units = "in", dpi = 300)

 # Sim 10
 output_file10 <- paste0(output_directory, "ewaa_fyear_sim10.png")
 ewaa_fyear9<-ggplot(ewaa_m3, aes(x = factor(real_year), y = ewaa_fyear9, color = factor(age),
                                  group = factor(age))) +
   geom_text(aes(label=round(age,2)), size = 4.5) +
   geom_line(alpha = 0.85) +ylim(0, 2.5)+
   theme_bw() +
   labs(x = "Year", y = "Weight") + annotate("text", x = 62, y = 2.5, parse = FALSE,label = "Sim 10",size=5)+
   theme(axis.title = element_text(size = 17),
         axis.text = element_text(size = 15),
         legend.title = element_text(size = 17),
         legend.text = element_text(size = 15),
         axis.text.x = element_text(angle = 90,hjust = 1,size=8),
         legend.position = "none",
         axis.title.x = element_blank(),
         axis.ticks=element_blank())
 # Save the plot to the specified directory
 ggsave(output_file10, plot = ewaa_fyear9, width = 12, height = 8, units = "in", dpi = 300)

 # Sim 11
 output_file11 <- paste0(output_directory, "ewaa_fyear_sim11.png")
 ewaa_fyear10<-ggplot(ewaa_m3, aes(x = factor(real_year), y = ewaa_fyear10, color = factor(age),
                                  group = factor(age))) +
   geom_text(aes(label=round(age,2)), size = 4.5) +
   geom_line(alpha = 0.85) +ylim(0, 2.5)+
   theme_bw() +
   labs(x = "Year", y = "Weight") + annotate("text", x = 62, y = 2.5, parse = FALSE,label = "Sim 11",size=5)+
   theme(axis.title = element_text(size = 17),
         axis.text = element_text(size = 15),
         legend.title = element_text(size = 17),
         legend.text = element_text(size = 15),
         axis.text.x = element_text(angle = 90,hjust = 1,size=8),
         legend.position = "none",
         axis.title.x = element_blank(),
         axis.ticks=element_blank())
 # Save the plot to the specified directory
 ggsave(output_file11, plot = ewaa_fyear9, width = 12, height = 8, units = "in", dpi = 300)

 ##########################################################################################
 #predict weight at age fcohort effect plot
 #
 #########################################################################################
 output_directory_fcoh <- "C:/GitHub/hake-assessment/Tables_figures_Sandra/ewaa_fcohort_sim/"
 #sim 1
 output_file_fcoh0 <- paste0(output_directory_fcoh, "ewaa_fcohort_sim1.png")
 ewaa_fcohort1<-ggplot(ewaa_m3, aes(x = factor(real_year), y = ewaa_fcohort, color = factor(age),
                                 group = factor(age))) +
   geom_text(aes(label=round(age,2)), size = 4.5) +
   geom_line(alpha = 0.85) +ylim(0, 2.5)+
   theme_bw() +
   labs(x = "Year", y = "Weight") + annotate("text", x = 62, y = 2.5, parse = FALSE,label = "Sim 1",size=5)+
   theme(axis.title = element_text(size = 17),
         axis.text = element_text(size = 15),
         legend.title = element_text(size = 17),
         legend.text = element_text(size = 15),
         axis.text.x = element_text(angle = 90,hjust = 1,size=8),
         legend.position = "none",
         axis.title.x = element_blank(),
         axis.ticks=element_blank())
 # Save the plot to the specified directory
 ggsave(output_file_fcoh0, plot = ewaa_fcohort1, width = 12, height = 8, units = "in", dpi = 300)


 #sim 2
 output_file_fcoh2 <- paste0(output_directory_fcoh, "ewaa_fcohort_sim2.png")
 ewaa_fcohort2<-ggplot(ewaa_m3, aes(x = factor(real_year), y = ewaa_fcohort1, color = factor(age),
                                    group = factor(age))) +
   geom_text(aes(label=round(age,2)), size = 4.5) +
   geom_line(alpha = 0.85) +ylim(0, 2.5)+
   theme_bw() +
   labs(x = "Year", y = "Weight") + annotate("text", x = 62, y = 2.5, parse = FALSE,label = "Sim 2",size=5)+
   theme(axis.title = element_text(size = 17),
         axis.text = element_text(size = 15),
         legend.title = element_text(size = 17),
         legend.text = element_text(size = 15),
         axis.text.x = element_text(angle = 90,hjust = 1,size=8),
         legend.position = "none",
         axis.title.x = element_blank(),
         axis.ticks=element_blank())
 # Save the plot to the specified directory
 ggsave(output_file_fcoh2, plot = ewaa_fcohort2, width = 12, height = 8, units = "in", dpi = 300)

 #sim 3
 output_file_fcoh3 <- paste0(output_directory_fcoh, "ewaa_fcohort_sim3.png")
 ewaa_fcohort3<-ggplot(ewaa_m3, aes(x = factor(real_year), y = ewaa_fcohort2, color = factor(age),
                                    group = factor(age))) +
   geom_text(aes(label=round(age,2)), size = 4.5) +
   geom_line(alpha = 0.85) +ylim(0, 2.5)+
   theme_bw() +
   labs(x = "Year", y = "Weight") + annotate("text", x = 62, y = 2.5, parse = FALSE,label = "Sim 3",size=5)+
   theme(axis.title = element_text(size = 17),
         axis.text = element_text(size = 15),
         legend.title = element_text(size = 17),
         legend.text = element_text(size = 15),
         axis.text.x = element_text(angle = 90,hjust = 1,size=8),
         legend.position = "none",
         axis.title.x = element_blank(),
         axis.ticks=element_blank())
 # Save the plot to the specified directory
 ggsave(output_file_fcoh3, plot = ewaa_fcohort3, width = 12, height = 8, units = "in", dpi = 300)


 #sim 4
 output_file_fcoh4 <- paste0(output_directory_fcoh, "ewaa_fcohort_sim4.png")
 ewaa_fcohort4<-ggplot(ewaa_m3, aes(x = factor(real_year), y = ewaa_fcohort3, color = factor(age),
                                    group = factor(age))) +
   geom_text(aes(label=round(age,2)), size = 4.5) +
   geom_line(alpha = 0.85) +ylim(0, 2.5)+
   theme_bw() +
   labs(x = "Year", y = "Weight") + annotate("text", x = 62, y = 2.5, parse = FALSE,label = "Sim 4",size=5)+
   theme(axis.title = element_text(size = 17),
         axis.text = element_text(size = 15),
         legend.title = element_text(size = 17),
         legend.text = element_text(size = 15),
         axis.text.x = element_text(angle = 90,hjust = 1,size=8),
         legend.position = "none",
         axis.title.x = element_blank(),
         axis.ticks=element_blank())
 # Save the plot to the specified directory
 ggsave(output_file_fcoh4, plot = ewaa_fcohort4, width = 12, height = 8, units = "in", dpi = 300)


 #sim 5
 output_file_fcoh5 <- paste0(output_directory_fcoh, "ewaa_fcohort_sim5.png")
 ewaa_fcohort5<-ggplot(ewaa_m3, aes(x = factor(real_year), y = ewaa_fcohort4, color = factor(age),
                                    group = factor(age))) +
   geom_text(aes(label=round(age,2)), size = 4.5) +
   geom_line(alpha = 0.85) +ylim(0, 2.5)+
   theme_bw() +
   labs(x = "Year", y = "Weight") + annotate("text", x = 62, y = 2.5, parse = FALSE,label = "Sim 5",size=5)+
   theme(axis.title = element_text(size = 17),
         axis.text = element_text(size = 15),
         legend.title = element_text(size = 17),
         legend.text = element_text(size = 15),
         axis.text.x = element_text(angle = 90,hjust = 1,size=8),
         legend.position = "none",
         axis.title.x = element_blank(),
         axis.ticks=element_blank())
 # Save the plot to the specified directory
 ggsave(output_file_fcoh5, plot = ewaa_fcohort5, width = 12, height = 8, units = "in", dpi = 300)


 #sim 6
 output_file_fcoh6 <- paste0(output_directory_fcoh, "ewaa_fcohort_sim6.png")
 ewaa_fcohort6<-ggplot(ewaa_m3, aes(x = factor(real_year), y = ewaa_fcohort5, color = factor(age),
                                    group = factor(age))) +
   geom_text(aes(label=round(age,2)), size = 4.5) +
   geom_line(alpha = 0.85) +ylim(0, 2.5)+
   theme_bw() +
   labs(x = "Year", y = "Weight") + annotate("text", x = 62, y = 2.5, parse = FALSE,label = "Sim 6",size=5)+
   theme(axis.title = element_text(size = 17),
         axis.text = element_text(size = 15),
         legend.title = element_text(size = 17),
         legend.text = element_text(size = 15),
         axis.text.x = element_text(angle = 90,hjust = 1,size=8),
         legend.position = "none",
         axis.title.x = element_blank(),
         axis.ticks=element_blank())
 # Save the plot to the specified directory
 ggsave(output_file_fcoh6, plot = ewaa_fcohort6, width = 12, height = 8, units = "in", dpi = 300)


 #sim 7
 output_file_fcoh7 <- paste0(output_directory_fcoh, "ewaa_fcohort_sim7.png")
 ewaa_fcohort7<-ggplot(ewaa_m3, aes(x = factor(real_year), y = ewaa_fcohort6, color = factor(age),
                                    group = factor(age))) +
   geom_text(aes(label=round(age,2)), size = 4.5) +
   geom_line(alpha = 0.85) +ylim(0, 2.5)+
   theme_bw() +
   labs(x = "Year", y = "Weight") + annotate("text", x = 62, y = 2.5, parse = FALSE,label = "Sim 7",size=5)+
   theme(axis.title = element_text(size = 17),
         axis.text = element_text(size = 15),
         legend.title = element_text(size = 17),
         legend.text = element_text(size = 15),
         axis.text.x = element_text(angle = 90,hjust = 1,size=8),
         legend.position = "none",
         axis.title.x = element_blank(),
         axis.ticks=element_blank())
 # Save the plot to the specified directory
 ggsave(output_file_fcoh7, plot = ewaa_fcohort7, width = 12, height = 8, units = "in", dpi = 300)

 #sim 8
 output_file_fcoh8 <- paste0(output_directory_fcoh, "ewaa_fcohort_sim8.png")
 ewaa_fcohort8<-ggplot(ewaa_m3, aes(x = factor(real_year), y = ewaa_fcohort7, color = factor(age),
                                    group = factor(age))) +
   geom_text(aes(label=round(age,2)), size = 4.5) +
   geom_line(alpha = 0.85) +ylim(0, 2.5)+
   theme_bw() +
   labs(x = "Year", y = "Weight") + annotate("text", x = 62, y = 2.5, parse = FALSE,label = "Sim 8",size=5)+
   theme(axis.title = element_text(size = 17),
         axis.text = element_text(size = 15),
         legend.title = element_text(size = 17),
         legend.text = element_text(size = 15),
         axis.text.x = element_text(angle = 90,hjust = 1,size=8),
         legend.position = "none",
         axis.title.x = element_blank(),
         axis.ticks=element_blank())
 # Save the plot to the specified directory
 ggsave(output_file_fcoh8, plot = ewaa_fcohort8, width = 12, height = 8, units = "in", dpi = 300)

 #sim 9
 output_file_fcoh9 <- paste0(output_directory_fcoh, "ewaa_fcohort_sim9.png")
 ewaa_fcohort9<-ggplot(ewaa_m3, aes(x = factor(real_year), y = ewaa_fcohort8, color = factor(age),
                                    group = factor(age))) +
   geom_text(aes(label=round(age,2)), size = 4.5) +
   geom_line(alpha = 0.85) +ylim(0, 2.5)+
   theme_bw() +
   labs(x = "Year", y = "Weight") + annotate("text", x = 62, y = 2.5, parse = FALSE,label = "Sim 9",size=5)+
   theme(axis.title = element_text(size = 17),
         axis.text = element_text(size = 15),
         legend.title = element_text(size = 17),
         legend.text = element_text(size = 15),
         axis.text.x = element_text(angle = 90,hjust = 1,size=8),
         legend.position = "none",
         axis.title.x = element_blank(),
         axis.ticks=element_blank())
 # Save the plot to the specified directory
 ggsave(output_file_fcoh9, plot = ewaa_fcohort9, width = 12, height = 8, units = "in", dpi = 300)

#Sim 10
 output_file_fcoh10 <- paste0(output_directory_fcoh, "ewaa_fcohort_sim10.png")
 ewaa_fcohort10<-ggplot(ewaa_m3, aes(x = factor(real_year), y = ewaa_fcohort9, color = factor(age),
                                    group = factor(age))) +
   geom_text(aes(label=round(age,2)), size = 4.5) +
   geom_line(alpha = 0.85) +ylim(0, 2.5)+
   theme_bw() +
   labs(x = "Year", y = "Weight") + annotate("text", x = 62, y = 2.5, parse = FALSE,label = "Sim 10",size=5)+
   theme(axis.title = element_text(size = 17),
         axis.text = element_text(size = 15),
         legend.title = element_text(size = 17),
         legend.text = element_text(size = 15),
         axis.text.x = element_text(angle = 90,hjust = 1,size=8),
         legend.position = "none",
         axis.title.x = element_blank(),
         axis.ticks=element_blank())
 # Save the plot to the specified directory
 ggsave(output_file_fcoh10, plot = ewaa_fcohort10, width = 12, height = 8, units = "in", dpi = 300)


 #Sim 11
 output_file_fcoh11 <- paste0(output_directory_fcoh, "ewaa_fcohort_sim11.png")
 ewaa_fcohort11<-ggplot(ewaa_m3, aes(x = factor(real_year), y = ewaa_fcohort10, color = factor(age),
                                     group = factor(age))) +
   geom_text(aes(label=round(age,2)), size = 4.5) +
   geom_line(alpha = 0.85) +ylim(0, 2.5)+
   theme_bw() +
   labs(x = "Year", y = "Weight") + annotate("text", x = 62, y = 2.5, parse = FALSE,label = "Sim 11",size=5)+
   theme(axis.title = element_text(size = 17),
         axis.text = element_text(size = 15),
         legend.title = element_text(size = 17),
         legend.text = element_text(size = 15),
         axis.text.x = element_text(angle = 90,hjust = 1,size=8),
         legend.position = "none",
         axis.title.x = element_blank(),
         axis.ticks=element_blank())
 # Save the plot to the specified directory
 ggsave(output_file_fcoh11, plot = ewaa_fcohort10, width = 12, height = 8, units = "in", dpi = 300)

 #############################################################################################
# predict weight at age fcohort and fyear effect plot
###############################################################################################
 #predict weight at age fcohort effect plot

 output_directory_fyear_fcoh <- "C:/GitHub/hake-assessment/Tables_figures_Sandra/ewaa_fyear_fcohort_sim/"
 #sim 1
 output_file_fyear_fcoh1 <- paste0(output_directory_fyear_fcoh, "ewaa_fyear_fcohort_sim1.png")
 ewaa_fyear_fcohort1<-ggplot(ewaa_m3, aes(x = factor(real_year), y = ewaa_fyear_fcoh, color = factor(age),
                                    group = factor(age))) +
   geom_text(aes(label=round(age,2)), size = 4.5) +
   geom_line(alpha = 0.85) +ylim(0, 2.5)+
   theme_bw() +
   labs(x = "Year", y = "Weight") + annotate("text", x = 62, y = 2.5, parse = FALSE,label = "Sim 1",size=5)+
   theme(axis.title = element_text(size = 17),
         axis.text = element_text(size = 15),
         legend.title = element_text(size = 17),
         legend.text = element_text(size = 15),
         axis.text.x = element_text(angle = 90,hjust = 1,size=8),
         legend.position = "none",
         axis.title.x = element_blank(),
         axis.ticks=element_blank())
 # Save the plot to the specified directory
 ggsave(output_file_fyear_fcoh1, plot = ewaa_fyear_fcohort1, width = 12, height = 8, units = "in", dpi = 300)


 #sim 2
 output_file_fyear_fcoh2 <- paste0(output_directory_fyear_fcoh, "ewaa_fyear_fcohort_sim2.png")
 ewaa_fyear_fcohort2<-ggplot(ewaa_m3, aes(x = factor(real_year), y = ewaa_fyear_fcoh1, color = factor(age),
                                          group = factor(age))) +
   geom_text(aes(label=round(age,2)), size = 4.5) +
   geom_line(alpha = 0.85) +ylim(0, 2.5)+
   theme_bw() +
   labs(x = "Year", y = "Weight") + annotate("text", x = 62, y = 2.5, parse = FALSE,label = "Sim 2",size=5)+
   theme(axis.title = element_text(size = 17),
         axis.text = element_text(size = 15),
         legend.title = element_text(size = 17),
         legend.text = element_text(size = 15),
         axis.text.x = element_text(angle = 90,hjust = 1,size=8),
         legend.position = "none",
         axis.title.x = element_blank(),
         axis.ticks=element_blank())
 # Save the plot to the specified directory
 ggsave(output_file_fyear_fcoh2, plot = ewaa_fyear_fcohort2, width = 12, height = 8, units = "in", dpi = 300)

 #sim 3
 output_file_fyear_fcoh3 <- paste0(output_directory_fyear_fcoh, "ewaa_fyear_fcohort_sim3.png")
 ewaa_fyear_fcohort3<-ggplot(ewaa_m3, aes(x = factor(real_year), y = ewaa_fyear_fcoh2, color = factor(age),
                                          group = factor(age))) +
   geom_text(aes(label=round(age,2)), size = 4.5) +
   geom_line(alpha = 0.85) +ylim(0, 2.5)+
   theme_bw() +
   labs(x = "Year", y = "Weight") + annotate("text", x = 62, y = 2.5, parse = FALSE,label = "Sim 3",size=5)+
   theme(axis.title = element_text(size = 17),
         axis.text = element_text(size = 15),
         legend.title = element_text(size = 17),
         legend.text = element_text(size = 15),
         axis.text.x = element_text(angle = 90,hjust = 1,size=8),
         legend.position = "none",
         axis.title.x = element_blank(),
         axis.ticks=element_blank())
 # Save the plot to the specified directory
 ggsave(output_file_fyear_fcoh3, plot = ewaa_fyear_fcohort3, width = 12, height = 8, units = "in", dpi = 300)

 #sim 4
 output_file_fyear_fcoh4 <- paste0(output_directory_fyear_fcoh, "ewaa_fyear_fcohort_sim4.png")
 ewaa_fyear_fcohort4<-ggplot(ewaa_m3, aes(x = factor(real_year), y = ewaa_fyear_fcoh3, color = factor(age),
                                          group = factor(age))) +
   geom_text(aes(label=round(age,2)), size = 4.5) +
   geom_line(alpha = 0.85) +ylim(0, 2.5)+
   theme_bw() +
   labs(x = "Year", y = "Weight") + annotate("text", x = 62, y = 2.5, parse = FALSE,label = "Sim 4",size=5)+
   theme(axis.title = element_text(size = 17),
         axis.text = element_text(size = 15),
         legend.title = element_text(size = 17),
         legend.text = element_text(size = 15),
         axis.text.x = element_text(angle = 90,hjust = 1,size=8),
         legend.position = "none",
         axis.title.x = element_blank(),
         axis.ticks=element_blank())
 # Save the plot to the specified directory
 ggsave(output_file_fyear_fcoh4, plot = ewaa_fyear_fcohort4, width = 12, height = 8, units = "in", dpi = 300)

 #sim 5
 output_file_fyear_fcoh5 <- paste0(output_directory_fyear_fcoh, "ewaa_fyear_fcohort_sim5.png")
 ewaa_fyear_fcohort5<-ggplot(ewaa_m3, aes(x = factor(real_year), y = ewaa_fyear_fcoh4, color = factor(age),
                                          group = factor(age))) +
   geom_text(aes(label=round(age,2)), size = 4.5) +
   geom_line(alpha = 0.85) +ylim(0, 2.5)+
   theme_bw() +
   labs(x = "Year", y = "Weight") + annotate("text", x = 62, y = 2.5, parse = FALSE,label = "Sim 5",size=5)+
   theme(axis.title = element_text(size = 17),
         axis.text = element_text(size = 15),
         legend.title = element_text(size = 17),
         legend.text = element_text(size = 15),
         axis.text.x = element_text(angle = 90,hjust = 1,size=8),
         legend.position = "none",
         axis.title.x = element_blank(),
         axis.ticks=element_blank())
 # Save the plot to the specified directory
 ggsave(output_file_fyear_fcoh5, plot = ewaa_fyear_fcohort5, width = 12, height = 8, units = "in", dpi = 300)


 #sim 6
 output_file_fyear_fcoh6 <- paste0(output_directory_fyear_fcoh, "ewaa_fyear_fcohort_sim6.png")
 ewaa_fyear_fcohort6<-ggplot(ewaa_m3, aes(x = factor(real_year), y = ewaa_fyear_fcoh5, color = factor(age),
                                          group = factor(age))) +
   geom_text(aes(label=round(age,2)), size = 4.5) +
   geom_line(alpha = 0.85) +ylim(0, 2.5)+
   theme_bw() +
   labs(x = "Year", y = "Weight") + annotate("text", x = 62, y = 2.5, parse = FALSE,label = "Sim 6",size=5)+
   theme(axis.title = element_text(size = 17),
         axis.text = element_text(size = 15),
         legend.title = element_text(size = 17),
         legend.text = element_text(size = 15),
         axis.text.x = element_text(angle = 90,hjust = 1,size=8),
         legend.position = "none",
         axis.title.x = element_blank(),
         axis.ticks=element_blank())
 # Save the plot to the specified directory
 ggsave(output_file_fyear_fcoh6, plot = ewaa_fyear_fcohort6, width = 12, height = 8, units = "in", dpi = 300)


 #sim 7
 output_file_fyear_fcoh7 <- paste0(output_directory_fyear_fcoh, "ewaa_fyear_fcohort_sim7.png")
 ewaa_fyear_fcohort7<-ggplot(ewaa_m3, aes(x = factor(real_year), y = ewaa_fyear_fcoh6, color = factor(age),
                                          group = factor(age))) +
   geom_text(aes(label=round(age,2)), size = 4.5) +
   geom_line(alpha = 0.85) +ylim(0, 2.5)+
   theme_bw() +
   labs(x = "Year", y = "Weight") + annotate("text", x = 62, y = 2.5, parse = FALSE,label = "Sim 7",size=5)+
   theme(axis.title = element_text(size = 17),
         axis.text = element_text(size = 15),
         legend.title = element_text(size = 17),
         legend.text = element_text(size = 15),
         axis.text.x = element_text(angle = 90,hjust = 1,size=8),
         legend.position = "none",
         axis.title.x = element_blank(),
         axis.ticks=element_blank())
 # Save the plot to the specified directory
 ggsave(output_file_fyear_fcoh7, plot = ewaa_fyear_fcohort7, width = 12, height = 8, units = "in", dpi = 300)

 #sim 8
 output_file_fyear_fcoh8 <- paste0(output_directory_fyear_fcoh, "ewaa_fyear_fcohort_sim8.png")
 ewaa_fyear_fcohort8<-ggplot(ewaa_m3, aes(x = factor(real_year), y = ewaa_fyear_fcoh7, color = factor(age),
                                          group = factor(age))) +
   geom_text(aes(label=round(age,2)), size = 4.5) +
   geom_line(alpha = 0.85) +ylim(0, 2.5)+
   theme_bw() +
   labs(x = "Year", y = "Weight") + annotate("text", x = 62, y = 2.5, parse = FALSE,label = "Sim 8",size=5)+
   theme(axis.title = element_text(size = 17),
         axis.text = element_text(size = 15),
         legend.title = element_text(size = 17),
         legend.text = element_text(size = 15),
         axis.text.x = element_text(angle = 90,hjust = 1,size=8),
         legend.position = "none",
         axis.title.x = element_blank(),
         axis.ticks=element_blank())
 # Save the plot to the specified directory
 ggsave(output_file_fyear_fcoh8, plot = ewaa_fyear_fcohort8, width = 12, height = 8, units = "in", dpi = 300)


 #sim 9
 output_file_fyear_fcoh9 <- paste0(output_directory_fyear_fcoh, "ewaa_fyear_fcohort_sim9.png")
 ewaa_fyear_fcohort9<-ggplot(ewaa_m3, aes(x = factor(real_year), y = ewaa_fyear_fcoh8, color = factor(age),
                                          group = factor(age))) +
   geom_text(aes(label=round(age,2)), size = 4.5) +
   geom_line(alpha = 0.85) +ylim(0, 2.5)+
   theme_bw() +
   labs(x = "Year", y = "Weight") + annotate("text", x = 62, y = 2.5, parse = FALSE,label = "Sim 9",size=5)+
   theme(axis.title = element_text(size = 17),
         axis.text = element_text(size = 15),
         legend.title = element_text(size = 17),
         legend.text = element_text(size = 15),
         axis.text.x = element_text(angle = 90,hjust = 1,size=8),
         legend.position = "none",
         axis.title.x = element_blank(),
         axis.ticks=element_blank())
 # Save the plot to the specified directory
 ggsave(output_file_fyear_fcoh9, plot = ewaa_fyear_fcohort9, width = 12, height = 8, units = "in", dpi = 300)


 #Sim 10
 output_file_fyear_fcoh10 <- paste0(output_directory_fyear_fcoh, "ewaa_fyear_fcohort_sim10.png")
 ewaa_fyear_fcohort10<-ggplot(ewaa_m3, aes(x = factor(real_year), y = ewaa_fyear_fcoh9, color = factor(age),
                                          group = factor(age))) +
   geom_text(aes(label=round(age,2)), size = 4.5) +
   geom_line(alpha = 0.85) +ylim(0, 2.5)+
   theme_bw() +
   labs(x = "Year", y = "Weight") + annotate("text", x = 62, y = 2.5, parse = FALSE,label = "Sim 10",size=5)+
   theme(axis.title = element_text(size = 17),
         axis.text = element_text(size = 15),
         legend.title = element_text(size = 17),
         legend.text = element_text(size = 15),
         axis.text.x = element_text(angle = 90,hjust = 1,size=8),
         legend.position = "none",
         axis.title.x = element_blank(),
         axis.ticks=element_blank())
 # Save the plot to the specified directory
 ggsave(output_file_fyear_fcoh10, plot = ewaa_fyear_fcohort10, width = 12, height = 8, units = "in", dpi = 300)

 #Sim 11
 output_file_fyear_fcoh11 <- paste0(output_directory_fyear_fcoh, "ewaa_fyear_fcohort_sim11.png")
 ewaa_fyear_fcohort11<-ggplot(ewaa_m3, aes(x = factor(real_year), y = ewaa_fyear_fcoh10, color = factor(age),
                                           group = factor(age))) +
   geom_text(aes(label=round(age,2)), size = 4.5) +
   geom_line(alpha = 0.85) +ylim(0, 2.5)+
   theme_bw() +
   labs(x = "Year", y = "Weight") + annotate("text", x = 62, y = 2.5, parse = FALSE,label = "Sim 10",size=5)+
   theme(axis.title = element_text(size = 17),
         axis.text = element_text(size = 15),
         legend.title = element_text(size = 17),
         legend.text = element_text(size = 15),
         axis.text.x = element_text(angle = 90,hjust = 1,size=8),
         legend.position = "none",
         axis.title.x = element_blank(),
         axis.ticks=element_blank())
 # Save the plot to the specified directory
 ggsave(output_file_fyear_fcoh11, plot = ewaa_fyear_fcohort10, width = 12, height = 8, units = "in", dpi = 300)


 #Real Predict weight from GAM # This in log escale to be comparable need to be in exponential
 output_file_fyear_fcoh12 <- paste0(output_directory_fyear_fcoh, "ewaa_GAM_03.png")
 head(ewaa_m3,n=10)
 ewaa_fyear_fcohort12<-ggplot(ewaa_m3, aes(x = factor(real_year), y = exp(pred_weight), color = factor(age),
                                           group = factor(age))) +
   geom_text(aes(label=round(age,2)), size = 4.5) +
   geom_line(alpha = 0.85) +ylim(0, 2.5)+
   theme_bw() +
   labs(x = "Year", y = "Weight") + annotate("text", x = 62, y = 2.5, parse = FALSE,label = "",size=5)+
   theme(axis.title = element_text(size = 17),
         axis.text = element_text(size = 15),
         legend.title = element_text(size = 17),
         legend.text = element_text(size = 15),
         axis.text.x = element_text(angle = 90,hjust = 1,size=8),
         legend.position = "none",
         axis.title.x = element_blank(),
         axis.ticks=element_blank())

 # Save the plot to the specified directory
 ggsave(output_file_fyear_fcoh12, plot = ewaa_fyear_fcohort12, width = 12, height = 8, units = "in", dpi = 300)

 #########################################################################################




