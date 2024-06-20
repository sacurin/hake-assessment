# Title: Empirical weight at age pacific hake MSE
# Purpose: To visualize the effect of weight at age, and to evaluate the cohort and year effect in pacific hake
# Author: Sandra Curin-Osorio, Kelli Johnson, Kristin Marshall, Aaron Berger
# method: Fit a non spatial temporal generalized linear mixed effects model (GLMM) with sdmTMB
# Project: Hake MSE
# Date: 4/24/2024

library(gfdata)
library(gfplot)
library(ggmagnify)
#library(hakedataUSA)
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
library(dplyr)
library(ggplot2)

#devtools::install("C:/GitHub/hake-assessment", dependencies = TRUE)
rm(list = c("plot_map"))

# Load the WAA matrix (only use fishery data), the data is in the hake-assessment folder
setwd("C:\\GitHub\\hake-assessment\\doc")
devtools::load_all()

# This code fits three statistical models to examine the effects of age, cohort, and sex on estimation of weight at age of pacific hake using the sdmTMB package

# Model 1:  w~s(age)+(1|fcohort)+(1|fyear)+sex
#(1|fcohort) indicates that we're allowing intercepts to vary randomly across different levels of cohort.
#(1|fyear) indicates that we're allowing intercepts to vary randomly across different levels of year, accounting for potential correlation among observations within the same year.
# Sex as a fixed effect in the model (Kelli, J)

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

# Calculate the variance-covariance matrix of random effects from M1
random_vcov <- vcov.sdmTMB(m1)

print(m1)
## Extracting parameters as a data frame:

# Extracting fixed effect from M1
fixed_eff_m1<-tidy(m1,conf.int=TRUE)
print(fixed_eff_m1)

#Dispersion parameters for random effect from M1
random_par_m1 <- tidy(m1,"ran_pars",conf.int=TRUE)
print(random_par_m1)
# term    estimate std.error conf.low conf.high
# <chr>      <dbl>     <dbl>    <dbl>     <dbl>
#   1 phi       0.261   0.000375   0.260      0.262 # A dispersion parameter for a distribution
# 2 sigma_G   0.199   0.0186     0.166      0.239   # IID random intercept variance (fcohort)
# 3 sigma_G   0.0986  0.0103     0.0803     0.121   # IID random intercept variance (fyear)

#random effect values by variable
random_eff_val_m1<-tidy(m1,"ran_vals",conf.int=TRUE)
print(random_eff_val_m1,n=110)


# Extracting dispersion parameter
dispersion_param_m1 <- 0.26
print(dispersion_param_m1)

# Extracting ML criterion at convergence
ML_criterion_m1 <- -154961.173  #this value is coming from the output
print(ML_criterion_m1)

# grab the internal parameter list at estimated values:
#pars <- sdmTMB::get_pars(m1)

#Perform several sanity checks from the M1
sanity(m1)

predictions_m1 <- predict(m1)
head(predictions_m1)
predictions_m1$resids <- residuals(m1) # randomized quantile residuals

#Figure 1 in phase 1 word document
#plots to visualize the residuals and predictions from the M1
ggplot(predictions_m1, aes(age, weight, col = resids)) + scale_colour_gradient2() +
  geom_point() + facet_wrap(~year)
hist(predictions_m1$resids)

#MODEL DIAGNOSTIC
residuals_m1 <- residuals(m1)
residuals_m1 <- residuals_m1[is.finite(residuals_m1)]
png(fs::path("C:/GitHub/hake-assessment/Tables_figures_Sandra/weight-at-age-qq_m1.png"))
qqnorm(residuals_m1); qqline(residuals_m1)
dev.off()
saveRDS(
  m1,
  file = fs::path("C:/GitHub/hake-assessment/Tables_figures_Sandra/weight-at-age-model_m1.rds")
)


# create prediction grid
pred_grid <- expand.grid(
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

# Get estimates
preds <- predict(m1, newdata = pred_grid)
preds <- preds |>
  dplyr::mutate(est_weight = exp(est))
saveRDS(preds,
  file = fs::path("C:/GitHub/hake-assessment/Tables_figures_Sandra/weight-at-age-preds_m1.rds")
)

##Figure 2 in phase 1 word document
# adjusted predictions for all 4 interaction terms
pr_m1<-ggplot(preds, aes(age, est_weight, colour = sex)) + geom_point() +
  facet_wrap(~sex)+labs(y="Predict weight")

#head(preds)
#ewa_m1<-ggplot(preds, aes(age, est_weight,colour=sex)) + geom_point()+labs(x = "Age", y = "Estimate weight")

# make preds# make EWAA
ewaa_m1 <- preds |>
  dplyr::group_by(year, age) |>
  dplyr::summarise(
    n = dplyr::n(),
    pred_weight = mean(exp(est)),
    upper=quantile(exp(est),probs = 0.975),
    lower = quantile(exp(est), probs = 0.025),
  ) |>
  as.data.frame() |>
  dplyr::select(year, age, pred_weight,upper,lower)

ewaa_complete <- pad_weight_at_age(ewaa_m1)
hakedataUSA:::write_wtatage_file(
  file = fs::path(bridge_dir, "70-tv-weight-at-age", "wtatage.ss"),
  data = ewaa_complete,
  maturity = hakedataUSA::maturity_at_age
)

#compute anomaly relative to the mean
# WAA_re_m1_all <- ewaa_m1 %>%
#   left_join(ewaa_m1, by = c("ages", "model")) %>%
#   mutate(anom = (exp(pred_weight) - mean(exp(est)) / mean(exp(est))))
# WAA_re_m1_all <- ewaa_m1 %>%mutate(anom = (exp(pred_weight) - pred_weight / pred_weight))%>%lwr_95 = exp(pred_weight - (1.96 * sd))


# Ggplot visualization of the predicted weight at age over the years, using M2: weight ~ 1 + s(age)+ (1|fyear) + cohort
# w~s(age)+(1|fcohort)+(1|fyear)+sex Plot
weight_OP1<-ggplot(ewaa_m1, aes(x = factor(year), y = pred_weight, color = factor(age),
           group = factor(age))) +
  geom_text(aes(label=round(age,2)), size = 4.5) +
  geom_line(alpha = 0.85) +ylim(0, 2.5)+
  theme_bw() +
  labs(x = "Year", y = "Predicted Weight") +
  geom_ribbon(aes(ymin = lower, ymax = upper,fill=factor(age)), alpha = 0.1)+
  theme(axis.title = element_text(size = 17),
        axis.text = element_text(size = 15),
        legend.title = element_text(size = 17),
        legend.text = element_text(size = 15),
        legend.position = "none",
        axis.title.x = element_blank(),
        axis.text.x=element_blank(),
        axis.ticks=element_blank())+ annotate("text", x = 40, y = 2.4, parse = FALSE,label = "weight ~ 1 + s(age) + (1|fcohort) +(1|fyear) + sex",size=5)

#function to help maps.

plot_map <- function(dat, column) {
  ggplot(dat, aes(age, weight, fill = {{ column }})) +
    geom_raster() +
    coord_fixed()
}

predict_plot<-plot_map(predictions_m1, exp(est)) +
  scale_fill_viridis_c(
    trans = "sqrt",
    # trim extreme high values to make spatial variation more visible
    na.value = "yellow", limits = c(0, quantile(exp(predictions_m1$est), 0.995))
  ) +
  facet_wrap(~year) +
  ggtitle("Prediction (fixed effects + all random effects)",
          subtitle = paste("maximum estimated biomass density =", round(max(exp(predictions_m1$est))))
  )

Predic_fix_effect<-plot_map(predictions_m1, exp(est)) +
  scale_fill_viridis_c(trans = "sqrt") +
  ggtitle("Prediction (fixed effects only)")

# summary of the Akaike Information Criterion (AIC) for the fitted model m1
AIC_m1<-AIC(m1)

##########################################################################
# Model 2- weight ~ s(age) + (1|fcohort) + sex
m2 <- sdmTMB::sdmTMB(
  data = weight_at_age_df %>%
    dplyr::filter(sex != "U") %>%
    dplyr::mutate(
      age = ifelse(age > 15, 15, age),
      cohort = year - age,
      fyear = as.factor(year),
      fcohort = as.factor(cohort)
    ),
  formula = weight ~ s(age) + (1|fcohort) + sex,
  family = lognormal(link = "log"),
  spatial = "off",
  spatiotemporal = "off",
  time = "year",
  control = sdmTMB::sdmTMBcontrol(newton_loops = 1)
)
print(m2)

## Extracting parameters as a data frame:

# Extracting fixed effect of M2
fixed_eff_m2<-tidy(m2,conf.int=TRUE)
print(fixed_eff_m2)
# term        estimate std.error conf.low conf.high
# <chr>          <dbl>     <dbl>    <dbl>     <dbl>
#   1 (Intercept)  -0.565    0.0336   -0.631    -0.499
# 2 sexM         -0.0948   0.00111  -0.0970   -0.0927

#Dispersion parameters for random effect of M2
random_par_m2 <- tidy(m2,"ran_pars",conf.int=TRUE)
print(random_par_m2)

# term    estimate std.error conf.low conf.high
# <chr>      <dbl>     <dbl>    <dbl>     <dbl>
# 1 phi        0.270  0.000388    0.270     0.271 # A dispersion parameter for a distribution
# 2 sigma_G    0.264  0.0241      0.221     0.316 # IID random intercept variance (fcohort)

#random effect values by variable
random_eff_val_m2<-tidy(m2,"ran_vals",conf.int=TRUE)
print(random_eff_val_m2,n=110)


# Extracting dispersion parameter of M2
dispersion_param_m2 <- 0.27
print(dispersion_param_m2)

# Extracting ML criterion at convergence of M2
ML_criterion_m2 <- -146632.857  #this value is coming from the output
print(ML_criterion_m2)

# grab the internal parameter list at estimated values:
#pars <- sdmTMB::get_pars(m1)


# Perform model assessment M2
sanity(m2)
predictions_m2 <- predict(m2)
head(predictions_m2)
predictions_m2$resids <- residuals(m2) # randomized quantile residuals

#Figure 3 in Phase 1 document (Residuals of weight at age by year for Pacific hake model 2.  )
ggplot(predictions_m2, aes(age, weight, col = resids)) + scale_colour_gradient2() +
  geom_point() + facet_wrap(~year)
hist(predictions_m1$resids)
residuals_m2 <- residuals(m2)
residuals_m2 <- residuals_m2[is.finite(residuals_m2)]
png(fs::path("", "srv", "hake", "other", "tv", "weight-at-age-qq_m2.png"))
qqnorm(residuals_m2); qqline(residuals_m2)
dev.off()

# Save the model
saveRDS(
  m2,
  file = fs::path("", "srv", "hake", "other", "tv", "weight-at-age-model.rds")
)

# Create prediction grid of M2
pred_grid_m2 <- expand.grid(
  year = unique(m2[["data"]][["year"]]),
  sex = unique(m2[["data"]][["sex"]]),
  age = 0:15
) %>%
  dplyr::mutate(cohort = year - age) %>%
  dplyr::filter(cohort %in% unique(m2[["data"]][["cohort"]])) %>%
  dplyr::mutate(
    fyear = as.factor(year),
    fcohort = as.factor(cohort)
  )

# Get estimates M2
preds_m2 <- predict(m2, newdata = pred_grid_m2)
preds_m2 <- preds_m2 %>%
  dplyr::mutate(est_weight = exp(est))
saveRDS(
  preds,
  file = fs::path("C:/GitHub/hake-assessment/Tables_figures_Sandra/weight-at-age-preds_m2.rds")
)
#visualize predictions
head(preds_m2)

# Make EWAA
ewaa_m2 <- preds_m2 %>%
  dplyr::group_by(year, age) %>%
  dplyr::summarise(
    n = dplyr::n(),
    pred_weight = mean(exp(est)),
    upper=quantile(exp(est),probs = 0.975),
    lower = quantile(exp(est), probs = 0.025),
  ) %>%
  as.data.frame() %>%
  dplyr::select(year, age, pred_weight,upper,lower)

ewaa_complete_m2 <- pad_weight_at_age(ewaa_m2)
hakedataUSA:::write_wtatage_file(
  file = fs::path(bridge_dir, "70-tv-weight-at-age", "wtatage.ss"),
  data = ewaa_complete,
  maturity = hakedataUSA::maturity_at_age
)

# Summary and visualization
summary(m2)
summary(ewaa_m2)

# Ggplot visualization of the predicted weight at age over the years, using M2: weight ~ 1 + s(age)+ (1|fyear) + cohort
weight_OP2<-ggplot(ewaa_m2, aes(x = factor(year), y = pred_weight, color = factor(age),
                                group = factor(age))) +
  geom_text(aes(label=round(age,2)), size = 4.5) +
  geom_line(alpha = 0.85) +
  theme_bw() + labs(x = "Year", y = "Predicted Weight") +ylim(0, 2.5)+geom_ribbon(aes(ymin = lower, ymax = upper,fill=factor(age)), alpha = 0.1)+
  theme(axis.title = element_text(size = 17),
        axis.text = element_text(size = 15),
        legend.title = element_text(size = 17),
        legend.text = element_text(size = 15),
        legend.position = "none",
        axis.title.x = element_blank(),
        axis.text.x=element_blank(),
        axis.ticks=element_blank())+ annotate("text", x = 42.5, y = 2.4, parse = FALSE,label = "weight ~ s(age) + (1|fcohort)  + sex",size=5)

# summary of the Akaike Information Criterion (AIC) for the fitted model m2
summary(AIC(m2))
AIC_m2<-AIC(m2)

##########################################################################################################################
# Model 3.	weight ~ s(age)+(1|fyear)+sex

m3 <- sdmTMB::sdmTMB(
  data = weight_at_age_df %>%
    dplyr::filter(sex != "U") %>%
    dplyr::mutate(
      age = ifelse(age > 15, 15, age),
      cohort = year - age,
      fyear = as.factor(year),
      fcohort = as.factor(cohort)
    ),
  formula = weight ~ s(age) + (1|fyear) + sex,
  family = lognormal(link = "log"),
  spatial = "off",
  spatiotemporal = "off",
  time = "year",
  control = sdmTMB::sdmTMBcontrol(newton_loops = 1)
)
print(m3)

# Extracting parameters as a data frame:

# Extracting fixed effect
fixed_eff_m3<-tidy(m3,conf.int=TRUE)
print(fixed_eff_m3)
# term        estimate std.error conf.low conf.high
# <chr>          <dbl>     <dbl>    <dbl>     <dbl>
#   1 (Intercept)  -0.565    0.0336   -0.631    -0.499
# 2 sexM         -0.0948   0.00111  -0.0970   -0.0927

#Dispersion parameters for random effect
random_par_m3 <- tidy(m3,"ran_pars",conf.int=TRUE)
print(random_par_m3)

# term    estimate std.error conf.low conf.high
# <chr>      <dbl>     <dbl>    <dbl>     <dbl>
#   1 phi        0.278  0.000399    0.277     0.279 # A dispersion parameter for a distribution
# 2 sigma_G    0.131  0.0132      0.107     0.159   # IID random intercept variance (fyear)

#random effect values by variable
random_eff_val_m3<-tidy(m3,"ran_vals",conf.int=TRUE)
print(random_eff_val_m3,n=110)


# Extracting dispersion parameter
dispersion_param_m3 <- 0.28
print(dispersion_param_m3)

# Extracting ML criterion at convergence
ML_criterion_m3 <- -139889.807  #this value is coming from the output
print(ML_criterion_m3)

# grab the internal parameter list at estimated values:
#pars <- sdmTMB::get_pars(m1)


# Perform model assessment
sanity(m3)
predictions_m3 <- predict(m3)
head(predictions_m3)
predictions_m3$resids <- residuals(m3) # randomized quantile residuals
#Figure 3 In the Phase 1 document.
ggplot(predictions_m3, aes(age, weight, col = resids)) + scale_colour_gradient2() +
  geom_point() + facet_wrap(~year)
hist(predictions_m1$resids)
residuals_m3 <- residuals(m3)
residuals_m3 <- residuals_m3[is.finite(residuals_m3)]
png(fs::path("", "srv", "hake", "other", "tv", "weight-at-age-qq_m3.png"))
qqnorm(residuals_m3); qqline(residuals_m3)
dev.off()
# Save the model
saveRDS(
  m3,
  file = fs::path("", "srv", "hake", "other", "tv", "weight-at-age-model_m3.rds")
)

# Create prediction grid for M3
pred_grid_m3 <- expand.grid(
  year = unique(m3[["data"]][["year"]]),
  sex = unique(m3[["data"]][["sex"]]),
  age = 0:15
) %>%
  dplyr::mutate(cohort = year - age) %>%
  dplyr::filter(cohort %in% unique(m3[["data"]][["cohort"]])) %>%
  dplyr::mutate(
    fyear = as.factor(year),
    fcohort = as.factor(cohort)
  )

# Get estimates fro M3
preds_m3 <- predict(m3, newdata = pred_grid_m3)
preds_m3 <- preds_m3 %>%
  dplyr::mutate(est_weight = exp(est))
saveRDS(
  preds,
  file = fs::path("C:/GitHub/hake-assessment/Tables_figures_Sandra/weight-at-age-preds_m3.rds")
)

head(preds_m3)

# Make EWAA
ewaa_m3 <- preds_m3 %>%
  dplyr::group_by(year, age) %>%
  dplyr::summarise(
    n = dplyr::n(),
    pred_weight = mean(exp(est)),
    upper=quantile(exp(est),probs = 0.975),
    lower = quantile(exp(est), probs = 0.025),
  ) %>%
  as.data.frame() %>%
  dplyr::select(year, age, pred_weight,upper,lower)

ewaa_complete_m3 <- pad_weight_at_age(ewaa_m3)
hakedataUSA:::write_wtatage_file(
  file = fs::path(bridge_dir, "70-tv-weight-at-age", "wtatage.ss"),
  data = ewaa_complete_m3,
  maturity = hakedataUSA::maturity_at_age
)

head(ewaa_m3)

# Summary and visualization
# # Ggplot visualization of the predicted weight at age over the years, using M3: weight ~ 1 + s(age)+ (1|fyear) + cohort
weight_OP3<-ggplot(ewaa_m3, aes(x = factor(year), y = pred_weight, color = factor(age),
                             group = factor(age))) +
  geom_text(aes(label=round(age,2)), size = 4.5) +
  geom_line(alpha = 0.85) + ylim(0, 2.5)+geom_ribbon(aes(ymin = lower, ymax = upper,fill=factor(age)), alpha = 0.1)+
  theme_bw() + labs(x = "Year", y = "Predicted Weight") +
  theme(axis.title = element_text(size = 17),
        axis.text = element_text(size = 15),
        legend.title = element_text(size = 17),
        legend.text = element_text(size = 15),
        axis.text.x = element_text(angle = 90,hjust = 1,size=8),
        legend.position = "none",
        axis.title.x=element_blank(),
        axis.ticks=element_blank())+ annotate("text", x = 43, y = 2.4, parse = FALSE,label = "weight ~ s(age)+ (1|fyear) + sex",size=5)

# summary of the Akaike Information Criterion (AIC) for the fitted model M3
summary(AIC(m3))
AIC_m3<-AIC(m3)

#######################################################################################################################################
# median of Predict weight by model- Summary plot
GAM_all_plot<-ggarrange(weight_OP1+ rremove("x.text"), weight_OP2+ rremove("x.text"), weight_OP3 ,
          ncol = 1, nrow = 3)

#######################################################################################################################################
#Visualization by model
# random parameters of each GAM

# data.frame
random_par_dat_m1<-data.frame(random_par_m1)
random_par_dat_m2<-data.frame(random_par_m2)
random_par_dat_m3<-data.frame(random_par_m3)

# Add model column
random_par_dat_m1$model <- "Age+(1|fcohort)+(1|fyear)"
random_par_dat_m2$model <- "Age+(1|fcohort)"
random_par_dat_m3$model <- "Age+(1|fyear)"

# Add model AIC
random_par_dat_m1$AIC <- AIC_m1
random_par_dat_m2$AIC <- AIC_m2
random_par_dat_m3$AIC <- AIC_m3

# Combine dataframes
combined_df <- rbind(random_par_dat_m1, random_par_dat_m2, random_par_dat_m3)

# Replace "sigma_G" with a symbol in Combined_df
combined_df <- combined_df %>%
  mutate(term = ifelse(term == "sigma_G", "\u03C3_G", term))

# Print the updated dataframe
print(combined_df)

par_plot <- ggplot(combined_df, aes(x = model, y = estimate,
                                   ymin = conf.low, ymax = conf.high,
                                   shape = factor(estimate),
                                   color = model)) +
  geom_pointrange(size = 1.1) +
  ggsci::scale_color_jco() +
  facet_wrap(~term, scales = "free_y", labeller = label_parsed, ncol = 1) +
  guides(color="none") +
  theme_bw() +
  theme(legend.position = 'none',
        axis.title = element_text(size = 17),
        axis.text = element_text(size = 15, color = "black"),
        strip.text = element_text(size = 17)) +
  labs(x = "Model", y = "Parameter Estimate")


# Visualize AIC across models
aic_plot <- ggplot(combined_df, aes(x = model, y = AIC, group = 1)) +
  geom_line(lty = 2, size = 1.3) +
  geom_point(aes(color = model),size = 5) +
  ggsci::scale_color_jco() +
  theme_bw() +
  theme(legend.position = "none",
        axis.title = element_text(size = 17),
        axis.text = element_text(size = 15, color = "black"),
        strip.text = element_text(size = 17)) + labs(x = "Model", y = "AIC")

plot_grid(par_plot, aic_plot, rel_widths = c(0.5, 0.5),
          align = "hv", axis = "bl", ncol = 2,
          labels = c("A", "B"), label_size = 23, hjust = -0.5)


####################################################################################################################
#Comparing weight estimated for the three model (M1-Age+fcohort, M2-Age+fcohort+fyear, M3- Age+fyear).
# Combine dataframes with 95% normal CIs

# Add model column
ewaa_m1$model <- "Age+(1|fcohort)+(1|fyear)"
ewaa_m2$model <- "Age+(1|fcohort)"
ewaa_m3$model <- "Age+(1|fyear)"

combined_ewaa_df <- rbind(ewaa_m1, ewaa_m2, ewaa_m3)
head(combined_ewaa_df)

ggplot(combined_ewaa_df,
       aes(x = year, y = pred_weight, color = factor(model), group = factor(model),
           fill = model)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.35, color = NA) +
  geom_line(alpha = 1, size = 1) +
  ggsci::scale_color_nejm( ) +
  ggsci::scale_fill_nejm( ) +
  theme_bw() +
  facet_wrap(~age, ncol = 5, scales = "free")+
  guides(color=guide_legend(ncol=3)) +
  labs(x = "Year", y = "Weight", color = "Model", fill = "Model") +
  theme(axis.title = element_text(size = 17),
        axis.text = element_text(size = 10, color = "black"),
        legend.title = element_text(size = 21),
        legend.text = element_text(size = 17),
        strip.text = element_text(size = 17),
        legend.position = "top",
        legend.background = element_blank(),
        legend.key.width = unit(0.75, "cm"))

#######################################################################################################
# OTHER GAM EFFECTS
#############################################################################
# 4.	weight ~  s(age)+ (1|fcohort)+ fyear + sex

# Fit the model with the new formula
m4 <- sdmTMB::sdmTMB(
  data = weight_at_age_df %>%
    dplyr::filter(sex != "U") %>%
    dplyr::mutate(
      age = ifelse(age > 15, 15, age),
      cohort = year - age,
      fyear = as.factor(year),
      fcohort = as.factor(cohort)
    ),
  formula = weight ~ s(age) + (1|fcohort) + fyear + sex,
  family = lognormal(link = "log"),
  spatial = "off",
  spatiotemporal = "off",
  time = "year",
  control = sdmTMB::sdmTMBcontrol(newton_loops = 1)
)
print(m4)

#Extract coefficients
tidy(m4)

#Confidence intervals on the fixed effects
tidy(m4,conf.int=TRUE)

#Confidence intervals on the random effects
tidy(m4,"ran_pars",conf.int=TRUE)

tidy(m4,"ran_vals")
print(tidy(m4,"ran_vals"),n=110)

# Perform model assessment
sanity(m4)
residuals_m4 <- residuals(m4)
residuals_m4 <- residuals_m4[is.finite(residuals_m4)]
png(fs::path("", "srv", "hake", "other", "tv", "weight-at-age-qq.png"))
qqnorm(residuals_m4); qqline(residuals_m4)
dev.off()

# Save the model
saveRDS(
  m4,
  file = fs::path("", "srv", "hake", "other", "tv", "weight-at-age-model_m4.rds")
)

# Create prediction grid
pred_grid_m4 <- expand.grid(
  year = unique(m4[["data"]][["year"]]),
  sex = unique(m4[["data"]][["sex"]]),
  age = 0:15
) %>%
  dplyr::mutate(cohort = year - age) %>%
  dplyr::filter(cohort %in% unique(m4[["data"]][["cohort"]])) %>%
  dplyr::mutate(
    fyear = as.factor(year),
    fcohort = as.factor(cohort)
  )

# Get estimates
preds_m4 <- predict(m4, newdata = pred_grid_m4)
preds_m4 <- preds_m4 %>%
  dplyr::mutate(est_weight = exp(est))
saveRDS(
  preds,
  file = fs::path("/srv/hake/other/tv", "weight-at-age-preds_m4.rds")
)

head(preds_m4)
ewaa_m4<-ggplot(preds_m4, aes(age, est_weight,colour=sex)) + geom_point()+labs(x = "Age", y = "Estimate weight")

# Make EWAA
ewaa_m4 <- preds_m4 %>%
  dplyr::group_by(year, age) %>%
  dplyr::summarise(
    n = dplyr::n(),
    pred_weight = mean(exp(est))
  ) %>%
  as.data.frame() %>%
  dplyr::select(year, age, pred_weight)

ewaa_complete_m4 <- pad_weight_at_age(ewaa_m4)
hakedataUSA:::write_wtatage_file(
  file = fs::path(bridge_dir, "70-tv-weight-at-age", "wtatage.ss"),
  data = ewaa_complete_m4,
  maturity = hakedataUSA::maturity_at_age
)

# Summary and visualization
summary(m4)
summary(ewaa_m4)
weight_OP4<-ggplot(ewaa_m4, aes(x = factor(year), y = pred_weight, color = factor(age),
                                group = factor(age))) +
  geom_text(aes(label=round(age,2)), size = 4.5) +
  geom_line(alpha = 0.85) + ylim(0, 2.5)+
  theme_bw() +
  labs(x = "Year", y = "Weight") +
  theme(axis.title = element_blank(),
        axis.text = element_text(size = 15),
        legend.title = element_text(size = 17),
        legend.text = element_text(size = 15),
        legend.position = "none",
        axis.text.x=element_blank(),
        axis.ticks=element_blank())+ annotate("text", x = 29, y = 2.5, parse = FALSE,label = "weight ~ s(age) + (1|fcohort) + fyear + sex ",size=5)


weight_OP2<-ggplot(ewaa_m2, aes(x = factor(year), y = pred_weight, color = factor(age),
                                group = factor(age))) +
  geom_text(aes(label=round(age,2)), size = 4.5) +
  geom_line(alpha = 0.85) +
  theme_bw() + labs(x = "Year", y = "Weight") +ylim(0, 2.5)+
  theme(axis.title = element_blank(),
        axis.text = element_text(size = 15),
        legend.title = element_text(size = 17),
        legend.text = element_text(size = 15),
        legend.position = "none",
        axis.text.x=element_blank(),
        axis.ticks=element_blank())+ annotate("text", x = 29, y = 2.5, parse = FALSE,label = "weight ~ s(age) + fcohort + (1|fyear) + sex",size=5)

summary(AIC(m4))

############################################################################
# 5-	weight ~ s(age)+ (1|fcohort) + fyear

# Fit the model with the new formula
m5 <- sdmTMB::sdmTMB(
  data = weight_at_age_df %>%
    dplyr::filter(sex != "U") %>%
    dplyr::mutate(
      age = ifelse(age > 15, 15, age),
      cohort = year - age,
      fyear = as.factor(year),
      fcohort = as.factor(cohort)
    ),
  formula = weight ~ s(age)+ (1|fcohort) + (fyear),
  family = lognormal(link = "log"),
  spatial = "off",
  spatiotemporal = "off",
  time = "year",
  control = sdmTMB::sdmTMBcontrol(newton_loops = 1)
)
print(m5)

#Extract coefficients
tidy(m5)

#Confidence intervals on the fixed effects
tidy(m5,conf.int=TRUE)

#Confidence intervals on the random effects
tidy(m5,"ran_pars",conf.int=TRUE)

tidy(m5,"ran_vals")
print(tidy(m5,"ran_vals"),n=110)

# Perform model assessment
sanity(m5)
residuals_m5 <- residuals(m5)
residuals_m5 <- residuals_m5[is.finite(residuals_m5)]
png(fs::path("", "srv", "hake", "other", "tv", "weight-at-age-qq.png"))
qqnorm(residuals_m5); qqline(residuals_m5)
dev.off()

# Save the model
saveRDS(
  m5,
  file = fs::path("", "srv", "hake", "other", "tv", "weight-at-age-model_m5.rds")
)

# Create prediction grid
pred_grid_m5 <- expand.grid(
  year = unique(m5[["data"]][["year"]]),
  age = 0:15
) %>%
  dplyr::mutate(cohort = year - age) %>%
  dplyr::filter(cohort %in% unique(m5[["data"]][["cohort"]])) %>%
  dplyr::mutate(
    fyear = as.factor(year),
    fcohort = as.factor(cohort)
  )

# Get estimates
preds_m5 <- predict(m5, newdata = pred_grid_m5)
preds_m5 <- preds_m5 %>%
  dplyr::mutate(est_weight = exp(est))
saveRDS(
  preds,
  file = fs::path("/srv/hake/other/tv", "weight-at-age-preds_m5.rds")
)

head(preds_m5)
ewaa_m5<-ggplot(preds_m5, aes(age, est_weight)) + geom_point()+labs(x = "Age", y = "Estimate weight")


# Make EWAA
ewaa_m5 <- preds_m5 %>%
  dplyr::group_by(year, age) %>%
  dplyr::summarise(
    n = dplyr::n(),
    pred_weight = mean(exp(est))
  ) %>%
  as.data.frame() %>%
  dplyr::select(year, age, pred_weight)

ewaa_complete_m5 <- pad_weight_at_age(ewaa_m5)
hakedataUSA:::write_wtatage_file(
  file = fs::path(bridge_dir, "70-tv-weight-at-age", "wtatage.ss"),
  data = ewaa_complete,
  maturity = hakedataUSA::maturity_at_age
)

# Summary and visualization
summary(m5)
summary(ewaa_m5)

weight_OP5<-ggplot(ewaa_m5, aes(x = factor(year), y = pred_weight, color = factor(age),
                                group = factor(age))) +
  geom_text(aes(label=round(age,2)), size = 4.5) +
  geom_line(alpha = 0.85) + ylim(0, 2.5)+
  theme_bw() +
  labs(x = "Year", y = "Weight") +
  theme(axis.title = element_text(size = 17),
        axis.text = element_text(size = 15),
        axis.text.x = element_text(angle = 90,hjust = 1,size=8),
        legend.title = element_text(size = 17),
        legend.text = element_text(size = 15),
        legend.position = "none")+ annotate("text", x = 32, y = 2.5, parse = FALSE,label = "weight ~ s(age) + (1|fcohort)+ fyear",size=5)

summary(AIC(m5))

windows()
jpeg("ewaa_scenarios.jpg", width = 350, height = "350")
ggarrange(weight_OP1+ rremove("x.text"), weight_OP2+ rremove("x.text"), weight_OP3 + rremove("x.text"),weight_OP4 + rremove("x.text"),weight_OP5 ,
          ncol = 2, nrow = 3)
dev.off()

##########################################################################################################################################
