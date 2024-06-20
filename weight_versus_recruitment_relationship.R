# Title: Empirical weight versus recruitment.
# Purpose 1: plotting weight estimated by fcohort and fyear vs year, only for discussion.
# Purpose 2: lineal regreation between random effect variable and recruitment deviations, Lorenzo request.
# Author: Sandra Curin-Osorio,
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
library(tidyr)
#devtools::install_github("pacific-hake/hake-assessment") Didn't work
#devtools::install("C:/GitHub/hake-assessment", dependencies = TRUE)

rm(list = c("plot_map"))
# Load the WAA matrix (only use fishery data)/where your hake-assessment/doc folder is located
setwd("C:\\GitHub\\hake-assessment\\doc")
devtools::load_all()

# The code used a sdmTMB to replicates the empirical weight at age used on the hake assessment. In the stock assessment
# the empirical weight at age is modelated using GAM. The code predict the weight given the age
# 1- w~s(age)+(1|fcohort)+(1|fyear)+sex, the same GAM of the stock assessment to be consistent by Kelli Johnson
#(1|fcohort) indicates that we're allowing intercepts to vary randomly across different levels of cohort.
#(1|fyear) indicates that we're allowing intercepts to vary randomly across different levels of year, accounting for potential correlation among observations within the same year.
# Sex as a fixed effect in the model

# exploring random effect of fcohort and fyear

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

#Print the model fit
print(m1)

## Extracting parameters as a data frame:
# Extracting fixed effect
fixed_eff_m1<-tidy(m1,conf.int=TRUE)
print(fixed_eff_m1)

#Dispersion parameters for random effect
random_par_m1 <- tidy(m1,"ran_pars",conf.int=TRUE)
print(random_par_m1)
# term    estimate std.error conf.low conf.high
# <chr>      <dbl>     <dbl>    <dbl>     <dbl>
#   1 phi       0.261   0.000375   0.260      0.262 # A dispersion parameter for a distribution
# 2 sigma_G   0.199   0.0186     0.166      0.239   # IID random intercept variance (fcohort)
# 3 sigma_G   0.0986  0.0103     0.0803     0.121   # IID random intercept variance (fyear)

#random effect values by variable
random_eff_val_m1<-tidy(m1,"ran_vals",conf.int=TRUE)
print(random_eff_val_m1,n=115)

# Separate the 'term' column into 'year' and 'term'
random_eff_val_m1 <- separate(random_eff_val_m1, term, into = c("term", "year"), sep = "_", extra = "merge")

# Gather the data into key-value pairs
random_eff_val_m1_long <- pivot_longer(random_eff_val_m1, cols = c("estimate", "std.error", "conf.low", "conf.high"), names_to = "metric", values_to = "value")

# Spread the 'metric' column into separate columns
random_eff_val_m1_final <- pivot_wider(random_eff_val_m1_long, names_from = "metric", values_from = "value")

# View the final table
print(random_eff_val_m1_final)

# Filter data for fcohort
fcohort_data <- random_eff_val_m1_final[grep("^fcohort", random_eff_val_m1_final$term), ]

# Plot for fcohort with expanded x-axis range
output_directory <- "C:/GitHub/hake-assessment/Tables_figures_Sandra/"
fcohort_val <- paste0(output_directory, "ewaa_fcohort_random_eff_val.png")
plot_fcohort <- ggplot(fcohort_data, aes(x = as.integer(year), y = estimate)) +
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.3, fill = "red") +
  geom_hline(yintercept = 0, linetype = "solid", color = "grey") +  # Add a line at y = 0
  labs(x = "Year", y = "Weight") +
  ggtitle("Weight estimate for (1|fcohort)") +
  scale_x_continuous(breaks = seq(1960, 2022, by = 3)) +  # Set breaks at every 5 years
  theme(axis.title.x = element_text(size = 14, face = "bold"),  # Customize x-axis label
        axis.title.y = element_text(size = 14, face = "bold"))  # Customize y-axis label
ggsave(fcohort_val, plot = plot_fcohort, width = 12, height = 8, units = "in", dpi = 300)


# Filter data for fyear
fyear_data <- random_eff_val_m1_final[grep("^fyear", random_eff_val_m1_final$term), ]

# Plot for fyear with expanded x-axis range
output_directory <- "C:/GitHub/hake-assessment/Tables_figures_Sandra/"
fyear_val <- paste0(output_directory, "ewaa_fyear_random_eff_val.png")
plot_fyear <- ggplot(fyear_data, aes(x = as.integer(year), y = estimate)) +
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.3, fill = "green") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +  # Add a line at y = 0
  labs(x = "Year", y = "Weight") +
  ggtitle("Weight estimate for (1|fyear)") +
  scale_x_continuous(breaks = seq(1960, 2023, by = 3)) +  # Set breaks at every 5 years
  theme(axis.title.x = element_text(size = 14, face = "bold"),  # Customize x-axis label
        axis.title.y = element_text(size = 14, face = "bold"))  # Customize y-axis label
ggsave(fyear_val, plot = plot_fyear, width = 12, height = 8, units = "in", dpi = 300)


##########################################################################################
#Relationships between random effect levels and recruitment devs
##########################################################################################


#Recruitment devs from 1970-2022 (From stock assessment) & fcohort
rec_dev<-c(1.90617,-0.238677,-0.758393,1.47529,-1.05657,	0.30881,-1.5827,-1.64343,-1.98697,0.123949,	2.68975,-1.37203,-1.34778,-0.770718,	2.49053,	-2.01555,	-1.67691,	1.74704,	0.697539,
           -1.90878,	1.35377,	0.192922,	-1.8685,	1.11944,	1.15696,	0.275784,	0.59495,	-0.0304116,	0.654935,	2.5325,	-0.953306,	0.1495,	-3.10058,	0.454618,	-3.09334,	0.977308,	0.729547,
           -3.53521,	1.71368,	0.318232,	2.70529,	-0.847241,	0.349782,	-0.982123,	1.9944,-3.24806,	1.62556	,0.274458,	-0.947068,	-1.31558,	1.38513	,2.08257,	0.700675) # Main_RecrDev is coming from stock assessment
fcohort<-c(0.0674,	0.0186,	-0.0911,	-0.0312,	-0.174,	-0.114,	-0.379,	-0.0108,	-0.351,	-0.231,	-0.113,	-0.399,	-0.411,	-0.396,	-0.0731,	-0.27,	-0.159,	-0.0553,	-0.0407,
           -0.0621,	-0.00717,	0.078,	-0.00856,	0.0644,	0.141,	0.0412,	0.014,	0.128,	0.047,	0.0227,	0.0716,	0.116,	0.123,	0.164,	0.105,	-0.0186,	-0.146,	-0.182,	-0.202,
           -0.139,	-0.171,	-0.0639,	-0.00992,	-0.0491,	-0.0481,	0.0837,	0.105,	0.186,	-0.0706,	0.0997,	0.193,	-0.0385,	-0.143) # is coming from random effect of the GAM output.

recr<-c(5.12E+06,	600537,	360687,	3.42E+06,	271253,	1.08E+06,	163661,	4.10E+06,	107392,	892200,	1.15E+07,	198145,	204249,	369849,	9.67E+06,	106920,	150098,
        4.62E+06,	1.62E+06,	118726,	3.10E+06,	966033,	122300,	2.39E+06,	2.47E+06,	1.02E+06,	1.39E+06,	745126,	1.46E+06,	9.34E+06,	289390,	901703,	35902.8,
        1.25E+06,	35695,	2.06E+06,	1.57E+06,	21385.3,	4.08E+06,	983435,	1.10E+07,	313297,	1.06E+06,	295406,	5.86E+06,	30405.4,	3.93E+06,	1.04E+06,	306810,
        209274,	3.10E+06,	6.12E+06,	1.53E+06)

# Recruitment deviations versus fcohort relationship
# Create a data frame
df_rec_dev <- data.frame(rec_dev = rec_dev, fcohort = fcohort)

# Fit linear regression model
model_rec_dev <- lm(rec_dev ~ fcohort, data = df_rec_dev)
summary(model_rec_dev)

# Extract p-value from the model with recruitment deviation summary
p_value_rec_dev <- summary(model_rec_dev)$coefficients[2, 4]

# Create the ggplot with regression line and p-value
output_directory <- "C:/GitHub/hake-assessment/Tables_figures_Sandra/"
fcohor_rec_dev <- paste0(output_directory, "ewaa_fcohor_rec_dev.png")
fcohor_RecDev <- ggplot(df_rec_dev, aes(x = fcohort, y = rec_dev)) +
  geom_point() +  # Scatter plot
  geom_smooth(method = "lm", se = TRUE, color = "blue") +  # Add linear regression line
  geom_text(x = max(df_rec_dev$fcohort), y = min(df_rec_dev$rec_dev), label = paste("p =", round(p_value_rec_dev, 4)),
            hjust = 1.2, vjust = -0.8, size = 6) +  # Add p-value to the plot
  labs(title = "Linear Relationship between fcohort and rec_dev") + labs(x="1|fcohort", y="Recruitment deviations")
# Save the plot to the specified directory
ggsave(fcohor_rec_dev, plot = fcohor_RecDev, width = 12, height = 8, units = "in", dpi = 300)

# Display the plot
print(fcohor_RecDev)

# Recruitment versus fcohort relationship
# Create a data frame
df_recr <- data.frame(recr = recr, fcohort = fcohort)

# Fit linear regression model
model_recr <- lm(recr ~ fcohort, data = df_recr)
summary(model_recr)

# Extract p-value from the model with recruitment deviation summary
p_value_recr <- summary(model_recr)$coefficients[2, 4]

# Create the ggplot with regression line and p-value
output_directory <- "C:/GitHub/hake-assessment/Tables_figures_Sandra/"
fcohor_recr <- paste0(output_directory, "ewaa_fcohor_recruitment.png")
fcohor_Recr <- ggplot(df_recr, aes(x = fcohort, y = recr)) +
  geom_point() +  # Scatter plot
  geom_smooth(method = "lm", se = TRUE, color = "blue") +  # Add linear regression line
  geom_text(x = max(df_recr$fcohort), y = min(df_recr$recr), label = paste("p =", round(p_value_recr, 4)),
            hjust = 1.2, vjust = -0.8, size = 6) +  # Add p-value to the plot
  labs(title = "Linear Relationship between fcohort and recruitment") + labs(x="1|fcohort", y="Recruitment")
# Save the plot to the specified directory
ggsave(fcohor_recr, plot = fcohor_Recr, width = 12, height = 8, units = "in", dpi = 300)

# Display the plot
print(fcohor_Recr)



# Fyear  versus Main recruitment deviations. For the fyear random effect only have data from 1977 to 2023
#Recruitment devs from 1977-2022 (From stock assessment)
rec_dev_fyear<-c(-1.64343,-1.98697,0.123949,	2.68975,-1.37203,-1.34778,-0.770718,	2.49053,	-2.01555,	-1.67691,	1.74704,	0.697539,	-1.90878,	1.35377,	0.192922,	-1.8685,	1.11944,	1.15696,	0.275784,	0.59495,	-0.0304116,	0.654935,	2.5325,	-0.953306,	0.1495,	-3.10058,	0.454618,	-3.09334,	0.977308,	0.729547,	-3.53521,	1.71368,	0.318232,	2.70529,	-0.847241,	0.349782,	-0.982123,	1.9944,-3.24806,	1.62556	,0.274458,	-0.947068,	-1.31558,	1.38513	,2.08257,	0.700675)
fyear<-c(0.214,	0.118,	0.211,	0.0327,	0.0217,	-0.0767,	-0.0769,	-0.0352,	0.0699,	0.0429,	-0.048,	-0.00144,	-0.0408,	-0.0326,	-0.0344,	-0.0317,	-0.197,	-0.115,	-0.0503,	-0.0979,	-0.0887,	-0.163,	-0.183,	0.0748,	0.108,	0.126,	0.0109,	-0.0744,	-0.0958,	-0.0829,	-0.146,	-0.0219,	-0.0419,	0.0169,	-0.00546,	0.00696,	0.0746,	0.156,	-0.0549,	-0.0548,	0.0357,	0.093,	-0.0232,	-0.00896,	0.0328,	0.0784)
recr_fyear<-c(4.10E+06,	107392,	892200,	1.15E+07,	198145,	204249,	369849,	9.67E+06,	106920,	150098,	4.62E+06,	1.62E+06,	118726,	3.10E+06,
              966033,	122300,	2.39E+06,	2.47E+06,	1.02E+06,	1.39E+06,	745126,	1.46E+06,	9.34E+06,	289390,	901703,	35902.8,	1.25E+06,	35695,	2.06E+06,
              1.57E+06,	21385.3,	4.08E+06,	983435,	1.10E+07,	313297,	1.06E+06,	295406,	5.86E+06,	30405.4,	3.93E+06,	1.04E+06,	306810,	209274,	3.10E+06,	6.12E+06,	1.53E+06)

# Create a data frame
dfyear <- data.frame(rec_dev_fyear = rec_dev_fyear, fyear = fyear)

# Fit linear regression model
model_dfyear <- lm(rec_dev_fyear ~ fyear, data = dfyear)
summary(model_dfyear)

# Extract p-value from the model summary
p_valuefyear <- summary(model_dfyear)$coefficients[2, 4]

# Create the ggplot with regression line and p-value
output_directory <- "C:/GitHub/hake-assessment/Tables_figures_Sandra/"
fyear_rec_dev <- paste0(output_directory, "ewaa_fyear_rec_dev.png")
fyear_RecDev <- ggplot(model_dfyear, aes(x = fyear, y = rec_dev_fyear)) +
  geom_point() +  # Scatter plot
  geom_smooth(method = "lm", se = TRUE, color = "blue") +  # Add linear regression line
  geom_text(x = 0.22, y = -3, label = paste("p =", round(p_valuefyear, 4)),
            hjust = 1.2, vjust = -0.8, size = 6) +  # Add p-value to the plot
  labs(title = "Linear Relationship between fyear and rec_dev") + labs(x="1|fyear", y="Recruitment deviations")
# Save the plot to the specified directory
ggsave(fyear_rec_dev, plot = fyear_RecDev, width = 12, height = 8, units = "in", dpi = 300)

## Fyear  versus Main recruitment. For the fyear random effect only have data from 1977 to 2023
#Recruitment devs from 1977-2022 (From stock assessment)
# Create a data frame
dfyear_recr <- data.frame(recr_fyear = recr_fyear, fyear = fyear)

# Fit linear regression model
model_dfyear_recr <- lm(recr_fyear ~ fyear, data = dfyear_recr)
summary(model_dfyear_recr)

# Extract p-value from the model summary
p_valuefyear_recr <- summary(model_dfyear_recr)$coefficients[2, 4]

# Create the ggplot with regression line and p-value
output_directory <- "C:/GitHub/hake-assessment/Tables_figures_Sandra/"
fyear_recr <- paste0(output_directory, "ewaa_fyear_recruitment.png")
fyear_Recru <- ggplot(model_dfyear_recr, aes(x = fyear, y = recr_fyear)) +
  geom_point() +  # Scatter plot
  geom_smooth(method = "lm", se = TRUE, color = "blue") +  # Add linear regression line
  geom_text(x = 0.22, y = -3, label = paste("p =", round(p_valuefyear_recr, 4)),
            hjust = 1.2, vjust = -0.8, size = 6) +  # Add p-value to the plot
  labs(title = "Linear Relationship between fyear and recruitment") + labs(x="1|fyear", y="Recruitment")
# Save the plot to the specified directory
ggsave(fyear_recr, plot = fyear_Recru, width = 12, height = 8, units = "in", dpi = 300)

