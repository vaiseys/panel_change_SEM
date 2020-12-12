# load packages
library(tidyverse)
library(broom)
library(lavaan)

# load, reshape, and nest data
load("gsslong.Rdata")
d <- gss.long %>% 
  pivot_longer(abany:xmovie, names_to = "variable", values_to = "y") %>% 
  pivot_wider(names_from = wave, values_from = y, names_prefix = "y") %>%
  arrange(variable, ds, idnum) %>%
  drop_na() %>% 
  group_by(variable, ds) %>% 
  nest()

# model formulas

## AUM1
aum1_mod <- "

  # structural part 
    y3 ~ alpha*1 + rho*y2 
    y2 ~ alpha*1 + rho*y1
    
  # measurement part
    U =~ 1*y3 + 1*y2

  # variances
    y3 ~~ v*y3
    y2 ~~ v*y2
    
  # covariance (constrained)
    U ~~ 0*y1
    
  # constraint var(U) = 0
    U ~~ 0*U
    
"

## AUM2
aum2_mod <- "

  # structural part 
    y3 ~ alpha3*1 + rho*y2 
    y2 ~ alpha2*1 + rho*y1
    
  # measurement part
    U =~ 1*y3 + 1*y2

  # variances
    y3 ~~ v*y3
    y2 ~~ v*y2
    
  # covariance (constrained)
    U ~~ 0*y1
    
  # constraint var(U) = 0
    U ~~ 0*U
    
"

## AUM3
aum3_mod <- "

  # structural part 
    y3 ~ alpha*1 + rho*y2 
    y2 ~ alpha*1 + rho*y1
    
  # measurement part
    U =~ 1*y3 + 1*y2

  # variances
    y3 ~~ v*y3
    y2 ~~ v*y2
    
  # covariance
    U ~~ tau*y1
    
  # var(U)
    U ~~ U
  
"

## AUM4
aum4_mod <- "

  # structural part 
    y3 ~ alpha3*1 + rho*y2 
    y2 ~ alpha2*1 + rho*y1
    
  # measurement part
    U =~ 1*y3 + 1*y2

  # variances
    y3 ~~ v*y3
    y2 ~~ v*y2
    
  # covariance (tau)
    U ~~ tau*y1
    
  # var(U)
    U ~~ U
  
"

## SDM1
sdm1_mod <- "

  # structural part 
    y3 ~ alpha*1 + rho*y2 
    y2 ~ alpha*1 + rho*y1
    
  # measurement part
    U =~ 1*y3 + 1*y2

  # variances
    y3 ~~ v*y3
    y2 ~~ v*y2
    
  # covariance (tau)
    U ~~ tau*y1
    
  # var(U)
    U ~~ U
    
  # constraint
    rho == 0
  
"

## SDM2
sdm2_mod <- "

  # structural part 
    y3 ~ alpha3*1 + rho*y2 
    y2 ~ alpha2*1 + rho*y1
    
  # measurement part
    U =~ 1*y3 + 1*y2

  # variances
    y3 ~~ v*y3
    y2 ~~ v*y2
    
  # covariance (tau)
    U ~~ tau*y1
    
  # var(U)
    U ~~ U
    
  # constraint
    rho == 0
  
"

## CFA
cfa_mod <- "

  # measurement part
    U =~ 1*y3 + 1*y2 + 1*y1

  # variances
    y3 ~~ v*y3
    y2 ~~ v*y2
    y1 ~~ v*y1
    
  # var(U)
    U ~~ U

  # intercepts
    y1 ~ int*1
    y2 ~ int*1
    y3 ~ int*1

"

# fitting functions
aum1_fit <- function (x) {
  sem(aum1_mod,
      sampling.weights = "wtpannr123",
      check.post = TRUE,
      se = "none",
      data = x)
}

aum2_fit <- function (x) {
  sem(aum2_mod,
      sampling.weights = "wtpannr123",
      check.post = TRUE,
      se = "none",
      data = x)
}

aum3_fit <- function (x) {
  sem(aum3_mod,
      sampling.weights = "wtpannr123",
      check.post = TRUE,
      se = "none",
      data = x)
}

aum4_fit <- function (x) {
  sem(aum4_mod,
      sampling.weights = "wtpannr123",
      check.post = TRUE,
      se = "none",
      data = x)
}

sdm1_fit <- function (x) {
  sem(sdm1_mod,
      sampling.weights = "wtpannr123",
      check.post = TRUE,
      se = "none",
      data = x)
}

sdm2_fit <- function (x) {
  sem(sdm2_mod,
      sampling.weights = "wtpannr123",
      check.post = TRUE,
      se = "none",
      data = x)
}

cfa_fit <- function (x) {
  sem(cfa_mod,
      sampling.weights = "wtpannr123",
      check.post = TRUE,
      se = "none",
      data = x)
}

p_aum1_fit <- possibly(aum1_fit, otherwise = "NOPE")
p_aum2_fit <- possibly(aum2_fit, otherwise = "NOPE")
p_aum3_fit <- possibly(aum3_fit, otherwise = "NOPE")
p_aum4_fit <- possibly(aum4_fit, otherwise = "NOPE")
p_sdm1_fit <- possibly(sdm1_fit, otherwise = "NOPE")
p_sdm2_fit <- possibly(sdm2_fit, otherwise = "NOPE")
p_cfa_fit <- possibly(cfa_fit, otherwise = "NOPE")


d <- d %>% 
  mutate(AUM1 = map(data, p_aum1_fit),
         AUM2 = map(data, p_aum2_fit),
         AUM3 = map(data, p_aum3_fit),
         AUM4 = map(data, p_aum4_fit),
         SDM1 = map(data, p_sdm1_fit),
         SDM2 = map(data, p_sdm2_fit),
         CFA = map(data, p_cfa_fit))


# Results
## glance function to extract
p_glance = possibly(glance, otherwise = "NOPE")

## add model summary stuff
results <- d %>%
  pivot_longer(AUM1:CFA, names_to = "mod_spec", values_to = "mod_object") %>% 
  mutate(glanced = map(mod_object, p_glance)) %>% 
  unnest(glanced) %>% 
  select(variable, ds, mod_spec, BIC, chisq, npar, rmsea, converged, nobs) %>% 
  mutate(type = case_when(
    mod_spec == "SDM1" | mod_spec == "SDM2" ~ "SDM",
    mod_spec == "CFA" ~ "CFA",
    TRUE ~ "AUM")) %>% 
  ungroup()

save(results, file = "results.Rdata")
