# USING INDEPENDENT COVARIATES IN EXPERIMENTAL DESIGNS:
# QUANTIFYING THE TRADE-OFF BETWEEN POWER BOOST AND TYPE I ERROR INFLATION

# BY WANG, SPARKS, GONZALES, HESS & LEDGERWOOD

# R CODE WRITTEN BY Y. ANDRE WANG & JOSEPH E. GONZALES
# LAST UPDATE: 07/16/2018

############################# INTRO ###############################

# The purpose of these simulations is to compare the changes in Type I error 
# rate and power from flexible covariate use. 
# We focus on two analytic practices that researchers might use in a 
# two-condition experiment where a promising, unanticipated covariate
# is measured before the experimental manipulation:

# Baseline Practice: Only test the effect of IV (X) on DV (Y)

# Flexible-Covariate Practice: Test the effect of X on Y, 
                # and the effect of X on Y after controlling for covariate (C)

# We created three functions that simulate the different practices
# we considered in the manuscript. Each can be used to estimate the trade-off
# between Type I error inflation and power boost.

# 1. SimData
  # See MAIN SIMULATIONS for details.
  # This function compares the trade-off produced by the baseline practice 
  # and the flexible-covariate practice, which is the trade-off we focus on.

# 2. SuppSimData
  # See ADDITIONAL SIMULATIONS for details.
  # This function compares the trade-off produced by the baseline practice,
  # add the flexible-covariate practice, in addition to other potential
  # practices for including a single, unanticipated covariate we considered.

# 3. SimData5
  # See MULTIPLE COVARIATES for details.
  # This function simulates the situation in which a researcher considers 
  # flexibly including five covariates in his analysis.

###################### MAIN SIMULATIONS ##########################

# The MASS package is used for the simulations
library(MASS)
# If you do not have the package installed, run install.packages('MASS') before
# running the line above.

########### Function for data simulation #############

# The function SimData specifies the simulation inputs. See below for 
# instructions on what inputs are available.
SimData <- function(Nsims, Ncell0, Ncell1 = Ncell0, diff.mu, rcy, SDY0, 
                    SDY1 = SDY0, SDC0, SDC1 = SDC0, alpha = 0.05, seed = 42) {
  
  ########### Data frame for storing simulated data ####
  
  # Create a list of variables
  varlist <- c("N", "N0", "N1", "SDY0", "SDY1", "SDC0", "SDC1", "Mudiff", 
               "Mdiff", "CohenD", "Rcy", "P.Rcy", "P.X.noC", "P.X.C", "P.X.C.I",
               "P.IntXC", "d.P.X.noC", "d.P.X.C")
  
  # Set a data frame for storing data
  temp <- as.data.frame(matrix(NA, nrow = Nsims, ncol = length(varlist)))
  colnames(temp) <- varlist
  
  # Set covariance of Y and C by condition (x = 0 vs. 1)
  Ccy0 <- SDY0 * SDC0 * rcy
  Ccy1 <- SDY1 * SDC1 * rcy
  
  ########### Data simulation and p-value collection loop ####
  
  # Set seed
  set.seed(seed)
  
  # Set conditions
  lapply(1:Nsims, function(i) {
    
    x0 <- as.data.frame(mvrnorm(Ncell0, c(0, 0), 
                                matrix(c(SDY0^2, Ccy0, Ccy0, SDC0^2), 
                                       nrow = 2, ncol = 2, byrow = TRUE), 
                                empirical = FALSE))
    colnames(x0) <- c("y", "c")
    
    x1 <- as.data.frame(mvrnorm(Ncell1, c(diff.mu, 0),
                                matrix(c(SDY1^2, Ccy1, Ccy1, SDC1^2), 
                                       nrow = 2, ncol = 2, byrow = TRUE), 
                                empirical = FALSE))
    colnames(x1) <- c("y", "c")
    dat <- rbind(x0, x1)
    dat$x <- c(rep(0, Ncell0), rep(1, Ncell1))
    
    temp[i, "N"] <<- dim(dat)[1]
    temp[i, "N0"] <<- Ncell0
    temp[i, "N1"] <<- Ncell1
    temp[i, "SDY0"] <<- ysd0 <- sd(dat[which(dat$x == 0), "y"])
    temp[i, "SDY1"] <<- ysd1 <- sd(dat[which(dat$x == 1), "y"])
    temp[i, "SDC0"] <<- csd0 <- sd(dat[which(dat$x == 0), "c"])
    temp[i, "SDC1"] <<- csd1 <- sd(dat[which(dat$x == 1), "c"])
    temp[i, "Mudiff"] <<- diff.mu
    temp[i, "Mdiff"] <<- diff.m <- abs(mean(dat[which(dat$x == 1), "y"]) - 
                                         mean(dat[which(dat$x == 0), "y"]))
    spool <- sqrt(((ysd0^2) * (Ncell0 - 1) + (ysd1^2) * (Ncell1 - 1)) / 
                    (Ncell0 + Ncell1 - 2))
    temp[i, "CohenD"] <<- diff.m/spool
    temp[i, "Rcy"] <<- cor(dat$y, dat$c)
    temp[i, "P.Rcy"] <<- cor.test(dat$y, dat$c)$p.value
    temp[i, "P.X.noC"] <<- coef(summary(lm(y ~ x, data = dat)))[2, 4]
    temp[i, "P.X.C"] <<- coef(summary(lm(y ~ x + c, data = dat)))[2, 4]
    temp[i, "P.X.C.I"] <<- coef(summary(lm(y ~ x + c + x * c, 
                                           data = dat)))[2, 4]
    temp[i, "P.IntXC"] <<- coef(summary(lm(y ~ x + c + x * c, 
                                           data = dat)))[4, 4]
    temp[i, "d.P.X.noC"] <<- coef(summary(lm(y ~ x, data = dat)))[2, 1]
    temp[i, "d.P.X.C"] <<- coef(summary(lm(y ~ x + c, data = dat)))[2, 1]
  })
  
  ######## Baseline Practice #########
  Baseline <- sum(temp$P.X.noC <= alpha)/dim(temp)[1]
  
  ######## Flexible-Covariate Practice ########
  Covariate <- sum(unlist(lapply(1:dim(temp)[1], function(i) {
    min(temp[i, c("P.X.noC", "P.X.C")])
  })) <= alpha)/dim(temp)[1]
  
  ######## Merge results #############
  Practice <- c("Baseline", "Covariate")
  if (diff.mu == 0) {
    TypeIError <- c(Baseline, Covariate)
    Results <- data.frame(Practice, TypeIError)
  } else {
    Power <- c(Baseline, Covariate)
    Results <- data.frame(Practice, Power)
  }
  return(list(Results, temp))
}

############################# Data Generation ###############################

# We ran 10,000 iterations for each condition and simulated data for a large
# number of conditions. Running the full syntax, depending on the speed of your
# CPU, would take many hours. To illustrate our method, the results below are
# generated for N = 160, d = 0.4, rho = .4. This will reproduce the corresponding
# values reported in the manuscript. You can reproduce any value reported in
# Table S1 and S2 (our full set of results) by changing corresponding parameter
# inputs.

##### *-- Power #####
N160_d0.4_rho0.4 <- SimData(Nsims = 10000, Ncell0 = 80, Ncell1 = 80, 
                            diff.mu = 0.4, rcy = 0.4, SDY0 = 1, SDY1 = 1, 
                            SDC0 = 1, SDC1 = 1, seed = 42)

# Nsims: Number of simulations
# Ncell0: Number of participants in Condition 0
# Ncell1: Number of participants in Condition 1 (by default Ncell1 = Ncell0)
# diff.mu: Effect size of x on y in the population (mean group difference)
# rcy: Correlation of c and y in the population
# SDY0: Standard deviation of y in Condition 0 in the population (set to 1)
# SDY1: SD of y in Condition 1 in the population (by default SDY1 = SDY0)
# SDC0: SD of c in Condition 0 in the population (set to 1)
# SDC1: SD of c in Condition 1 in the population (by default SDC1 = SDC0)

# Optionally, you can also set the following inputs:
# alpha: Alpha level (default = 0.05)
# seed: seed value you wish to use (default seed = 42, used in the manuscript)

# See generated data
View(N160_d0.4_rho0.4[[2]])

##### *-- Type I Error Rate #####
N160_d0_rho0.4 <- SimData(Nsims = 10000, Ncell0 = 80, Ncell1 = 80, diff.mu = 0, 
                          rcy = 0.4, SDY0 = 1, SDY1 = 1, SDC0 = 1, SDC1 = 1, 
                          seed = 42)
# We simply set d = 0 to obtain Type I error rate

# See generated data
View(N160_d0_rho0.4[[2]])

############################# Results #######################################

##### *-- Power #####
N160_d0.4_rho0.4[1]
# You should get 0.7131 for Baseline and 0.8111 for Covariate
# Power boost = 0.098

##### *-- Type I Error Rate #####
N160_d0_rho0.4[1]
# You should get 0.0509 for Baseline and 0.0675 for Covariate
# Inflation = 0.0166

##### *-- Number of additional N needed for comparable power boost ####
# Reported in Table 1
# Calculations conducted in G*Power using t-test/Means: Difference
# between two independent means (two groups)/Type of Power Analysis: A priori

##### *-- Accuracy and Precision of the Two Practices ####
# Running the syntax below will reproduce the corresponding values reported
# in Table S3.

dat_N160_d0.4_rho0.4 <- as.data.frame(N160_d0.4_rho0.4[2])
dat_N160_d0_rho0.4 <- as.data.frame(N160_d0_rho0.4[2])

# Mean of the sample effect size estimates by the Baseline Practice
round(mean(dat_N160_d0.4_rho0.4$d.P.X.noC), 3) # 0.402

# Mean of the sample effect size estimates by the Covariate Practice
round(mean(dat_N160_d0.4_rho0.4$d.P.X.C), 3) # 0.402

# SD of the sample effect size estimates by the Baseline Practice
round(sd(dat_N160_d0.4_rho0.4$d.P.X.noC), 3) # 0.159

# SD of the sample effect size estimates by the Covariate Practice
round(sd(dat_N160_d0.4_rho0.4$d.P.X.C), 3) # 0.146


######################## ADDITIONAL SIMULATIONS ##########################

# In this section, the function SuppSimData below allows us to simulate 
# practices in addition to the two practices described above when a researcher
# is considering a single unanticipated, promising covariate. These additional 
# practices are described in Section IV of the Supplemental Materials. 
# For situations in which a researcher is considering multiple unanticipated 
# covariates, see the section "MULTIPLE COVARIATES" below.

# Covariate Practice described by Simmons et al. (2011): 
  # Test (1) the effect of X on Y, 
  # (2) the effect of X on Y controlling for C,
  # (3) the interactive effect of X and C on Y, 
  # and (4) the effect of X on Y controlling for both C
  # and the interaction between C and Y.

# Modified Covariate Practice (A): 
  # Test the effect of X on Y, 
  # and the effect of X on Y controlling for C...
  # only if r(C, Y) >= .3 in the sample.

# Modified Covariate Practice (B):
  # Test the effect of X on Y, 
  # and the effect of X on Y controlling for C...
  # only if r(C, Y) is significantly different from 0 in the sample.

# Refer to the function SimData() above for annotations on our procedure.

SuppSimData <- function(Nsims, Ncell0, Ncell1 = Ncell0, diff.mu, rcy, SDY0, 
                        SDY1 = SDY0, SDC0, SDC1 = SDC0, alpha = 0.05, seed = 42) 
  {
  varlist <- c("N", "N0", "N1", "SDY0", "SDY1", "SDC0", "SDC1", "Mudiff", 
               "Mdiff", "CohenD", "Rcy", "P.Rcy", "P.X.noC", "P.X.C", "P.X.C.I",
               "P.IntXC", "d.P.X.noC", "d.P.X.C")
  temp <- as.data.frame(matrix(NA, nrow = Nsims, ncol = length(varlist)))
  colnames(temp) <- varlist
  Ccy0 <- SDY0 * SDC0 * rcy
  Ccy1 <- SDY1 * SDC1 * rcy
  set.seed(seed)
  lapply(1:Nsims, function(i) {
    x0 <- as.data.frame(mvrnorm(Ncell0, c(0, 0), 
                                matrix(c(SDY0^2, Ccy0, Ccy0, SDC0^2), 
                                       nrow = 2, ncol = 2, byrow = TRUE), 
                                empirical = FALSE))
    colnames(x0) <- c("y", "c")
    x1 <- as.data.frame(mvrnorm(Ncell1, c(diff.mu, 0), 
                                matrix(c(SDY1^2, Ccy1, Ccy1, SDC1^2), 
                                       nrow = 2, ncol = 2, byrow = TRUE), 
                                empirical = FALSE))
    colnames(x1) <- c("y", "c")
    dat <- rbind(x0, x1)
    dat$x <- c(rep(0, Ncell0), rep(1, Ncell1))
    
    temp[i, "N"] <<- dim(dat)[1]
    temp[i, "N0"] <<- Ncell0
    temp[i, "N1"] <<- Ncell1
    temp[i, "SDY0"] <<- ysd0 <- sd(dat[which(dat$x == 0), "y"])
    temp[i, "SDY1"] <<- ysd1 <- sd(dat[which(dat$x == 1), "y"])
    temp[i, "SDC0"] <<- csd0 <- sd(dat[which(dat$x == 0), "c"])
    temp[i, "SDC1"] <<- csd1 <- sd(dat[which(dat$x == 1), "c"])
    temp[i, "Mudiff"] <<- diff.mu
    temp[i, "Mdiff"] <<- diff.m <- abs(mean(dat[which(dat$x == 1), "y"]) - 
                                         mean(dat[which(dat$x == 0), "y"]))
    spool <- sqrt(((ysd0^2) * (Ncell0 - 1) + (ysd1^2) * (Ncell1 - 1)) / 
                    (Ncell0 + Ncell1 - 2))
    temp[i, "CohenD"] <<- diff.m/spool
    temp[i, "Rcy"] <<- cor(dat$y, dat$c)
    temp[i, "P.Rcy"] <<- cor.test(dat$y, dat$c)$p.value
    temp[i, "P.X.noC"] <<- coef(summary(lm(y ~ x, data = dat)))[2, 4]
    temp[i, "P.X.C"] <<- coef(summary(lm(y ~ x + c, data = dat)))[2, 4]
    temp[i, "P.X.C.I"] <<- coef(summary(lm(y ~ x + c + x * c, 
                                           data = dat)))[2, 4]
    temp[i, "P.IntXC"] <<- coef(summary(lm(y ~ x + c + x * c, 
                                           data = dat)))[4, 4]
    temp[i, "d.P.X.noC"] <<- coef(summary(lm(y ~ x, data = dat)))[2, 1]
    temp[i, "d.P.X.C"] <<- coef(summary(lm(y ~ x + c, data = dat)))[2, 1]
  })
  
  ######## Baseline Practice #########
  Baseline <- sum(temp$P.X.noC <= alpha)/dim(temp)[1]
  
  ######## Flexible-Covariate Practice ########
  Covariate <- sum(unlist(lapply(1:dim(temp)[1], function(i) {
    min(temp[i, c("P.X.noC", "P.X.C")])
  })) <= alpha)/dim(temp)[1]
  
  ######## Simmons et al. (2011) Practice #####
  Simmons <- sum(unlist(lapply(1:dim(temp)[1], function(i) {
    if (temp[i, "P.IntXC"] <= alpha) {
      temp[i, "P.IntXC"]
    } else {
      min(temp[i, c("P.X.noC", "P.X.C", "P.X.C.I")])
    }
  })) <= alpha)/dim(temp)[1]
  
  ######## Modified Covariate Practice (a) ####
  Modified.A <- sum(unlist(lapply(1:dim(temp)[1], function(i) {
    if (temp[i, "Rcy"] >= 0.3) {
      min(temp[i, c("P.X.noC", "P.X.C")])
    } else {
      temp[i, "P.X.noC"]
    }
  })) <= alpha)/dim(temp)[1]
  
  ######## Modified Covariate Practice (b) ####
  Modified.B <- sum(unlist(lapply(1:dim(temp)[1], function(i) {
    if (temp[i, "P.Rcy"] <= 0.05) {
      min(temp[i, c("P.X.noC", "P.X.C")])
    } else {
      temp[i, "P.X.noC"]
    }
  })) <= alpha)/dim(temp)[1]
  
  ######## Merge results ########
  Practice <- c("Baseline", "Covariate", "Simmons", 
                "Modified.A", "Modified.B")
  if (diff.mu == 0) {
    TypeIError <- c(Baseline, Covariate, Simmons, Modified.A, Modified.B)
    Results <- data.frame(Practice, TypeIError)
  } else {
    Power <- c(Baseline, Covariate, Simmons, Modified.A, Modified.B)
    Results <- data.frame(Practice, Power)
  }
  return(list(Results, temp))
}

############################# Data Generation ###############################

# Data Generation for N = 160, rho = .4
# Type I error rate
Supp_N160_d0_rho0.4 <- SuppSimData(10000, 80, 80, 0, 0.4, 1, 1, 1, 1)
# Power (d = 0.4)
Supp_N160_d0.4_rho0.4 <- SuppSimData(10000, 80, 80, 0.4, 0.4, 1, 1, 1, 1)

# Data Generation for N = 160, rho = .2
# Type I error rate
Supp_N160_d0_rho0.2 <- SuppSimData(10000, 80, 80, 0, 0.2, 1, 1, 1, 1)
# Power (d = 0.4)
Supp_N160_d0.4_rho0.2 <- SuppSimData(10000, 80, 80, 0.4, 0.2, 1, 1, 1, 1)


############################# Results #######################################

#### *-- Flexibly Including Both C and the X*C Interaction ####

# Power
Supp_N160_d0.4_rho0.4[1]
# You should get 0.8111 for the Covariate Practice, and
# 0.8223 for the Simmons approach (change = 0.0112)

# Type I error rate
Supp_N160_d0_rho0.4[1]
# You should get 0.0675 for the Covariate Practice, and
# 0.1182 for the Simmons approach (change = 0.0507)

#### *-- Using the Observed Correlation in the Sample ####

# We explored two variations of the Covariate Practice that use the observed
# correlation between C and Y in the sample:
# Modified Covariate Practice (A), and Modified Covariate Practice (B). 
# See the function SuppSimData() above for details.

# When rho < .3: Power
Supp_N160_d0.4_rho0.2[1]
# You should get 0.7166 for Modified (A), and 0.7413 for Modified (B).
# Modified (A) performs similarly to Baseline.
# Modified (B) performs similarly to Covariate.

# When rho < .3: Type I error rate
Supp_N160_d0_rho0.2[1]
# You should get 0.0499 for Modified (A), and 0.0575 for Modified (B).
# Modified (A) performs similarly to Baseline. 
# Modified (B) performs similarly to Covariate.

# When rho > .3: Power
Supp_N160_d0.4_rho0.4[1]
# You should get 0.8033 for Modified (A), and 0.8110 for Modified (B).
# Both Modified (A) and Modified (B) perform similarly to Covariate.

# When rho > .3: Type I error rate
Supp_N160_d0_rho0.4[1]
# You should get 0.0668 for Modified (A), and 0.0675 for Modified (B).
# Both Modified (A) and Modified (B) perform similarly to Covariate.

# As you can see, neither Modified (A) nor Modified (B) offers a more favorable
# trade-off than the Baseline or the Covariate practice.
# Using the observed correlation in the sample does not offer protection against
# Type I error inflation. It can undercut power boost.

#### *-- Adjusting Alpha Level ####

dat_Supp_N160_d0_rho0.4 <- as.data.frame(Supp_N160_d0_rho0.4[2])

# Find the alpha level that keeps Type I error rate at .05
dat_Supp_N160_d0_rho0.4$cov.pvalue <- pmin(dat_Supp_N160_d0_rho0.4[, 13], 
                                           dat_Supp_N160_d0_rho0.4[, 14])

adjusted.alpha <- sort(dat_Supp_N160_d0_rho0.4$cov.pvalue, 
                       partial = dim(dat_Supp_N160_d0_rho0.4)[1] * 
                         0.05)[dim(dat_Supp_N160_d0_rho0.4)[1] * 0.05]

adj1.N160_d0.4_rho0.4 <- SimData(10000, 80, 80, 0.4, 0.4, 1, 1, 1, 1, 
                                 alpha = adjusted.alpha)
adj1.N160_d0_rho0.4 <- SimData(10000, 80, 80, 0, 0.4, 1, 1, 1, 1, 
                               alpha = adjusted.alpha)

# Type I Error Rate
adj1.N160_d0_rho0.4[1]
# As expected, Type I error rate is exactly .05 now for the Covariate Practice.

# Power
adj1.N160_d0.4_rho0.4[1]
# You should get 0.7716 for the Covariate Practice.
# Note that since we used the alpha adjustment,
# the baseline power is not the baseline output produced here;
# it should be the one prior to adjustment (0.7131).
# Power boost = 0.7716 - 0.7131 = 0.0585.

# Alternatively, we can directly use a smaller alpha value (e.g., alpha = .033)
adj2.N160_d0_rho0.4 <- SimData(10000, 80, 80, 0, 0.4, 1, 1, 1, 1, 
                               alpha = 0.033)
adj2.N160_d0.4_rho0.4 <- SimData(10000, 80, 80, 0.4, 0.4, 1, 1, 1, 1, 
                                 alpha = 0.033)

# Type I Error Rate
adj2.N160_d0_rho0.4[1]
# Type I error rate is now 0.0457.

# Power
adj2.N160_d0.4_rho0.4[1]
# You should get 0.7546 for the Covariate Practice.
# Power boost = 0.7546 - 0.7131 = 0.0415.

# What happens when rho is smaller (e.g., rho = .1)?
adj2.N160_d0.4_rho0.1 <- SimData(10000, 80, 80, 0.4, 0.1, 1, 1, 1, 1, 
                                 alpha = 0.033)
adj2.N160_d0_rho0.1 <- SimData(10000, 80, 80, 0, 0.1, 1, 1, 1, 1, 
                               alpha = 0.033)

# Type I Error Rate
adj2.N160_d0_rho0.1[1]
# Type I error rate is now 0.0367.

# Power
adj2.N160_d0.4_rho0.1[1]
# You should get 0.6663 for the Covariate Practice.
# We actually have a power loss of 0.7131 - 0.6663 = 0.0468 now.


######################## MULTIPLE COVARIATES ##########################

# In this section, we use the function SimData5 to simulate a situation
# in which a researcher considers flexibly including five covariates in 
# his analysis. Refer to SimData for the core parameters. Additional
# inputs are explained in the "Data Generation" section below.

SimData5 <- function(Nsims, Ncell0, Ncell1 = Ncell0, diff.mu, rc1y, rc2y, rc3y,
                     rc4y, rc5y, rc1c2, rc1c3, rc1c4, rc1c5, rc2c3, rc2c4, 
                     rc2c5, rc3c4, rc3c5, rc4c5, SDY0, SDY1 = SDY0, SDC1_0, 
                     SDC1_1 = SDC1_0, SDC2_0, SDC2_1 = SDC2_0, SDC3_0, 
                     SDC3_1 = SDC3_0, SDC4_0, SDC4_1 = SDC4_0, SDC5_0, 
                     SDC5_1 = SDC5_0, alpha = 0.05, seed = 42) {
  
  ########### Data frame for storing simulated data ####
  
  # Create a list of variables
  varlist <- c("N", "N0", "N1", "SDY0", "SDY1", "SDC1_0", "SDC1_1", "SDC2_0", 
               "SDC2_1", "SDC3_0", "SDC3_1", "SDC4_0", "SDC4_1", "SDC5_0", 
               "SDC5_1","Mudiff", "Mdiff", "CohenD", "Rc1y", "Rc2y", "Rc3y",
               "Rc4y", "Rc5y", "P.Rc1y", "P.Rc2y", "P.Rc3y", "P.Rc4y", "P.Rc5y", 
               "P.X.noC", "P.X.C1", "P.X.C2", "P.X.C3", "P.X.C4", "P.X.C5",
               "P.X.C1.I", "P.X.C2.I", "P.X.C3.I", "P.X.C4.I", "P.X.C5.I",
               "P.IntXC1", "P.IntXC2", "P.IntXC3", "P.IntXC4", "P.IntXC5",
               "d.P.X.noC", "d.P.X.C1", "d.P.X.C2", "d.P.X.C3", "d.P.X.C4",
               "d.P.X.C5")
  
  # Set a dataframe for storing data
  temp <- as.data.frame(matrix(NA, nrow = Nsims, ncol = length(varlist)))
  colnames(temp) <- varlist
  
  # Set covariances of Y and Cs by condition (x = 0 vs. 1)
  Cc1y0 <- SDY0 * SDC1_0 * rc1y
  Cc1y1 <- SDY1 * SDC1_1 * rc1y
  
  Cc2y0 <- SDY0 * SDC2_0 * rc2y
  Cc2y1 <- SDY1 * SDC2_1 * rc2y
  
  Cc3y0 <- SDY0 * SDC3_0 * rc3y
  Cc3y1 <- SDY1 * SDC3_1 * rc3y
  
  Cc4y0 <- SDY0 * SDC4_0 * rc4y
  Cc4y1 <- SDY1 * SDC4_1 * rc4y
  
  Cc5y0 <- SDY0 * SDC5_0 * rc5y
  Cc5y1 <- SDY1 * SDC5_1 * rc5y
  
  # Set covariances among Cs by condition (x = 0 vs. 1)
  
  C0c1c2 <- SDC1_0 * SDC2_0 * rc1c2
  C0c1c3 <- SDC1_0 * SDC3_0 * rc1c3
  C0c1c4 <- SDC1_0 * SDC4_0 * rc1c4
  C0c1c5 <- SDC1_0 * SDC5_0 * rc1c5
  C0c2c3 <- SDC2_0 * SDC3_0 * rc2c3
  C0c2c4 <- SDC2_0 * SDC4_0 * rc2c4
  C0c2c5 <- SDC2_0 * SDC5_0 * rc2c5
  C0c3c4 <- SDC3_0 * SDC4_0 * rc3c4
  C0c3c5 <- SDC3_0 * SDC5_0 * rc3c5
  C0c4c5 <- SDC4_0 * SDC5_0 * rc4c5
  
  C1c1c2 <- SDC1_1 * SDC2_1 * rc1c2
  C1c1c3 <- SDC1_1 * SDC3_1 * rc1c3
  C1c1c4 <- SDC1_1 * SDC4_1 * rc1c4
  C1c1c5 <- SDC1_1 * SDC5_1 * rc1c5
  C1c2c3 <- SDC2_1 * SDC3_1 * rc2c3
  C1c2c4 <- SDC2_1 * SDC4_1 * rc2c4
  C1c2c5 <- SDC2_1 * SDC5_1 * rc2c5
  C1c3c4 <- SDC3_1 * SDC4_1 * rc3c4
  C1c3c5 <- SDC3_1 * SDC5_1 * rc3c5
  C1c4c5 <- SDC4_1 * SDC5_1 * rc4c5
  
  ########### Data simulation and p-value collection loop ####
  
  # Set seed
  set.seed(seed)
  
  # Set conditions
  lapply(1:Nsims, function(i) {
    
    x0 <- as.data.frame(
      mvrnorm(Ncell0, c(rep(0, 6)), 
              matrix(c(SDY0^2, Cc1y0, Cc2y0, Cc3y0, Cc4y0, Cc5y0,
                       Cc1y0, SDC1_0^2, C0c1c2, C0c1c3, C0c1c4, C0c1c5,
                       Cc2y0, C0c1c2, SDC2_0^2, C0c2c3, C0c2c4, C0c2c5,
                       Cc3y0, C0c1c3, C0c2c3, SDC3_0^2, C0c3c4, C0c3c5,
                       Cc4y0, C0c1c4, C0c2c4, C0c3c4, SDC4_0^2, C0c4c5,
                       Cc5y0, C0c1c5, C0c2c5, C0c3c5, C0c4c5, SDC5_0^2), 
                     nrow = 6, ncol = 6, byrow = TRUE), empirical = FALSE))
    colnames(x0) <- c("y", "c1", "c2", "c3", "c4", "c5")
    
    x1 <- as.data.frame(
      mvrnorm(Ncell1, c(diff.mu, rep(0, 5)), 
              matrix(c(SDY1^2, Cc1y1, Cc2y1, Cc3y1, Cc4y1, Cc5y1,
                       Cc1y1, SDC1_1^2, C1c1c2, C1c1c3, C1c1c4, C1c1c5,
                       Cc2y1, C1c1c2, SDC2_1^2, C1c2c3, C1c2c4, C1c2c5,
                       Cc3y1, C1c1c3, C1c2c3, SDC3_1^2, C1c3c4, C1c3c5,
                       Cc4y1, C1c1c4, C1c2c4, C1c3c4, SDC4_1^2, C1c4c5,
                       Cc5y1, C1c1c5, C1c2c5, C1c3c5, C1c4c5, SDC5_1^2), 
                     nrow = 6, ncol = 6, byrow = TRUE),
              empirical = FALSE))
    colnames(x1) <- c("y", "c1", "c2", "c3", "c4", "c5")
    
    dat <- rbind(x0, x1)
    dat$x <- c(rep(0, Ncell0), rep(1, Ncell1))
    
    temp[i, "N"] <<- dim(dat)[1]
    temp[i, "N0"] <<- Ncell0
    temp[i, "N1"] <<- Ncell1
    temp[i, "SDY0"] <<- ysd0 <- sd(dat[which(dat$x == 0), "y"])
    temp[i, "SDY1"] <<- ysd1 <- sd(dat[which(dat$x == 1), "y"])
    temp[i, "SDC1_0"] <<- sd(dat[which(dat$x == 0), "c1"])
    temp[i, "SDC1_1"] <<- sd(dat[which(dat$x == 1), "c1"])
    temp[i, "SDC2_0"] <<- sd(dat[which(dat$x == 0), "c2"])
    temp[i, "SDC2_1"] <<- sd(dat[which(dat$x == 1), "c2"])
    temp[i, "SDC3_0"] <<- sd(dat[which(dat$x == 0), "c3"])
    temp[i, "SDC3_1"] <<- sd(dat[which(dat$x == 1), "c3"])
    temp[i, "SDC4_0"] <<- sd(dat[which(dat$x == 0), "c4"])
    temp[i, "SDC4_1"] <<- sd(dat[which(dat$x == 1), "c4"])
    temp[i, "SDC5_0"] <<- sd(dat[which(dat$x == 0), "c5"])
    temp[i, "SDC5_1"] <<- sd(dat[which(dat$x == 1), "c5"])
    temp[i, "Mudiff"] <<- diff.mu
    temp[i, "Mdiff"] <<- diff.m <- abs(mean(dat[which(dat$x == 1), "y"]) - 
                                         mean(dat[which(dat$x == 0), "y"]))
    spool <- sqrt(((ysd0^2) * (Ncell0 - 1) + (ysd1^2) * (Ncell1 - 1)) / 
                    (Ncell0 + Ncell1 - 2))
    temp[i, "CohenD"] <<- diff.m/spool
    temp[i, "Rc1y"] <<- cor(dat$y, dat$c1)
    temp[i, "Rc2y"] <<- cor(dat$y, dat$c2)
    temp[i, "Rc3y"] <<- cor(dat$y, dat$c3)
    temp[i, "Rc4y"] <<- cor(dat$y, dat$c4)
    temp[i, "Rc5y"]  <<- cor(dat$y, dat$c5)
    temp[i, "P.Rc1y"] <<- cor.test(dat$y, dat$c1)$p.value
    temp[i, "P.Rc2y"] <<- cor.test(dat$y, dat$c2)$p.value
    temp[i, "P.Rc3y"] <<- cor.test(dat$y, dat$c3)$p.value
    temp[i, "P.Rc4y"] <<- cor.test(dat$y, dat$c4)$p.value
    temp[i, "P.Rc5y"] <<- cor.test(dat$y, dat$c5)$p.value
    temp[i, "P.X.noC"] <<- coef(summary(lm(y ~ x, data = dat)))[2, 4]
    temp[i, "P.X.C1"] <<- coef(summary(lm(y ~ x + c1, data = dat)))[2, 4]
    temp[i, "P.X.C2"] <<- coef(summary(lm(y ~ x + c2, data = dat)))[2, 4]
    temp[i, "P.X.C3"] <<- coef(summary(lm(y ~ x + c3, data = dat)))[2, 4]
    temp[i, "P.X.C4"] <<- coef(summary(lm(y ~ x + c4, data = dat)))[2, 4]
    temp[i, "P.X.C5"] <<- coef(summary(lm(y ~ x + c5, data = dat)))[2, 4]
    temp[i, "P.X.C1.I"] <<- coef(summary(lm(y ~ x + c1 + x * c1, 
                                            data = dat)))[2, 4]
    temp[i, "P.X.C2.I"] <<- coef(summary(lm(y ~ x + c2 + x * c2, 
                                            data = dat)))[2, 4]
    temp[i, "P.X.C3.I"] <<- coef(summary(lm(y ~ x + c3 + x * c3, 
                                            data = dat)))[2, 4]
    temp[i, "P.X.C4.I"] <<- coef(summary(lm(y ~ x + c4 + x * c4, 
                                            data = dat)))[2, 4]
    temp[i, "P.X.C5.I"] <<- coef(summary(lm(y ~ x + c5 + x * c5, 
                                            data = dat)))[2, 4]
    temp[i, "P.IntXC1"] <<- coef(summary(lm(y ~ x + c1 + x * c1, 
                                            data = dat)))[4, 4]
    temp[i, "P.IntXC2"] <<- coef(summary(lm(y ~ x + c2 + x * c2, 
                                            data = dat)))[4, 4]
    temp[i, "P.IntXC3"] <<- coef(summary(lm(y ~ x + c3 + x * c3, 
                                            data = dat)))[4, 4]
    temp[i, "P.IntXC4"] <<- coef(summary(lm(y ~ x + c4 + x * c4, 
                                            data = dat)))[4, 4]
    temp[i, "P.IntXC5"] <<- coef(summary(lm(y ~ x + c5 + x * c5, 
                                            data = dat)))[4, 4]
    temp[i, "d.P.X.noC"] <<- coef(summary(lm(y ~ x, data = dat)))[2, 1]
    temp[i, "d.P.X.C1"] <<- coef(summary(lm(y ~ x + c1, data = dat)))[2, 1]
    temp[i, "d.P.X.C2"] <<- coef(summary(lm(y ~ x + c2, data = dat)))[2, 1]
    temp[i, "d.P.X.C3"] <<- coef(summary(lm(y ~ x + c3, data = dat)))[2, 1]
    temp[i, "d.P.X.C4"] <<- coef(summary(lm(y ~ x + c4, data = dat)))[2, 1]
    temp[i, "d.P.X.C5"] <<- coef(summary(lm(y ~ x + c5, data = dat)))[2, 1]
  })
  
  ######## Baseline Practice #########
  Baseline <- sum(temp$P.X.noC <= alpha)/dim(temp)[1]
  
  ######## Covariate Practice ########
  Covariate <- sum(unlist(lapply(1:dim(temp)[1], function(i) {
    min(temp[i, c("P.X.noC", "P.X.C1", "P.X.C2", "P.X.C3", "P.X.C4", "P.X.C5")])
  })) <= alpha)/dim(temp)[1]
  
  ######## Merge results #############
  Practice <- c("Baseline", "Covariate")
  if (diff.mu == 0) {
    TypeIError <- c(Baseline, Covariate)
    Results <- data.frame(Practice, TypeIError)
  } else {
    Power <- c(Baseline, Covariate)
    Results <- data.frame(Practice, Power)
  }
  return(list(Results, temp))
}


############################# Data Generation ###############################

# *-- N = 200, rho = .4 (C1, C2) and .2 (C3 - C5) ####
# Type I error rate, with intercorrelations similar to those among Big Five
Sim5_N200_d0_rho44222_b5 <- SimData5(
  Nsims = 10000, Ncell0 = 100, diff.mu = 0, 
  rc1y = 0.4, rc2y = 0.4, rc3y = 0.2, rc4y = 0.2, rc5y = 0.2,
  rc1c2 = .20, rc1c3 = .43, rc1c4 = .21, rc1c5 = -.17, rc2c3 = .29, 
  rc2c4 = .43, rc2c5 = -.43, rc3c4 = .26, rc3c5 = -.36, rc4c5 = -.36,
  SDY0 = 1, SDC1_0 = 1, SDC2_0 = 1, SDC3_0 = 1, SDC4_0 = 1, SDC5_0 = 1)

# rc1y, rc2y, rc3y, rc4y, rc5y: population correlation between each C and Y
# rc1c2, rc1c3, rc1c4, rc1c5, rc2c3, rc2c4, rc2c5, rc3c4, rc3c5, rc4c5:
# population intercorrelations among Cs
# SDY0, SDC1_0, SDC2_0, SDC3_0, SDC4_0, SDC5_0: SDs of Y and Cs
# By default SDs of all variables are set to be equal across conditions
# But you can change them by specifying:
# SDY1, SDC1_1, SDC2_1, SDC3_1, SDC4_1, SDC5_1

# Power (d = 0.4), with intercorrelations similar to those among Big Five
Sim5_N200_d0.4_rho44222_b5 <- SimData5(
  Nsims = 10000, Ncell0 = 100, diff.mu = 0.4, 
  rc1y = 0.4, rc2y = 0.4, rc3y = 0.2, rc4y = 0.2, rc5y = 0.2, 
  rc1c2 = .20, rc1c3 = .43, rc1c4 = .21, rc1c5 = -.17, rc2c3 = .29, 
  rc2c4 = .43, rc2c5 = -.43, rc3c4 = .26, rc3c5 = -.36, rc4c5 = -.36, 
  SDY0 = 1, SDC1_0 = 1, SDC2_0 = 1, SDC3_0 = 1, SDC4_0 = 1, SDC5_0 = 1)

# Type I error rate, with intercorrelations = 0
Sim5_N200_d0_rho44222_in <- SimData5(
  Nsims = 10000, Ncell0 = 100, diff.mu = 0, 
  rc1y = 0.4, rc2y = 0.4, rc3y = 0.2, rc4y = 0.2, rc5y = 0.2,
  rc1c2 = 0, rc1c3 = 0, rc1c4 = 0, rc1c5 = 0, rc2c3 = 0, 
  rc2c4 = 0, rc2c5 = 0, rc3c4 = 0, rc3c5 = 0, rc4c5 = 0,
  SDY0 = 1, SDC1_0 = 1, SDC2_0 = 1, SDC3_0 = 1, SDC4_0 = 1, SDC5_0 = 1)

# Power (d = 0.4), with intercorrelations = 0
Sim5_N200_d0.4_rho44222_in <- SimData5(
  Nsims = 10000, Ncell0 = 100, diff.mu = 0.4, 
  rc1y = 0.4, rc2y = 0.4, rc3y = 0.2, rc4y = 0.2, rc5y = 0.2,
  rc1c2 = 0, rc1c3 = 0, rc1c4 = 0, rc1c5 = 0, rc2c3 = 0, 
  rc2c4 = 0, rc2c5 = 0, rc3c4 = 0, rc3c5 = 0, rc4c5 = 0,
  SDY0 = 1, SDC1_0 = 1, SDC2_0 = 1, SDC3_0 = 1, SDC4_0 = 1, SDC5_0 = 1)

# *-- N = 200, rhos = .5 ####
# Type I error rate, with intercorrelations = .1
Sim5_N200_d0_rho55555_in <- SimData5(
  Nsims = 10000, Ncell0 = 100, diff.mu = 0, 
  rc1y = 0.5, rc2y = 0.5, rc3y = 0.5, rc4y = 0.5, rc5y = 0.5,
  rc1c2 = 0.1, rc1c3 = 0.1, rc1c4 = 0.1, rc1c5 = 0.1, rc2c3 = 0.1, 
  rc2c4 = 0.1, rc2c5 = 0.1, rc3c4 = 0.1, rc3c5 = 0.1, rc4c5 = 0.1,
  SDY0 = 1, SDC1_0 = 1, SDC2_0 = 1, SDC3_0 = 1, SDC4_0 = 1, SDC5_0 = 1)

# Power (d = 0.4), with intercorrelations = .1
Sim5_N200_d0.4_rho55555_in <- SimData5(
  Nsims = 10000, Ncell0 = 100, diff.mu = 0.4, 
  rc1y = 0.5, rc2y = 0.5, rc3y = 0.5, rc4y = 0.5, rc5y = 0.5,
  rc1c2 = 0.1, rc1c3 = 0.1, rc1c4 = 0.1, rc1c5 = 0.1, rc2c3 = 0.1, 
  rc2c4 = 0.1, rc2c5 = 0.1, rc3c4 = 0.1, rc3c5 = 0.1, rc4c5 = 0.1,
  SDY0 = 1, SDC1_0 = 1, SDC2_0 = 1, SDC3_0 = 1, SDC4_0 = 1, SDC5_0 = 1)

############################# Results #######################################

# *-- N = 200, rho = .4 (C1, C2) and .2 (C3 - C5) ####
# Type I error rate (Big Five intercorrelations)
Sim5_N200_d0_rho44222_b5[1]
# You should get 0.0497 for Baseline and 0.1010 for Covariate
# Inflation = 0.1010 - 0.0497 = 0.0513

# Power (d = 0.4, Big Five intercorrelations)
Sim5_N200_d0.4_rho44222_b5[1]
# You should get 0.8020 for Baseline and 0.9197 for Covariate
# Power boost = 0.9197 - 0.8020 = 0.1177

# Type I error rate (covariates independent from each other)
Sim5_N200_d0_rho44222_in[1]
# You should get 0.0504 for Baseline and 0.1020 for Covariate
# Inflation = 0.1020 - 0.0504 = 0.0516

# Power (d = 0.4, covariates independent from each other)
Sim5_N200_d0.4_rho44222_in[1]
# You should get 0.9245 for Baseline and 0.7966 for Covariate
# Power boost = 0.9245 - 0.7966 = 0.1279

# *-- N = 200, rhos = .5 ####
# Type I error rate (with intercorrelations = .1)
Sim5_N200_d0_rho55555_in[1]
# You should get 0.0531 for Baseline and 0.1545 for Covariate
# Inflation = 0.1545 - 0.0531 = 0.1014

# Power (d = 0.4, with intercorrelations = .1)
Sim5_N200_d0.4_rho55555_in[1]
# You should get 0.7977 for Baseline and 0.9816 for Covariate
# Power boost = 0.9816 - 0.7977 = 0.1839
