# USING INDEPENDENT COVARIATES IN EXPERIMENTAL DESIGNS:
# QUANTIFYING THE TRADE-OFF BETWEEN TYPE I ERROR INFLATION AND POWER BOOST

# BY WANG, SPARKS, GONZALES, HESS & LEDGERWOOD

# R CODE WRITTEN BY Y. ANDRE WANG & JOSEPH E. GONZALES
# LAST UPDATE: 11/29/2016

############################# INTRO ###############################

# The purpose of these simulations is to compare the changes in
# Type I error rate and power from flexible covariate use. 
# We focus on two analytic strategies that researchers might use 
# in a two-condition experiment where a covariate is measured 
# before the experimental manipulation:

# Baseline Strategy: Only test the effect of IV (x) on DV (y)

# Covariate Strategy: Test the effect of x on y, 
                # and the effect of x on y after controlling for covariate (c)

###################### SIMULATIONS SETUP ##########################

# The MASS package is used for the simulations
library(MASS)
# If you do not have the package installed, run install.packages('MASS') before
# running the line above.

########### Function for data simulation #############

# The function below specifies the simulation inputs. See below for instructions
# on what inputs are available.
SimData <- function(Nsims, Ncell0, Ncell1 = Ncell0, diff.mu, rcy, SDY0, 
                    SDY1 = SDY0, SDC0, SDC1 = SDC0, alpha = 0.05, seed = 42) {
  
  ########### Data frame for storing simulated data ####
  
  # Create a list of variables
  varlist <- c("N", "N0", "N1", "SDY0", "SDY1", "SDC0", "SDC1", "Mudiff", 
               "Mdiff", "CohenD", "Rcy", "P.Rcy", "P.X.noC", "P.X.C", "P.X.C.I",
               "P.IntXC", "d.P.X.noC", "d.P.X.C")
  
  # Set a dataframe for storing data
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
    
    x0 <- as.data.frame(mvrnorm(Ncell0, c(0, 0), matrix(c(SDY0^2, Ccy0, Ccy0, 
                SDC0^2), nrow = 2, ncol = 2, byrow = TRUE), empirical = FALSE))
    colnames(x0) <- c("y", "c")
    
    x1 <- as.data.frame(mvrnorm(Ncell1, c(diff.mu, 0), matrix(c(SDY1^2, Ccy1, 
          Ccy1, SDC1^2), nrow = 2, ncol = 2, byrow = TRUE), empirical = FALSE))
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
    temp[i, "P.X.C.I"] <<- coef(summary(lm(y ~ x + c + x * c, data = dat)))[2, 4]
    temp[i, "P.IntXC"] <<- coef(summary(lm(y ~ x + c + x * c, data = dat)))[4, 4]
    temp[i, "d.P.X.noC"] <<- coef(summary(lm(y ~ x, data = dat)))[2, 1]
    temp[i, "d.P.X.C"] <<- coef(summary(lm(y ~ x + c, data = dat)))[2, 1]
  })
  
  ######## Baseline Strategy #########
  Baseline <- sum(temp$P.X.noC <= alpha)/dim(temp)[1]
  
  ######## Covariate Strategy ########
  Covariate <- sum(unlist(lapply(1:dim(temp)[1], function(i) {
    min(temp[i, c("P.X.noC", "P.X.C")])
  })) <= alpha)/dim(temp)[1]
  
  ######## Merge results #############
  Strategy <- c("Baseline", "Covariate")
  if (diff.mu == 0) {
    TypeIError <- c(Baseline, Covariate)
    Results <- data.frame(Strategy, TypeIError)
  } else {
    Power <- c(Baseline, Covariate)
    Results <- data.frame(Strategy, Power)
  }
  return(list(Results, temp))
}

############################# DATA GENERATION ###############################

# We ran 10,000 iterations for each condition and simulated data for a large
# number of conditions. Running the full syntax, depending on the speed of your
# CPU, would take many hours. To illustrate our method, the results below are
# generated for N = 160, d = .4, rho = .4. This will reproduce the corresponding
# values reported in the manuscript. You can reproduce any value reported in
# Table S1 and S2 (our full set of results) by changing corresponding parameter
# inputs.

##### Power #####
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
# seed: seed value you wish to use (default seed = 42, as used in the manuscript)

# See generated data
View(N160_d0.4_rho0.4[2])

##### Type I Error Rate #####
N160_d0_rho0.4 <- SimData(Nsims = 10000, Ncell0 = 80, Ncell1 = 80, diff.mu = 0, 
                          rcy = 0.4, SDY0 = 1, SDY1 = 1, SDC0 = 1, SDC1 = 1, 
                          seed = 42)
# We simply set d = 0 to obtain Type I error rate

# See generated data
View(N160_d0_rho0.4[2])

############################# RESULTS ###############################

##### Power #####
N160_d0.4_rho0.4[1]
# You should get 0.7131 for Baseline and 0.8111 for Covariate
# Power boost = 0.098

##### Type I Error Rate #####
N160_d0_rho0.4[1]
# You should get 0.0509 for Baseline and 0.0675 for Covariate
# Inflation = 0.0166

##### Number of additional participants needed for comparable power boost #####
# Reported in Table 1
# Calculations conducted in G*Power using t-test/Means: Difference
# between two independent means (two groups)/Type of Power Analysis: A priori

##### Accuracy and Precision of the Two Strategies ####
# Running the syntax below will reproduce the corresponding values reported
# in Table S3.

dat_N160_d0.4_rho0.4 <- as.data.frame(N160_d0.4_rho0.4[2])
dat_N160_d0_rho0.4 <- as.data.frame(N160_d0_rho0.4[2])

# Mean of the sample effect size estimates by the Baseline Strategy
round(mean(dat_N160_d0.4_rho0.4$d.P.X.noC), 3)

# Mean of the sample effect size estimates by the Covariate Strategy
round(mean(dat_N160_d0.4_rho0.4$d.P.X.C), 3)

# SD of the sample effect size estimates by the Baseline Strategy
round(sd(dat_N160_d0.4_rho0.4$d.P.X.noC), 3)

# SD of the sample effect size estimates by the Covariate Strategy
round(sd(dat_N160_d0.4_rho0.4$d.P.X.C), 3)


######################## ADDITIONAL SIMULATIONS ##########################

# The expanded function below allows us to simulate strategies in addition to the
# two strategies described above. These additional strategies are described in
# Section IV of the Supplemental Materials.

# Kitchen-Sink Approach to Covariate Use: 
# Test the effect of X on Y, 
# the effect of X on Y controlling for C,
# the interactive effect of X and C on Y, 
# and the effect of X on Y controlling for both C
# and the interaction between C and Y.

# Qualified Covariate Strategy (A): 
# Test the effect of X on Y, 
# and the effect of X on Y controlling for C only if r(C, Y) >= .3 in the sample.

# Qualified Covariate Strategy (B):
# Test the effect of X on Y, 
# and the effect of X on Y controlling for C only if r(C, Y) is 
# significantly different from 0 in the sample.

# Refer to the original function SimData() above for annotations on our procedure.

SuppSimData <- function(Nsims, Ncell0, Ncell1 = Ncell0, diff.mu, rcy, SDY0, 
                        SDY1 = SDY0, SDC0, SDC1 = SDC0, alpha = 0.05, seed = 42) {
  varlist <- c("N", "N0", "N1", "SDY0", "SDY1", "SDC0", "SDC1", "Mudiff", "Mdiff", 
               "CohenD", "Rcy", "P.Rcy", "P.X.noC", "P.X.C", "P.X.C.I", "P.IntXC",
               "d.P.X.noC", "d.P.X.C")
  temp <- as.data.frame(matrix(NA, nrow = Nsims, ncol = length(varlist)))
  colnames(temp) <- varlist
  Ccy0 <- SDY0 * SDC0 * rcy
  Ccy1 <- SDY1 * SDC1 * rcy
  set.seed(seed)
  lapply(1:Nsims, function(i) {
    x0 <- as.data.frame(mvrnorm(Ncell0, c(0, 0), matrix(c(SDY0^2, Ccy0, Ccy0, 
             SDC0^2), nrow = 2, ncol = 2, byrow = TRUE), empirical = FALSE))
    colnames(x0) <- c("y", "c")
    x1 <- as.data.frame(mvrnorm(Ncell1, c(diff.mu, 0), matrix(c(SDY1^2, Ccy1, 
            Ccy1, SDC1^2), nrow = 2, ncol = 2, byrow = TRUE), empirical = FALSE))
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
    temp[i, "P.X.C.I"] <<- coef(summary(lm(y ~ x + c + x * c, data = dat)))[2, 4]
    temp[i, "P.IntXC"] <<- coef(summary(lm(y ~ x + c + x * c, data = dat)))[4, 4]
    temp[i, "d.P.X.noC"] <<- coef(summary(lm(y ~ x, data = dat)))[2, 1]
    temp[i, "d.P.X.C"] <<- coef(summary(lm(y ~ x + c, data = dat)))[2, 1]
  })
  
  ######## Baseline Strategy #########
  Baseline <- sum(temp$P.X.noC <= alpha)/dim(temp)[1]
  
  ######## Covariate Strategy ########
  Covariate <- sum(unlist(lapply(1:dim(temp)[1], function(i) {
    min(temp[i, c("P.X.noC", "P.X.C")])
  })) <= alpha)/dim(temp)[1]
  
  ######## Kitchen-Sink Strategy #########
  Kitchen_Sink <- sum(unlist(lapply(1:dim(temp)[1], function(i) {
    if (temp[i, "P.IntXC"] <= alpha) {
      temp[i, "P.IntXC"]
    } else {
      min(temp[i, c("P.X.noC", "P.X.C", "P.X.C.I")])
    }
  })) <= alpha)/dim(temp)[1]
  
  ######## Qualified Covariate Strategy (a) ####
  Qualified.A <- sum(unlist(lapply(1:dim(temp)[1], function(i) {
    if (temp[i, "Rcy"] >= 0.3) {
      min(temp[i, c("P.X.noC", "P.X.C")])
    } else {
      temp[i, "P.X.noC"]
    }
  })) <= alpha)/dim(temp)[1]
  
  ######## Qualified Covariate Strategy (b) ####
  Qualified.B <- sum(unlist(lapply(1:dim(temp)[1], function(i) {
    if (temp[i, "P.Rcy"] <= 0.05) {
      min(temp[i, c("P.X.noC", "P.X.C")])
    } else {
      temp[i, "P.X.noC"]
    }
  })) <= alpha)/dim(temp)[1]
  
  ######## Merge results ########
  Strategy <- c("Baseline", "Covariate", "Kitchen_Sink", 
                "Qualified.A", "Qualified.B")
  if (diff.mu == 0) {
    TypeIError <- c(Baseline, Covariate, Kitchen_Sink, Qualified.A, Qualified.B)
    Results <- data.frame(Strategy, TypeIError)
  } else {
    Power <- c(Baseline, Covariate, Kitchen_Sink, Qualified.A, Qualified.B)
    Results <- data.frame(Strategy, Power)
  }
  return(list(Results, temp))
}


# Data Generation for N = 160, d = 0.4, rho = .4
Supp_N160_d0.4_rho0.4 <- SuppSimData(10000, 80, 80, 0.4, 0.4, 1, 1, 1, 1)
Supp_N160_d0_rho0.4 <- SuppSimData(10000, 80, 80, 0, 0.4, 1, 1, 1, 1)

# Data Generation for N = 160, d = 0.4, rho = .2
Supp_N160_d0.4_rho0.2 <- SuppSimData(10000, 80, 80, 0.4, 0.2, 1, 1, 1, 1)
Supp_N160_d0_rho0.2 <- SuppSimData(10000, 80, 80, 0, 0.2, 1, 1, 1, 1)

#### Allowing C to Interact with Y ####

# Power
Supp_N160_d0.4_rho0.4[1]
# You should get 0.8111 for the Covariate Strategy, and
# 0.8223 for the Kitchen Sink approach (change = 0.0112)

# Type I error rate
Supp_N160_d0_rho0.4[1]
# You should get 0.0675 for the Covariate Strategy, and
# 0.1182 for the Kitchen Sink approach (change = 0.0507)

#### Using the Observed Correlation in the Sample ####

# We explored two variations of the Covariate Strategy that use the observed
# correlation between C and Y in the sample:
# Qualified Covariate Strategy (A), and Qualified Covariate Strategy (B). 
# See the function SuppSimData() above for details.

# When rho < .3: Power
Supp_N160_d0.4_rho0.2[1]
# You should get 0.7166 for Qualified (A), and 0.7413 for Qualified (B).
# Qualified (A) performs similarly to Baseline.
# Qualified (B) performs similarly to Covariate.

# When rho < .3: Type I error rate
Supp_N160_d0_rho0.2[1]
# You should get 0.0499 for Qualified (A), and 0.0575 for Qualified (B).
# Qualified (A) performs similarly to Baseline. 
# Qualified (B) performs similarly to Covariate.

# When rho > .3: Power
Supp_N160_d0.4_rho0.4[1]
# You should get 0.8033 for Qualified (A), and 0.8110 for Qualified (B).
# Both Qualified (A) and Qualified (B) perform similarly to Covariate.

# When rho > .3: Type I error rate
Supp_N160_d0_rho0.4[1]
# You should get 0.0668 for Qualified (A), and 0.0675 for Qualified (B).
# Both Qualified (A) and Qualified (B) perform similarly to Covariate.

# As you can see, neither Qualified (A) nor Qualified (B) offers a more favorable
# tradeoff than the Baseline or the Covariate.
# Using the observed correlation in the sample does not offer protection against
# Type I error inflation. It can undercut power boost.

#### Adjusting Alpha Level ####

dat_Supp_N160_d0_rho0.4 <- as.data.frame(Supp_N160_d0_rho0.4[2])

# Find the alpha level that controls Type I error rate at .05
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
# As expected, Type I error rate is exactly .05 now for the Covariate Strategy.

# Power
adj1.N160_d0.4_rho0.4[1]
# You should get 0.7716 for the Covariate Strategy.
# Note that since we used the alpha adjustment,
# the baseline is not the baseline output produced here;
# it should be the one prior to adjustment (0.7131).
# Power boost = 0.7716 - 0.7131 = 0.0585.

# Alternatively, we can directly use a smaller alpha value (e.g., alpha = .033)
adj2.N160_d0.4_rho0.4 <- SimData(10000, 80, 80, 0.4, 0.4, 1, 1, 1, 1, 
                                 alpha = 0.033)
adj2.N160_d0_rho0.4 <- SimData(10000, 80, 80, 0, 0.4, 1, 1, 1, 1, 
                               alpha = 0.033)

# Type I Error Rate
adj2.N160_d0_rho0.4[1]
# Type I error rate is now 0.0457.

# Power
adj2.N160_d0.4_rho0.4[1]
# You should get 0.7546 for the Covariate Strategy.
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
# You should get 0.6663 for the Covariate Strategy.
# We actually have a power loss of 0.7131 - 0.6663 = 0.0468 now.
