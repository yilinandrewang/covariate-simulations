# Covariate Simulations of Type I and II errors
# Written by Y. Andre Wang
# Modified by Joseph E. Gonzales
# Last Update: 05/03/16

############################# INTRO ###############################

# The purpose of these simulations is to compare Type I and II errors that result from using different analysis strategies
# We specifically examine the effect of analysis strategies that flexibly include covariates
# We consider four strategies that researchers might use...
#
# Strategy 1. Only look at the effect of IV (x) on DV (y)
#
# Strategy 2. Look at the effect of x on y, 
#               and the effect of x on y after controlling for covariate (c)
#
# Strategy 3. Look at the effect of x on y, 
#               and the effect of x on y after controlling for c, 
#               but only if the correlation between y and c is significantly greater than 0
#
# Strategy 4. Look at the effect of x on y,
#               and the effect of x on y after controlling for c,
#               and the main effect of x on y by treating both x and c as IVs,
#               and the interaction of x and c on y by treating both a and c as IVs.
#
# ...and compare the Type I and II errors of using these analysis strategies...
# ...if researchers only report the most favorable result (i.e., smallest p value observed)


###################### SIMULATIONS SETUP ##########################

####### The following packages are used for the simulations
library(MASS);library(MBESS)

####### Function for data simulation ####### 
PowerSimData <- function(Nsims, Ncell0, Ncell1 = Ncell0, diff.mu, rcy, 
				 SDY0, SDY1 = SDY0, SDC0, SDC1 = SDC0, alpha=0.05, seed=42) { 
  
  ###### The function above specifies what inputs we can use for simulations
  # Nsims: Number of Simulations,
  # Ncell0: Number of participants in Condition 0
  # Ncell1: Number of participants in Condition 1 (by default Ncell1 = Ncell0)
  # diff.mu: Effect size of x on y in the population (mean group difference)
  # rcy: Correlation of c and y in the population
  # SDY0: Standard deviation of y in Condition 0 in the population
  # SDY1: Standard deviation of y in Condition 1 in the population (by default SDY1 = SDY0)
  # SDC0: Standard deviation of c in Condition 0 in the population
  # SDC1: Standard deviation of c in Condition 1 in the population (by default SDC1 = SDC0)
  # alpha: Alpha level (default = 0.05)
  # seed: seed value you wish to use; default is no seed

####### Data frame for simulated data #######
  
  # List of variables
  varlist <- c("N","N0","N1","SDY0","SDY1","SDC0","SDC1","Mudiff","Mdiff",
	"CohenD","Rcy","P.Rcy","P.X.noC","P.X.C","P.X.C.I","P.IntXC")
  
  # Set dataframe
  temp <- as.data.frame(matrix(NA, nrow = Nsims, ncol = length(varlist)))
  colnames(temp) <- varlist

  # Covariance of Y and C by x condition (x = 0 vs. 1)
  Ccy0 <- SDY0*SDC0*rcy
  Ccy1 <- SDY1*SDC1*rcy

####### Data simulation and p-value collection loop #######
  set.seed(seed)
  lapply(1:Nsims, function(i){ 

    x0 <- as.data.frame(mvrnorm(Ncell0, c(0,0), 
		matrix(c(SDY0^2, Ccy0, Ccy0, SDC0^2),nrow=2,ncol=2,byrow=TRUE), 
		empirical = FALSE))
    colnames(x0) <- c("y","c")

    x1 <- as.data.frame(mvrnorm(Ncell1, c(diff.mu,0), 
		matrix(c(SDY1^2, Ccy1, Ccy1, SDC1^2),nrow=2,ncol=2,byrow=TRUE), 
		empirical = FALSE))
    colnames(x1) <- c("y","c")
    dat <- rbind(x0,x1)
    dat$x <- c(rep(0,Ncell0),rep(1,Ncell1))

    temp[i,"Rcy"] <<- cor(dat$y, dat$c)
    temp[i,"P.Rcy"] <<- cor.test(dat$y, dat$c)$p.value
    temp[i,"P.X.noC"] <<- coef(summary(lm(y ~ x, data = dat)))[2,4]
    temp[i,"P.X.C"] <<- coef(summary(lm(y ~ x + c, data = dat)))[2,4]
    temp[i,"P.X.C.I"] <<- coef(summary(lm(y ~ x + c + x*c, data = dat)))[2,4]
    temp[i,"P.IntXC"] <<- coef(summary(lm(y ~ x + c + x*c, data = dat)))[4,4]    
    temp[i,"N"] <<- dim(dat)[1]
    temp[i,"N0"] <<- Ncell0
    temp[i,"N1"] <<- Ncell1
    temp[i,"SDY0"] <<- ysd0 <- sd(dat[which(dat$x == 0),"y"])
    temp[i,"SDY1"] <<- ysd1 <- sd(dat[which(dat$x == 1),"y"])
    temp[i,"SDC0"] <<- csd0 <- sd(dat[which(dat$x == 0),"c"])
    temp[i,"SDC1"] <<- csd1 <- sd(dat[which(dat$x == 1),"c"])
    temp[i,"Mudiff"] <<- diff.mu
    temp[i,"Mdiff"] <<- diff.m <- abs(mean(dat[which(dat$x == 1),"y"]) - 
						  mean(dat[which(dat$x == 0),"y"]))
	spool <- sqrt(((ysd0^2)*(Ncell0-1) + (ysd1^2)*(Ncell1-1))/(Ncell0 + Ncell1 - 2))
    temp[i,"CohenD"] <<- diff.m/spool
	})

######## Strategy 1 #########
    st1 <- sum(temp$P.X.noC <= alpha)/dim(temp)[1]

######## Strategy 2 ######## 
    st2 <- sum(unlist(lapply(1:dim(temp)[1], function(i){
	min(temp[i,c("P.X.noC","P.X.C")])})) <= alpha)/dim(temp)[1]

######## Strategy 3 ########
    st3 <- sum(unlist(lapply(1:dim(temp)[1], function(i){
	if(temp[i,"P.Rcy"] <= alpha){min(temp[i,c("P.X.noC","P.X.C")])}else{temp[i,"P.X.noC"]}
	})) <= alpha)/dim(temp)[1]

######## Strategy 4 ########
    st4 <- sum(unlist(lapply(1:dim(temp)[1], function(i){
	if(temp[i,"P.IntXC"] <= alpha){temp[i,"P.IntXC"]}else{
		min(temp[i,c("P.X.noC","P.X.C","P.X.C.I")])
		}})) <= alpha)/dim(temp)[1]
  
######## Strategy 6 #######
    st6 <- sum(unlist(lapply(1:dim(temp)[1], function(i){
    if(temp[i,"Rcy"] >= 0.3){min(temp[i,c("P.X.noC","P.X.C")])}else{temp[i,"P.X.noC"]}
  })) <= alpha)/dim(temp)[1]

######## Merge Results ########
	Strategy <- c("ST1","ST2","ST3","ST4","ST6")
    if(diff.mu == 0){
	TypeIError <- c(st1, st2, st3, st4, st6)
	Results <- data.frame(Strategy,TypeIError)
		}else{
	Power <- c(st1, st2, st3, st4, st6)
	Results <- data.frame(Strategy,Power)
		}
return(list(Results,temp))
}

######### Alpha Level Adjustment ####

n100_d0_rcy0.75 <- read.csv("n100_d0_rcy0.75.csv")


ff<-list(n20_d0_rcy0,n20_d0_rcy0.25,n20_d0_rcy0.5,n20_d0_rcy0.75,
         n50_d0_rcy0,n50_d0_rcy0.25,n50_d0_rcy0.5,n50_d0_rcy0.75,
         n100_d0_rcy0,n100_d0_rcy0.25,n100_d0_rcy0.5,n100_d0_rcy0.75)

lapply(ff, function(x){   
  # ST1
  x$ST1.p <- pmin(x[,13])
  alpha.ST1 <- sort(x$ST1.p,partial=dim(x)[1]*0.05)[dim(x)[1]*0.05]
  
  # ST2
  x$ST2.p <- pmin(x[,13], x[,14])
  alpha.ST2 <- sort(x$ST2.p,partial=dim(x)[1]*0.05)[dim(x)[1]*0.05]
  
  # ST3
  x$ST3.p <- unlist(lapply(1:dim(x)[1], function(i){
    if(x[i,"P.Rcy"] <= 0.05){min(x[i,c("P.X.noC","P.X.C")])}
    else{x[i,"P.X.noC"]}}))
  alpha.ST3 <- sort(x$ST3.p,partial=dim(x)[1]*0.05)[dim(x)[1]*0.05]
  
  # ST4
  x$ST4.p <- unlist(lapply(1:dim(x)[1], function(i){
    if(x[i,"P.IntXC"] <= 0.05){x[i,"P.IntXC"]}
    else{min(x[i,c("P.X.noC","P.X.C","P.X.C.I")])}}))
  alpha.ST4 <- sort(x$ST4.p,partial=dim(x)[1]*0.05)[dim(x)[1]*0.05]
  
  # ST6
  x$ST6.p <- unlist(lapply(1:dim(x)[1], function(i){
    if(x[i,"Rcy"] >= 0.3){min(x[i,c("P.X.noC","P.X.C")])}else{x[i,"P.X.noC"]}}))
  alpha.ST6 <- sort(x$ST6.p,partial=dim(x)[1]*0.05)[dim(x)[1]*0.05]
  
  # Summary
  AlphaAdj <- c(alpha.ST1, alpha.ST2, alpha.ST3, alpha.ST4, alpha.ST6)
  AlphaAdj
}
)

######### Power (diff.mu > 0) ####

# d = 0.2
n20_d0.2_rcy0[[1]];n20_d0.2_rcy0.25[[1]];n20_d0.2_rcy0.5[[1]];n20_d0.2_rcy0.75[[1]]
n50_d0.2_rcy0[[1]];n50_d0.2_rcy0.25[[1]];n50_d0.2_rcy0.5[[1]];n50_d0.2_rcy0.75[[1]]
n100_d0.2_rcy0[[1]];n100_d0.2_rcy0.25[[1]];n100_d0.2_rcy0.5[[1]];n100_d0.2_rcy0.75[[1]]

# d = 0.3
n20_d0.3_rcy0[[1]];n20_d0.3_rcy0.25[[1]];n20_d0.3_rcy0.5[[1]];n20_d0.3_rcy0.75[[1]]
n50_d0.3_rcy0[[1]];n50_d0.3_rcy0.25[[1]];n50_d0.3_rcy0.5[[1]];n50_d0.3_rcy0.75[[1]]
n100_d0.3_rcy0[[1]];n100_d0.3_rcy0.25[[1]];n100_d0.3_rcy0.5[[1]];n100_d0.3_rcy0.75[[1]]

# d = 0.4
n20_d0.4_rcy0[[1]];n20_d0.4_rcy0.25[[1]];n20_d0.4_rcy0.5[[1]];n20_d0.4_rcy0.75[[1]]
n50_d0.4_rcy0[[1]];n50_d0.4_rcy0.25[[1]];n50_d0.4_rcy0.5[[1]];n50_d0.4_rcy0.75[[1]]
n100_d0.4_rcy0[[1]];n100_d0.4_rcy0.25[[1]];n100_d0.4_rcy0.5[[1]];n100_d0.4_rcy0.75[[1]]

# d = 0.5
n20_d0.5_rcy0[[1]];n20_d0.5_rcy0.25[[1]];n20_d0.5_rcy0.5[[1]];n20_d0.5_rcy0.75[[1]]
n50_d0.5_rcy0[[1]];n50_d0.5_rcy0.25[[1]];n50_d0.5_rcy0.5[[1]];n50_d0.5_rcy0.75[[1]]
n100_d0.5_rcy0[[1]];n100_d0.5_rcy0.25[[1]];n100_d0.5_rcy0.5[[1]];n100_d0.5_rcy0.75[[1]]

# d = 0.8
n20_d0.8_rcy0[[1]];n20_d0.8_rcy0.25[[1]];n20_d0.8_rcy0.5[[1]];n20_d0.8_rcy0.75[[1]]
n50_d0.8_rcy0[[1]];n50_d0.8_rcy0.25[[1]];n50_d0.8_rcy0.5[[1]];n50_d0.8_rcy0.75[[1]]
n100_d0.8_rcy0[[1]];n100_d0.8_rcy0.25[[1]];n100_d0.8_rcy0.5[[1]];n100_d0.8_rcy0.75[[1]]

# Power after alpha adjustment ####

# n = 50 ####

n50_d0.8_rcy0.25 <- read.csv("n50_d0.8_rcy0.25.csv")
sum(unlist(lapply(1:dim(n50_d0.8_rcy0.25)[1], function(i){
  min(n50_d0.8_rcy0.25[i,c("P.X.noC","P.X.C")])})) <= 0.033)/dim(n50_d0.8_rcy0.25)[1]

n50_d0_rcy0 <- read.csv("n50_d0_rcy0.csv")
sum(unlist(lapply(1:dim(n50_d0_rcy0)[1], function(i){
  min(n50_d0_rcy0[i,c("P.X.noC","P.X.C")])})) <= 0.033)/dim(n50_d0_rcy0)[1]

n50_d0.2_rcy0 <- read.csv("n50_d0.2_rcy0.csv")
sum(unlist(lapply(1:dim(n50_d0.2_rcy0)[1], function(i){
  min(n50_d0.2_rcy0[i,c("P.X.noC","P.X.C")])})) <= 0.033)/dim(n50_d0.2_rcy0)[1]

n50_d0.3_rcy0 <- read.csv("n50_d0.3_rcy0.csv")
sum(unlist(lapply(1:dim(n50_d0.3_rcy0)[1], function(i){
  min(n50_d0.3_rcy0[i,c("P.X.noC","P.X.C")])})) <= 0.033)/dim(n50_d0.3_rcy0)[1]

n50_d0.4_rcy0 <- read.csv("n50_d0.4_rcy0.csv")
sum(unlist(lapply(1:dim(n50_d0.4_rcy0)[1], function(i){
  min(n50_d0.4_rcy0[i,c("P.X.noC","P.X.C")])})) <= 0.033)/dim(n50_d0.4_rcy0)[1]

n50_d0.5_rcy0 <- read.csv("n50_d0.5_rcy0.csv")
sum(unlist(lapply(1:dim(n50_d0.5_rcy0)[1], function(i){
  min(n50_d0.5_rcy0[i,c("P.X.noC","P.X.C")])})) <= 0.033)/dim(n50_d0.5_rcy0)[1]

n50_d0.8_rcy0 <- read.csv("n50_d0.8_rcy0.csv")
sum(unlist(lapply(1:dim(n50_d0.8_rcy0)[1], function(i){
  min(n50_d0.8_rcy0[i,c("P.X.noC","P.X.C")])})) <= 0.033)/dim(n50_d0.8_rcy0)[1]

n50_d0_rcy0.75 <- read.csv("n50_d0_rcy0.75.csv")
sum(unlist(lapply(1:dim(n50_d0_rcy0.75)[1], function(i){
  min(n50_d0_rcy0.75[i,c("P.X.noC","P.X.C")])})) <= 0.033)/dim(n50_d0_rcy0.75)[1]

n50_d0.2_rcy0.75 <- read.csv("n50_d0.2_rcy0.75.csv")
sum(unlist(lapply(1:dim(n50_d0.2_rcy0.75)[1], function(i){
  min(n50_d0.2_rcy0.75[i,c("P.X.noC","P.X.C")])})) <= 0.033)/dim(n50_d0.2_rcy0.75)[1]

n50_d0.3_rcy0.75 <- read.csv("n50_d0.3_rcy0.75.csv")
sum(unlist(lapply(1:dim(n50_d0.3_rcy0.75)[1], function(i){
  min(n50_d0.3_rcy0.75[i,c("P.X.noC","P.X.C")])})) <= 0.033)/dim(n50_d0.3_rcy0.75)[1]

n50_d0.4_rcy0.75 <- read.csv("n50_d0.4_rcy0.75.csv")
sum(unlist(lapply(1:dim(n50_d0.4_rcy0.75)[1], function(i){
  min(n50_d0.4_rcy0.75[i,c("P.X.noC","P.X.C")])})) <= 0.033)/dim(n50_d0.4_rcy0.75)[1]

n50_d0.5_rcy0.75 <- read.csv("n50_d0.5_rcy0.75.csv")
sum(unlist(lapply(1:dim(n50_d0.5_rcy0.75)[1], function(i){
  min(n50_d0.5_rcy0.75[i,c("P.X.noC","P.X.C")])})) <= 0.033)/dim(n50_d0.5_rcy0.75)[1]

n50_d0.8_rcy0.75 <- read.csv("n50_d0.8_rcy0.75.csv")
sum(unlist(lapply(1:dim(n50_d0.8_rcy0.75)[1], function(i){
  min(n50_d0.8_rcy0.75[i,c("P.X.noC","P.X.C")])})) <= 0.033)/dim(n50_d0.8_rcy0.75)[1]

# n = 20 ####

# rcy = 0.75 ####

n20_d0_rcy0.75 <- read.csv("n20_d0_rcy0.75.csv")
sum(unlist(lapply(1:dim(n20_d0_rcy0.75)[1], function(i){
  min(n20_d0_rcy0.75[i,c("P.X.noC","P.X.C")])})) <= 0.033)/dim(n20_d0_rcy0.75)[1]

n20_d0.2_rcy0.75 <- read.csv("n20_d0.2_rcy0.75.csv")
sum(unlist(lapply(1:dim(n20_d0.2_rcy0.75)[1], function(i){
  min(n20_d0.2_rcy0.75[i,c("P.X.noC","P.X.C")])})) <= 0.033)/dim(n20_d0.2_rcy0.75)[1]

n20_d0.3_rcy0.75 <- read.csv("n20_d0.3_rcy0.75.csv")
sum(unlist(lapply(1:dim(n20_d0.3_rcy0.75)[1], function(i){
  min(n20_d0.3_rcy0.75[i,c("P.X.noC","P.X.C")])})) <= 0.033)/dim(n20_d0.3_rcy0.75)[1]

n20_d0.4_rcy0.75 <- read.csv("n20_d0.4_rcy0.75.csv")
sum(unlist(lapply(1:dim(n20_d0.4_rcy0.75)[1], function(i){
  min(n20_d0.4_rcy0.75[i,c("P.X.noC","P.X.C")])})) <= 0.033)/dim(n20_d0.4_rcy0.75)[1]

n20_d0.5_rcy0.75 <- read.csv("n20_d0.5_rcy0.75.csv")
sum(unlist(lapply(1:dim(n20_d0.5_rcy0.75)[1], function(i){
  min(n20_d0.5_rcy0.75[i,c("P.X.noC","P.X.C")])})) <= 0.033)/dim(n20_d0.5_rcy0.75)[1]

n20_d0.8_rcy0.75 <- read.csv("n20_d0.8_rcy0.75.csv")
sum(unlist(lapply(1:dim(n20_d0.8_rcy0.75)[1], function(i){
  min(n20_d0.8_rcy0.75[i,c("P.X.noC","P.X.C")])})) <= 0.033)/dim(n20_d0.8_rcy0.75)[1]

# rcy = 0.5 ####
n20_d0_rcy0.5 <- read.csv("n20_d0_rcy0.5.csv")
sum(unlist(lapply(1:dim(n20_d0_rcy0.5)[1], function(i){
  min(n20_d0_rcy0.5[i,c("P.X.noC","P.X.C")])})) <= 0.033)/dim(n20_d0_rcy0.5)[1]

n20_d0.2_rcy0.5 <- read.csv("n20_d0.2_rcy0.5.csv")
sum(unlist(lapply(1:dim(n20_d0.2_rcy0.5)[1], function(i){
  min(n20_d0.2_rcy0.5[i,c("P.X.noC","P.X.C")])})) <= 0.033)/dim(n20_d0.2_rcy0.5)[1]

n20_d0.3_rcy0.5 <- read.csv("n20_d0.3_rcy0.5.csv")
sum(unlist(lapply(1:dim(n20_d0.3_rcy0.5)[1], function(i){
  min(n20_d0.3_rcy0.5[i,c("P.X.noC","P.X.C")])})) <= 0.033)/dim(n20_d0.3_rcy0.5)[1]

n20_d0.4_rcy0.5 <- read.csv("n20_d0.4_rcy0.5.csv")
sum(unlist(lapply(1:dim(n20_d0.4_rcy0.5)[1], function(i){
  min(n20_d0.4_rcy0.5[i,c("P.X.noC","P.X.C")])})) <= 0.033)/dim(n20_d0.4_rcy0.5)[1]

n20_d0.5_rcy0.5 <- read.csv("n20_d0.5_rcy0.5.csv")
sum(unlist(lapply(1:dim(n20_d0.5_rcy0.5)[1], function(i){
  min(n20_d0.5_rcy0.5[i,c("P.X.noC","P.X.C")])})) <= 0.033)/dim(n20_d0.5_rcy0.5)[1]

n20_d0.8_rcy0.5 <- read.csv("n20_d0.8_rcy0.5.csv")
sum(unlist(lapply(1:dim(n20_d0.8_rcy0.5)[1], function(i){
  min(n20_d0.8_rcy0.5[i,c("P.X.noC","P.X.C")])})) <= 0.033)/dim(n20_d0.8_rcy0.5)[1]

# rcy = 0.25 ####
n20_d0_rcy0.25 <- read.csv("n20_d0_rcy0.25.csv")
sum(unlist(lapply(1:dim(n20_d0_rcy0.25)[1], function(i){
  min(n20_d0_rcy0.25[i,c("P.X.noC","P.X.C")])})) <= 0.033)/dim(n20_d0_rcy0.25)[1]

n20_d0.2_rcy0.25 <- read.csv("n20_d0.2_rcy0.25.csv")
sum(unlist(lapply(1:dim(n20_d0.2_rcy0.25)[1], function(i){
  min(n20_d0.2_rcy0.25[i,c("P.X.noC","P.X.C")])})) <= 0.033)/dim(n20_d0.2_rcy0.25)[1]

n20_d0.3_rcy0.25 <- read.csv("n20_d0.3_rcy0.25.csv")
sum(unlist(lapply(1:dim(n20_d0.3_rcy0.25)[1], function(i){
  min(n20_d0.3_rcy0.25[i,c("P.X.noC","P.X.C")])})) <= 0.033)/dim(n20_d0.3_rcy0.25)[1]

n20_d0.4_rcy0.25 <- read.csv("n20_d0.4_rcy0.25.csv")
sum(unlist(lapply(1:dim(n20_d0.4_rcy0.25)[1], function(i){
  min(n20_d0.4_rcy0.25[i,c("P.X.noC","P.X.C")])})) <= 0.033)/dim(n20_d0.4_rcy0.25)[1]

n20_d0.5_rcy0.25 <- read.csv("n20_d0.5_rcy0.25.csv")
sum(unlist(lapply(1:dim(n20_d0.5_rcy0.25)[1], function(i){
  min(n20_d0.5_rcy0.25[i,c("P.X.noC","P.X.C")])})) <= 0.033)/dim(n20_d0.5_rcy0.25)[1]

n20_d0.8_rcy0.25 <- read.csv("n20_d0.8_rcy0.25.csv")
sum(unlist(lapply(1:dim(n20_d0.8_rcy0.25)[1], function(i){
  min(n20_d0.8_rcy0.25[i,c("P.X.noC","P.X.C")])})) <= 0.033)/dim(n20_d0.8_rcy0.25)[1]


# rcy = 0 ####
n20_d0_rcy0 <- read.csv("n20_d0_rcy0.csv")
sum(unlist(lapply(1:dim(n20_d0_rcy0)[1], function(i){
  min(n20_d0_rcy0[i,c("P.X.noC","P.X.C")])})) <= 0.033)/dim(n20_d0_rcy0)[1]

n20_d0.2_rcy0 <- read.csv("n20_d0.2_rcy0.csv")
sum(unlist(lapply(1:dim(n20_d0.2_rcy0)[1], function(i){
  min(n20_d0.2_rcy0[i,c("P.X.noC","P.X.C")])})) <= 0.033)/dim(n20_d0.2_rcy0)[1]

n20_d0.3_rcy0 <- read.csv("n20_d0.3_rcy0.csv")
sum(unlist(lapply(1:dim(n20_d0.3_rcy0)[1], function(i){
  min(n20_d0.3_rcy0[i,c("P.X.noC","P.X.C")])})) <= 0.033)/dim(n20_d0.3_rcy0)[1]

n20_d0.4_rcy0 <- read.csv("n20_d0.4_rcy0.csv")
sum(unlist(lapply(1:dim(n20_d0.4_rcy0)[1], function(i){
  min(n20_d0.4_rcy0[i,c("P.X.noC","P.X.C")])})) <= 0.033)/dim(n20_d0.4_rcy0)[1]

n20_d0.5_rcy0 <- read.csv("n20_d0.5_rcy0.csv")
sum(unlist(lapply(1:dim(n20_d0.5_rcy0)[1], function(i){
  min(n20_d0.5_rcy0[i,c("P.X.noC","P.X.C")])})) <= 0.033)/dim(n20_d0.5_rcy0)[1]

n20_d0.8_rcy0 <- read.csv("n20_d0.8_rcy0.csv")
sum(unlist(lapply(1:dim(n20_d0.8_rcy0)[1], function(i){
  min(n20_d0.8_rcy0[i,c("P.X.noC","P.X.C")])})) <= 0.033)/dim(n20_d0.8_rcy0)[1]


# n = 100 ####

# rcy = 0.75 ####

n100_d0_rcy0.75 <- read.csv("n100_d0_rcy0.75.csv")
sum(unlist(lapply(1:dim(n100_d0_rcy0.75)[1], function(i){
  min(n100_d0_rcy0.75[i,c("P.X.noC","P.X.C")])})) <= 0.033)/dim(n100_d0_rcy0.75)[1]

n100_d0.2_rcy0.75 <- read.csv("n100_d0.2_rcy0.75.csv")
sum(unlist(lapply(1:dim(n100_d0.2_rcy0.75)[1], function(i){
  min(n100_d0.2_rcy0.75[i,c("P.X.noC","P.X.C")])})) <= 0.033)/dim(n100_d0.2_rcy0.75)[1]

n100_d0.3_rcy0.75 <- read.csv("n100_d0.3_rcy0.75.csv")
sum(unlist(lapply(1:dim(n100_d0.3_rcy0.75)[1], function(i){
  min(n100_d0.3_rcy0.75[i,c("P.X.noC","P.X.C")])})) <= 0.033)/dim(n100_d0.3_rcy0.75)[1]

n100_d0.4_rcy0.75 <- read.csv("n100_d0.4_rcy0.75.csv")
sum(unlist(lapply(1:dim(n100_d0.4_rcy0.75)[1], function(i){
  min(n100_d0.4_rcy0.75[i,c("P.X.noC","P.X.C")])})) <= 0.033)/dim(n100_d0.4_rcy0.75)[1]

n100_d0.5_rcy0.75 <- read.csv("n100_d0.5_rcy0.75.csv")
sum(unlist(lapply(1:dim(n100_d0.5_rcy0.75)[1], function(i){
  min(n100_d0.5_rcy0.75[i,c("P.X.noC","P.X.C")])})) <= 0.033)/dim(n100_d0.5_rcy0.75)[1]

n100_d0.8_rcy0.75 <- read.csv("n100_d0.8_rcy0.75.csv")
sum(unlist(lapply(1:dim(n100_d0.8_rcy0.75)[1], function(i){
  min(n100_d0.8_rcy0.75[i,c("P.X.noC","P.X.C")])})) <= 0.033)/dim(n100_d0.8_rcy0.75)[1]

# rcy = 0.5 ####
n100_d0_rcy0.5 <- read.csv("n100_d0_rcy0.5.csv")
sum(unlist(lapply(1:dim(n100_d0_rcy0.5)[1], function(i){
  min(n100_d0_rcy0.5[i,c("P.X.noC","P.X.C")])})) <= 0.033)/dim(n100_d0_rcy0.5)[1]

n100_d0.2_rcy0.5 <- read.csv("n100_d0.2_rcy0.5.csv")
sum(unlist(lapply(1:dim(n100_d0.2_rcy0.5)[1], function(i){
  min(n100_d0.2_rcy0.5[i,c("P.X.noC","P.X.C")])})) <= 0.033)/dim(n100_d0.2_rcy0.5)[1]

n100_d0.3_rcy0.5 <- read.csv("n100_d0.3_rcy0.5.csv")
sum(unlist(lapply(1:dim(n100_d0.3_rcy0.5)[1], function(i){
  min(n100_d0.3_rcy0.5[i,c("P.X.noC","P.X.C")])})) <= 0.033)/dim(n100_d0.3_rcy0.5)[1]

n100_d0.4_rcy0.5 <- read.csv("n100_d0.4_rcy0.5.csv")
sum(unlist(lapply(1:dim(n100_d0.4_rcy0.5)[1], function(i){
  min(n100_d0.4_rcy0.5[i,c("P.X.noC","P.X.C")])})) <= 0.033)/dim(n100_d0.4_rcy0.5)[1]

n100_d0.5_rcy0.5 <- read.csv("n100_d0.5_rcy0.5.csv")
sum(unlist(lapply(1:dim(n100_d0.5_rcy0.5)[1], function(i){
  min(n100_d0.5_rcy0.5[i,c("P.X.noC","P.X.C")])})) <= 0.033)/dim(n100_d0.5_rcy0.5)[1]

n100_d0.8_rcy0.5 <- read.csv("n100_d0.8_rcy0.5.csv")
sum(unlist(lapply(1:dim(n100_d0.8_rcy0.5)[1], function(i){
  min(n100_d0.8_rcy0.5[i,c("P.X.noC","P.X.C")])})) <= 0.033)/dim(n100_d0.8_rcy0.5)[1]

# rcy = 0.25 ####
n100_d0_rcy0.25 <- read.csv("n100_d0_rcy0.25.csv")
sum(unlist(lapply(1:dim(n100_d0_rcy0.25)[1], function(i){
  min(n100_d0_rcy0.25[i,c("P.X.noC","P.X.C")])})) <= 0.033)/dim(n100_d0_rcy0.25)[1]

n100_d0.2_rcy0.25 <- read.csv("n100_d0.2_rcy0.25.csv")
sum(unlist(lapply(1:dim(n100_d0.2_rcy0.25)[1], function(i){
  min(n100_d0.2_rcy0.25[i,c("P.X.noC","P.X.C")])})) <= 0.033)/dim(n100_d0.2_rcy0.25)[1]

n100_d0.3_rcy0.25 <- read.csv("n100_d0.3_rcy0.25.csv")
sum(unlist(lapply(1:dim(n100_d0.3_rcy0.25)[1], function(i){
  min(n100_d0.3_rcy0.25[i,c("P.X.noC","P.X.C")])})) <= 0.033)/dim(n100_d0.3_rcy0.25)[1]

n100_d0.4_rcy0.25 <- read.csv("n100_d0.4_rcy0.25.csv")
sum(unlist(lapply(1:dim(n100_d0.4_rcy0.25)[1], function(i){
  min(n100_d0.4_rcy0.25[i,c("P.X.noC","P.X.C")])})) <= 0.033)/dim(n100_d0.4_rcy0.25)[1]

n100_d0.5_rcy0.25 <- read.csv("n100_d0.5_rcy0.25.csv")
sum(unlist(lapply(1:dim(n100_d0.5_rcy0.25)[1], function(i){
  min(n100_d0.5_rcy0.25[i,c("P.X.noC","P.X.C")])})) <= 0.033)/dim(n100_d0.5_rcy0.25)[1]

n100_d0.8_rcy0.25 <- read.csv("n100_d0.8_rcy0.25.csv")
sum(unlist(lapply(1:dim(n100_d0.8_rcy0.25)[1], function(i){
  min(n100_d0.8_rcy0.25[i,c("P.X.noC","P.X.C")])})) <= 0.033)/dim(n100_d0.8_rcy0.25)[1]


# rcy = 0 ####
n100_d0_rcy0 <- read.csv("n100_d0_rcy0.csv")
sum(unlist(lapply(1:dim(n100_d0_rcy0)[1], function(i){
  min(n100_d0_rcy0[i,c("P.X.noC","P.X.C")])})) <= 0.033)/dim(n100_d0_rcy0)[1]

n100_d0.2_rcy0 <- read.csv("n100_d0.2_rcy0.csv")
sum(unlist(lapply(1:dim(n100_d0.2_rcy0)[1], function(i){
  min(n100_d0.2_rcy0[i,c("P.X.noC","P.X.C")])})) <= 0.033)/dim(n100_d0.2_rcy0)[1]

n100_d0.3_rcy0 <- read.csv("n100_d0.3_rcy0.csv")
sum(unlist(lapply(1:dim(n100_d0.3_rcy0)[1], function(i){
  min(n100_d0.3_rcy0[i,c("P.X.noC","P.X.C")])})) <= 0.033)/dim(n100_d0.3_rcy0)[1]

n100_d0.4_rcy0 <- read.csv("n100_d0.4_rcy0.csv")
sum(unlist(lapply(1:dim(n100_d0.4_rcy0)[1], function(i){
  min(n100_d0.4_rcy0[i,c("P.X.noC","P.X.C")])})) <= 0.033)/dim(n100_d0.4_rcy0)[1]

n100_d0.5_rcy0 <- read.csv("n100_d0.5_rcy0.csv")
sum(unlist(lapply(1:dim(n100_d0.5_rcy0)[1], function(i){
  min(n100_d0.5_rcy0[i,c("P.X.noC","P.X.C")])})) <= 0.033)/dim(n100_d0.5_rcy0)[1]

n100_d0.8_rcy0 <- read.csv("n100_d0.8_rcy0.csv")
sum(unlist(lapply(1:dim(n100_d0.8_rcy0)[1], function(i){
  min(n100_d0.8_rcy0[i,c("P.X.noC","P.X.C")])})) <= 0.033)/dim(n100_d0.8_rcy0)[1]


# Strategy 6 #####

# n = 50
n50_d0.4_rcy0 <- read.csv("n50_d0.4_rcy0.csv")
sum(unlist(lapply(1:dim(n50_d0.4_rcy0)[1], function(i){
  if(n50_d0.4_rcy0[i,"Rcy"] >= 0.3){min(n50_d0.4_rcy0[i,c("P.X.noC","P.X.C")])}else{n50_d0.4_rcy0[i,"P.X.noC"]}
})) <= 0.05)/dim(n50_d0.4_rcy0)[1]

n50_d0.4_rcy0.25 <- read.csv("n50_d0.4_rcy0.25.csv")
sum(unlist(lapply(1:dim(n50_d0.4_rcy0.25)[1], function(i){
  if(n50_d0.4_rcy0.25[i,"Rcy"] >= 0.3){min(n50_d0.4_rcy0.25[i,c("P.X.noC","P.X.C")])}else{n50_d0.4_rcy0.25[i,"P.X.noC"]}
})) <= 0.05)/dim(n50_d0.4_rcy0.25)[1]

n50_d0.4_rcy0.5 <- read.csv("n50_d0.4_rcy0.5.csv")
sum(unlist(lapply(1:dim(n50_d0.4_rcy0.5)[1], function(i){
  if(n50_d0.4_rcy0.5[i,"Rcy"] >= 0.3){min(n50_d0.4_rcy0.5[i,c("P.X.noC","P.X.C")])}else{n50_d0.4_rcy0.5[i,"P.X.noC"]}
})) <= 0.05)/dim(n50_d0.4_rcy0.5)[1]

n50_d0.4_rcy0.75 <- read.csv("n50_d0.4_rcy0.75.csv")
sum(unlist(lapply(1:dim(n50_d0.4_rcy0.75)[1], function(i){
  if(n50_d0.4_rcy0.75[i,"Rcy"] >= 0.3){min(n50_d0.4_rcy0.75[i,c("P.X.noC","P.X.C")])}else{n50_d0.4_rcy0.75[i,"P.X.noC"]}
})) <= 0.05)/dim(n50_d0.4_rcy0.75)[1]

n50_d0_rcy0 <- read.csv("n50_d0_rcy0.csv")
sum(unlist(lapply(1:dim(n50_d0_rcy0)[1], function(i){
  if(n50_d0_rcy0[i,"Rcy"] >= 0.3){min(n50_d0_rcy0[i,c("P.X.noC","P.X.C")])}else{n50_d0_rcy0[i,"P.X.noC"]}
})) <= 0.05)/dim(n50_d0_rcy0)[1]

n50_d0_rcy0.25 <- read.csv("n50_d0_rcy0.25.csv")
sum(unlist(lapply(1:dim(n50_d0_rcy0.25)[1], function(i){
  if(n50_d0_rcy0.25[i,"Rcy"] >= 0.3){min(n50_d0_rcy0.25[i,c("P.X.noC","P.X.C")])}else{n50_d0_rcy0.25[i,"P.X.noC"]}
})) <= 0.05)/dim(n50_d0_rcy0.25)[1]

n50_d0_rcy0.5 <- read.csv("n50_d0_rcy0.5.csv")
sum(unlist(lapply(1:dim(n50_d0_rcy0.5)[1], function(i){
  if(n50_d0_rcy0.5[i,"Rcy"] >= 0.3){min(n50_d0_rcy0.5[i,c("P.X.noC","P.X.C")])}else{n50_d0_rcy0.5[i,"P.X.noC"]}
})) <= 0.05)/dim(n50_d0_rcy0.5)[1]

n50_d0_rcy0.75 <- read.csv("n50_d0_rcy0.75.csv")
sum(unlist(lapply(1:dim(n50_d0_rcy0.75)[1], function(i){
  if(n50_d0_rcy0.75[i,"Rcy"] >= 0.3){min(n50_d0_rcy0.75[i,c("P.X.noC","P.X.C")])}else{n50_d0_rcy0.75[i,"P.X.noC"]}
})) <= 0.05)/dim(n50_d0_rcy0.75)[1]

n50_d0.2_rcy0 <- read.csv("n50_d0.2_rcy0.csv")
sum(unlist(lapply(1:dim(n50_d0.2_rcy0)[1], function(i){
  if(n50_d0.2_rcy0[i,"Rcy"] >= 0.3){min(n50_d0.2_rcy0[i,c("P.X.noC","P.X.C")])}else{n50_d0.2_rcy0[i,"P.X.noC"]}
})) <= 0.05)/dim(n50_d0.2_rcy0)[1]

n50_d0.2_rcy0.25 <- read.csv("n50_d0.2_rcy0.25.csv")
sum(unlist(lapply(1:dim(n50_d0.2_rcy0.25)[1], function(i){
  if(n50_d0.2_rcy0.25[i,"Rcy"] >= 0.3){min(n50_d0.2_rcy0.25[i,c("P.X.noC","P.X.C")])}else{n50_d0.2_rcy0.25[i,"P.X.noC"]}
})) <= 0.05)/dim(n50_d0.2_rcy0.25)[1]

n50_d0.2_rcy0.5 <- read.csv("n50_d0.2_rcy0.5.csv")
sum(unlist(lapply(1:dim(n50_d0.2_rcy0.5)[1], function(i){
  if(n50_d0.2_rcy0.5[i,"Rcy"] >= 0.3){min(n50_d0.2_rcy0.5[i,c("P.X.noC","P.X.C")])}else{n50_d0.2_rcy0.5[i,"P.X.noC"]}
})) <= 0.05)/dim(n50_d0.2_rcy0.5)[1]

n50_d0.2_rcy0.75 <- read.csv("n50_d0.2_rcy0.75.csv")
sum(unlist(lapply(1:dim(n50_d0.2_rcy0.75)[1], function(i){
  if(n50_d0.2_rcy0.75[i,"Rcy"] >= 0.3){min(n50_d0.2_rcy0.75[i,c("P.X.noC","P.X.C")])}else{n50_d0.2_rcy0.75[i,"P.X.noC"]}
})) <= 0.05)/dim(n50_d0.2_rcy0.75)[1]

n50_d0.8_rcy0 <- read.csv("n50_d0.8_rcy0.csv")
sum(unlist(lapply(1:dim(n50_d0.8_rcy0)[1], function(i){
  if(n50_d0.8_rcy0[i,"Rcy"] >= 0.3){min(n50_d0.8_rcy0[i,c("P.X.noC","P.X.C")])}else{n50_d0.8_rcy0[i,"P.X.noC"]}
})) <= 0.05)/dim(n50_d0.8_rcy0)[1]

n50_d0.8_rcy0.25 <- read.csv("n50_d0.8_rcy0.25.csv")
sum(unlist(lapply(1:dim(n50_d0.8_rcy0.25)[1], function(i){
  if(n50_d0.8_rcy0.25[i,"Rcy"] >= 0.3){min(n50_d0.8_rcy0.25[i,c("P.X.noC","P.X.C")])}else{n50_d0.8_rcy0.25[i,"P.X.noC"]}
})) <= 0.05)/dim(n50_d0.8_rcy0.25)[1]

n50_d0.8_rcy0.5 <- read.csv("n50_d0.8_rcy0.5.csv")
sum(unlist(lapply(1:dim(n50_d0.8_rcy0.5)[1], function(i){
  if(n50_d0.8_rcy0.5[i,"Rcy"] >= 0.3){min(n50_d0.8_rcy0.5[i,c("P.X.noC","P.X.C")])}else{n50_d0.8_rcy0.5[i,"P.X.noC"]}
})) <= 0.05)/dim(n50_d0.8_rcy0.5)[1]

n50_d0.8_rcy0.75 <- read.csv("n50_d0.8_rcy0.75.csv")
sum(unlist(lapply(1:dim(n50_d0.8_rcy0.75)[1], function(i){
  if(n50_d0.8_rcy0.75[i,"Rcy"] >= 0.3){min(n50_d0.8_rcy0.75[i,c("P.X.noC","P.X.C")])}else{n50_d0.8_rcy0.75[i,"P.X.noC"]}
})) <= 0.05)/dim(n50_d0.8_rcy0.75)[1]

n50_d0.5_rcy0 <- read.csv("n50_d0.5_rcy0.csv")
sum(unlist(lapply(1:dim(n50_d0.5_rcy0)[1], function(i){
  if(n50_d0.5_rcy0[i,"Rcy"] >= 0.3){min(n50_d0.5_rcy0[i,c("P.X.noC","P.X.C")])}else{n50_d0.5_rcy0[i,"P.X.noC"]}
})) <= 0.05)/dim(n50_d0.5_rcy0)[1]

n50_d0.5_rcy0.25 <- read.csv("n50_d0.5_rcy0.25.csv")
sum(unlist(lapply(1:dim(n50_d0.5_rcy0.25)[1], function(i){
  if(n50_d0.5_rcy0.25[i,"Rcy"] >= 0.3){min(n50_d0.5_rcy0.25[i,c("P.X.noC","P.X.C")])}else{n50_d0.5_rcy0.25[i,"P.X.noC"]}
})) <= 0.05)/dim(n50_d0.5_rcy0.25)[1]

n50_d0.5_rcy0.5 <- read.csv("n50_d0.5_rcy0.5.csv")
sum(unlist(lapply(1:dim(n50_d0.5_rcy0.5)[1], function(i){
  if(n50_d0.5_rcy0.5[i,"Rcy"] >= 0.3){min(n50_d0.5_rcy0.5[i,c("P.X.noC","P.X.C")])}else{n50_d0.5_rcy0.5[i,"P.X.noC"]}
})) <= 0.05)/dim(n50_d0.5_rcy0.5)[1]

n50_d0.5_rcy0.75 <- read.csv("n50_d0.5_rcy0.75.csv")
sum(unlist(lapply(1:dim(n50_d0.5_rcy0.75)[1], function(i){
  if(n50_d0.5_rcy0.75[i,"Rcy"] >= 0.3){min(n50_d0.5_rcy0.75[i,c("P.X.noC","P.X.C")])}else{n50_d0.5_rcy0.75[i,"P.X.noC"]}
})) <= 0.05)/dim(n50_d0.5_rcy0.75)[1]

n50_d0.3_rcy0 <- read.csv("n50_d0.3_rcy0.csv")
sum(unlist(lapply(1:dim(n50_d0.3_rcy0)[1], function(i){
  if(n50_d0.3_rcy0[i,"Rcy"] >= 0.3){min(n50_d0.3_rcy0[i,c("P.X.noC","P.X.C")])}else{n50_d0.3_rcy0[i,"P.X.noC"]}
})) <= 0.05)/dim(n50_d0.3_rcy0)[1]

n50_d0.3_rcy0.25 <- read.csv("n50_d0.3_rcy0.25.csv")
sum(unlist(lapply(1:dim(n50_d0.3_rcy0.25)[1], function(i){
  if(n50_d0.3_rcy0.25[i,"Rcy"] >= 0.3){min(n50_d0.3_rcy0.25[i,c("P.X.noC","P.X.C")])}else{n50_d0.3_rcy0.25[i,"P.X.noC"]}
})) <= 0.05)/dim(n50_d0.3_rcy0.25)[1]

n50_d0.3_rcy0.5 <- read.csv("n50_d0.3_rcy0.5.csv")
sum(unlist(lapply(1:dim(n50_d0.3_rcy0.5)[1], function(i){
  if(n50_d0.3_rcy0.5[i,"Rcy"] >= 0.3){min(n50_d0.3_rcy0.5[i,c("P.X.noC","P.X.C")])}else{n50_d0.3_rcy0.5[i,"P.X.noC"]}
})) <= 0.05)/dim(n50_d0.3_rcy0.5)[1]

n50_d0.3_rcy0.75 <- read.csv("n50_d0.3_rcy0.75.csv")
sum(unlist(lapply(1:dim(n50_d0.3_rcy0.75)[1], function(i){
  if(n50_d0.3_rcy0.75[i,"Rcy"] >= 0.3){min(n50_d0.3_rcy0.75[i,c("P.X.noC","P.X.C")])}else{n50_d0.3_rcy0.75[i,"P.X.noC"]}
})) <= 0.05)/dim(n50_d0.3_rcy0.75)[1]



# n = 100

n100_d0.4_rcy0 <- read.csv("n100_d0.4_rcy0.csv")
sum(unlist(lapply(1:dim(n100_d0.4_rcy0)[1], function(i){
  if(n100_d0.4_rcy0[i,"Rcy"] >= 0.3){min(n100_d0.4_rcy0[i,c("P.X.noC","P.X.C")])}else{n100_d0.4_rcy0[i,"P.X.noC"]}
})) <= 0.05)/dim(n100_d0.4_rcy0)[1]

n100_d0.4_rcy0.25 <- read.csv("n100_d0.4_rcy0.25.csv")
sum(unlist(lapply(1:dim(n100_d0.4_rcy0.25)[1], function(i){
  if(n100_d0.4_rcy0.25[i,"Rcy"] >= 0.3){min(n100_d0.4_rcy0.25[i,c("P.X.noC","P.X.C")])}else{n100_d0.4_rcy0.25[i,"P.X.noC"]}
})) <= 0.05)/dim(n100_d0.4_rcy0.25)[1]

n100_d0.4_rcy0.5 <- read.csv("n100_d0.4_rcy0.5.csv")
sum(unlist(lapply(1:dim(n100_d0.4_rcy0.5)[1], function(i){
  if(n100_d0.4_rcy0.5[i,"Rcy"] >= 0.3){min(n100_d0.4_rcy0.5[i,c("P.X.noC","P.X.C")])}else{n100_d0.4_rcy0.5[i,"P.X.noC"]}
})) <= 0.05)/dim(n100_d0.4_rcy0.5)[1]

n100_d0.4_rcy0.75 <- read.csv("n100_d0.4_rcy0.75.csv")
sum(unlist(lapply(1:dim(n100_d0.4_rcy0.75)[1], function(i){
  if(n100_d0.4_rcy0.75[i,"Rcy"] >= 0.3){min(n100_d0.4_rcy0.75[i,c("P.X.noC","P.X.C")])}else{n100_d0.4_rcy0.75[i,"P.X.noC"]}
})) <= 0.05)/dim(n100_d0.4_rcy0.75)[1]

n100_d0_rcy0 <- read.csv("n100_d0_rcy0.csv")
sum(unlist(lapply(1:dim(n100_d0_rcy0)[1], function(i){
  if(n100_d0_rcy0[i,"Rcy"] >= 0.3){min(n100_d0_rcy0[i,c("P.X.noC","P.X.C")])}else{n100_d0_rcy0[i,"P.X.noC"]}
})) <= 0.05)/dim(n100_d0_rcy0)[1]

n100_d0_rcy0.25 <- read.csv("n100_d0_rcy0.25.csv")
sum(unlist(lapply(1:dim(n100_d0_rcy0.25)[1], function(i){
  if(n100_d0_rcy0.25[i,"Rcy"] >= 0.3){min(n100_d0_rcy0.25[i,c("P.X.noC","P.X.C")])}else{n100_d0_rcy0.25[i,"P.X.noC"]}
})) <= 0.05)/dim(n100_d0_rcy0.25)[1]

n100_d0_rcy0.5 <- read.csv("n100_d0_rcy0.5.csv")
sum(unlist(lapply(1:dim(n100_d0_rcy0.5)[1], function(i){
  if(n100_d0_rcy0.5[i,"Rcy"] >= 0.3){min(n100_d0_rcy0.5[i,c("P.X.noC","P.X.C")])}else{n100_d0_rcy0.5[i,"P.X.noC"]}
})) <= 0.05)/dim(n100_d0_rcy0.5)[1]

n100_d0_rcy0.75 <- read.csv("n100_d0_rcy0.75.csv")
sum(unlist(lapply(1:dim(n100_d0_rcy0.75)[1], function(i){
  if(n100_d0_rcy0.75[i,"Rcy"] >= 0.3){min(n100_d0_rcy0.75[i,c("P.X.noC","P.X.C")])}else{n100_d0_rcy0.75[i,"P.X.noC"]}
})) <= 0.05)/dim(n100_d0_rcy0.75)[1]

n100_d0.2_rcy0 <- read.csv("n100_d0.2_rcy0.csv")
sum(unlist(lapply(1:dim(n100_d0.2_rcy0)[1], function(i){
  if(n100_d0.2_rcy0[i,"Rcy"] >= 0.3){min(n100_d0.2_rcy0[i,c("P.X.noC","P.X.C")])}else{n100_d0.2_rcy0[i,"P.X.noC"]}
})) <= 0.05)/dim(n100_d0.2_rcy0)[1]

n100_d0.2_rcy0.25 <- read.csv("n100_d0.2_rcy0.25.csv")
sum(unlist(lapply(1:dim(n100_d0.2_rcy0.25)[1], function(i){
  if(n100_d0.2_rcy0.25[i,"Rcy"] >= 0.3){min(n100_d0.2_rcy0.25[i,c("P.X.noC","P.X.C")])}else{n100_d0.2_rcy0.25[i,"P.X.noC"]}
})) <= 0.05)/dim(n100_d0.2_rcy0.25)[1]

n100_d0.2_rcy0.5 <- read.csv("n100_d0.2_rcy0.5.csv")
sum(unlist(lapply(1:dim(n100_d0.2_rcy0.5)[1], function(i){
  if(n100_d0.2_rcy0.5[i,"Rcy"] >= 0.3){min(n100_d0.2_rcy0.5[i,c("P.X.noC","P.X.C")])}else{n100_d0.2_rcy0.5[i,"P.X.noC"]}
})) <= 0.05)/dim(n100_d0.2_rcy0.5)[1]

n100_d0.2_rcy0.75 <- read.csv("n100_d0.2_rcy0.75.csv")
sum(unlist(lapply(1:dim(n100_d0.2_rcy0.75)[1], function(i){
  if(n100_d0.2_rcy0.75[i,"Rcy"] >= 0.3){min(n100_d0.2_rcy0.75[i,c("P.X.noC","P.X.C")])}else{n100_d0.2_rcy0.75[i,"P.X.noC"]}
})) <= 0.05)/dim(n100_d0.2_rcy0.75)[1]

n100_d0.8_rcy0 <- read.csv("n100_d0.8_rcy0.csv")
sum(unlist(lapply(1:dim(n100_d0.8_rcy0)[1], function(i){
  if(n100_d0.8_rcy0[i,"Rcy"] >= 0.3){min(n100_d0.8_rcy0[i,c("P.X.noC","P.X.C")])}else{n100_d0.8_rcy0[i,"P.X.noC"]}
})) <= 0.05)/dim(n100_d0.8_rcy0)[1]

n100_d0.8_rcy0.25 <- read.csv("n100_d0.8_rcy0.25.csv")
sum(unlist(lapply(1:dim(n100_d0.8_rcy0.25)[1], function(i){
  if(n100_d0.8_rcy0.25[i,"Rcy"] >= 0.3){min(n100_d0.8_rcy0.25[i,c("P.X.noC","P.X.C")])}else{n100_d0.8_rcy0.25[i,"P.X.noC"]}
})) <= 0.05)/dim(n100_d0.8_rcy0.25)[1]

n100_d0.8_rcy0.5 <- read.csv("n100_d0.8_rcy0.5.csv")
sum(unlist(lapply(1:dim(n100_d0.8_rcy0.5)[1], function(i){
  if(n100_d0.8_rcy0.5[i,"Rcy"] >= 0.3){min(n100_d0.8_rcy0.5[i,c("P.X.noC","P.X.C")])}else{n100_d0.8_rcy0.5[i,"P.X.noC"]}
})) <= 0.05)/dim(n100_d0.8_rcy0.5)[1]

n100_d0.8_rcy0.75 <- read.csv("n100_d0.8_rcy0.75.csv")
sum(unlist(lapply(1:dim(n100_d0.8_rcy0.75)[1], function(i){
  if(n100_d0.8_rcy0.75[i,"Rcy"] >= 0.3){min(n100_d0.8_rcy0.75[i,c("P.X.noC","P.X.C")])}else{n100_d0.8_rcy0.75[i,"P.X.noC"]}
})) <= 0.05)/dim(n100_d0.8_rcy0.75)[1]

n100_d0.5_rcy0 <- read.csv("n100_d0.5_rcy0.csv")
sum(unlist(lapply(1:dim(n100_d0.5_rcy0)[1], function(i){
  if(n100_d0.5_rcy0[i,"Rcy"] >= 0.3){min(n100_d0.5_rcy0[i,c("P.X.noC","P.X.C")])}else{n100_d0.5_rcy0[i,"P.X.noC"]}
})) <= 0.05)/dim(n100_d0.5_rcy0)[1]

n100_d0.5_rcy0.25 <- read.csv("n100_d0.5_rcy0.25.csv")
sum(unlist(lapply(1:dim(n100_d0.5_rcy0.25)[1], function(i){
  if(n100_d0.5_rcy0.25[i,"Rcy"] >= 0.3){min(n100_d0.5_rcy0.25[i,c("P.X.noC","P.X.C")])}else{n100_d0.5_rcy0.25[i,"P.X.noC"]}
})) <= 0.05)/dim(n100_d0.5_rcy0.25)[1]

n100_d0.5_rcy0.5 <- read.csv("n100_d0.5_rcy0.5.csv")
sum(unlist(lapply(1:dim(n100_d0.5_rcy0.5)[1], function(i){
  if(n100_d0.5_rcy0.5[i,"Rcy"] >= 0.3){min(n100_d0.5_rcy0.5[i,c("P.X.noC","P.X.C")])}else{n100_d0.5_rcy0.5[i,"P.X.noC"]}
})) <= 0.05)/dim(n100_d0.5_rcy0.5)[1]

n100_d0.5_rcy0.75 <- read.csv("n100_d0.5_rcy0.75.csv")
sum(unlist(lapply(1:dim(n100_d0.5_rcy0.75)[1], function(i){
  if(n100_d0.5_rcy0.75[i,"Rcy"] >= 0.3){min(n100_d0.5_rcy0.75[i,c("P.X.noC","P.X.C")])}else{n100_d0.5_rcy0.75[i,"P.X.noC"]}
})) <= 0.05)/dim(n100_d0.5_rcy0.75)[1]

n100_d0.3_rcy0 <- read.csv("n100_d0.3_rcy0.csv")
sum(unlist(lapply(1:dim(n100_d0.3_rcy0)[1], function(i){
  if(n100_d0.3_rcy0[i,"Rcy"] >= 0.3){min(n100_d0.3_rcy0[i,c("P.X.noC","P.X.C")])}else{n100_d0.3_rcy0[i,"P.X.noC"]}
})) <= 0.05)/dim(n100_d0.3_rcy0)[1]

n100_d0.3_rcy0.25 <- read.csv("n100_d0.3_rcy0.25.csv")
sum(unlist(lapply(1:dim(n100_d0.3_rcy0.25)[1], function(i){
  if(n100_d0.3_rcy0.25[i,"Rcy"] >= 0.3){min(n100_d0.3_rcy0.25[i,c("P.X.noC","P.X.C")])}else{n100_d0.3_rcy0.25[i,"P.X.noC"]}
})) <= 0.05)/dim(n100_d0.3_rcy0.25)[1]

n100_d0.3_rcy0.5 <- read.csv("n100_d0.3_rcy0.5.csv")
sum(unlist(lapply(1:dim(n100_d0.3_rcy0.5)[1], function(i){
  if(n100_d0.3_rcy0.5[i,"Rcy"] >= 0.3){min(n100_d0.3_rcy0.5[i,c("P.X.noC","P.X.C")])}else{n100_d0.3_rcy0.5[i,"P.X.noC"]}
})) <= 0.05)/dim(n100_d0.3_rcy0.5)[1]

n100_d0.3_rcy0.75 <- read.csv("n100_d0.3_rcy0.75.csv")
sum(unlist(lapply(1:dim(n100_d0.3_rcy0.75)[1], function(i){
  if(n100_d0.3_rcy0.75[i,"Rcy"] >= 0.3){min(n100_d0.3_rcy0.75[i,c("P.X.noC","P.X.C")])}else{n100_d0.3_rcy0.75[i,"P.X.noC"]}
})) <= 0.05)/dim(n100_d0.3_rcy0.75)[1]



# n = 20

n20_d0.4_rcy0 <- read.csv("n20_d0.4_rcy0.csv")
sum(unlist(lapply(1:dim(n20_d0.4_rcy0)[1], function(i){
  if(n20_d0.4_rcy0[i,"Rcy"] >= 0.3){min(n20_d0.4_rcy0[i,c("P.X.noC","P.X.C")])}else{n20_d0.4_rcy0[i,"P.X.noC"]}
})) <= 0.05)/dim(n20_d0.4_rcy0)[1]

n20_d0.4_rcy0.25 <- read.csv("n20_d0.4_rcy0.25.csv")
sum(unlist(lapply(1:dim(n20_d0.4_rcy0.25)[1], function(i){
  if(n20_d0.4_rcy0.25[i,"Rcy"] >= 0.3){min(n20_d0.4_rcy0.25[i,c("P.X.noC","P.X.C")])}else{n20_d0.4_rcy0.25[i,"P.X.noC"]}
})) <= 0.05)/dim(n20_d0.4_rcy0.25)[1]

n20_d0.4_rcy0.5 <- read.csv("n20_d0.4_rcy0.5.csv")
sum(unlist(lapply(1:dim(n20_d0.4_rcy0.5)[1], function(i){
  if(n20_d0.4_rcy0.5[i,"Rcy"] >= 0.3){min(n20_d0.4_rcy0.5[i,c("P.X.noC","P.X.C")])}else{n20_d0.4_rcy0.5[i,"P.X.noC"]}
})) <= 0.05)/dim(n20_d0.4_rcy0.5)[1]

n20_d0.4_rcy0.75 <- read.csv("n20_d0.4_rcy0.75.csv")
sum(unlist(lapply(1:dim(n20_d0.4_rcy0.75)[1], function(i){
  if(n20_d0.4_rcy0.75[i,"Rcy"] >= 0.3){min(n20_d0.4_rcy0.75[i,c("P.X.noC","P.X.C")])}else{n20_d0.4_rcy0.75[i,"P.X.noC"]}
})) <= 0.05)/dim(n20_d0.4_rcy0.75)[1]

n20_d0_rcy0 <- read.csv("n20_d0_rcy0.csv")
sum(unlist(lapply(1:dim(n20_d0_rcy0)[1], function(i){
  if(n20_d0_rcy0[i,"Rcy"] >= 0.3){min(n20_d0_rcy0[i,c("P.X.noC","P.X.C")])}else{n20_d0_rcy0[i,"P.X.noC"]}
})) <= 0.05)/dim(n20_d0_rcy0)[1]

n20_d0_rcy0.25 <- read.csv("n20_d0_rcy0.25.csv")
sum(unlist(lapply(1:dim(n20_d0_rcy0.25)[1], function(i){
  if(n20_d0_rcy0.25[i,"Rcy"] >= 0.3){min(n20_d0_rcy0.25[i,c("P.X.noC","P.X.C")])}else{n20_d0_rcy0.25[i,"P.X.noC"]}
})) <= 0.05)/dim(n20_d0_rcy0.25)[1]

n20_d0_rcy0.5 <- read.csv("n20_d0_rcy0.5.csv")
sum(unlist(lapply(1:dim(n20_d0_rcy0.5)[1], function(i){
  if(n20_d0_rcy0.5[i,"Rcy"] >= 0.3){min(n20_d0_rcy0.5[i,c("P.X.noC","P.X.C")])}else{n20_d0_rcy0.5[i,"P.X.noC"]}
})) <= 0.05)/dim(n20_d0_rcy0.5)[1]

n20_d0_rcy0.75 <- read.csv("n20_d0_rcy0.75.csv")
sum(unlist(lapply(1:dim(n20_d0_rcy0.75)[1], function(i){
  if(n20_d0_rcy0.75[i,"Rcy"] >= 0.3){min(n20_d0_rcy0.75[i,c("P.X.noC","P.X.C")])}else{n20_d0_rcy0.75[i,"P.X.noC"]}
})) <= 0.05)/dim(n20_d0_rcy0.75)[1]

n20_d0.2_rcy0 <- read.csv("n20_d0.2_rcy0.csv")
sum(unlist(lapply(1:dim(n20_d0.2_rcy0)[1], function(i){
  if(n20_d0.2_rcy0[i,"Rcy"] >= 0.3){min(n20_d0.2_rcy0[i,c("P.X.noC","P.X.C")])}else{n20_d0.2_rcy0[i,"P.X.noC"]}
})) <= 0.05)/dim(n20_d0.2_rcy0)[1]

n20_d0.2_rcy0.25 <- read.csv("n20_d0.2_rcy0.25.csv")
sum(unlist(lapply(1:dim(n20_d0.2_rcy0.25)[1], function(i){
  if(n20_d0.2_rcy0.25[i,"Rcy"] >= 0.3){min(n20_d0.2_rcy0.25[i,c("P.X.noC","P.X.C")])}else{n20_d0.2_rcy0.25[i,"P.X.noC"]}
})) <= 0.05)/dim(n20_d0.2_rcy0.25)[1]

n20_d0.2_rcy0.5 <- read.csv("n20_d0.2_rcy0.5.csv")
sum(unlist(lapply(1:dim(n20_d0.2_rcy0.5)[1], function(i){
  if(n20_d0.2_rcy0.5[i,"Rcy"] >= 0.3){min(n20_d0.2_rcy0.5[i,c("P.X.noC","P.X.C")])}else{n20_d0.2_rcy0.5[i,"P.X.noC"]}
})) <= 0.05)/dim(n20_d0.2_rcy0.5)[1]

n20_d0.2_rcy0.75 <- read.csv("n20_d0.2_rcy0.75.csv")
sum(unlist(lapply(1:dim(n20_d0.2_rcy0.75)[1], function(i){
  if(n20_d0.2_rcy0.75[i,"Rcy"] >= 0.3){min(n20_d0.2_rcy0.75[i,c("P.X.noC","P.X.C")])}else{n20_d0.2_rcy0.75[i,"P.X.noC"]}
})) <= 0.05)/dim(n20_d0.2_rcy0.75)[1]

n20_d0.8_rcy0 <- read.csv("n20_d0.8_rcy0.csv")
sum(unlist(lapply(1:dim(n20_d0.8_rcy0)[1], function(i){
  if(n20_d0.8_rcy0[i,"Rcy"] >= 0.3){min(n20_d0.8_rcy0[i,c("P.X.noC","P.X.C")])}else{n20_d0.8_rcy0[i,"P.X.noC"]}
})) <= 0.05)/dim(n20_d0.8_rcy0)[1]

n20_d0.8_rcy0.25 <- read.csv("n20_d0.8_rcy0.25.csv")
sum(unlist(lapply(1:dim(n20_d0.8_rcy0.25)[1], function(i){
  if(n20_d0.8_rcy0.25[i,"Rcy"] >= 0.3){min(n20_d0.8_rcy0.25[i,c("P.X.noC","P.X.C")])}else{n20_d0.8_rcy0.25[i,"P.X.noC"]}
})) <= 0.05)/dim(n20_d0.8_rcy0.25)[1]

n20_d0.8_rcy0.5 <- read.csv("n20_d0.8_rcy0.5.csv")
sum(unlist(lapply(1:dim(n20_d0.8_rcy0.5)[1], function(i){
  if(n20_d0.8_rcy0.5[i,"Rcy"] >= 0.3){min(n20_d0.8_rcy0.5[i,c("P.X.noC","P.X.C")])}else{n20_d0.8_rcy0.5[i,"P.X.noC"]}
})) <= 0.05)/dim(n20_d0.8_rcy0.5)[1]

n20_d0.8_rcy0.75 <- read.csv("n20_d0.8_rcy0.75.csv")
sum(unlist(lapply(1:dim(n20_d0.8_rcy0.75)[1], function(i){
  if(n20_d0.8_rcy0.75[i,"Rcy"] >= 0.3){min(n20_d0.8_rcy0.75[i,c("P.X.noC","P.X.C")])}else{n20_d0.8_rcy0.75[i,"P.X.noC"]}
})) <= 0.05)/dim(n20_d0.8_rcy0.75)[1]

n20_d0.5_rcy0 <- read.csv("n20_d0.5_rcy0.csv")
sum(unlist(lapply(1:dim(n20_d0.5_rcy0)[1], function(i){
  if(n20_d0.5_rcy0[i,"Rcy"] >= 0.3){min(n20_d0.5_rcy0[i,c("P.X.noC","P.X.C")])}else{n20_d0.5_rcy0[i,"P.X.noC"]}
})) <= 0.05)/dim(n20_d0.5_rcy0)[1]

n20_d0.5_rcy0.25 <- read.csv("n20_d0.5_rcy0.25.csv")
sum(unlist(lapply(1:dim(n20_d0.5_rcy0.25)[1], function(i){
  if(n20_d0.5_rcy0.25[i,"Rcy"] >= 0.3){min(n20_d0.5_rcy0.25[i,c("P.X.noC","P.X.C")])}else{n20_d0.5_rcy0.25[i,"P.X.noC"]}
})) <= 0.05)/dim(n20_d0.5_rcy0.25)[1]

n20_d0.5_rcy0.5 <- read.csv("n20_d0.5_rcy0.5.csv")
sum(unlist(lapply(1:dim(n20_d0.5_rcy0.5)[1], function(i){
  if(n20_d0.5_rcy0.5[i,"Rcy"] >= 0.3){min(n20_d0.5_rcy0.5[i,c("P.X.noC","P.X.C")])}else{n20_d0.5_rcy0.5[i,"P.X.noC"]}
})) <= 0.05)/dim(n20_d0.5_rcy0.5)[1]

n20_d0.5_rcy0.75 <- read.csv("n20_d0.5_rcy0.75.csv")
sum(unlist(lapply(1:dim(n20_d0.5_rcy0.75)[1], function(i){
  if(n20_d0.5_rcy0.75[i,"Rcy"] >= 0.3){min(n20_d0.5_rcy0.75[i,c("P.X.noC","P.X.C")])}else{n20_d0.5_rcy0.75[i,"P.X.noC"]}
})) <= 0.05)/dim(n20_d0.5_rcy0.75)[1]

n20_d0.3_rcy0 <- read.csv("n20_d0.3_rcy0.csv")
sum(unlist(lapply(1:dim(n20_d0.3_rcy0)[1], function(i){
  if(n20_d0.3_rcy0[i,"Rcy"] >= 0.3){min(n20_d0.3_rcy0[i,c("P.X.noC","P.X.C")])}else{n20_d0.3_rcy0[i,"P.X.noC"]}
})) <= 0.05)/dim(n20_d0.3_rcy0)[1]

n20_d0.3_rcy0.25 <- read.csv("n20_d0.3_rcy0.25.csv")
sum(unlist(lapply(1:dim(n20_d0.3_rcy0.25)[1], function(i){
  if(n20_d0.3_rcy0.25[i,"Rcy"] >= 0.3){min(n20_d0.3_rcy0.25[i,c("P.X.noC","P.X.C")])}else{n20_d0.3_rcy0.25[i,"P.X.noC"]}
})) <= 0.05)/dim(n20_d0.3_rcy0.25)[1]

n20_d0.3_rcy0.5 <- read.csv("n20_d0.3_rcy0.5.csv")
sum(unlist(lapply(1:dim(n20_d0.3_rcy0.5)[1], function(i){
  if(n20_d0.3_rcy0.5[i,"Rcy"] >= 0.3){min(n20_d0.3_rcy0.5[i,c("P.X.noC","P.X.C")])}else{n20_d0.3_rcy0.5[i,"P.X.noC"]}
})) <= 0.05)/dim(n20_d0.3_rcy0.5)[1]

n20_d0.3_rcy0.75 <- read.csv("n20_d0.3_rcy0.75.csv")
sum(unlist(lapply(1:dim(n20_d0.3_rcy0.75)[1], function(i){
  if(n20_d0.3_rcy0.75[i,"Rcy"] >= 0.3){min(n20_d0.3_rcy0.75[i,c("P.X.noC","P.X.C")])}else{n20_d0.3_rcy0.75[i,"P.X.noC"]}
})) <= 0.05)/dim(n20_d0.3_rcy0.75)[1]

###################### SIMULATION RESULTS ##########################

######## Input Format for PowerSimData ######## 
# Number of Simulations,
# Number of participants in Condition 0,
# Number of participants in Condition 1,
# Effect size of x on y in the population,
# Correlation of c and y in the population,
# Standard deviation of y in Condition 0 in the population,
# Standard deviation of y in Condition 1 in the population,
# Standard deviation of c in Condition 0 in the population,
# Standard deviation of c in Condition 1 in the population

# n = 20, rcy = 0 ####
n20_d0_rcy0 <- PowerSimData(10000,20,20,0,0,1,1,1,1)
n20_d0.2_rcy0 <- PowerSimData(10000,20,20,0.2,0,1,1,1,1)
n20_d0.3_rcy0 <- PowerSimData(10000,20,20,0.3,0,1,1,1,1)
n20_d0.4_rcy0 <- PowerSimData(10000,20,20,0.4,0,1,1,1,1)
n20_d0.5_rcy0 <- PowerSimData(10000,20,20,0.5,0,1,1,1,1)
n20_d0.8_rcy0 <- PowerSimData(10000,20,20,0.8,0,1,1,1,1)

# n = 20, rcy = 0.25 ####
n20_d0_rcy0.25 <- PowerSimData(10000,20,20,0,0.25,1,1,1,1)
n20_d0.2_rcy0.25 <- PowerSimData(10000,20,20,0.2,0.25,1,1,1,1)
n20_d0.3_rcy0.25 <- PowerSimData(10000,20,20,0.3,0.25,1,1,1,1)
n20_d0.4_rcy0.25 <- PowerSimData(10000,20,20,0.4,0.25,1,1,1,1)
n20_d0.5_rcy0.25 <- PowerSimData(10000,20,20,0.5,0.25,1,1,1,1)
n20_d0.8_rcy0.25 <- PowerSimData(10000,20,20,0.8,0.25,1,1,1,1)

# n = 20, rcy = 0.5 ####
n20_d0_rcy0.5 <- PowerSimData(10000,20,20,0,0.5,1,1,1,1)
n20_d0.2_rcy0.5 <- PowerSimData(10000,20,20,0.2,0.5,1,1,1,1)
n20_d0.3_rcy0.5 <- PowerSimData(10000,20,20,0.3,0.5,1,1,1,1)
n20_d0.4_rcy0.5 <- PowerSimData(10000,20,20,0.4,0.5,1,1,1,1)
n20_d0.5_rcy0.5 <- PowerSimData(10000,20,20,0.5,0.5,1,1,1,1)
n20_d0.8_rcy0.5 <- PowerSimData(10000,20,20,0.8,0.5,1,1,1,1)

# n = 20, rcy = 0.75 ####
n20_d0_rcy0.75 <- PowerSimData(10000,20,20,0,0.75,1,1,1,1)
n20_d0.2_rcy0.75 <- PowerSimData(10000,20,20,0.2,0.75,1,1,1,1)
n20_d0.3_rcy0.75 <- PowerSimData(10000,20,20,0.3,0.75,1,1,1,1)
n20_d0.4_rcy0.75 <- PowerSimData(10000,20,20,0.4,0.75,1,1,1,1)
n20_d0.5_rcy0.75 <- PowerSimData(10000,20,20,0.5,0.75,1,1,1,1)
n20_d0.8_rcy0.75 <- PowerSimData(10000,20,20,0.8,0.75,1,1,1,1)

# n = 50, rcy = 0 ####
n50_d0_rcy0 <- PowerSimData(10000,50,50,0,0,1,1,1,1)
n50_d0.2_rcy0 <- PowerSimData(10000,50,50,0.2,0,1,1,1,1)
n50_d0.3_rcy0 <- PowerSimData(10000,50,50,0.3,0,1,1,1,1)
n50_d0.4_rcy0 <- PowerSimData(10000,50,50,0.4,0,1,1,1,1)
n50_d0.5_rcy0 <- PowerSimData(10000,50,50,0.5,0,1,1,1,1)
n50_d0.8_rcy0 <- PowerSimData(10000,50,50,0.8,0,1,1,1,1)

# n = 50, rcy = 0.25 ####
n50_d0_rcy0.25 <- PowerSimData(10000,50,50,0,0.25,1,1,1,1)
n50_d0.2_rcy0.25 <- PowerSimData(10000,50,50,0.2,0.25,1,1,1,1)
n50_d0.3_rcy0.25 <- PowerSimData(10000,50,50,0.3,0.25,1,1,1,1)
n50_d0.4_rcy0.25 <- PowerSimData(10000,50,50,0.4,0.25,1,1,1,1)
n50_d0.5_rcy0.25 <- PowerSimData(10000,50,50,0.5,0.25,1,1,1,1)
n50_d0.8_rcy0.25 <- PowerSimData(10000,50,50,0.8,0.25,1,1,1,1)

# n = 50, rcy = 0.5 ####
n50_d0_rcy0.5 <- PowerSimData(10000,50,50,0,0.5,1,1,1,1)
n50_d0.2_rcy0.5 <- PowerSimData(10000,50,50,0.2,0.5,1,1,1,1)
n50_d0.3_rcy0.5 <- PowerSimData(10000,50,50,0.3,0.5,1,1,1,1)
n50_d0.4_rcy0.5 <- PowerSimData(10000,50,50,0.4,0.5,1,1,1,1)
n50_d0.5_rcy0.5 <- PowerSimData(10000,50,50,0.5,0.5,1,1,1,1)
n50_d0.8_rcy0.5 <- PowerSimData(10000,50,50,0.8,0.5,1,1,1,1)

# n = 50, rcy = 0.75 ####
n50_d0_rcy0.75 <- PowerSimData(10000,50,50,0,0.75,1,1,1,1)
n50_d0.2_rcy0.75 <- PowerSimData(10000,50,50,0.2,0.75,1,1,1,1)
n50_d0.3_rcy0.75 <- PowerSimData(10000,50,50,0.3,0.75,1,1,1,1)
n50_d0.4_rcy0.75 <- PowerSimData(10000,50,50,0.4,0.75,1,1,1,1)
n50_d0.5_rcy0.75 <- PowerSimData(10000,50,50,0.5,0.75,1,1,1,1)
n50_d0.8_rcy0.75 <- PowerSimData(10000,50,50,0.8,0.75,1,1,1,1)

# n = 100, rcy = 0 ####
n100_d0_rcy0 <- PowerSimData(10000,100,100,0,0,1,1,1,1)
n100_d0.2_rcy0 <- PowerSimData(10000,100,100,0.2,0,1,1,1,1)
n100_d0.3_rcy0 <- PowerSimData(10000,100,100,0.3,0,1,1,1,1)
n100_d0.4_rcy0 <- PowerSimData(10000,100,100,0.4,0,1,1,1,1)
n100_d0.5_rcy0 <- PowerSimData(10000,100,100,0.5,0,1,1,1,1)
n100_d0.8_rcy0 <- PowerSimData(10000,100,100,0.8,0,1,1,1,1)

# n = 100, rcy = 0.25 ####
n100_d0_rcy0.25 <- PowerSimData(10000,100,100,0,0.25,1,1,1,1)
n100_d0.2_rcy0.25 <- PowerSimData(10000,100,100,0.2,0.25,1,1,1,1)
n100_d0.3_rcy0.25 <- PowerSimData(10000,100,100,0.3,0.25,1,1,1,1)
n100_d0.4_rcy0.25 <- PowerSimData(10000,100,100,0.4,0.25,1,1,1,1)
n100_d0.5_rcy0.25 <- PowerSimData(10000,100,100,0.5,0.25,1,1,1,1)
n100_d0.8_rcy0.25 <- PowerSimData(10000,100,100,0.8,0.25,1,1,1,1)

# n = 100, rcy = 0.5 ####
n100_d0_rcy0.5 <- PowerSimData(10000,100,100,0,0.5,1,1,1,1)
n100_d0.2_rcy0.5 <- PowerSimData(10000,100,100,0.2,0.5,1,1,1,1)
n100_d0.3_rcy0.5 <- PowerSimData(10000,100,100,0.3,0.5,1,1,1,1)
n100_d0.4_rcy0.5 <- PowerSimData(10000,100,100,0.4,0.5,1,1,1,1)
n100_d0.5_rcy0.5 <- PowerSimData(10000,100,100,0.5,0.5,1,1,1,1)
n100_d0.8_rcy0.5 <- PowerSimData(10000,100,100,0.8,0.5,1,1,1,1)

# n = 100, rcy = 0.75 ####
n100_d0_rcy0.75 <- PowerSimData(10000,100,100,0,0.75,1,1,1,1)
n100_d0.2_rcy0.75 <- PowerSimData(10000,100,100,0.2,0.75,1,1,1,1)
n100_d0.3_rcy0.75 <- PowerSimData(10000,100,100,0.3,0.75,1,1,1,1)
n100_d0.4_rcy0.75 <- PowerSimData(10000,100,100,0.4,0.75,1,1,1,1)
n100_d0.5_rcy0.75 <- PowerSimData(10000,100,100,0.5,0.75,1,1,1,1)
n100_d0.8_rcy0.75 <- PowerSimData(10000,100,100,0.8,0.75,1,1,1,1)

# Save output for future use ####
write.table(n20_d0.2_rcy0[[2]], "n20_d0.2_rcy0.csv",append = TRUE,sep = ",",row.names = FALSE,col.names = TRUE)
write.table(n20_d0.3_rcy0[[2]], "n20_d0.3_rcy0.csv",append = TRUE,sep = ",",row.names = FALSE,col.names = TRUE)
write.table(n20_d0.4_rcy0[[2]], "n20_d0.4_rcy0.csv",append = TRUE,sep = ",",row.names = FALSE,col.names = TRUE)
write.table(n20_d0.5_rcy0[[2]], "n20_d0.5_rcy0.csv",append = TRUE,sep = ",",row.names = FALSE,col.names = TRUE)
write.table(n20_d0.8_rcy0[[2]], "n20_d0.8_rcy0.csv",append = TRUE,sep = ",",row.names = FALSE,col.names = TRUE)
write.table(n20_d0_rcy0.25[[2]], "n20_d0_rcy0.25.csv",append = TRUE,sep = ",",row.names = FALSE,col.names = TRUE)
write.table(n20_d0.2_rcy0.25[[2]], "n20_d0.2_rcy0.25.csv",append = TRUE,sep = ",",row.names = FALSE,col.names = TRUE)
write.table(n20_d0.3_rcy0.25[[2]], "n20_d0.3_rcy0.25.csv",append = TRUE,sep = ",",row.names = FALSE,col.names = TRUE)
write.table(n20_d0.4_rcy0.25[[2]], "n20_d0.4_rcy0.25.csv",append = TRUE,sep = ",",row.names = FALSE,col.names = TRUE)
write.table(n20_d0.5_rcy0.25[[2]], "n20_d0.5_rcy0.25.csv",append = TRUE,sep = ",",row.names = FALSE,col.names = TRUE)
write.table(n20_d0.8_rcy0.25[[2]], "n20_d0.8_rcy0.25.csv",append = TRUE,sep = ",",row.names = FALSE,col.names = TRUE)
write.table(n20_d0_rcy0.5[[2]], "n20_d0_rcy0.5.csv",append = TRUE,sep = ",",row.names = FALSE,col.names = TRUE)
write.table(n20_d0.2_rcy0.5[[2]], "n20_d0.2_rcy0.5.csv",append = TRUE,sep = ",",row.names = FALSE,col.names = TRUE)
write.table(n20_d0.3_rcy0.5[[2]], "n20_d0.3_rcy0.5.csv",append = TRUE,sep = ",",row.names = FALSE,col.names = TRUE)
write.table(n20_d0.4_rcy0.5[[2]], "n20_d0.4_rcy0.5.csv",append = TRUE,sep = ",",row.names = FALSE,col.names = TRUE)
write.table(n20_d0.5_rcy0.5[[2]], "n20_d0.5_rcy0.5.csv",append = TRUE,sep = ",",row.names = FALSE,col.names = TRUE)
write.table(n20_d0.8_rcy0.5[[2]], "n20_d0.8_rcy0.5.csv",append = TRUE,sep = ",",row.names = FALSE,col.names = TRUE)
write.table(n20_d0_rcy0.75[[2]], "n20_d0_rcy0.75.csv",append = TRUE,sep = ",",row.names = FALSE,col.names = TRUE)
write.table(n20_d0.2_rcy0.75[[2]], "n20_d0.2_rcy0.75.csv",append = TRUE,sep = ",",row.names = FALSE,col.names = TRUE)
write.table(n20_d0.3_rcy0.75[[2]], "n20_d0.3_rcy0.75.csv",append = TRUE,sep = ",",row.names = FALSE,col.names = TRUE)
write.table(n20_d0.4_rcy0.75[[2]], "n20_d0.4_rcy0.75.csv",append = TRUE,sep = ",",row.names = FALSE,col.names = TRUE)
write.table(n20_d0.5_rcy0.75[[2]], "n20_d0.5_rcy0.75.csv",append = TRUE,sep = ",",row.names = FALSE,col.names = TRUE)
write.table(n20_d0.8_rcy0.75[[2]], "n20_d0.8_rcy0.75.csv",append = TRUE,sep = ",",row.names = FALSE,col.names = TRUE)
write.table(n50_d0_rcy0[[2]], "n50_d0_rcy0.csv",append = TRUE,sep = ",",row.names = FALSE,col.names = TRUE)
write.table(n50_d0.2_rcy0[[2]], "n50_d0.2_rcy0.csv",append = TRUE,sep = ",",row.names = FALSE,col.names = TRUE)
write.table(n50_d0.3_rcy0[[2]], "n50_d0.3_rcy0.csv",append = TRUE,sep = ",",row.names = FALSE,col.names = TRUE)
write.table(n50_d0.4_rcy0[[2]], "n50_d0.4_rcy0.csv",append = TRUE,sep = ",",row.names = FALSE,col.names = TRUE)
write.table(n50_d0.5_rcy0[[2]], "n50_d0.5_rcy0.csv",append = TRUE,sep = ",",row.names = FALSE,col.names = TRUE)
write.table(n50_d0.8_rcy0[[2]], "n50_d0.8_rcy0.csv",append = TRUE,sep = ",",row.names = FALSE,col.names = TRUE)
write.table(n50_d0_rcy0.25[[2]], "n50_d0_rcy0.25.csv",append = TRUE,sep = ",",row.names = FALSE,col.names = TRUE)
write.table(n50_d0.2_rcy0.25[[2]], "n50_d0.2_rcy0.25.csv",append = TRUE,sep = ",",row.names = FALSE,col.names = TRUE)
write.table(n50_d0.3_rcy0.25[[2]], "n50_d0.3_rcy0.25.csv",append = TRUE,sep = ",",row.names = FALSE,col.names = TRUE)
write.table(n50_d0.4_rcy0.25[[2]], "n50_d0.4_rcy0.25.csv",append = TRUE,sep = ",",row.names = FALSE,col.names = TRUE)
write.table(n50_d0.5_rcy0.25[[2]], "n50_d0.5_rcy0.25.csv",append = TRUE,sep = ",",row.names = FALSE,col.names = TRUE)
write.table(n50_d0.8_rcy0.25[[2]], "n50_d0.8_rcy0.25.csv",append = TRUE,sep = ",",row.names = FALSE,col.names = TRUE)
write.table(n50_d0_rcy0.5[[2]], "n50_d0_rcy0.5.csv",append = TRUE,sep = ",",row.names = FALSE,col.names = TRUE)
write.table(n50_d0.2_rcy0.5[[2]], "n50_d0.2_rcy0.5.csv",append = TRUE,sep = ",",row.names = FALSE,col.names = TRUE)
write.table(n50_d0.3_rcy0.5[[2]], "n50_d0.3_rcy0.5.csv",append = TRUE,sep = ",",row.names = FALSE,col.names = TRUE)
write.table(n50_d0.4_rcy0.5[[2]], "n50_d0.4_rcy0.5.csv",append = TRUE,sep = ",",row.names = FALSE,col.names = TRUE)
write.table(n50_d0.5_rcy0.5[[2]], "n50_d0.5_rcy0.5.csv",append = TRUE,sep = ",",row.names = FALSE,col.names = TRUE)
write.table(n50_d0.8_rcy0.5[[2]], "n50_d0.8_rcy0.5.csv",append = TRUE,sep = ",",row.names = FALSE,col.names = TRUE)
write.table(n50_d0_rcy0.75[[2]], "n50_d0_rcy0.75.csv",append = TRUE,sep = ",",row.names = FALSE,col.names = TRUE)
write.table(n50_d0.2_rcy0.75[[2]], "n50_d0.2_rcy0.75.csv",append = TRUE,sep = ",",row.names = FALSE,col.names = TRUE)
write.table(n50_d0.3_rcy0.75[[2]], "n50_d0.3_rcy0.75.csv",append = TRUE,sep = ",",row.names = FALSE,col.names = TRUE)
write.table(n50_d0.4_rcy0.75[[2]], "n50_d0.4_rcy0.75.csv",append = TRUE,sep = ",",row.names = FALSE,col.names = TRUE)
write.table(n50_d0.5_rcy0.75[[2]], "n50_d0.5_rcy0.75.csv",append = TRUE,sep = ",",row.names = FALSE,col.names = TRUE)
write.table(n50_d0.8_rcy0.75[[2]], "n50_d0.8_rcy0.75.csv",append = TRUE,sep = ",",row.names = FALSE,col.names = TRUE)
write.table(n100_d0_rcy0[[2]], "n100_d0_rcy0.csv",append = TRUE,sep = ",",row.names = FALSE,col.names = TRUE)
write.table(n100_d0.2_rcy0[[2]], "n100_d0.2_rcy0.csv",append = TRUE,sep = ",",row.names = FALSE,col.names = TRUE)
write.table(n100_d0.3_rcy0[[2]], "n100_d0.3_rcy0.csv",append = TRUE,sep = ",",row.names = FALSE,col.names = TRUE)
write.table(n100_d0.4_rcy0[[2]], "n100_d0.4_rcy0.csv",append = TRUE,sep = ",",row.names = FALSE,col.names = TRUE)
write.table(n100_d0.5_rcy0[[2]], "n100_d0.5_rcy0.csv",append = TRUE,sep = ",",row.names = FALSE,col.names = TRUE)
write.table(n100_d0.8_rcy0[[2]], "n100_d0.8_rcy0.csv",append = TRUE,sep = ",",row.names = FALSE,col.names = TRUE)
write.table(n100_d0_rcy0.25[[2]], "n100_d0_rcy0.25.csv",append = TRUE,sep = ",",row.names = FALSE,col.names = TRUE)
write.table(n100_d0.2_rcy0.25[[2]], "n100_d0.2_rcy0.25.csv",append = TRUE,sep = ",",row.names = FALSE,col.names = TRUE)
write.table(n100_d0.3_rcy0.25[[2]], "n100_d0.3_rcy0.25.csv",append = TRUE,sep = ",",row.names = FALSE,col.names = TRUE)
write.table(n100_d0.4_rcy0.25[[2]], "n100_d0.4_rcy0.25.csv",append = TRUE,sep = ",",row.names = FALSE,col.names = TRUE)
write.table(n100_d0.5_rcy0.25[[2]], "n100_d0.5_rcy0.25.csv",append = TRUE,sep = ",",row.names = FALSE,col.names = TRUE)
write.table(n100_d0.8_rcy0.25[[2]], "n100_d0.8_rcy0.25.csv",append = TRUE,sep = ",",row.names = FALSE,col.names = TRUE)
write.table(n100_d0_rcy0.5[[2]], "n100_d0_rcy0.5.csv",append = TRUE,sep = ",",row.names = FALSE,col.names = TRUE)
write.table(n100_d0.2_rcy0.5[[2]], "n100_d0.2_rcy0.5.csv",append = TRUE,sep = ",",row.names = FALSE,col.names = TRUE)
write.table(n100_d0.3_rcy0.5[[2]], "n100_d0.3_rcy0.5.csv",append = TRUE,sep = ",",row.names = FALSE,col.names = TRUE)
write.table(n100_d0.4_rcy0.5[[2]], "n100_d0.4_rcy0.5.csv",append = TRUE,sep = ",",row.names = FALSE,col.names = TRUE)
write.table(n100_d0.5_rcy0.5[[2]], "n100_d0.5_rcy0.5.csv",append = TRUE,sep = ",",row.names = FALSE,col.names = TRUE)
write.table(n100_d0.8_rcy0.5[[2]], "n100_d0.8_rcy0.5.csv",append = TRUE,sep = ",",row.names = FALSE,col.names = TRUE)
write.table(n100_d0_rcy0.75[[2]], "n100_d0_rcy0.75.csv",append = TRUE,sep = ",",row.names = FALSE,col.names = TRUE)
write.table(n100_d0.2_rcy0.75[[2]], "n100_d0.2_rcy0.75.csv",append = TRUE,sep = ",",row.names = FALSE,col.names = TRUE)
write.table(n100_d0.3_rcy0.75[[2]], "n100_d0.3_rcy0.75.csv",append = TRUE,sep = ",",row.names = FALSE,col.names = TRUE)
write.table(n100_d0.4_rcy0.75[[2]], "n100_d0.4_rcy0.75.csv",append = TRUE,sep = ",",row.names = FALSE,col.names = TRUE)
write.table(n100_d0.5_rcy0.75[[2]], "n100_d0.5_rcy0.75.csv",append = TRUE,sep = ",",row.names = FALSE,col.names = TRUE)
write.table(n100_d0.8_rcy0.75[[2]], "n100_d0.8_rcy0.75.csv",append = TRUE,sep = ",",row.names = FALSE,col.names = TRUE)

###################### END ######################

#### Sandbox ####
sb1 <- PowerSimData(100,20,20,0,0,1,1,1,1)
sb2 <- PowerSimData(100,20,20,0.5,0,1,1,1,1)
head(sb1[[1]])
head(sb2[[1]])
head(sb1[[2]])
head(sb2[[2]])

x <- rnorm(100)
y <- rnorm(100)
c <- rnorm(100)
dat <- as.data.frame(c(x,y,c))
coef(summary(lm(y ~ x + c + x*c, data = dat)))[,4] 
summary(lm(y ~ x + c + x*c, data = dat))