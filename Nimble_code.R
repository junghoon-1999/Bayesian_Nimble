################################################################################
# Question 4
################################################################################

library(arm)

# Define the number of samples for the prior distribution for Se and Sp
p = 1000000

# Sample from the prior distributions
se <- rbeta(p, 21, 19)
sp <- rbeta(p, 82, 2)

# Define the distribution for prevalence and sample from it
prev <- rnorm(p, 0, 10)
prev <- invlogit(prev)

# The distribution of prevalence
hist(prev)

# Use the previously defined samles to obtain an estimate distribution for the PPV
ppv <- prev*se/(prev*se + (1-prev)*(1-sp))
hist(ppv)
print(c("The estimated mean value of PPV is:", mean(ppv)))

################################################################################
# Question 5
################################################################################

library(nimble)
library(igraph)
library(coda)
library(R6)

#######################################model1###################################
# Normalize the data in the X_hp and X_hi columns
ZX_hp <- (dataQ5$X_hp-mean(dataQ5$X_hp))/sd(dataQ5$X_hp)
ZX_hi <- (dataQ5$X_hi-mean(dataQ5$X_hi))/sd(dataQ5$X_hi)
y <- dataQ5$ncarriers

# Specify the statistical model
BinNormCode <- nimbleCode({
  
  # Specify the likelihood
  for(i in 1 : n) {
    y[i] ~ dbinom(p[i],ni[i])
    logitp[i] <- beta0 + beta1*X_hp[i] + beta2*X_hi[i]
    p[i] <- exp(logitp[i])/(1 + exp(logitp[i]))
    
  }
  
  # Specify the prior
  beta0 ~ dnorm(0, 5)
  beta1 ~ dnorm(0, 5)
  beta2 ~ dnorm(0, 5)
})

# Define the constant values in the model 
BinNormConsts <- list(n = 15)

# Define the data values
BinNormData <- list(y = y, X_hp = ZX_hp, X_hi = ZX_hi, ni = dataQ5$ni)

# Set initial values               
BinNormInits <- list(beta0=0.1, beta1=0.1, beta2=0.1) 

# Build the model 
BinNorm <- nimbleModel(code = BinNormCode, name = "BinNorm", constants = BinNormConsts,
                       data = BinNormData, inits<-BinNormInits)

# Compile the model 
CBinNorm<- compileNimble(BinNorm)

# Set up quantities to monitor 
BinNormConf <- configureMCMC(BinNorm, enableWAIC = TRUE, monitors = c('beta0','beta1', 'beta2'), print = TRUE) #M1 

# Build the MCMC algorithm
BinNormMCMC <- buildMCMC(BinNormConf)

# Compile the MCMC chain 
CBinNormMCMC <- compileNimble(BinNormMCMC, project = BinNorm)

set.seed(15)

# Set up initial values 
BinNormInits <- list(list(beta0=0, beta1=0, beta2=0), 
                     list(beta0=1, beta1=1, beta2=1))

# Run the MCMC to obtain points in the posterior
posterior <- runMCMC(CBinNormMCMC, niter = 200000, thin=5, nburnin=100000, 
                     summary = TRUE, samples = TRUE, nchains=2, 
                     samplesAsCodaMCMC=TRUE, inits = BinNormInits) 

# Combine the two chains
combinedchains <- mcmc.list(posterior$samples$chain1, posterior$samples$chain2)

# Plot the combined chains for all the parameter values
plot(combinedchains[,c('beta0','beta1','beta2')]) 

# obtain an ACF plot for the two chains
autocorr.plot(posterior$samples$chain1)
autocorr.plot(posterior$samples$chain2)

# Obtain the BRG value and the gelman plot for the combined chains
gelman.diag(combinedchains)
gelman.plot(combinedchains)

# Obtain the effective sample size and the WAIC for the model 
ESS1 <- effectiveSize(combinedchains)
WAIC1 <- calculateWAIC(CBinNormMCMC)

# Use the summary option to obtain the standard deviation of the posterior
post_sum <- posterior$summary$all.chains

# Use the summary function to obtain the Monte Carlo Error for the generated values
summary <- summary(combinedchains)

# Check if the Monte Carlo Error/posterior summary is smaller than 0.05
MCError <- summary$statistics[, 4]/ESS1[2]
print(MCError/post_sum[,3] < 0.05)
print('If all are True, the MCE/posteriorSD is samller than 0.05')

###############################Model 2##########################################

# Specify the statistical model
BinNormCode <- nimbleCode({
  
  # Specify the likelihood
  for(i in 1 : n) {
    y[i] ~ dbinom(p[i],ni[i])
    logitp[i] <- beta0  + beta1*X_hp[i] # Model 2
    p[i] <- exp(logitp[i])/(1 + exp(logitp[i]))
    
  }
  
  # Specify the prior
  beta0 ~ dnorm(0, 5)
  beta1 ~ dnorm(0, 5)
})

# Define the constant values in the model
BinNormConsts <- list(n = 15)

# Define the data values
BinNormData <- list(y = y, X_hp = ZX_hp, ni = dataQ5$ni)

# Set initial values                
BinNormInits <- list(beta0=0.1, beta1=0.1) 

# Build the model 
BinNorm2 <- nimbleModel(code = BinNormCode, name = "BinNorm2", constants = BinNormConsts,
                        data = BinNormData, inits<-BinNormInits)

# Compile the model 
CBinNorm2<- compileNimble(BinNorm2)

# Set up quantities to monitor
BinNormConf2 <- configureMCMC(BinNorm2, enableWAIC = TRUE, monitors = c('beta0','beta1'), print = TRUE) #M1 

# Build the MCMC algorithm
BinNormMCMC2 <- buildMCMC(BinNormConf2)

# Compile the MCMC chain 
CBinNormMCMC2 <- compileNimble(BinNormMCMC2, project = BinNorm2)

set.seed(15)

# Set up initial values 
BinNormInits2 <- list(list(beta0=0, beta1=0), 
                      list(beta0=1, beta1=1))

# Run the MCMC to obtain points in the posterior
posterior2 <- runMCMC(CBinNormMCMC2, niter = 200000, thin=5, nburnin=100000, 
                      summary = TRUE, samples = TRUE, nchains=2, 
                      samplesAsCodaMCMC=TRUE, inits = BinNormInits2) 

# Combine the two chains
combinedchains2 <- mcmc.list(posterior2$samples$chain1, posterior2$samples$chain2)

# Plot the combined chains for all the parameter values
plot(combinedchains2[,c('beta0','beta1')]) 

# obtain an ACF plot for the two chains
autocorr.plot(posterior2$samples$chain1)
autocorr.plot(posterior2$samples$chain2)

# Obtain the BRG value and the gelman plot for the combined chains
gelman.diag(combinedchains2)
gelman.plot(combinedchains2)

# Obtain the effective sample size and the WAIC for the model 
ESS2 <- effectiveSize(combinedchains2)
WAIC2 <- calculateWAIC(CBinNormMCMC2)

# Use the summary option to obtain the standard deviation of the posterior
post_sum2 <- posterior2$summary$all.chains

# Use the summary function to obtain the Monte Carlo Error for the generated values
summary2 <- summary(combinedchains2)

# Check if the Monte Carlo Error/posterior summary is smaller than 0.05
MCError2 <- summary2$statistics[, 4]/ESS2[2]
print(unlist(MCError2/post_sum2[,3] < 0.05))
print('If all are True, the MCE/posteriorSD is samller than 0.05')

####################################Model3######################################
# Specify the statistical model
BinNormCode <- nimbleCode({
  
  # Specify the likelihood
  for(i in 1 : n) {
    y[i] ~ dbinom(p[i],ni[i])
    logitp[i] <- beta0  + beta2*X_hi[i] # Model 2
    p[i] <- exp(logitp[i])/(1 + exp(logitp[i]))
    
  }
  
  # Specify the prior
  beta0 ~ dnorm(0, 5)
  beta2 ~ dnorm(0, 5)
})

# Define the constant values in the model 
BinNormConsts <- list(n = 15)

# Define the data values
BinNormData <- list(y = y, X_hi = ZX_hi, ni = dataQ5$ni)

# Set initial values                
BinNormInits <- list(beta0=0.1, beta2=0.1) 

# Build the model 
BinNorm3 <- nimbleModel(code = BinNormCode, name = "BinNorm3", constants = BinNormConsts,
                        data = BinNormData, inits<-BinNormInits)

# Compile the model 
CBinNorm3<- compileNimble(BinNorm3)

# Set up quantities to monitor
BinNormConf3 <- configureMCMC(BinNorm3, enableWAIC = TRUE, monitors = c('beta0','beta2'), print = TRUE) #M1 

# build the MCMC algorithm
BinNormMCMC3 <- buildMCMC(BinNormConf3)

# compile the MCMC chain 
CBinNormMCMC3 <- compileNimble(BinNormMCMC3, project = BinNorm3)

set.seed(15)

# Set up initial values 
BinNormInits3 <- list(list(beta0=0, beta2=0), 
                      list(beta0=1, beta2=1))

# Run the MCMC to obtain points in the posterior
posterior3 <- runMCMC(CBinNormMCMC3, niter = 200000, thin=5, nburnin=100000, 
                      summary = TRUE, samples = TRUE, nchains=2, 
                      samplesAsCodaMCMC=TRUE, inits = BinNormInits3) 

# Combine the two chains
combinedchains3 <- mcmc.list(posterior3$samples$chain1, posterior3$samples$chain2)

# Plot the combined chains for all the parameter values
plot(combinedchains3[,c('beta0','beta2')]) 

# obtain an ACF plot for the two chains
autocorr.plot(posterior3$samples$chain1)
autocorr.plot(posterior3$samples$chain2)

# Obtain the BRG value and the gelman plot for the combined chains
gelman.diag(combinedchains3)
gelman.plot(combinedchains3)

# Obtain the effective sample size and the WAIC for the model 
ESS3 <-effectiveSize(combinedchains3)
WAIC3 <- calculateWAIC(CBinNormMCMC3)

# Use the summary option to obtain the standard deviation of the posterior
post_sum3 <- posterior3$summary$all.chains

# Use the summary function to obtain the Monte Carlo Error for the generated values
summary3 <- summary(combinedchains3)

# Check if the Monte Carlo Error/posterior summary is smaller than 0.05
MCError3 <- summary3$statistics[, 4]/ESS3[2]
print(MCError3/post_sum3[,3] < 0.05)
print('If all are True, the MCE/posteriorSD is samller than 0.05')

###############################################################################
library(arm)

# Save the samples from the chain after burn-in
beta_1 <- posterior$samples$chain1
beta_2 <- posterior$samples$chain2

# Combine the two chains together
post_samp <- rbind(invlogit(beta_1[,1] + beta_1[,2]*mean(ZX_hp) + beta_1[,3]*mean(ZX_hi)), 
                   invlogit(beta_2[,1] + beta_2[,2]*mean(ZX_hp) + beta_2[,3]*mean(ZX_hi)))
# Plot the posterior samples
hist(post_samp)

# Obtain samples from Se and Sp
p = length(post_samp)
se <- rbeta(p, 21, 19)
sp <- rbeta(p, 82, 2)

# Use the previously obtained samples to get samples for the PPV
prev <- post_samp
ppv <- prev*se/(prev*se + (1-prev)*(1-sp))

# Plot the PPV and obtain the mean from the estimated distribution
hist(ppv)
mean(ppv)

################################################################################
# Question 6
################################################################################
# Part 1

# Generate 10000 values for each parameter
set.seed(15)
beta0_sample <- sample(c(beta_1[,1], beta_2[,1]), 10000, replace = TRUE)
beta1_sample <- sample(c(beta_1[,2], beta_2[,2]), 10000, replace = TRUE)
beta2_sample <- sample(c(beta_1[,3], beta_2[,3]), 10000, replace = TRUE)

# Define an empty data frame to store results 
q6_1data <- data.frame(matrix(ncol = 15, nrow = 10000))
colnames(q6_1data) <- c('Hospital 1', 'Hospital 2', 'Hospital 3','Hospital 4', 'Hospital 5',
                        'Hospital 6', 'Hospital 7', 'Hospital 8', 'Hospital 9', 'Hospital 10',
                        'Hospital 11', 'Hospital 12', 'Hospital 13', 'Hospital 14', 'Hospital 15')

# Fit the data to the previously generated parameter values
for ( i in 1:15){
  pred_samp <- dataQ5$ni[i]*invlogit(beta0_sample+beta1_sample*ZX_hp[i]+beta2_sample*ZX_hi[i])
  q6_1data[i] <- pred_samp
}

head(q6_1data)

# Part 2

# Create an empty data frame to store the results
q6_2data <- data.frame(matrix(ncol = 3, nrow = 10000))
colnames(q6_2data) <- c("mean","median","max")

# Store summary values for the each row in q6_1
for ( i in 1:10000){
  q6_2data[i,1] <- mean(as.numeric(q6_1data[i,]))
  q6_2data[i,2] <- median(as.numeric(q6_1data[i,]))
  q6_2data[i,3] <- max(as.numeric(q6_1data[i,]))
}

# Produce the distribution for each summary value
hist(q6_2data$mean)
hist(q6_2data$median)
hist(q6_2data$max)

# Part 3
# Calculate the actual summary values from the original data set 
true_mean = mean(dataQ5$ncarriers)
true_median = median(dataQ5$ncarriers)
true_max = max(dataQ5$ncarriers)

# Calculate the p-value for each summary value
mean_pvalue = min(sum(q6_2data[,1] < true_mean)/10000,(sum(q6_2data[,1] > true_mean))/10000)*2
median_pvalue =min(sum(q6_2data[,2] < true_median)/10000,(sum(q6_2data[,2] > true_median))/10000)*2
max_pvalue = min(sum(q6_2data[,3] < true_max)/10000,(sum(q6_2data[,3] > true_max))/10000)*2
