#basic occupancy model:
#ecological model: fixed effects for al, pb and ddt and random intercept for OTU and Lake
#observation model: fixed effects for year, and random intercepts for OTU and Lake

# setwd("C:/Users/diana.bowler/OneDrive - NINA/eDNA")

# Specify model in BUGS language
sink("eDNA_basic_ST.txt")
cat("
    model {
    
    # Priors
    
    #observational process
    mean.p ~ dunif(0, 1)         # Detection intercept on prob. scale
    alpha0 <- logit(mean.p)      # Detection intercept
    alpha1 ~ dunif(-20, 20)      # Detection slope on year 

    #ecological processes
    mean.psi ~ dunif(0, 1)       # Occupancy intercept on prob. scale
    beta0 <- logit(mean.psi)     # Occupancy intercept
    beta1 ~ dunif(-20, 20)       # Occupancy slope on Al
    beta2 ~ dunif(-20, 20)       # Occupancy slope on Pb
    beta3 ~ dunif(-20, 20)       # Occupancy slope on ddt

    #random effects:

    #random intercept on occupancy
    random.otu.beta.sd ~ dunif(0,10)
    random.otu.beta.tau <- pow(random.otu.beta.sd,-2)
    for(i in 1:nOTU){
      random.beta.otu[i]~dnorm(0,random.otu.beta.tau)
    }

    #random intercept on deection
    random.otu.alpha.sd ~ dunif(0,10)
    random.otu.alpha.tau <- pow(random.otu.alpha.sd,-2)
    for(i in 1:nOTU){
      random.alpha.otu[i]~dnorm(0,random.otu.alpha.tau)
    }

    #random slopes on the covariates
    #beta1
    for(i in 1:nOTU){
      random.beta1.otu[i]~dnorm(0,beta1.tau)
    }
    beta1.sd ~ dunif(0,10)
    beta1.tau <- beta1.sd,-2)

    #beta2
    for(i in 1:nOTU){
      random.beta2.otu[i]~dnorm(0,beta2.tau)
    }
    beta2.sd ~ dunif(0,10)
    beta2.tau <- beta2.sd,-2)

    #beta1
    for(i in 1:nOTU){
      random.beta3.otu[i]~dnorm(0,beta3.tau)
    }
    beta3.sd ~ dunif(0,10)
    beta3.tau <- beta3.sd,-2)


    # Likelihood
    for (i in 1:M) {
        # Ecological process
        z[i] ~ dbern(psi[i])      # True occupancy z at site i
        logit(psi[i]) <- beta0 + beta1 * Al[i] + beta2 * Pb[i] + beta3 * ddt[i] + #fixed effects
                          random.beta.otu[OTU[i]]+ #random OTU intercept
                          Al[i] * random.beta1.otu[OTU[i]]+ #random OTU effect of AL
                          Pb[i] * random.beta2.otu[OTU[i]]+ #random OTU effect of Pb
                          ddt[i] * random.beta3.otu[OTU[i]]+ #random OTU effect of ddt

    for (j in 1:J) {
    # Observation model for the actual observations -  observation process
      y[i,j] ~ dbern(z[i] * p[i,j])    # Detection-nondetection at i and j
      logit(p[i,j]) <- alpha0 + alpha1 * Year[i] +  random.alpha.otu[OTU[i]]
    }

    }
    
    ## Derived quantities
    N.occ <- sum(z[])       # Number of occupied sites among sample of M
    psi.fs <- N.occ/M
    
    }
    ",fill = TRUE)
sink()

# load the data
load("bugs.data.ST.Rdata")

# Initial values: must give for same quantities as priors given !
zst <- as.numeric(apply(bugs.data.ST$y, 1, max))        
inits <- function(){ list (z=zst) }

# Parameters monitored
params <- c("mean.p","mean.psi","alpha0", "alpha1", "beta0", "beta1", "beta2","beta3","beta1.sd","beta2.sd","beta3.sd") 

# MCMC settings
# ni <- 250000   ;   nt <- 100   ;   nb <- 50000   ;   nc <- 3
ni <- 250   ;   nt <- 10   ;   nb <- 50   ;   nc <- 3

# Call JAGS and summarize posteriors
library(jagsUI)
out1_ST <- jags(bugs.data.ST, inits=inits, params, "eDNA_basic_ST.txt", n.chains = nc,
              n.thin = nt, n.iter = ni, n.burnin = nb)
print(out1_ST, dig = 3)

#to retrieve a different parameter at a later stage
out1 <- update(out1_ST,parameters.to.monitor="random.beta.otu",n.iter=1000)
