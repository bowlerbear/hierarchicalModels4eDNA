#basic occupancy model:
#ecological model: fixed effects for al, pb and ddt and random intercept for OTU and Lake
#observation model: fixed effects for year, and random intercepts for OTU and Lake

setwd("C:/Users/diana.bowler/OneDrive - NINA/eDNA")

# Specify model in BUGS language
sink("eDNA_basic.txt")
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

    #random OTU/Lake effects priors
    random.lake.beta.sd ~ dunif(0,10)
    random.lake.beta.tau <- pow(random.lake.beta.sd,-2)
    for(i in 1:nLake){
      random.beta.lake[i]~dnorm(0,random.lake.beta.tau)
    }

    random.otu.beta.sd ~ dunif(0,10)
    random.otu.beta.tau <- pow(random.otu.beta.sd,-2)
    for(i in 1:nOTU){
    random.beta.otu[i]~dnorm(0,random.otu.beta.tau)
    }

    random.lake.alpha.sd ~ dunif(0,10)
    random.lake.alpha.tau <- pow(random.lake.alpha.sd,-2)
    for(i in 1:nLake){
    random.alpha.lake[i]~dnorm(0,random.lake.alpha.tau)
    }
    
    random.otu.alpha.sd ~ dunif(0,10)
    random.otu.alpha.tau <- pow(random.otu.alpha.sd,-2)
    for(i in 1:nOTU){
    random.alpha.otu[i]~dnorm(0,random.otu.alpha.tau)
    }


    # Likelihood
    for (i in 1:M) {
        # Ecological process
        z[i] ~ dbern(psi[i])      # True occupancy z at site i
        logit(psi[i]) <- beta0 + beta1 * Al[i] + beta2 * Pb[i] + beta3 * ddt[i] + random.beta.otu[OTU[i]] + random.beta.lake[Lake[i]]

    for (j in 1:J) {
    # Observation model for the actual observations -  observation process
      y[i,j] ~ dbern(z[i] * p[i,j])    # Detection-nondetection at i and j
      logit(p[i,j]) <- alpha0 + alpha1 * Year[i] +  random.alpha.otu[OTU[i]] + random.alpha.lake[Lake[i]]
    }

    }
    
    ## Derived quantities
    N.occ <- sum(z[])       # Number of occupied sites among sample of M
    psi.fs <- N.occ/M
    
    }
    ",fill = TRUE)
sink()

# Initial values: must give for same quantities as priors given !
zst <- as.numeric(apply(bugs.data$y, 1, max))        
inits <- function(){ list (z=zst) }

# Parameters monitored
params <- c("alpha0", "alpha1", "beta0", "beta1", "beta2","beta3","alpha1") 

# MCMC settings
ni <- 2500   ;   nt <- 10   ;   nb <- 200   ;   nc <- 3

# Call JAGS and summarize posteriors
library(jagsUI)
out1 <- jags(bugs.data, inits=inits, params, "eDNA_basic.txt", n.chains = nc,
              n.thin = nt, n.iter = ni, n.burnin = nb)
print(out1, dig = 3)
