sink("BUGSmodel.txt")
cat("
    model {
    
    # Priors for fixed effects
    mean.p ~ dunif(0, 1)         # Detection intercept on prob. scale
    alpha0 <- logit(mean.p)      # Detection intercept
    alpha.depth ~ dnorm(0,0.001)  # Detection slope on depth
    mean.psi ~ dunif(0, 1)       # Occupancy intercept on prob. scale
    beta0 <- logit(mean.psi)     # Occupancy intercept
    beta.time ~ dnorm(0,0.001)  # Occupancy slope on time
    alpha.ts ~ dnorm(0,0.001)

    #Priors for random effects    

    #random lake effect
    lake.sd ~ dunif(0,10)
    lake.tau <- pow(lake.sd,-2)
    for (i in 1:n.Lakes){
      random.lake[i] ~ dnorm(0,lake.tau)
    }

    #OTU random effects on occupancy  
    state.sd ~ dunif(0,10)
    state.tau <- pow(state.sd,-2)
    for (i in 1:n.OTU){
      OTU.state[i] ~ dnorm(0,state.tau)
    }

    #OTU random effects on detection 
    det.sd ~ dunif(0,10)
    det.tau <- pow(det.sd,-2)
    for (i in 1:n.OTU){
      OTU.det[i] ~ dnorm(0,det.tau)
    }
  
    # Likelihood
    for (i in 1:M) {
    # True state model for the partially observed true state
    z[i] ~ dbern(psi[i])      # True occupancy z at site i
    logit(psi[i]) <- beta0 + beta.time * time[i] + random.lake[site[i]] + OTU.state[otu[i]]
    
    for (j in 1:J) {
    # Observation model for the actual observations
    y[i,j] ~ dbern(p.eff[i,j])    # Detection-nondetection at i and j
    p.eff[i,j] <- z[i] * p[i,j]   
    logit(p[i,j]) <- alpha0 + alpha.depth * depth[i] + alpha.ts * trophicStatus[i] + OTU.det[otu[i]]
    }
    }
    
    }
    ",fill = TRUE)
sink()