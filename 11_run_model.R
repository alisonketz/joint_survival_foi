

############################################################################################
############################################################################################
############################################################################################
###
### Model Statement
###
############################################################################################
############################################################################################
############################################################################################

modelcode <- nimbleCode({


  ##############################
  ### Priors
  ##############################
  
  beta_sex ~ dnorm(0, .01)

  ##############################
  ### Force of infection model
  ##############################
  mprec  ~ dgamma(1, 1)
  fprec  ~ dgamma(1, 1)
  mprec1 <- 0.0000001 * mprec
  fprec1 <- 0.0000001 * fprec
  m_age_foi[1] ~ dnorm(0, mprec1)
  f_age_foi[1] ~ dnorm(0, fprec1)
  m_age_foi[2] ~ dnorm(0, mprec1)
  f_age_foi[2] ~ dnorm(0, fprec1)
  for (i in 3:n_agem) {
    m_age_foi[i]~dnorm(2 * m_age_foi[i-1] - m_age_foi[i-2], mprec)
  }
  for (i in 3:n_agef) {
    f_age_foi[i]~dnorm(2 * f_age_foi[i-1] - f_age_foi[i-2], fprec)
  }
  m_age_foi_mu <- mean(m_age_foi[1:n_agem])
  f_age_foi_mu <- mean(f_age_foi[1:n_agef])

  # Period effects
  tmprec  ~ dgamma(1, 1)
  tfprec  ~ dgamma(1, 1)

  f_period_foi[1:n_period] ~ dcar_normal(adj = adj_period[1:n_adj_period],
                                  weights = weights_period[1:n_adj_period],
                                  num = num_period[1:n_period],
                                  tau = tfprec,
                                  zero_mean = 1)
  
  m_period_foi[1:n_period] ~ dcar_normal(adj = adj_period[1:n_adj_period],
                                  weights = weights_period[1:n_adj_period],
                                  num = num_period[1:n_period],
                                  tau = tmprec,
                                  zero_mean = 1)

  #Beseg-York-Mollie-2 Spatial model
  # phi[1:n_sect] ~ dcar_normal(sadj[1:n_adj_sp],
  #                             sweights[1:n_adj_sp],
  #                             snum[1:n_sect],
  #                             1,
  #                             zero_mean = 1) #sum to zero constraint
  # hprec ~ dgamma(1, 1)
  # # hsd ~ T(dnorm(1, 1), 0, 10)
  # hprec1 <- hprec * 0.0000001
  # sprec ~ dgamma(1, 1)  
  # rho ~ dbeta(1, 1)# mixing parameter
  # for (j in 1:n_sect){
  #   hetero[j]~dnorm(0, hprec1)
  #   space[j] <- (1 / sqrt(sprec)) * (sqrt((1 - rho)) * hetero[j] +
  #                sqrt(rho / scale) * phi[j])
  # }

  ############################################################
  ############################################################
  ### Age/period Survival Model
  ############################################################
  ############################################################

  ##############################
  ### Susceptibles
  ##############################

  #Priors for intercept and covariate
  beta0_sus_temp ~ dnorm(0, .01)
  sus_mix ~ dunif(-1, 1)
  beta0_sus <- beta0_sus_temp * sus_mix

  #Priors for Age and Period effects
  #Age effects
  for (k in 1:nknots_age) {
    b_age[k] ~ dnorm(0, tau_age)
  }
  tau_age ~ dgamma(1, 1)

  for (t in 1:nT_age_surv) {
    age_effect_survival[t] <- inprod(b_age[1:nknots_age],
                            Z_age[t, 1:nknots_age])
  }

  #Period effects
  for (k in 1:nknots_period) {
    b_period[k] ~ dnorm(0, tau_period)
  }
  tau_period ~ dgamma(1, 1)
  for (t in 1:nT_period_surv) {
    period_effect_surv[t] <- inprod(b_period[1:nknots_period],
                               Z_period[t, 1:nknots_period])
  }

  ##################################
  ## Infected survival intercept
  ##################################

  #Priors for intercept and covariate
  beta0_inf_temp ~ dnorm(0, .01)
  inf_mix ~ dunif(-1, 1)
  beta0_inf <- beta0_inf_temp * inf_mix

  ############################################################
  ## setting up period effects across full study timeline
  ############################################################

  period_effect_survival[(nT_period_presurv + 1):nT_overall_period] <- period_effect_surv[1:nT_period_surv]

  #######################################################################
  #######################################################################
  ## Likelihoods of Joint Model
  #######################################################################
  #######################################################################

  #######################################################################
  ###
  ###   User defined distribution for likelihood for
  ###   infected harvest deer
  ###
  ###   d_fit_hunt_pos
  ###   Overleaf Equation 3
  ###
  #######################################################################

  for (i in 1:nInfHarvest) {

    y_hunt_pos[i] ~ dInfHarvest(
                  a = hunt_pos_ageweeks[i], #age (weeks) at harvest
                  sex = hunt_pos_sex[i],
                  age2date = hunt_pos_age2date[i],
                  beta_sex = beta_sex,
                  beta0_sus = beta0_sus,
                  beta0_inf = beta0_inf,
                  age_effect_surv = age_effect_survival[1:nT_age_surv],
                  period_effect_surv = period_effect_survival[1:nT_overall_period],
                  f_age_foi = f_age_foi[1:n_agef],
                  m_age_foi = m_age_foi[1:n_agem],
                  age_lookup_f = age_lookup_f[1:n_age_lookup_f],
                  age_lookup_m = age_lookup_m[1:n_age_lookup_m],
                  period_lookup = period_lookup[1:nT_overall_period],
                  f_period_foi = f_period_foi[1:n_period],
                  m_period_foi = m_period_foi[1:n_period],
                  space = space[sect_hunt_pos[i]]
                  )
  }

})#end model statement

#######################################
### Data for Model Fitting
#######################################

nimData <- list(Z_period = Z_period,
                Z_age = Z_age,
                num_period = num_period,
                adj_period = adj_period,
                weights_period = weights_period,
                age_lookup_f = age_lookup_f,
                age_lookup_m = age_lookup_m,
                period_effect_survival = period_effect_survival,
                y_hunt_pos = rep(1, nrow(d_fit_hunt_pos)),
                hunt_pos_ageweeks = d_fit_hunt_pos$ageweeks,
                hunt_pos_sex = d_fit_hunt_pos$sex,
                hunt_pos_age2date = d_fit_hunt_pos$birthweek - 1
                )


#######################################
### Constants for MCMC
#######################################

nimConsts <- list(
                  nT_overall_period = nT_overall,
                  nT_period_presurv = nT_period_presurv,
                  nknots_age = nknots_age,
                  nknots_period = nknots_period,
                  nT_age_surv = nT_age_surv,
                  nT_period_surv = nT_period_surv,
                  n_age = n_age,
                  n_agef = n_agef,
                  n_agem = n_agem,
                  n_period = n_period,
                  n_adj_period = n_adj_period,
                  # n_adj_sp = n_adj_sp,
                  n_sect = n_sect,
                  # scale = scale,
                  space = rep(0, n_sect),
                  n_age_lookup_f = length(age_lookup_f),
                  n_age_lookup_m = length(age_lookup_m),
                  period_lookup = period_lookup,
                  nInfHarvest = nrow(d_fit_hunt_pos),
                  sect_hunt_pos = sect_hunt_pos
                  )


#######################################
### Initial Values for MCMC
#######################################
initsFun <- function()list(
                          beta_sex = rnorm(1, 0, .01),
                          beta0_sus_temp = rnorm(1, -6, 0.0001),
                          sus_mix = 1,
                          beta0_inf_temp = rnorm(1, -3, 0.0001),
                          inf_mix = 1,
                          b_age = rnorm(nknots_age) * 10^-4,
                          b_period = rnorm(nknots_period) * 10^-4,
                          tau_period = runif(1, .1, 1),
                          tau_age = runif(1, .1, .4),
                          # hprec = runif(1, .01, 1),
                          # sprec = runif(1, 3, 4),
                          mprec = runif(1, 1.5, 1.7),
                          fprec = runif(1, 2.7, 4.2),
                          # tfsd = runif(1,1,2),
                          # tmsd = runif(1,3,4),
                          tmprec = runif(1, 4.2, 6.8),
                          tfprec = runif(1, 2.27, 3.44),
                          # phi = seq(.0001, .0005, length = n_sect),
                          # hetero = seq(.0001, .0005, length = n_sect),
                          # rho = rbeta(1, 1, 1),
                          m_period_foi = seq(-2, 2, length = n_period),
                          f_period_foi = seq(-2, 2, length = n_period),
                          m_age_foi = seq(-6, -4, length = n_agem),
                          f_age_foi = seq(-7, -5, length = n_agef)
                          )
nimInits <- initsFun()

Rmodel <- nimbleModel(code = modelcode,
                      constants = nimConsts,
                      data = nimData,
                      inits = initsFun(),
                      calculate = FALSE,
                      check = FALSE
                      )
Rmodel$initializeInfo()

#######################################
### Parameters to trace in MCMC
#######################################

parameters <- c(
              "beta_sex",
              "beta0_sus",
              "beta0_inf",
              "age_effect_survival",
              "period_effect_surv",
              "tau_period",
              "tau_age",
              # "hprec",
              # "sprec",
             "f_age_foi",
              "m_age_foi",
              "mprec",
              "fprec",
              # "m_age_mu",
              # "f_age_mu",
              "f_period_foi",
              "m_period_foi",
              "tmprec",
              # "tmsd",
              "tfprec"
              # "tfsd",
              # "space",
              # "hetero",
              # "phi",
              # "rho"
               )

starttime <- Sys.time()
confMCMC <- configureMCMC(Rmodel,
                         monitors = parameters, 
                         thin = 1,
                         # enableWAIC = TRUE,
                         useConjugacy = FALSE)
nimMCMC <- buildMCMC(confMCMC)
Cnim <- compileNimble(Rmodel)
CnimMCMC <- compileNimble(nimMCMC,
                         project = Rmodel)
mcmcout <- runMCMC(CnimMCMC,
                  niter = 200,
                  nburnin = 100,
                  nchains = 1,
                  inits = initsFun,
                  samplesAsCodaMCMC = TRUE,
                  summary = TRUE
                  )

runtime <- difftime(Sys.time(),
                    starttime,
                    units = "min")


# reps  <- 2000
# bin <- reps * .5
# n_thin <- 1
# n_chains <- 3
# starttime <- Sys.time()
# cl <- makeCluster(n_chains)
# clusterExport(cl, c("modelcode",
#                     "initsFun",
#                     "nimData",
#                     "nimConsts",
#                     "parameters",
#                     "reps",
#                     "bin",
#                     "n_thin",
#                     "dSurvival",
#                     "dFOIcollar",
#                     "dFOIinf",
#                     "dFOIhunt",
#                     # "dSurvival_icap",
#                     "calc_age_inf_rec",
#                     "calc_age_inf_idead",
#                     "calc_age_inf_icap"
#                     ))
# for (j in seq_along(cl)) {
#   set.seed(j+1000)
#   init <- initsFun()
#   clusterExport(cl[j], "init")
# }
# mcmcout <-  mcmc.list(clusterEvalQ(cl, {
#   library(nimble)
#   library(coda)

#   Ccalc_age_inf_rec <- compileNimble(calc_age_inf_rec)
#   Ccalc_age_inf_idead <- compileNimble(calc_age_inf_idead)
#   Ccalc_age_inf_icap <- compileNimble(calc_age_inf_icap)
#   # Cpo_indx_fun <- compileNimble(po_indx_fun)
#   # Cpr_indx_fun <- compileNimble(pr_indx_fun)

#   ##############################################################
#   ###
#   ### Execute MCMC
#   ###
#   ##############################################################

#   Rmodel <- nimbleModel(code = modelcode,
#                         name = "modelcode",
#                         constants = nimConsts,
#                         data = nimData,
#                         inits = init)
#   confMCMC <- configureMCMC(Rmodel,
#                             monitors = parameters,
#                             thin = n_thin,
#                             useConjugacy = FALSE)
#   nimMCMC <- buildMCMC(confMCMC)
#   Cnim <- compileNimble(Rmodel)
#   CnimMCMC <- compileNimble(nimMCMC,
#                             project = Rmodel)
#   mcmcout <- runMCMC(CnimMCMC,
#                      niter = reps,
#                      nburnin = bin,
#                      samplesAsCodaMCMC = TRUE)

#   return(mcmcout)
# }))
# runtime <- difftime(Sys.time(), starttime, units = "min")
# stopCluster(cl)

save(mcmcout, file = "mcmcout.Rdata")
save(runtime, file = "runtime.Rdata")

#not calculating waic, because too many params would need to be traced
# posteriorSamplesMatrix <- rbind(mcmcout[[1]], mcmcout[[2]], mcmcout[[3]])
# CnimMCMC$run(5)   ## non-zero number of iterations
# nimble:::matrix2mv(posteriorSamplesMatrix, CnimMCMC$mvSamples)
# # CnimMCMC$enableWAIC <- TRUE
# waic_spline <- calculateWAIC(posteriorSamplesMatrix,Rmodel)


# waic_spline_covs <- mcmcout$WAIC
# save(waic_spline, file = "waic_spline.Rdata")



###
### save model run
###

# save(runtime,file="results/runtime.Rdata")
# save(mcmcout,file="results/mcmcout.Rdata")