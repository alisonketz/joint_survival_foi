

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
  ### Force of infection model
  ##############################
    mprec  ~ dgamma(1, 1)
    fprec  ~ dgamma(1, 1)
    mprec1 <- 0.0000001 * mprec
    fprec1 <- 0.0000001 * fprec
    m_age[1] ~ dnorm(0, mprec1)
    f_age[1] ~ dnorm(0, fprec1)
    m_age[2] ~ dnorm(0, mprec1)
    f_age[2] ~ dnorm(0, fprec1)
    for (i in 3:n_agem) {
      m_age[i]~dnorm(2 * m_age[i-1] - m_age[i-2], mprec)
    }
    for (i in 3:n_agef) {
      f_age[i]~dnorm(2 * f_age[i-1] - f_age[i-2], fprec)
    }
    m_age_mu <- mean(m_age[1:n_agem])
    f_age_mu <- mean(f_age[1:n_agef])

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
  sus_beta0_temp ~ dnorm(0, .01)
  sus_mix ~ dunif(-1, 1)
  sus_beta0 <- sus_beta0_temp * sus_mix
  sus_beta_sex ~ dnorm(0, .01)

  #Priors for Age and Period effects
  #Age effects
  for (k in 1:nknots_age) {
    b_age[k] ~ dnorm(0, tau_age)
  }
  tau_age ~ dgamma(1, 1)

  for (t in 1:nT_age_surv) {
    age_effect_surv[t] <- inprod(b_age[1:nknots_age],
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
  inf_beta0_temp ~ dnorm(0, .01)
  inf_mix ~ dunif(-1, 1)
  inf_beta0 <- inf_beta0_temp * inf_mix  
  # inf_beta0 ~ T(dnorm(0,.01),0,Inf)
  inf_beta_sex ~ dnorm(0, .01)

  ##########################################################################################
  ###
  ### Infection hazard for different age classes, years(periods), and sexes 
  ###
  ##########################################################################################

  # for(j in 1:n_period_foi){
  #   for(i in 1:n_agef_foi){
  #     llambda_foi[i,j,2] <- f_age_foi[i] + f_period_foi[j]
  #   }
  #   for(i in 1:n_agem_foi){
  #     llambda_foi[i,j,1] <- m_age_foi[i] + m_period_foi[j]
  #   }
  # }

  # for(i in 1:nT_age_surv) {
  #   for(j in 1:nT_period_surv) {
  #     llambda_sus[i, j, 1] <- age_effect_surv[t] + period_effect_surv[t] + sus_beta0
  #     llambda_sus[i, j, 2] <- age_effect_surv[t] + period_effect_surv[t] + sus_beta0 + beta_sex
  #     llambda_inf[i, j, 1] <- age_effect_surv[t] + period_effect_surv[t] + inf_beta0
  #     llambda_inf[i, j, 2] <- age_effect_surv[t] + period_effect_surv[t] + inf_beta0 + inf_beta_sex
  #   }
  # }



})#end model statement



#######################################
### Data for Model Fitting
#######################################

nimData <- list(sus_censor = d_fit_sus$censor,
                Z_period = Z_period,
                Z_age = Z_age,
                sus_left_age = d_fit_sus$left_age,
                sus_right_age = d_fit_sus$right_age,
                sus_sex = d_fit_sus$sex,
                # inf_Z_period =  inf_Z_period,
                # inf_Z_age = inf_Z_age,
                icap_left_age = d_fit_icap$left_age,
                icap_right_age = d_fit_icap$right_age,
                icap_sex = d_fit_icap$sex,
                icap_censor = d_fit_icap$censor,
                rec_cens_censor = d_fit_rec_cens$censor,
                rec_cens_left_age = d_fit_rec_cens$left_age,
                rec_cens_left_period = d_fit_rec_cens$left_period,
                rec_cens_recap_age = d_fit_rec_cens$ageweek_recap,
                rec_cens_recap_period = d_fit_rec_cens$periodweek_recap,
                rec_cens_sex = d_fit_rec_cens$sex,
                rec_censor = d_fit_rec$censor,
                rec_left_age = d_fit_rec$left_age,
                rec_right_age = d_fit_rec$right_age,
                rec_left_period = d_fit_rec$left_period,
                rec_sex = d_fit_rec$sex,
                rec_censonly_censor = d_fit_rec_censonly$censor,
                rec_censonly_censor_inf = d_fit_rec_censonly$censor,
                rec_censonly_left_age = d_fit_rec_censonly$left_age,
                rec_censonly_right_age = d_fit_rec_censonly$right_age,
                rec_censonly_left_period = d_fit_rec_censonly$left_period,
                rec_censonly_recap_age = d_fit_rec_censonly$ageweek_recap,
                rec_censonly_recap_period = d_fit_rec_censonly$periodweek_recap,
                rec_censonly_sex = d_fit_rec_censonly$sex,
                idead_cens_censor = d_fit_idead_cens$censor,
                idead_cens_left_age = d_fit_idead_cens$left_age,
                idead_cens_right_age = d_fit_idead_cens$right_age,
                idead_cens_sex = d_fit_idead_cens$sex,
                idead_cens_left_period = d_fit_idead_cens$left_period,
                idead_cens_right_period = d_fit_idead_cens$right_period,
                idead_censor = d_fit_idead$censor,
                idead_left_age = idead_left_age,
                idead_right_age = d_fit_idead$right_age,
                idead_left_period = d_fit_idead$left_period,
                idead_sex = d_fit_idead$sex,
                teststatus_col = d_fit_col$censor,
                left_age_col = d_fit_col$left_age,
                right_age_col = d_fit_col$right_age,
                sexfoi_col = d_fit_col$sex,
                teststatus_inf = d_fit_inf$censor,
                left_age_inf = d_fit_inf$left_age,
                sexfoi_inf = d_fit_inf$sex,
                num_period = num_period,
                adj_period = adj_period,
                weights_period = weights_period,
                # snum = num_sp,
                # sadj = adj_sp,
                # sweights = weights_sp
                teststatus_hunt = cwd_df$teststatus,
                ageweeks_hunt = cwd_df$ageweeks,
                sexfoi_hunt = cwd_df$sex
                )


#######################################
### Constants for MCMC
#######################################

nimConsts <- list(sus_records = n_fit_sus,
                  icap_records = n_fit_icap,
                  nknots_age = nknots_age,
                  nknots_period = nknots_period,
                  sus_age2date = sus_age2date,
                  # inf_nknots_age = inf_nknots_age,
                  # inf_nknots_period = inf_nknots_period,
                  icap_age2date = icap_age2date,
                  n_period_lookup_icap =n_period_lookup_icap,
                  period_lookup_icap = period_lookup_icap,
                  n_icap_cens_post = n_icap_cens_post,
                  n_icap_post = n_icap_post,
                  n_icap_cens_pre = n_icap_cens_pre,
                  icap_post_la_indx = icap_post_la_indx,
                  icap_cens_post_age2date = icap_cens_post_age2date,
                  age_pre_indx = age_pre_indx,
                  rec_la_indx = rec_la_indx,
                  n_fit_rec = n_fit_rec,
                  n_fit_rec_cens = n_fit_rec_cens,
                  rec_cens_age2date = rec_cens_age2date,
                  rec_censonly_age2date = rec_censonly_age2date,
                  idead_cens_age2date = idead_cens_age2date,
                  n_fit_idead_cens = n_fit_idead_cens,
                  idead_indx = idead_indx,
                  n_fit_idead = n_fit_idead,
                  age_week_indx = age_week_indx,
                  period_week_indx = period_week_indx_col,
                  nT_age = nT_age,
                  nT_period = nT_period,
                  n_age = n_age,
                  n_agef = n_agef,
                  n_agem = n_agem,
                  n_period = n_period,
                  n_fit_col = n_fit_col,
                  n_fit_inf = n_fit_inf,
                  n_fit_hunt = n_fit_hunt,
                  n_age_lookup_col = n_age_lookup_col,
                  age_lookup_col_f = age_lookup_col_f,
                  age_lookup_col_m = age_lookup_col_m,
                  period_lookup_col = period_lookup_col,
                  n_period_lookup_col = n_period_lookup_col,
                  n_age_lookup_col_inf = n_age_lookup_col_inf,
                  age_lookup_col_inf_m = age_lookup_col_inf_m,
                  age_lookup_col_inf_f = age_lookup_col_inf_f,
                  period_lookup_col_inf = period_lookup_col_inf,
                  n_period_lookup_col_inf = n_period_lookup_col_inf,
                  n_adj_period = n_adj_period,
                  # n_adj_sp = n_adj_sp,
                  n_sect = n_sect,
                  sect_col = d_fit_col$sect,
                  age2date_col = age2date_col,
                  sect_inf_col = d_fit_inf$sect,
                  birth_week_inf = d_fit_inf$birth_week,
                  # scale = scale,
                  space = rep(0,n_sect),
                  n_age_lookup = n_age_lookup,
                  age_lookup_f = age_lookup_f,
                  age_lookup_m = age_lookup_m,
                  period_lookup_hunt = period_lookup_hunt,
                  n_period_lookup_hunt = n_period_lookup_hunt,
                  birthweek_hunt = cwd_df$birthweek,
                  sect_hunt = sect
                  )


#######################################
### Initial Values for MCMC
#######################################
initsFun <- function()list(tau_period = runif(1, .1, 1),
                          tau_age = runif(1, .1, .4),
                          sus_beta0_temp = rnorm(1, -6, 0.0001),
                          sus_mix = 1,
                          sus_beta_sex = rnorm(1,0,.1),
                          b_period = rnorm(nknots_period) * 10^-4,
                          b_age = rnorm(nknots_age) * 10^-4,
                          # beta_sex = rnorm(1, 0, .01),
                          # inf_tau_period = runif(1, .01, 1),
                          # inf_tau_age = runif(1, .1, 1),
                          inf_beta0_temp = rnorm(1, -3, 0.0001),
                          inf_mix = 1,
                          # inf_b_period = rnorm(inf_nknots_period) * 10^-4,
                          # inf_b_age = rnorm(inf_nknots_age) * 10^-4,
                          inf_beta_sex = rnorm(1, 0, .01),
                          rec_cens_age_add = rec_cens_age_add_init,
                          rec_censonly_age_add = rec_censonly_age_add_init,
                          idead_age_add = idead_age_add_init,
                          idead_left_age = idead_left_age_init,
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
                          m_period = seq(-2, 2, length = n_period),
                          f_period = seq(-2, 2, length = n_period),
                          m_age = seq(-6, -4, length = n_agem),
                          f_age = seq(-7, -5, length = n_agef),
                          llambda_foi = array(0,c(n_age,n_period,2))
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
              # "sus_beta0",
              # "sus_beta_sex",
              # "age_effect",
              # "period_effect",
              # "tau_period",
              # "tau_age",
              # "inf_beta0",
              # "inf_beta_sex",
              # "inf_age_effect",
              # "inf_period_effect",
              # "inf_tau_period",
              # "inf_tau_age",
              "m_age_mu",
              "f_age_mu",
              # "hprec",
              # "sprec",
              "mprec",
              "fprec",
              "tmprec",
              # "tmsd",
              "tfprec",
              # "tfsd",
              "f_age",
              "m_age",
              "f_period",
              "m_period"#,
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
                  niter = 2000,
                  nburnin = 1000,
                  nchains = 3,
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