#######################################################################
###
###   Likelihoods for each kind of data listed below
###
#######################################################################

# d_fit_hunt_neg
# d_fit_hunt_pos
# d_fit_sus_cens_posttest
# d_fit_sus_cens_postno
# d_fit_sus_mort_posttest
# d_fit_sus_mort_postno 
# d_fit_icap_cens
# d_fit_icap_mort
# d_fit_rec_neg_mort
# d_fit_rec_neg_cens
# d_fit_rec_pos_mort
# d_fit_rec_pos_cens
# d_fit_idead
# d_fit_endlive

#######################################################################
###
###   User defined distribution for likelihood for
###   infected harvest deer
###
###   d_fit_hunt_pos
###   Overleaf Equation 3
###
#######################################################################


# dInfHarvest <- nimble::nimbleFunction(
#     run = function(
#         ### argument type declarations
#         x = integer(0),
#         a = integer(0), #age (weeks) at harvest
#         sex = integer(0),
#         age2date = integer(0),
#         beta_sex = double(0),
#         beta0_sus = double(0),
#         beta0_inf = double(0),
#         age_effect_surv = double(1),
#         period_effect_surv = double(1),
#         f_age_foi = double(1),
#         m_age_foi = double(1),
#         age_lookup_f = double(1),
#         age_lookup_m = double(1),
#         period_lookup = double(1),
#         f_period_foi = double(1),
#         m_period_foi = double(1),
#         space = double(0),
#         log = double(0)
#         ) {

#     lam_foi <- nimNumeric(a)
#     lam_sus <- nimNumeric(a - 1)
#     lam_inf <- nimNumeric(a)
#     lik_temp <- nimNumeric(a)
#     indx_foi_age_f <- nimNumeric(a)
#     indx_foi_age_m <- nimNumeric(a)
#     indx_foi_period <- nimNumeric(a)
    
#     #############################################
#     # preliminary hazards for the likelihood
#     #############################################
#     indx_foi_age_f[1:a] <- age_lookup_f[1:a]
#     indx_foi_age_m[1:a] <- age_lookup_m[1:a]
#     indx_foi_period[1:a] <- period_lookup[(1 + age2date):(a + age2date)]

#     lam_foi[1:a] <- exp(rep(space, a) +
#                         sex * (f_age_foi[indx_foi_age_f[1:a]] + 
#                                f_period_foi[indx_foi_period[1:a]])+
#                         (1 - sex) * (m_age_foi[indx_foi_age_m[1:a]] + 
#                                      m_period_foi[indx_foi_period[1:a]])
#                     )
#     lam_sus[1:(a - 1)] <- exp(rep(beta0_sus, (a - 1)) +
#                             age_effect_surv[1:(a - 1)] +
#                             period_effect_surv[(1 + age2date):(a - 1 + age2date)] +
#                             rep(beta_sex * sex, (a - 1))
#                            )
#     lam_inf[1:a] <- exp(rep(beta0_inf, a) +
#                         age_effect_surv[1:a] +
#                         period_effect_surv[(1 + age2date):(a + age2date)] +
#                         rep(beta_sex * sex, a)
#                         )

#     #######################################
#     ### calculating the joint likelihood
#     #######################################

#     lik_temp[1] <- lam_foi[1] * exp(-sum(lam_inf[1:(a - 1)]))

#     for (k in 2:(a - 1)) {
#         lik_temp[k] <- lam_foi[k] *
#                exp(-sum(lam_sus[1:(k - 1)] + lam_foi[1:(k - 1)])) *
#                exp(-sum(lam_inf[k:(a - 1)]))
#     }
#     lik_temp[a] <- lam_foi[a] * exp(-sum(lam_foi[1:(a - 1)] +
#                                          lam_sus[1:(a - 1)]))
#     lik <- lam_inf[a] * sum(lik_temp[1:a])
#     llik <- log(lik)
#     returnType(double(0))
#     if(log) return(llik) else return(exp(llik))    ## return log-likelihood
#   })


# nimble::registerDistributions(list(
#     dInfHarvest = list(
#         BUGSdist = 'dInfHarvest(a,sex,age2date,beta_sex,beta0_sus,beta0_inf,age_effect_surv,period_effect_surv,f_age_foi,m_age_foi,age_lookup_f,age_lookup_m,f_period_foi,m_period_foi,period_lookup,space)',
#         types = c("a = integer(0)",
#                     "sex = integer(0)",
#                     "age2date = integer(0)",
#                     "beta_sex = double(0)",
#                     "beta0_sus = double(0)",
#                     "beta0_inf = double(0)",
#                     "age_effect_surv = double(1)",
#                     "period_effect_surv = double(1)",
#                     "f_age_foi = double(1)",
#                     "m_age_foi = double(1)",
#                     "age_lookup_f = double(1)",
#                     "age_lookup_m = double(1)",
#                     "period_lookup = double(1)",
#                     "f_period_foi = double(1)",
#                     "m_period_foi = double(1)",
#                     "space = double(0)",
#                     "log = double(0)"
#                   ),
#         discrete = TRUE
#     )
# ))

# # for a user-defined distribution
# assign('dInfHarvest', dInfHarvest, envir = .GlobalEnv)

# dInfHarvest(
#         x = 1,
#         a = d_fit_hunt_pos$ageweeks[1], #age (weeks) at harvest
#         sex = d_fit_hunt_pos$sex[1],
#         age2date = d_fit_hunt_pos$birthweek[1]-1,
#         beta_sex = beta_sex,
#         beta0_sus = beta0_sus,
#         beta0_inf = beta0_inf,
#         age_effect_surv = age_effect_survival_test,
#         period_effect_surv = period_effect_survival_test,
#         f_age_foi = f_age_foi,
#         m_age_foi = m_age_foi,
#         age_lookup_f = age_lookup_f,
#         age_lookup_m = age_lookup_m,
#         period_lookup = period_lookup,
#         f_period_foi = f_period_foi,
#         m_period_foi = m_period_foi,
#         space = space_mn[1],
#         log = TRUE
#         )

# test=c()
#   for (i in 1:nrow(d_fit_hunt_pos)) {
#     test[i] <- dInfHarvest(x = 1,
#                   a = d_fit_hunt_pos$ageweeks[i], #age (weeks) at harvest
#                   sex = d_fit_hunt_pos$sex[i],
#                   age2date = d_fit_hunt_pos$birthweek[i],
#                   beta_sex = beta_sex,
#                   beta0_sus = beta0_sus,
#                   beta0_inf = beta0_inf,
#                   age_effect_surv = age_effect_survival_test,
#                   period_effect_surv = period_effect_survival_test,
#                   f_age_foi = f_age_foi,
#                   m_age_foi = m_age_foi,
#                   age_lookup_f = age_lookup_f,
#                   age_lookup_m = age_lookup_m,
#                   period_lookup = period_lookup,
#                   f_period_foi = f_period_foi,
#                   m_period_foi = m_period_foi,
#                   space = 0,
#                   log = TRUE
#                   )
#   }
# which(is.na(test))

#######################################################################
###
###   User defined distribution for likelihood for
###   infected harvest deer, revised for multiple deer
###
###   d_fit_hunt_pos
###   Overleaf Equation 3
###
#######################################################################


dInfHarvest <- nimble::nimbleFunction(
    run = function(
        ### argument type declarations
        x = double(0),
        n_samples = integer(0), # number of samples in dataset
		a = double(1), #age (weeks) at harvest
        sex = double(1),
        age2date = double(1),
        beta_sex = double(0),
        beta0_sus = double(0),
        beta0_inf = double(0),
        age_effect_surv = double(1),
        period_effect_surv = double(1),
        f_age_foi = double(1),
        m_age_foi = double(1),
        age_lookup_f = double(1),
        age_lookup_m = double(1),
        period_lookup = double(1),
        f_period_foi = double(1),
        m_period_foi = double(1),
        space = double(1),
        sect = double(1),
        log = double(0)
        ) {

	# start the loop through individuals
	sumllik <- 0
	for(i in 1:n_samples) {
		#intitialize vectors
		lam_foi <- nimNumeric(a[i])
		lam_sus <- nimNumeric(a[i] - 1)
		lam_inf <- nimNumeric(a[i])
		lik_temp <- nimNumeric(a[i])
		indx_foi_age_f <- nimNumeric(a[i])
		indx_foi_age_m <- nimNumeric(a[i])
		indx_foi_period <- nimNumeric(a[i])
		
		#############################################
		# preliminary hazards for the likelihood
		#############################################
		indx_foi_age_f[1:a[i]] <- age_lookup_f[1:a[i]]
		indx_foi_age_m[1:a[i]] <- age_lookup_m[1:a[i]]
		indx_foi_period[1:a[i]] <- period_lookup[(1 + age2date[i]):(a[i] + age2date[i])]
	
		lam_foi[1:a[i]] <- exp(rep(space[sect[i]], a[i]) +
							sex[i] * (f_age_foi[indx_foi_age_f[1:a[i]]] + 
								f_period_foi[indx_foi_period[1:a[i]]])+
					  (1 - sex[i]) * (m_age_foi[indx_foi_age_m[1:a[i]]] + 
								m_period_foi[indx_foi_period[1:a[i]]])
						)
		lam_sus[1:(a[i] - 1)] <- exp(rep(beta0_sus, (a[i] - 1)) +
								age_effect_surv[1:(a[i] - 1)] +
								period_effect_surv[(1 + age2date[i]):(a[i] - 1 + age2date[i])] +
								rep(beta_sex * sex[i], (a[i] - 1))
							)
		lam_inf[1:a[i]] <- exp(rep(beta0_inf, a[i]) +
							age_effect_surv[1:a[i]] +
							period_effect_surv[(1 + age2date[i]):(a[i] + age2date[i])] +
							rep(beta_sex * sex[i], a[i])
							)
	
		#######################################
		### calculating the joint likelihood
		#######################################
	
		lik_temp[1] <- lam_foi[1] * exp(-sum(lam_inf[1:(a[i] - 1)]))
	
		for (j in 2:(a[i] - 1)) {
			lik_temp[j] <- lam_foi[j] *
				exp(-sum(lam_sus[1:(j - 1)] + lam_foi[1:(j - 1)])) *
				exp(-sum(lam_inf[j:(a[i] - 1)]))
		}
		lik_temp[a[i]] <- lam_foi[a[i]] * exp(-sum(lam_foi[1:(a[i] - 1)] +	lam_sus[1:(a[i] - 1)]))
		
		lik <- lam_inf[a[i]] * sum(lik_temp[1:a[i]])
		sumllik <- sumllik + log(lik)
	}
    returnType(double(0))
    if(log) return(sumllik) else return(exp(sumllik))    ## return log-likelihood
  })


nimble::registerDistributions(list(
    dInfHarvest = list(
        BUGSdist = 'dInfHarvest(n_samples,a,sex,age2date,beta_sex,beta0_sus,beta0_inf,age_effect_surv,period_effect_surv,f_age_foi,m_age_foi,age_lookup_f,age_lookup_m,f_period_foi,m_period_foi,period_lookup,space,sect)',
        types = c("value=double(0)",
					"n_samples = integer(0)",
				    "a = double(1)",
                    "sex = double(1)",
                    "age2date = double(1)",
                    "beta_sex = double(0)",
                    "beta0_sus = double(0)",
                    "beta0_inf = double(0)",
                    "age_effect_surv = double(1)",
                    "period_effect_surv = double(1)",
                    "f_age_foi = double(1)",
                    "m_age_foi = double(1)",
                    "age_lookup_f = double(1)",
                    "age_lookup_m = double(1)",
                    "period_lookup = double(1)",
                    "f_period_foi = double(1)",
                    "m_period_foi = double(1)",
                    "space = double(1)",
                    "sect = double(1)",
                    "log = double(0)"
                  ),
        discrete = TRUE
    )
))

# for a user-defined distribution
assign('dInfHarvest', dInfHarvest, envir = .GlobalEnv)

start <- Sys.time()
test <- dInfHarvest(
        x = rep(1,nrow(d_fit_hunt_pos)),
		  n_samples = nrow(d_fit_hunt_pos),
        a = d_fit_hunt_pos$ageweeks, #age (weeks) at harvest
        sex = d_fit_hunt_pos$sex,
        age2date = d_fit_hunt_pos$birthweek-1,
        beta_sex = beta_sex,
        beta0_sus = beta0_sus,
        beta0_inf = beta0_inf,
        age_effect_surv = age_effect_survival_test,
        period_effect_surv = period_effect_survival_test,
        f_age_foi = f_age_foi,
        m_age_foi = m_age_foi,
        age_lookup_f = age_lookup_f,
        age_lookup_m = age_lookup_m,
        period_lookup = period_lookup,
        f_period_foi = f_period_foi,
        m_period_foi = m_period_foi,
        space = rep(0,nrow(d_fit_hunt_pos)),
        sect = d_fit_hunt_pos$sect,
        log = TRUE
        )
(end<- Sys.time()-start)
test



#######################################################################
###
###   User defined distribution for likelihood for
###   uninfected harvest deer
###
###   d_fit_hunt_neg
###
###   Overleaf Equation (5)
###
#######################################################################

# dSusHarvest <- nimble::nimbleFunction(
#     run = function(
#         ### argument type declarations
#         x = integer(0),
#         a = integer(0),
#         sex = integer(0),
#         age2date = integer(0),
#         beta_sex = double(0),
#         beta0_sus = double(0),
#         age_effect_surv = double(1),
#         period_effect_surv = double(1),
#         f_age_foi = double(1),
#         m_age_foi = double(1),
#         age_lookup_f = double(1),
#         age_lookup_m = double(1),
#         period_lookup=double(1),
#         f_period_foi=double(1),
#         m_period_foi=double(1),
#         space = double(0),
#         log = double(0)
#         ) {
#     lik <- 0 #intialize likelihood
#     llik <- 0 #intialize log-likelihood
#     lam_foi <- nimNumeric(a)
#     lam_sus <- nimNumeric(a)
#     lam_temp <- nimNumeric(a)
#     indx_foi_age_f <- nimNumeric(a)
#     indx_foi_age_m <- nimNumeric(a)
#     indx_foi_period <- nimNumeric(a)
#     # lam_temp <- nimNumeric(a-1)
    
#     #############################################
#     # preliminary hazards for the likelihood
#     #############################################
#     indx_foi_age_f <- age_lookup_f[1:a]
#     indx_foi_age_m <- age_lookup_m[1:a]
#     indx_foi_period <- period_lookup[(1 + age2date):(a + age2date)]

#     lam_foi[1:a] <- exp(rep(space, a) +
#                         sex * (f_age_foi[indx_foi_age_f[1:a]] +
#                                f_period_foi[indx_foi_period[1:a]]) +
#                         (1 - sex) * (m_age_foi[indx_foi_age_m[1:a]] +
#                                      m_period_foi[indx_foi_period[1:a]])
#             )

#     lam_sus[1:a] <- exp(rep(beta0_sus, a) +
#                             age_effect_surv[1:a] +
#                             period_effect_surv[(1 + age2date):(a + age2date)] +
#                             rep(beta_sex * sex, a)
#                             )

#     #######################################
#     ###
#     ### calculating the joint likelihood
#     ###
#     #######################################

#     lik <- exp(-(sum(lam_foi[1:(a - 1)] + 
#                      lam_sus[1:(a - 1)]) - lam_foi[a])) * lam_sus[a]

#     llik <- log(lik)

#     returnType(double(0))
#     if(log) return(llik) else return(exp(llik))    ## return log-likelihood
#   })


# nimble::registerDistributions(list(
#     dSusHarvest = list(
#         BUGSdist = 'dSusHarvest(a,sex,age2date,beta_sex,beta0_sus,age_effect_surv,period_effect_surv,f_age_foi,m_age_foi,age_lookup_f,age_lookup_m,f_period_foi,m_period_foi,period_lookup,space)',
#         types = c("a = integer(0)",
#                   "sex = integer(0)",
#                   "age2date = integer(0)",
#                   "beta_sex = double(0)",
#                   "beta0_sus = double(0)",
#                   "age_effect_surv = double(1)",
#                   "period_effect_surv = double(1)",
#                   "f_age_foi = double(1)",
#                   "m_age_foi = double(1)",
#                   "age_lookup_f = double(1)",
#                   "age_lookup_m = double(1)",
#                   "period_lookup=double(1)",
#                   "f_period_foi=double(1)",
#                   "m_period_foi=double(1)",
#                   "space = double(0)",
#                   "log = double(0)"
#                   ),
#         discrete = TRUE
#     )
# ))

# # for a user-defined distribution
# assign('dSusHarvest', dSusHarvest, envir = .GlobalEnv)

# i=1
# dSusHarvest(
#         x = 1,
#         a = d_fit_hunt_neg$ageweeks[i], #age (weeks) at harvest
#         sex = d_fit_hunt_neg$sex[i],
#         age2date = d_fit_hunt_neg$birthweek[i]-1,
#         beta_sex = beta_sex,
#         beta0_sus = beta0_sus,
#         age_effect_surv = age_effect_survival_test,
#         period_effect_surv = period_effect_survival_test,
#         f_age_foi = f_age_foi,
#         m_age_foi = m_age_foi,
#         age_lookup_f = age_lookup_f,
#         age_lookup_m = age_lookup_m,
#         period_lookup = period_lookup,
#         f_period_foi = f_period_foi,
#         m_period_foi = m_period_foi,
#         space = space_mn[1],
#         log = TRUE
#         )

# test=c()
# for (i in 1:nrow(d_fit_hunt_neg)) {
#         test[i] <- dSusHarvest(
#         x = 1,
#         a = d_fit_hunt_neg$ageweeks[i], #age (weeks) at harvest
#         sex = d_fit_hunt_neg$sex[i],
#         age2date = d_fit_hunt_neg$birthweek[i] - 1,
#         beta_sex = beta_sex,
#         beta0_sus = beta0_sus,
#         age_effect_surv = age_effect_survival_test,
#         period_effect_surv = period_effect_survival_test,
#         f_age_foi = f_age_foi,
#         m_age_foi = m_age_foi,
#         age_lookup_f = age_lookup_f,
#         age_lookup_m = age_lookup_m,
#         period_lookup = period_lookup,
#         f_period_foi = f_period_foi,
#         m_period_foi = m_period_foi,
#         space = space_mn[1],
#         log = TRUE
#         )
# }
# which(is.na(test))


#######################################################################
###
###   User defined distribution for likelihood for
###   uninfected harvest deer
###
###   d_fit_hunt_neg
###
###   Overleaf Equation (5)
###
#######################################################################

dSusHarvest <- nimble::nimbleFunction(
    run = function(
        ### argument type declarations
        x = double(0),
        n_samples = double(0), # number of samples in dataset
		a = double(1), #age (weeks) at harvest
        sex = double(1),
        age2date = double(1),
        beta_sex = double(0),
        beta0_sus = double(0),
        age_effect_surv = double(1),
        period_effect_surv = double(1),
        f_age_foi = double(1),
        m_age_foi = double(1),
        age_lookup_f = double(1),
        age_lookup_m = double(1),
        period_lookup=double(1),
        f_period_foi=double(1),
        m_period_foi=double(1),
        space = double(1),
        sect = double(1),
        log = double(0)
        ) {

	# start the loop through individuals	
	sumllik <- 0
    # sumlliktemp <- nimNumeric(n_samples)
	for(i in 1:n_samples){
		#intitialize vectors
		lam_foi <- nimNumeric(a[i])
		lam_sus <- nimNumeric(a[i])
		lik_temp <- nimNumeric(a[i])
		indx_foi_age_f <- nimNumeric(a[i])
		indx_foi_age_m <- nimNumeric(a[i])
		indx_foi_period <- nimNumeric(a[i])
		
		#############################################
		# preliminary hazards for the likelihood
		#############################################
		indx_foi_age_f[1:a[i]] <- age_lookup_f[1:a[i]]
		indx_foi_age_m[1:a[i]] <- age_lookup_m[1:a[i]]
		indx_foi_period[1:a[i]] <- period_lookup[(1 + age2date[i]):(a[i] + age2date[i])]

		lam_foi[1:a[i]] <- exp(rep(space[sect[i]], a[i]) +
							sex[i] * (f_age_foi[indx_foi_age_f[1:a[i]]] +
								f_period_foi[indx_foi_period[1:a[i]]]) +
							(1 - sex[i]) * (m_age_foi[indx_foi_age_m[1:a[i]]] +
										m_period_foi[indx_foi_period[1:a[i]]])
				)
	
		lam_sus[1:a[i]] <- exp(rep(beta0_sus, a[i]) +
								age_effect_surv[1:a[i]] +
								period_effect_surv[(1 + age2date[i]):(a[i] + age2date[i])] +
								rep(beta_sex * sex[i], a[i])
								)

		#######################################
		###
		### calculating the joint likelihood
		###
		#######################################
	
		lik <- exp(-(sum(lam_foi[1:(a[i] - 1)] + 
						lam_sus[1:(a[i] - 1)]) - lam_foi[a[i]])) * lam_sus[a[i]]
	
		sumllik <- sumllik + log(lik)
        # sumlliktemp[i] <- sumllik + log(lik)
	}
   
	returnType(double(0))
	if(log) return(sumllik) else return(exp(sumllik))    ## return log-likelihood
  })


nimble::registerDistributions(list(
    dSusHarvest = list(
        BUGSdist = 'dSusHarvest(n_samples,a,sex,age2date,beta_sex,beta0_sus,age_effect_surv,period_effect_surv,f_age_foi,m_age_foi,age_lookup_f,age_lookup_m,f_period_foi,m_period_foi,period_lookup,space,sect)',
        types = c("value=double(0)",
				  "a = double(1)",
				  "n_samples = double(0)",
                  "sex = double(1)",
                  "age2date = double(1)",
                  "beta_sex = double(0)",
                  "beta0_sus = double(0)",
                  "age_effect_surv = double(1)",
                  "period_effect_surv = double(1)",
                  "f_age_foi = double(1)",
                  "m_age_foi = double(1)",
                  "age_lookup_f = double(1)",
                  "age_lookup_m = double(1)",
                  "period_lookup=double(1)",
                  "f_period_foi=double(1)",
                  "m_period_foi=double(1)",
                  "space = double(1)",
                  "sect = double(1)",
                  "log = double(0)"
                  ),
        discrete = TRUE
    )
))

assign('dSusHarvest', dSusHarvest, envir = .GlobalEnv)

#  dSusHarvest(
#         x = rep(1,nrow(d_fit_hunt_neg)),
# 		n_samples = nrow(d_fit_hunt_neg),
#         a = d_fit_hunt_neg$ageweeks, #age (weeks) at harvest
#         sex = d_fit_hunt_neg$sex,
#         age2date = d_fit_hunt_neg$birthweek-1,
#         beta_sex = beta_sex,
#         beta0_sus = beta0_sus,
#         age_effect_surv = age_effect_survival_test,
#         period_effect_surv = period_effect_survival_test,
#         f_age_foi = f_age_foi,
#         m_age_foi = m_age_foi,
#         age_lookup_f = age_lookup_f,
#         age_lookup_m = age_lookup_m,
#         period_lookup = period_lookup,
#         f_period_foi = f_period_foi,
#         m_period_foi = m_period_foi,
#         space = rep(0,n_sect),
#         sect = sect_hunt_neg,
#         log = TRUE
#         )

#######################################################################
#######################################################################
#######################################################################
###
### Likelihoods for radiomarked deer
###
#######################################################################
#######################################################################
#######################################################################

#######################################################################
###
###   User defined distribution for likelihood for
###   Uninfected radio-marked deer right censor:
###   Test neg at cap and censoring
###
###   d_fit_sus_cens_posttest
###   Overleaf Equation 7
###
#######################################################################

dSusCensTest <- nimble::nimbleFunction(
    run = function(
       ### argument type declarations
        x = integer(0),
        e = integer(0), #e, age of entry
        r = integer(0), #r, age of last known alive
        sex = integer(0),
        age2date = integer(0),
        beta_sex = double(0),
        beta0_sus = double(0),
        age_effect_surv = double(1),
        period_effect_surv = double(1),
        f_age_foi = double(1),
        m_age_foi = double(1),
        age_lookup_f = double(1),
        age_lookup_m = double(1),
        period_lookup = double(1),
        f_period_foi = double(1),
        m_period_foi = double(1),
        space = double(0),
        log = double()
        ) {

    lik <- 0 #intialize likelihood
    llik <- 0 #intialize log-likelihood
    lam_foi <- nimNumeric(r - 1)
    lam_sus <- nimNumeric(r - 1)
    lam_inf <- nimNumeric(r - 1)
    lik_temp <- nimNumeric(r - 1)
    indx_sus_period <- nimNumeric(r - 1)

    #############################################
    # preliminary hazards for the likelihood
    #############################################
    indx_sus_age <- e:(r - 1)
    n_indx_sus <- r - e
    indx_sus_period <- (e + age2date):(r - 1 + age2date)

    #survival hazard for susceptible deer
    lam_sus[e:(r - 1)] <- exp(
            rep(beta0_sus, n_indx_sus) +
            age_effect_surv[e:(r - 1)] +
            period_effect_surv[indx_sus_period[1:n_indx_sus]] +
            rep(beta_sex * sex, n_indx_sus)
    )

    #force of infection infection hazard
    indx_foi_age_f <- age_lookup_f[1:(r - 1)]
    indx_foi_age_m <- age_lookup_m[1:(r - 1)]
    indx_foi_period <- period_lookup[(1 + age2date):(r - 1 + age2date)]

    lam_foi[1:(r - 1)] <- exp(rep(space, r - 1) +
            sex * (f_age_foi[indx_foi_age_f[1:(r - 1)]] +
                   f_period_foi[indx_foi_period[1:(r - 1)]]) +
            (1 - sex) * (m_age_foi[indx_foi_age_m[1:(r - 1)]] +
                         m_period_foi[indx_foi_period[1:(r - 1)]])
    )
    #######################################
    ### calculating the joint likelihood
    #######################################

    lik <- exp(-sum(lam_sus[e:(r - 1)])) *
            exp(-sum(lam_foi[1:(r - 1)]))

    llik <- log(lik)

    returnType(double(0))
    if(log) return(llik) else return(exp(llik))    ## return log-likelihood
  })


nimble::registerDistributions(list(
    dSusCensTest = list(
        BUGSdist = 'dSusCensTest(e,r,sex,age2date,beta_sex,beta0_sus,age_effect_surv,period_effect_surv,f_age_foi,m_age_foi,age_lookup_f,age_lookup_m,f_period_foi,m_period_foi,period_lookup,space)',
        types = c("e = integer(0)",
                  "r = integer(0)",
                  "sex = integer(0)",
                  "age2date = integer(0)",
                  "beta_sex = double(0)",
                  "beta0_sus = double(0)",
                  "age_effect_surv = double(1)",
                  "period_effect_surv = double(1)",
                  "f_age_foi = double(1)",
                  "m_age_foi = double(1)",
                  "age_lookup_f = double(1)",
                  "age_lookup_m = double(1)",
                  "period_lookup=double(1)",
                  "f_period_foi=double(1)",
                  "m_period_foi=double(1)",
                  "space = double(0)",
                  "log = double(0)"
                  ),
        discrete = TRUE
    )
))


###for a user-defined distribution
assign('dSusCensTest', dSusCensTest, envir = .GlobalEnv)

# i = 1
# dSusCensTest(
#         x = 1,
#         e = d_fit_sus_cens_posttest$left_age_e[i],
#         r = d_fit_sus_cens_posttest$right_age_r[i],
#         sex = d_fit_sus_cens_posttest$sex[i],
#         age2date = sus_cens_posttest_age2date[i],
#         beta_sex = beta_sex,
#         beta0_sus = beta0_sus,
#         age_effect_surv = age_effect_survival_test,
#         period_effect_surv = period_effect_survival_test,
#         f_age_foi = f_age_foi,
#         m_age_foi = m_age_foi,
#         age_lookup_f = age_lookup_col_f,
#         age_lookup_m = age_lookup_col_m,
#         period_lookup = period_lookup,
#         f_period_foi = f_period_foi,
#         m_period_foi = m_period_foi,
#         space = 0,
#         log = TRUE
#         )
#######################################################################
###
###   User defined distribution for likelihood for
###   Uninfected radio-marked deer right censored:
###   Test neg at cap and censoring
###
###   d_fit_sus_cens_postno
###   d_fit_endlive
###
###   Overleaf Equation (9)
###
#######################################################################

dSusCensNo <- nimble::nimbleFunction(
    run = function(
       ### argument type declarations
        x = double(0),
        n_samples = integer(0), # number of samples in dataset
        e = double(1), #e, age of entry
        r = double(1), #r, age of last known alive
        sex = double(1),
		age2date = double(1),
        beta_sex = double(0),
        beta0_sus = double(0),
        beta0_inf = double(0),
        age_effect_surv = double(1),
        period_effect_surv = double(1),
        f_age_foi = double(1),
        m_age_foi = double(1),
        age_lookup_f = double(1),
        age_lookup_m = double(1),
        period_lookup = double(1),
        f_period_foi = double(1),
        m_period_foi = double(1),
        space = double(1),
        sect = double(1),
        log = double()
        ) {

    sumllik <- 0 #intialize log-likelihood
	for(i in 1:n_samples) {
		lam_foi <- nimNumeric(r[i] - 1)
		lam_sus <- nimNumeric(r[i] - 1)
		lam_inf <- nimNumeric(r[i] - 1)
		lik_temp <- nimNumeric(r[i] - 1)
		indx_sus_period <- nimNumeric(r[i] - 1)
		
		#############################################
		# preliminary hazards for the likelihood
		#############################################
		
		indx_sus_age <- e[i]:(r[i] - 1)
		n_indx_sus <- r[i] - e[i]
		indx_sus_period <- (e[i] + age2date[i]):(r[i] - 1 + age2date[i])
	
		#survival hazard for susceptible deer
		lam_sus[indx_sus_age] <- exp(
			rep(beta0_sus, n_indx_sus) +
				age_effect_surv[indx_sus_age] +
				period_effect_surv[indx_sus_period[1:n_indx_sus]] +
				rep(beta_sex * sex[i], n_indx_sus)
		)
		#survival hazard while infected
		lam_inf[indx_sus_age] <- exp(rep(beta0_inf, n_indx_sus) +
			age_effect_surv[indx_sus_age] +
			period_effect_surv[indx_sus_period] +
			rep(beta_sex * sex[i], n_indx_sus)
			)
	
		#force of infection infection hazard
		indx_foi_age_f <- age_lookup_f[1:(r[i] - 1)]
		indx_foi_age_m <- age_lookup_m[1:(r[i] - 1)]
		indx_foi_period <- period_lookup[(1 + age2date[i]):(r[i] - 1 + age2date[i])]
	
		lam_foi[1:(r[i] - 1)] <- exp(rep(space[sect[i]], r[i] - 1) +
				sex[i] * (f_age_foi[indx_foi_age_f[1:(r[i] - 1)]] +
					f_period_foi[indx_foi_period[1:(r[i] - 1)]]) +
				(1 - sex[i]) * (m_age_foi[indx_foi_age_m[1:(r[i] - 1)]] +
							m_period_foi[indx_foi_period[1:(r[i] - 1)]])
		)
		#######################################
		### calculating the joint likelihood
		#######################################
	
		lik_temp[(e[i] + 1)] <- lam_foi[(e[i] + 1)] * exp(-sum(lam_inf[(e[i] + 1):(r[i] - 1)]))
	    if(r[i]-e[i]>2){
			for(k in (e[i] + 2):(r[i] - 1)){
				lik_temp[k] <- lam_foi[k] *
							exp(-sum(lam_sus[(e[i] + 1):(k - 1)])) *
							exp(-sum(lam_foi[(e[i] + 1):(k - 1)])) *
							exp(-sum(lam_inf[k:(r[i] - 1)]))
			}
		}
		lik <- exp(-lam_sus[e[i]]) *
			exp(-sum(lam_foi[1:e[i]])) *
			sum(lik_temp[(e[i] + 1):(r[i] - 1)]) +
			exp(-sum(lam_foi[1:(r[i] - 1)])) *
			exp(-sum(lam_sus[e[i]:(r[i] - 1)]))

		sumllik <- sumllik + log(lik)
	}
    returnType(double(0))
    if(log) return(sumllik) else return(exp(sumllik))    ## return log-likelihood
  })

nimble::registerDistributions(list(
    dSusCensNo = list(
        BUGSdist = 'dSusCensNo(n_samples,e,r,sex,age2date,beta_sex,beta0_sus,beta0_inf,age_effect_surv,period_effect_surv,f_age_foi,m_age_foi,age_lookup_f,age_lookup_m,f_period_foi,m_period_foi,period_lookup,space,sect)',
        types = c("value=double(0)",
				  "n_samples = integer(0)",
                  "e = double(1)",
                  "r = double(1)",
                  "sex = double(1)",
                  "age2date = double(1)",
                  "beta_sex = double(0)",
                  "beta0_sus = double(0)",
                  "beta0_inf = double(0)",
                  "age_effect_surv = double(1)",
                  "period_effect_surv = double(1)",
                  "f_age_foi = double(1)",
                  "m_age_foi = double(1)",
                  "age_lookup_f = double(1)",
                  "age_lookup_m = double(1)",
                  "period_lookup=double(1)",
                  "f_period_foi=double(1)",
                  "m_period_foi=double(1)",
                  "space = double(1)",
                  "sect = double(1)",
                  "log = double()"
                  ),
        discrete = TRUE
    )
))

###for a user-defined distribution
assign('dSusCensNo', dSusCensNo, envir = .GlobalEnv)

start <- Sys.time()
test <- dSusCensNo(
        x = rep(1,nrow(d_fit_sus_cens_postno)),
		  n_samples = nrow(d_fit_sus_cens_postno),
        e = d_fit_sus_cens_postno$left_age_e,
        r = d_fit_sus_cens_postno$right_age_r,
        sex = d_fit_sus_cens_postno$sex,
        age2date = sus_cens_postno_age2date,
        beta_sex = beta_sex,
        beta0_sus = beta0_sus,
        beta0_inf = beta0_inf,
        age_effect_surv = age_effect_survival_test,
        period_effect_surv = period_effect_survival_test,
        f_age_foi = f_age_foi,
        m_age_foi = m_age_foi,
        age_lookup_f = age_lookup_col_f,
        age_lookup_m = age_lookup_col_m,
        period_lookup = period_lookup,
        f_period_foi = f_period_foi,
        m_period_foi = m_period_foi,
        space = rep(0,n_sect),
		  sect = d_fit_sus_cens_postno$sect,
        log = TRUE
        )
(end<- Sys.time()-start)
test

start <- Sys.time()
test2 <- dSusCensNo(
        x = rep(1,nrow(d_fit_endlive)),
		  n_samples = nrow(d_fit_endlive),
        e = d_fit_endlive$left_age_e,
        r = d_fit_endlive$right_age_r,
        sex = d_fit_endlive$sex,
        age2date = endlive_age2date,
        beta_sex = beta_sex,
        beta0_sus = beta0_sus,
        beta0_inf = beta0_inf,
        age_effect_surv = age_effect_survival_test,
        period_effect_surv = period_effect_survival_test,
        f_age_foi = f_age_foi,
        m_age_foi = m_age_foi,
        age_lookup_f = age_lookup_col_f,
        age_lookup_m = age_lookup_col_m,
        period_lookup = period_lookup,
        f_period_foi = f_period_foi,
        m_period_foi = m_period_foi,
        space = rep(0,n_sect),
		  sect = d_fit_endlive$sect,
        log = TRUE
        )
(end<- Sys.time()-start)
test2

#######################################################################
###
###   User defined distribution for likelihood for
###   uninfected radio-marked deer mortalities:
###   test neg at cap and tested mort
###
###   d_fit_sus_mort_posttest
###
###   Overleaf Equation (11)
###
#######################################################################


dSusMortTest <- nimble::nimbleFunction(
    run = function(
       ### argument type declarations
        x = double(0),
        n_samples = integer(0), # number of samples in dataset
        e = double(1), #e, age of entry
        r = double(1), #r, age of last known alive
        s = double(1), #s, age of known mortality
        sex = double(1),
        fast = double(1),
		age2date = double(1),
        beta_sex = double(0),
        beta0_sus = double(0),
        age_effect_surv = double(1),
        period_effect_surv = double(1),
        f_age_foi = double(1),
        m_age_foi = double(1),
        age_lookup_f = double(1),
        age_lookup_m = double(1),
        period_lookup = double(1),
        f_period_foi = double(1),
        m_period_foi = double(1),
        space = double(1),
        sect = double(1),
        log = double()
        ) {

    sumllik <- 0 #intialize log-likelihood
	for(i in 1:n_samples){
		lam_foi <- nimNumeric(s[i] - 1)
		lam_sus <- nimNumeric(s[i] - 1)
		indx_foi_age_f <- nimNumeric(s[i] - 1)
		indx_foi_age_m <- nimNumeric(s[i] - 1)
		indx_foi_period <- nimNumeric(s[i] - 1)
	
		#############################################
		# preliminary hazards for the likelihood
		#############################################
		n_indx_sus <- s[i] - e[i]
		#survival hazard for susceptible deer
		lam_sus[e[i]:(s[i] - 1)] <- exp(rep(beta0_sus, n_indx_sus) +
				age_effect_surv[e[i]:(s[i] - 1)] +
				period_effect_surv[(e[i] + age2date[i]):(s[i] - 1 + age2date[i])] +
				rep(beta_sex * sex[i], n_indx_sus)
		)
		#force of infection infection hazard
		indx_foi_age_f[1:(s[i] - 1)] <- age_lookup_f[1:(s[i] - 1)]
		indx_foi_age_m[1:(s[i] - 1)] <- age_lookup_m[1:(s[i] - 1)]
		indx_foi_period[1:(s[i] - 1)] <- period_lookup[(1 + age2date[i]):(s[i] - 1 + age2date[i])]
	
		lam_foi[1:(s[i] - 1)] <- exp(rep(space[sect[i]], s[i] - 1) +
				sex[i] * (f_age_foi[indx_foi_age_f[1:(s[i] - 1)]] +
					f_period_foi[indx_foi_period[1:(s[i] - 1)]]) +
				(1 - sex[i]) * (m_age_foi[indx_foi_age_m[1:(s[i] - 1)]] +
							m_period_foi[indx_foi_period[1:(s[i] - 1)]])
		)
	
		#######################################
		### calculating the joint likelihood
		#######################################
 #   	lik <- (1 - exp(-sum(lam_sus[r:(s - 1)]))) *
 #   	       exp(-sum(lam_sus[e:(r - 1)])) *
 #   	       exp(-sum(lam_foi[1:(s - 1)]))
#		if(e[i]<r[i]){
#				lik <- (1 - exp(-sum(lam_sus[r[i]:(s[i] - 1)]))) *
#					exp(-sum(lam_sus[e[i]:(r[i] - 1)])) *
#					exp(-sum(lam_foi[1:(s[i] - 1)]))
#		}else{
#				lik <- (1 - exp(-sum(lam_sus[r[i]:(s[i] - 1)]))) *
#					exp(-sum(lam_foi[1:(s[i] - 1)]))
#		}
			lik <- (1 - fast[i]) * ((1 - exp(-sum(lam_sus[r[i]:(s[i] - 1)]))) *
                 exp(-sum(lam_sus[e[i]:(r[i] - 1)])) *
                 exp(-sum(lam_foi[1:(s[i] - 1)]))) +
            fast[i] * ((1 - exp(-sum(lam_sus[r[i]:(s[i] - 1)]))) *
            exp(-sum(lam_foi[1:(s[i] - 1)])))
	
		sumllik <- sumllik + log(lik)
	}
	returnType(double(0))
	if(log) return(sumllik) else return(exp(sumllik))    ## return log-likelihood
  })


nimble::registerDistributions(list(
    dSusMortTest = list(
        BUGSdist = 'dSusMortTest(n_samples,e,r,s,sex,fast,age2date,beta_sex,beta0_sus,age_effect_surv,period_effect_surv,f_age_foi,m_age_foi,age_lookup_f,age_lookup_m,f_period_foi,m_period_foi,period_lookup,space,sect)',
        types = c("value=double(0)",
				  "a = double(1)",
				  "n_samples = integer(0)",
                  "e = double(1)",
                  "r = double(1)",
                  "s = double(1)",
                  "sex = double(1)",
				  "fast = double(1)",
                  "age2date = double(1)",
                  "beta_sex = double(0)",
                  "beta0_sus = double(0)",
                  "age_effect_surv = double(1)",
                  "period_effect_surv = double(1)",
                  "f_age_foi = double(1)",
                  "m_age_foi = double(1)",
                  "age_lookup_f = double(1)",
                  "age_lookup_m = double(1)",
                  "period_lookup=double(1)",
                  "f_period_foi=double(1)",
                  "m_period_foi=double(1)",
                  "space = double(1)",
                  "sect = double(1)",
                  "log = double(0)"
                  ),
        discrete = TRUE
    )
))

# ###for a user-defined distribution
assign('dSusMortTest', dSusMortTest, envir = .GlobalEnv)

# start <- Sys.time()
# test <- dSusMortTest(
#         x = rep(1,nrow(d_fit_sus_mort_posttest)),
#		  n_samples = nrow(d_fit_sus_mort_posttest),
#         e = d_fit_sus_mort_posttest$left_age_e,
#         r = d_fit_sus_mort_posttest$right_age_r,
#         s = d_fit_sus_mort_posttest$right_age_s,
#         sex = d_fit_sus_mort_posttest$sex,
#		  fast = d_fit_sus_mort_posttest$fast,
#         age2date = sus_mort_posttest_age2date,
#         beta_sex = beta_sex,
#         beta0_sus = beta0_sus,
#         age_effect_surv = age_effect_survival_test,
#         period_effect_surv = period_effect_survival_test,
#         f_age_foi = f_age_foi,
#         m_age_foi = m_age_foi,
#         age_lookup_f = age_lookup_col_f,
#         age_lookup_m = age_lookup_col_m,
#         period_lookup = period_lookup,
#         f_period_foi = f_period_foi,
#         m_period_foi = m_period_foi,
#         space = rep(0,n_sect),
#		  sect = sect_sus_mort_posttest,
#         log = TRUE
#         )
# (end<- Sys.time()-start)
# test

#######################################################################
###
###   User defined distribution for likelihood for
###   uninfected radio-marked deer mortalities:
###   test neg at cap and no test at mortality, no recap
###
###   d_fit_sus_mort_postno
###
###   Overleaf Equation (13)
###
#######################################################################

dSusMortNoTest <- nimble::nimbleFunction(
    run = function(
        ### argument type declarations
        x = integer(0),
        e = integer(0), #e, age of entry
        r = integer(0), #r, age of last known alive
        s = integer(0), #s, age of known mortality
        dn1 = integer(0), #right of last test negative
        sex = integer(0),
        age2date = integer(0),
        beta_sex = double(0),
        beta0_sus = double(0),
        beta0_inf = double(0),
        age_effect_surv = double(1),
        period_effect_surv = double(1),
        f_age_foi = double(1),
        m_age_foi = double(1),
        age_lookup_f = double(1),
        age_lookup_m = double(1),
        period_lookup=double(1),
        f_period_foi=double(1),
        m_period_foi=double(1),
        space = double(0),
        log = double(0)
        ) {
    
    lam_foi <- nimNumeric(s - 1)
    lam_sus <- nimNumeric(s - 1)
    lam_inf <- nimNumeric(s - 1)
    lik_temp <- nimNumeric(s - 1)

    #############################################
    # preliminary hazards for the likelihood
    #############################################

    #survival hazard for susceptible deer
    n_indx_sus <- length(e:(s - 1))
    lam_sus[e:(s - 1)] <- exp(rep(beta0_sus, n_indx_sus) +
                    age_effect_surv[e:(s - 1)] +
                    period_effect_surv[(e + age2date):(s - 1 + age2date)] +
                    rep(beta_sex * sex, n_indx_sus))

    #survival hazard while infected
    n_indx_inf <- length(e:(s - 1))
    lam_inf[e:(s - 1)] <- exp(rep(beta0_inf, n_indx_inf) +
        age_effect_surv[e:(s - 1)] +
        period_effect_surv[(e + age2date):(s - 1 + age2date)] +
        rep(beta_sex * sex, n_indx_inf)
        )

    #force of infection infection hazard
    lam_foi[1:(s-1)] <- exp(rep(space, s - 1) +
        sex * (f_age_foi[age_lookup_f[1:(s - 1)]] +
            f_period_foi[period_lookup[(1 + age2date):(s - 1 + age2date)]]) +
        (1 - sex) * (m_age_foi[age_lookup_m[1:(s - 1)]] +
            m_period_foi[period_lookup[(1 + age2date):(s - 1 + age2date)]])
        )

    #total probability of getting infected and dying before end of the study
    lik_temp[dn1 + 1] <- lam_foi[dn1 + 1] *
                      exp(-sum(lam_inf[(dn1 + 1):(r - 1)])) *
                      (1 - exp(-sum(lam_inf[r:(s - 1)])))

    if((r - dn1)>3){
        for (k in (dn1 + 2):(r - 2)) {
            lik_temp[k] <- lam_foi[k] *
                      exp(-sum(lam_foi[(dn1 + 1):(k - 1)])) *
                      exp(-sum(lam_sus[(dn1 + 1):(k - 1)])) *
                      exp(-sum(lam_inf[k:(r - 1)])) *
                      (1 - exp(-sum(lam_inf[r:(s - 1)])))
        }
    }
    for (k in (r - 1):(s - 2)) {
        lik_temp[k] <- lam_foi[k] *
                      exp(-sum(lam_foi[(dn1 + 1):(k - 1)])) *
                      exp(-sum(lam_sus[(dn1 + 1):(k - 1)])) *
                      (1 - exp(-sum(lam_inf[k:(s - 1)])))
    }
    lik_temp[(s - 1)] <- lam_foi[(s - 1)] *
                        exp(-sum(lam_foi[(dn1 + 1):((s - 1) - 1)])) *
                        exp(-sum(lam_sus[(dn1 + 1):((s - 1) - 1)])) *
                        (lam_inf[(s - 1)])
    lik <- exp(-sum(lam_sus[e:dn1])) *
           exp(-sum(lam_foi[1:dn1])) *
           sum(lik_temp[(dn1 + 1):(s - 1)]) +
           exp(-sum(lam_foi[1:(s - 1)])) *
           exp(-sum(lam_sus[e:(r - 1)])) *
           (1 - exp(-sum(lam_sus[r:(s - 1)])))

    llik <- log(lik)

    returnType(double(0))
    if(log) return(llik) else return(exp(llik))
  })

nimble::registerDistributions(list(
    dSusMortNoTest = list(
        BUGSdist = 'dSusMortNoTest(e,r,s,dn1,sex,age2date,beta_sex,beta0_sus,beta0_inf,age_effect_surv,period_effect_surv,f_age_foi,m_age_foi,age_lookup_f,age_lookup_m,f_period_foi,m_period_foi,period_lookup,space)',
        types = c("e = integer(0)",
                  "r = integer(0)",
                  "s = integer(0)",
                  "dn1 = integer(0)",
                  "sex = integer(0)",
                  "age2date = integer(0)",
                  "beta_sex = double(0)",
                  "beta0_inf = double(0)",
                  "beta0_sus = double(0)",
                  "age_effect_surv = double(1)",
                  "period_effect_surv = double(1)",
                  "f_age_foi = double(1)",
                  "m_age_foi = double(1)",
                  "age_lookup_f = double(1)",
                  "age_lookup_m = double(1)",
                  "period_lookup=double(1)",
                  "f_period_foi=double(1)",
                  "m_period_foi=double(1)",
                  "space = double(0)",
                  "log = double(0)"
                  ),
        discrete = TRUE
    )
))

###for a user-defined distribution
assign('dSusMortNoTest', dSusMortNoTest, envir = .GlobalEnv)

# i=6
# dSusMortNoTest(
#         x = 1,
#         e = d_fit_sus_mort_postno$left_age_e[i],
#         r = d_fit_sus_mort_postno$right_age_r[i],
#         s = d_fit_sus_mort_postno$right_age_s[i],
#         dn1 = d_fit_sus_mort_postno$left_age_e[i],
#         sex = d_fit_sus_mort_postno$sex[i],
#         age2date = sus_mort_postno_age2date[i],
#         beta_sex = beta_sex,
#         beta0_sus = beta0_sus,
#         beta0_inf = beta0_inf,
#         age_effect_surv = age_effect_survival_test,
#         period_effect_surv = period_effect_survival_test,
#         f_age_foi = f_age_foi,
#         m_age_foi = m_age_foi,
#         age_lookup_f = age_lookup_col_f,
#         age_lookup_m = age_lookup_col_m,
#         period_lookup = period_lookup,
#         f_period_foi = f_period_foi,
#         m_period_foi = m_period_foi,
#         space = 0,
#         log = TRUE
#         )

# test = c()
# for(i in 1:nrow(d_fit_sus_mort_postno)) {
#     test[i] <-  dSusMortNoTest(
#         x = 1,
#         e = d_fit_sus_mort_postno$left_age_e[i],
#         r = d_fit_sus_mort_postno$right_age_r[i],
#         s = d_fit_sus_mort_postno$right_age_s[i],
#         dn1 = d_fit_sus_mort_postno$left_age_e[i],
#         sex = d_fit_sus_mort_postno$sex[i],
#         age2date = sus_mort_postno_age2date[i],
#         beta_sex = beta_sex,
#         beta0_sus = beta0_sus,
#         beta0_inf = beta0_inf,
#         age_effect_surv = age_effect_survival_test,
#         period_effect_surv = period_effect_survival_test,
#         f_age_foi = f_age_foi,
#         m_age_foi = m_age_foi,
#         age_lookup_f = age_lookup_col_f,
#         age_lookup_m = age_lookup_col_m,
#         period_lookup = period_lookup,
#         f_period_foi = f_period_foi,
#         m_period_foi = m_period_foi,
#         space = 0,
#         log = TRUE
#         )
# }
# test

#######################################################################
###
###   User defined distribution for likelihood for
###   infected deer mortalities for radio marked deer that
###   enter the study as test positive at capture
###
###   d_fit_icap_cens
###
###   Overleaf Equation (15)
###
#######################################################################

dIcapCens <- nimble::nimbleFunction(
    run = function(
        ### argument type declarations
        x = integer(0),
        e = integer(0), #e, age of entry
        r = integer(0), #r, age of last known alive
        sex = integer(0),
        age2date = integer(0),
        beta_sex = double(0),
        beta0_sus = double(0),
        beta0_inf = double(0),
        age_effect_surv = double(1),
        period_effect_surv = double(1),
        f_age_foi = double(1),
        m_age_foi = double(1),
        age_lookup_f = double(1),
        age_lookup_m = double(1),
        period_lookup=double(1),
        f_period_foi=double(1),
        m_period_foi=double(1),
        space = double(0),
        log = double()
        ) {

    llik <- 0 #intialize log-likelihood
    lam_inf <- nimNumeric(r - 1)
    lam_foi <- nimNumeric(r - 1)
    lam_sus <- nimNumeric(r - 1)
    lik_temp <- nimNumeric(e - 1)
    
    #############################################
    # preliminary hazards for the likelihood
    #############################################
    n_indx_sus <- length(1:(e - 2))
    #survival hazard for susceptible deer
    lam_sus[1:(e - 2)] <- exp(rep(beta0_sus, n_indx_sus) +
                    age_effect_surv[1:(e - 2)] +
                    period_effect_surv[(1 + age2date):(e - 2 + age2date)] +
                    rep(beta_sex * sex, n_indx_sus))

    #survival hazard while infected
    n_indx_inf <- length(1:(r - 1))
    lam_inf[1:(r - 1)] <- exp(rep(beta0_inf, n_indx_inf) +
        age_effect_surv[1:(r - 1)] +
        period_effect_surv[(1 + age2date):(r - 1 + age2date)] +
        rep(beta_sex * sex, n_indx_inf)
        )

    #force of infection infection hazard
    lam_foi[1:(e-1)] <- exp(rep(space, e - 1) +
        sex * (f_age_foi[age_lookup_f[1:(e - 1)]] +
            f_period_foi[period_lookup[(1 + age2date):(e - 1 + age2date)]]) +
        (1 - sex) * (m_age_foi[age_lookup_m[1:(e - 1)]] +
            m_period_foi[period_lookup[(1 + age2date):(e - 1 + age2date)]])
        )
    
    #######################################
    ### calculating the joint likelihood
    #######################################
    lik_temp[1] <- lam_foi[1] * exp(-sum(lam_inf[1:(e - 1)]))
    for(k in 2:(e - 1)){
     lik_temp[k] <- lam_foi[k] *
               exp(-sum(lam_sus[1:(k - 1)])) *
               exp(-sum(lam_foi[1:(k - 1)])) *
               exp(-sum(lam_inf[k:(e - 1)]))
    }
       lik <- exp(-sum(lam_inf[e:(r - 1)])) * sum(lik_temp[1:(e - 1)])
    llik <- log(lik)
    returnType(double(0))
    if(log) return(llik) else return(exp(llik))    ## return log-likelihood
  })


nimble::registerDistributions(list(
    dIcapCens = list(
        BUGSdist = 'dIcapCens(e,r,sex,age2date,beta_sex,beta0_inf,beta0_sus,age_effect_surv,period_effect_surv,f_age_foi,m_age_foi,age_lookup_f,age_lookup_m,f_period_foi,m_period_foi,period_lookup,space)',
        types = c("e = integer(0)",
                  "r = integer(0)",
                  "sex = integer(0)",
                  "age2date = integer(0)",
                  "beta_sex = double(0)",
                  "beta0_inf = double(0)",
                  "beta0_sus = double(0)",
                  "age_effect_surv = double(1)",
                  "period_effect_surv = double(1)",
                  "f_age_foi = double(1)",
                  "m_age_foi = double(1)",
                  "age_lookup_f = double(1)",
                  "age_lookup_m = double(1)",
                  "period_lookup=double(1)",
                  "f_period_foi=double(1)",
                  "m_period_foi=double(1)",
                  "space = double(0)",
                  "log = double()"
                  ),
        discrete = TRUE
    )
))

# for a user-defined distribution
assign('dIcapCens', dIcapCens, envir = .GlobalEnv)

i <- 1
dIcapCens(
        x = 1,
        e = d_fit_icap_cens$left_age_e[i],
        r = d_fit_icap_cens$right_age_r[i],
        sex = d_fit_icap_cens$sex[i],
        age2date = icap_cens_age2date[i],
        beta_sex = beta_sex,
        beta0_inf = beta0_inf,
        beta0_sus = beta0_sus,
        age_effect_surv = age_effect_survival_test,
        period_effect_surv = period_effect_survival_test,
        f_age_foi = f_age_foi,
        m_age_foi = m_age_foi,
        age_lookup_f = age_lookup_col_f,
        age_lookup_m = age_lookup_col_m,
        period_lookup = period_lookup,
        f_period_foi = f_period_foi,
        m_period_foi = m_period_foi,
        space = 0,
        log = TRUE
        )
# test=c()
# for(i in 1:nrow(d_fit_icap_cens)){
#     test[i] <-  dIcapCens(
#         x = 1,
#         e = d_fit_icap_cens$left_age_e[i],
#         r = d_fit_icap_cens$right_age_r[i],
#         sex = d_fit_icap_cens$sex[i],
#         age2date = icap_cens_age2date[i],
#         beta_sex = beta_sex,
#         beta0_sus = beta0_sus,
#         beta0_inf = beta0_inf,
#         age_effect_surv = age_effect_survival_test,
#         period_effect_surv = period_effect_survival_test,
#         f_age_foi = f_age_foi,
#         m_age_foi = m_age_foi,
#         age_lookup_f = age_lookup_col_f,
#         age_lookup_m = age_lookup_col_m,
#         period_lookup = period_lookup,
#         f_period_foi = f_period_foi,
#         m_period_foi = m_period_foi,
#         space = 0,
#         log = TRUE
#         )
# }
# test

#######################################################################
###
###   User defined distribution for likelihood for
###   infected deer mortalities for radio marked deer that
###   enter the study as test positive at capture
###
###   d_fit_icap_mort
###
###   Overleaf Equation (17)
###
#######################################################################

dIcapMort <- nimble::nimbleFunction(
    run = function(
        ### argument type declarations
        x = integer(0),
        e = integer(0), #e, age of entry
        r = integer(0), #r, age of last known alive
        s = integer(0), #s, age of known mortality
        sex = integer(0),
        fast = integer(0),
        age2date = integer(0),
        beta_sex = double(0),
        beta0_sus = double(0),
        beta0_inf = double(0),
        age_effect_surv = double(1),
        period_effect_surv = double(1),
        f_age_foi = double(1),
        m_age_foi = double(1),
        age_lookup_f = double(1),
        age_lookup_m = double(1),
        period_lookup = double(1),
        f_period_foi = double(1),
        m_period_foi = double(1),
        space = double(0),
        log = double()
        ) {

    llik <- 0 #intialize log-likelihood
    lam_inf <- nimNumeric(s - 1)
    lam_foi <- nimNumeric(s - 1)
    lam_sus <- nimNumeric(s - 1)
    lik_temp <- nimNumeric(e - 1)

    #############################################
    # preliminary hazards for the likelihood
    #############################################

    n_indx_sus <- length(1:(e - 2))
    #survival hazard for susceptible deer
    lam_sus[1:(e - 2)] <- exp(rep(beta0_sus, n_indx_sus) +
                    age_effect_surv[1:(e - 2)] +
                    period_effect_surv[(1 + age2date):(e - 2 + age2date)] +
                    rep(beta_sex * sex, n_indx_sus))

    #survival hazard while infected
    n_indx_inf <- length(1:(s - 1))
    lam_inf[1:(s - 1)] <- exp(rep(beta0_inf, n_indx_inf) +
        age_effect_surv[1:(s - 1)] +
        period_effect_surv[(1 + age2date):(s - 1 + age2date)] +
        rep(beta_sex * sex, n_indx_inf)
        )

    #force of infection infection hazard
    lam_foi[1:(e - 1)] <- exp(rep(space, e - 1) +
        sex * (f_age_foi[age_lookup_f[1:(e - 1)]] +
            f_period_foi[period_lookup[(1 + age2date):(e - 1 + age2date)]]) +
        (1 - sex) * (m_age_foi[age_lookup_m[1:(e - 1)]] +
            m_period_foi[period_lookup[(1 + age2date):(e - 1 + age2date)]])
        )

    #######################################
    ### calculating the joint likelihood
    #######################################
    lik_temp[1] <- lam_foi[1] * exp(-sum(lam_inf[1:(e - 1)]))

    for(k in 2:(e - 1)){
     lik_temp[k] <- lam_foi[k] *
               exp(-sum(lam_sus[1:(k - 1)])) *
               exp(-sum(lam_foi[1:(k - 1)])) *
               exp(-sum(lam_inf[k:(e - 1)]))
    }
    lik <- (1 - fast) * (exp(-sum(lam_inf[e:(r - 1)])) *
           (1 - exp(-sum(lam_inf[r:(s - 1)]))) *
           sum(lik_temp[1:(e - 1)])) +
           fast * ((1 - exp(-sum(lam_inf[r:(s - 1)]))) *
           sum(lik_temp[1:(e - 1)]))
    llik <- log(lik)
    returnType(double(0))
    if(log) return(llik) else return(exp(llik))    ## return log-likelihood
  })


nimble::registerDistributions(list(
    dIcapMort = list(
        BUGSdist = 'dIcapMort(e,r,s,sex,fast,age2date,beta_sex,beta0_sus,beta0_inf,age_effect_surv,period_effect_surv,f_age_foi,m_age_foi,age_lookup_f,age_lookup_m,f_period_foi,m_period_foi,period_lookup,space)',
        types = c("e = integer(0)",
                    "r = integer(0)",
                    "s = integer(0)",
                    "sex = integer(0)",
                    "fast = integer(0)",
                    "age2date = integer(0)",
                    "beta_sex = double(0)",
                    "beta0_sus = double(0)",
                    "beta0_inf = double(0)",
                    "age_effect_surv = double(1)",
                    "period_effect_surv = double(1)",
                    "f_age_foi = double(1)",
                    "m_age_foi = double(1)",
                    "age_lookup_f = double(1)",
                    "age_lookup_m = double(1)",
                    "period_lookup=double(1)",
                    "f_period_foi=double(1)",
                    "m_period_foi=double(1)",
                    "space = double(0)",
                    "log = double()"
                  ),
        discrete = TRUE
    )
))

### for a user-defined distribution
assign('dIcapMort', dIcapMort, envir = .GlobalEnv)

# i=2
# dIcapMort(
#         x = 1,
#         e = d_fit_icap_mort$left_age_e[i],
#         r = d_fit_icap_mort$right_age_r[i],
#         s = d_fit_icap_mort$right_age_s[i],
#         sex = d_fit_icap_mort$sex[i],
#         age2date = icap_mort_age2date[i],
#         beta_sex = beta_sex,
#         beta0_sus = beta0_sus,
#         beta0_inf = beta0_inf,
#         age_effect_surv = age_effect_survival_test,
#         period_effect_surv = period_effect_survival_test,
#         f_age_foi = f_age_foi,
#         m_age_foi = m_age_foi,
#         age_lookup_f = age_lookup_col_f,
#         age_lookup_m = age_lookup_col_m,
#         period_lookup = period_lookup,
#         f_period_foi = f_period_foi,
#         m_period_foi = m_period_foi,
#         space = 0,
#         log = TRUE
#         )

# test <- c()
# for(i in 1:nrow(d_fit_icap_mort)){
# test[i] <- dIcapMort(
#         x = 1,
#         e = d_fit_icap_mort$left_age_e[i],
#         r = d_fit_icap_mort$right_age_r[i],
#         s = d_fit_icap_mort$right_age_s[i],
#         sex = d_fit_icap_mort$sex[i],
#         age2date = icap_mort_age2date[i],
#         beta_sex = beta_sex,
#         beta0_sus = beta0_sus,
#         beta0_inf = beta0_inf,
#         age_effect_surv = age_effect_survival_test,
#         period_effect_surv = period_effect_survival_test,
#         f_age_foi = f_age_foi,
#         m_age_foi = m_age_foi,
#         age_lookup_f = age_lookup_col_f,
#         age_lookup_m = age_lookup_col_m,
#         period_lookup = period_lookup,
#         f_period_foi = f_period_foi,
#         m_period_foi = m_period_foi,
#         space = 0,
#         log = TRUE
#         )
#  }
# test

#######################################################################
###
###   User defined distribution for likelihood for
###   uninfected deer that were test neg at capture,
###   then test negative at recap, that are right censored, 
###   and have been tested post censoring
###
###   d_fit_rec_neg_cens_posttest
###
###   Overleaf Equation (19)
###
#######################################################################

dRecNegCensTest <- nimble::nimbleFunction(
    run = function(
        ### argument type declarations
        x = integer(0),
        e = integer(0), #e, age of entry
        r = integer(0), #r, age of last known alive
        sex = integer(0),
        age2date = integer(0),
        beta_sex = double(0),
        beta0_sus = double(0),
        beta0_inf = double(0),
        age_effect_surv = double(1),
        period_effect_surv = double(1),
        f_age_foi = double(1),
        m_age_foi = double(1),
        age_lookup_f = double(1),
        age_lookup_m = double(1),
        period_lookup=double(1),
        f_period_foi=double(1),
        m_period_foi=double(1),
        space = double(0),
        log = double()
        ) {

    llik <- 0 #intialize log-likelihood
    lam_foi <- nimNumeric(r - 1)
    lam_sus <- nimNumeric(r - 1)
    lam_inf <- nimNumeric(r - 1)

    #############################################
    # preliminary hazards for the likelihood
    #############################################

    indx_period_surv <- (e + age2date):(r - 1 + age2date)

    #survival hazard for susceptible deer
    lam_sus[e:(r - 1)] <- exp(rep(beta0_sus, r - e) +
                          age_effect_surv[e:(r - 1)] +
                          period_effect_surv[indx_period_surv[e:(r - 1)]] +
                          rep(beta_sex * sex, r - e)
                          )

    #force of infection infection hazard
    indx_foi_age_f <- age_lookup_f[1:(r - 1)]
    indx_foi_age_m <- age_lookup_m[1:(r - 1)]
    indx_foi_period <- period_lookup[(1 + age2date):(r - 1 + age2date)]

    lam_foi[1:(r - 1)] <- exp(rep(space, r - 1) +
            sex * (f_age_foi[indx_foi_age_f[1:(r - 1)]] +
                   f_period_foi[indx_foi_period[1:(r - 1)]]) +
            (1 - sex) * (m_age_foi[indx_foi_age_m[1:(r - 1)]] +
                         m_period_foi[indx_foi_period[1:(r - 1)]])
    )

    #######################################
    ### calculating the joint likelihood
    #######################################

    lik <- exp(-sum(lam_sus[e:(r - 1)])) *
           exp(-sum(lam_foi[1:(r - 1)]))
    llik <- log(lik)
    returnType(double(0))
    if(log) return(llik) else return(exp(llik))    ## return log-likelihood
  })

nimble::registerDistributions(list(
    dRecNegCensTest = list(
        BUGSdist = 'dRecNegCensTest(e,r,sex,age2date,beta_sex,beta0_sus,beta0_inf,age_effect_surv,period_effect_surv,f_age_foi,m_age_foi,age_lookup_f,age_lookup_m,f_period_foi,m_period_foi,period_lookup,space)',
        types = c("e = integer(0)",
                    "r = integer(0)",
                    "sex = integer(0)",
                    "age2date = integer(0)",
                    "beta_sex = double(0)",
                    "beta0_sus = double(0)",
                    "beta0_inf = double(0)",
                    "age_effect_surv = double(1)",
                    "period_effect_surv = double(1)",
                    "f_age_foi = double(1)",
                    "m_age_foi = double(1)",
                    "age_lookup_f = double(1)",
                    "age_lookup_m = double(1)",
                    "period_lookup=double(1)",
                    "f_period_foi=double(1)",
                    "m_period_foi=double(1)",
                    "space = double(0)",
                    "log = double(0)"
                  ),
        discrete = TRUE
    )
))

##for a user-defined distribution
assign('dRecNegCensTest', dRecNegCensTest, envir = .GlobalEnv)

# i=1
# dRecNegCensTest(
#         x = 1,
#         e = d_fit_rec_neg_cens_posttest$left_age_e[i],
#         r = d_fit_rec_neg_cens_posttest$right_age_r[i],
#         dn1 = d_fit_rec_neg_cens_posttest$ageweek_recap[i],
#         sex = d_fit_rec_neg_cens_posttest$sex[i],
#         age2date = rec_neg_cens_posttest_age2date[i],
#         beta_sex = beta_sex,
#         beta0_sus = beta0_sus,
#         beta0_inf = beta0_inf,
#         age_effect_surv = age_effect_survival_test,
#         period_effect_surv = period_effect_survival_test,
#         f_age_foi = f_age_foi,
#         m_age_foi = m_age_foi,
#         age_lookup_f = age_lookup_col_f,
#         age_lookup_m = age_lookup_col_m,
#         period_lookup = period_lookup,
#         f_period_foi = f_period_foi,
#         m_period_foi = m_period_foi,
#         space = 0,
#         log = TRUE
#         )

# test <- c()
# for(i in 1:nrow(d_fit_rec_neg_cens_posttest)){
# test[i] <- dRecNegCensTest(
#         x = 1,
#         e = d_fit_rec_neg_cens_posttest$left_age_e[i],
#         r = d_fit_rec_neg_cens_posttest$right_age_r[i],
#         dn1 = d_fit_rec_neg_cens_posttest$ageweek_recap[i],
#         sex = d_fit_rec_neg_cens_posttest$sex[i],
#         age2date = rec_neg_cens_posttest_age2date[i],
#         beta_sex = beta_sex,
#         beta0_sus = beta0_sus,
#         beta0_inf = beta0_inf,
#         age_effect_surv = age_effect_survival_test,
#         period_effect_surv = period_effect_survival_test,
#         f_age_foi = f_age_foi,
#         m_age_foi = m_age_foi,
#         age_lookup_f = age_lookup_col_f,
#         age_lookup_m = age_lookup_col_m,
#         period_lookup = period_lookup,
#         f_period_foi = f_period_foi,
#         m_period_foi = m_period_foi,
#         space = 0,
#         log = TRUE
#         )  
#  }
# test

#######################################################################
###
###   User defined distribution for likelihood for
###   uninfected deer that were test neg at capture,
###   then test negative at recap, that are right censored, 
###   and have been tested post censoring
###
###   d_fit_rec_neg_cens_postno
###
###   Overleaf Equation (21)
###
#######################################################################

dRecNegCensPostNo <- nimble::nimbleFunction(
    run = function(
        ### argument type declarations
        x = integer(0),
        e = integer(0), #e, age of entry
        r = integer(0), #r, age of last known alive
        dn1 = integer(0), #interval of last test negative
        sex = integer(0),
        age2date = integer(0),
        beta_sex = double(0),
        beta0_sus = double(0),
        beta0_inf = double(0),
        age_effect_surv = double(1),
        period_effect_surv = double(1),
        f_age_foi = double(1),
        m_age_foi = double(1),
        age_lookup_f = double(1),
        age_lookup_m = double(1),
        period_lookup=double(1),
        f_period_foi=double(1),
        m_period_foi=double(1),
        space = double(0),
        log = double()
        ) {

    llik <- 0 #intialize log-likelihood
    lam_foi <- nimNumeric(r)
    lam_sus <- nimNumeric(r)
    lam_inf <- nimNumeric(r)
    lik_temp <- nimNumeric(r-1)

    #############################################
    indx_sus_age <- e:(r - 1)
    n_indx_sus <- r - e
    indx_sus_period <- (e + age2date):(r - 1 + age2date)

    #survival hazard for susceptible deer
    lam_sus[e:(r - 1)] <- exp(
        rep(beta0_sus, n_indx_sus) +
            age_effect_surv[e:(r - 1)] +
            period_effect_surv[indx_sus_period[1:n_indx_sus]] +
            rep(beta_sex * sex, n_indx_sus)
    )

    #survival hazard while infected
    lam_inf[1:(r - 1)] <- exp(rep(beta0_inf, (r - 1)) +
        age_effect_surv[1:(r - 1)] +
        period_effect_surv[(1 + age2date):(r - 1 + age2date)] +
        rep(beta_sex * sex, (r - 1))
        )

    #force of infection infection hazard
    indx_foi_age_f <- age_lookup_f[1:(r - 1)]
    indx_foi_age_m <- age_lookup_m[1:(r - 1)]
    indx_foi_period <- period_lookup[(1 + age2date):(r - 1 + age2date)]

    lam_foi[1:(r - 1)] <- exp(rep(space, r - 1) +
            sex * (f_age_foi[indx_foi_age_f[1:(r - 1)]] +
                   f_period_foi[indx_foi_period[1:(r - 1)]]) +
            (1 - sex) * (m_age_foi[indx_foi_age_m[1:(r - 1)]] +
                         m_period_foi[indx_foi_period[1:(r - 1)]])
    )

    #######################################
    ### calculating the joint likelihood
    #######################################

    lik_temp[(dn1+1)] <- lam_foi[(dn1 + 1)] * exp(-sum(lam_inf[(dn1 + 1):(r - 1)]))

    for(k in (dn1 + 2):(r - 1)){
        lik_temp[k] <- lam_foi[k] *
                      exp(-sum(lam_sus[dn1:(k - 1)])) *
                      exp(-sum(lam_inf[k:(r - 1)]))
    }

    lik <- exp(-sum(lam_sus[e:dn1])) *
           exp(-sum(lam_foi[1:dn1])) *
           sum(lik_temp[(e + 1):(r - 1)]) +
           exp(-sum(lam_foi[1:(r - 1)])) *
           exp(-sum(lam_sus[e:(r - 1)]))

    llik <- log(lik)
    returnType(double(0))
    if(log) return(llik) else return(exp(llik))    ## return log-likelihood
  })

nimble::registerDistributions(list(
    dRecNegCensPostNo = list(
        BUGSdist = 'dRecNegCensPostNo(e,r,dn1,sex,age2date,beta_sex,beta0_inf,beta0_sus,age_effect_surv,period_effect_surv,f_age_foi,m_age_foi,age_lookup_f,age_lookup_m,f_period_foi,m_period_foi,period_lookup,space)',
        types = c("e = integer(0)",
                    "r = integer(0)",
                    "dn1 = integer(0)",
                    "sex = integer(0)",
                    "age2date = integer(0)",
                    "beta_sex = double(0)",
                    "beta0_inf = double(0)",
                    "beta0_sus = double(0)",
                    "age_effect_surv = double(1)",
                    "period_effect_surv = double(1)",
                    "f_age_foi = double(1)",
                    "m_age_foi = double(1)",
                    "age_lookup_f = double(1)",
                    "age_lookup_m = double(1)",
                    "period_lookup=double(1)",
                    "f_period_foi=double(1)",
                    "m_period_foi=double(1)",
                    "space = double(0)",
                    "log = double(0)"
                  ),
        discrete = TRUE
    )
))

###for a user-defined distribution
assign('dRecNegCensPostNo', dRecNegCensPostNo, envir = .GlobalEnv)

# i=1
# dRecNegCensPostNo(
#         x = 1,
#         e = d_fit_rec_neg_cens_postno$left_age_e[i],
#         r = d_fit_rec_neg_cens_postno$right_age_r[i],
#         dn1 = d_fit_rec_neg_cens_postno$ageweek_recap[i],
#         sex = d_fit_rec_neg_cens_postno$sex[i],
#         age2date = rec_neg_cens_postno_age2date[i],
#         beta_sex = beta_sex,
#         beta0_sus = beta0_sus,
#         beta0_inf = beta0_inf,
#         age_effect_surv = age_effect_survival_test,
#         period_effect_surv = period_effect_survival_test,
#         f_age_foi = f_age_foi,
#         m_age_foi = m_age_foi,
#         age_lookup_f = age_lookup_col_f,
#         age_lookup_m = age_lookup_col_m,
#         period_lookup = period_lookup,
#         f_period_foi = f_period_foi,
#         m_period_foi = m_period_foi,
#         space = 0,
#         log = TRUE
#         )


# test <- c()
# for(i in 1:nrow(d_fit_rec_neg_cens_postno)){
# test[i] <- dRecNegCensPostNo(
#         x = 1,
#         e = d_fit_rec_neg_cens_postno$left_age_e[i],
#         r = d_fit_rec_neg_cens_postno$right_age_r[i],
#         dn1 = d_fit_rec_neg_cens_postno$ageweek_recap[i],
#         sex = d_fit_rec_neg_cens_postno$sex[i],
#         age2date = rec_neg_cens_postno_age2date[i],
#         beta_sex = beta_sex,
#         beta0_sus = beta0_sus,
#         beta0_inf = beta0_inf,
#         age_effect_surv = age_effect_survival_test,
#         period_effect_surv = period_effect_survival_test,
#         f_age_foi = f_age_foi,
#         m_age_foi = m_age_foi,
#         age_lookup_f = age_lookup_col_f,
#         age_lookup_m = age_lookup_col_m,
#         period_lookup = period_lookup,
#         f_period_foi = f_period_foi,
#         m_period_foi = m_period_foi,
#         space = 0,
#         log = TRUE
#         )
#  }
# test

#######################################################################
###
###   User defined distribution for likelihood for
###   uninfected deer that were test neg at capture,
###   then test negative at recap,
###   that die
###
###   d_fit_rec_neg_mort
###
###   Overleaf Equation (23)
###
#######################################################################

dRecNegMort <- nimble::nimbleFunction(
    run = function(
        ### argument type declarations
        x = integer(0),
        e = integer(0), #e, age of entry
        r = integer(0), #r, age of last known alive
        s = integer(0), #s, age of mortality
        sex = integer(0),
        age2date = integer(0),
        beta_sex = double(0),
        beta0_sus = double(0),
        age_effect_surv = double(1),
        period_effect_surv = double(1),
        f_age_foi = double(1),
        m_age_foi = double(1),
        age_lookup_f = double(1),
        age_lookup_m = double(1),
        period_lookup=double(1),
        f_period_foi=double(1),
        m_period_foi=double(1),
        space = double(0),
        log = double(0)
        ) {

    llik <- 0 #intialize log-likelihood
    lam_foi <- nimNumeric(s-1)
    lam_sus <- nimNumeric(s-1)

    #############################################
    # preliminary hazards for the likelihood
    #############################################

    #survival hazard for susceptible deer
    lam_sus[e:(s-1)] <- exp(rep(beta0_sus, s - e) +
                          age_effect_surv[e:(s - 1)] +
                          period_effect_surv[(e + age2date):(s - 1 + age2date)] +
                          rep(beta_sex * sex, s - e)
                      )
    #force of infection infection hazard
    indx_foi_age_f <- age_lookup_f[1:(s - 1)]
    indx_foi_age_m <- age_lookup_m[1:(s - 1)]
    indx_foi_period <- period_lookup[(1 + age2date):(s - 1 + age2date)]

    lam_foi[1:(s - 1)] <- exp(rep(space, s - 1) +
            sex * (f_age_foi[indx_foi_age_f[1:(s - 1)]] +
                   f_period_foi[indx_foi_period[1:(s - 1)]]) +
            (1 - sex) * (m_age_foi[indx_foi_age_m[1:(s - 1)]] +
                         m_period_foi[indx_foi_period[1:(s - 1)]])
    )

    #######################################
    ### calculating the joint likelihood
    #######################################

    lik <- exp(-sum(lam_sus[e:(r - 1)])) *
           (1 - exp(-sum(lam_sus[(r:(s - 1))]))) *
           exp(-sum(lam_foi[1:(s - 1)]))

    llik <- log(lik)

    returnType(double(0))
    if(log) return(llik) else return(exp(llik))    ## return log-likelihood
  })

nimble::registerDistributions(list(
    dRecNegMort = list(
        BUGSdist = 'dRecNegMort(e,r,s,sex,age2date,beta_sex,beta0_sus,age_effect_surv,period_effect_surv,f_age_foi,m_age_foi,age_lookup_f,age_lookup_m,f_period_foi,m_period_foi,period_lookup,space)',
        types = c("e = integer(0)",
                  "r = integer(0)",
                  "s = integer(0)",
                  "sex = integer(0)",
                  "age2date = integer(0)",
                  "beta_sex = double(0)",
                  "beta0_sus = double(0)",
                  "age_effect_surv = double(1)",
                  "period_effect_surv = double(1)",
                  "f_age_foi = double(1)",
                  "m_age_foi = double(1)",
                  "age_lookup_f = double(1)",
                  "age_lookup_m = double(1)",
                  "period_lookup = double(1)",
                  "f_period_foi = double(1)",
                  "m_period_foi = double(1)",
                  "space = double(0)",
                  "log = double(0)"
                  ),
        discrete = TRUE
    )
))

# ###for a user-defined distribution
assign('dRecNegMort', dRecNegMort, envir = .GlobalEnv)

# i=1
# dRecNegMort(
#         x = 1,
#         e = d_fit_rec_neg_mort$left_age_e[i],
#         r = d_fit_rec_neg_mort$right_age_r[i],
#         s = d_fit_rec_neg_mort$right_age_s[i],
#         sex = d_fit_rec_neg_mort$sex[i],
#         age2date = rec_neg_mort_age2date[i],
#         beta_sex = beta_sex,
#         beta0_sus = beta0_sus,
#         age_effect_surv = age_effect_survival_test,
#         period_effect_surv = period_effect_survival_test,
#         f_age_foi = f_age_foi,
#         m_age_foi = m_age_foi,
#         age_lookup_f = age_lookup_col_f,
#         age_lookup_m = age_lookup_col_m,
#         period_lookup = period_lookup,
#         f_period_foi = f_period_foi,
#         m_period_foi = m_period_foi,
#         space = 0,
#         log = TRUE
#         )

# test <- c()
# for(i in 1:nrow(d_fit_rec_neg_mort)){
# test[i] <- dRecNegMort(
#         x = 1,
#         e = d_fit_rec_neg_mort$left_age_e[i],
#         r = d_fit_rec_neg_mort$right_age_r[i],
#         s = d_fit_rec_neg_mort$right_age_s[i],
#         sex = d_fit_rec_neg_mort$sex[i],
#         age2date = rec_neg_mort_age2date[i],
#         beta_sex = beta_sex,
#         beta0_sus = beta0_sus,
#         age_effect_surv = age_effect_survival_test,
#         period_effect_surv = period_effect_survival_test,
#         f_age_foi = f_age_foi,
#         m_age_foi = m_age_foi,
#         age_lookup_f = age_lookup_col_f,
#         age_lookup_m = age_lookup_col_m,
#         period_lookup = period_lookup,
#         f_period_foi = f_period_foi,
#         m_period_foi = m_period_foi,
#         space = 0,
#         log = TRUE
#         )
#  }
# test

#######################################################################
###
###   User defined distribution for likelihood for
###   deer that were test neg at capture,
###   then test positive at recap,
###   than die
###
###   d_fit_rec_pos_mort
###
###   Overleaf Equation (25)
###
#######################################################################

dRecPosMort <- nimble::nimbleFunction(
    run = function(
        ### argument type declarations
        x = integer(0),
        e = integer(0), #e, age of entry
        r = integer(0), #r, age of last known alive
        s = integer(0), #s, age of last known alive
        dn1 = integer(0), #interval of last test negative
        dn = integer(0), #int tested positive
        sex = integer(0),
        age2date = integer(0),
        beta_sex = double(0),
        beta0_sus = double(0),
        beta0_inf = double(0),
        age_effect_surv = double(1),
        period_effect_surv = double(1),
        f_age_foi = double(1),
        m_age_foi = double(1),
        age_lookup_f = double(1),
        age_lookup_m = double(1),
        period_lookup=double(1),
        f_period_foi=double(1),
        m_period_foi=double(1),
        space = double(0),
        log = double()
        ) {

    lik <- 0 #intialize log-likelihood
    llik <- 0 #intialize log-likelihood
    lam_foi <- nimNumeric(s - 1)
    lam_sus <- nimNumeric(s - 1)
    lam_inf <- nimNumeric(s - 1)
    lik_temp  <- nimNumeric(s - 1)

    #############################################
    # preliminary hazards for the likelihood
    #############################################

    #survival hazard for susceptible deer
    n_indx_sus <- length(e:(dn - 1))
    lam_sus[e:(dn - 1)] <- exp(rep(beta0_sus, n_indx_sus)  +
                    age_effect_surv[e:(dn - 1)] +
                    period_effect_surv[(e + age2date):(dn - 1 + age2date)] +
                    rep(beta_sex * sex, n_indx_sus))

    #survival hazard while infected
    n_indx_inf <- length(dn1:(s - 1))
    lam_inf[dn1:(s - 1)] <- exp(rep(beta0_inf,n_indx_inf) +
                                age_effect_surv[dn1:(s - 1)] +
                                period_effect_surv[(dn1 + age2date):
                                                   (s - 1 + age2date)] +
                                rep(beta_sex * sex,n_indx_inf)
                            )
    #force of infection infection hazard
    lam_foi[1:(dn - 1)] <- exp(rep(space,(dn - 1)) +
                  sex * (f_age_foi[age_lookup_f[1:(dn - 1)]] +
                            f_period_foi[period_lookup[(1 + age2date):((dn - 1) + age2date)]])
                             +
                  (1 - sex) * (m_age_foi[age_lookup_m[1:(dn - 1)]] +
                                  m_period_foi[period_lookup[(1 + age2date):((dn - 1) + age2date)]])
                  )
    #######################################
    ### calculating the joint likelihood
    #######################################
      lik_temp[(dn1 + 1)] <- lam_foi[(dn1 + 1)] *
               exp(-sum(lam_foi[1:dn1])) *
               exp(-sum(lam_inf[(dn1 + 1):(r - 1)]))

    for(k in (dn1+2):(dn-1)) {
      lik_temp[k] <- lam_foi[k] *
               exp(-sum(lam_sus[e:(k - 1)])) *
               exp(-sum(lam_foi[1:(k - 1)])) *
               exp(-sum(lam_inf[k:(r - 1)]))
    }

    lik <- (1 - exp(-sum(lam_inf[r:(s - 1)]))) *
            sum(lik_temp[dn1:(dn - 1)])
    llik <- log(lik)

    returnType(double(0))
    if(log) return(llik) else return(exp(llik))    ## return log-likelihood
  })

nimble::registerDistributions(list(
    dRecPosMort = list(
        BUGSdist = 'dRecPosMort(e,r,s,dn1,dn,sex,age2date,beta_sex,beta0_sus,beta0_inf,age_effect_surv,period_effect_surv,f_age_foi,m_age_foi,age_lookup_f,age_lookup_m,f_period_foi,m_period_foi,period_lookup,space)',
        types = c("e = integer(0)",
                    "r = integer(0)",
                    "s = integer(0)",
                    "dn1 = integer(0)",
                    "dn = integer(0)",
                    "sex = integer(0)",
                    "age2date = integer(0)",
                    "beta_sex = double(0)",
                    "beta0_sus = double(0)",
                    "beta0_inf = double(0)",
                    "age_effect_surv = double(1)",
                    "period_effect_surv = double(1)",
                    "f_age_foi = double(1)",
                    "m_age_foi = double(1)",
                    "age_lookup_f = double(1)",
                    "age_lookup_m = double(1)",
                    "period_lookup=double(1)",
                    "f_period_foi=double(1)",
                    "m_period_foi=double(1)",
                    "space = double(0)",
                    "log = double(0)"
                  ),
        discrete = TRUE
    )
))

###Global Declaration so Nimble can access
assign('dRecPosMort', dRecPosMort, envir = .GlobalEnv)

# i=2
# dRecPosMort(
#         x = 1,
#         e = d_fit_rec_pos_mort$left_age_e[i],
#         r = d_fit_rec_pos_mort$right_age_r[i],
#         s = d_fit_rec_pos_mort$right_age_s[i],
#         dn1 = d_fit_rec_pos_mort$left_age_e[i],
#         dn = d_fit_rec_pos_mort$ageweek_recap[i],
#         sex = d_fit_rec_pos_mort$sex[i],
#         age2date = rec_pos_mort_age2date[i],
#         beta_sex = beta_sex,
#         beta0_sus = beta0_sus,
#         beta0_inf = beta0_inf,
#         age_effect_surv = age_effect_survival_test,
#         period_effect_surv = period_effect_survival_test,
#         f_age_foi = f_age_foi,
#         m_age_foi = m_age_foi,
#         age_lookup_f = age_lookup_col_f,
#         age_lookup_m = age_lookup_col_m,
#         period_lookup = period_lookup,
#         f_period_foi = f_period_foi,
#         m_period_foi = m_period_foi,
#         space = 0,
#         log = TRUE
#         )

# test <- c()
# for(i in 1:nrow(d_fit_rec_pos_mort)){
# test[i] <- dRecPosMort(
#         x = 1,
#         e = d_fit_rec_pos_mort$left_age_e[i],
#         r = d_fit_rec_pos_mort$right_age_r[i],
#         s = d_fit_rec_pos_mort$right_age_s[i],
#         dn1 = d_fit_rec_pos_mort$left_age_e[i],
#         dn = d_fit_rec_pos_mort$ageweek_recap[i],
#         sex = d_fit_rec_pos_mort$sex[i],
#         age2date = rec_pos_mort_age2date[i],
#         beta_sex = beta_sex,
#         beta0_sus = beta0_sus,
#         beta0_inf = beta0_inf,
#         age_effect_surv = age_effect_survival_test,
#         period_effect_surv = period_effect_survival_test,
#         f_age_foi = f_age_foi,
#         m_age_foi = m_age_foi,
#         age_lookup_f = age_lookup_col_f,
#         age_lookup_m = age_lookup_col_m,
#         period_lookup = period_lookup,
#         f_period_foi = f_period_foi,
#         m_period_foi = m_period_foi,
#         space = 0,
#         log = TRUE
#         )
#  }
# test

#######################################################################
###
###   User defined distribution for likelihood for
###   infected deer that were test neg at capture,
###   then test positive at recap,
###   than these were right censored
###
###   d_fit_rec_pos_cens
###
###   Overleaf Equation (27)
###
#######################################################################

dRecPosCens <- nimble::nimbleFunction(
    run = function(
        ### argument type declarations
        x = integer(0),
        e = integer(0), #e, age of entry
        r = integer(0), #r, age of last known alive
        dn1 = integer(0), #interval of last test negative
        dn = integer(0), #int tested positive
        sex = integer(0),
        age2date = integer(0),
        beta_sex = double(0),
        beta0_sus = double(0),
        beta0_inf = double(0),
        age_effect_surv = double(1),
        period_effect_surv = double(1),
        f_age_foi = double(1),
        m_age_foi = double(1),
        age_lookup_f = double(1),
        age_lookup_m = double(1),
        period_lookup=double(1),
        f_period_foi=double(1),
        m_period_foi=double(1),
        space = double(0),
        log = double()
        ) {

    llik<-0 #intialize log-likelihood
    lam_foi <- nimNumeric(r-1)
    lam_sus <- nimNumeric(r-1)
    lam_inf <- nimNumeric(r-1)
    lik_temp <- nimNumeric(r-1)

    #############################################
    # preliminary hazards for the likelihood
    #############################################

    #survival hazard for susceptible deer
    n_indx_sus <- length(e:(dn - 1))
    lam_sus[e:(dn - 1)] <- exp(rep(beta0_sus, n_indx_sus)  +
                    age_effect_surv[e:(dn - 1)] +
                    period_effect_surv[(e + age2date):((dn - 1) + age2date)] +
                    rep(beta_sex * sex, n_indx_sus))
    
    #survival hazard while infected
    n_indx_inf <- length(dn1:(r - 1))
    
    lam_inf[dn1:(r - 1)] <- exp(rep(beta0_inf,n_indx_inf) +
        age_effect_surv[dn1:(r - 1)] +
        period_effect_surv[(dn1 + age2date):(r - 1 + age2date)] +
        rep(beta_sex * sex,n_indx_inf)
        )

    #force of infection infection hazard
    lam_foi[1:(dn - 1)] <- exp(rep(space, (dn - 1)) +
        sex * (f_age_foi[age_lookup_f[1:(dn - 1)]] +
            f_period_foi[period_lookup[(1 + age2date):((dn - 1) + age2date)]]) +
        (1 - sex) * (m_age_foi[age_lookup_m[1:(dn - 1)]] +
            m_period_foi[period_lookup[(1 + age2date):((dn - 1) + age2date)]])
        )

    #######################################
    ### calculating the joint likelihood
    #######################################

    lik_temp[(dn1 + 1)] <- lam_foi[(dn1 + 1)] *
                    exp(-sum(lam_foi[1:dn1])) *
                    exp(-sum(lam_inf[(dn1 + 1):(r - 1)]))

    for(k in (dn1 + 2):(dn - 1)) {
      lik_temp[k] <- lam_foi[k] *
                    exp(-sum(lam_sus[(dn1 + 1):(k - 1)])) *
                    exp(-sum(lam_foi[1:(k-1)])) *
                    exp(-sum(lam_inf[k:(r-1)]))
    }
    lik <- exp(-sum(lam_sus[e:dn1])) *
           exp(-sum(lam_foi[1:dn1])) *
           sum(lik_temp[(dn1 + 1):(dn - 1)])
    llik <- log(lik)
    returnType(double(0))
    if(log) return(llik) else return(exp(llik))    ## return log-likelihood
  })

nimble::registerDistributions(list(
    dRecPosCens = list(
        BUGSdist = 'dRecPosCens(e,r,dn1,dn,sex,age2date,beta_sex,beta0_sus,beta0_inf,age_effect_surv,period_effect_surv,f_age_foi,m_age_foi,age_lookup_f,age_lookup_m,f_period_foi,m_period_foi,period_lookup,space)',
        types = c("e = integer(0)",
                    "r = integer(0)",
                    "dn1 = integer(0)",
                    "dn = integer(0)",
                    "sex = integer(0)",
                    "age2date = integer(0)",
                    "beta_sex = double(0)",
                    "beta0_inf = double(0)",
                    "beta0_sus = double(0)",
                    "age_effect_surv = double(1)",
                    "period_effect_surv = double(1)",
                    "f_age_foi = double(1)",
                    "m_age_foi = double(1)",
                    "age_lookup_f = double(1)",
                    "age_lookup_m = double(1)",
                    "period_lookup=double(1)",
                    "f_period_foi=double(1)",
                    "m_period_foi=double(1)",
                    "space = double(0)",
                    "log = double(0)"
                  ),
        discrete = TRUE
    )
))

###for a user-defined distribution
assign('dRecPosCens', dRecPosCens, envir = .GlobalEnv)

# i=1
# dRecPosCens(
#         x = 1,
#         e = d_fit_rec_pos_cens$left_age_e[i],
#         r = d_fit_rec_pos_cens$right_age_r[i],
#         dn1 = d_fit_rec_pos_cens$left_age_e[i],
#         dn = d_fit_rec_pos_cens$ageweek_recap[i],
#         sex = d_fit_rec_pos_cens$sex[i],
#         age2date = rec_pos_cens_age2date[i],
#         beta_sex = beta_sex,
#         beta0_sus = beta0_sus,
#         beta0_inf = beta0_inf,
#         age_effect_surv = age_effect_survival_test,
#         period_effect_surv = period_effect_survival_test,
#         f_age_foi = f_age_foi,
#         m_age_foi = m_age_foi,
#         age_lookup_f = age_lookup_col_f,
#         age_lookup_m = age_lookup_col_m,
#         period_lookup = period_lookup,
#         f_period_foi = f_period_foi,
#         m_period_foi = m_period_foi,
#         space = 0,
#         log = TRUE
#         )


# test <- c()
# for(i in 1:nrow(d_fit_rec_pos_cens)){
# test[i] <- dRecPosCens(
#         x = 1,
#         e = d_fit_rec_pos_cens$left_age_e[i],
#         r = d_fit_rec_pos_cens$right_age_r[i],
#         dn1 = d_fit_rec_pos_cens$left_age_e[i],
#         dn = d_fit_rec_pos_cens$ageweek_recap[i],
#         sex = d_fit_rec_pos_cens$sex[i],
#         age2date = rec_pos_cens_age2date[i],
#         beta_sex = beta_sex,
#         beta0_sus = beta0_sus,
#         beta0_inf = beta0_inf,
#         age_effect_surv = age_effect_survival_test,
#         period_effect_surv = period_effect_survival_test,
#         f_age_foi = f_age_foi,
#         m_age_foi = m_age_foi,
#         age_lookup_f = age_lookup_col_f,
#         age_lookup_m = age_lookup_col_m,
#         period_lookup = period_lookup,
#         f_period_foi = f_period_foi,
#         m_period_foi = m_period_foi,
#         space = 0,
#         log = TRUE
#         )
#  }
# test


#######################################################################
###
###   User defined distribution for likelihood for
###   infected deer mortalities for radio marked deer that
###   enter the study as test negative at capture
###
###   d_fit_idead
###
###   Overleaf Equation (29)
###
#######################################################################

dNegCapPosMort <- nimble::nimbleFunction(
    run = function(
        ### argument type declarations
        x = integer(0),
        n_samples = integer(0), # number of samples in dataset
        e = double(1), #e, age of entry
        r = double(1), #r, age of last known alive
        s = double(1), #s, age of known mortality
        dn1 = double(1), #last interval of test negative
        dn = double(1), #interval of test positive (could be interval of mortality)
        sex = double(1),
        age2date = double(1),
        beta_sex = double(0),
        beta0_sus = double(0),
        beta0_inf = double(0),
        age_effect_surv = double(1),
        period_effect_surv = double(1),
        f_age_foi = double(1),
        m_age_foi = double(1),
        age_lookup_f = double(1),
        age_lookup_m = double(1),
        period_lookup=double(1),
        f_period_foi=double(1),
        m_period_foi=double(1),
		sect = double(1),
        space = double(1),
        log = double()
        ) {

    sumllik <- 0 #intialize log-likelihood
	for(i in 1:n_samples){
		lam_foi <- nimNumeric(s[i] - 1,init=FALSE)
		lam_sus <- nimNumeric(s[i] - 1,init=FALSE)
		lam_inf <- nimNumeric(s[i] - 1,init=FALSE)
		lik_temp  <- 0
	
		#############################################
		# preliminary hazards for the likelihood
		#############################################
		if(sex[i] ==1){ # females
			#survival hazard for susceptible deer
			lam_sus[e[i]:(dn[i] - 1)] <- exp(beta0_sus  + 
					age_effect_surv[e[i]:(dn[i] - 1)] +
					period_effect_surv[(e[i] + age2date[i]):(dn[i] - 1 + age2date[i])] +
					beta_sex)
			#survival hazard while infected
			lam_inf[dn1[i]:(s[i] - 1)] <- exp(beta0_inf +
					age_effect_surv[dn1[i]:(s[i] - 1)] +
					period_effect_surv[(dn1[i] + age2date[i]):(s[i] - 1 + age2date[i])] +
					beta_sex)
			#force of infection infection hazard
			lam_foi[1:(dn[i] - 1)] <- exp(space[sect[i]] +
					f_age_foi[age_lookup_f[1:(dn[i] - 1)]] +
					f_period_foi[period_lookup[(1 + age2date[i]):(dn[i] - 1 + age2date[i])]])
		}else{ # males
			#survival hazard for susceptible deer
			lam_sus[e[i]:(dn[i] - 1)] <- exp(beta0_sus  + 
					age_effect_surv[e[i]:(dn[i] - 1)] +
					period_effect_surv[(e[i] + age2date[i]):(dn[i] - 1 + age2date[i])])
			#survival hazard while infected
			lam_inf[dn1[i]:(s[i] - 1)] <- exp(beta0_inf +
					age_effect_surv[dn1[i]:(s[i] - 1)] +
					period_effect_surv[(dn1[i] + age2date[i]):(s[i] - 1 + age2date[i])])
			#force of infection infection hazard
			lam_foi[1:(dn[i] - 1)] <- exp(space[sect[i]] +
					m_age_foi[age_lookup_m[1:(dn[i] - 1)]] +
					m_period_foi[period_lookup[(1 + age2date[i]):(dn[i] - 1 + age2date[i])]])
		}

		#######################################
		### calculating the joint likelihood
		#######################################
	
		if((r[i] - dn1[i]) > 1){
			lik_temp[dn1[i] + 1] <- lam_foi[dn1[i] + 1] *
						exp(-sum(lam_inf[(dn1[i] + 1):(r[i] - 1)])) *
						(1 - exp(-sum(lam_inf[r[i]:(s[i] - 1)])))
			lik_temp[(s[i] - 1)] <- lam_foi[(s[i] - 1)] *
							exp(-sum(lam_foi[(dn1[i] + 1):((s[i] - 1) - 1)])) *
							#exp(-sum(lam_foi[(dn1[i] + 1):(s[i] - 1)])) *
							exp(-sum(lam_sus[(dn1[i] + 1):((s[i] - 1) - 1)])) *
							lam_inf[(s[i] - 1)]
		}
	
		if((r[i] - dn1[i]) == 1){
			lik_temp[(s[i] - 1)] <- lam_foi[(s[i] - 1)] *
							lam_inf[(s[i] - 1)]
		}
		
		if((r[i] - dn1[i]) > 3){
		for (k in (dn1[i] + 2):(r[i] - 2)) {
			lik_temp[k] <- lam_foi[k] *
						exp(-sum(lam_foi[(dn1[i] + 1):(k - 1)])) *
						#exp(-sum(lam_foi[(dn1[i] + 1):k])) *
						exp(-sum(lam_sus[(dn1[i] + 1):(k - 1)])) *
						exp(-sum(lam_inf[k:(r[i] - 1)])) *
						(1 - exp(-sum(lam_inf[r[i]:(s[i] - 1)])))
		}
		}
		if((r[i] - dn1[i]) > 2){
		for (k in (r[i] - 1):(s[i] - 2)) {
			lik_temp[k] <- lam_foi[k] *
						exp(-sum(lam_foi[(dn1[i] + 1):(k - 1)])) *
						#exp(-sum(lam_foi[(dn1[i] + 1):k])) *
						exp(-sum(lam_sus[(dn1[i] + 1):(k - 1)])) *
						(1 - exp(-sum(lam_inf[k:(s[i] - 1)])))
		}
		}
							
		sumllik <- sumllik -sum(lam_sus[e[i]:dn1[i]]) -
			sum(lam_foi[1:dn1[i]]) +
			log(sum(lik_temp[(dn1[i] + 1):(s[i] - 1)]))
	
	}
    returnType(double(0))
    if(log) return(sumllik) else return(exp(sumllik))    ## return log-likelihood
  })


nimble::registerDistributions(list(
    dNegCapPosMort = list(
        BUGSdist = 'dNegCapPosMort(n_samples,e,r,s,dn1,dn,sex,age2date,beta_sex,beta0_sus,beta0_inf,age_effect_surv,period_effect_surv,f_age_foi,m_age_foi,age_lookup_f,age_lookup_m,f_period_foi,m_period_foi,period_lookup,sect,space)',
        types = c("value = integer(0)",
				  "n_samples = integer(0)",
				  "e = double(1)",
                  "r = double(1)",
                  "s = double(1)",
                  "dn1 = double(1)",
				  "dn = double(0)",
                  "sex = double(1)",
                  "age2date = double(1)",
                  "beta_sex = double(0)",
                  "beta0_sus = double(0)",
                  "beta0_inf = double(0)",
                  "age_effect_surv = double(1)",
                  "period_effect_surv = double(1)",
                  "f_age_foi = double(1)",
                  "m_age_foi = double(1)",
                  "age_lookup_f = double(1)",
                  "age_lookup_m = double(1)",
                  "period_lookup = double(1)",
                  "f_period_foi = double(1)",
                  "m_period_foi = double(1)",
				  "sect = double(1)",
                  "space = double(0)",
                  "log = double()"
                  ),
        discrete = TRUE
    )
))

### for a user-defined distribution
assign('dNegCapPosMort', dNegCapPosMort, envir = .GlobalEnv)

start <- Sys.time()
test <-  dNegCapPosMort(
        x = 1,
		  n_samples = nrow(d_fit_idead),
        e = d_fit_idead$left_age_e,
        r = d_fit_idead$right_age_r,
        s = d_fit_idead$right_age_s,
        dn1 = d_fit_idead$left_age_e,
        dn = d_fit_idead$right_age_s,
        sex = d_fit_idead$sex,
        age2date = idead_age2date,
        beta_sex = beta_sex,
        beta0_sus = beta0_sus,
        beta0_inf = beta0_inf,
        age_effect_surv = age_effect_survival_test,
        period_effect_surv = period_effect_survival_test,
        f_age_foi = f_age_foi,
        m_age_foi = m_age_foi,
        age_lookup_f = age_lookup_col_f,
        age_lookup_m = age_lookup_col_m,
        period_lookup = period_lookup,
        f_period_foi = f_period_foi,
        m_period_foi = m_period_foi,
		  sect = d_fit_idead$sect,
        space = rep(0,n_sect),
        log = TRUE
        )
(end<- Sys.time()-start)
test
