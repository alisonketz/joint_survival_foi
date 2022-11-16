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
###
#######################################################################

dInfHarvest <- nimble::nimbleFunction(
    run = function(
        ### argument type declarations
        x=double(0),
        a = double(0), #age (weeks) at harvest
        age2date_surv = double(0),
        beta0_sus = double(0),
        beta0_inf = double(0),
        period_effect_surv = double(0),
        age_effect_surv = double(0),
        sex = double(0),
        beta_sex = double(0),
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
    lam_foi <- nimNumeric(a)
    lam_sus <- nimNumeric(a - 1)
    lam_inf <- nimNumeric(a)
    lik_temp <- nimNumeric(a)

    #############################################
    # preliminary hazards for the likelihood
    #############################################
    indx_foi_age_f <- age_lookup_f[1:a]
    indx_foi_age_m <- age_lookup_m[1:a]
    indx_foi_period <- period_lookup[(1 + age2date_foi):(a + age2date_foi)]

    lam_foi[1:a] <- exp(rep(space, a) +
                        sex * (f_age_foi[indx_foi_age_f[1:a]] + 
                               f_period_foi[indx_foi_period[1:a]])
                        (1 - sex) * (m_age_foi[indx_foi_age_m[1:a]] + 
                                     m_period_foi[indx_foi_period[1:a]])
            )

    indx_sus_period <- (1 + age2date_surv):(a - 1 + age2date_surv)
    lam_sus[1:(a - 1)] <- exp(rep(beta0_sus, (a - 1)) +
                            age_effect_surv[1:(a - 1)] +
                            period_effect_surv[indx_sus_period[1:( a - 1)]] +
                            rep(beta_sex * sex, (a - 1))
                            )

    indx_inf_period <- (1 + age2date_surv):(a + age2date_surv)
    lam_inf[1:a] <- exp(rep(beta0_inf, a) +
                        age_effect_surv[1:a] +
                        period_effect_surv[indx_inf_period[1:a]] +
                        rep(beta_sex * sex, a)
                        )
    #######################################
    ### calculating the joint likelihood
    #######################################

    lik_temp[1] <- lam_foi[1] * exp(-sum(lam_inf[1:(a - 1)]))

    for(j in 2:(a-1)){
        lik_temp[j] <- lam_foi[j] *
               exp(-sum(lam_sus[1:(j - 1)] + lam_foi[1:(j - 1)])) *
               exp(-sum(lam_inf[j:(a - 1)]))
    }
    lik_temp[a] <- lam_foi[a] * exp(-sum(lam_foi[1:(a - 1)] + lam_sus[1:(a - 1)]))

    lik <- lam_inf[a] * sum(lik_temp[1:a])
    llik <- log(lik)

    returnType(double(0))
    if(log) return(llik) else return(exp(llik))    ## return log-likelihood
  })


nimble::registerDistributions(list(
    dInfHarvest = list(
        BUGSdist = 'dInfHarvest(a,age2date_surv,beta0_inf,beta0_sus,period_effect_surv,age_effect_surv,sex,beta_sex,f_age_foi,m_age_foi,age_lookup_f,age_lookup_m,f_period_foi,m_period_foi,period_lookup,space)',
        types = c("a = double(0)",
                    "age2date_surv = double(0)",
                    "beta0_inf = double(0)",
                    "beta0_sus = double(0)",
                    "period_effect_surv = double(1)",
                    "age_effect_surv = double(1)",
                    "sex = double(0)",
                    "beta_sex = double(0)",
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

# for a user-defined distribution
assign('dInfHarvest', dInfHarvest, envir = .GlobalEnv)


#######################################################################
###
###   User defined distribution for likelihood for 
###   uninfected harvest deer
###
###   d_fit_hunt_neg
###
#######################################################################

dSusHarvest <- nimble::nimbleFunction(
    run = function(
        ### argument type declarations
        x = double(0),
        a = double(0),
        age2date_surv = double(0),
        beta0_sus = double(0),
        period_effect_surv = double(0),
        age_effect_surv = double(0),
        sex = double(0),
        beta_sex = double(0),
        age2date_foi = double(0),
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

    llik<-0 #intialize log-likelihood

    lam_foi <- nimNumeric(a)
    lam_sus <- nimNumeric(a)
    lam_temp <- nimNumeric(a)

    #############################################
    # preliminary hazards for the likelihood
    #############################################
    indx_foi_age_f <- age_lookup_f[1:a]
    indx_foi_age_m <- age_lookup_m[1:a]
    indx_foi_period <- period_lookup[(1 + age2date_foi):(a + age2date_foi)]

    lam_foi[1:a] <- exp(rep(space, a) +
                        sex * (f_age_foi[indx_foi_age_f[1:a]] +
                               f_period_foi[indx_foi_period[1:a]])
                        (1 - sex) * (m_age_foi[indx_foi_age_m[1:a]] +
                                     m_period_foi[indx_foi_period[1:a]])
            )

    indx_sus_period <- (1 + age2date_surv):(a - 1 + age2date_surv)
    lam_sus[1:(a - 1)] <- exp(rep(beta0_sus, (a - 1)) +
                            age_effect_surv[1:(a - 1)] +
                            period_effect_surv[indx_sus_period[1:(a - 1)]] +
                            rep(beta_sex * sex, (a - 1))
                            )

    lam_temp[1:(a - 1)] <- lam_foi[1:(a - 1)] + lam_sus[1:(a - 1)]


    #######################################
    ###
    ### calculating the joint likelihood
    ###
    #######################################

    lik <- exp(-(sum(lam_temp[1:(a - 1)]) - lam_foi[a])) * lam_sus[a]
    llik <- log(lik)

    returnType(double(0))
    if(log) return(llik) else return(exp(llik))    ## return log-likelihood
  })


nimble::registerDistributions(list(
    dSusHarvest = list(
        BUGSdist = 'dSusHarvest(a,age2date_surv,beta0_sus,period_effect_surv,age_effect_surv,sex,beta_sex,f_age_foi,m_age_foi,age_lookup_f,age_lookup_m,f_period_foi,m_period_foi,period_lookup,space)',
        types = c("a = double(0)",
                  "age2date_surv = double(0)",
                  "beta0_sus = double(0)",
                  "period_effect_surv = double(1)",
                  "age_effect_surv = double(1)",
                  "sex = double(0)",
                  "beta_sex = double(0)",
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

# for a user-defined distribution
assign('dSusHarvest', dSusHarvest, envir = .GlobalEnv)


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
###   uninfected radio-marked deer interval censored:
###   test neg at cap and tested at interval censoring
###
###   d_fit_sus_cens_posttest
###
###   test neg at cap and Not tested at interval censoring
###
###   d_fit_sus_cens_postno
###
#######################################################################



#######################################################################
###
###   User defined distribution for likelihood for
###   Uninfected radio-marked deer interval censor:
###   Test neg at cap and censoring
###
###   d_fit_sus_cens_postno
###
#######################################################################

dSusCensNo <- nimble::nimbleFunction(
    run = function(
       ### argument type declarations
        x = double(),
        e = double(0), #e, age of entry
        r = double(0), #r, age of last known alive
        age2date_surv = double(0),
        beta0_sus = double(0),
        beta0_inf = double(0),
        age_effect_surv = double(1),
        period_effect_surv = double(1),
        sex = double(0),
        beta_sex = double(0),
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

    lik<-0 #intialize likelihood
    llik<-0 #intialize log-likelihood
    lam_foi <- nimNumeric(r - 1)
    lam_sus <- nimNumeric(r - 1)
    lam_inf <- nimNumeric(r - 1)

    #############################################
    # preliminary hazards for the likelihood
    #############################################
    indx_sus_age <- e:(r - 1)
    n_indx_sus <- r - e
    indx_sus_period <- (e + age2date_surv):(r - 1 + age2date_surv)

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
        period_effect_surv[(1 + age2date_surv):(r - 1 + age2date_surv)] +
        rep(beta_sex * sex, (r - 1))
        )

    #force of infection infection hazard
    indx_foi_age_f <- age_lookup_f[1:(r - 1)]
    indx_foi_age_m <- age_lookup_m[1:(r - 1)]
    indx_foi_period <- period_lookup[(1 + age2date_foi):(s - 1 + age2date_foi)]

    lam_foi[1:(r - 1)] <- exp(rep(space, r - 1) +
            sex * (f_age_foi[indx_foi_age_f[1:(r - 1)]] +
                   f_period_foi[indx_foi_period[1:(r - 1)]])
            (1 - sex) * (m_age_foi[indx_foi_age_m[1:(r - 1)]] +
                         m_period_foi[indx_foi_period[1:(r - 1)]])
    )
    #######################################
    ### calculating the joint likelihood
    #######################################

    liktemp[(e + 1)] <- lam_foi[(e + 1)] * exp(-sum(lam_inf[(e + 1):(r - 1)]))

    for(k in (e + 2):(r - 1)){
        liktemp[k] <- lam_foi[k] *
                      exp(-sum(lam_sus[e:(k-1)])) *
                      exp(-sum(lam_inf[k:(r - 1)]))
    }

    lik <- exp(-lam_sus[e]) *
           exp(-sum(lam_foi[1:e])) *
           sum(liktemp[(e + 1):(r - 1)]) +
           exp(-sum(lam_foi[1:(r - 1)])) * exp(-sum(lam_sus[e:(r - 1)]))

    llik <- log(lik)

    returnType(double(0))
    if(log) return(llik) else return(exp(llik))    ## return log-likelihood
  })


nimble::registerDistributions(list(
    dSusCensNo = list(
        BUGSdist = 'dSusCensNo(e,r,age2date_surv,beta0_inf,beta0_sus,period_effect_surv,age_effect_surv,sex,beta_sex,f_age_foi,m_age_foi,age_lookup_f,age_lookup_m,f_period_foi,m_period_foi,period_lookup,space)',
        types = c("x = double(0)",
                  "e = double(0)",
                  "r = double(0)",
                  "age2date_surv = double(0)",
                  "beta0_inf = double(0)",
                  "beta0_sus = double(0)",
                  "period_effect_surv = double(1)",
                  "age_effect_surv = double(1)",
                  "sex = double(0)",
                  "beta_sex = double(0)",
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
assign('dSusCensNo', dSusCensNo, envir = .GlobalEnv)



#######################################################################
###
###   User defined distribution for likelihood for
###   Uninfected radio-marked deer interval censor:
###   Test neg at cap and censoring
###
###   d_fit_sus_cens_posttest
###
#######################################################################

dSusCensTest <- nimble::nimbleFunction(
    run = function(
       ### argument type declarations
        x = double(),
        e = double(0), #e, age of entry
        r = double(0), #r, age of last known alive
        age2date_surv = double(0),
        beta0_sus = double(0),
        age_effect_surv = double(1),
        period_effect_surv = double(1),
        sex = double(0),
        beta_sex = double(0),
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

    lik<-0 #intialize likelihood
    llik<-0 #intialize log-likelihood
    lam_foi <- nimNumeric(r - 1)
    lam_sus <- nimNumeric(r - 1)

    #############################################
    # preliminary hazards for the likelihood
    #############################################
    indx_sus_age <- e:(r - 1)
    n_indx_sus <- r - e
    indx_sus_period <- (e + age2date_surv):(r - 1 + age2date_surv)

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
    indx_foi_period <- period_lookup[(1 + age2date_foi):(s - 1 + age2date_foi)]

    lam_foi[1:(s - 1)] <- exp(rep(space, s - 1) +
            sex * (f_age_foi[indx_foi_age_f[1:(r - 1)]] +
                   f_period_foi[indx_foi_period[1:(r - 1)]])
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
    dSusCens = list(
        BUGSdist = 'dSusCensTest(e,r,age2date_surv,sex,beta0_sus,period_effect_surv,age_effect_surv,,beta_sex,f_age_foi,m_age_foi,age_lookup_f,age_lookup_m,f_period_foi,m_period_foi,period_lookup,space)',
        types = c("x = double(0)",
                  "e = double(0)",
                  "r = double(0)",
                  "age2date_surv = double(0)",
                  "beta0_sus = double(0)",
                  "period_effect_surv = double(1)",
                  "age_effect_surv = double(1)",
                  "sex = double(0)",
                  "beta_sex = double(0)",
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


#######################################################################
###
###   User defined distribution for likelihood for
###   uninfected radio-marked deer mortalities:
###   test neg at cap and tested mort
###
###   d_fit_sus_mort_posttest
###
#######################################################################


dSusMortTest <- nimble::nimbleFunction(
    run = function(
       ### argument type declarations
        x = double(),
        e = double(0), #e, age of entry
        r = double(0), #r, age of last known alive
        s = double(0), #s, age of known mortality
        age2date_surv = double(0),
        beta0_sus = double(0),
        age_effect_surv = double(1),
        period_effect_surv = double(1),
        sex = double(0),
        beta_sex = double(0),
        age2date_foi = double(0),
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
    lam_foi <- nimNumeric(s - 1)
    lam_sus <- nimNumeric(s - 1)

    #############################################
    # preliminary hazards for the likelihood
    #############################################
    n_indx_sus <- s - e
    #survival hazard for susceptible deer
    lam_sus[e:(s - 1)] <- exp(rep(beta0_sus, n_indx_sus) +
            age_effect_surv[e:(s - 1)] +
            period_effect_surv[(e + age2date_surv):(s - 1 + age2date_surv)] +
            rep(beta_sex * sex, n_indx_sus)
    )

    #force of infection infection hazard
    indx_foi_age_f <- age_lookup_f[1:(s - 1)]
    indx_foi_age_m <- age_lookup_m[1:(s - 1)]
    indx_foi_period <- period_lookup[(1 + age2date_foi):(s - 1 + age2date_foi)]

    lam_foi[1:(s - 1)] <- exp(rep(space, s - 1) +
            sex * (f_age_foi[indx_foi_age_f[1:(s - 1)]] +
                   f_period_foi[indx_foi_period[1:(s - 1)]])
            (1 - sex) * (m_age_foi[indx_foi_age_m[1:(s - 1)]] +
                         m_period_foi[indx_foi_period[1:(s - 1)]])
    )

    #######################################
    ### calculating the joint likelihood
    #######################################

    lik <- (1 - exp(-sum(lam_sus[r:(s - 1)]))) *
            exp(-sum(lam_sus[e:(r - 1)])) *
            exp(-sum(lam_foi[1:(s - 1)]))

    llik <- log(lik)

    returnType(double(0))
    if(log) return(llik) else return(exp(llik))    ## return log-likelihood
  })


nimble::registerDistributions(list(
    dSusMortTest = list(
        BUGSdist = 'dSusMortTest(e,r,s,age2date_surv,beta0_inf,beta0_sus,period_effect_surv,age_effect_surv,sex,beta_sex,age2date_foi,f_age_foi,m_age_foi,age_lookup_f,age_lookup_m,f_period_foi,m_period_foi,period_lookup,space)',
        types = c("x = double(0)",
                  "e = double(0)",
                  "r = double(0)",
                  "s = double(0)",
                  "age2date_surv = double(0)",
                  "beta0_inf = double(0)",
                  "beta0_sus = double(0)",
                  "period_effect_surv = double(1)",
                  "age_effect_surv = double(1)",
                  "sex = double(0)",
                  "beta_sex = double(0)",
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
assign('dSusMortTest', dSusMortTest, envir = .GlobalEnv)



#######################################################################
###
###   User defined distribution for likelihood for 
###   uninfected radio-marked deer mortalities:
###   test neg at cap and tested mort
###
###   d_fit_sus_mort_postno
###
#######################################################################


dSusMortNoTest <- nimble::nimbleFunction(
    run = function(
        ### argument type declarations
        x = double(),
        e = double(0), #e, age of entry
        r = double(0), #r, age of last known alive
        s = double(0), #s, age of known mortality
        dn1 = double(0), #interval of last test negative
        beta0_sus = double(0),
        age_effect_surv = double(1),
        period_effect_surv = double(1),
        sex = double(0),
        beta_sex = double(0),
        age2date_foi = double(0),
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
    lam_foi <- nimNumeric(s - 1)
    lam_sus <- nimNumeric(s - 1)
    lam_inf <- nimNumeric(s - 1)
    
    ############################################################
    ###
    ### overall likelihood uses the following components
    ### exp(-sum(lam_sus[e:(s-1)]))
    ### exp(-sum(lam_foi[1:(s - 1)]))
    ### exp(-sum(lam_inf[x:(r-1)]))
    ###
    ############################################################

    #############################################
    # preliminary hazards for the likelihood
    #############################################

    n_indx_sus <- length(e:(s - 1))
    #survival hazard for susceptible deer
    lam_sus[e:(s - 1)] <- exp(rep(beta0_sus, n_indx_sus) +
                    age_effect_surv[e:(s - 1)] +
                    period_effect_surv[(e + age2date_surv):(s - 1 + age2date_surv)] +
                    rep(beta_sex * sex, n_indx_sus))

    #survival hazard while infected
    n_indx_inf <- length(e:(s - 1))
    lam_inf[e:(s - 1)] <- exp(rep(beta0_inf, n_indx_inf) +
        age_effect_surv[e:(s - 1)] +
        period_effect_surv[(e + age2date_surv):(s - 1 + age2date_surv)] +
        rep(beta_sex * sex, n_indx_inf)
        )

    #force of infection infection hazard
    lam_foi[1:(s-1)] <- exp(rep(space, s - 1) +
        sex * (f_age_foi[age_lookup_f[1:(s - 1)]] +
            f_period_foi[period_lookup[(1 + age2date_foi):(s - 1 + age2date_foi)]]) +
        (1 - sex) * (m_age_foi[age_lookup_m[1:(e - 1)]] +
            m_period_foi[period_lookup[(1 + age2date_foi):(s - 1 + age2date_foi)]])
        )

    #total probability of getting infected and dying before end of the study
    liktemp[dn1+1] <- lam_foi[dn1+1] * 
                      exp(-sum(lam_inf[(dn1 + 1):(s - 1)])) 

    for(k in (dn1+2):(s - 1)){
        liktemp[k] <- lam_foi[k] *
                      exp(-sum(lam_sus[dn1:(k - 1)]))*
                      exp(-sum(lam_inf[(k:(r - 1))]))
    }

    lik <- exp(-sum(lam_sus[e:dn1])) *
           exp(-sum(lam_foi[1:dn1])) *
           (1 - exp(-lam_inf[r:(s - 1)])) *
           sum(liktemp[(dn1 + 1):(s - 1)]) +
           exp(-sum(lam_foi[1:(s - 1)])) *
           exp(-sum(lam_sus[1:(r - 1)])) *
           (1 - exp(-sum(lam_sus[r:(s - 1)])))

    llik <- log(lik)

    returnType(double(0))
    if(log) return(llik) else return(exp(llik))
  })

nimble::registerDistributions(list(
    dSusMortNoTest = list(
        BUGSdist = 'dSusMortNoTest(e,r,s,dn1,age2date_surv,beta0_inf,beta0_sus,period_effect_surv,age_effect_surv,sex,beta_sex,age2date_foi,f_age_foi,m_age_foi,age_lookup_f,age_lookup_m,f_period_foi,m_period_foi,period_lookup,space)',
        types = c("x = double(0)",
                    "e = double(0)",
                    "r = double(0)",
                    "dn1 = double(0)",
                    "age2date_surv = double(0)",
                    "beta0_inf = double(0)",
                    "beta0_sus = double(0)",
                    "period_effect_surv = double(1)",
                    "age_effect_surv = double(1)",
                    "sex = double(0)",
                    "beta_sex = double(0)",
                    "age2date_foi = double(0)",
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



#######################################################################
###
###   User defined distribution for likelihood for
###   infected deer mortalities for radio marked deer that
###   enter the study as test positive at capture
###
###   d_fit_icap_cens
###
#######################################################################

dIcapCens <- nimble::nimbleFunction(
    run = function(
        ### argument type declarations
        x = double(),
        e = double(0), #e, age of entry
        r = double(0), #r, age of last known alive
        dn1 = double(0), #interval of last test negative
        age2date_surv = double(0),
        beta0_sus = double(0),
        beta0_inf = double(0),
        age_effect_surv = double(1),
        period_effect_surv = double(1),
        sex = double(0),
        beta_sex = double(0),
        age2date_foi = double(0),
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
    lam_inf <- nimNumeric(r)
    lam_foi <- nimNumeric(r)
    lam_sus <- nimNumeric(r)
    
    #############################################
    # preliminary hazards for the likelihood
    #############################################

    n_indx_sus <- length(1:(e - 2))
    #survival hazard for susceptible deer
    lam_sus[1:(e - 2)] <- exp(rep(beta0_sus, n_indx_sus) +
                    age_effect_surv[1:(e - 2)] +
                    period_effect_surv[(e + age2date_surv):(e - 2 + age2date_surv)] +
                    rep(beta_sex * sex, n_indx_sus))

    #survival hazard while infected
    n_indx_inf <- length(1:(r - 1))
    lam_inf[1:(r - 1)] <- exp(rep(beta0_inf, n_indx_inf) +
        age_effect_surv[1:(r - 1)] +
        period_effect_surv[(1 + age2date_surv):(r - 1 + age2date_surv)] +
        rep(beta_sex * sex, n_indx_inf)
        )

    #force of infection infection hazard
    lam_foi[1:(e-1)] <- exp(rep(space, e - 1) +
        sex * (f_age_foi[age_lookup_f[1:(e - 1)]] +
            f_period_foi[period_lookup[(1 + age2date_foi):(e - 1 + age2date_foi)]]) +
        (1 - sex) * (m_age_foi[age_lookup_m[1:(e - 1)]] +
            m_period_foi[period_lookup[(1 + age2date_foi):(e - 1 + age2date_foi)]])
        )
    
    #######################################
    ### calculating the joint likelihood
    #######################################

    for(k in 1:(e - 1)) {
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
        BUGSdist = 'dIcapCens(e,r,age2date_surv,beta0_inf,beta0_sus,period_effect_surv,age_effect_surv,sex,beta_sex,age2date_foi,f_age_foi,m_age_foi,age_lookup_f,age_lookup_m,f_period_foi,m_period_foi,period_lookup,space)',
        types = c("x = double(0)",
                    "e = double(0)", 
                    "r = double(0)",
                    "age2date_surv = double(0)",
                    "beta0_inf = double(0)",
                    "beta0_sus = double(0)",
                    "period_effect_surv = double(1)",
                    "age_effect_surv = double(1)",
                    "sex = double(0)",
                    "beta_sex = double(0)",
                    "age2date_foi = double(0)", 
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

#######################################################################
###
###   User defined distribution for likelihood for
###   infected deer mortalities for radio marked deer that
###   enter the study as test positive at capture
###
###   d_fit_icap_mort
###
#######################################################################

dIcapMort <- nimble::nimbleFunction(
    run = function(
        ### argument type declarations
        x = double(0),
        e = double(0), #e, age of entry
        r = double(0), #r, age of last known alive
        s = double(0), #s, age of known mortality
        age2date_surv = double(0),
        beta0_sus = double(0),
        beta0_inf = double(0),
        age_effect_surv = double(1),
        period_effect_surv = double(1),
        sex = double(0), 
        beta_sex = double(0),
        age2date_foi = double(0),
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
    lam_inf <- nimNumeric(s)
    lam_foi <- nimNumeric(s)
    lam_sus <- nimNumeric(s)
    
    #############################################
    # preliminary hazards for the likelihood
    #############################################

    n_indx_sus <- length(1:(e-2))
    #survival hazard for susceptible deer
    lam_sus[1:(e - 2)] <- exp(rep(beta0_sus, n_indx_sus)  + 
                    age_effect_surv[1:(e - 2)] +
                    period_effect_surv[(e + age2date_surv):(e - 2 + age2date_surv)] +
                    rep(beta_sex * sex, n_indx_sus))

    #survival hazard while infected
    n_indx_inf <- length(1:(s - 1))
    lam_inf[1:(s - 1)] <- exp(rep(beta0_inf, n_indx_inf) +
        age_effect_surv[1:(s - 1)] +
        period_effect_surv[(1 + age2date_surv):(s - 1 + age2date_surv)] +
        rep(beta_sex * sex, n_indx_inf)
        )

    #force of infection infection hazard
    lam_foi[1:(e-1)] <- exp(rep(space, e - 1) +
        sex * (f_age_foi[age_lookup_f[1:(e - 1)]] +
            f_period_foi[period_lookup[(1 + age2date_foi):(e - 1 + age2date_foi)]]) +
        (1 - sex) * (m_age_foi[age_lookup_m[1:(e - 1)]] +
            m_period_foi[period_lookup[(1 + age2date_foi):(e - 1 + age2date_foi)]])
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

    lik <- exp(-sum(lam_inf[e:(r - 1)])) *
           (1 - exp(-sum(lam_inf[r:(s - 1)]))) *
           sum(lik_temp[1:(e - 1)])
    llik <- log(lik)

    returnType(double(0))
    if(log) return(llik) else return(exp(llik))    ## return log-likelihood
  })


nimble::registerDistributions(list(
    dIcapMort = list(
        BUGSdist = 'dIcapMort(e,r,s,age2date_surv,beta0_inf,beta0_sus,period_effect_surv,age_effect_surv,sex,beta_sex,age2date_foi,f_age_foi,m_age_foi,age_lookup_f,age_lookup_m,f_period_foi,m_period_foi,period_lookup,space)',
        types = c("x = double(0)",
                    "e = double(0)", 
                    "r = double(0)", 
                    "s = double(0)", 
                    "age2date_surv = double(0)",
                    "beta0_inf = double(0)",
                    "beta0_sus = double(0)",
                    "period_effect_surv = double(1)",
                    "age_effect_surv = double(1)",
                    "sex = double(0)",
                    "beta_sex = double(0)",
                    "age2date_foi = double(0)",
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
assign('dIcapMort', dIcapMort, envir = .GlobalEnv)


#######################################################################
###
###   User defined distribution for likelihood for
###   infected deer that were test neg at capture,
###   then test negative at recap,
###   that die
###
###   d_fit_rec_neg_mort
###
#######################################################################

dRecNegMort <- nimble::nimbleFunction(
    run = function(
        ### argument type declarations
        x = double(0),
        e = double(0), #e, age of entry
        r = double(0), #r, age of last known alive
        s = double(0), #s, age of mortality
        dn1 = double(0), #interval of last test negative
        age2date_surv = double(0),
        sex = double(0),
        beta_sex = double(0),
        age2date_foi = double(0),
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
    lam_foi <- nimNumeric(s)
    lam_sus <- nimNumeric(s)
    lam_inf <- nimNumeric(s)

    #############################################
    # preliminary hazards for the likelihood
    #############################################

    inx_period_surv <- (e + age2date_surv):(s - 1 + age2date_surv)
    #survival hazard for susceptible deer
    lam_sus[1:(s - e)] <- exp(rep(beta0_sus, s - e) +
                          age_effect_surv[1:(s - e)] +
                          period_effect_surv[indx_period_surv[1:(s - e)]] +
                          rep(beta_sex * sex, s - e)
                      )
    #force of infection infection hazard
    indx_foi_age_f <- age_lookup_f[1:(s - 1)]
    indx_foi_age_m <- age_lookup_m[1:(s - 1)]
    indx_foi_period <- period_lookup[(1 + age2date_foi):(r - 1 + age2date_foi)]

    lam_foi[1:(s - 1)] <- exp(rep(space, s - 1) +
            sex * (f_age_foi[indx_foi_age_f[1:(s - 1)]] +
                   f_period_foi[indx_foi_period[1:(s - 1)]])
            (1 - sex) * (m_age_foi[indx_foi_age_m[1:(s - 1)]] +
                         m_period_foi[indx_foi_period[1:(s - 1)]])
    )
    #######################################
    ### calculating the joint likelihood
    #######################################
    for (k in 1:(dn1 - 1)) {
        liktemp[k] <- lambda_foi[k] * exp(-sum(lam_sus[1:k]))
    }
    lik <- exp(-sum(lam_sus[e:(r - 1)])) *
           (1 - exp(sum(lam_sus[(r:(s - 1))]))) *
           exp(-sum(lambda_foi[dn1:(s - 1)])) *
           sum(liktemp[1:(dn1 - 1)])

    llik <- log(lik)

    returnType(double(0))
    if(log) return(llik) else return(exp(llik))    ## return log-likelihood
  })

nimble::registerDistributions(list(
    dRecNegMort = list(
        BUGSdist = 'dRecNegMort(e,r,s,dn1,age2date_surv,beta0_inf,beta0_sus,period_effect_surv,age_effect_surv,sex,beta_sex,age2date_foi,f_age_foi,m_age_foi,age_lookup_f,age_lookup_m,f_period_foi,m_period_foi,period_lookup,space)',
        types = c("x = double(0)",
                  "e = double(0)",
                  "r = double(0)",
                  "dn = double(0)",
                  "age2date_surv = double(0)",
                  "beta0_inf = double(0)",
                  "beta0_sus = double(0)",
                  "period_effect_surv = double(1)",
                  "age_effect_surv = double(1)",
                  "sex = double(0)",
                  "beta_sex = double(0)",
                  "age2date_foi = double(0)",
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
assign('dRecNegMort', dRecNegMort, envir = .GlobalEnv)


#######################################################################
###
###   User defined distribution for likelihood for
###   infected deer that were test neg at capture,
###   then test negative at recap, that are interval censored, 
###   and have been tested post censoring
###
###   d_fit_rec_neg_cens
###
#######################################################################

dRecNegCens <- nimble::nimbleFunction(
    run = function(
        ### argument type declarations
        x = double(0),
        e = double(0), #e, age of entry
        r = double(0), #r, age of last known alive
        dn1 = double(0), #interval of last test negative
        age2date_surv = double(0),
        sex = double(0),
        beta_sex = double(0),
        age2date_foi = double(0),
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

    #############################################
    # preliminary hazards for the likelihood
    #############################################

    inx_period_surv <- (e + age2date_surv):(r - 1 + age2date_surv)

    #survival hazard for susceptible deer
    lam_sus[1:(r - e)] <- exp(rep(beta0_sus, r - e) +
                          age_effect_surv[1:(r - e)] +
                          period_effect_surv[indx_period_surv[1:(r - e)]] +
                          rep(beta_sex * sex, r - e)
                          )

    #force of infection infection hazard
    indx_foi_age_f <- age_lookup_f[1:(r - 1)]
    indx_foi_age_m <- age_lookup_m[1:(r - 1)]
    indx_foi_period <- period_lookup[(1 + age2date_foi):(r - 1 + age2date_foi)]

    lam_foi[1:(r - 1)] <- exp(rep(space, r - 1) +
            sex * (f_age_foi[indx_foi_age_f[1:(r - 1)]] +
                   f_period_foi[indx_foi_period[1:(r - 1)]])
            (1 - sex) * (m_age_foi[indx_foi_age_m[1:(r - 1)]] +
                         m_period_foi[indx_foi_period[1:(r - 1)]])
    )

    #######################################
    ### calculating the joint likelihood
    #######################################

    for (k in 1:(dn1 - 1)) {
        liktemp[k] <- lambda_foi[k] * exp(-sum(lam_sus[1:k]))
    }
    lik <- sum(liktemp[1:(dn - 1)]) *
           exp(-sum(lam_sus[e:(r - 1)])) *
           exp(-sum(lam_foi[1:(r - 1)]))
    llik <- log(lik)

    returnType(double(0))
    if(log) return(llik) else return(exp(llik))    ## return log-likelihood
  })

nimble::registerDistributions(list(
    dRecNegCens = list(
        BUGSdist = 'dRecNegCens(e,r,dn1,age2date_surv,beta0_inf,beta0_sus,period_effect_surv,age_effect_surv,sex,beta_sex,age2date_foi,f_age_foi,m_age_foi,age_lookup_f,age_lookup_m,f_period_foi,m_period_foi,period_lookup,space)',
        types = c("x = double(0)",
                    "e = double(0)",
                    "r = double(0)",
                    "dn1 = double(0)",
                    "age2date_surv = double(0)",
                    "beta0_inf = double(0)",
                    "beta0_sus = double(0)",
                    "period_effect_surv = double(1)",
                    "age_effect_surv = double(1)",
                    "sex = double(0)",
                    "beta_sex = double(0)",
                    "age2date_foi = double(0)",
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
assign('dRecNegCens', dRecNegCens, envir = .GlobalEnv)

#######################################################################
###
###   User defined distribution for likelihood for
###   infected deer that were test neg at capture,
###   then test positive at recap,
###   than die
###
###   d_fit_rec_pos_mort
###
#######################################################################


dRecPosMort <- nimble::nimbleFunction(
    run = function(
        ### argument type declarations
        x = double(0),
        e = double(0), #e, age of entry
        r = double(0), #r, age of last known alive
        s = double(0), #s, age of last known alive
        dn1 = double(0), #interval of last test negative
        dn = double(0), #int tested positive
        age2date_surv = double(0),
        sex = double(0),
        beta_sex = double(0),
        age2date_foi = double(0),
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
    lam_foi <- nimNumeric(s)
    lam_sus <- nimNumeric(s)
    lam_inf <- nimNumeric(s)

    #############################################
    # preliminary hazards for the likelihood
    #############################################

    n_indx_sus <- length(e:(dn - 1))

    #survival hazard for susceptible deer
    lam_sus[e:(dn - 1)] <- exp(rep(beta0_sus, n_indx_sus)  + 
                    age_effect_surv[e:(dn - 1)] +
                    period_effect_surv[(e + age2date_surv):(dn - 1 + age2date_surv)] +
                    rep(beta_sex * sex, n_indx_sus))
    
    #survival hazard while infected
    n_indx_inf <- length(dn1:(s - 1))
    
    lam_inf[dn1:(s - 1)] <- exp(rep(beta0_inf,n_indx_inf) +
                                age_effect_surv[dn1:(s - 1)] +
                                period_effect_surv[(dn1 + age2date_surv):
                                                   (s - 1 + age2date_surv)] +
                                rep(beta_sex * sex,n_indx_inf)
                            )
    #force of infection infection hazard
    lam_foi[1:dn] <- exp(rep(space,dn) +
                  sex * (f_age_foi[age_lookup_f[1:dn]] +
                            f_period_foi[period_lookup[(1 + age2date_foi):(dn + age2date_foi)]]) +
                  (1 - sex) * (m_age_foi[age_lookup_m[1:dn]] +
                                  m_period_foi[period_lookup[(1 + age2date_foi):(dn + age2date_foi)]])
                  )

    #######################################
    ### calculating the joint likelihood
    #######################################

    for(k in dn1:dn) {
      lik_temp[x] <- lam_foi[k] *
               exp(-sum(lam_sus[e:(k - 1)])) *
               exp(-sum(lam_foi[1:(k - 1)])) *
               exp(-sum(lam_inf[x:(r - 1)]))
    }

    lik <- (1 - exp(-sum(lam_inf[r:(s - 1)]))) * sum(lik_temp[dn1:dn])
    llik <- log(lik)

    returnType(double(0))
    if(log) return(llik) else return(exp(llik))    ## return log-likelihood
  })

nimble::registerDistributions(list(
    dRecPosMort = list(
        BUGSdist = 'dRecPosMort(e,r,s,dn1,dn,age2date_surv,beta0_inf,beta0_sus,period_effect_surv,age_effect_surv,sex,beta_sex,age2date_foi,f_age_foi,m_age_foi,age_lookup_f,age_lookup_m,f_period_foi,m_period_foi,period_lookup,space)',
        types = c("x = double(0)",
                    "e = double(0)",
                    "r = double(0)",
                    "s = double(0)",
                    "dn1 = double(0)",
                    "dn = double(0)",
                    "age2date_surv = double(0)",
                    "beta0_inf = double(0)",
                    "beta0_sus = double(0)",
                    "period_effect_surv = double(1)",
                    "age_effect_surv = double(1)",
                    "sex = double(0)",
                    "beta_sex = double(0)",
                    "age2date_foi = double(0)",
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


#######################################################################
###
###   User defined distribution for likelihood for 
###   infected deer that were test neg at capture,
###   then test positive at recap,
###   that are interval censored
###
###   d_fit_rec_pos_cens
###
#######################################################################

dRecPosCens <- nimble::nimbleFunction(
    run = function(
        ### argument type declarations
        x = double(),
        e = double(0), #e, age of entry
        r = double(0), #r, age of last known alive
        dn1 = double(0), #interval of last test negative
        dn = double(0), #int tested positive
        age2date_surv = double(0),
        sex = double(0),
        beta_sex = double(0),
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
    lam_foi <- nimNumeric(r)
    lam_sus <- nimNumeric(r)
    lam_inf <- nimNumeric(r)

    #############################################
    # preliminary hazards for the likelihood
    #############################################

    n_indx_sus <- length(e:(dn - 1))

    #survival hazard for susceptible deer
    lam_sus[e:(dn - 1)] <- exp(rep(beta0_sus, n_indx_sus)  + 
                    age_effect_surv[e:(dn - 1)] +
                    period_effect_surv[(e + age2date_surv):(dn - 1 + age2date_surv)] +
                    rep(beta_sex * sex, n_indx_sus))
    
    #survival hazard while infected
    n_indx_inf <- length(dn1:(r - 1))
    
    lam_inf[dn1:(r - 1)] <- exp(rep(beta0_inf,n_indx_inf) +
        age_effect_surv[dn1:(r - 1)] +
        period_effect_surv[(dn1 + age2date_surv):(r - 1 + age2date_surv)] +
        rep(beta_sex * sex,n_indx_inf)
        )

    #force of infection infection hazard
    lam_foi[1:(dn-1)] <- exp(rep(space, dn - 1) +
        sex * (f_age_foi[age_lookup_f[1:(dn - 1)]] +
            f_period_foi[period_lookup[(1 + age2date_foi):(dn - 1 + age2date_foi)]]) +
        (1 - sex) * (m_age_foi[age_lookup_m[1:(dn - 1)]] +
            m_period_foi[period_lookup[(1 + age2date_foi):(dn - 1 + age2date_foi)]])
        )

    #######################################
    ### calculating the joint likelihood
    #######################################

    for(k in dn1:(dn - 1)) {
      lik_temp[k] <- lam_foi[k] *
               exp(-sum(lam_sus[e:(k-1)])) *
               exp(-sum(lam_foi[1:(k-1)])) *
               exp(-sum(lam_foi[k:(dn-1)]))     
    }

    lik <- exp(-sum(lam_inf[dn:(r-1)])) * sum(lik_temp[dn1:(dn-1)])
    llik <- log(lik)

    returnType(double(0))
    if(log) return(llik) else return(exp(llik))    ## return log-likelihood
  })

nimble::registerDistributions(list(
    dRecPosCens = list(
        BUGSdist = 'dRecPosCens(e,r,dn1,dn,age2date_surv,beta0_inf,beta0_sus,period_effect_surv,age_effect_surv,sex,beta_sex,f_age_foi,m_age_foi,age_lookup_f,age_lookup_m,f_period_foi,m_period_foi,period_lookup,space)',
        types = c("x = double(0)",
                    "e = double(0)", 
                    "r = double(0)", 
                    "dn1 = double(0)", 
                    "dn = double(0)", 
                    "age2date_surv = double(0)",
                    "beta0_inf = double(0)",
                    "beta0_sus = double(0)",
                    "period_effect_surv = double(1)",
                    "age_effect_surv = double(1)",
                    "sex = double(0)",
                    "beta_sex = double(0)",
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





#######################################################################
###
###   User defined distribution for likelihood for
###   infected deer mortalities for radio marked deer that
###   enter the study as test negative at capture
###
###   d_fit_idead
###
#######################################################################

dNegCapPosMort <- nimble::nimbleFunction(
    run = function(
        ### argument type declarations
        x = double(),
        e = double(0), #e, age of entry
        r = double(0), #r, age of last known alive
        s = double(0), #s, age of known mortality
        dn1 = double(0), #last interval of test negative
        dn = double(0), #interval of test positive (could be interval of mortality)
        age2date_surv = double(0),
        sex = double(0),
        beta_sex = double(0),
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
    lam_inf <- nimNumeric(s)
    lam_foi <- nimNumeric(s)
    lam_sus <- nimNumeric(s)
    
    #############################################
    # preliminary hazards for the likelihood
    #############################################
    

    n_indx_sus <- length(e:(dn - 1))

    #survival hazard for susceptible deer
    lam_sus[e:(dn - 1)] <- exp(rep(beta0_sus, n_indx_sus)  + 
                    age_effect_surv[e:(dn - 1)] +
                    period_effect_surv[(e + age2date_surv):(dn - 1 + age2date_surv)] +
                    rep(beta_sex * sex, n_indx_sus))
    
    #survival hazard while infected
    n_indx_inf <- length(dn1:(s - 1))
    lam_inf[dn1:(s - 1)] <- exp(rep(beta0_inf, n_indx_inf) +
        age_effect_surv[dn1:(s - 1)] +
        period_effect_surv[(dn1 + age2date_surv):(s - 1 + age2date_surv)] +
        rep(beta_sex * sex, n_indx_inf)
        )

    #force of infection infection hazard
    lam_foi[1:(dn-1)] <- exp(rep(space, dn - 1) +
        sex * (f_age_foi[age_lookup_f[1:(dn - 1)]] +
            f_period_foi[period_lookup[(1 + age2date_foi):(dn - 1 + age2date_foi)]]) +
        (1 - sex) * (m_age_foi[age_lookup_m[1:(dn - 1)]] +
            m_period_foi[period_lookup[(1 + age2date_foi):(dn - 1 + age2date_foi)]])
        )

    #######################################
    ### calculating the joint likelihood
    #######################################

    for(k in dn1:r){
     lik_temp[k] <- lam_foi[k] *
               exp(-sum(lam_sus[e:(k-1)])) *
               exp(-sum(lam_foi[1:(k-1)])) *
               exp(-sum(lam_inf[k:(r-1)]))
    }
    lik_temp[dn] <- lam_foi[dn] *
               exp(-sum(lam_sus[e:(dn-1)])) *
               exp(-sum(lam_foi[1:(dn-1)])) *

    lik <- (1 - exp(-sum(lam_inf[r:(s-1)])))*
            sum(lik_temp[dn1:dn])
    llik <- log(lik)

    returnType(double(0))
    if(log) return(llik) else return(exp(llik))    ## return log-likelihood
  })


nimble::registerDistributions(list(
    dNegCapPosMort = list(
        BUGSdist = 'dNegCapPosMort(e,r,s,dn1,dn,age2date_surv,beta0_inf,beta0_sus,period_effect_surv,age_effect_surv,sex,beta_sex,f_age_foi,m_age_foi,age_lookup_f,age_lookup_m,f_period_foi,m_period_foi,period_lookup,space)',
        types = c("e = double(0)", 
                    "r = double(0)", 
                    "s = double(0)", 
                    "dn1 = double(0)", 
                    "dn = double(0)", 
                    "age2date_surv = double(0)",
                    "beta0_inf = double(0)",
                    "beta0_sus = double(0)",
                    "period_effect_surv = double(1)",
                    "age_effect_surv = double(1)",
                    "sex = double(0)",
                    "beta_sex = double(0)",
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
assign('dNegCapPosMort', dNegCapPosMort, envir = .GlobalEnv)


