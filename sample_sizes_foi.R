sample_size_surveillance <- cwd_df %>% count(yearkill)

sample_size_surveillance$year <- 2002:2022

sample_size_surveillance
sample_size_surveillance$n[20] <- sample_size_surveillance$n[20]+sample_size_surveillance$n[21]
sample_size_surveillance <- sample_size_surveillance[1:20,]

col_n <- d_col_temp %>% count(year_cap)
col_n <- c(col_n$n,0)
compare_samples <- data.frame(year = sample_size_surveillance$year[sample_size_surveillance$year>2016],
                              col_n=col_n,
                              surv_n = sample_size_surveillance$n[sample_size_surveillance$year>2016]
                                )


#run FOI model without surveillance data

names(d_fit_col)
length(unique(d_fit_col$year_cap))
n_period  <- 5


dFOIcollar <- nimble::nimbleFunction(
    run = function(
        ### argument type declarations
        x=double(),
        left = double(0),
        right = double(0),
        sexfoi = double(0),
        f_age = double(1),
        m_age = double(1),
        age_lookup = double(1),
        space = double(0),
        log = double()
        ) {

    logL<-0 #intialize log-likelihood
    gam <-nimNumeric(right)
    for (k in left:(right-1)) {
        gam[k] <- space +
                  sexfoi * (f_age[age_lookup[k]]) +
                  (1 - sexfoi) * (m_age[age_lookup[k]])
    }
    #total probability of getting infected
    p <- 1 - exp(-sum(exp(gam[left:(right-1)])))
    
    logL <- dbinom(x,1,p,log=TRUE)

    returnType(double(0))
    if(log) return(logL) else return(exp(logL))    ## return log-likelihood
  })


nimble::registerDistributions(list(
    dFOIcollar = list(
        BUGSdist = 'dFOIcollar(left,right,sexfoi,f_age,m_age,age_lookup,space)',
        types = c('p = double(0)',
                  'left=double(0)',
                  'right=double(0)',
                  'sexfoi=double(0)',
                  'f_age=double(1)',
                  'm_age=double(1)',
                  'age_lookup=double(1)',
                  'space=double(0)'
                  ),
        discrete = TRUE
    )
))

# for a user-defined distribution
assign('dFOIcollar', dFOIcollar, envir = .GlobalEnv)


##################################################################################################################################
###
###   User defined distribution for FOI from GPS with antemortem tests + at capture
###
##################################################################################################################################


dFOIinf <- nimble::nimbleFunction(
    run = function(
        ### argument type declarations
        x=double(0),
        left = double(0),
        sexfoi = double(0),
        f_age = double(1),
        m_age = double(1),
        age_lookup = double(1),
        space = double(0),
        log = double()
        ) {

    logL<-0 #intialize log-likelihood
    gam <-nimNumeric(left-1)

    for (k in 1:(left-1)) {
        gam[k] <- space +
                  sexfoi * (f_age[age_lookup[k]]) +
                  (1 - sexfoi) * (m_age[age_lookup[k]])
    }
    #total probability of getting infected
    p <- 1 - exp(-sum(exp(gam[1:(left-1)])))
    
    logL <- dbinom(x,1,p,log=TRUE)

    returnType(double(0))
    if(log) return(logL) else return(exp(logL))    ## return log-likelihood
  })


nimble::registerDistributions(list(
    dFOIinf = list(
        BUGSdist = 'dFOIinf(left,sexfoi,f_age,m_age,age_lookup,space)',
        types = c('p = double(0)',
                  'gam = double(1)',
                  'left=double(0)',
                  'sexfoi=double(0)',
                  'f_age=double(1)',
                  'm_age=double(1)',
                  'age_lookup=double(1)',
                  'space=double(0)'
                  ),
        discrete = TRUE
    )
))

assign('dFOIinf', dFOIinf, envir = .GlobalEnv)




#########################
###
### Model code
###
##########################

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

    #Likelihood for Collared individuals that were test negative at capture
    for (i in 1:n_fit_col) {
      teststatus_col[i] ~ dFOIcollar(left = left_age_col[i],
                  right = right_age_col[i],
                  sexfoi = sexfoi_col[i],
                  f_age = f_age[1:n_agef],
                  m_age = m_age[1:n_agem],
                  age_lookup = age_lookup_col[1:n_age_lookup_col],
                  space = space[sect_col[i]]
                  )
    }

    #Likelihood for individuals infected at capture
    for (i in 1:n_fit_inf) {
      teststatus_inf[i] ~ dFOIinf(left = left_age_inf[i],
                  sexfoi = sexfoi_inf[i],
                  f_age = f_age[1:n_agef],
                  m_age = m_age[1:n_agem],
                  age_lookup = age_lookup_col_inf[1:n_age_lookup_col_inf],
                  space = space[sect_inf_col[i]]
                  )
    }

})#end model statement



#######################################
### Data for Model Fitting
#######################################

nimData <- list(
                teststatus_col = d_fit_col$censor,
                left_age_col = d_fit_col$left_age,
                right_age_col = d_fit_col$right_age,
                sexfoi_col = d_fit_col$sex,
                teststatus_inf = d_fit_inf$censor,
                left_age_inf = d_fit_inf$left_age,
                sexfoi_inf = d_fit_inf$sex
                # snum = num_sp,
                # sadj = adj_sp,
                # sweights = weights_sp
                # teststatus_hunt = cwd_df$teststatus,
                # ageweeks_hunt = cwd_df$ageweeks,
                # sexfoi_hunt = cwd_df$sex
                )


#######################################
### Constants for MCMC
#######################################

nimConsts <- list(n_agem = n_agem,
                  n_agef = n_agef,
                  n_period = n_period,
                  n_age = n_age,
                  n_fit_col = n_fit_col,
                  n_fit_inf = n_fit_inf,
                  n_age_lookup_col = n_age_lookup_col,
                  age_lookup_col = age_lookup_col,
                  n_age_lookup_col_inf = n_age_lookup_col_inf,
                  age_lookup_col_inf = age_lookup_col_inf,
                  n_sect = n_sect,
                  sect_col = d_fit_col$sect,
                  sect_inf_col = d_fit_inf$sect,
                  space = rep(0,n_sect)
                  )


#######################################
### Initial Values for MCMC
#######################################
initsFun <- function()list(
                          mprec = runif(1, 1.5, 1.7),
                          fprec = runif(1, 2.7, 4.2),
                        #   m_period = seq(-2, 2, length = n_period),
                        #   f_period = seq(-2, 2, length = n_period),
                          m_age = seq(-6, -4, length = n_agem),
                          f_age = seq(-7, -5, length = n_agef)
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
              "m_age_mu",
              "f_age_mu",
              "mprec",
              "fprec",
              "f_age",
              "m_age"
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
                  niter = 10000,
                  nburnin = 5000,
                  nchains = 3,
                  inits = initsFun,
                  samplesAsCodaMCMC = TRUE,
                  summary = TRUE
                  )

runtime <- difftime(Sys.time(),
                    starttime,
                    units = "min")

mcmcout_infonly <- mcmcout
save(mcmcout_infonly,"results/mcmcout_infonly.Rdata")

out <- mcmcout$samples
fit_sum <- mcmcout$summary$all.chains

gd <- gelman.diag(out[,c(grep("beta",rownames(fit_sum)),
                        grep("prec",rownames(fit_sum)),
                        grep("tau",rownames(fit_sum)))],
                  multivariate = FALSE)
gd
es <- effectiveSize(out[,c(grep("beta",rownames(fit_sum)),
                        grep("prec",rownames(fit_sum)),
                           grep("tau",rownames(fit_sum)))])
es
pdf("figures/traceplots_infonly.pdf")
traceplot(out[,"mprec"],ylab="mprec")
traceplot(out[,"fprec"],ylab="fprec")
for(i in 1:n_age){
  traceplot(out[,paste0("f_age[",i,"]")],ylab=paste0("f_age[",i,"]"))
}
for(i in 1:n_age){
  traceplot(out[,paste0("m_age[",i,"]")],ylab=paste0("m_age[",i,"]"))
}
dev.off()



Agegroups=1:n_agef
haz.mean = fit_sum[grep("f_age",rownames(fit_sum)),1][1:n_agef]
haz.lower = fit_sum[grep("f_age",rownames(fit_sum)),4][1:n_agef]
haz.upper = fit_sum[grep("f_age",rownames(fit_sum)),5][1:n_agef]

temp = data.frame(Agegroups,haz.mean,haz.lower,haz.upper)
temp <- rbind(temp,temp[nrow(temp),])
temp$Agegroups[8] <- 8

fage.plot=ggplot(data=temp)+geom_step(aes(x=Agegroups,y=haz.mean),direction="hv")+
  geom_step(aes(x=Agegroups,y=haz.lower),direction="hv",linetype=3)+
  geom_step(aes(x=Agegroups,y=haz.upper),direction="hv",linetype=3)+
  theme_bw()+ylim(-31,-2.5)+ggtitle("Female")+
  scale_x_continuous(breaks=(1:n_agef)+.5,labels=c("Fawn","1.5","2.5","3.5","4.5-5.5","6.5-8.5","9.5+"))

fage.plot


Agegroups=1:n_agem
haz.mean = fit_sum[grep("m_age",rownames(fit_sum)),1][1:n_agem]
haz.lower = fit_sum[grep("m_age",rownames(fit_sum)),4][1:n_agem]
haz.upper = fit_sum[grep("m_age",rownames(fit_sum)),5][1:n_agem]

temp = data.frame(Agegroups,haz.mean,haz.lower,haz.upper)
temp <- rbind(temp,temp[nrow(temp),])
temp$Agegroups[n_agem+1] <- n_agem+1

mage.plot=ggplot(data=temp)+geom_step(aes(x=Agegroups,y=haz.mean),direction="hv")+
  geom_step(aes(x=Agegroups,y=haz.lower),direction="hv",linetype=3)+
  geom_step(aes(x=Agegroups,y=haz.upper),direction="hv",linetype=3)+
  theme_bw()+
  ylim(-31,-2.5)+ggtitle("Male")+
  scale_x_continuous(breaks = (1:n_agem)+.5,labels =c("Fawn","1.5","2.5","3.5","4.5-5.5","6.5+"))

mage.plot

combo_infonly <- grid.arrange(mage.plot,fage.plot,ncol=2)

ggsave("figures/combo_infonly.png",combo_infonly,height=4,width=8)


mcmcout_inf <- mcmcout
save(mcmcout_inf, file = "mcmcout_inf.Rdata")



########################################################################
########################################################################
########################################################################



dFOIhunt <- nimble::nimbleFunction(
    run = function(
        ### argument type declarations
        x=double(),
        agew = double(0),
        sexfoi = double(0),
        f_age = double(1),
        m_age = double(1),
        age_lookup = double(1),
        space = double(0),
        log = double()
        ) {
    
    logL<-0 #intialize log-likelihood
    gam <-nimNumeric(agew)
    for (k in 1:agew) {
        gam[k] <- space +
                  sexfoi * (f_age[age_lookup[k]]) +
                  (1 - sexfoi) * (m_age[age_lookup[k]])
    }
    p <- 1 - exp(-sum(exp(gam[1:agew])))
    
    logL <- dbinom(x,1,p,log=TRUE)

    returnType(double(0))
    if(log) return(logL) else return(exp(logL))    ## return log-likelihood
  })


nimble::registerDistributions(list(
    dFOIhunt = list(
        BUGSdist = 'dFOIhunt(agew,sexfoi,f_age,m_age,age_lookup,space)',
        types = c('p = double(0)',
                  'gam = double(1)',
                  'agew=double(0)',
                  'sexfoi=double(0)',
                  'f_age=double(1)',
                  'm_age=double(1)',
                  'age_lookup=double(1)',
                  'space=double(0)'
                  ),
        discrete = TRUE
    )
))


assign('dFOIhunt', dFOIhunt, envir = .GlobalEnv)



##################################################
### Model code
###################################################



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

    #Likelihood for Collared individuals that were test negative at capture
    # for (i in 1:n_fit_col) {
    #   teststatus_col[i] ~ dFOIcollar(left = left_age_col[i],
    #               right = right_age_col[i],
    #               sexfoi = sexfoi_col[i],
    #               f_age = f_age[1:n_agef],
    #               m_age = m_age[1:n_agem],
    #               age_lookup = age_lookup_col[1:n_age_lookup_col],
    #               space = space[sect_col[i]]
    #               )
    # }

    # #Likelihood for individuals infected at capture
    # for (i in 1:n_fit_inf) {
    #   teststatus_inf[i] ~ dFOIinf(left = left_age_inf[i],
    #               sexfoi = sexfoi_inf[i],
    #               f_age = f_age[1:n_agef],
    #               m_age = m_age[1:n_agem],
    #               age_lookup = age_lookup_col_inf[1:n_age_lookup_col_inf],
    #               space = space[sect_inf_col[i]]
    #               )
    # }

    #Likelihood for hunter harvest surveillence data
    for (i in 1:n_fit_hunt) {
      teststatus_hunt[i] ~ dFOIhunt(agew = ageweeks_hunt[i],
                  sexfoi = sexfoi_hunt[i],
                  f_age = f_age[1:n_agef],
                  m_age = m_age[1:n_agem],
                  age_lookup = age_lookup[1:n_age_lookup],
                  space = space[sect_hunt[i]]
                  )
    }

})#end model statement



#######################################
### Data for Model Fitting
#######################################

cwd_df$kill_year[cwd_df$kill_year==2022] <- 2021
dim(cwd_df)
cwd_df <- cwd_df[cwd_df$kill_year>2016,]
n_fit_hunt <- nrow(cwd_df)

nimData <- list(
                # teststatus_col = d_fit_col$censor,
                # left_age_col = d_fit_col$left_age,
                # right_age_col = d_fit_col$right_age,
                # sexfoi_col = d_fit_col$sex,
                # teststatus_inf = d_fit_inf$censor,
                # left_age_inf = d_fit_inf$left_age,
                # sexfoi_inf = d_fit_inf$sex
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

nimConsts <- list(n_agem = n_agem,
                  n_agef = n_agef,
                  n_period = n_period,
                  n_age = n_age,
                  n_fit_col = n_fit_col,
                  n_fit_inf = n_fit_inf,
                  n_age_lookup_col = n_age_lookup_col,
                  age_lookup_col = age_lookup_col,
                  n_age_lookup_col_inf = n_age_lookup_col_inf,
                  age_lookup_col_inf = age_lookup_col_inf,
                  n_sect = n_sect,
                  sect_col = d_fit_col$sect,
                  sect_inf_col = d_fit_inf$sect,
                  n_age_lookup = n_age_lookup,
                  age_lookup = age_lookup,
                  period_lookup_hunt = period_lookup_hunt,
                  n_period_lookup_hunt = n_period_lookup_hunt,
                  birthweek_hunt = cwd_df$birthweek,
                  space = rep(0,n_sect),
                  sect_hunt = sect
                  )


#######################################
### Initial Values for MCMC
#######################################
initsFun <- function()list(
                          mprec = runif(1, 1.5, 1.7),
                          fprec = runif(1, 2.7, 4.2),
                        #   m_period = seq(-2, 2, length = n_period),
                        #   f_period = seq(-2, 2, length = n_period),
                          m_age = seq(-6, -4, length = n_agem),
                          f_age = seq(-7, -5, length = n_agef)
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
              "m_age_mu",
              "f_age_mu",
              "mprec",
              "fprec",
              "f_age",
              "m_age"
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
                  niter = 10000,
                  nburnin = 5000,
                  nchains = 3,
                  inits = initsFun,
                  samplesAsCodaMCMC = TRUE,
                  summary = TRUE
                  )

runtime <- difftime(Sys.time(),
                    starttime,
                    units = "min")
mcmcout_survonly <- mcmcout
save(mcmcout_survonly,file="results/mcmcout_survonly.Rdata")

out <- mcmcout$samples
fit_sum <- mcmcout$summary$all.chains

gd <- gelman.diag(out[,c(grep("beta",rownames(fit_sum)),
                        grep("prec",rownames(fit_sum)),
                        grep("tau",rownames(fit_sum)))],
                  multivariate = FALSE)
gd
es <- effectiveSize(out[,c(grep("beta",rownames(fit_sum)),
                        grep("prec",rownames(fit_sum)),
                           grep("tau",rownames(fit_sum)))])
es
pdf("figures/traceplots_survonly.pdf")
traceplot(out[,"mprec"],ylab="mprec")
traceplot(out[,"fprec"],ylab="fprec")
for(i in 1:n_age){
  traceplot(out[,paste0("f_age[",i,"]")],ylab=paste0("f_age[",i,"]"))
}
for(i in 1:n_age){
  traceplot(out[,paste0("m_age[",i,"]")],ylab=paste0("m_age[",i,"]"))
}
dev.off()


Agegroups=1:n_agef
haz.mean = fit_sum[grep("f_age",rownames(fit_sum)),1][1:n_agef]
haz.lower = fit_sum[grep("f_age",rownames(fit_sum)),4][1:n_agef]
haz.upper = fit_sum[grep("f_age",rownames(fit_sum)),5][1:n_agef]

temp <- data.frame(Agegroups,haz.mean,haz.lower,haz.upper)
temp <- rbind(temp,temp[nrow(temp),])
temp$Agegroups[8] <- 8

fage.plot <- ggplot(data=temp)+geom_step(aes(x=Agegroups,y=haz.mean),direction="hv")+
  geom_step(aes(x=Agegroups,y=haz.lower),direction="hv",linetype=3)+
  geom_step(aes(x=Agegroups,y=haz.upper),direction="hv",linetype=3)+
  theme_bw()+ylim(-31,-2.5)+ggtitle("Female")+
  scale_x_continuous(breaks=(1:n_agef)+.5,labels=c("Fawn","1.5","2.5","3.5","4.5-5.5","6.5-8.5","9.5+"))
fage.plot


Agegroups=1:n_agem
haz.mean = fit_sum[grep("m_age",rownames(fit_sum)),1][1:n_agem]
haz.lower = fit_sum[grep("m_age",rownames(fit_sum)),4][1:n_agem]
haz.upper = fit_sum[grep("m_age",rownames(fit_sum)),5][1:n_agem]

temp <- data.frame(Agegroups,haz.mean,haz.lower,haz.upper)
temp <- rbind(temp,temp[nrow(temp),])
temp$Agegroups[n_agem+1] <- n_agem+1

mage.plot=ggplot(data=temp)+geom_step(aes(x=Agegroups,y=haz.mean),direction="hv")+
  geom_step(aes(x=Agegroups,y=haz.lower),direction="hv",linetype=3)+
  geom_step(aes(x=Agegroups,y=haz.upper),direction="hv",linetype=3)+
  scale_x_continuous(breaks = (1:n_agem)+.5,labels =c("Fawn","1.5","2.5","3.5","4.5-5.5","6.5+"))+
  theme_bw()+
  ylim(-31,-2.5)+ggtitle("Male")
  
mage.plot

combo_survonly <- grid.arrange(mage.plot,fage.plot,ncol=2)

ggsave("figures/combo_survonly.png",combo_survonly,height=4,width=8)

mcmcout_surv <- mcmcout
save(mcmcout_surv, file = "mcmcout_surv.Rdata")


########################################################################
########################################################################
########################################################################


##################################################
### Model code
###################################################

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

    #Likelihood for Collared individuals that were test negative at capture
    for (i in 1:n_fit_col) {
      teststatus_col[i] ~ dFOIcollar(left = left_age_col[i],
                  right = right_age_col[i],
                  sexfoi = sexfoi_col[i],
                  f_age = f_age[1:n_agef],
                  m_age = m_age[1:n_agem],
                  age_lookup = age_lookup_col[1:n_age_lookup_col],
                  space = space[sect_col[i]]
                  )
    }

    #Likelihood for individuals infected at capture
    for (i in 1:n_fit_inf) {
      teststatus_inf[i] ~ dFOIinf(left = left_age_inf[i],
                  sexfoi = sexfoi_inf[i],
                  f_age = f_age[1:n_agef],
                  m_age = m_age[1:n_agem],
                  age_lookup = age_lookup_col_inf[1:n_age_lookup_col_inf],
                  space = space[sect_inf_col[i]]
                  )
    }

    #Likelihood for hunter harvest surveillence data
    for (i in 1:n_fit_hunt) {
      teststatus_hunt[i] ~ dFOIhunt(agew = ageweeks_hunt[i],
                  sexfoi = sexfoi_hunt[i],
                  f_age = f_age[1:n_agef],
                  m_age = m_age[1:n_agem],
                  age_lookup = age_lookup[1:n_age_lookup],
                  space = space[sect_hunt[i]]
                  )
    }

})#end model statement

#######################################
### Data for Model Fitting
#######################################

# cwd_df$kill_year[cwd_df$kill_year==2022] <- 2021
# dim(cwd_df)
# cwd_df <- cwd_df[cwd_df$kill_year>2016,]
# n_fit_hunt <- nrow(cwd_df)

nimData <- list(
                teststatus_col = d_fit_col$censor,
                left_age_col = d_fit_col$left_age,
                right_age_col = d_fit_col$right_age,
                sexfoi_col = d_fit_col$sex,
                teststatus_inf = d_fit_inf$censor,
                left_age_inf = d_fit_inf$left_age,
                sexfoi_inf = d_fit_inf$sex,
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

nimConsts <- list(n_agem = n_agem,
                  n_agef = n_agef,
                  n_period = n_period,
                  n_age = n_age,
                  n_fit_col = n_fit_col,
                  n_fit_inf = n_fit_inf,
                  n_age_lookup_col = n_age_lookup_col,
                  age_lookup_col = age_lookup_col,
                  n_age_lookup_col_inf = n_age_lookup_col_inf,
                  age_lookup_col_inf = age_lookup_col_inf,
                  n_sect = n_sect,
                  sect_col = d_fit_col$sect,
                  sect_inf_col = d_fit_inf$sect,
                  n_age_lookup = n_age_lookup,
                  age_lookup = age_lookup,
                  period_lookup_hunt = period_lookup_hunt,
                  n_period_lookup_hunt = n_period_lookup_hunt,
                  birthweek_hunt = cwd_df$birthweek,
                  space = rep(0,n_sect),
                  sect_hunt = sect
                  )


#######################################
### Initial Values for MCMC
#######################################
initsFun <- function()list(
                          mprec = runif(1, 1.5, 1.7),
                          fprec = runif(1, 2.7, 4.2),
                        #   m_period = seq(-2, 2, length = n_period),
                        #   f_period = seq(-2, 2, length = n_period),
                          m_age = seq(-6, -4, length = n_agem),
                          f_age = seq(-7, -5, length = n_agef)
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
              "m_age_mu",
              "f_age_mu",
              "mprec",
              "fprec",
              "f_age",
              "m_age"
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
                  niter = 10000,
                  nburnin = 5000,
                  nchains = 3,
                  inits = initsFun,
                  samplesAsCodaMCMC = TRUE,
                  summary = TRUE
                  )

runtime <- difftime(Sys.time(),
                    starttime,
                    units = "min")

out <- mcmcout$samples
fit_sum <- mcmcout$summary$all.chains

gd <- gelman.diag(out[,c(grep("beta",rownames(fit_sum)),
                        grep("prec",rownames(fit_sum)),
                        grep("tau",rownames(fit_sum)))],
                  multivariate = FALSE)
gd
es <- effectiveSize(out[,c(grep("beta",rownames(fit_sum)),
                        grep("prec",rownames(fit_sum)),
                           grep("tau",rownames(fit_sum)))])
es
pdf("figures/traceplots_surv_inf.pdf")
traceplot(out[,"mprec"],ylab="mprec")
traceplot(out[,"fprec"],ylab="fprec")
for(i in 1:n_age){
  traceplot(out[,paste0("f_age[",i,"]")],ylab=paste0("f_age[",i,"]"))
}
for(i in 1:n_age){
  traceplot(out[,paste0("m_age[",i,"]")],ylab=paste0("m_age[",i,"]"))
}
dev.off()


Agegroups=1:n_agef
haz.mean = fit_sum[grep("f_age",rownames(fit_sum)),1][1:n_agef]
haz.lower = fit_sum[grep("f_age",rownames(fit_sum)),4][1:n_agef]
haz.upper = fit_sum[grep("f_age",rownames(fit_sum)),5][1:n_agef]

temp = data.frame(Agegroups,haz.mean,haz.lower,haz.upper)
temp <- rbind(temp,temp[nrow(temp),])
temp$Agegroups[8] <- 8

fage.plot=ggplot(data=temp)+geom_step(aes(x=Agegroups,y=haz.mean),direction="hv")+
  geom_step(aes(x=Agegroups,y=haz.lower),direction="hv",linetype=3)+
  geom_step(aes(x=Agegroups,y=haz.upper),direction="hv",linetype=3)+
  theme_bw()+ylim(-31,-2.5)+ggtitle("Female")+
  scale_x_continuous(breaks=(1:n_agef)+.5,labels=c("Fawn","1.5","2.5","3.5","4.5-5.5","6.5-8.5","9.5+"))
fage.plot


Agegroups=1:n_agem
haz.mean = fit_sum[grep("m_age",rownames(fit_sum)),1][1:n_agem]
haz.lower = fit_sum[grep("m_age",rownames(fit_sum)),4][1:n_agem]
haz.upper = fit_sum[grep("m_age",rownames(fit_sum)),5][1:n_agem]

temp = data.frame(Agegroups,haz.mean,haz.lower,haz.upper)
temp <- rbind(temp,temp[nrow(temp),])
temp$Agegroups[n_agem+1] <- n_agem+1

mage.plot=ggplot(data=temp)+geom_step(aes(x=Agegroups,y=haz.mean),direction="hv")+
  geom_step(aes(x=Agegroups,y=haz.lower),direction="hv",linetype=3)+
  geom_step(aes(x=Agegroups,y=haz.upper),direction="hv",linetype=3)+
  theme_bw()+
  ylim(-31,-2.5)+ggtitle("Male")+
  scale_x_continuous(breaks=(1:n_agem)+.5,labels=c("Fawn","1.5","2.5","3.5","4.5-5.5","6.5+"))
mage.plot

combo_surv_inf <- grid.arrange(mage.plot,fage.plot,ncol=2)

ggsave("figures/combo_surv_inf.png",combo_surv_inf,height=4,width=8)

mcmcout_surv_inf <- mcmcout
save(mcmcout_surv_inf, file = "mcmcout_surv_inf.Rdata")




######################################################
######################################################
######################################################
###
### Aggregating numeric table
###
######################################################
######################################################
######################################################


load("mcmcout_inf.Rdata")
load("mcmcout_surv.Rdata")
load("mcmcout_surv_inf.Rdata")


fit_sum_inf <- mcmcout_inf$summary$all.chains
fit_sum_surv <- mcmcout_surv$summary$all.chains
fit_sum_surv_inf <- mcmcout_surv_inf$summary$all.chains


Agegroups=rep(1:n_agef,3)
haz.mean <- c(fit_sum_inf[grep("f_age",rownames(fit_sum_inf)),1][1:n_agef],
              fit_sum_surv[grep("f_age",rownames(fit_sum_surv)),1][1:n_agef],
              fit_sum_surv_inf[grep("f_age",rownames(fit_sum_surv_inf)),1][1:n_agef])

haz.lower <- c(fit_sum_inf[grep("f_age",rownames(fit_sum_inf)),4][1:n_agef],
               fit_sum_surv[grep("f_age",rownames(fit_sum_surv)),4][1:n_agef],
               fit_sum_surv_inf[grep("f_age",rownames(fit_sum_surv_inf)),4][1:n_agef])

haz.upper <- c(fit_sum_surv_inf[grep("f_age",rownames(fit_sum_surv_inf)),5][1:n_agef],
               fit_sum_surv_inf[grep("f_age",rownames(fit_sum_surv_inf)),5][1:n_agef],
               fit_sum_surv_inf[grep("f_age",rownames(fit_sum_surv_inf)),5][1:n_agef])
sex=rep("Female",3*n_agef)
model=c(rep("Collar Only",n_agef),rep("Surveillance Only",n_agef),rep("Collar & Surveillance",n_agef))
tempf <- data.frame(model,sex,Agegroups,haz.mean,haz.lower,haz.upper)

Agegroups=rep(1:n_agem,3)
haz.mean <- c(fit_sum_inf[grep("m_age",rownames(fit_sum_inf)),1][1:n_agem],
              fit_sum_surv[grep("m_age",rownames(fit_sum_surv)),1][1:n_agem],
              fit_sum_surv_inf[grep("m_age",rownames(fit_sum_surv_inf)),1][1:n_agem])

haz.lower <- c(fit_sum_inf[grep("m_age",rownames(fit_sum_inf)),4][1:n_agem],
               fit_sum_surv[grep("m_age",rownames(fit_sum_surv)),4][1:n_agem],
               fit_sum_surv_inf[grep("m_age",rownames(fit_sum_surv_inf)),4][1:n_agem])

haz.upper <- c(fit_sum_surv_inf[grep("m_age",rownames(fit_sum_surv_inf)),5][1:n_agem],
               fit_sum_surv_inf[grep("m_age",rownames(fit_sum_surv_inf)),5][1:n_agem],
               fit_sum_surv_inf[grep("m_age",rownames(fit_sum_surv_inf)),5][1:n_agem])

sex=rep("Male",3*n_agem)
model=c(rep("Collar Only",n_agem),
        rep("Surveillance Only",n_agem),
        rep("Collar & Surveillance",n_agem))

tempm <- data.frame(model,sex,Agegroups,haz.mean,haz.lower,haz.upper)

output <- rbind(tempf,tempm)
write.csv(output,file="output_sample_sizes_foi.csv",row.names=FALSE)
