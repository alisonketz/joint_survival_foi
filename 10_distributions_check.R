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

####################################################################
###
### Data and paramters needed to calculate likelihoods
###
###
####################################################################

beta0_sus <- -5
beta0_inf <- -3.2
# floor(as.duration(ymd("1992-09-01") %--% ymd("2021-09-01"))/dweeks(1))
# floor(as.duration(ymd("1992-09-01") %--% ymd("2022-05-31"))/dweeks(1))
# floor(as.duration(ymd("1992-05-15") %--% ymd("1992-09-01"))/dweeks(1))
# floor(as.duration(ymd("1992-09-01") %--% ymd("1993-01-01"))/dweeks(1))
# floor(as.duration(ymd("1993-01-01") %--% ymd("1993-09-01"))/dweeks(1))
# floor(as.duration(ymd("1993-09-01") %--% ymd("1994-01-01"))/dweeks(1))
# floor(as.duration(ymd("1993-01-01") %--% ymd("2002-01-01")))
# floor(as.duration(ymd("2002-03-01") %--% ymd("2002-09-01"))/dweeks(1)) - 1
# floor(as.duration(ymd("2003-01-01") %--% ymd("2002-09-01"))/dweeks(1)) - 1
# floor(as.duration(ymd("2003-01-01") %--% ymd("2003-09-01"))/dweeks(1)) - 1
# floor(as.duration(ymd("2004-01-01") %--% ymd("2003-09-01"))/dweeks(1)) - 1
# floor(as.duration(ymd("2022-01-01") %--% ymd("2022-05-31"))/dweeks(1)) - 1

period_effect_survival <- c(rep(-6,15),#may1992-sep1992
                            rep(c(rep(-4.5,18),rep(-6,34)),29),#sep1992 - sep2021
                            rep(-4.5,19),#sep2021-jan2022
                            rep(-6,20))#jan2022-May15,2022

age_effect_survival  <- exp(-.01*seq(1:962))
age_effect_survival <- age_effect_survival - mean(age_effect_survival)
beta_sex <- -.5

#load FOI parameters based on run w/o survival integrated
load("~/Documents/Transmission/Transmission_v3/27_urw2_bym2_timeRW1_2021/fit_sum.Rdata")
fit_sum_foi <- fit_sum
f_age_foi <- fit_sum_foi[grep("f_age",rownames(fit_sum_foi)),][1:7,1]
m_age_foi <- fit_sum_foi[grep("m_age",rownames(fit_sum_foi)),][1:6,1]
f_period_foi <- c(fit_sum_foi[grep("f_time",rownames(fit_sum_foi)),1])
m_period_foi <- c(fit_sum_foi[grep("m_time",rownames(fit_sum_foi)),1])
space_mn <- rep(0,n_sect) #leave 0 for now
# space_mn  <- c(fit_sum_foi[grep("space",rownames(fit_sum_foi)),1])[sect_num]



