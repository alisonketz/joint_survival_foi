##########################################################################
###
### week time blocks for harvest
###
###
#########################################################################

# #year 1, 2017
# non_harvest_survival_end_yr1 <- week(ymd('2017-09-16')) - start_week
# bow_start_yr1 <- week(ymd('2017-09-16')) - start_week + 1
# bow_end_yr1 <- week(ymd('2018-01-07'))+ 52 - start_week + 1
# gun_start_yr1 <- week(ymd('2017-11-18'))- start_week + 1
# gun_end_yr1 <- week(ymd('2017-11-26'))- start_week + 1
# gun_holiday_start_yr1 <- week(ymd('2017-12-24'))- start_week + 1
# gun_holiday_end_yr1 <- week(ymd('2018-01-01')) + 52 - start_week + 1
# study_end_yr1 <- week(ymd('2018-01-07')) + 52 - start_week + 1

# #year 2,2018
# non_harvest_survival_end_yr2 <- week(ymd('2018-09-15'))-start_week
# bow_start_yr2 <- week(ymd('2018-09-15'))-start_week+1
# bow_end_yr2 <- week(ymd('2019-01-06'))-start_week+1+nT_yr1
# gun_start_yr2 <- week(ymd('2018-11-17'))-start_week+1
# gun_end_yr2 <- week(ymd('2018-11-25'))-start_week+1
# gun_holiday_start_yr2 <- week(ymd('2017-12-24'))-start_week+1
# gun_holiday_end_yr2 <- week(ymd('2019-01-01'))-start_week+1+nT_yr1
# study_end_yr2 <- bow_end_yr2

# #year 3, 2019
# non_harvest_survival_end_yr3 <- week(ymd('2019-09-14')) - start_week#+nT_yr1+nT_yr2
# # Bow & Crossbow: Sept. 14 – Jan_ 5, 2020
# bow_start_yr3 <- week(ymd('2019-09-14')) - start_week + 1#+nT_yr1+nT_yr2
# bow_end_yr3 <- week(ymd('2020-01-05')) - start_week + 1 + nT_yr2
# # 9-day Gun Hunt: Nov. 23 – Dec. 1
# gun_start_yr3 <- week(ymd('2019-11-21')) - start_week + 1#+nT_yr1+nT_yr2
# gun_end_yr3 <- week(ymd('2019-12-01')) - start_week + 1#+nT_yr1+nT_yr2
# # December Holiday Hunt: Dec. 24 – Jan_ 1, 2020
# gun_holiday_start_yr3 <- week(ymd('2019-12-24')) - start_week + 1#+nT_yr1+nT_yr2
# gun_holiday_end_yr3 <- week(ymd('2020-01-01')) - start_week + 1 + nT_yr2
# study_end_yr3 <- bow_end_yr3

# #year 4, 2020
# non_harvest_survival_end_yr4 <- week(ymd('2020-09-12')) - start_week#+nT_yr1+nT_yr2+nT_yr3
# # Bow & Crossbow: Sept. 12 – Jan_ 3, 2021
# bow_start_yr4 <- week(ymd('2020-09-12'))-start_week+1#+nT_yr3
# bow_end_yr4 <- week(ymd('2021-01-03'))-start_week+1+nT_yr3
# # 9-day Gun Hunt: Nov. 21 – 29
# gun_start_yr4 <- week(ymd('2020-11-21'))-start_week+1#+nT_yr1+nT_yr2+nT_yr3
# gun_end_yr4 <- week(ymd('2020-12-01'))-start_week+1#+nT_yr1+nT_yr2+nT_yr3
# # December Holiday Hunt: Dec. 24 – Jan_ 1, 2020
# gun_holiday_start_yr4 <- week(ymd('2020-12-24'))-start_week+1#+nT_yr1+nT_yr2+nT_yr3
# gun_holiday_end_yr4 <- week(ymd('2021-01-01'))-start_week+1+nT_yr3
# study_end_yr4 <- bow_end_yr4


# #year 5, 2020
# non_harvest_survival_end_yr5 <- week(ymd('2021-09-12')) - start_week#+nT_yr1+nT_yr2+nT_yr3
# # Bow & Crossbow: Sept. 12 – Jan_ 3, 2021
# bow_start_yr5 <- week(ymd('2021-09-12'))-start_week+1#+nT_yr3
# bow_end_yr5 <- week(ymd('2022-01-03'))-start_week+1+nT_yr4
# # 9-day Gun Hunt: Nov. 21 – 29
# gun_start_yr5 <- week(ymd('2021-11-20'))-start_week+1#+nT_yr1+nT_yr2+nT_yr3
# gun_end_yr5 <- week(ymd('2021-11-28'))-start_week+1#+nT_yr1+nT_yr2+nT_yr3
# # December Holiday Hunt: Dec. 24 – Jan_ 1, 2020
# gun_holiday_start_yr5 <- week(ymd('2021-12-24')) - start_week + 1#+nT_yr1+nT_yr2+nT_yr3
# gun_holiday_end_yr5 <- week(ymd('2022-01-01')) - start_week + 1 + nT_yr4
# study_end_yr5 <- bow_end_yr5


# # likely bow harvest
# # which(d_mort$weapon!="Rifle" & d_mort$cause1 == "Hunter harvest")
# d_mort$bow <- 0
# d_mort$bow[which(d_mort$weapon!="Rifle" & d_mort$cause1 == "Hunter harvest")[-1]] <- 1

# #setting harvest_season for derived parameters
# bow_harvest_hazard_yr1 <- c(rep(0,non_harvest_survival_end_yr1),rep(1,bow_end_yr1-bow_start_yr1+1))
# bow_harvest_hazard_yr2 <- c(rep(0,non_harvest_survival_end_yr2),rep(1,bow_end_yr2-bow_start_yr2+1))
# bow_harvest_hazard_yr3 <- c(rep(0,non_harvest_survival_end_yr3),rep(1,bow_end_yr3-bow_start_yr3+1))
# bow_harvest_hazard_yr4 <- c(rep(0,non_harvest_survival_end_yr4),rep(1,bow_end_yr4-bow_start_yr4+1))
# bow_harvest_hazard_yr5 <- c(rep(0,non_harvest_survival_end_yr4),rep(1,bow_end_yr4-bow_start_yr5+1))

# bow_harvest_haz <- c(bow_harvest_hazard_yr1,bow_harvest_hazard_yr2,bow_harvest_hazard_yr3,bow_harvest_hazard_yr4)

# gun_harvest_haz_yr1 <- c(rep(0,gun_start_yr1-1),rep(1,gun_end_yr1-gun_start_yr1+1),rep(0,bow_end_yr1-gun_end_yr1))
# gun_harvest_haz_yr2 <- c(rep(0,gun_start_yr2-1),rep(1,gun_end_yr2-gun_start_yr2+1),rep(0,bow_end_yr2-gun_end_yr2))
# gun_harvest_haz_yr3 <- c(rep(0,gun_start_yr3-1),rep(1,gun_end_yr3-gun_start_yr3+1),rep(0,bow_end_yr3-gun_end_yr3))
# gun_harvest_haz_yr4 <- c(rep(0,gun_start_yr4-1),rep(1,gun_end_yr4-gun_start_yr4+1),rep(0,bow_end_yr4-gun_end_yr4))

# gun_harvest_haz <- c(gun_harvest_haz_yr1,gun_harvest_haz_yr2,gun_harvest_haz_yr3,gun_harvest_haz_yr4)







#notes on cleaning Tooth data
###
### For all individuals where an interval of ages is returned from the
### cementum annuli data, we are going to use the minimum age from that interval
### i.e. 1-2 == 1, or 10-11 == 10

#this individual was called a big 1.5 yr old buck, so I will denote him as 1 year old
# d_tooth[d_tooth$age == "0-1",]
# d_cap[d_cap$lowtag==d_tooth$lowtag[d_tooth$age == "0-1"],]

#there were 4 with 1-2 year old designations
#three of these were 20mo and 1 was >2yrs,
#these will all be treated as 1 year olds
# d_tooth[d_tooth$age=="1-2",1]
# d_cap[d_cap$lowtag %in% d_tooth$lowtag[d_tooth$age == "1-2"],c(27,57)]
# names(d_cap)
# d_tooth$age[d_tooth$age == "1-2"] <- "1.5"
#there's 5 that were 2-3 yrs old, one was 20mth age class, but since tooth pulled
#should denote it at 2 yrs old
# d_tooth[d_tooth$age=="2-3",]
# d_cap[d_cap$lowtag %in% d_tooth$lowtag[d_tooth$age == "2-3"],c(27,57)]
# d_tooth$age[d_tooth$age == "2-3"] <- "2.5"


#there are 7 for which we have no exact ages
# check<- d_cap[d_cap$lowtag %in% d_tooth$lowtag[d_tooth$age=="NOAGE"],]
# check$observationsnotes

# check[,c(4,5,23,27,29,57)]
# check[1,c(4,5,23,27,29,57)]#20mo in ageclass
# check[2,c(4,5,23,27,29,57)]#tooth removal broke tooth,
# check[3,c(4,5,23,27,29,57)]#tooth removed, broke tooth
# check[4,c(4,5,23,27,29,57)]#tooth removed, broke tooth
# check[5,c(4,5,23,27,29,57)]#5-6 years old in comments
# check[6,c(4,5,23,27,29,57)]#tooth removed, broke tooth
# check[7,c(4,5,23,27,29,57)]#tooth removed, broke tooth


