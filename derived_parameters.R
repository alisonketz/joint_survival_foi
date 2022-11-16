##########################
  ### Derived parameters
  ##########################

  # #Males
  # for (t in 1:nT_age) {
  #   llambda_age[t, 1] <- beta0 + age_effect[t]
  #   UCH0_age[t, 1] <- exp(llambda_age[t, 1])
  #   S0_age[t, 1] <- exp(-sum(UCH0_age[1:t, 1]))
  # }
  # for (t in 1:nT_period) {
  #   llambda_period[t, 1] <- beta0 + period_effect[t]
  #   UCH0_period[t, 1] <- exp(llambda_period[t, 1])
  #   S0_period[t, 1] <- exp(-sum(UCH0_period[1:t, 1]))
  # }

  # #Females
  # for (t in 1:nT_age) {
  #   llambda_age[t, 2] <- beta0 + age_effect[t] + beta_sex
  #   UCH0_age[t, 2] <- exp(llambda_age[t, 2])
  #   S0_age[t, 2] <- exp(-sum(UCH0_age[1:t, 2]))
  # }
  # for (t in 1:nT_period) {
  #   llambda_period[t, 2] <- beta0 + period_effect[t] + beta_sex
  #   UCH0_period[t, 2] <- exp(llambda_period[t, 2])
  #   S0_period[t, 2] <- exp(-sum(UCH0_period[1:t, 2]))
  # }

  ###
  ### Full hazard surface
  ### calculating the cause-specific hazards based on the constant hazards model
  ###

  # for(i in 1:nT_age){
  #   for(j in 1:nT_period){
  #     llambda_all[i,j,1] <- beta0 + age_effect[i] + period_effect[j]
  #     llambda_all[i,j,2] <- beta0 + age_effect[i] + period_effect[j] + beta_sex
  #     UCH0_all[i,j,1] <- exp(llambda_all[i,j,1])
  #     UCH0_all[i,j,2] <- exp(llambda_all[i,j,2])
  #     S0_all[i,j,1] <- exp(-sum(UCH0_all[1:i,1:j,1]))
  #     S0_all[i,j,2] <- exp(-sum(UCH0_all[1:i,1:j,2]))
  #     for(k in 1:n_causes){
  #       cause_haz[i,j,1] <- exp(llambda_all[i,j,1])*p_cause[k]
  #       cause_haz[i,j,2] <- exp(llambda_all[i,j,2])*p_cause[k]
  #     }

  #   }
  # }