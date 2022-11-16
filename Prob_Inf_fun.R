
########################################
### Function to calculate the 
### probability of infection in a given 
### interval for a single individual
########################################


prob_inf <- function(

        nT_age,
        left_age,
        nT_period,
        left_period,
        age_effect,
        period_effect

){

  maxtimes <- nT_period - left_period

  ########################################
  ### Log hazard
  ########################################
  hazard <- rep(NA, nT_age)
  hazard <- c(rep(0,left_age - 1),
              beta0 + age_effect[left_age:(left_age + maxtimes - 1)] +
                      period_effect[left_period:(left_period + maxtimes - 1)],
                      rep(0, nT_age - (left_age - 1 + maxtimes)))
 ########################################
  ### Probability of mortality
  ########################################

  test_stat <- rep(0, nT_age)
  for (j in left_age:(left_age + maxtimes - 1)) {
      if (j==left_age) {
        test_stat[j] <- (1 - exp(-exp(hazard[j])))
      } else {
        test_stat[j] <-
          (1 - exp(-exp(hazard[j])))*exp(-sum(exp(hazard[left_age:(j - 1)])))
      }
    
      test_stat[left_age + maxtimes] <-
          exp(-sum(exp(hazard[left_age:(left_age + maxtimes - 1)])))
  }


}