####################################################################
############ Code to reproduce the main analysis reported in #######
############ Ranzani, Hitchings, et al, BMJ 2021 10.1136/bmj.n2015
############

library(tidyverse)
library(tidylog)
library(survival)
library(rms)
library(broom)
library(broom.helpers)
library(broom.mixed)


#### loading data ####
case_control_data_primary <- read_csv("PATH/CoronaVac_BMJ_CaseControlPairs.csv")

#### adjusting factors / exposure #####
case_control_data_primary <-
  case_control_data_primary %>% 
  mutate(
    vcat_combined = factor(vcat_combined,
                                  levels = c("Unvaccinated", "Dose 1 0-13 days",
                                             "Dose 1 >=14 days", "Dose 2 0-13 days",
                                             "Dose 2 >=14 days")),
    
    vcat_combined_bias = factor(vcat_combined_bias,
                                levels = c("Unvaccinated", "Dose 1 0-6 days","Dose 1 7-13 days",
                                           "Dose 1 >=14 days", "Dose 2 0-13 days",
                                           "Dose 2 >=14 days"))
  )

#### Main Analysis, Symptomatic ####


model_crude_tbl <-
  clogit(case ~ vcat_combined + strata(stratum_no), data = case_control_data_primary) %>% 
  broom.helpers::tidy_and_attach(exponentiate = TRUE, conf.int = TRUE)


model_adj_tbl <-
  clogit(case ~ vcat_combined + age + n_comorb_cat + strata(stratum_no), data = case_control_data_primary) %>% 
  broom.helpers::tidy_and_attach(exponentiate = TRUE,conf.int = TRUE)


models <- bind_rows(model_crude_tbl %>% mutate(model = "Crude"), 
                    model_adj_tbl   %>% mutate(model = "Adj"))

models <-
  models %>% 
  filter(str_detect(term, "Dose")) %>% 
  mutate(OR = paste(
    round(estimate,3), " (",
    round(conf.low,3), "-",
    round(conf.high,3), ")", sep = ""),
    p_value = as.character(round(p.value,3)),
    ve = 
        paste(
          round(1-estimate,3)*100, "% (",
          round(1-conf.high,3)*100, "-", 
          round(1-conf.low,3)*100, ")",
          sep = ""))

#### Main Analysis, Hospitalization / Death ####

model_crude_hosp <-
  clogit(case ~ vcat_combined + strata(stratum_no), data = case_control_data_primary %>% filter(hosp_pairs==1)) %>% 
  broom.helpers::tidy_and_attach(exponentiate = TRUE,conf.int = TRUE)

model_adj_hosp <-
  clogit(case ~ vcat_combined + rcs(age,3) + n_comorb_cat + strata(stratum_no), data = case_control_data_primary %>% filter(hosp_pairs==1)) %>% 
  broom.helpers::tidy_and_attach(exponentiate = TRUE,conf.int = TRUE)



model_crude_death <-
  clogit(case ~ vcat_combined + strata(stratum_no), data = case_control_data_primary %>% filter(death_pairs==1)) %>% 
  broom.helpers::tidy_and_attach(exponentiate = TRUE,conf.int = TRUE)

model_adj_death <-
  clogit(case ~ vcat_combined + rcs(age,5) + n_comorb_cat + strata(stratum_no), data = case_control_data_primary %>% filter(death_pairs==1)) %>% 
  broom.helpers::tidy_and_attach(exponentiate = TRUE,conf.int = TRUE)



models_severity <- bind_rows(
  
  model_crude_hosp   %>% mutate(model = "Crude", outcome = "Hospitalization"), 
  model_adj_hosp     %>% mutate(model = "Adj"  , outcome = "Hospitalization"),
  
  model_crude_death  %>% mutate(model = "Crude", outcome = "Death"),
  model_adj_death    %>% mutate(model = "Adj"  , outcome = "Death"),
  
)

models_severity <-
  models_severity %>% 
  filter(str_detect(term, "Dose")) %>% 
  mutate(OR = paste(
    round(estimate,3), " (",
    round(conf.low,3), "-",
    round(conf.high,3), ")", sep = ""),
    ve = 
        paste(
          round(1-estimate,3)*100, "% (",
          round(1-conf.high,3)*100, "-", 
          round(1-conf.low,3)*100, ")",
          sep = ""),
    p_value = as.character(round(p.value,3)),
  )

#### Sensitivity analysis, Symptomatic ####


model_adj_time_tbl <-
  clogit(case ~ vcat_combined + age + n_comorb_cat + strata(stratum_no) + day_of_year_anonymized, data = case_control_data_primary) %>%
  broom.helpers::tidy_and_attach(exponentiate = TRUE,conf.int = TRUE)


model_adj_time <-
  model_adj_time_tbl %>% 
  filter(str_detect(term, "Dose")) %>% 
  mutate(OR = paste(
    round(estimate,3), " (",
    round(conf.low,3), "-",
    round(conf.high,3), ")", sep = ""),
    p_value = as.character(round(p.value,3)),
    ve = 
        paste(
          round(1-estimate,3)*100, "% (",
          round(1-conf.high,3)*100, "-", 
          round(1-conf.low,3)*100, ")",
          sep = ""))


#### Sensitivity analysis, Symptomatic, bias indicator 0-6 ####


model_adj_tbl_bias <-
  clogit(case ~ vcat_combined_bias + age + n_comorb_cat + strata(stratum_no), data = case_control_data_primary) %>% 
  broom.helpers::tidy_and_attach(exponentiate = TRUE,conf.int = TRUE)

model_adj_hosp_bias <-
  clogit(case ~ vcat_combined_bias + rcs(age,3) + n_comorb_cat + strata(stratum_no), data = case_control_data_primary %>% filter(hosp_pairs==1)) %>% 
  broom.helpers::tidy_and_attach(exponentiate = TRUE,conf.int = TRUE)


model_adj_death_bias <-
  clogit(case ~ vcat_combined_bias + rcs(age,5) + n_comorb_cat + strata(stratum_no), data = case_control_data_primary %>% filter(death_pairs==1)) %>% 
  broom.helpers::tidy_and_attach(exponentiate = TRUE,conf.int = TRUE)



models_bias <- 
  bind_rows(
    model_adj_tbl_bias   %>% mutate(model = "Adj", outcome = "Case"),
    model_adj_hosp_bias  %>% mutate(model = "Adj", outcome = "Hosp"),
    model_adj_death_bias %>% mutate(model = "Adj", outcome = "Death")
  )



models_bias <-
  models_bias %>% 
  filter(str_detect(term, "Dose")) %>% 
  mutate(OR = paste(
    round(estimate,3), " (",
    round(conf.low,3), "-",
    round(conf.high,3), ")", sep = ""),
    p_value = as.character(round(p.value,3)),
    ve = 
        paste(
          round(1-estimate,3)*100, "% (",
          round(1-conf.high,3)*100, "-", 
          round(1-conf.low,3)*100, ")",
          sep = ""))

#### Interactions with Age group, main analysis, Symptomatic, Hosp, Death ####

model_withdummy_agep  <- clogit(case ~ 
                                  factor(vcat_combined=="Dose 1 0-13 days") + 
                                  factor(vcat_combined=="Dose 1 >=14 days") + 
                                  factor(vcat_combined=="Dose 2 0-13 days") + 
                                  factor(vcat_combined=="Dose 2 >=14 days") + 
                                  age + factor(age_plot) + n_comorb_cat + strata(stratum_no), data = case_control_data_primary)

model_withdummy_agep_2dose  <- clogit(case ~ 
                                        factor(vcat_combined=="Dose 1 0-13 days") + 
                                        factor(vcat_combined=="Dose 1 >=14 days") + 
                                        factor(vcat_combined=="Dose 2 0-13 days") + 
                                        factor(vcat_combined=="Dose 2 >=14 days") + 
                                        relevel(factor(vcat_combined=="Dose 2 >=14 days"),ref="FALSE"):factor(age_plot) + 
                                        age + n_comorb_cat + strata(stratum_no), data = case_control_data_primary)

p_int_agep <- anova(model_withdummy_agep, model_withdummy_agep_2dose)


model_int_age_1  <- clogit(case ~ vcat_combined * relevel(factor(age_plot), ref = 1) + age +  n_comorb_cat + strata(stratum_no), data = case_control_data_primary)
model_int_age_2  <- clogit(case ~ vcat_combined * relevel(factor(age_plot), ref = 2) + age +  n_comorb_cat + strata(stratum_no), data = case_control_data_primary)
model_int_age_3  <- clogit(case ~ vcat_combined * relevel(factor(age_plot), ref = 3) + age +  n_comorb_cat + strata(stratum_no), data = case_control_data_primary)



# hosp by age
model_withdummy_agep_hosp  <- clogit(case ~ 
                                       factor(vcat_combined=="Dose 1 0-13 days") + 
                                       factor(vcat_combined=="Dose 1 >=14 days") + 
                                       factor(vcat_combined=="Dose 2 0-13 days") + 
                                       factor(vcat_combined=="Dose 2 >=14 days") + 
                                       rcs(age,3) + factor(age_plot) + n_comorb_cat + strata(stratum_no), 
                                     data = case_control_data_primary %>% filter(hosp_pairs == 1))

model_withdummy_agep_2dose_hosp <- clogit(case ~ 
                                            factor(vcat_combined=="Dose 1 0-13 days") + 
                                            factor(vcat_combined=="Dose 1 >=14 days") + 
                                            factor(vcat_combined=="Dose 2 0-13 days") + 
                                            factor(vcat_combined=="Dose 2 >=14 days") + 
                                            relevel(factor(vcat_combined=="Dose 2 >=14 days"),ref="FALSE"):factor(age_plot) + 
                                            rcs(age,3) + n_comorb_cat + strata(stratum_no), 
                                          data = case_control_data_primary %>% filter(hosp_pairs==1))

p_int_agep_hosp <- anova(model_withdummy_agep_hosp, model_withdummy_agep_2dose_hosp)


model_int_age_1_hosp  <- clogit(case ~ vcat_combined * relevel(factor(age_plot), ref = 1) + rcs(age,3) +  n_comorb_cat + strata(stratum_no), data = case_control_data_primary %>% filter(hosp_pairs==1))
model_int_age_2_hosp  <- clogit(case ~ vcat_combined * relevel(factor(age_plot), ref = 2) + rcs(age,3) +  n_comorb_cat + strata(stratum_no), data = case_control_data_primary %>% filter(hosp_pairs==1))
model_int_age_3_hosp  <- clogit(case ~ vcat_combined * relevel(factor(age_plot), ref = 3) + rcs(age,3) +  n_comorb_cat + strata(stratum_no), data = case_control_data_primary %>% filter(hosp_pairs==1))

# death

model_withdummy_agep_death  <- clogit(case ~ 
                                        factor(vcat_combined=="Dose 1 0-13 days") + 
                                        factor(vcat_combined=="Dose 1 >=14 days") + 
                                        factor(vcat_combined=="Dose 2 0-13 days") + 
                                        factor(vcat_combined=="Dose 2 >=14 days") + 
                                        rcs(age,5) + factor(age_plot) + n_comorb_cat + strata(stratum_no), 
                                      data = case_control_data_primary %>% filter(death_pairs == 1))

model_withdummy_agep_2dose_death <- clogit(case ~ 
                                             factor(vcat_combined=="Dose 1 0-13 days") + 
                                             factor(vcat_combined=="Dose 1 >=14 days") + 
                                             factor(vcat_combined=="Dose 2 0-13 days") + 
                                             factor(vcat_combined=="Dose 2 >=14 days") + 
                                             relevel(factor(vcat_combined=="Dose 2 >=14 days"),ref="FALSE"):factor(age_plot) + 
                                             rcs(age,5) + n_comorb_cat + strata(stratum_no), 
                                           data = case_control_data_primary %>% filter(death_pairs==1))

p_int_agep_death <- anova(model_withdummy_agep_death, model_withdummy_agep_2dose_death)


model_int_age_1_death  <- clogit(case ~ vcat_combined * relevel(factor(age_plot), ref = 1) + rcs(age,5) +  n_comorb_cat + strata(stratum_no), data = case_control_data_primary %>% filter(death_pairs==1))
model_int_age_2_death  <- clogit(case ~ vcat_combined * relevel(factor(age_plot), ref = 2) + rcs(age,5) +  n_comorb_cat + strata(stratum_no), data = case_control_data_primary %>% filter(death_pairs==1))
model_int_age_3_death  <- clogit(case ~ vcat_combined * relevel(factor(age_plot), ref = 3) + rcs(age,5) +  n_comorb_cat + strata(stratum_no), data = case_control_data_primary %>% filter(death_pairs==1))



## tidy
model_int_age_1 <- model_int_age_1 %>% broom.helpers::tidy_and_attach(exponentiate = TRUE,conf.int = TRUE)
model_int_age_2 <- model_int_age_2 %>% broom.helpers::tidy_and_attach(exponentiate = TRUE,conf.int = TRUE)
model_int_age_3 <- model_int_age_3 %>% broom.helpers::tidy_and_attach(exponentiate = TRUE,conf.int = TRUE)

model_int_age_1_hosp <- model_int_age_1_hosp %>% broom.helpers::tidy_and_attach(exponentiate = TRUE,conf.int = TRUE)
model_int_age_2_hosp <- model_int_age_2_hosp %>% broom.helpers::tidy_and_attach(exponentiate = TRUE,conf.int = TRUE)
model_int_age_3_hosp <- model_int_age_3_hosp %>% broom.helpers::tidy_and_attach(exponentiate = TRUE,conf.int = TRUE)

model_int_age_1_death <- model_int_age_1_death %>% broom.helpers::tidy_and_attach(exponentiate = TRUE,conf.int = TRUE)
model_int_age_2_death <- model_int_age_2_death %>% broom.helpers::tidy_and_attach(exponentiate = TRUE,conf.int = TRUE)
model_int_age_3_death <- model_int_age_3_death %>% broom.helpers::tidy_and_attach(exponentiate = TRUE,conf.int = TRUE)

# p int
p_age        <- p_int_agep$`P(>|Chi|)`
p_age_hosp   <- p_int_agep_hosp$`P(>|Chi|)`
p_age_death  <- p_int_agep_death$`P(>|Chi|)`



models_int_main <- 
  bind_rows(
    model_int_age_1 %>% mutate(model = "70-74", outcome = "Symptomatic", p_int = p_age[2]),
    model_int_age_2 %>% mutate(model = "75-79", outcome = "Symptomatic", p_int = NA_real_),
    model_int_age_3 %>% mutate(model = "80+"  , outcome = "Symptomatic", p_int = NA_real_),
    
    model_int_age_1_hosp %>% mutate(model = "70-74", outcome = "Hospitalized", p_int = p_age_hosp[2]),
    model_int_age_2_hosp %>% mutate(model = "75-79", outcome = "Hospitalized", p_int = NA_real_),
    model_int_age_3_hosp %>% mutate(model = "80+"  , outcome = "Hospitalized", p_int = NA_real_),
    
    model_int_age_1_death %>% mutate(model = "70-74", outcome = "Deaths", p_int = p_age_death[2]),
    model_int_age_2_death %>% mutate(model = "75-79", outcome = "Deaths", p_int = NA_real_),
    model_int_age_3_death %>% mutate(model = "80+"  , outcome = "Deaths", p_int = NA_real_),
    
  )

models_int_main <-
  models_int_main %>% 
  filter(str_detect(term, "Dose 2 >=14 days") & 
           !str_detect(term, "relevel"))  %>% 
  mutate(OR = paste(
    round(estimate,3), " (",
    round(conf.low,3), "-",
    round(conf.high,3), ")", sep = ""),
    ve = 
        paste(
          format(round(1-estimate,3) *100, nsmall=1, trim = TRUE, scientific = FALSE), "% (",
          format(round(1-conf.high,3)*100, nsmall=1, trim = TRUE, scientific = FALSE), "-", 
          format(round(1-conf.low,3) *100, nsmall=1, trim = TRUE, scientific = FALSE), ")",
          sep = ""),
    p_int = as.character(round(p_int,3)),
  )

#### Interactions other subgroups, main analysis, Symptomatic ####


# Sex
model_withdummy_sex  <- clogit(case ~ 
                                 factor(vcat_combined=="Dose 1 0-13 days") + 
                                 factor(vcat_combined=="Dose 1 >=14 days") + 
                                 factor(vcat_combined=="Dose 2 0-13 days") + 
                                 factor(vcat_combined=="Dose 2 >=14 days") + 
                                 age + n_comorb_cat + sex + strata(stratum_no), data = case_control_data_primary)

model_withdummy_sex_2dose  <- clogit(case ~ 
                                       factor(vcat_combined=="Dose 1 0-13 days") + 
                                       factor(vcat_combined=="Dose 1 >=14 days") + 
                                       factor(vcat_combined=="Dose 2 0-13 days") + 
                                       factor(vcat_combined=="Dose 2 >=14 days") + 
                                       relevel(factor(vcat_combined=="Dose 2 >=14 days"),ref="FALSE"):sex + 
                                       age + n_comorb_cat + strata(stratum_no), data = case_control_data_primary)

p_int_sex <- anova(model_withdummy_sex, model_withdummy_sex_2dose)


model_int_sex_1  <- clogit(case ~ vcat_combined * relevel(factor(sex), ref = 1) + age + n_comorb_cat + strata(stratum_no), data = case_control_data_primary)
model_int_sex_2  <- clogit(case ~ vcat_combined * relevel(factor(sex), ref = 2) + age + n_comorb_cat + strata(stratum_no), data = case_control_data_primary)

#Comorb


model_withdummy_comor  <- clogit(case ~ factor(vcat_combined=="Dose 1 0-13 days") + 
                                   factor(vcat_combined=="Dose 1 >=14 days") + 
                                   factor(vcat_combined=="Dose 2 0-13 days") + 
                                   factor(vcat_combined=="Dose 2 >=14 days") + 
                                   age  + anycomorb + strata(stratum_no), data = case_control_data_primary)

model_withdummy_comor_2dose  <- clogit(case ~ factor(vcat_combined=="Dose 1 0-13 days") + 
                                         factor(vcat_combined=="Dose 1 >=14 days") + 
                                         factor(vcat_combined=="Dose 2 0-13 days") + 
                                         factor(vcat_combined=="Dose 2 >=14 days") + 
                                         relevel(factor(vcat_combined=="Dose 2 >=14 days"),ref="FALSE"):anycomorb + 
                                         age  + strata(stratum_no), data = case_control_data_primary)

p_int_comor <- anova(model_withdummy_comor, model_withdummy_comor_2dose)


model_int_co_1   <- clogit(case ~ vcat_combined * relevel(factor(anycomorb), ref = 1) + age  + strata(stratum_no), data = case_control_data_primary)
model_int_co_2   <- clogit(case ~ vcat_combined * relevel(factor(anycomorb), ref = 2) + age  + strata(stratum_no), data = case_control_data_primary)

# Cardio

model_withdummy_cvd  <- clogit(case ~ 
                                 factor(vcat_combined=="Dose 1 0-13 days") + 
                                 factor(vcat_combined=="Dose 1 >=14 days") + 
                                 factor(vcat_combined=="Dose 2 0-13 days") + 
                                 factor(vcat_combined=="Dose 2 >=14 days") + 
                                 age + n_comorb_cat + CARDIOPATI_m + strata(stratum_no), data = case_control_data_primary)

model_withdummy_cvd_2dose  <- clogit(case ~ 
                                       factor(vcat_combined=="Dose 1 0-13 days") + 
                                       factor(vcat_combined=="Dose 1 >=14 days") + 
                                       factor(vcat_combined=="Dose 2 0-13 days") + 
                                       factor(vcat_combined=="Dose 2 >=14 days") + 
                                       relevel(factor(vcat_combined=="Dose 2 >=14 days"),ref="FALSE"):CARDIOPATI_m + 
                                       age + n_comorb_cat + strata(stratum_no), data = case_control_data_primary)

p_int_cvd <- anova(model_withdummy_cvd, model_withdummy_cvd_2dose)


model_int_cardio_1   <- clogit(case ~ vcat_combined * relevel(factor(CARDIOPATI_m), ref = 1) + age  + n_comorb_cat+strata(stratum_no), data = case_control_data_primary)
model_int_cardio_2   <- clogit(case ~ vcat_combined * relevel(factor(CARDIOPATI_m), ref = 2) + age  + n_comorb_cat+strata(stratum_no), data = case_control_data_primary)

# Diabetes

model_withdummy_dm  <- clogit(case ~ 
                                factor(vcat_combined=="Dose 1 0-13 days") + 
                                factor(vcat_combined=="Dose 1 >=14 days") + 
                                factor(vcat_combined=="Dose 2 0-13 days") + 
                                factor(vcat_combined=="Dose 2 >=14 days") + 
                                age + n_comorb_cat + DIABETES_m + strata(stratum_no), data = case_control_data_primary)

model_withdummy_dm_2dose  <- clogit(case ~ 
                                      factor(vcat_combined=="Dose 1 0-13 days") + 
                                      factor(vcat_combined=="Dose 1 >=14 days") + 
                                      factor(vcat_combined=="Dose 2 0-13 days") + 
                                      factor(vcat_combined=="Dose 2 >=14 days") + 
                                      relevel(factor(vcat_combined=="Dose 2 >=14 days"),ref="FALSE"):DIABETES_m + 
                                      age + n_comorb_cat + strata(stratum_no), data = case_control_data_primary)

p_int_dm <- anova(model_withdummy_dm, model_withdummy_dm_2dose)


model_int_dm_1   <- clogit(case ~ vcat_combined * relevel(factor(DIABETES_m), ref = 1) + age  + n_comorb_cat+strata(stratum_no), data = case_control_data_primary)
model_int_dm_2   <- clogit(case ~ vcat_combined * relevel(factor(DIABETES_m), ref = 2) + age  + n_comorb_cat+strata(stratum_no), data = case_control_data_primary)


# Metro


model_withdummy_sp  <- clogit(case ~ 
                                factor(vcat_combined=="Dose 1 0-13 days") + 
                                factor(vcat_combined=="Dose 1 >=14 days") + 
                                factor(vcat_combined=="Dose 2 0-13 days") + 
                                factor(vcat_combined=="Dose 2 >=14 days") + 
                                age + n_comorb_cat + metro_sp + strata(stratum_no), data = case_control_data_primary)

model_withdummy_sp_2dose  <- clogit(case ~ 
                                      factor(vcat_combined=="Dose 1 0-13 days") + 
                                      factor(vcat_combined=="Dose 1 >=14 days") + 
                                      factor(vcat_combined=="Dose 2 0-13 days") + 
                                      factor(vcat_combined=="Dose 2 >=14 days") + 
                                      relevel(factor(vcat_combined=="Dose 2 >=14 days"),ref="FALSE"):metro_sp + 
                                      age + n_comorb_cat + strata(stratum_no), data = case_control_data_primary)

p_int_sp <- anova(model_withdummy_sp, model_withdummy_sp_2dose)



model_int_sp_1  <- clogit(case ~ vcat_combined * relevel(factor(metro_sp), ref = 1) + age + n_comorb_cat + strata(stratum_no), data = case_control_data_primary)
model_int_sp_2  <- clogit(case ~ vcat_combined * relevel(factor(metro_sp), ref = 2) + age + n_comorb_cat + strata(stratum_no), data = case_control_data_primary)


model_int_sex_1 <- model_int_sex_1 %>% broom.helpers::tidy_and_attach(exponentiate = TRUE,conf.int = TRUE)
model_int_sex_2 <- model_int_sex_2 %>% broom.helpers::tidy_and_attach(exponentiate = TRUE,conf.int = TRUE)
model_int_co_1  <- model_int_co_1  %>% broom.helpers::tidy_and_attach(exponentiate = TRUE,conf.int = TRUE)
model_int_co_2  <- model_int_co_2  %>% broom.helpers::tidy_and_attach(exponentiate = TRUE,conf.int = TRUE)
model_int_cardio_1  <- model_int_cardio_1  %>% broom.helpers::tidy_and_attach(exponentiate = TRUE,conf.int = TRUE)
model_int_cardio_2  <- model_int_cardio_2  %>% broom.helpers::tidy_and_attach(exponentiate = TRUE,conf.int = TRUE)
model_int_dm_1  <- model_int_dm_1  %>% broom.helpers::tidy_and_attach(exponentiate = TRUE,conf.int = TRUE)
model_int_dm_2  <- model_int_dm_2  %>% broom.helpers::tidy_and_attach(exponentiate = TRUE,conf.int = TRUE)
model_int_sp_1  <- model_int_sp_1  %>% broom.helpers::tidy_and_attach(exponentiate = TRUE,conf.int = TRUE)
model_int_sp_2  <- model_int_sp_2  %>% broom.helpers::tidy_and_attach(exponentiate = TRUE,conf.int = TRUE)

p_sex <- p_int_sex$`P(>|Chi|)`
p_cor <- p_int_comor$`P(>|Chi|)`
p_cvd <- p_int_cvd$`P(>|Chi|)`
p_dm  <- p_int_dm$`P(>|Chi|)`
p_sp  <- p_int_sp$`P(>|Chi|)`


models_int_main_groups <- 
  bind_rows(
    
    model_int_sex_1 %>% mutate(model = "Female"   , category = "Sex", p_int = p_sex[2]),
    model_int_sex_2 %>% mutate(model = "Male"     , category = "Sex", p_int = NA_real_),
    
    model_int_co_1 %>% mutate(model = "Not reported",        category = "Comor", p_int = p_cor[2],),
    model_int_co_2 %>% mutate(model = "Reported", category = "Comor",p_int = NA_real_),
    
    model_int_cardio_1 %>% mutate(model = "Not reported",         category = "CVD",  p_int = p_cvd[2]),
    model_int_cardio_2 %>% mutate(model = "Reported",  category = "CVD",  p_int = NA_real_),
    
    model_int_dm_1 %>% mutate(model = "Not reported",        category = "DM", p_int = p_dm[2]),
    model_int_dm_2 %>% mutate(model = "Reported", category = "DM", p_int = NA_real_),
    
    model_int_sp_1 %>% mutate(model = "Within Grande São Paulo Health Region", category = "Region",p_int = p_sp[2]),
    model_int_sp_2 %>% mutate(model = "Outside Grande São Paulo Health Region",category = "Region",p_int =NA_real_)
    
  )

models_int_main_groups <-
  models_int_main_groups %>% 
  filter(str_detect(term, "Dose 2 >=14 days") & 
           !str_detect(term, "relevel"))  %>% 
  mutate(OR = paste(
    round(estimate,3), " (",
    round(conf.low,3), "-",
    round(conf.high,3), ")", sep = ""),
    p_int = as.character(round(p_int,3)),
    ve = 
        paste(
          format(round(1-estimate,3) *100, nsmall=1, trim = TRUE), "% (",
          format(round(1-conf.high,3)*100, nsmall=1, trim = TRUE), "-", 
          format(round(1-conf.low,3) *100, nsmall=1, trim = TRUE), ")",
          sep = "")
  )

