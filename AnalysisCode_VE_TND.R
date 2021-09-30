
library(tidyverse)
library(tidylog)
library(lubridate)
library(tableone)
library(survival)
library(broom)
library(broom.helpers)
library(reshape2)

#### Matched data
case_control_data_primary <- read_csv('ChAdOx1_CaseControlPairs.csv')
case_control_data_primary = case_control_data_primary %>% mutate(
  vcat_combined = factor(vcat_combined,levels=c("Unvaccinated","Dose 1 0-13 days","Dose 1 14-27 days","Dose 1 >=28 days","Dose 2 0-13 days","Dose 2 >=14 days")),
  vcat_combined_biweek   = factor(ifelse((!is.na(DATA_APLICACAO_1_VACINA) & DATA_APLICACAO_1_VACINA>pcrcollection_date) | is.na(DATA_APLICACAO_1_VACINA),'Unvaccinated',
                                         ifelse((!is.na(DATA_APLICACAO_1_VACINA) & DATA_APLICACAO_1_VACINA>pcrcollection_date-14) &
                                                  (is.na(DATA_APLICACAO_2_VACINA) | (!is.na(DATA_APLICACAO_2_VACINA) & DATA_APLICACAO_2_VACINA>pcrcollection_date)),'Dose 1 0-13 days',
                                                ifelse((!is.na(DATA_APLICACAO_1_VACINA) & DATA_APLICACAO_1_VACINA>pcrcollection_date-28) &
                                                         (is.na(DATA_APLICACAO_2_VACINA) | (!is.na(DATA_APLICACAO_2_VACINA) & DATA_APLICACAO_2_VACINA>pcrcollection_date)),'Dose 1 14-27 days',
                                                       ifelse((!is.na(DATA_APLICACAO_1_VACINA) & DATA_APLICACAO_1_VACINA>pcrcollection_date-42) &
                                                                (is.na(DATA_APLICACAO_2_VACINA) | (!is.na(DATA_APLICACAO_2_VACINA) & DATA_APLICACAO_2_VACINA>pcrcollection_date)),'Dose 1 28-41 days',
                                                              ifelse((!is.na(DATA_APLICACAO_1_VACINA) & DATA_APLICACAO_1_VACINA>pcrcollection_date-56) &
                                                                       (is.na(DATA_APLICACAO_2_VACINA) | (!is.na(DATA_APLICACAO_2_VACINA) & DATA_APLICACAO_2_VACINA>pcrcollection_date)),'Dose 1 42-55 days',
                                                                     ifelse((!is.na(DATA_APLICACAO_1_VACINA) & DATA_APLICACAO_1_VACINA<=pcrcollection_date-56) &
                                                                              (is.na(DATA_APLICACAO_2_VACINA) | (!is.na(DATA_APLICACAO_2_VACINA) & DATA_APLICACAO_2_VACINA>pcrcollection_date)),'Dose 1 >=56 days',
                                                                            ifelse(!is.na(DATA_APLICACAO_2_VACINA) & DATA_APLICACAO_2_VACINA>pcrcollection_date-14,'Dose 2 0-13 days','Dose 2 >=14 days'
                                                                            ))))))))
  )
case_control_data_primary = case_control_data_primary %>% mutate(
  vcat_combined_biweek = factor(vcat_combined_biweek,
                                levels = c("Unvaccinated", "Dose 1 0-13 days",
                                           "Dose 1 14-27 days", "Dose 1 28-41 days",
                                           "Dose 1 42-55 days","Dose 1 >=56 days",
                                           "Dose 2 0-13 days","Dose 2 >=14 days"))
)

## Vector of variables to summarize
myVars <- c(
            "age_5y", "sex", "race","metro_sp",
            "n_comorb_cat", "CARDIOPATI_m", "DIABETES_m",
            "visits_cat","previous_inf", "vcat_combined",
            "hosp_yn","icu","resp_supp_imv","death_yn"
)


## Vector of categorical variables that need transformation
catVars <- c(
  "age_5y", "sex", "race","metro_sp",
  "n_comorb_cat", "CARDIOPATI_m", "DIABETES_m",
  "visits_cat","previous_inf", "vcat_combined",
  "hosp_yn","icu","resp_supp_imv","death_yn"
)



tab_match <- tableone::CreateTableOne(
  vars = myVars, factorVars = catVars,
  strata = "case", includeNA = TRUE,
  smd = FALSE, addOverall = FALSE,
  data = case_control_data_primary)

tab_match_export <- print(tab_match,
                       quote = FALSE, noSpaces = TRUE, printToggle = FALSE, 
                       smd = FALSE, test = FALSE,
                       showAllLevels = FALSE)

#### Running analysis ####
model_crude_tbl <-
  clogit(case ~ vcat_combined + strata(stratum_no), data = case_control_data_primary) %>% 
  broom.helpers::tidy_and_attach(exponentiate = TRUE,conf.int = TRUE)

model_adj_tbl <-
  clogit(case ~ vcat_combined + n_comorb_cat + previous_inf + strata(stratum_no), data = case_control_data_primary) %>% 
  broom.helpers::tidy_and_attach(exponentiate = TRUE,conf.int = TRUE)


models <- bind_rows(model_crude_tbl %>% mutate(model = "Crude"), 
                    model_adj_tbl   %>% mutate(model = "Adj"))

# Primary models
models <-
  models %>% 
  mutate(OR = paste(
         round(estimate,2), " (",
         round(conf.low,2), "-",
         round(conf.high,2), ")", sep = ""),
         p_value = as.character(round(p.value,3)),
         ve = case_when(
           str_detect(term, "Dose") ~
           paste(
           round(1-estimate,3)*100, "% (",
           round(1-conf.high,3)*100, "-", 
           round(1-conf.low,3)*100, ")",
           sep = "")))

#### Effectiveness by bi-week following first dose (data for Figure 3)
model_crude_tbl <-
  clogit(case ~ vcat_combined_biweek + strata(stratum_no), data = case_control_data_primary) %>% 
  broom.helpers::tidy_and_attach(exponentiate = TRUE,conf.int = TRUE)

model_adj_tbl <-
  clogit(case ~ vcat_combined_biweek +  n_comorb_cat + previous_inf + strata(stratum_no), data = case_control_data_primary) %>% 
  broom.helpers::tidy_and_attach(exponentiate = TRUE,conf.int = TRUE)

model_adj_death_week <-
  clogit(case ~ vcat_combined_biweek +  n_comorb_cat + previous_inf + strata(stratum_no), data = case_control_data_primary %>% filter(death_pairs==1)) %>% 
  broom.helpers::tidy_and_attach(exponentiate = TRUE,conf.int = TRUE)


models <- bind_rows(model_crude_tbl %>% mutate(model = "Crude"), 
                    model_adj_tbl   %>% mutate(model = "Adj"),
                    model_adj_hosp_week %>% mutate(model="Hosp Adj"),
                    model_adj_death_week %>% mutate(model="Death Adj")
)

models <-
  models %>% 
  mutate(OR = paste(
    round(estimate,2), " (",
    round(conf.low,2), "-",
    round(conf.high,2), ")", sep = ""),
    p_value = as.character(round(p.value,3)),
    ve = case_when(
      str_detect(term, "Dose") ~
        paste(
          sprintf('%.1f',round(1-estimate,3)*100), "% (",
          sprintf('%.1f',round(1-conf.high,3)*100), "-", 
          sprintf('%.1f',round(1-conf.low,3)*100), ")",
          sep = "")),
    ve_plot  = (1-estimate),
    ve_plot_low   = (1-conf.high),
    ve_plot_high  = (1-conf.low))

#### Effectiveness against severe outcomes

# separating stratas
strata_cases_hosp <- 
  case_control_data_primary %>% 
  filter(hospitalization==1 & case == 1) %>% 
  select(stratum_no)

strata_cases_icu <- 
  case_control_data_primary %>% 
  filter(icu==1 & case == 1) %>% 
  select(stratum_no)

strata_cases_imv <- 
  case_control_data_primary %>% 
  filter(resp_supp_imv==1 & case == 1) %>% 
  select(stratum_no)

strata_cases_death <- 
  case_control_data_primary %>% 
  filter(death_2==1 & case == 1) %>% 
  select(stratum_no)

case_control_data_primary <-
  case_control_data_primary %>% 
  mutate(
    hosp_pairs = case_when(
      stratum_no %in% strata_cases_hosp$stratum_no   ~ 1),
    icu_pairs = case_when(
      stratum_no %in% strata_cases_icu$stratum_no    ~ 1),
    imv_pairs = case_when(
      stratum_no %in% strata_cases_imv$stratum_no    ~ 1),
    death_pairs = case_when(
      stratum_no %in% strata_cases_death$stratum_no  ~ 1),
  )

case_control_data_primary %>% filter(hosp_pairs==1) %>% distinct(id) %>% count()
case_control_data_primary %>% filter(death_pairs==1) %>% distinct(id) %>% count()


## main models

#### Running analysis ####

# Hospitalization
model_crude_hosp <-
  clogit(case ~ vcat_combined + strata(stratum_no), data = case_control_data_primary %>% filter(hosp_pairs==1)) %>% 
  broom.helpers::tidy_and_attach(exponentiate = TRUE,conf.int = TRUE)

model_adj_hosp <-
  clogit(case ~ vcat_combined +  n_comorb_cat + previous_inf + strata(stratum_no), data = case_control_data_primary %>% filter(hosp_pairs==1)) %>% 
  broom.helpers::tidy_and_attach(exponentiate = TRUE,conf.int = TRUE)

# ICU admission
model_crude_icu <-
  clogit(case ~ vcat_combined + strata(stratum_no), data = case_control_data_primary %>% filter(icu_pairs==1)) %>% 
  broom.helpers::tidy_and_attach(exponentiate = TRUE,conf.int = TRUE)

model_adj_icu <-
  clogit(case ~ vcat_combined +  n_comorb_cat + previous_inf + strata(stratum_no), data = case_control_data_primary %>% filter(icu_pairs==1)) %>% 
  broom.helpers::tidy_and_attach(exponentiate = TRUE,conf.int = TRUE)

# Mechanical ventilation
model_crude_imv <-
  clogit(case ~ vcat_combined + strata(stratum_no), data = case_control_data_primary %>% filter(imv_pairs==1)) %>% 
  broom.helpers::tidy_and_attach(exponentiate = TRUE,conf.int = TRUE)

model_adj_imv <-
  clogit(case ~ vcat_combined +  n_comorb_cat + previous_inf + strata(stratum_no), data = case_control_data_primary %>% filter(imv_pairs==1)) %>% 
  broom.helpers::tidy_and_attach(exponentiate = TRUE,conf.int = TRUE)

# Death
model_crude_death <-
  clogit(case ~ vcat_combined + strata(stratum_no), data = case_control_data_primary %>% filter(death_pairs==1)) %>% 
  broom.helpers::tidy_and_attach(exponentiate = TRUE,conf.int = TRUE)

model_adj_death <-
  clogit(case ~ vcat_combined +  n_comorb_cat + previous_inf + strata(stratum_no), data = case_control_data_primary %>% filter(death_pairs==1)) %>% 
  broom.helpers::tidy_and_attach(exponentiate = TRUE,conf.int = TRUE)


models_severity <- bind_rows(
  
  model_crude_hosp   %>% mutate(model = "Crude", outcome = "Hospitalization", n = nrow(strata_cases_hosp)), 
  model_adj_hosp     %>% mutate(model = "Adj"  , outcome = "Hospitalization", n = nrow(strata_cases_hosp)),
  
  model_crude_icu    %>% mutate(model = "Crude", outcome = "ICU", n = nrow(strata_cases_icu)), 
  model_adj_icu      %>% mutate(model = "Adj"  , outcome = "ICU", n = nrow(strata_cases_icu)),
  
  model_crude_imv    %>% mutate(model = "Crude", outcome = "iMV", n = nrow(strata_cases_imv)), 
  model_adj_imv      %>% mutate(model = "Adj"  , outcome = "iMV", n = nrow(strata_cases_imv)),
  
  model_crude_death  %>% mutate(model = "Crude", outcome = "Death", n = nrow(strata_cases_death)),
  model_adj_death    %>% mutate(model = "Adj"  , outcome = "Death", n = nrow(strata_cases_death)),
  
)

# Models for severe outcomes (primary results)
models_severity <-
  models_severity %>% 
  mutate(OR = paste(
    round(estimate,2), " (",
    round(conf.low,2), "-",
    round(conf.high,2), ")", sep = ""),
    p_value = as.character(round(p.value,3)),
    ve = case_when(
      str_detect(term, "Dose") ~
        paste(
          round(1-estimate,3)*100, "% (",
          round(1-conf.high,3)*100, "-", 
          round(1-conf.low,3)*100, ")",
          sep = "")))


### Subgroup analyses for 1st dose ###
# Age
model_withdummy_agep_1dose  <- clogit(case ~ 
                                        factor(vcat_combined=="Dose 1 0-13 days") + 
                                        factor(vcat_combined=="Dose 1 14-27 days") + 
                                        factor(vcat_combined=="Dose 1 >=28 days") + 
                                        factor(vcat_combined=="Dose 2 0-13 days") + 
                                        factor(vcat_combined=="Dose 2 >=14 days") + 
                                        relevel(factor(vcat_combined=="Dose 1 >=28 days"),ref="FALSE"):factor(age_70) + 
                                         n_comorb_cat + previous_inf + strata(stratum_no), data = case_control_data_primary)

p_int_agep <- anova(model_withdummy_agep, model_withdummy_agep_1dose)


# Sex
model_withdummy_sex_1dose  <- clogit(case ~ 
                                       factor(vcat_combined=="Dose 1 0-13 days") + 
                                       factor(vcat_combined=="Dose 1 14-27 days") + 
                                       factor(vcat_combined=="Dose 1 >=28 days") + 
                                       factor(vcat_combined=="Dose 2 0-13 days") + 
                                       factor(vcat_combined=="Dose 2 >=14 days") + 
                                       relevel(factor(vcat_combined=="Dose 1 >=28 days"),ref="FALSE"):sex + 
                                        n_comorb_cat + previous_inf + strata(stratum_no), data = case_control_data_primary)

p_int_sex <- anova(model_withdummy_sex, model_withdummy_sex_1dose)


#Comorb
model_withdummy_comor_1dose  <- clogit(case ~ factor(vcat_combined=="Dose 1 0-13 days") + 
                                         factor(vcat_combined=="Dose 1 14-27 days") + 
                                         factor(vcat_combined=="Dose 1 >=28 days") + 
                                         factor(vcat_combined=="Dose 2 0-13 days") + 
                                         factor(vcat_combined=="Dose 2 >=14 days") + 
                                         relevel(factor(vcat_combined=="Dose 1 >=28 days"),ref="FALSE"):anycomorb + 
                                          previous_inf + strata(stratum_no), data = case_control_data_primary)

p_int_comor <- anova(model_withdummy_comor, model_withdummy_comor_1dose)


# Cardio
model_withdummy_cvd_1dose  <- clogit(case ~ 
                                       factor(vcat_combined=="Dose 1 0-13 days") + 
                                       factor(vcat_combined=="Dose 1 14-27 days") + 
                                       factor(vcat_combined=="Dose 1 >=28 days") + 
                                       factor(vcat_combined=="Dose 2 0-13 days") + 
                                       factor(vcat_combined=="Dose 2 >=14 days") + 
                                       relevel(factor(vcat_combined=="Dose 1 >=28 days"),ref="FALSE"):CARDIOPATI_m + 
                                        n_comorb_cat + previous_inf + strata(stratum_no), data = case_control_data_primary)

p_int_cvd <- anova(model_withdummy_cvd, model_withdummy_cvd_1dose)

# Diabetes
model_withdummy_dm_1dose  <- clogit(case ~ 
                                      factor(vcat_combined=="Dose 1 0-13 days") + 
                                      factor(vcat_combined=="Dose 1 14-27 days") + 
                                      factor(vcat_combined=="Dose 1 >=28 days") + 
                                      factor(vcat_combined=="Dose 2 0-13 days") + 
                                      factor(vcat_combined=="Dose 2 >=14 days") + 
                                      relevel(factor(vcat_combined=="Dose 1 >=28 days"),ref="FALSE"):DIABETES_m + 
                                       n_comorb_cat + previous_inf + strata(stratum_no), data = case_control_data_primary)

p_int_dm <- anova(model_withdummy_dm, model_withdummy_dm_1dose)

# Metro
model_withdummy_sp_1dose  <- clogit(case ~ 
                                      factor(vcat_combined=="Dose 1 0-13 days") + 
                                      factor(vcat_combined=="Dose 1 14-27 days") + 
                                      factor(vcat_combined=="Dose 1 >=28 days") + 
                                      factor(vcat_combined=="Dose 2 0-13 days") + 
                                      factor(vcat_combined=="Dose 2 >=14 days") + 
                                      relevel(factor(vcat_combined=="Dose 1 >=28 days"),ref="FALSE"):metro_sp + 
                                       n_comorb_cat + previous_inf + strata(stratum_no), data = case_control_data_primary)

p_int_sp <- anova(model_withdummy_sp, model_withdummy_sp_1dose)

p_age <- p_int_agep$`P(>|Chi|)`
p_sex <- p_int_sex$`P(>|Chi|)`
p_cor <- p_int_comor$`P(>|Chi|)`
p_cvd <- p_int_cvd$`P(>|Chi|)`
p_dm  <- p_int_dm$`P(>|Chi|)`
p_sp  <- p_int_sp$`P(>|Chi|)`


models_int <- 
  bind_rows(
    model_int_age_1 %>% mutate(model = "60-69", category = "Age", p_int = p_age[2], n = nage1),
    model_int_age_2 %>% mutate(model = ">=70", category = "Age", p_int = NA_real_, n = nage2),
    
    model_int_sex_1 %>% mutate(model = "Female"   , category = "Sex", p_int = p_sex[2], n = nsex1),
    model_int_sex_2 %>% mutate(model = "Male"     , category = "Sex", p_int = NA_real_, n = nsex2),
    
    model_int_co_1 %>% mutate(model = "Not reported",        category = "Comor", p_int = p_cor[2], n = nco1),
    model_int_co_2 %>% mutate(model = "Reported", category = "Comor",p_int = NA_real_, n = nco2),
    
    model_int_cardio_1 %>% mutate(model = "Not reported",         category = "CVD",  p_int = p_cvd[2], n = ncvd1),
    model_int_cardio_2 %>% mutate(model = "Reported",  category = "CVD",  p_int = NA_real_, n = ncvd2),
    
    model_int_dm_1 %>% mutate(model = "Not reported",        category = "DM", p_int = p_dm[2], n = ndm1),
    model_int_dm_2 %>% mutate(model = "Reported", category = "DM", p_int = NA_real_, n = ndm2),
    
    model_int_sp_1 %>% mutate(model = "Within Grande São Paulo Health Region", category = "Region",p_int = p_sp[2], n = nsp1),
    model_int_sp_2 %>% mutate(model = "Outside Grande São Paulo Health Region",category = "Region",p_int =NA_real_, n = nsp2)
    
  )

# Models for subgroup analyses (Supplementary Table 10)
models_int <-
  models_int %>% 
  mutate(OR = paste(
    round(estimate,2), " (",
    round(conf.low,2), "-",
    round(conf.high,2), ")", sep = ""),
    p_value = as.character(round(p.value,3)),
    p_int = as.character(round(p_int,2)),
    ve = case_when(
      str_detect(term, "Dose") ~
        paste(
          format(round(1-estimate,3) *100, nsmall=1, trim = TRUE), "% (",
          format(round(1-conf.high,3)*100, nsmall=1, trim = TRUE), "-", 
          format(round(1-conf.low,3) *100, nsmall=1, trim = TRUE), ")",
          sep = "")),
    ve_plot  = (1-estimate),
    ve_plot_low   = (1-conf.high),
    ve_plot_high  = (1-conf.low),
  )


models_int_plot <-
  models_int %>% 
  filter(
    term == "vcat_combinedDose 1 >=28 days" 
  ) %>% 
  mutate(term = "Dose 1 >=28 days"
  ) %>% 
  select(term, model, category, ve, ve_plot, ve_plot_low, ve_plot_high, p_int, n)


# Effectiveness against severe outcomes by age, single dose and two doses

# hosp by age
model_withdummy_agep_hosp  <- clogit(case ~ 
                                         factor(vcat_combined=="Dose 1 0-13 days") + 
                                         factor(vcat_combined=="Dose 1 14-27 days") + 
                                       factor(vcat_combined=="Dose 1 >=28 days") + 
                                       factor(vcat_combined=="Dose 2 0-13 days") + 
                                         factor(vcat_combined=="Dose 2 >=14 days") + 
                                          factor(age_70) + n_comorb_cat  + strata(stratum_no), 
                                       data = case_control_data_primary %>% filter(hosp_pairs == 1))

model_withdummy_agep_1dose_hosp <- clogit(case ~ 
                                            factor(vcat_combined=="Dose 1 0-13 days") + 
                                            factor(vcat_combined=="Dose 1 14-27 days") + 
                                            factor(vcat_combined=="Dose 1 >=28 days") + 
                                            factor(vcat_combined=="Dose 2 0-13 days") + 
                                            factor(vcat_combined=="Dose 2 >=14 days") + 
                                            relevel(factor(vcat_combined=="Dose 1 >=28 days"),ref="FALSE"):factor(age_70) + 
                                             n_comorb_cat  + strata(stratum_no), 
                                          data = case_control_data_primary %>% filter(hosp_pairs==1))

model_withdummy_agep_2dose_hosp <- clogit(case ~ 
                                            factor(vcat_combined=="Dose 1 0-13 days") + 
                                            factor(vcat_combined=="Dose 1 14-27 days") + 
                                            factor(vcat_combined=="Dose 1 >=28 days") + 
                                            factor(vcat_combined=="Dose 2 0-13 days") + 
                                            factor(vcat_combined=="Dose 2 >=14 days") + 
                                            relevel(factor(vcat_combined=="Dose 2 >=14 days"),ref="FALSE"):factor(age_70) + 
                                             n_comorb_cat  + strata(stratum_no), 
                                          data = case_control_data_primary %>% filter(hosp_pairs==1))

p_int_agep_hosp_1dose <- anova(model_withdummy_agep_hosp, model_withdummy_agep_1dose_hosp)
p_int_agep_hosp_2dose <- anova(model_withdummy_agep_hosp, model_withdummy_agep_2dose_hosp)


model_int_age_1_hosp  <- clogit(case ~ vcat_combined * relevel(factor(age_70), ref = 1) +   n_comorb_cat  + strata(stratum_no), data = case_control_data_primary %>% filter(hosp_pairs==1))
model_int_age_2_hosp  <- clogit(case ~ vcat_combined * relevel(factor(age_70), ref = 2) +   n_comorb_cat  + strata(stratum_no), data = case_control_data_primary %>% filter(hosp_pairs==1))


model_withdummy_agep_death  <- clogit(case ~ 
                                        factor(vcat_combined=="Dose 1 0-13 days") + 
                                        factor(vcat_combined=="Dose 1 14-27 days") + 
                                        factor(vcat_combined=="Dose 1 >=28 days") + 
                                        factor(vcat_combined=="Dose 2 0-13 days") + 
                                        factor(vcat_combined=="Dose 2 >=14 days") + 
                                         factor(age_70) + n_comorb_cat  + strata(stratum_no), 
                                    data = case_control_data_primary %>% filter(death_pairs == 1))

model_withdummy_agep_1dose_death <- clogit(case ~ 
                                             factor(vcat_combined=="Dose 1 0-13 days") + 
                                             factor(vcat_combined=="Dose 1 14-27 days") + 
                                             factor(vcat_combined=="Dose 1 >=28 days") + 
                                             factor(vcat_combined=="Dose 2 0-13 days") + 
                                             factor(vcat_combined=="Dose 2 >=14 days") + 
                                             relevel(factor(vcat_combined=="Dose 1 >=28 days"),ref="FALSE"):factor(age_70) + 
                                            n_comorb_cat  + strata(stratum_no), 
                                         data = case_control_data_primary %>% filter(death_pairs==1))

p_int_agep_death_1dose <- anova(model_withdummy_agep_death, model_withdummy_agep_1dose_death)

model_withdummy_agep_2dose_death <- clogit(case ~ 
                                             factor(vcat_combined=="Dose 1 0-13 days") + 
                                             factor(vcat_combined=="Dose 1 14-27 days") + 
                                             factor(vcat_combined=="Dose 1 >=28 days") + 
                                             factor(vcat_combined=="Dose 2 0-13 days") + 
                                             factor(vcat_combined=="Dose 2 >=14 days") + 
                                             relevel(factor(vcat_combined=="Dose 2 >=14 days"),ref="FALSE"):factor(age_70) + 
                                              n_comorb_cat  + strata(stratum_no), 
                                           data = case_control_data_primary %>% filter(death_pairs==1))

p_int_agep_death_2dose <- anova(model_withdummy_agep_death, model_withdummy_agep_2dose_death)


model_int_age_1_death  <- clogit(case ~ vcat_combined * relevel(factor(age_70), ref = 1) +   n_comorb_cat + strata(stratum_no), data = case_control_data_primary %>% filter(death_pairs==1))
model_int_age_2_death  <- clogit(case ~ vcat_combined * relevel(factor(age_70), ref = 2) +   n_comorb_cat + strata(stratum_no), data = case_control_data_primary %>% filter(death_pairs==1))

## tidying
model_int_age_1_hosp <- model_int_age_1_hosp %>% broom.helpers::tidy_and_attach(exponentiate = TRUE,conf.int = TRUE)
model_int_age_2_hosp <- model_int_age_2_hosp %>% broom.helpers::tidy_and_attach(exponentiate = TRUE,conf.int = TRUE)

model_int_age_1_death <- model_int_age_1_death %>% broom.helpers::tidy_and_attach(exponentiate = TRUE,conf.int = TRUE)
model_int_age_2_death <- model_int_age_2_death %>% broom.helpers::tidy_and_attach(exponentiate = TRUE,conf.int = TRUE)


p_age_hosp   <- p_int_agep_hosp_1dose$`P(>|Chi|)`
p_age_death  <- p_int_agep_death_1dose$`P(>|Chi|)`

nage1_hosp <- case_control_data_primary %>% filter(hosp_pairs==1) %>% count(age_70) %>% filter(age_70==0) %>% select(n) %>% as.numeric()
nage2_hosp <- case_control_data_primary %>% filter(hosp_pairs==1) %>% count(age_70) %>% filter(age_70==1) %>% select(n) %>% as.numeric()

nage1_death <- case_control_data_primary %>% filter(death_pairs==1) %>% count(age_70) %>% filter(age_70==0) %>% select(n) %>% as.numeric()
nage2_death <- case_control_data_primary %>% filter(death_pairs==1) %>% count(age_70) %>% filter(age_70==1) %>% select(n) %>% as.numeric()



models_int_severity <- 
  bind_rows(
    model_int_age_1_hosp %>% mutate(model = "60-69", outcome = "Hospitalized", p_int = p_age_hosp[2],   n = nage1_hosp),
    model_int_age_2_hosp %>% mutate(model = "70+", outcome = "Hospitalized", p_int = NA_real_,        n = nage2_hosp),
    
    model_int_age_1_death %>% mutate(model = "60-69", outcome = "Deaths", p_int = p_age_death[2],  n = nage1_death),
    model_int_age_2_death %>% mutate(model = "70+", outcome = "Deaths", p_int = NA_real_,        n = nage2_death)
    )

# Subgroup analyses by age for severe outcomes, single dose (Supplementary Table 11)
models_int_severity <-
  models_int_severity %>% 
  mutate(OR = paste(
    round(estimate,2), " (",
    round(conf.low,2), "-",
    round(conf.high,2), ")", sep = ""),
    p_value = as.character(round(p.value,3)),
    p_int = as.character(round(p_int,2)),
    ve = case_when(
      str_detect(term, "Dose") ~
        paste(
          format(round(1-estimate,3) *100, nsmall=1, trim = TRUE, scientific = FALSE), "% (",
          format(round(1-conf.high,3)*100, nsmall=1, trim = TRUE, scientific = FALSE), "-", 
          format(round(1-conf.low,3) *100, nsmall=1, trim = TRUE, scientific = FALSE), ")",
          sep = "")),
    ve_plot  = (1-estimate),
    ve_plot_low   = (1-conf.high),
    ve_plot_high  = (1-conf.low),
  )


models_int_severity_plot <-
  models_int_severity %>% 
  filter(
    term == "vcat_combinedDose 1 >=28 days" 
  ) %>% 
  mutate(term = "Dose 1 >=28 days"
  ) %>% 
  select(term, model, outcome, ve, ve_plot, ve_plot_low, ve_plot_high, p_int, n)


p_age_hosp   <- p_int_agep_hosp_2dose$`P(>|Chi|)`
p_age_icu    <- p_int_agep_icu_2dose$`P(>|Chi|)`
p_age_imv    <- p_int_agep_imv_2dose$`P(>|Chi|)`
p_age_death  <- p_int_agep_death_2dose$`P(>|Chi|)`

models_int_severity <- 
  bind_rows(
    model_int_age_1_hosp %>% mutate(model = "60-69", outcome = "Hospitalized", p_int = p_age_hosp[2],   n = nage1_hosp),
    model_int_age_2_hosp %>% mutate(model = "70+", outcome = "Hospitalized", p_int = NA_real_,        n = nage2_hosp),
    
    model_int_age_1_icu %>% mutate(model = "60-69", outcome = "Admitted to the ICU", p_int = p_age_icu[2],    n = nage1_icu),
    model_int_age_2_icu %>% mutate(model = "70+", outcome = "Admitted to the ICU", p_int = NA_real_,        n = nage2_icu),
    
    model_int_age_1_imv %>% mutate(model = "60-69", outcome = "Mechanically ventilated", p_int = p_age_imv[2],    n = nage1_imv),
    model_int_age_2_imv %>% mutate(model = "70+", outcome = "Mechanically ventilated", p_int = NA_real_,        n = nage2_imv),
    
    model_int_age_1_death %>% mutate(model = "60-69", outcome = "Deaths", p_int = p_age_death[2],  n = nage1_death),
    model_int_age_2_death %>% mutate(model = "70+", outcome = "Deaths", p_int = NA_real_,        n = nage2_death)
  )

# Effectiveness against severe outcomes by age, two doses (Supplementary Table 11)
models_int_severity <-
  models_int_severity %>% 
  mutate(OR = paste(
    round(estimate,2), " (",
    round(conf.low,2), "-",
    round(conf.high,2), ")", sep = ""),
    p_value = as.character(round(p.value,3)),
    p_int = as.character(round(p_int,2)),
    ve = case_when(
      str_detect(term, "Dose") ~
        paste(
          format(round(1-estimate,3) *100, nsmall=1, trim = TRUE, scientific = FALSE), "% (",
          format(round(1-conf.high,3)*100, nsmall=1, trim = TRUE, scientific = FALSE), "-", 
          format(round(1-conf.low,3) *100, nsmall=1, trim = TRUE, scientific = FALSE), ")",
          sep = "")),
    ve_plot  = (1-estimate),
    ve_plot_low   = (1-conf.high),
    ve_plot_high  = (1-conf.low),
  )

models_int_severity_plot <-
  models_int_severity %>% 
  filter(
    term == "vcat_combinedDose 2 >=14 days" 
  ) %>% 
  mutate(term = "Dose 2 >=14 days"
  ) %>% 
  select(term, model, outcome, ve, ve_plot, ve_plot_low, ve_plot_high, p_int, n)


# Make supplementary figures
dps=case_control_data_primary %>% mutate(week=epiweek(pcrcollection_date)) %>% dplyr::group_by(stratum_no) %>% 
  dplyr::summarise(dp1v0=as.numeric(sum(vcat_combined=="Unvaccinated")==1 & sum(vcat_combined=="Dose 1 0-13 days")==1)
                   ,dp2v0=as.numeric(sum(vcat_combined=="Unvaccinated")==1 & sum(vcat_combined=="Dose 1 14-27 days")==1)
                   ,dp3v0=as.numeric(sum(vcat_combined=="Unvaccinated")==1 & sum(vcat_combined=="Dose 1 >=28 days")==1)
                   ,dp4v0=as.numeric(sum(vcat_combined=="Unvaccinated")==1 & sum(vcat_combined=="Dose 2 0-13 days")==1)
                   ,dp5v0=as.numeric(sum(vcat_combined=="Unvaccinated")==1 & sum(vcat_combined=="Dose 2 >=14 days")==1)
                   ,dp2v1=as.numeric(sum(vcat_combined=="Dose 1 0-13 days")==1 & sum(vcat_combined=="Dose 1 14-27 days")==1)
                   ,dp3v1=as.numeric(sum(vcat_combined=="Dose 1 0-13 days")==1 & sum(vcat_combined=="Dose 1 >=28 days")==1)
                   ,dp4v1=as.numeric(sum(vcat_combined=="Dose 1 0-13 days")==1 & sum(vcat_combined=="Dose 2 0-13 days")==1)
                   ,dp5v1=as.numeric(sum(vcat_combined=="Dose 1 0-13 days")==1 & sum(vcat_combined=="Dose 2 >=14 days")==1)
                   ,dp3v2=as.numeric(sum(vcat_combined=="Dose 1 14-27 days")==1 & sum(vcat_combined=="Dose 1 >=28 days")==1)
                   ,dp4v2=as.numeric(sum(vcat_combined=="Dose 1 14-27 days")==1 & sum(vcat_combined=="Dose 2 0-13 days")==1)
                   ,dp5v2=as.numeric(sum(vcat_combined=="Dose 1 14-27 days")==1 & sum(vcat_combined=="Dose 2 >=14 days")==1)
                   ,dp4v3=as.numeric(sum(vcat_combined=="Dose 1 >=28 days")==1 & sum(vcat_combined=="Dose 2 0-13 days")==1)
                   ,dp5v3=as.numeric(sum(vcat_combined=="Dose 1 >=28 days")==1 & sum(vcat_combined=="Dose 2 >=14 days")==1)
                   ,dp5v4=as.numeric(sum(vcat_combined=="Dose 2 0-13 days")==1 & sum(vcat_combined=="Dose 2 >=14 days")==1)
                   ,d=max(pcrcollection_date)
                   ,w=max(week)
  )
dps$dp_any=dps$dp1v0+dps$dp2v0+dps$dp3v0+dps$dp4v0+dps$dp2v1+dps$dp3v1+dps$dp4v1+dps$dp3v2+dps$dp4v2+dps$dp4v3


# Table of discordant pairs (Supplementary Table 1)
discpairtable=sapply(c("Unvaccinated","Dose 1 0-13 days","Dose 1 14-27 days","Dose 1 >=28 days","Dose 2 0-13 days","Dose 2 >=14 days"),
                     function(x) sapply(c("Unvaccinated","Dose 1 0-13 days","Dose 1 14-27 days","Dose 1 >=28 days","Dose 2 0-13 days","Dose 2 >=14 days"),
                                        function(y) sum(case_control_data_primary %>% dplyr::group_by(stratum_no) %>%
                                                          dplyr::summarise(s=as.numeric(sum(vcat_combined==x & case==1)==1 & sum(vcat_combined==y & case==0)==1)) %>% select(s))))

# Table of discordant hospitalized pairs (Supplementary Table 2)
discpairtable_hosp=sapply(c("Unvaccinated","Dose 1 0-13 days","Dose 1 14-27 days","Dose 1 >=28 days","Dose 2 0-13 days","Dose 2 >=14 days"),
                     function(x) sapply(c("Unvaccinated","Dose 1 0-13 days","Dose 1 14-27 days","Dose 1 >=28 days","Dose 2 0-13 days","Dose 2 >=14 days"),
                                        function(y) sum(case_control_data_primary %>% filter(hosp_pairs==1) %>% dplyr::group_by(stratum_no) %>%
                                                          dplyr::summarise(s=as.numeric(sum(vcat_combined==x & case==1)==1 & sum(vcat_combined==y & case==0)==1)) %>% select(s))))

# Table of discordant died pairs (Supplementary Table 3)
discpairtable_death=sapply(c("Unvaccinated","Dose 1 0-13 days","Dose 1 14-27 days","Dose 1 >=28 days","Dose 2 0-13 days","Dose 2 >=14 days"),
                          function(x) sapply(c("Unvaccinated","Dose 1 0-13 days","Dose 1 14-27 days","Dose 1 >=28 days","Dose 2 0-13 days","Dose 2 >=14 days"),
                                             function(y) sum(case_control_data_primary %>% filter(death_pairs==1) %>% dplyr::group_by(stratum_no) %>%
                                                               dplyr::summarise(s=as.numeric(sum(vcat_combined==x & case==1)==1 & sum(vcat_combined==y & case==0)==1)) %>% select(s))))


theme_set(theme_classic())
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
cols = gg_color_hue(15)
# Supplementary Figure 1
ggplot(data = melt(dps %>% dplyr::group_by(w) %>% 
                     dplyr::summarise(d01=sum(dp1v0),d02=sum(dp2v0),d03=sum(dp3v0),d04=sum(dp4v0),d05=sum(dp5v0),
                                      d12=sum(dp2v1),d13=sum(dp3v1),d14=sum(dp4v1),dp15=sum(dp5v1),
                                      d23=sum(dp3v2),d24=sum(dp4v2),d25=sum(dp5v2),
                                      d34=sum(dp4v3),d35=sum(dp5v3),
                                      d45=sum(dp5v4)),
                   id='w'),
       aes(x=w,y=value,fill=variable))+geom_bar(position='stack',stat='identity',colour='black')+
         ylab('Number of discordant pairs')+ylim(c(0,400))+
         scale_x_continuous(breaks=seq(1,25,2),labels=seq(1,25,2),name = 'Weeks from 17 January 2021')+
         scale_fill_manual(labels=c('1 dose 0-13 days vs. Unvaccinated',
                                    '1 dose 14-27 days vs. Unvaccinated',
                                    '1 dose 28+ days vs. Unvaccinated',
                                    '2 doses 0-13 days vs. Unvaccinated',
                                    '2 doses 14+ days vs. Unvaccinated',
                                    '1 dose 14-27 days vs. 1 dose 0-13 days',
                                    '1 dose 28+ days vs. 1 dose 0-13 days',
                                    '2 doses 0-13 days vs. 1 dose 0-13 days',
                                    '2 doses 14+ days vs. 1 dose 0-13 days',
                                    '1 dose 28+ days vs. 1 dose 14-27 days',
                                    '2 doses 0-13 days vs. 1 dose 14-27 days',
                                    '2 doses 14+ days vs. 1 dose 14-27 days',
                                    '2 doses 0-13 days vs. 1 dose 28+ days',
                                    '2 doses 14+ days vs. 1 dose 28+ days',
                                    '2 doses 14+ days vs. 2 doses 0-13 days'),
                           values=cols,name='Discordant pair type')+
         theme(legend.text = element_text(size=6),
               legend.title = element_text(size=6))

data_for_figs2 = case_control_data_primary %>% select(case,pcrcollection_date,DATA_APLICACAO_1_VACINA,
                                                      DATA_APLICACAO_2_VACINA)
# Supplementary Figure 2
ggplot(data_for_figs2 %>% filter(!is.na(DATA_APLICACAO_1_VACINA)) %>% 
         mutate(td = as.numeric(pcrcollection_date-DATA_APLICACAO_1_VACINA),
                cc=ifelse(case==1,"Cases","Controls"))) + 
  geom_histogram(aes(x=td),fill='white',col='black') + 
         facet_wrap(~cc,ncol=1) + 
         theme_classic(base_size = 6.5)+
         labs(x = "Days from 1st vaccine dose to sample collection", y = "Frequency")

ggplot(data_for_figs2 %>% mutate(td = as.numeric(pcrcollection_date-DATA_APLICACAO_2_VACINA),
                                 cc=ifelse(case==1,"Cases","Controls"))) + 
         geom_histogram(aes(x=td),fill='white',col='black') + 
         facet_wrap(~cc,ncol=1) + 
         theme_classic(base_size = 6.5)+
         labs(x = "Days from 2nd vaccine dose to sample collection", y = "Frequency")

