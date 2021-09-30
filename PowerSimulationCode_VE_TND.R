# Supporting code for "ChAdOx1 vaccine effectiveness in older adults during SARS-CoV-2 Gamma variant transmission in Sao Paulo State"
# Matt Hitchings
# No analytical power formula applies to 1:1 matched case-control study with categorical exposure
# Code is to estimate the power of a given case-control data set with specific distribution of discordant pairs through simulation

require(dplyr)
require(survival)

# Inputs: matched case control data set
# Outline: randomly assign vaccination data based on assumed odds ratios and estimate VE
assumed_ors_primary = c(1,1,0.8,0.6,0.6,0.5)

# Make case-control data
run_mod = function(data,exposure,exposure_cat) {
  
  formula_unadj = formula(paste0('case~factor(',exposure,')+strata(stratum_no)'))
  
  mod_unadj = clogit(formula_unadj,data=data)
  
  ests_cis_unadj = cbind(exp(mod_unadj$coefficients),exp(confint(mod_unadj)))
  
  pvals = sapply(1:length(mod_unadj$coefficients),function(x) round(2*pnorm(-abs(mod_unadj$coefficients[x]/sqrt(vcov(mod_unadj)[x,x])),lower.tail=T),2))
  
  return(list(ests_cis_unadj,pvals))
  
}

# Read in data
# 1:1 matched case control data set, exposures are vcat_beforetest (at least one dose) and vcat_combined (1 and 2 doses)
# Variables: stratum_no, case (0,1), vcat_combined (Unvaccinated, Dose 1 0-13 days, Dose 1 14-27 days, Dose 1 >=28 days, Dose 2 0-13 days, Dose 2 >=14 days)
case_control_data_primary = read_csv('ChAdOx1_CaseControlPairs.csv')

case_control_data_primary$vcat_combined_permute = case_control_data_primary$vcat_combined

case_control_data_primary = case_control_data_primary %>% select(stratum_no,case,vcat_combined,vcat_combined_permute) %>% arrange(stratum_no,case)

# Label each pair as discordant
case_control_data_primary$disc_primary = sapply(case_control_data_primary$stratum_no,
                                                function(x) 
                                                  paste(sum(case_control_data_primary$vcat_combined[case_control_data_primary$stratum_no==x]=="Unvaccinated"),
                                                        sum(case_control_data_primary$vcat_combined[case_control_data_primary$stratum_no==x]=="Dose 1 0-13 days"),
                                                        sum(case_control_data_primary$vcat_combined[case_control_data_primary$stratum_no==x]=="Dose 1 >=14 days"),
                                                        sum(case_control_data_primary$vcat_combined[case_control_data_primary$stratum_no==x]=="Dose 2 0-13 days"),
                                                        sum(case_control_data_primary$vcat_combined[case_control_data_primary$stratum_no==x]=="Dose 2 >=14 days"),
                                                        sep='_'))


ves_primary = rep(NA,1000)
ves_uci_primary = rep(NA,1000)
ps_primary = rep(NA,1000)

for (i in 1:1000) {
  
  # Primary analysis
  vcat_labels = c("Unvaccinated","Dose 1 0-13 days","Dose 1 14-27 days", "Dose 1 >=28 days","Dose 2 0-13 days","Dose 2 >=14 days")
  for (p1 in 1:5) {
    for (p2 in (p1+1):6) {
      pair = rep(0,6)
      pair[c(p1,p2)]=1
      
      pair = paste0(pair,collapse='_')
      case_control_data_primary$vcat_combined_permute[case_control_data_primary$disc_primary==pair & case_control_data_primary$case==1] = 
        sapply(1:length(case_control_data_primary$vcat_combined[case_control_data_primary$disc_primary==pair & case_control_data_primary$case==1]),function(x) 
          vcat_labels[c(p1,p2)][1+rbinom(1,1,assumed_ors_primary[p2]/assumed_ors_primary[p1]/(1+assumed_ors_primary[p2]/assumed_ors_primary[p1]))])
      case_control_data_primary$vcat_combined_permute[case_control_data_primary$disc_primary==pair & case_control_data_primary$case==0] = 
        sapply(case_control_data_primary$stratum_no[case_control_data_primary$disc_primary==pair & case_control_data_primary$case==0],
               function(x) setdiff(vcat_labels[c(p1,p2)],case_control_data_primary$vcat_combined_permute[case_control_data_primary$stratum_no==x & case_control_data_primary$case==1])
        )
    }
  }
  
  case_control_data_primary$vcat_combined_permute = relevel(factor(case_control_data_primary$vcat_combined_permute),ref='Unvaccinated')
  
  t=run_mod(case_control_data_primary,'vcat_combined_permute')
  ves_primary_onedose[i] =  unname(t[[1]][1,1])
  ves_uci_primary_onedose[i] =  unname(t[[1]][1,3])
  ps_primary_onedose[i] =  unname(t[[2]][1])
  
  ves_primary_twodose[i] =  unname(t[[1]][4,1])
  ves_uci_primary_twodose[i] =  unname(t[[1]][4,3])
  ps_primary_twodose[i] =  unname(t[[2]][4])
  
}

summary(ves_primary)
summary(ves_uci_primary)

# Estimated power
mean(ps_primary<=0.05)

mean(ves_uci_primary<=0.8)
mean(ves_uci_primary<=0.7)