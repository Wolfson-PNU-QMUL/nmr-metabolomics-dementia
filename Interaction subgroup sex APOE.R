## Interaction and subgroup analyses by sex and APOE e4 carrier status: NMR metabolites and incident dementia (ACD, AD, VaD)
## Corresponds to Methods (Interaction and subgroup analysis) and Supplementary Tables 6-8
## Part 1: metabolite x sex, metabolite x APOE e4 and metabolite x sex x APOE e4 interactions
## Part 2: subgroup analyses stratified by sex, APOE e4 carrier status and their combinations


## Set directories---------------------------------------------------------------
data_dir <- "path/to/data_derived/"   # folder containing input data
out_dir  <- "path/to/output/"         # folder for results

## Load packages-----------------------------------------------------------------
library(dplyr)
library(survival)
library(data.table)

### Load data--------------------------------------------------------------------
# nmr_dem.csv: analysis dataset (n = 269,284) with NMR metabolites, dementia outcomes and covariates
# Covariates were imputed using KNN (k = 10). Participants with prevalent or early-onset (<60 years) dementia,
# >20% missing metabolite values, missing APOE genotype, or diagnoses of both AD and VaD were excluded
nmr_dem <- fread(file.path(data_dir, "nmr_dem.csv"))
nmr_dem <- as.data.frame(nmr_dem)

### Recode covariates------------------------------------------------------------
nmr_dem$townsend <- as.numeric(nmr_dem$townsend)
nmr_dem$centre_0 <- as.factor(nmr_dem$centre_0)
# Define the mapping for codes to countries
nmr_dem$centre_0 <- ifelse(nmr_dem$centre_0 %in% c(11012, 11021, 11011, 11008,
                                                   11024, 11020, 11018,
                                                   11010, 11016, 11001, 11017, 11009,
                                                   11013, 11002, 11007, 11014, 11006), "England",
                           ifelse(nmr_dem$centre_0 %in% c(11003, 11022, 11023), "Wales",
                                  ifelse(nmr_dem$centre_0 %in% c(11004, 11005), "Scotland", NA)))

nmr_dem <- as.data.frame(nmr_dem)

# APOE e4 carrier status (1 or 2 alleles = carrier)
nmr_dem$apoe[nmr_dem$apoe != 0] <- 1

# Metabolite names
metab_names <- nmr_dem %>%
  dplyr::select(bOHbutyrate:Omega_6_pct_PUFA) %>%
  colnames()

## Part 1: Interaction models----------------------------------------------------
cox_results <- data.frame()


for (i in c("dem","AD","VaD")){
  print(i)
  col_name <- paste0("int_", i)
  # For AD and VaD, exclude participants with other types of dementia
  nmr_dem_new <- nmr_dem[!(nmr_dem$int_dem == 1 &nmr_dem[[col_name]] == 0), ]
  
  # Rank-based inverse normal transformation of metabolites
  for(j in metab_names){
    nmr_dem_new[,j] <- qnorm((rank(nmr_dem_new[,j], na.last='keep')-0.5)/sum(!is.na(nmr_dem_new[,j])))
  }
  
  for (f in metab_names){
    print(f)
    
    cox_model_sex <- coxph(Surv(nmr_dem_new[,paste0("time_",i)], nmr_dem_new[,paste0("int_",i)]) ~ scale(nmr_dem_new[,f])*sex + centre_0 + ethn + townsend + age + age2 + apoe, data = nmr_dem_new)
    cox_model_apoe <- coxph(Surv(nmr_dem_new[,paste0("time_",i)], nmr_dem_new[,paste0("int_",i)]) ~ scale(nmr_dem_new[,f])*apoe + centre_0 + ethn + townsend + age + age2 + sex + age * sex + age2*sex, data = nmr_dem_new)
    cox_model_apoesex <- coxph(Surv(nmr_dem_new[,paste0("time_",i)], nmr_dem_new[,paste0("int_",i)]) ~ scale(nmr_dem_new[,f])*sex*apoe + centre_0 + ethn + townsend + age + age2, data = nmr_dem_new)
    
    # Row 12: metabolite x sex interaction term
    cox_results[f,paste0("HR_sex.",i)] <- summary(cox_model_sex)$coefficients[12, "exp(coef)"] # hazard ratio
    cox_results[f,paste0("cilower_sex.",i)] <- summary(cox_model_sex)$conf.int[12, "lower .95"]      # Lower CI
    cox_results[f,paste0("ciupper_sex.",i)] <- summary(cox_model_sex)$conf.int[12, "upper .95"]      # Upper CI
    cox_results[f,paste0("p_sex.",i)] <- summary(cox_model_sex)$coefficients[12, "Pr(>|z|)"]
    
    # Row 12: metabolite x APOE e4 interaction term
    cox_results[f,paste0("HR_apoe.",i)] <- summary(cox_model_apoe)$coefficients[12, "exp(coef)"] # hazard ratio
    cox_results[f,paste0("cilower_apoe.",i)] <- summary(cox_model_apoe)$conf.int[12, "lower .95"]      # Lower CI
    cox_results[f,paste0("ciupper_apoe.",i)] <- summary(cox_model_apoe)$conf.int[12, "upper .95"]      # Upper CI
    cox_results[f,paste0("p_apoe.",i)] <- summary(cox_model_apoe)$coefficients[12, "Pr(>|z|)"]
    
    # Row 15: metabolite x sex x APOE e4 interaction term
    cox_results[f,paste0("HR_apoesexmetab.",i)] <- summary(cox_model_apoesex)$coefficients[15, "exp(coef)"] # hazard ratio
    cox_results[f,paste0("cilower_apoesexmetab.",i)] <- summary(cox_model_apoesex)$conf.int[15, "lower .95"]      # Lower CI
    cox_results[f,paste0("ciupper_apoesexmetab.",i)] <- summary(cox_model_apoesex)$conf.int[15, "upper .95"]      # Upper CI
    cox_results[f,paste0("p_apoesexmetab.",i)] <- summary(cox_model_apoesex)$coefficients[15, "Pr(>|z|)"]
    
  }
  
  cox_results[ ,paste0("FDR_sex.",i)] <- p.adjust(cox_results[,paste0("p_sex.",i)],method = "BH") # FDR adjust p value
  cox_results[ ,paste0("FDR_apoe.",i)] <- p.adjust(cox_results[,paste0("p_apoe.",i)],method = "BH") # FDR adjust p value
  cox_results[ ,paste0("FDR_apoesexmetab.",i)] <- p.adjust(cox_results[,paste0("p_apoesexmetab.",i)],method = "BH") # FDR adjust p value
  
}


## Part 2: Subgroup analyses-----------------------------------------------------
cox_results_subgroups <- data.frame()

for (i in c("dem","AD","VaD")){
  print(i)
  col_name <- paste0("int_", i)
  # For AD and VaD, exclude participants with other types of dementia
  nmr_dem_new <- nmr_dem[!(nmr_dem$int_dem == 1 &nmr_dem[[col_name]] == 0), ]
  
  # Subgroups by sex, APOE e4 carrier status and their combinations
  nmr_dem_male <- nmr_dem_new[nmr_dem_new$sex== "Male",]
  nmr_dem_female <- nmr_dem_new[nmr_dem_new$sex== "Female",]
  nmr_dem_noapoe <- nmr_dem_new[nmr_dem_new$apoe == 0,]
  nmr_dem_apoe <- nmr_dem_new[nmr_dem_new$apoe != 0,]
  nmr_dem_male_noapoe <- nmr_dem_new[nmr_dem_new$sex=="Male" & nmr_dem_new$apoe == 0, ]
  nmr_dem_female_noapoe <- nmr_dem_new[nmr_dem_new$sex=="Female" & nmr_dem_new$apoe == 0, ]
  nmr_dem_male_apoe <- nmr_dem_new[nmr_dem_new$sex=="Male" & nmr_dem_new$apoe != 0, ]
  nmr_dem_female_apoe <- nmr_dem_new[nmr_dem_new$sex=="Female" & nmr_dem_new$apoe != 0, ]
  
  df_list <- list(nmr_dem_female, nmr_dem_male,nmr_dem_noapoe,nmr_dem_apoe,nmr_dem_female_noapoe,nmr_dem_male_noapoe,nmr_dem_female_apoe,nmr_dem_male_apoe)
  
  # Rank-based inverse normal transformation of metabolites within each subgroup
  for (k in 1:length(df_list)) {
    for (j in metab_names) {
      df_list[[k]][[j]] <- qnorm((rank(df_list[[k]][[j]], na.last = 'keep') - 0.5) / sum(!is.na(df_list[[k]][[j]])))
    }
  }
  nmr_dem_female <- df_list[[1]]
  nmr_dem_male <- df_list[[2]]
  nmr_dem_noapoe <- df_list[[3]]
  nmr_dem_apoe <- df_list[[4]]
  nmr_dem_female_noapoe <- df_list[[5]]
  nmr_dem_male_noapoe <- df_list[[6]]
  nmr_dem_female_apoe <- df_list[[7]]
  nmr_dem_male_apoe <- df_list[[8]]
  
  for (f in metab_names){
    print(f)
    
    # Sex subgroups: adjusted for APOE e4; APOE subgroups: adjusted for sex, age*sex, age2*sex
    cox_model_male <- coxph(Surv(nmr_dem_male[,paste0("time_",i)], nmr_dem_male[,paste0("int_",i)]) ~ scale(nmr_dem_male[,f]) + centre_0 + ethn + townsend + age + age2 + apoe, data = nmr_dem_male)
    cox_model_female <- coxph(Surv(nmr_dem_female[,paste0("time_",i)], nmr_dem_female[,paste0("int_",i)]) ~ scale(nmr_dem_female[,f]) + centre_0 + ethn + townsend + age + age2 + apoe, data = nmr_dem_female)
    
    cox_model_apoe <- coxph(Surv(nmr_dem_apoe[,paste0("time_",i)], nmr_dem_apoe[,paste0("int_",i)]) ~ scale(nmr_dem_apoe[,f]) + centre_0 + ethn + townsend + age + age2 + age *sex + age2*sex, data = nmr_dem_apoe)
    cox_model_noapoe <- coxph(Surv(nmr_dem_noapoe[,paste0("time_",i)], nmr_dem_noapoe[,paste0("int_",i)]) ~ scale(nmr_dem_noapoe[,f]) + centre_0 + ethn + townsend + age + age2 + age *sex + age2*sex, data = nmr_dem_noapoe)
    
    cox_model_male_noapoe <- coxph(Surv(nmr_dem_male_noapoe[,paste0("time_",i)], nmr_dem_male_noapoe[,paste0("int_",i)]) ~ scale(nmr_dem_male_noapoe[,f]) + centre_0 + ethn + townsend + age + age2, data = nmr_dem_male_noapoe)
    cox_model_female_noapoe <- coxph(Surv(nmr_dem_female_noapoe[,paste0("time_",i)], nmr_dem_female_noapoe[,paste0("int_",i)]) ~ scale(nmr_dem_female_noapoe[,f]) + centre_0 + ethn + townsend + age + age2, data = nmr_dem_female_noapoe)
    cox_model_male_apoe <- coxph(Surv(nmr_dem_male_apoe[,paste0("time_",i)], nmr_dem_male_apoe[,paste0("int_",i)]) ~ scale(nmr_dem_male_apoe[,f]) + centre_0 + ethn + townsend + age + age2, data = nmr_dem_male_apoe)
    cox_model_female_apoe <- coxph(Surv(nmr_dem_female_apoe[,paste0("time_",i)], nmr_dem_female_apoe[,paste0("int_",i)]) ~ scale(nmr_dem_female_apoe[,f]) + centre_0 + ethn + townsend + age + age2, data = nmr_dem_female_apoe)
    
    
    cox_results_subgroups[f,paste0("HR_female.",i)] <- summary(cox_model_female)$coefficients[1, "exp(coef)"] # hazard ratio
    cox_results_subgroups[f,paste0("cilower_female.",i)] <- summary(cox_model_female)$conf.int[1, "lower .95"]      # Lower CI
    cox_results_subgroups[f,paste0("ciupper_female.",i)] <- summary(cox_model_female)$conf.int[1, "upper .95"]      # Upper CI
    cox_results_subgroups[f,paste0("p_female.",i)] <- summary(cox_model_female)$coefficients[1, "Pr(>|z|)"]
    cox_results_subgroups[f,paste0("HR_male.",i)] <- summary(cox_model_male)$coefficients[1, "exp(coef)"] # hazard ratio
    cox_results_subgroups[f,paste0("cilower_male.",i)] <- summary(cox_model_male)$conf.int[1, "lower .95"]      # Lower CI
    cox_results_subgroups[f,paste0("ciupper_male.",i)] <- summary(cox_model_male)$conf.int[1, "upper .95"]      # Upper CI
    cox_results_subgroups[f,paste0("p_male.",i)] <- summary(cox_model_male)$coefficients[1, "Pr(>|z|)"]
    
    
    cox_results_subgroups[f,paste0("HR_noapoe.",i)] <- summary(cox_model_noapoe)$coefficients[1, "exp(coef)"] # hazard ratio
    cox_results_subgroups[f,paste0("cilower_noapoe.",i)] <- summary(cox_model_noapoe)$conf.int[1, "lower .95"]      # Lower CI
    cox_results_subgroups[f,paste0("ciupper_noapoe.",i)] <- summary(cox_model_noapoe)$conf.int[1, "upper .95"]      # Upper CI
    cox_results_subgroups[f,paste0("p_noapoe.",i)] <- summary(cox_model_noapoe)$coefficients[1, "Pr(>|z|)"]
    cox_results_subgroups[f,paste0("HR_apoe.",i)] <- summary(cox_model_apoe)$coefficients[1, "exp(coef)"] # hazard ratio
    cox_results_subgroups[f,paste0("cilower_apoe.",i)] <- summary(cox_model_apoe)$conf.int[1, "lower .95"]      # Lower CI
    cox_results_subgroups[f,paste0("ciupper_apoe.",i)] <- summary(cox_model_apoe)$conf.int[1, "upper .95"]      # Upper CI
    cox_results_subgroups[f,paste0("p_apoe.",i)] <- summary(cox_model_apoe)$coefficients[1, "Pr(>|z|)"]
    
    cox_results_subgroups[f,paste0("HR_female_noapoe.",i)] <- summary(cox_model_female_noapoe)$coefficients[1, "exp(coef)"] # hazard ratio
    cox_results_subgroups[f,paste0("cilower_female_noapoe.",i)] <- summary(cox_model_female_noapoe)$conf.int[1, "lower .95"]      # Lower CI
    cox_results_subgroups[f,paste0("ciupper_female_noapoe.",i)] <- summary(cox_model_female_noapoe)$conf.int[1, "upper .95"]      # Upper CI
    cox_results_subgroups[f,paste0("p_female_noapoe.",i)] <- summary(cox_model_female_noapoe)$coefficients[1, "Pr(>|z|)"]
    cox_results_subgroups[f,paste0("HR_male_noapoe.",i)] <- summary(cox_model_male_noapoe)$coefficients[1, "exp(coef)"] # hazard ratio
    cox_results_subgroups[f,paste0("cilower_male_noapoe.",i)] <- summary(cox_model_male_noapoe)$conf.int[1, "lower .95"]      # Lower CI
    cox_results_subgroups[f,paste0("ciupper_male_noapoe.",i)] <- summary(cox_model_male_noapoe)$conf.int[1, "upper .95"]      # Upper CI
    cox_results_subgroups[f,paste0("p_male_noapoe.",i)] <- summary(cox_model_male_noapoe)$coefficients[1, "Pr(>|z|)"]
    cox_results_subgroups[f,paste0("HR_female_apoe.",i)] <- summary(cox_model_female_apoe)$coefficients[1, "exp(coef)"] # hazard ratio
    cox_results_subgroups[f,paste0("cilower_female_apoe.",i)] <- summary(cox_model_female_apoe)$conf.int[1, "lower .95"]      # Lower CI
    cox_results_subgroups[f,paste0("ciupper_female_apoe.",i)] <- summary(cox_model_female_apoe)$conf.int[1, "upper .95"]      # Upper CI
    cox_results_subgroups[f,paste0("p_female_apoe.",i)] <- summary(cox_model_female_apoe)$coefficients[1, "Pr(>|z|)"]
    cox_results_subgroups[f,paste0("HR_male_apoe.",i)] <- summary(cox_model_male_apoe)$coefficients[1, "exp(coef)"] # hazard ratio
    cox_results_subgroups[f,paste0("cilower_male_apoe.",i)] <- summary(cox_model_male_apoe)$conf.int[1, "lower .95"]      # Lower CI
    cox_results_subgroups[f,paste0("ciupper_male_apoe.",i)] <- summary(cox_model_male_apoe)$conf.int[1, "upper .95"]      # Upper CI
    cox_results_subgroups[f,paste0("p_male_apoe.",i)] <- summary(cox_model_male_apoe)$coefficients[1, "Pr(>|z|)"]
    
    
    
  }
  
  
}


## Save results (summary statistics only) --------------------------------------
write.csv(cox_results, file.path(out_dir, "cox_results_apoesex.csv"))
write.csv(cox_results_subgroups, file.path(out_dir, "cox_results_subgroups.csv"))
