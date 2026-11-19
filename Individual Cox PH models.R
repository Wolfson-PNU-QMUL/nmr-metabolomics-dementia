## Individual Cox regressions: NMR metabolites and incident dementia (ACD, AD, VaD)
## Model 1: age, age2, sex, age*sex, age2*sex, centre, ethnicity, Townsend index
## Model 2: Model 1 + APOE e4; Model 3: Model 2 + education, smoking, alcohol, obesity, diabetes, hypertension

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

# Metabolite names
metab_names <- nmr_dem %>%
  dplyr::select(bOHbutyrate:Omega_6_pct_PUFA) %>%
  colnames()

## Model 1 ---------------------------------------------------------------------
cox_results_model1 <- data.frame()


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
    
    cox_model <- coxph(Surv(nmr_dem_new[,paste0("time_",i)], nmr_dem_new[,paste0("int_",i)]) ~ scale(nmr_dem_new[,f]) + centre_0 + ethn + townsend + age + age2 + sex + age*sex + age2*sex, data = nmr_dem_new)
    
    cox_results_model1[f,paste0("HR.",i)] <- summary(cox_model)$coefficients[1, "exp(coef)"] # hazard ratio
    cox_results_model1[f,paste0("ci_lower.",i)] <- summary(cox_model)$conf.int[1, "lower .95"]      # Lower CI
    cox_results_model1[f,paste0("ci_upper.",i)] <- summary(cox_model)$conf.int[1, "upper .95"]      # Upper CI
    cox_results_model1[f,paste0("p.",i)] <- summary(cox_model)$coefficients[1, "Pr(>|z|)"]
    
    # Proportional hazards assumption (Schoenfeld residuals)
    cox_model1 <- coxph(Surv(nmr_dem_new[,paste0("time_",i)], nmr_dem_new[,paste0("int_",i)]) ~ scale(nmr_dem_new[,f]) + centre_0 + ethn + townsend + age + sex, data = nmr_dem_new)
    ph <- cox.zph(cox_model1)
    cox_results_model1[f,paste0("ph.",i)] <- ph$table[, "p"][1]
  }
  
  cox_results_model1[ ,paste0("FDR.",i)] <- p.adjust(cox_results_model1[,paste0("p.",i)],method = "BH")# FDR adjust p value
  
}


## Model 2 ---------------------------------------------------------------------
cox_results_model2 <- data.frame()


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
    
    cox_model <- coxph(Surv(nmr_dem_new[,paste0("time_",i)], nmr_dem_new[,paste0("int_",i)]) ~ scale(nmr_dem_new[,f]) + centre_0 + ethn + townsend + age + age2 + sex + age*sex + age2*sex + apoe, data = nmr_dem_new)
    
    cox_results_model2[f,paste0("HR.",i)] <- summary(cox_model)$coefficients[1, "exp(coef)"] # hazard ratio
    cox_results_model2[f,paste0("ci_lower.",i)] <- summary(cox_model)$conf.int[1, "lower .95"]      # Lower CI
    cox_results_model2[f,paste0("ci_upper.",i)] <- summary(cox_model)$conf.int[1, "upper .95"]      # Upper CI
    cox_results_model2[f,paste0("p.",i)] <- summary(cox_model)$coefficients[1, "Pr(>|z|)"]
    
    # Proportional hazards assumption (Schoenfeld residuals)
    cox_model2 <- coxph(Surv(nmr_dem_new[,paste0("time_",i)], nmr_dem_new[,paste0("int_",i)]) ~ scale(nmr_dem_new[,f]) + centre_0 + ethn + townsend + age + sex + apoe, data = nmr_dem_new)
    ph <- cox.zph(cox_model2)
    cox_results_model2[f,paste0("ph.",i)] <- ph$table[, "p"][1]
    
  }
  
  cox_results_model2[ ,paste0("FDR.",i)] <- p.adjust(cox_results_model2[,paste0("p.",i)],method = "BH")# FDR adjust p value
  
}


## Model 3 ---------------------------------------------------------------------
cox_results_model3 <- data.frame()


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
    
    cox_model <- coxph(Surv(nmr_dem_new[,paste0("time_",i)], nmr_dem_new[,paste0("int_",i)]) ~ scale(nmr_dem_new[,f]) + centre_0 + ethn + townsend + age + age2 + sex + age*sex + age2*sex + apoe + edu + smoking + alcohol + obesity + diabetes + hypertension, data = nmr_dem_new)
    
    cox_results_model3[f,paste0("HR.",i)] <- summary(cox_model)$coefficients[1, "exp(coef)"] # hazard ratio
    cox_results_model3[f,paste0("ci_lower.",i)] <- summary(cox_model)$conf.int[1, "lower .95"]      # Lower CI
    cox_results_model3[f,paste0("ci_upper.",i)] <- summary(cox_model)$conf.int[1, "upper .95"]      # Upper CI
    cox_results_model3[f,paste0("p.",i)] <- summary(cox_model)$coefficients[1, "Pr(>|z|)"]
    
    # Proportional hazards assumption (Schoenfeld residuals)
    cox_model3 <- coxph(Surv(nmr_dem_new[,paste0("time_",i)], nmr_dem_new[,paste0("int_",i)]) ~ scale(nmr_dem_new[,f]) + centre_0 + ethn + townsend + age + sex + apoe + edu + smoking + alcohol + obesity + diabetes + hypertension, data = nmr_dem_new)
    ph <- cox.zph(cox_model3)
    cox_results_model3[f,paste0("ph.",i)] <- ph$table[, "p"][1]
    
  }
  
  cox_results_model3[ ,paste0("FDR.",i)] <- p.adjust(cox_results_model3[,paste0("p.",i)],method = "BH")# FDR adjust p value
  
}


## Save results (summary statistics only) --------------------------------------
write.csv(cox_results_model1, file.path(out_dir, "res_cox_model1.csv"))
write.csv(cox_results_model2, file.path(out_dir, "res_cox_model2.csv"))
write.csv(cox_results_model3, file.path(out_dir, "res_cox_model3.csv"))