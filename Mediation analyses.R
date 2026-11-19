## Mediation analyses: key metabolites as mediators of the associations between modifiable risk factors and incident dementia
## Product-of-coefficients approach (Baron and Kenny):
##   path a: linear regression of metabolite on risk factor
##   path b: Cox regression of dementia on metabolite, adjusted for the risk factor
##   path c': direct effect of the risk factor, from the same Cox model as path b
##   path c: total effect of the risk factor, from a Cox model without the metabolite
## Indirect effect = a * b; proportion mediated = indirect effect / total effect; 95% CIs from 1,000 bootstrap samples
## Indirect effects are calculated only when both path a and path b are significant
## Model 1: age, age2, sex, age*sex, age2*sex, assessment centre, ethnicity, Townsend index (as in the Cox models);
## Model 2: Model 1 + APOE e4; Model 3: Model 2 + all remaining modifiable risk factors

## Set directories---------------------------------------------------------------
data_dir <- "path/to/data_derived/"   # folder containing input data
out_dir  <- "path/to/output/"         # folder for results

## Load packages-----------------------------------------------------------------
library(data.table)
library(survival)
library(dplyr)
library(boot)

### Load data--------------------------------------------------------------------
# nmr_dem.csv: analysis dataset (n = 269,284) with NMR metabolites, dementia outcomes and covariates
nmr_dem <- fread(file.path(data_dir, "nmr_dem.csv"))
nmr_dem <- as.data.frame(nmr_dem)

### Prepare variables------------------------------------------------------------
# Modifiable risk factors (0/1); smoking and alcohol: 0 = never, 1 = former or current
nmr_dem$edu <- ifelse(nmr_dem$edu=="Secondary Education",1,0)
nmr_dem$smoking <- ifelse(nmr_dem$smoking=="Never",0,1)
nmr_dem$alcohol <- ifelse(nmr_dem$alcohol=="Never",0,1)
nmr_dem$obesity <- ifelse(nmr_dem$obesity=="Yes",1,0)
nmr_dem$diabetes <- ifelse(nmr_dem$diabetes=="Yes",1,0)
nmr_dem$hypertension <- ifelse(nmr_dem$hypertension=="Yes",1,0)

# Covariates
nmr_dem$townsend <- as.numeric(nmr_dem$townsend)
nmr_dem$centre_0 <- as.factor(nmr_dem$centre_0)
# Define the mapping for codes to countries
nmr_dem$centre_0 <- ifelse(nmr_dem$centre_0 %in% c(11012, 11021, 11011, 11008,
                                                   11024, 11020, 11018,
                                                   11010, 11016, 11001, 11017, 11009,
                                                   11013, 11002, 11007, 11014, 11006), "England",
                           ifelse(nmr_dem$centre_0 %in% c(11003, 11022, 11023), "Wales",
                                  ifelse(nmr_dem$centre_0 %in% c(11004, 11005), "Scotland", NA)))

# Rank-based inverse normal transformation and scaling of metabolites
metab_names <- nmr_dem%>%
  dplyr::select(bOHbutyrate:Omega_6_pct_PUFA) %>%
  colnames()

for(i in metab_names){
  nmr_dem[,i] <- qnorm((rank(nmr_dem[,i], na.last='keep')-0.5)/sum(!is.na(nmr_dem[,i])))
  nmr_dem[,i] <- scale(nmr_dem[,i])
}

### Key metabolites from the LASSO-Cox models-------------------------------------
# Top 15 (ACD, AD) or 10 (VaD) metabolites ranked by mean regularised coefficient across 100 iterations
get_top_metab <- function(res, n_top) {
  cox_lasso_coef_m1 <- lapply(res, function(x) x$cox_lasso_coef_m1)
  coef_df <- do.call(rbind, cox_lasso_coef_m1)
  df_summary <- coef_df %>%
    dplyr::group_by(metab) %>%
    dplyr::summarise(mean_s0 = sum(s0)/100) %>%
    dplyr::arrange(desc(abs(mean_s0)))
  return(df_summary[1:n_top,]$metab)
}

key_metab <- unique(c(get_top_metab(readRDS(file.path(out_dir, "lasso_dem_full_3splits.rds")), 15),
                      get_top_metab(readRDS(file.path(out_dir, "lasso_AD_full_3splits.rds")), 15),
                      get_top_metab(readRDS(file.path(out_dir, "lasso_VaD_full_3splits.rds")), 10)))

### Exposures, outcomes and models-----------------------------------------------
rf_names <- c("edu", "smoking", "alcohol", "obesity", "diabetes", "hypertension")
outcomes <- c("dem", "AD", "VaD")

cov_m1 <- "age + age2 + sex + age*sex + age2*sex + centre_0 + ethn + townsend"
cov_m2 <- paste(cov_m1, "+ apoe")
# Model 3: Model 2 + all modifiable risk factors other than the exposure
get_covs <- function(model, rf) {
  if (model == "model1") return(cov_m1)
  if (model == "model2") return(cov_m2)
  if (model == "model3") return(paste(cov_m2, "+", paste(setdiff(rf_names, rf), collapse = " + ")))
}

### Paths a, b, c' and c-----------------------------------------------------------
results_list <- list()

for (k in outcomes){
  print(k)
  col_name <- paste0("int_", k)
  # For AD and VaD, exclude participants with other types of dementia
  nmr_dem_new <- nmr_dem[!(nmr_dem$int_dem == 1 &nmr_dem[[col_name]] == 0), ]
  nmr_dem_new$time_tmp <- nmr_dem_new[,paste0("time_",k)]
  nmr_dem_new$int_tmp <- nmr_dem_new[,paste0("int_",k)]
  
  for (m in c("model1", "model2", "model3")){
    # APOE e4 as an additional exposure in Model 1 only (adjusted for in Models 2 and 3)
    exposures <- if (m == "model1") c("apoe", rf_names) else rf_names
    for (rf in exposures){
      covs <- get_covs(m, rf)
      nmr_dem_new$rf_tmp <- nmr_dem_new[,rf]   # exposure (risk factor)
      
      # Path c (total effect: exposure, without the metabolite)
      cox_model_c <- coxph(as.formula(paste("Surv(time_tmp, int_tmp) ~ scale(rf_tmp) +", covs)), data = nmr_dem_new)
      coef_c <- summary(cox_model_c)$coefficients
      
      for (i in key_metab){
        print(paste(m, rf, i))
        nmr_dem_new$met_tmp <- nmr_dem_new[,i]   # mediator (metabolite)
        
        # Path a (linear association between exposure and mediator)
        linear_a <- lm(as.formula(paste("scale(met_tmp) ~ scale(rf_tmp) +", covs)), data = nmr_dem_new)
        coef_a <- summary(linear_a)$coefficients
        conf_a <- confint(linear_a, level = 0.95)
        
        # Path b (mediator) and path c' (direct effect of exposure), from the same Cox model
        cox_model_b <- coxph(as.formula(paste("Surv(time_tmp, int_tmp) ~ scale(met_tmp) + scale(rf_tmp) +", covs)), data = nmr_dem_new)
        coef_b <- summary(cox_model_b)$coefficients
        conf_b <- confint(cox_model_b, level = 0.95)
        
        results_list[[paste(k, m, rf, i, sep = "_")]] <- data.frame(
          Outcome = k, Model = m, Risk_factor = rf, Metabolite = i,
          a = coef_a["scale(rf_tmp)", 1], a_lower = conf_a["scale(rf_tmp)", 1], a_upper = conf_a["scale(rf_tmp)", 2], p_a = coef_a["scale(rf_tmp)", 4],
          b = coef_b["scale(met_tmp)", 1], b_lower = conf_b["scale(met_tmp)", 1], b_upper = conf_b["scale(met_tmp)", 2], p_b = coef_b["scale(met_tmp)", 5],
          direct = coef_b["scale(rf_tmp)", 1], direct_lower = conf_b["scale(rf_tmp)", 1], direct_upper = conf_b["scale(rf_tmp)", 2], p_direct = coef_b["scale(rf_tmp)", 5],
          total = coef_c["scale(rf_tmp)", 1], p_total = coef_c["scale(rf_tmp)", 5]
        )
      }
    }
  }
}

# Combine all results into a single data frame
merged_results_paths <- do.call(rbind, results_list)
write.csv(merged_results_paths, file.path(out_dir, "res_mediation_paths.csv"), row.names=F)

### Indirect effect and proportion mediated (bootstrap)----------------------------
# Only pathways where both path a and path b are significant
paths_sig <- merged_results_paths[merged_results_paths$p_a < 0.05 & merged_results_paths$p_b < 0.05, ]

# Bootstrap function
mediation_function <- function(data_used, idx, covs) {
  data_temp <- data_used[idx, ]
  a <- lm(as.formula(paste("scale(met_tmp) ~ scale(rf_tmp) +", covs)), data = data_temp)$coefficients["scale(rf_tmp)"]
  b <- coxph(as.formula(paste("Surv(time_tmp, int_tmp) ~ scale(met_tmp) + scale(rf_tmp) +", covs)), data = data_temp)$coefficients["scale(met_tmp)"]
  c <- coxph(as.formula(paste("Surv(time_tmp, int_tmp) ~ scale(rf_tmp) +", covs)), data = data_temp)$coefficients["scale(rf_tmp)"]
  
  # Calculate indirect effect and proportion of mediation effect
  ACME <- a * b # Indirect effect
  prop.med <- ACME/c # Proportion of mediation effect
  return(c(ACME, prop.med))
}

res_mediation <- list()

for (r in seq_len(nrow(paths_sig))){
  k <- paths_sig$Outcome[r]; m <- paths_sig$Model[r]; rf <- paths_sig$Risk_factor[r]; i <- paths_sig$Metabolite[r]
  print(paste(k, m, rf, i))
  covs <- get_covs(m, rf)
  
  col_name <- paste0("int_", k)
  nmr_dem_new <- nmr_dem[!(nmr_dem$int_dem == 1 &nmr_dem[[col_name]] == 0), ]
  nmr_dem_new$time_tmp <- nmr_dem_new[,paste0("time_",k)]
  nmr_dem_new$int_tmp <- nmr_dem_new[,paste0("int_",k)]
  nmr_dem_new$rf_tmp <- nmr_dem_new[,rf]
  nmr_dem_new$met_tmp <- nmr_dem_new[,i]
  
  # Run bootstrapping
  set.seed(777)
  boot_mediation <- boot(nmr_dem_new, mediation_function, R = 1000, covs = covs)
  
  # Confidence intervals for indirect effect (ACME) and proportion mediated
  boot_ci_ACME <- boot.ci(boot.out = boot_mediation, index = 1, type = c("norm", "perc"))
  boot_ci_prop.med <- boot.ci(boot.out = boot_mediation, index = 2, type = c("norm", "perc"))
  
  res_mediation[[paste(k, m, rf, i, sep = "_")]] <- data.frame(
    Outcome = k, Model = m, Risk_factor = rf, Metabolite = i,
    ACME = boot_mediation$t0[1],
    SE_ACME = sd(boot_mediation$t[, 1]),
    ciLower_norm_ACME = boot_ci_ACME$normal[2],
    ciUpper_norm_ACME = boot_ci_ACME$normal[3],
    ciLower_perc_ACME = boot_ci_ACME$percent[4],
    ciUpper_perc_ACME = boot_ci_ACME$percent[5],
    prop.med = boot_mediation$t0[2],
    SE_prop.med = sd(boot_mediation$t[, 2]),
    ciLower_norm_prop.med = boot_ci_prop.med$normal[2],
    ciUpper_norm_prop.med = boot_ci_prop.med$normal[3],
    ciLower_perc_prop.med = boot_ci_prop.med$percent[4],
    ciUpper_perc_prop.med = boot_ci_prop.med$percent[5]
  )
}

# Combine all results into a single data frame
merged_results_mediation <- do.call(rbind, res_mediation)
write.csv(merged_results_mediation, file.path(out_dir, "res_mediation_bootstrap.csv"), row.names = F)