## Associations of key metabolites and metabolomic risk scores (MetRS) with global and regional brain imaging measures
## Corresponds to Methods (Associations with Brain Imaging-Derived Phenotypes), Figure 4 and Figure 5
## MetRS: weighted sum of the top 15 (ACD, AD) or 10 (VaD) metabolites, using mean regularised coefficients
## from the LASSO-Cox models (LASSO_Cox_MetRS_*.R)
## Linear regression with three levels of adjustment (consistent with the individual Cox models), each additionally
## adjusted for imaging centre, total intracranial volume and time between baseline and imaging visit
## Model 1: age, age2, sex, age*sex, age2*sex, assessment centre, ethnicity, Townsend index
## Model 2: Model 1 + APOE e4; Model 3: Model 2 + education, smoking, alcohol, obesity, diabetes, hypertension
## Figures 4 and 5 show Model 2 results

## Set directories---------------------------------------------------------------
data_dir <- "path/to/data_derived/"   # folder containing input data
out_dir  <- "path/to/output/"         # folder for results

## Load packages-----------------------------------------------------------------
library(tibble)
library(dplyr)
library(tidyr)
library(data.table)

### Load data--------------------------------------------------------------------
# nmr_img.csv: participants with NMR metabolites and brain imaging data (n = 28,068), including covariates (KNN imputed),
# brain age, global and regional (Desikan-Killiany and ASEG) imaging measures, TIV, imaging centre (centre_2)
# and time between baseline and imaging visit (time_int)
nmr_img <- fread(file.path(data_dir, "nmr_img.csv"))
nmr_img <- as.data.frame(nmr_img)

# Names of metabolites and imaging measures
nmr_info <- read.csv(file.path(data_dir, "nmr_info_newgroup.csv"))
metab_names <- nmr_info$Biomarker
vol_names_global <- tibble(nmr_img) %>%
  dplyr::select(brainage,brainage_Delta, "Cortex_L": "Deep WMHs")%>%
  colnames()
vol_names_regional <- tibble(nmr_img) %>%
  dplyr::select("x3rd-ventricle": "rh_transversetemporal")%>%
  colnames()

# Rank-based inverse normal transformation of metabolites
for(i in metab_names){
  nmr_img[,i] <- qnorm((rank(nmr_img[,i], na.last='keep')-0.5)/sum(!is.na(nmr_img[,i])))
}

### Metabolomic risk scores (MetRS)-----------------------------------------------
# Results of LASSO-Cox models (100 iterations) for each dementia outcome
res_dem <- readRDS(file.path(out_dir, "lasso_dem_full_3splits.rds"))
res_AD <- readRDS(file.path(out_dir, "lasso_AD_full_3splits.rds"))
res_VaD <- readRDS(file.path(out_dir, "lasso_VaD_full_3splits.rds"))

# Function to get the top metabolites and calculate the MetRS
calc_metrs <- function(res, n_top, score_name) {
  # Mean regularised coefficient across 100 iterations (0 when not selected)
  cox_lasso_coef_m1 <- lapply(res, function(x) x$cox_lasso_coef_m1)
  coef_df <- do.call(rbind, cox_lasso_coef_m1)
  
  df_summary <- coef_df %>%
    dplyr::group_by(metab) %>%
    dplyr::summarise(mean_s0 = sum(s0)/100) %>%
    dplyr::arrange(desc(abs(mean_s0)))
  
  df_summary_top <- df_summary[1:n_top,]
  
  # reshape individual metabolomic data
  img_df <- nmr_img[,c("eid",df_summary_top$metab)] #only keep eid and metabolites data
  img_df_reshaped <- img_df %>% gather(key = "metab", value = "value",-eid) #reshape individual metabolomic data to 3 columns("eid",metab","value")
  
  # merge reshaped individual metabolomic data with mean coefficients
  nmr_df_coef <- inner_join(img_df_reshaped, df_summary_top, by = "metab")
  
  # In each individual, calculate the sum of metabolite values multiplied by the coefficients of corresponding metabolites
  mrs <- nmr_df_coef %>%
    dplyr::group_by(eid) %>%
    dplyr::summarise(score = sum(value * mean_s0, na.rm = TRUE))
  colnames(mrs) <- c("eid", score_name)
  
  return(list(top_metab = df_summary_top$metab, mrs = mrs))
}

metrs_dem <- calc_metrs(res_dem, 15, "MetRS_ACD")
metrs_AD <- calc_metrs(res_AD, 15, "MetRS_AD")
metrs_VaD <- calc_metrs(res_VaD, 10, "MetRS_VaD")

topmetab_dem <- metrs_dem$top_metab
topmetab_AD <- metrs_AD$top_metab
topmetab_VaD <- metrs_VaD$top_metab

nmr_mrs_img <- nmr_img %>%
  inner_join(., metrs_dem$mrs,by="eid") %>%
  inner_join(., metrs_AD$mrs,by="eid") %>%
  inner_join(., metrs_VaD$mrs,by="eid")

### Prepare covariates------------------------------------------------------------
nmr_mrs_img$townsend <- as.numeric(nmr_mrs_img$townsend)
nmr_mrs_img$centre_0 <- as.factor(nmr_mrs_img$centre_0)
nmr_mrs_img$centre_2 <- as.factor(nmr_mrs_img$centre_2)
# Define the mapping for codes to countries
nmr_mrs_img$centre_0 <- ifelse(nmr_mrs_img$centre_0 %in% c(11012, 11021, 11011, 11008,
                                                           11024, 11020, 11018,
                                                           11010, 11016, 11001, 11017, 11009,
                                                           11013, 11002, 11007, 11014, 11006), "England",
                               ifelse(nmr_mrs_img$centre_0 %in% c(11003, 11022, 11023), "Wales",
                                      ifelse(nmr_mrs_img$centre_0 %in% c(11004, 11005), "Scotland", NA)))

nmr_mrs_img <- as.data.frame(nmr_mrs_img)

# Key metabolites from the three dementia-specific signatures, and MetRS
metab_name_nomrs <- unique(c(topmetab_dem, topmetab_AD, topmetab_VaD))
metab_all_names <- c(metab_name_nomrs, "MetRS_ACD", "MetRS_AD", "MetRS_VaD")

# Covariates for the three models
cov_m1 <- "TIV + time_int + centre_2 + age + age2 + sex + age*sex + age2*sex + centre_0 + ethn + townsend"
cov_m2 <- paste(cov_m1, "+ apoe")
cov_m3 <- paste(cov_m2, "+ edu + smoking + alcohol + obesity + diabetes + hypertension")
cov_models <- list(model1 = cov_m1, model2 = cov_m2, model3 = cov_m3)

# Function to run linear regressions of each imaging measure on each metabolite / MetRS
run_lm <- function(img_names, predictor_names, covs) {
  lm_results <- data.frame()
  for (i in img_names){
    print(i)
    for (f in predictor_names){
      nmr_mrs_img$y_tmp <- nmr_mrs_img[,i]   # imaging measure
      nmr_mrs_img$x_tmp <- nmr_mrs_img[,f]   # metabolite or MetRS
      lm_model <- lm(as.formula(paste("scale(y_tmp) ~ scale(x_tmp) +", covs)), data=nmr_mrs_img)
      lm_results[f,paste0("beta.",i)]<-summary(lm_model)$coefficients[2,1]# beta value
      lm_results[f,paste0("p.",i)]<-summary(lm_model)$coefficients[2,4]# p value
    }
    lm_results[ ,paste0("FDR.",i)] <- p.adjust(lm_results[,paste0("p.",i)],method = "BH")# FDR adjust p value
  }
  return(lm_results)
}

# 1) Global neuroimaging----------------------------------------------------------
lm_global <- list()
for (m in names(cov_models)){
  lm_global[[m]] <- run_lm(vol_names_global, metab_all_names, cov_models[[m]])
  write.csv(lm_global[[m]], file.path(out_dir, paste0("lm_results_global_img_", m, ".csv")))
}

# 2) Regional neuroimaging--------------------------------------------------------
lm_regional <- list()
for (m in names(cov_models)){
  lm_regional[[m]] <- run_lm(vol_names_regional, metab_all_names, cov_models[[m]])
  write.csv(lm_regional[[m]], file.path(out_dir, paste0("lm_results_regional_img_", m, ".csv")))
}