## LASSO-regularised Cox regression for all-cause dementia (ACD): metabolite selection, metabolomic risk score (MetRS) and predictive performance
## Participants with neuroimaging data are excluded (reserved for imaging analyses)
## Data split: training 50%, validation 25%, test 25%, stratified by outcome; metabolites inverse rank normalised within each subset
## Part 1: 10-fold cross-validation (C-index) in the training set, repeated over 100 seeds
## Part 2: number of metabolites selected in the validation set
## Part 3: MetRS (top 15 metabolites) evaluated in the test set against Models 1-3 

## Set directories---------------------------------------------------------------
data_dir <- "path/to/data_derived/"   # folder containing input data
out_dir  <- "path/to/output/"         # folder for results

## Load packages-----------------------------------------------------------------
library(data.table)
library(caret)
library(glmnet)
library(survival)
library(boot)
library(dplyr)
library(tidyr)
library(parallel)

## Prepare data------------------------------------------------------------------
# nmr_dem.csv: analysis dataset (n = 269,284) with NMR metabolites, dementia outcomes and covariates
# Covariates were imputed using KNN (k = 10). Participants with prevalent or early-onset (<60 years) dementia,
# >20% missing metabolite values, missing APOE genotype, or diagnoses of both AD and VaD were excluded
nmr_dem <- fread(file.path(data_dir, "nmr_dem.csv"))
nmr_dem <- as.data.frame(nmr_dem)

# Participants with neuroimaging data (Desikan-Killiany and ASEG)
img <- fread(file.path(data_dir, "img_dk_aseg.csv"))
img <- as.data.frame(img)

# Exclude participants with neuroimaging data
cox_lasso_data <- subset(nmr_dem, !(eid %in% img$eid))

# Get metabnames
metab_names <- nmr_dem%>%
  select(bOHbutyrate:length(.)) %>%
  colnames()

## Split the data to train, validation and test----------------------------------
set.seed(777)

# index 1
index1 <- createDataPartition(cox_lasso_data$int_dem, p = 0.5, list = FALSE)

train_data <- cox_lasso_data[index1, ]
opti_data <- cox_lasso_data[-index1, ]

# index 2
index2 <- createDataPartition(opti_data$int_dem, p = 0.5, list = FALSE)

valid_data <- opti_data[index2, ]
test_data <- opti_data[-index2, ]

# List of datasets to iterate over
datasets <- list(train_data = train_data, valid_data = valid_data, test_data = test_data)

# Function to apply rank-based normal transformation and scaling
transform_metabolites <- function(df, metab_names) {
  df[, metab_names] <- apply(df[, metab_names], 2, function(x) {
    x_transformed <- qnorm((rank(x, na.last = 'keep') - 0.5) / sum(!is.na(x)))
    return(scale(x_transformed))
  })
  return(df)
}

# Apply the transformation to each dataset separately
datasets <- lapply(datasets, transform_metabolites, metab_names = metab_names)

train_data <- datasets$train_data
valid_data <- datasets$valid_data
test_data <- datasets$test_data

# Train data
y_train <- as.matrix(Surv(train_data$time_dem, train_data$int_dem))

metab_train <- train_data %>%
  select(bOHbutyrate:length(.)) %>%
  as.matrix()

## Part 1: LASSO-Cox with cross-validation in the training set--------------------
# Function to run glmnet model
run_iteration <- function(X) {
  set.seed(X)
  ## Cross validation for lambda_model metabolites
  cvfit_metab <- cv.glmnet(metab_train, y_train, family = 'cox', alpha = 1, type.measure = 'C')
  
  # Fit the final model using the optimal lambda
  lasso.model.all <- glmnet(metab_train, y_train, alpha =1, lambda = cvfit_metab$lambda.1se, family = "cox")
  
  # Extract coefficients
  cox_lasso_coef <- as.matrix(coef(lasso.model.all))
  cox_lasso_coef <- as.data.frame(cox_lasso_coef)
  cox_lasso_coef$metab <- rownames(cox_lasso_coef)
  # Filter non-zero coefficients
  cox_lasso_coef_m1 <- cox_lasso_coef[cox_lasso_coef$s0 != 0, ]
  cox_lasso_coef_m1$seed <- X
  
  return(list(
    cvfit_metab = cvfit_metab,
    cox_lasso_coef_m1 = cox_lasso_coef_m1
  ))
  
}

# Repeat over 100 random seeds (parallel processing)
set.seed(777)
seed_iterations <- sample(1:5000, 100)
res <- mclapply(seed_iterations, run_iteration, mc.cores = 20)
names(res) <- seed_iterations

# Save results
saveRDS(res, file.path(out_dir, "lasso_dem_full_3splits.rds"))




# Function to calculate median, SE, 95% CI and p value of bootstrapped values
get_bootstrap_summary <- function(boot_result, metric_names) {
  # Calculate the median of each bootstrapped metric, ignoring NA values
  median_vals <- apply(boot_result$t, 2, function(x) median(x, na.rm = TRUE))
  # Calculate the standard error (SE) for each metric, ignoring NA values
  se_vals <- apply(boot_result$t, 2, function(x) sd(x, na.rm = TRUE))
  # Calculate 95% confidence intervals for each metric, ignoring NA values
  ci_vals <- apply(boot_result$t, 2, function(x) quantile(x, c(0.025, 0.975), na.rm = TRUE))
  # Calculate p-values based on Z-scores
  p_estimate <- 2 * pnorm(-abs(median_vals / se_vals))
  # Combine results into a data frame
  summary_df <- data.frame(
    Metric = metric_names,
    Median = median_vals,
    SE = se_vals,
    CI_Lower = ci_vals[1, ],  # 2.5% percentile
    CI_Upper = ci_vals[2, ],  # 97.5% percentile
    p_value = p_estimate
  )
  return(summary_df)
}


## Part 2: Select the number of metabolites in the validation set----------------
# Summarise selected metabolites across iterations
# Mean regularised coefficient across 100 iterations (0 when not selected)
cox_lasso_coef_m1<- lapply(res, function(x) x$cox_lasso_coef_m1)
coef_df <- do.call(rbind, cox_lasso_coef_m1)

df_summary <- coef_df %>%
  dplyr::group_by(metab) %>%
  dplyr::summarise(mean_s0 = sum(s0)/100) %>%
  dplyr::arrange(desc(abs(mean_s0)))

valid_mrs <- valid_data
## Calculate metabolomic risk scores (MetRS) adding top-ranked metabolites one at a time
for (i in 1:length(df_summary$metab)) {
  top_metab <- df_summary[1:i,]
  
  #reshape individual metabolomic data
  valid_df <- valid_data[,c("eid",top_metab$metab)] #only keep eid and metabolites data
  valid_df_reshaped <- valid_df%>% gather(key = "metab", value = "value",-eid) #reshape individual metabolomic data to 3 columns("eid",metab","coef”)
  
  ##metabolomic risk score
  nmr_df_coef <- inner_join(valid_df_reshaped, top_metab, by = "metab") #merge reshaped individual metabolomic data with glmnet coefficients
  
  ##In each individual,calculate the sum of original metabolite value multiply coefficients of corresponding metabolites
  mrs_scores <- nmr_df_coef %>%
    dplyr::group_by(eid) %>%
    dplyr::summarise(mrs = sum(value * mean_s0, na.rm = TRUE))
  
  colnames(mrs_scores) <- c("eid",paste0("mrs_",i))
  
  valid_mrs <- valid_mrs %>%
    dplyr::inner_join(., mrs_scores,by="eid")
  
}
write.csv(valid_mrs, file = file.path(out_dir, "lasso_valid_mrs_dem_3splits.csv"), row.names=FALSE)

## C-index of each MetRS
for (i in 1:length(df_summary$metab)){
  cox_model_valid <- coxph(Surv(valid_mrs$time_dem, valid_mrs$int_dem) ~ valid_mrs[, paste0("mrs_", i)])
  df_summary[i,"Cindex_valid"] <- concordance(cox_model_valid)$concordance
}

write.csv(df_summary, file.path(out_dir, "lasso_coef_cindex_dem_3splits.csv"))

## Bootstrap comparison of C-indices between sets of metabolites
# Define a bootstrap function to calculate the difference in C-index
bootstrap_c_index_nmetab <- function(data, indices) {
  # Subset the data based on bootstrap indices
  d <- data[indices, ]
  # Initialize lists to store models, C-index values, and differences
  cox_boot <- vector("list", length(df_summary$metab))
  C <- numeric(length(df_summary$metab))
  Cdiff <- numeric(length(df_summary$metab) - 1)
  # Loop through the metabolites
  for (i in 1:length(df_summary$metab)) {
    # Dynamically construct the formula for each model
    formula <- as.formula(paste0("Surv(time_dem, int_dem) ~ mrs_", i))
    cox_boot[[i]] <- coxph(formula, data = d)
    
    # Calculate the C-index for each model
    C[i] <- concordance(cox_boot[[i]])$concordance
  }
  
  
  # Calculate specific differences requested
  Cdiff1 <- C[10] - C[5]
  Cdiff2 <- C[15] - C[10]
  Cdiff3 <- C[20] - C[15]
  Cdiff4 <- C[25] - C[20]
  Cdiff5 <- C[30] - C[25]
  Cdiff6 <- C[35] - C[30]
  
  # Return the C-indexes, their differences, and specified additional differences
  return(c(C, Cdiff1 = Cdiff1, Cdiff2 = Cdiff2, Cdiff3 = Cdiff3,
           Cdiff4 = Cdiff4, Cdiff5 = Cdiff5, Cdiff6 = Cdiff6))
}

# Perform the bootstrapping on each dataset and get summary tables
set.seed(777)
boot_results_valid <- boot(data = valid_mrs, statistic = bootstrap_c_index_nmetab, R = 1000, parallel = "multicore")

# Save the bootstrap results
saveRDS(boot_results_valid, file.path(out_dir, "lasso_boot_valid_Cindex1000_dem.rds"))

# Get summary tables for each bootstrap result
valid_summary <- get_bootstrap_summary(boot_results_valid,
                                       metric_names = c(paste0("C", 1:length(df_summary$metab)),
                                                        "Cdiff1", "Cdiff2", "Cdiff3", "Cdiff4", "Cdiff5", "Cdiff6"))

write.csv(valid_summary, file.path(out_dir, "lasso_res_Cindex1000_datasets_dem.csv"), row.names = F)

df_summary$Cindex_median_valid <- valid_summary$Median[1:length(df_summary$metab)]

write.csv(df_summary, file.path(out_dir, "lasso_coef_cindex1000_withmedian_dem_3splits.csv"))

## Part 3: Predictive performance of the MetRS in the test set--------------------
# Top 15 metabolites (selected in Part 2)
df_summary_top <- df_summary[1:15,]

test_mrs <- test_data
## Calculate metabolomic risk score (MetRS)

#reshape individual metabolomic data
test_df <- test_data[,c("eid",df_summary_top$metab)] #only keep eid and metabolites data
test_df_reshaped <- test_df%>% gather(key = "metab", value = "value",-eid) #reshape individual metabolomic data to 3 columns("eid",metab","coef”)

##metabolomic risk score
nmr_df_coef <- inner_join(test_df_reshaped, df_summary_top, by = "metab") #merge reshaped individual metabolomic data with glmnet coefficients

##In each individual,calculate the sum of original metabolite value multiply coefficients of corresponding metabolites
mrs_scores <- nmr_df_coef %>%
  dplyr::group_by(eid) %>%
  dplyr::summarise(top_mrs = sum(value * mean_s0, na.rm = TRUE))

test_mrs <- test_mrs %>%
  dplyr::inner_join(., mrs_scores,by="eid")

write.csv(test_mrs, file = file.path(out_dir, "test_topmrs_dem_3splits.csv"), row.names=FALSE)

test_mrs$townsend <- as.numeric(test_mrs$townsend)
test_mrs$centre_0 <- as.factor(test_mrs$centre_0)
# Define the mapping for codes to countries
test_mrs$centre_0 <- ifelse(test_mrs$centre_0 %in% c(11012, 11021, 11011, 11008,
                                                     11024, 11020, 11018,
                                                     11010, 11016, 11001, 11017, 11009,
                                                     11013, 11002, 11007, 11014, 11006), "England",
                            ifelse(test_mrs$centre_0 %in% c(11003, 11022, 11023), "Wales",
                                   ifelse(test_mrs$centre_0 %in% c(11004, 11005), "Scotland", NA)))

## Bootstrap comparison of C-indices with and without MetRS----------------------
# cox1: MetRS only; cox2/cox3: Model 1 without/with MetRS; cox4/cox5: Model 2 without/with MetRS; cox6/cox7: Model 3 without/with MetRS
# Model 1: age, age2, sex, age*sex, age2*sex, centre, ethnicity, Townsend index; Model 2: Model 1 + APOE e4;
# Model 3: Model 2 + education, smoking, alcohol, obesity, diabetes, hypertension
# Define a bootstrap function to calculate the difference in C-index
bootstrap_c_index_models <- function(data, indices) {
  # Subset the data based on bootstrap indices
  d <- data[indices, ]
  
  # Fit Cox models on the subsetted data
  cox1 <- coxph(Surv(time_dem, int_dem) ~ top_mrs, data=d)
  cox2 <- coxph(Surv(time_dem, int_dem) ~ centre_0 + ethn + townsend + age + age2 + sex + age*sex + age2*sex, data=d)
  cox3 <- coxph(Surv(time_dem, int_dem) ~ centre_0 + ethn + townsend + age + age2 + sex + age*sex + age2*sex + top_mrs, data=d)
  cox4 <- coxph(Surv(time_dem, int_dem) ~ centre_0 + ethn + townsend + age + age2 + sex + age*sex + age2*sex + apoe, data=d)
  cox5 <- coxph(Surv(time_dem, int_dem) ~ centre_0 + ethn + townsend + age + age2 + sex + age*sex + age2*sex + apoe + top_mrs, data=d)
  cox6 <- coxph(Surv(time_dem, int_dem) ~ centre_0 + ethn + townsend + age + age2 + sex + age*sex + age2*sex + apoe + edu + smoking + alcohol + obesity + diabetes + hypertension, data=d)
  cox7 <- coxph(Surv(time_dem, int_dem) ~ centre_0 + ethn + townsend + age + age2 + sex + age*sex + age2*sex + apoe + edu + smoking + alcohol + obesity + diabetes + hypertension + top_mrs, data=d)
  
  # Calculate concordance for each model (with error handling)
  C1 <- tryCatch(concordance(cox1)$concordance, error = function(e) NA)
  C2 <- tryCatch(concordance(cox2)$concordance, error = function(e) NA)
  C3 <- tryCatch(concordance(cox3)$concordance, error = function(e) NA)
  C4 <- tryCatch(concordance(cox4)$concordance, error = function(e) NA)
  C5 <- tryCatch(concordance(cox5)$concordance, error = function(e) NA)
  C6 <- tryCatch(concordance(cox6)$concordance, error = function(e) NA)
  C7 <- tryCatch(concordance(cox7)$concordance, error = function(e) NA)
  
  # Compute differences
  Cdiff1 <- C3 - C2
  Cdiff2 <- C5 - C4
  Cdiff3 <- C7 - C6
  
  # Return the results
  return(c(C1=C1, C2=C2, C3=C3, C4=C4, C5=C5, C6=C6, C7=C7, Cdiff1=Cdiff1, Cdiff2=Cdiff2, Cdiff3=Cdiff3))
}

# Perform the bootstrapping
set.seed(777)
boot_results_test <- boot(data = test_mrs, statistic = bootstrap_c_index_models, R = 2000, parallel = "multicore")

# Save the bootstrap results
saveRDS(boot_results_test, file.path(out_dir, "boot_Cindex_test_dem.rds"))

# Get summary tables for each bootstrap result
test_summary <- get_bootstrap_summary(boot_results_test,
                                      metric_names = c("C1", "C2", "C3", "C4", "C5", "C6", "C7",
                                                       "Cdiff1", "Cdiff2", "Cdiff3"))

write.csv(test_summary, file.path(out_dir, "lasso_res_Cindex_test_dem.csv"), row.names = F)