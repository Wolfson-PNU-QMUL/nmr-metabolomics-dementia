## Example: Bidirectional two-sample MR: key metabolites <-> dementia (ACD, AD, VaD) and WMH phenotypes
## Part 1: instrument selection, harmonisation and IVW-MR 
## Part 2: FDR correction, instrument strength, and sensitivity analyses
##         (weighted median, weighted mode, MR-Egger bootstrap, Cochran's Q,
##          Egger intercept, leave-one-out, MR-PRESSO)

## Load packages
library(TwoSampleMR)
library(ieugwasr)
library(data.table)
library(dplyr)
library(LDlinkR)
library(tidyverse)
library(MRPRESSO)

## Paths (edit to your local setup)
metab_gwas_dir <- "GWAS_Sumstat/NMR_Karjalainen2024/"   # metabolite GWAS: <accession>_withrsid.tsv
dem_gwas_dir   <- "GWAS_Sumstat/dementia/"              # outcome GWAS: <outcome>_munged.tsv
ld_ref         <- "LD_reference/EUR/EUR"                # 1000 Genomes Phase 3 EUR, PLINK format
plink_bin      <- "LD_reference/EUR/plink"
out_dir        <- "results/"
ldlink_token   <- Sys.getenv("LDLINK_TOKEN")            # personal LDlink token, set in .Renviron

dir.create(paste0(out_dir, "metab/harmonised"), recursive = TRUE, showWarnings = FALSE)
dir.create(paste0(out_dir, "dem/harmonised"), recursive = TRUE, showWarnings = FALSE)

outcomes <- c("ACD", "AD", "VaD", "WMHs", "WMHs_atrophy")
dementia_traits <- c("ACD", "AD", "VaD")   # APOE region excluded from these GWAS

#### Functions--------------------------------------------------------------------

# Replace missing outcome SNPs with LD proxies (r2 >= 0.8; highest r2, then closest distance)
munge_proxies <- function(LDLink_file, outcome, outcome_clump){
  LDLink_file_path <- LDLink_file
  proxy_snps <- read_tsv(LDLink_file_path, skip = 1, col_names = F) %>%
    dplyr::rename(id = X1, func = X2, proxy_snp = X3, coord = X4, alleles = X5, maf = X6,
                  distance = X7, dprime = X8, rsq = X9, correlated_alleles = X10, FORGEdb = X11, RegulomeDB = X12) %>%
    separate(coord, c('chr', 'pos'), sep = ":") %>%
    mutate(snp = ifelse(id == 1, proxy_snp, NA),
           chr = str_replace(chr, 'chr', ""),
           chr = as.numeric(chr),
           pos = as.numeric(pos)) %>%
    fill(snp, .direction = 'down') %>%
    relocate(snp, .before = proxy_snp) %>%
    dplyr::select(-id, -func, -FORGEdb, -RegulomeDB) %>%
    filter(rsq >= 0.8)
  
  # Munge proxy snp and outcome data
  proxy_outcome <- left_join(
    proxy_snps, outcome, by = c("proxy_snp" = "SNP")
  ) %>%
    separate(correlated_alleles, c("target_a1.outcome", "proxy_a1.outcome",
                                   "target_a2.outcome", "proxy_a2.outcome"), sep = ",|=") %>%
    filter(!is.na(chr.outcome)) %>%
    arrange(snp, -rsq, abs(distance)) %>%
    group_by(snp) %>%
    slice(1) %>%
    ungroup() %>%
    mutate(
      proxy.outcome = TRUE,
      target_snp.outcome = snp,
      proxy_snp.outcome = proxy_snp,
    ) %>%
    mutate(
      new_effect_allele.outcome = case_when(
        proxy_a1.outcome == effect_allele.outcome & proxy_a2.outcome == other_allele.outcome ~ target_a1.outcome,
        proxy_a2.outcome == effect_allele.outcome & proxy_a1.outcome == other_allele.outcome ~ target_a2.outcome,
        TRUE ~ NA_character_
      ),
      new_other_allele.outcome = case_when(
        proxy_a1.outcome == effect_allele.outcome & proxy_a2.outcome == other_allele.outcome ~ target_a2.outcome,
        proxy_a2.outcome == effect_allele.outcome & proxy_a1.outcome == other_allele.outcome ~ target_a1.outcome,
        TRUE ~ NA_character_
      ),
      effect_allele.outcome = new_effect_allele.outcome,
      other_allele.outcome = new_other_allele.outcome
    ) %>%
    dplyr::select(-proxy_snp, -chr, -pos, -alleles, -maf, -distance, -rsq, -dprime,
                  -new_effect_allele.outcome, -new_other_allele.outcome) %>%
    relocate(target_a1.outcome, proxy_a1.outcome, target_a2.outcome, proxy_a2.outcome, .after = proxy_snp.outcome) %>%
    dplyr::rename(SNP = snp) %>%
    relocate(SNP, .after = samplesize.outcome)
  
  # Merge outcome and proxy outcomes
  outcome_dat <- bind_rows(
    outcome_clump, proxy_outcome
  ) %>%
    arrange(chr.outcome, pos.outcome)
  
  outcome_dat
}

# Wrapper function for read_exposure_data
read_exposure_data_with_formula <- function(filename) {
  read_exposure_data(
    filename = filename,
    sep = "\t",
    snp_col = "SNP",
    beta_col = "BETA",
    se_col = "SE",
    effect_allele_col = "A2",
    other_allele_col = "A1",
    eaf_col = "FRQ",
    pval_col = "P",
    samplesize_col = "N",
    min_pval = 1e-200,
    chr_col = "CHR",
    pos_col = "BP",
    log_pval = FALSE
  )
}

# Wrapper function for read_outcome_data (snps = NULL reads all SNPs)
read_outcome_data_with_formula <- function(filename, snps = NULL) {
  read_outcome_data(
    filename = filename,
    snps = snps,
    sep = "\t",
    snp_col = "SNP",
    beta_col = "BETA",
    se_col = "SE",
    effect_allele_col = "A2",
    other_allele_col = "A1",
    eaf_col = "FRQ",
    pval_col = "P",
    samplesize_col = "N",
    min_pval = 1e-200,
    chr_col = "CHR",
    pos_col = "BP",
    log_pval = FALSE
  )
}

# Restrict to MAF > 0.01, clump, apply p < 5e-8 (or p < 5e-6 if fewer than 5 IVs),
# remove weak instruments (F < 10); returns NULL if fewer than 5 IVs remain
select_ivs <- function(exp_dat, exposure_name) {
  exp_dat$exposure <- exposure_name
  exp_dat$rsid <- exp_dat$SNP
  exp_dat$pval <- exp_dat$pval.exposure
  exp_dat$id <- exp_dat$exposure
  
  exp_dat <- exp_dat[exp_dat$pval.exposure < 5e-6 &
                       exp_dat$eaf.exposure > 0.01 & exp_dat$eaf.exposure < 0.99, ]
  exp_dat <- ld_clump(
    dat = exp_dat,
    clump_kb = 10000,
    clump_r2 = 0.001,
    bfile = ld_ref,
    plink_bin = plink_bin
  )
  
  if (nrow(exp_dat[exp_dat$pval.exposure < 5e-8, ]) >= 5) {
    exp_dat <- exp_dat[exp_dat$pval.exposure < 5e-8, ]
  }
  exp_dat <- exp_dat[(exp_dat$beta.exposure^2) / (exp_dat$se.exposure^2) >= 10, ]
  if (nrow(exp_dat) < 5) return(NULL)
  exp_dat
}

# Read outcome data for the IVs; use LD proxies for IVs missing in the outcome
get_outcome <- function(outcome_file, outcome_name, exp_dat) {
  out_dat <- read_outcome_data_with_formula(outcome_file, snps = exp_dat$SNP)
  out_dat$outcome <- outcome_name
  
  exp_snps_wo <- anti_join(exp_dat, out_dat, by = "SNP") # Exposure SNPs not present in outcome
  
  if (nrow(exp_snps_wo) > 0) {
    out_all <- read_outcome_data_with_formula(outcome_file)
    out_all$outcome <- outcome_name
    outcome_clump <- semi_join(out_all, exp_dat, by = "SNP") # exposure SNPs present in outcome
    
    LDproxy_batch(exp_snps_wo$SNP,
                  pop = "EUR",
                  r2d = "r2",
                  token = ldlink_token,
                  append = T,
                  genome_build = "grch37")
    out_dat <- munge_proxies("combined_query_snp_list_grch37.txt", out_all, outcome_clump)
    file.remove("combined_query_snp_list_grch37.txt")
  }
  out_dat
}

# Harmonise; exclude palindromic SNPs with MAF > 0.42 (action = 2) and,
# for analyses involving dementia GWAS (ACD, AD, VaD), the APOE region (chr19:45,020,859-45,844,508, GRCh37)
harmonise_mr <- function(exp_dat, out_dat, dementia_gwas) {
  harmonise_data(exp_dat, out_dat, action = 2) %>%
    as_tibble() %>%
    mutate(
      apoe_region = chr.outcome == 19 & between(pos.outcome, 45020859, 45844508),
      mr_keep = ifelse(mr_keep == FALSE | (dementia_gwas & apoe_region), FALSE, TRUE)
    )
}

#### Part 1. Primary IVW-MR-------------------------------------------------------

keymetab_list <- fread("keymetab_list.csv", header = T)   # columns: metab, accession

### 1a. metabolites >> dementia / WMHs
MRres_metab_dem <- data.frame()

for (f in 1:nrow(keymetab_list)) {
  exp_dat <- read_exposure_data_with_formula(
    paste0(metab_gwas_dir, keymetab_list$accession[f], "_withrsid.tsv")
  )
  exp_dat <- select_ivs(exp_dat, keymetab_list$metab[f])
  if (is.null(exp_dat)) next   # fewer than 5 IVs at p < 5e-6
  
  for (i in outcomes) {
    out_dat <- get_outcome(paste0(dem_gwas_dir, i, "_munged.tsv"), i, exp_dat)
    dat <- harmonise_mr(exp_dat, out_dat, dementia_gwas = i %in% dementia_traits)
    write.csv(dat, paste0(out_dir, "metab/harmonised/harmonised_", keymetab_list$metab[f], "_", i, ".csv"), row.names = F)
    
    res <- mr(dat, method_list = c("mr_ivw_mre"))
    
    MRres_metab_dem <- rbind(
      MRres_metab_dem,
      data.frame(exp = keymetab_list$metab[f], out = i, path = "metab>>dementia",
                 coef = res$b, se = res$se, p = res$pval)
    )
  }
}

### 1b. dementia / WMHs >> metabolites
MRres_dem_metab <- data.frame()

for (f in outcomes) {
  exp_dat <- read_exposure_data_with_formula(paste0(dem_gwas_dir, f, "_munged.tsv"))
  exp_dat <- select_ivs(exp_dat, f)
  if (is.null(exp_dat)) next
  
  for (j in 1:nrow(keymetab_list)) {
    out_dat <- get_outcome(paste0(metab_gwas_dir, keymetab_list$accession[j], "_withrsid.tsv"),
                           keymetab_list$metab[j], exp_dat)
    dat <- harmonise_mr(exp_dat, out_dat, dementia_gwas = f %in% dementia_traits)
    write.csv(dat, paste0(out_dir, "dem/harmonised/harmonised_", f, "_", keymetab_list$metab[j], ".csv"), row.names = F)
    
    res <- mr(dat, method_list = c("mr_ivw_mre"))
    
    MRres_dem_metab <- rbind(
      MRres_dem_metab,
      data.frame(exp = f, out = keymetab_list$metab[j], path = "dem>>metab",
                 coef = res$b, se = res$se, p = res$pval)
    )
  }
}

#### Part 2. FDR, instrument strength and sensitivity analyses---------------------

res_IVW <- rbind(MRres_metab_dem, MRres_dem_metab)
res_IVW$FDR <- p.adjust(res_IVW$p, method = "BH")

# Read harmonised data and keep the SNPs used in the MR
read_harmonised <- function(path, exp, out) {
  folder <- ifelse(path == "metab>>dementia", "metab", "dem")
  dat <- read.csv(paste0(out_dir, folder, "/harmonised/harmonised_", exp, "_", out, ".csv"))
  dat[dat$mr_keep == TRUE, ]
}

### 2a. Number of IVs, F-statistics and I2
for (f in 1:nrow(res_IVW)) {
  dat <- read_harmonised(res_IVW$path[f], res_IVW$exp[f], res_IVW$out[f])
  F_stat <- (dat$beta.exposure^2) / (dat$se.exposure^2)
  
  res_IVW[f, "method"] <- "Inverse variance weighted (multiplicative random effects)"
  res_IVW[f, "p_threshold"] <- ifelse(all(dat$pval.exposure < 5e-8), 5e-8, 5e-6)
  res_IVW[f, "nsnp"] <- nrow(dat)
  res_IVW[f, "Min_F"] <- min(F_stat)
  res_IVW[f, "Max_F"] <- max(F_stat)
  res_IVW[f, "Mean_F"] <- mean(F_stat)
  res_IVW[f, "I_squared"] <- Isq(dat$beta.exposure, dat$se.exposure)
}

write.csv(res_IVW, paste0(out_dir, "res_IVW_keymetab_bidir.csv"), row.names = F)

### Sensitivity analyses on results with pFDR < 0.1
res_IVW_sig <- res_IVW[res_IVW$FDR < 0.1, ]

res_otherMR <- data.frame()
HPres <- data.frame()
loo_results <- data.frame()
Presso_res <- data.frame()

for (f in 1:nrow(res_IVW_sig)) {
  dat <- read_harmonised(res_IVW_sig$path[f], res_IVW_sig$exp[f], res_IVW_sig$out[f])
  
  # 2b. Weighted median, weighted mode, MR-Egger (bootstrap SE)
  set.seed(777)
  mr_res <- mr(dat, method_list = c("mr_weighted_median", "mr_weighted_mode", "mr_egger_regression_bootstrap"))
  wm  <- mr_res[mr_res$method == "Weighted median", ]
  wmo <- mr_res[mr_res$method == "Weighted mode", ]
  eg  <- mr_res[mr_res$method == "MR Egger (bootstrap)", ]
  
  res_otherMR[f, "path"] <- res_IVW_sig$path[f]
  res_otherMR[f, "exp"] <- res_IVW_sig$exp[f]
  res_otherMR[f, "out"] <- res_IVW_sig$out[f]
  res_otherMR[f, "beta_Wmedian"] <- wm$b
  res_otherMR[f, "se_Wmedian"] <- wm$se
  res_otherMR[f, "pval_Wmedian"] <- wm$pval
  res_otherMR[f, "beta_Wmode"] <- wmo$b
  res_otherMR[f, "se_Wmode"] <- wmo$se
  res_otherMR[f, "pval_Wmode"] <- wmo$pval
  res_otherMR[f, "beta_Egger"] <- eg$b
  res_otherMR[f, "se_Egger"] <- eg$se
  res_otherMR[f, "pval_Egger"] <- eg$pval
  
  # 2c. Cochran's Q (IVW) and MR-Egger intercept
  heterogeneity_results <- mr_heterogeneity(dat)
  het_ivw <- heterogeneity_results[heterogeneity_results$method == "Inverse variance weighted", ]
  pleiotropy <- mr_pleiotropy_test(dat)
  
  HPres[f, "path"] <- res_IVW_sig$path[f]
  HPres[f, "exp"] <- res_IVW_sig$exp[f]
  HPres[f, "out"] <- res_IVW_sig$out[f]
  HPres[f, "Q"] <- het_ivw$Q
  HPres[f, "Q_df"] <- het_ivw$Q_df
  HPres[f, "Q_pval"] <- het_ivw$Q_pval
  HPres[f, "egger_intercept"] <- pleiotropy$egger_intercept
  HPres[f, "egger_intercept_se"] <- pleiotropy$se
  HPres[f, "egger_intercept_pval"] <- pleiotropy$pval
  
  # 2d. Leave-one-out (IVW, multiplicative random effects)
  res_loo <- mr_leaveoneout(dat, method = mr_ivw_mre)
  loo_results <- rbind(loo_results,
                       data.frame(path = res_IVW_sig$path[f],
                                  exp = res_IVW_sig$exp[f],
                                  out = res_IVW_sig$out[f],
                                  excluded_SNP = res_loo$SNP,
                                  beta = res_loo$b,
                                  se = res_loo$se,
                                  pval = res_loo$p))
  
  # 2e. MR-PRESSO (global, outlier and distortion tests)
  set.seed(777)
  mr_presso_results <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure",
                                 SdOutcome = "se.outcome", SdExposure = "se.exposure",
                                 OUTLIERtest = TRUE, DISTORTIONtest = TRUE,
                                 data = dat, NbDistribution = 1000, SignifThreshold = 0.05)
  presso <- mr_presso_results$`MR-PRESSO results`
  
  Presso_res[f, "path"] <- res_IVW_sig$path[f]
  Presso_res[f, "exp"] <- res_IVW_sig$exp[f]
  Presso_res[f, "out"] <- res_IVW_sig$out[f]
  Presso_res[f, "beta_corrected"] <- mr_presso_results$`Main MR results`$`Causal Estimate`[2]
  Presso_res[f, "p_corrected"] <- mr_presso_results$`Main MR results`$`P-value`[2]
  Presso_res[f, "RSSobs"] <- if (is.null(presso$`Global Test`$RSSobs)) NA else presso$`Global Test`$RSSobs
  Presso_res[f, "p_global"] <- if (is.null(presso$`Global Test`$Pvalue)) NA else presso$`Global Test`$Pvalue
  Presso_res[f, "beta_distortion"] <- if (is.null(presso$`Distortion Test`$`Distortion Coefficient`)) NA else presso$`Distortion Test`$`Distortion Coefficient`
  Presso_res[f, "p_distortion"] <- if (is.null(presso$`Distortion Test`$Pvalue)) NA else presso$`Distortion Test`$Pvalue
}

write.csv(res_otherMR, paste0(out_dir, "res_otherMR_bidir.csv"), row.names = F)
write.csv(HPres, paste0(out_dir, "HPres_bidir.csv"), row.names = F)
write.csv(loo_results, paste0(out_dir, "loo_res_bidir.csv"), row.names = F)
write.csv(Presso_res, paste0(out_dir, "Presso_res_bidir.csv"), row.names = F)