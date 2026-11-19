## Mediation MR: risk factor (RF) >> metabolite >> dementia
## Pairs are taken forward when both path a (RF >> metabolite) and path b (metabolite >> dementia)
## were significant in univariable MR and passed sensitivity analyses.
## Path a and total effect (c): univariable IVW-MR; path b and direct effect: IVW-MVMR (RF + metabolite,
## MVMR package), with conditional F-statistics and Q-statistic for instrument strength and pleiotropy.
## Indirect effect = a x b; proportion mediated = a x b / c; 95% CI by the delta method.

## Load packages
library(data.table)
library(tidyverse)
library(TwoSampleMR)
library(LDlinkR)
library(ieugwasr)
library(MVMR)

## Paths (edit to your local setup)
rf_gwas_dir    <- "GWAS_Sumstat/RF/"                    # <RF>_munged.tsv
metab_gwas_dir <- "GWAS_Sumstat/NMR_Karjalainen2024/"   # <accession>_withrsid.tsv
dem_gwas_dir   <- "GWAS_Sumstat/dementia/"              # <outcome>_munged.tsv
ld_ref         <- "LD_reference/EUR/EUR"
plink_bin      <- "LD_reference/EUR/plink"
out_dir        <- "results/mediation/"
ldlink_token   <- Sys.getenv("LDLINK_TOKEN")

dir.create(paste0(out_dir, "harmonised"), recursive = TRUE, showWarnings = FALSE)

#### Functions--------------------------------------------------------------------

# Replace missing SNPs with LD proxies (r2 >= 0.8; highest r2, then closest distance)
# type = "outcome" or "exposure"
munge_proxies <- function(LDLink_file, dat, dat_clump, type = "outcome"){
  proxy_snps <- read_tsv(LDLink_file, skip = 1, col_names = F) %>%
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
  
  # Rename columns to a common suffix, munge, then rename back
  sfx <- paste0(".", type)
  names(dat) <- sub(sfx, ".x_", names(dat), fixed = TRUE)
  
  proxy_dat <- left_join(proxy_snps, dat, by = c("proxy_snp" = "SNP")) %>%
    separate(correlated_alleles, c("target_a1", "proxy_a1", "target_a2", "proxy_a2"), sep = ",|=") %>%
    filter(!is.na(chr.x_)) %>%
    arrange(snp, -rsq, abs(distance)) %>%
    group_by(snp) %>%
    slice(1) %>%
    ungroup() %>%
    mutate(
      proxy.x_ = TRUE,
      target_snp.x_ = snp,
      proxy_snp.x_ = proxy_snp,
      new_effect_allele = case_when(
        proxy_a1 == effect_allele.x_ & proxy_a2 == other_allele.x_ ~ target_a1,
        proxy_a2 == effect_allele.x_ & proxy_a1 == other_allele.x_ ~ target_a2,
        TRUE ~ NA_character_
      ),
      new_other_allele = case_when(
        proxy_a1 == effect_allele.x_ & proxy_a2 == other_allele.x_ ~ target_a2,
        proxy_a2 == effect_allele.x_ & proxy_a1 == other_allele.x_ ~ target_a1,
        TRUE ~ NA_character_
      ),
      effect_allele.x_ = new_effect_allele,
      other_allele.x_ = new_other_allele
    ) %>%
    dplyr::select(-proxy_snp, -chr, -pos, -alleles, -maf, -distance, -rsq, -dprime,
                  -new_effect_allele, -new_other_allele,
                  -target_a1, -proxy_a1, -target_a2, -proxy_a2) %>%
    dplyr::rename(SNP = snp)
  
  names(proxy_dat) <- sub(".x_", sfx, names(proxy_dat), fixed = TRUE)
  
  bind_rows(dat_clump, proxy_dat)
}

# Keep the SNPs in snps; use LD proxies for SNPs missing in dat
add_proxies <- function(dat, snps, type = "outcome") {
  dat_overlap <- dat[dat$SNP %in% snps, ]
  snps_wo <- setdiff(snps, dat$SNP)
  if (length(snps_wo) == 0) return(dat_overlap)
  
  if (file.exists("combined_query_snp_list_grch37.txt")) file.remove("combined_query_snp_list_grch37.txt")  # avoid appending to an old file
  LDproxy_batch(snps_wo,
                pop = "EUR",
                r2d = "r2",
                token = ldlink_token,
                append = T,
                genome_build = "grch37")
  dat_final <- munge_proxies("combined_query_snp_list_grch37.txt", dat, dat_overlap, type)
  file.remove("combined_query_snp_list_grch37.txt")
  dat_final
}

# Functions for reading GWAS summary statistics
read_exposure_data_with_formula <- function(filename) {
  read_exposure_data(
    filename = filename, sep = "\t",
    snp_col = "SNP", beta_col = "BETA", se_col = "SE",
    effect_allele_col = "A2", other_allele_col = "A1", eaf_col = "FRQ",
    pval_col = "P", samplesize_col = "N", min_pval = 1e-200,
    chr_col = "CHR", pos_col = "BP", log_pval = FALSE
  )
}

read_outcome_data_with_formula <- function(filename) {
  read_outcome_data(
    filename = filename, sep = "\t",
    snp_col = "SNP", beta_col = "BETA", se_col = "SE",
    effect_allele_col = "A2", other_allele_col = "A1", eaf_col = "FRQ",
    pval_col = "P", samplesize_col = "N", min_pval = 1e-200,
    chr_col = "CHR", pos_col = "BP", log_pval = FALSE
  )
}

# IV selection (as in the primary MR): MAF > 0.01, clumping (r2 < 0.001, 10,000 kb),
# p < 5e-8 (or p < 5e-6 if fewer than 5 IVs), F >= 10, APOE region excluded (outcome is dementia);
# returns NULL if fewer than 5 IVs remain
select_ivs <- function(exp_dat) {
  exp_dat$rsid <- exp_dat$SNP
  exp_dat$pval <- exp_dat$pval.exposure
  exp_dat$id <- exp_dat$exposure
  
  exp_dat <- exp_dat[which(exp_dat$pval.exposure < 5e-6 &
                             exp_dat$eaf.exposure > 0.01 & exp_dat$eaf.exposure < 0.99), ]
  exp_dat <- ld_clump(dat = exp_dat, clump_kb = 10000, clump_r2 = 0.001,
                      bfile = ld_ref, plink_bin = plink_bin)
  exp_dat <- exp_dat %>%
    filter(!(chr.exposure == 19 & between(pos.exposure, 45020859, 45844508)))
  
  if (nrow(exp_dat[exp_dat$pval.exposure < 5e-8, ]) >= 5) {
    exp_dat <- exp_dat[exp_dat$pval.exposure < 5e-8, ]
  }
  exp_dat <- exp_dat[(exp_dat$beta.exposure^2) / (exp_dat$se.exposure^2) >= 10, ]
  if (nrow(exp_dat) < 5) return(NULL)
  exp_dat
}

#### Select RF-metabolite-dementia pairs-------------------------------------------

mr_res_RF <- read.csv("mr_res_RF.csv")
mr_res_metab <- read.csv("mr_res_metab.csv")

patha_uniMR <- mr_res_RF %>% filter(sensitivity == "PASS" & path == "RF>>metab")
pathb_uniMR <- mr_res_metab %>% filter(sensitivity == "PASS" & path == "metab>>dementia")

# After merging: exp = RF, out = metabolite, out.y = dementia outcome
sig_pairs <- merge(patha_uniMR, pathb_uniMR, by.x = "out", by.y = "exp")

keymetab_list <- fread("keymetab_list.csv")   # columns: metab, accession
sig_pairs <- merge(sig_pairs, keymetab_list, by.x = "out", by.y = "metab")

#### Mediation MR-----------------------------------------------------------------

res_paths <- data.frame()
res_mediation <- data.frame()
res_mvmr_diag <- data.frame()

for (f in 1:nrow(sig_pairs)) {
  rf    <- sig_pairs$exp[f]
  metab <- sig_pairs$out[f]
  dem   <- sig_pairs$out.y[f]
  
  # Read data
  exposure1 <- read_exposure_data_with_formula(paste0(rf_gwas_dir, rf, "_munged.tsv"))
  exposure1$exposure <- rf
  
  exposure2 <- read_exposure_data_with_formula(paste0(metab_gwas_dir, sig_pairs$accession[f], "_withrsid.tsv"))
  exposure2$exposure <- metab
  
  outcome <- read_outcome_data_with_formula(paste0(dem_gwas_dir, dem, "_munged.tsv"))
  outcome$outcome <- dem
  
  exp2_out <- read_outcome_data_with_formula(paste0(metab_gwas_dir, sig_pairs$accession[f], "_withrsid.tsv"))
  exp2_out$outcome <- metab
  
  # IV selection for each exposure
  exposure1_dat <- select_ivs(exposure1)
  exposure2_dat <- select_ivs(exposure2)
  if (is.null(exposure1_dat) || is.null(exposure2_dat)) next
  
  ## 1. Path a (RF >> metabolite) and total effect c (RF >> dementia): univariable IVW-MR
  patha_mrdat <- harmonise_data(exposure1_dat, add_proxies(exp2_out, exposure1_dat$SNP, "outcome"), action = 2)
  res_a <- mr(patha_mrdat, method_list = c("mr_ivw_mre"))
  
  pathc_mrdat <- harmonise_data(exposure1_dat, add_proxies(outcome, exposure1_dat$SNP, "outcome"), action = 2)
  res_c <- mr(pathc_mrdat, method_list = c("mr_ivw_mre"))
  
  ## 2. MVMR (RF + metabolite >> dementia)
  # Union of IVs from both exposures
  mvmr_snps <- unique(c(exposure1_dat$SNP, exposure2_dat$SNP))
  
  # Extract the SNPs from both exposures, using proxies where missing
  exp1_final <- add_proxies(exposure1, mvmr_snps, "exposure")
  exp2_final <- add_proxies(exposure2, mvmr_snps, "exposure")
  comb_exp <- bind_rows(exp1_final, exp2_final)
  
  # Second round of clumping on the combined SNP list (minimum p across the two exposures)
  comb_exp_wide <- comb_exp %>%
    group_by(SNP) %>%
    summarise(pval = min(pval.exposure)) %>%
    mutate(rsid = SNP, id = "mvmr")
  clump_comb_exp <- ld_clump(dat = comb_exp_wide, clump_kb = 10000, clump_r2 = 0.001,
                             bfile = ld_ref, plink_bin = plink_bin)
  comb_exp_independent <- comb_exp %>% filter(SNP %in% clump_comb_exp$rsid)
  
  # Outcome data for the MVMR SNPs
  out_final <- add_proxies(outcome, unique(comb_exp_independent$SNP), "outcome")
  
  mvdat <- mv_harmonise_data(comb_exp_independent, out_final)
  saveRDS(mvdat, paste0(out_dir, "harmonised/harmonised_mvdat_", rf, "_", metab, "_", dem, ".rds"))
  
  # Format for the MVMR package; columns of exposure_beta are exposure ids
  r_input <- format_mvmr(BXGs = mvdat$exposure_beta,
                         BYG = mvdat$outcome_beta,
                         seBXGs = mvdat$exposure_se,
                         seBYG = mvdat$outcome_se,
                         RSID = rownames(mvdat$exposure_beta))
  exp_names <- mvdat$expname$exposure[match(colnames(mvdat$exposure_beta), mvdat$expname$id.exposure)]
  
  # IVW-MVMR
  mvmr_fit <- ivw_mvmr(r_input, gencov = 0)
  mvmr_res <- data.frame(exposure = exp_names,
                         b = mvmr_fit[, 1],
                         se = mvmr_fit[, 2],
                         pval = mvmr_fit[, 4])
  
  # Conditional F-statistics (instrument strength) and Q-statistic (heterogeneity/pleiotropy)
  cond_F <- strength_mvmr(r_input, gencov = 0)
  q_res <- pleiotropy_mvmr(r_input, gencov = 0)
  
  ## 3. Effects
  # Path a: RF >> metabolite
  a <- res_a$b;  se_a <- res_a$se;  pval_a <- res_a$pval
  # Path b: metabolite >> dementia, adjusted for RF
  b <- mvmr_res$b[mvmr_res$exposure == metab]
  se_b <- mvmr_res$se[mvmr_res$exposure == metab]
  pval_b <- mvmr_res$pval[mvmr_res$exposure == metab]
  # Direct effect: RF >> dementia, adjusted for metabolite
  direct <- mvmr_res$b[mvmr_res$exposure == rf]
  se_direct <- mvmr_res$se[mvmr_res$exposure == rf]
  pval_direct <- mvmr_res$pval[mvmr_res$exposure == rf]
  # Total effect: RF >> dementia
  total <- res_c$b;  se_c <- res_c$se;  pval_c <- res_c$pval
  
  # Indirect effect and proportion mediated (delta method)
  indirect <- a * b
  se_indirect <- sqrt((b^2 * se_a^2) + (a^2 * se_b^2))
  
  prop_med <- (a * b) / total
  se_prop_med <- sqrt((b^2 * se_a^2) / total^2 + ((a^2 * se_b^2) / total^2) + (((a * b)^2 * se_c^2) / total^4))
  
  res_paths <- rbind(
    res_paths,
    data.frame(Factor = rf, Metabolite = metab, Outcome = dem,
               Path = c("path_a", "path_b", "path_direct", "path_total"),
               Coef = c(a, b, direct, total),
               SE = c(se_a, se_b, se_direct, se_c),
               P_Value = c(pval_a, pval_b, pval_direct, pval_c))
  )
  
  res_mvmr_diag <- rbind(
    res_mvmr_diag,
    data.frame(Factor = rf, Metabolite = metab, Outcome = dem,
               nsnp = nrow(r_input),
               condF_RF = as.numeric(cond_F[1, which(exp_names == rf)]),
               condF_metab = as.numeric(cond_F[1, which(exp_names == metab)]),
               Q = q_res$Qstat,
               Q_pval = q_res$Qpval)
  )
  
  res_mediation <- rbind(
    res_mediation,
    data.frame(Factor = rf, Metabolite = metab, Outcome = dem,
               Path = c("indirect", "prop.med"),
               Coef = c(indirect, prop_med),
               SE = c(se_indirect, se_prop_med),
               Lower_95 = c(indirect - 1.96 * se_indirect, prop_med - 1.96 * se_prop_med),
               Upper_95 = c(indirect + 1.96 * se_indirect, prop_med + 1.96 * se_prop_med))
  )
}

write.csv(res_paths, paste0(out_dir, "path_MRmediation.csv"), row.names = F)
write.csv(res_mediation, paste0(out_dir, "res_MRmediation.csv"), row.names = F)
write.csv(res_mvmr_diag, paste0(out_dir, "MVMR_condF_Q.csv"), row.names = F)