## Statistical colocalisation: Glutamine and AD at rs2939302 (chr12:56,860,020, GRCh37)
## 1) metabolite vs outcome in a +/-250 kb window
## 2) blood eQTLs (eQTLGen) for candidate genes vs the metabolite and vs the outcome
## coloc.abf with default priors (p1 = p2 = 1e-4, p12 = 1e-5)

library(data.table)
library(dplyr)
library(coloc)

## ---- Paths ---------------------------------------------------------------------
out_dir            <- "results/coloc/"
gwas_exposure_file <- "GWAS_Sumstat/coloc/Glutamine_munged.tsv"   # SNP, CHR, BP, A1, A2, BETA, SE, P, FRQ, N
gwas_outcome_file  <- "GWAS_Sumstat/dementia/AD_munged.tsv"
eqtl_dir           <- "eQTLGen/"                                  # <gene>_munged.tsv: SNP, CHR, BP, A1, A2, P, FRQ, N

## ---- 1. Region -----------------------------------------------------------------
chr_of_interest <- 12
snp_bp          <- 56860020
region_start    <- snp_bp - 250000
region_end      <- snp_bp + 250000
region_label    <- sprintf("Chr%d:%d-%d", chr_of_interest, region_start, region_end)

## ---- 2. Load GWAS summary statistics for the region ----------------------------
GWAS_exposure <- fread(gwas_exposure_file) %>%
  filter(CHR == chr_of_interest, between(BP, region_start, region_end)) %>%
  mutate(MAF = ifelse(FRQ < 0.5, FRQ, 1 - FRQ)) %>%
  distinct(SNP, .keep_all = TRUE)

GWAS_outcome <- fread(gwas_outcome_file) %>%
  filter(CHR == chr_of_interest, between(BP, region_start, region_end)) %>%
  mutate(MAF = ifelse(FRQ < 0.5, FRQ, 1 - FRQ)) %>%
  distinct(SNP, .keep_all = TRUE)

# Align alleles of dataset y to dataset x (A2 = effect allele); flip beta_y if swapped
align_alleles <- function(df, effect_x, other_x, effect_y, other_y, beta_y = NULL) {
  same    <- df[[effect_x]] == df[[effect_y]] & df[[other_x]] == df[[other_y]]
  swapped <- df[[effect_x]] == df[[other_y]]  & df[[other_x]] == df[[effect_y]]
  if (!is.null(beta_y)) df[[beta_y]][swapped] <- -df[[beta_y]][swapped]
  df[same | swapped, ]
}

## ---- 3. Metabolite vs outcome ---------------------------------------------------
input <- merge(GWAS_exposure, GWAS_outcome, by = "SNP", suffixes = c("_exposure", "_outcome"))
input <- as.data.frame(input)
input <- align_alleles(input, "A2_exposure", "A1_exposure", "A2_outcome", "A1_outcome", "BETA_outcome")
nrow(input)

s_outcome <- 0.1   # proportion of cases in the AD GWAS

res_coloc <- coloc.abf(
  dataset1 = list(snp     = input$SNP,
                  type    = "quant",
                  beta    = input$BETA_exposure,
                  varbeta = input$SE_exposure^2,
                  MAF     = input$MAF_exposure,
                  N       = input$N_exposure),
  dataset2 = list(snp     = input$SNP,
                  type    = "cc",
                  beta    = input$BETA_outcome,
                  varbeta = input$SE_outcome^2,
                  MAF     = input$MAF_outcome,
                  N       = input$N_outcome,
                  s       = s_outcome)
)

coloc_main <- data.frame(region = region_label, t(res_coloc$summary))
write.csv(coloc_main, paste0(out_dir, "coloc_Glutamine_AD.csv"), row.names = FALSE)

## ---- 4. eQTL vs metabolite and eQTL vs outcome ---------------------------------
candidate_genes <- c("SPRYD4", "GLS2", "RBMS2", "STAT2", "SUOX")   # genes annotated to credible SNPs

eqtl_coloc <- data.frame()

for (gene in candidate_genes) {
  
  Gene <- fread(paste0(eqtl_dir, gene, "_munged.tsv")) %>%
    filter(CHR == chr_of_interest, between(BP, region_start, region_end)) %>%
    mutate(MAF = ifelse(FRQ < 0.5, FRQ, 1 - FRQ)) %>%
    distinct(SNP, .keep_all = TRUE)
  
  ## gene expression vs metabolite
  input_gene_exposure <- as.data.frame(merge(Gene, GWAS_exposure, by = "SNP", suffixes = c("_Gene", "_exposure")))
  input_gene_exposure <- align_alleles(input_gene_exposure, "A2_Gene", "A1_Gene", "A2_exposure", "A1_exposure", "BETA_exposure")
  
  result_exposure <- coloc.abf(
    dataset1 = list(snp     = input_gene_exposure$SNP,
                    type    = "quant",
                    pvalues = input_gene_exposure$P_Gene,
                    MAF     = input_gene_exposure$MAF_Gene,
                    N       = input_gene_exposure$N_Gene),
    dataset2 = list(snp     = input_gene_exposure$SNP,
                    type    = "quant",
                    beta    = input_gene_exposure$BETA_exposure,
                    varbeta = input_gene_exposure$SE_exposure^2,
                    MAF     = input_gene_exposure$MAF_exposure,
                    N       = input_gene_exposure$N_exposure)
  )
  
  ## gene expression vs outcome
  input_gene_outcome <- as.data.frame(merge(Gene, GWAS_outcome, by = "SNP", suffixes = c("_Gene", "_outcome")))
  input_gene_outcome <- align_alleles(input_gene_outcome, "A2_Gene", "A1_Gene", "A2_outcome", "A1_outcome", "BETA_outcome")
  
  result_outcome <- coloc.abf(
    dataset1 = list(snp     = input_gene_outcome$SNP,
                    type    = "quant",
                    pvalues = input_gene_outcome$P_Gene,
                    MAF     = input_gene_outcome$MAF_Gene,
                    N       = input_gene_outcome$N_Gene),
    dataset2 = list(snp     = input_gene_outcome$SNP,
                    type    = "cc",
                    beta    = input_gene_outcome$BETA_outcome,
                    varbeta = input_gene_outcome$SE_outcome^2,
                    MAF     = input_gene_outcome$MAF_outcome,
                    N       = input_gene_outcome$N_outcome,
                    s       = s_outcome)
  )
  
  eqtl_coloc <- rbind(
    eqtl_coloc,
    data.frame(gene = gene, region = region_label, trait = "Glutamine", t(result_exposure$summary)),
    data.frame(gene = gene, region = region_label, trait = "AD", t(result_outcome$summary))
  )
}

write.csv(eqtl_coloc, paste0(out_dir, "eqtl_coloc_Glutamine_AD.csv"), row.names = FALSE)