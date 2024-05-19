# load necessary packages
source("../scripts/load_R_packages.R")

# read in the report comparing genotype missingness for individual snps sequenced in batches 1 and 2 of the vitamin D data
batch_effect_allele_missingness_report <-
  read.table(
    file = "vitamin_D_merged_plink_case_control_to_check_for_missingness.missing",
    stringsAsFactors = FALSE,
    header = TRUE
    )

# obtain a list of SNP ids that have a P-value of less than 0.05 
# for the null hypothesis that there is no difference in genotype missingness for the 
# individual SNP sequenced in batch 1 & 2 based on Fisher's exact test
SNPsWithSignificantDifferencesInGenotypeMissingness <-
batch_effect_allele_missingness_report %>%
  filter(
    .data = .,
    as.numeric(x = P) < 0.05
    ) %>%
  select(
    .data = .,
    SNP
    )
#save the file to the file system
write_tsv(
  x = SNPsWithSignificantDifferencesInGenotypeMissingness,
  file = "SNPsWithSignificantDifferencesInGenotypeMissingness.txt",
  col_names = FALSE
  )