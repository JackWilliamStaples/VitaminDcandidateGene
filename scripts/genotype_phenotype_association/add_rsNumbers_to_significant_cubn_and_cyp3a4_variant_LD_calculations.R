# This is a script for adding rs-numbers to the variants in the significant_cubn_and_cyp3a4_variant_LD_calculations.txt file
# with the rs-ID lookup tale

# load necessary packages
source("../scripts/load_R_packages.R")

# load the significant_cubn_and_cyp3a4_variant_LD_calculations.txt file
significant_cubn_and_cyp3a4_variant_LD_calculations <-
  read.table(
    file = "./significant_cubn_and_cyp3a4_variant_LD_calculations.txt",
    header = TRUE
    )

# load the dbSNP rsID lookup table
rsIDlookupTable <-
  read.delim(
    file = "./PLINKvariantID_dbSNP_rsID_lookupTable_commonVariants.txt",
    header = TRUE
      )

significant_cubn_and_cyp3a4_variant_LD_calculations <-
# add rs-numbers for the SNP_A column
significant_cubn_and_cyp3a4_variant_LD_calculations %>%
  left_join(
    x = .,
    y = rsIDlookupTable %>%
        select(
          .data = .,
          SNP,
          existing_variant_VEP
        ),
    by = c("SNP_A" = "SNP")
      ) %>%
  rename(
    .data = .,
    "rs_A" = existing_variant_VEP
    ) %>%
# add rs-numbers for the SNP_B column
  left_join(
    x = .,
    y = rsIDlookupTable %>%
        select(
          .data = .,
          SNP,
          existing_variant_VEP
        ),
    by = c("SNP_B" = "SNP")
      ) %>%
  rename(
    .data = .,
    "rs_B" = existing_variant_VEP
    ) %>%
  # reorder the columns
  select(
    .data = .,
    CHR_A:SNP_A,
    rs_A,
    CHR_B:SNP_B,
    rs_B,
    R2
    )

# save the resulting file to the file system
significant_cubn_and_cyp3a4_variant_LD_calculations %>%
  write_tsv(
    x = .,
    file = "./significant_cubn_and_cyp3a4_variant_LD_calculations.txt"
    )


