# load necessary packages
source("./scripts/load_R_packages.R")

# read in the out.diff.sites_in_files detailing differences and similarities between batch 1 and batch 2 VCF files
Site_differences <-
read.table(
  file = "./out.diff.sites_in_files",
  header = TRUE
  )

# create a file with the sites that are only contained in batch 1 and the non-matching overlapping sites
Site_differences %>%
  # filter to sites in file 1 only
  filter(
    .data = .,
    (IN_FILE == "1") | (IN_FILE == "O")
    ) %>%
  # create a column of variant IDs with the format with the format CHROM:POS:REF:ALT
  # using the columns that denote the variants in batch 1 (e.g., POS1 and REF1)
  mutate(
    .data = .,
    VariantID = paste0(
                        CHROM,":",
                        POS1,":",
                        REF1,":",
                        ALT1
                      )
    ) %>%
  # select the variant ID column
  select(
    .data = .,
    VariantID
    ) %>%
  # save the resulting file to the file system without a header
  write_tsv(
    x = .,
    file = "batch1_sites_for_filtering.txt",
    col_names = FALSE
  )

# create a file with the sites that are only contained in batch 2 and the non-matching overlapping sites
Site_differences %>%
  # filter to sites in file 2 only
  filter(
    .data = .,
    (IN_FILE == "2") | (IN_FILE == "O")
  ) %>%
  # create a column of variant IDs with the format with the format CHROM:POS:REF:ALT
  # using the columns that denote the variants in batch 2 (e.g., POS2 and REF2 )
  mutate(
    .data = .,
    VariantID = paste0(
      CHROM,":",
      POS2,":",
      REF2,":",
      ALT2
          )
  ) %>%
  # select the variant ID column
  select(
    .data = .,
    VariantID
  ) %>%
  # save the resulting file to the file system without a header
  write_tsv(
    x = .,
    file = "batch2_sites_for_filtering.txt",
    col_names = FALSE
  )

# create a file with variants that are in BOTH batches 1 and 2
Site_differences %>%
  filter(
    .data = .,
    IN_FILE == "B"
    ) %>%
  # create a variant ID column with the format CHROM:POS:REF:ALT
  # using the columns that denote batch1 in case the POS1,REF1,ALT1 are not the exact same as the POS2,REF2,ALT2
  mutate(
    .data = .,
    VariantID = paste0(
                        CHROM,":",
                        POS1,":",
                        REF1,":",
                        ALT1
                      )
    ) %>%
  # select the variant ID column
  select(
    .data = .,
    VariantID
  ) %>%
  # save the resulting file to the file system without a header
  write_tsv(
    x = .,
    file = "shared_sites_between_batch1_and_2_batch1_variantIDs.txt",
    col_names = FALSE
  )

# create a file with variants that are in BOTH batches 1 and 2
Site_differences %>%
  filter(
    .data = .,
    IN_FILE == "B"
    ) %>%
  # create a variant ID column with the format CHROM:POS:REF:ALT
  # using either the columns that denote batch2 in case the POS2,REF2,ALT2 are not the exact same as the POS1,REF1,ALT1
  mutate(
    .data = .,
    VariantID = paste0(
                        CHROM,":",
                        POS2,":",
                        REF2,":",
                        ALT2
                      )
    ) %>%
  # select the variant ID column
  select(
    .data = .,
    VariantID
  ) %>%
  # save the resulting file to the file system without a header
  write_tsv(
    x = .,
    file = "shared_sites_between_batch1_and_2_batch2_variantIDs.txt",
    col_names = FALSE
  )

print(x = "batch1_sites_for_filtering.txt file created")
print(x = "batch2_sites_for_filtering.txt file created")
print(x = "shared_sites_between_batch1_and_2_batch1_variantIDs.txt file created")
print(x = "shared_sites_between_batch1_and_2_batch2_variantIDs.txt file created")
