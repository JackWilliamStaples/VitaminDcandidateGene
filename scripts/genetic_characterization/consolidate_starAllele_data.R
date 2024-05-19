# this is a script for consolidating the star allele haplotype, diplotype, and phenotype calling results 
# from starGazer, PGxPOP, and custom star allele calling into tables for manuscripts

# load necessary packages
source("./scripts/load_R_packages.R")

# load star allele data from starGazer software for CYP2R1 and UGT1A4
CYP2R1haplotypeData <-
  read.delim(
    file = "./starGazer/cyp2r1_HaplotypeFrequencyTable.tsv",
    header = TRUE
      )

CYP2R1diplotypeData <-
  read.delim(
    file = "./starGazer/cyp2r1_DiplotypeFrequencyTable.tsv",
    header = TRUE
    )

UGT1A4haplotypeData <-
  read.delim(
    file = "./starGazer/ugt1a4_HaplotypeFrequencyTable.tsv",
    header = TRUE
      )

UGT1A4diplotypeData <-
  read.delim(
    file = "./starGazer/ugt1a4_DiplotypeFrequencyTable.tsv",
    header = TRUE
    )

# load star allele data for CYP3A4, called using a custom method
CYP3A4haplotypeData <-
  read.delim(
    file = "./starGazer/custom_cyp3A4starAllele_calls/HaplotypeFrequencyTable.tsv",
    header = TRUE
    )

CYP3A4diplotypeData <-
  read.delim(
    file = "./starGazer/custom_cyp3A4starAllele_calls/DiplotypeFrequencyTable.tsv",
    header = TRUE
    )

# load star allele data for UGT1A1, called with PGx-POP software
UGT1A1haplotypeData <-
  read.delim(
    file = "./pgx_pop_annotation/ugt1a1_HaplotypeFrequencyTable.tsv",
    header = TRUE
    )

UGT1A1diplotypeData <-
  read.delim(
    file = "./pgx_pop_annotation/ugt1a1_DiplotypeFrequencyTable.tsv",
    header = TRUE
    )

# bind all haplotype tables together by row
StarAlleleHaplotypeTableToSave <-
  list(
    CYP2R1haplotypeData,
    CYP3A4haplotypeData,
    UGT1A1haplotypeData,
    UGT1A4haplotypeData
  ) %>%
  do.call(
    what = "rbind",
    args = .
    ) %>%
  # update frequency column name
  select(
    .data = .,
    everything(),
    "Frequency (N = 938)" = Haplotype.Frequency..N...938.
    )

# bind all diplotype tables together by row
StarAlleleDiplotypeTableToSave <-
  list(
    CYP2R1diplotypeData,
    CYP3A4diplotypeData,
    UGT1A1diplotypeData,
    UGT1A4diplotypeData
  ) %>%
  do.call(
    what = "rbind",
    args = .
    ) %>%
  # update the frequency column name
  select(
    .data = .,
    everything(),
    "Frequency (N = 469)" = Frequency..N...469.
  )

# make a directory for storing the results if it does not exist
if(!dir.exists(paths = "./starAlleleCallingResults/"))
{
  dir.create(path = "./starAlleleCallingResults/")
}

# save to the file system
StarAlleleHaplotypeTableToSave %>%
  write_tsv(
    x = .,
    file = "./starAlleleCallingResults/StarAlleleHaplotypeTableToSave.tsv"
      )

StarAlleleDiplotypeTableToSave %>%
  write_tsv(
    x = .,
    file = "./starAlleleCallingResults/StarAlleleDiplotypeTableToSave.tsv"
    )

print(x = "################## Star allele haplotype and diplotype tables saved ################################3")

