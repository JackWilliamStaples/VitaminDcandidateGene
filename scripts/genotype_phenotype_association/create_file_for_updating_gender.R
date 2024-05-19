# this is a script for updating the gender of the study participants in the PLINK genotype file
# load necessary packages
source("../scripts/load_R_packages.R")


# read in the sample lookup tables from batch1 and batch2 of sequencing
batch1_lookup_table <-
read_excel(path = "../sample_info_batch_1/lookup_woodahl_grc_1_updated.xlsx")
batch2_lookup_table <-
  read_excel(path = "../sample_info_batch_2/lookup_woodahl_grc_custom_2.xlsx")

LookupTableOfAllSamples <-
# bind the lookup tables from batch 1 and batch 2 by row
rbind(
  batch1_lookup_table,
  batch2_lookup_table
  ) %>%
  # select relevant columns
  select(
    .data = .,
    "QualityControlNote" = Note,
    "VCFfileID" = `NWGC LIMS ID`,
    "VisitCode" = `UDF/Investigator Sample Name`,
    "Gender" = `UDF/Gender`
    ) %>%
  # filter to samples that pass QC, i.e., a quality control note of "Released"
  filter(
    .data = .,
    QualityControlNote == "Released"
    ) %>%
  # remove the duplicate sample (visit code A-313) that was removed during sequencing quality control
  filter(
    .data = .,
    VisitCode != "A-313"
    )

# create a file that can be used to update the gender of the study participants in the PLINK files
ParticipantGenderPLINKfile <-
  LookupTableOfAllSamples %>%
  # convert the gender column to 1 for male and 2 for female
  mutate(
    .data = .,
    Gender = if_else(
                condition = Gender == "Male",
                true = 1,
                false = 2
                  )
    ) %>%
  select(
    .data = .,
    "FID" = VCFfileID,
    "IID" = VCFfileID,
    Gender
    )

# save the participant gender plink file to the file system
write_tsv(
  x = ParticipantGenderPLINKfile,
  file = "./ParticipantGenderPLINKfile.txt",
  col_names = FALSE
  )