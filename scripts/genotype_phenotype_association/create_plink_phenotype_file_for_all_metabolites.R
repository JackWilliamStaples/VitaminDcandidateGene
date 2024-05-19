# This is a script for creating phenotype files for association of all vitamin D candidate gene snps with vitamin D and vitamin D metabolite phenotypes
# This script will also create a missing value report and identify participants with metabolite quantitation duplicates
# This will also create a file with demographic and seasonal covariate data for adjusting regressions
# This will also create files of binary and continuous phenotypes for rare variant analysis
# This will create histograms and box and whisker plots for each metabolite to assess deviations from normality and identify outlying observations
# This will also perform Shapiro-Wilke's tests for deviations from normality on untransformed and log10() transformed metabolite values
# This will create correlation plots of all pairwise D3 and D2 metabolite correlations
# This will also test for heteroscedasticity for appropriate D3 and D2 metabolite correlations

# load necessary packages
source("../scripts/load_R_packages.R")

# load genetic-association analysis functions
source("../scripts/genotype_phenotype_association/genetic_association_analysis_functions.R")

# read in the metabolite data from the first round of metabolite quantitation
metaboliteData <-
  read.delim(
    file = "../metabolite_data/Vit_D_Master_clean.csv",
    sep = ","
  ) %>%
  # change the name of the participant ID column to VisitCode since this is a column of VisitCodes
  select(
    .data = .,
    "VisitCode" = ParticipantID,
    everything()
    ) %>%
  # remove any starting or trailing whitespace from each column in the table,
  # if the entry in the dataframe is an empty string, convert it to an NA
  RemoveWhiteSpaceAndChangeEmptyStringsToNA(
    DataFrameToUse = .
  )

# read in the metabolite data from the second round of metabolite quantitation
metaboliteData_round2 <-
  read.delim(
    file = "../metabolite_data/UM_Vit_D_data_summary_A417_toA560_clean.csv",
    sep = ","
    ) %>%
  # change the name of the participant ID column to VisitCode since this is a column of VisitCodes
  select(
    .data = .,
    "VisitCode" = ParticipantID,
    everything()
    ) %>%
  # remove any starting or trailing whitespace from each column in the table,
  # if the entry in the dataframe is an empty string, convert it to an NA
  RemoveWhiteSpaceAndChangeEmptyStringsToNA(
    DataFrameToUse = .
  ) %>%
  # create a column filled with N/As for the 1,25(OH)2D2 metabolite since it was not quantitated in this round
  mutate(
    .data = .,
    active_1alpha25OH2D2_ngperml = as.character(x = "N/A")
    )

# read in the de-identified demographic data with information of study date, age, gender, height, and weight
EncountersData <-
read.csv(
  file = "../covariate_and_demographic_data_deidentified/Encounters-Grid view.csv",
  header = TRUE
    ) %>%
  # remove any starting or trailing whitespace from each column in the table,
  # if the entry in the dataframe is an empty string, convert it to an NA
  RemoveWhiteSpaceAndChangeEmptyStringsToNA(
    DataFrameToUse = .
    )

ParticipantData <-
  # load de-identified participant data
  read.csv(
    file = "../covariate_and_demographic_data_deidentified/Participant info-Grid view_deidentified.csv",
    header = TRUE 
      ) %>%
  # remove any starting or trailing whitespace from each column in the table,
  # if the entry in the dataframe is an empty string, convert it to an NA
  RemoveWhiteSpaceAndChangeEmptyStringsToNA(
    DataFrameToUse = .
  )

# read in the sample lookup tables from batch1 and batch2 of sequencing, 
  # remove any leading/trailing whitespace for each entry in the tables and convert empty string to missing
batch1_lookup_table <-
  read_excel(path = "../sample_info_batch_1/lookup_woodahl_grc_1_updated.xlsx") %>%
  RemoveWhiteSpaceAndChangeEmptyStringsToNA(DataFrameToUse = .)

batch2_lookup_table <-
  read_excel(path = "../sample_info_batch_2/lookup_woodahl_grc_custom_2.xlsx") %>%
  RemoveWhiteSpaceAndChangeEmptyStringsToNA(DataFrameToUse = .)

# create a table of visit code and VCF sample IDs to convert the visit codes in the metabolite data spreadsheet to sequencing sample IDs
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
    "VCFfileID" = NWGC.LIMS.ID,
    "VisitCode" = UDF.Investigator.Sample.Name,
    "Gender" = UDF.Gender
  ) %>%
  # filter to samples that pass QC ---- a quality control note of "Released"
  filter(
    .data = .,
    QualityControlNote == "Released"
  ) %>%
  # remove the duplicate sample (visit code A-313) that was removed during sequencing quality control
  filter(
    .data = .,
    VisitCode != "A-313"
  ) %>%
  select(
    .data = .,
    VCFfileID,
    VisitCode
    )

################################ code for cleaning the metabolite data ##############################################

# save the names of the columns from round 1 of metabolite quantitation with metabolite concentrations 
# as an array by identifying the nanogram per mL and picogram per mL suffixes
ColumnsWithMetaboliteConcentrations <-
  metaboliteData %>%
  names(x = .) %>%
  grepl(pattern = "(ngperml|pgperml)",x = .) %>%
  which(x = .) %>%
  names(x = metaboliteData)[.]

# create the final PLINK phenotype file from the first round of metabolite quantitation 
VitaminDplinkPhenotypeFile_round1 <-
  metaboliteData %>%
  # clean up all of the metabolite data columns:
    # change all of the BLQs to NA,
    # change the NDs to NA
    # change any values with a "<" symbol to NA
    # change anything with a N/A to NA
  mutate_at(
    .tbl = .,
    .vars = ColumnsWithMetaboliteConcentrations,
    .funs = DataCleaningFunctionToApplyToAllMetaboliteColumns
  ) %>%
  # convert the study date column to a character vector
  mutate(
    .data = .,
    StudyDate = StudyDate %>% as.character(x = .)
  ) %>%
  # convert the units of the columns that have concentrations in picograms per ml
  # to concentrations of nanogram per ml
  mutate_at(
    .tbl = .,
    .vars = vars(contains(match = "pgperml")),
    .funs = list(ConvertPicoGramsToNanoGrams = ConvertPicoGramsToNanoGrams)
  ) %>%
  # change the suffix of the column names that have been converted to units of nanograms per mL to ngperml
  # and deselect the old pgperml columns
  select(
    .data = .,
    everything(),
    -beta_4beta25OH2D3_pgperml,
    -active_1alpha25OH2D3_pgperml,
    -VitaminD2_pgperml,
    -primary_25OHD2_pgperml,
    -active_1alpha25OH2D2_pgperml,
    "beta_4beta25OH2D3_ngperml" = beta_4beta25OH2D3_pgperml_ConvertPicoGramsToNanoGrams,
    "active_1alpha25OH2D3_ngperml" = active_1alpha25OH2D3_pgperml_ConvertPicoGramsToNanoGrams,
    "VitaminD2_ngperml" = VitaminD2_pgperml_ConvertPicoGramsToNanoGrams,
    "primary_25OHD2_ngperml" = primary_25OHD2_pgperml_ConvertPicoGramsToNanoGrams,
    "active_1alpha25OH2D2_ngperml" = active_1alpha25OH2D2_pgperml_ConvertPicoGramsToNanoGrams
  ) %>%
  # create a column where the 25OHD3 and 25OHD2 are summed together
  # to get a measure of the patients' overall 25(OH)D status:
      # 25(OH)D with no subscript is 25(OH)D3 + 25(OH)D2.
  # A patient’s overall 25(OH)D status is the sum of 25(OH)D3 and 25(OH)D2:
      # (A) 25(OH)D deficiency is < 12 ng/mL,
      # (B) insufficiency is between 12 ng/mL and 20 ng/mL,
      # (C) and sufficiency is 20 ng/mL based on Institute of Medicine cutoffs,
          # but there is disagreement on sufficiency… some groups say 30 ng/mL.
  mutate(
    .data = .,
    Overall_primary_25OHD_ngperml = mapply(
                                           FUN = function(currentD3value,currentD2value)
                                           {
                                             ValueToReturn <-
                                               c(
                                                 currentD3value,
                                                 currentD2value
                                                 ) %>%
                                               sum(
                                                 .,
                                                 na.rm = TRUE
                                                 )
                                             return(ValueToReturn)
                                           },
                                           primary_25OHD3_ngperml,
                                           primary_25OHD2_ngperml,
                                           SIMPLIFY = FALSE
                                           ) %>%
                                          unlist(x = .)
  ) %>%
  # determine whether the patient is "sufficient", "insufficient", or "deficient"
  # based on the cutoffs defined above in (A) - (C)
  mutate(
    .data = .,
    OverallVitDstatus = if_else(
                                # if the value is less than 12, code as deficient
                                condition = Overall_primary_25OHD_ngperml < 12,
                                true = "Deficient (< 12 ng/mL)",
                                false = if_else(
                                                # if the value is between 12 and 20, code as insufficient
                                                condition = Overall_primary_25OHD_ngperml >= 12 & Overall_primary_25OHD_ngperml <= 20,
                                                true = "Insufficient (between 12 & 20 ng/mL)",
                                                false = if_else(
                                                                # if the value is greater than 50, code as high
                                                                condition = Overall_primary_25OHD_ngperml > 50,
                                                                true = "High (> 50 ng/mL)",
                                                                # otherwise code as sufficient
                                                                false = "Sufficient (between 20 & 50 ng/mL)"
                                                                )
                                                )
                                )
        ) %>%
    # create a column with a 1=yes for sufficiency or a 0=no for below sufficiency 
    # based on the 20ng/ml cutoff if the Overall_primary_25OHD_ngperml is non-missing
        # this phenotype can be used in logistic regressions to test variants associated with sufficiency/insufficiency
    mutate(
      .data = .,
      VitDsufficiency = if_else(
                                condition = Overall_primary_25OHD_ngperml > 20,
                                true = as.numeric(x = 1),
                                false = if_else(
                                                condition = is.na(x = Overall_primary_25OHD_ngperml),
                                                true = as.numeric(x = NA),
                                                false = as.numeric(x = 0)
                                                )
                                ),
      # create a column with a 1=yes for deficiency or a 0=no for above deficiency 
      # based on a 12ng/ml cutoff if the Overall_primary_25OHD_ngperml is non-missing
      VitDdeficiency = if_else(
                                condition = Overall_primary_25OHD_ngperml < 12,
                                true = as.numeric(x = 1),
                                false = if_else(
                                                condition = is.na(x = Overall_primary_25OHD_ngperml),
                                                true = as.numeric(x = NA),
                                                false = as.numeric(x = 0)
                                              )
                              )
    ) %>%
    # create columns with metabolite ratios that are relevant to the vitamin D metabolic pathway
    mutate(
      .data = .,
      primary_25OHD3_VitaminD3_ratio = primary_25OHD3_ngperml/VitaminD3_ngperml,
      active_1alpha25OH2D3_primary_25OHD3_ratio = active_1alpha25OH2D3_ngperml/primary_25OHD3_ngperml,
      R_24R25OH2D3_primary_25OHD3_ratio = R_24R25OH2D3_ngperml/primary_25OHD3_ngperml,
      beta_4beta25OH2D3_primary_25OHD3_ratio = beta_4beta25OH2D3_ngperml/primary_25OHD3_ngperml,
      primary_25OHD2_VitaminD2_ratio = primary_25OHD2_ngperml/VitaminD2_ngperml,
      active_1alpha25OH2D2_primary_25OHD2_ratio = active_1alpha25OH2D2_ngperml/primary_25OHD2_ngperml
      ) %>%
    # replace all of the NA's in every row with a -9 (the missing phenotype code for PLINK)
    mutate_all(
      .tbl = .,
      .funs = ~ tidyr::replace_na(
                    data = .x,
                    replace = -9
                    )
      ) %>%
    # deselect unnecessary columns
    select(
        .data = .,
        -AltParticipantID,
        -StudyDate,
        -AnalysisDate
        )

# save the names of the columns with metabolite concentrations from the second round of metabolite quantitation 
# as an array by identifying the nanogram per mL and picogram per mL suffixes
ColumnsWithMetaboliteConcentrations <-
  metaboliteData_round2 %>%
  names(x = .) %>%
  grepl(pattern = "(ngperml|pgperml)",x = .) %>%
  which(x = .) %>%
  names(x = metaboliteData_round2)[.]

# create a vitamin D phenotype file from round 2 of metabolite quantitation
VitaminDplinkPhenotypeFile_round2 <-
  metaboliteData_round2 %>%
  # clean up all of the metabolite data columns:
    # change all of the BLQs to NA,
    # change the NDs to NA
    # change any less than (<) signs to an NA
    # change anything with a N/A to NA
  mutate_at(
    .tbl = .,
    .vars = ColumnsWithMetaboliteConcentrations,
    .funs = DataCleaningFunctionToApplyToAllMetaboliteColumns
  ) %>%
  # create a column where the 25OHD3 and 25OHD2 are summed together
  # to get a measure of the patients' overall 25(OH)D status
  # 25(OH)D with no subscript is 25(OH)D3 + 25(OH)D2.
  # A patient’s overall 25(OH)D status is the sum of 25(OH)D3 and 25(OH)D2:
      # (A) 25(OH)D deficiency is < 12 ng/mL,
      # (B) insufficiency is between 12 ng/mL and 20 ng/mL,
      # (C) and sufficiency is 20 ng/mL based on Institute of Medicine cutoffs,
          # but there is disagreement on sufficiency… some groups say 30 ng/mL.
  mutate(
    .data = .,
    Overall_primary_25OHD_ngperml = mapply(
                                           FUN = function(currentD3value,currentD2value)
                                           {
                                             ValueToReturn <-
                                               c(
                                                 currentD3value,
                                                 currentD2value
                                                 ) %>%
                                               sum(
                                                 .,
                                                 na.rm = TRUE
                                                 )
                                             return(ValueToReturn)
                                           },
                                           primary_25OHD3_ngperml,
                                           primary_25OHD2_ngperml,
                                           SIMPLIFY = FALSE
                                           ) %>%
                                          unlist(x = .)
  ) %>%
  # determine whether the patient is "sufficient", "insufficient", or "deficient"
  # based on the cutoffs defined above in (A) - (C)
  mutate(
    .data = .,
    OverallVitDstatus = if_else(
                                # if the value is less than 12, code as deficient
                                condition = Overall_primary_25OHD_ngperml < 12,
                                true = "Deficient (< 12 ng/mL)",
                                false = if_else(
                                                # if the value is between 12 and 20, code as insufficient
                                                condition = Overall_primary_25OHD_ngperml >= 12 & Overall_primary_25OHD_ngperml <= 20,
                                                true = "Insufficient (between 12 & 20 ng/mL)",
                                                false = if_else(
                                                                # if the value is greater than 50, code as high
                                                                condition = Overall_primary_25OHD_ngperml > 50,
                                                                true = "High (> 50 ng/mL)",
                                                                # otherwise code as sufficient
                                                                false = "Sufficient (between 20 & 50 ng/mL)"
                                                                )
                                                )
                                )
        ) %>%
    # create a column with a 1=yes for sufficiency or a 0=no for below sufficiency based on the 20ng/ml cutoff if the Overall_primary_25OHD_ngperml is non-missing
        # this phenotype can be used in logistic regressions to test variants associated with sufficiency/insufficiency
    mutate(
      .data = .,
      VitDsufficiency = if_else(
                                condition = Overall_primary_25OHD_ngperml > 20,
                                true = as.numeric(x = 1),
                                false = if_else(
                                                condition = is.na(x = Overall_primary_25OHD_ngperml),
                                                true = as.numeric(x = NA),
                                                false = as.numeric(x = 0)
                                                )
                                ),
      # create a column with a 1=yes for deficiency or a 0=no for above deficiency based on a 12ng/ml cutoff if the Overall_primary_25OHD_ngperml is non-missing
      VitDdeficiency = if_else(
                                condition = Overall_primary_25OHD_ngperml < 12,
                                true = as.numeric(x = 1),
                                false = if_else(
                                                condition = is.na(x = Overall_primary_25OHD_ngperml),
                                                true = as.numeric(x = NA),
                                                false = as.numeric(x = 0)
                                              )
                              )
    ) %>%
    # create columns with metabolite ratios that are relevant to the vitamin D metabolic pathway
    mutate(
      .data = .,
      primary_25OHD3_VitaminD3_ratio = primary_25OHD3_ngperml/VitaminD3_ngperml,
      active_1alpha25OH2D3_primary_25OHD3_ratio = active_1alpha25OH2D3_ngperml/primary_25OHD3_ngperml,
      R_24R25OH2D3_primary_25OHD3_ratio = R_24R25OH2D3_ngperml/primary_25OHD3_ngperml,
      beta_4beta25OH2D3_primary_25OHD3_ratio = beta_4beta25OH2D3_ngperml/primary_25OHD3_ngperml,
      primary_25OHD2_VitaminD2_ratio = primary_25OHD2_ngperml/VitaminD2_ngperml,
      active_1alpha25OH2D2_primary_25OHD2_ratio = active_1alpha25OH2D2_ngperml/primary_25OHD2_ngperml
      ) %>%
    # replace all of the NA's in every row with a -9, the missing phenotype code for PLINK
    mutate_all(
      .tbl = .,
      .funs = ~ tidyr::replace_na(
                    data = .x,
                    replace = -9
                    )
      ) %>%
   # deselect unneccesary columns
   select(
     .data = .,
     -AssayDate,
     -AssayNumber,
     -volume_assayed_ml
     )

# bind the vitamin D plink phenotype files from round 1 and round 2 of metabolite quantitation into one file by row
# with quantitation duplicates included
VitaminDplinkPhenotypeFileWithDuplicates <-
  rbind(
    VitaminDplinkPhenotypeFile_round1,
    VitaminDplinkPhenotypeFile_round2
  ) 

### identify the duplicate samples in the VitaminDplinkPhenotypeFileWithDuplicates based on duplicate participant IDs
ParticipantIDsWithMetaboliteQuantitationDuplicates <-
  VitaminDplinkPhenotypeFileWithDuplicates %>%
  # join the VitaminDplinkPhenotypeFileWithDuplicates with the encounters data to obtain participants IDs
  left_join(
    x = .,
    y = EncountersData,
    by = c("VisitCode" = "Visit.code")
  ) %>%
  # select the Participant.ID column and convert to an array
  pull(
    .data = .,
    Participant.ID
    ) %>%
  # identify duplicated participant IDs
  duplicated(x = .) %>%
  # convert the duplicated participant IDs to row indexes
  which(x = .) %>%
  # subset the rows of the VitaminDplinkPhenotypeFileWithDuplicates with the duplicated participant IDs
  VitaminDplinkPhenotypeFileWithDuplicates[.,] %>%
  # join the visit codes that have a duplicate with the EncountersData to obtain participant IDs corresponding to
  # the visit code
  left_join(
    x = .,
    y = EncountersData,
    by = c("VisitCode" = "Visit.code")
  ) %>%
  pull(
    .data = .,
    Participant.ID
    ) %>%
  unique(x = .)

# print a message about the number of participants who participated in the study more than once
paste(
  "There are",
  length(x = ParticipantIDsWithMetaboliteQuantitationDuplicates),
  "unique participants that participated in the vitamin D study more than once out of",
  nrow(x = VitaminDplinkPhenotypeFileWithDuplicates),
  "total metabolite samples that were quantitated"
  ) %>%
print(x = .)

DuplicateMetaboliteQuantitationEncounters <-
# filter the encounters data to the participant IDs that are ParticipantIDsWithMetaboliteQuantitationDuplicates
EncountersData %>%
  filter(
    .data = .,
    Participant.ID %in% ParticipantIDsWithMetaboliteQuantitationDuplicates
    ) %>%
  # arrange the data by participant ID so the duplicates are next to each other
  arrange(
    .data = .,
    Participant.ID
    ) %>%
  # join the duplicate encounters with the metabolite quantitation data based on visit code
  left_join(
    x = .,
    y = VitaminDplinkPhenotypeFileWithDuplicates,
    by = c("Visit.code" = "VisitCode")
    ) %>%
  # filter out samples that were not quantitated at all
  filter(
    .data = .,
    Study != "Low draw P01 - no usable sample"
    )

# print a message with the number of duplicate encounters
paste(
  "There were a total of ",
  nrow(x = DuplicateMetaboliteQuantitationEncounters),
  "duplicate vitamin D quantitation encounters for",
  length(x = ParticipantIDsWithMetaboliteQuantitationDuplicates),
  "unique participants that were quantitated more than once"
  ) %>%
  print(x = .)

# save the duplicate metabolite quantitation encounters to the file system
if(!dir.exists("../DuplicateMetaboliteQuantitationEncounters/"))
{
  dir.create("../DuplicateMetaboliteQuantitationEncounters/")
  
  DuplicateMetaboliteQuantitationEncounters %>%
    write_tsv(
      x = .,
      file = "../DuplicateMetaboliteQuantitationEncounters/DuplicateMetaboliteQuantitationEncounters.txt"
      )
}

#### identify the participant IDs and the visit codes of the DuplicateMetaboliteQuantitationEncounters that were quantitated on the latest
#### date and earliest dates for each duplicate participant
DuplicateMetaboliteQuantitationEncountersStudyCodes <-
  DuplicateMetaboliteQuantitationEncounters %>%
  # convert the study date column to an R date data type with the format month/day/year
  mutate(
    .data = .,
    Study.Date = Study.Date %>%
                 as.Date(
                   x = .,
                   '%m/%d/%y'
                   )
    ) %>%
  # group by the participant ID
  group_by(
    .data = .,
    Participant.ID
    ) %>%
  # split into a list of dataframes by participant ID
  group_split(.tbl = .) %>%
  # loop through each dataframe
  lapply(
    X = .,
    FUN = function(dataFrameForCurrentParticipantID)
    {
      # identify the earliest and latest visit code based on the date for the current participant ID
        DataToReturn <-
            dataFrameForCurrentParticipantID %>%
            # filter to the the latest study date (the max date) for the participant ID
            filter(
              .data = .,
              Study.Date == max(Study.Date,na.rm = TRUE)
              ) %>%
            # change the name of the visit code column to Visit.code_laterDate and select the Participant.ID
            select(
              .data = .,
              "Visit.code__laterDate" = Visit.code,
              Participant.ID
            ) %>%
            # add a column with the visit code from the earliest date
            mutate(
              .data = .,
              Visit.code__EarliestDate = dataFrameForCurrentParticipantID %>%
                                          # filter to the the earliest study date (the min date) for the participant ID
                                          filter(
                                            .data = .,
                                            Study.Date == min(Study.Date,na.rm = TRUE)
                                          ) %>%
                                          # pull out the visit code for the earliest date
                                          pull(
                                            .data = .,
                                            Visit.code
                                            )
              )
        
        return(DataToReturn)
    }
      ) %>%
  # bind the list of dataframes back together by row
  do.call(
    what = "rbind",
    args = .
    )

#########  identify the specific missing value symbols (ND, BLQ, blq, N/A, and values with a "<" symbol) for each metabolite in the quantitation data 
#########  from round 1 and round 2 of quantitation and compute the number of appearances of each missing symbol for all the samples that were quantitated (N = 532)

  MissingValueSymbolSummaries <-
      # loop through the metabolite quantitation data from rounds 1 and 2
      list(
        "round1" = metaboliteData,
        "round2" = metaboliteData_round2
      ) %>%
      lapply(
        X = .,
        FUN = function(currentDataSet)
        {
          DataToReturn <-
            currentDataSet %>%
            # select metabolite value columns only based on the ngperml and pgperml suffixes
            names(x = .) %>%
            grepl(
              pattern = "(ngperml|pgperml)",
              x = .
              ) %>%
            which(x = .) %>%
            names(x = currentDataSet)[.] %>%
            sort(x = .) %>%
            currentDataSet[.] %>%
            # loop through each metabolite column
            lapply(
              X = .,
              FUN = function(currentMetaboliteValues)
                {

                   # count the number of occurances of each missing value symbol (ND,blq,BLQ,N/A, or a < sign)
                   MissingValueSummary <-
                      # create a one column dataframe with the metabolite values
                      data.frame(
                        "Metabolite_value" = currentMetaboliteValues
                      ) %>%
                      # filter to rows that contain missing values (indicated by ND,blq,BLQ,N/A, or a < sign)
                      filter(
                        .data = .,
                        grepl(
                          pattern = "(BLQ|blq|ND|N/A|<)",
                          x = Metabolite_value
                            )
                        ) %>%
                      # if the metabolite value contains a "blq", convert it to upper-case,
                        # if the metabolite value contains a "<" change the value to a "<" symbol only
                      mutate(
                        .data = .,
                        Metabolite_value = if_else(
                                                   condition = Metabolite_value == "blq",
                                                   true = toupper(x = Metabolite_value),
                                                   false = if_else(
                                                                  condition = grepl(pattern = "<",x = Metabolite_value),
                                                                  true = "<",
                                                                  false = Metabolite_value
                                                                    )
                                                     )
                        ) %>%
                       # group the results by the missing metabolite value symbol
                      group_by(
                        .data = .,
                        Metabolite_value
                        ) %>%
                       # count the number of occurances of each missing metabolite value symbol
                      summarise(
                        .data = .,
                        Symbol_Count = n()
                        )

                   return(MissingValueSummary)

                }
              )

          # remove the units from the names of the metabolites since they do not matter in this case
          names(x = DataToReturn) <-
            DataToReturn %>%
            names(x = .) %>%
            lapply(
              X = .,
              FUN = function(currentMetabolite)
              {
                ValueToReturn <-
                  currentMetabolite %>%
                    gsub(
                      pattern = "(_ngperml|_pgperml)",
                      replacement = "",
                      x = .
                      )
                return(ValueToReturn)
              }
              ) %>%
            unlist(x = .)

          return(DataToReturn)
        }
          )
  
 # combine the missing value symbol summaries together from round 1 and 2 of quantitation for each metabolite
  # and tally the total number of observations with each missing value symbol
  mapply(
    FUN = function(
                   currentMissingValueSummaryRound1,
                   currentMissingValueSummaryRound2,
                   currentMetaboliteName
                   )
    {
      DataToReturn <-
        # bind the missing value summaries for the current metabolite from each quantitation round by row
        list(
          currentMissingValueSummaryRound1,
          currentMissingValueSummaryRound2
        ) %>%
        do.call(
          what = "rbind",
          args = .
          ) %>%
        # group by the Metabolite_value symbol and sum the missing symbol count
        group_by(
          .data = .,
          Metabolite_value
          ) %>%
        mutate(
          .data = .,
          Symbol_Count = sum(Symbol_Count)
          ) %>%
        # ungroup the result
        ungroup(x = .) %>%
        # obtain unique rows
        unique(x = .) %>%
        # pivot the table so the missing symbols are the columns and the symbol counts are the values
        pivot_wider(
          data = .,
          names_from = Metabolite_value,
          values_from = Symbol_Count
          ) %>%
         # compute the total number of missing values by summing across the columns
         mutate(
           .data = .,
           Total_Missing = rowSums(across(where(is.numeric)))
         ) %>%
         # create a column with the current metabolite
        mutate(
          .data = .,
          Metabolite = currentMetaboliteName
          ) %>%
         # move the metabolite column to the front of the table
        select(
          .data = .,
          Metabolite,
          everything()
          )

        # save the final missing value symbol summary to the file system labeled with the current metabolite
        if(!dir.exists(paths = "./MissingValueSymbolSummaries/"))
        {
          dir.create(path = "./MissingValueSymbolSummaries/")
        }
      
      DataToReturn %>%
        write_tsv(
          x = .,
          file = glue("./MissingValueSymbolSummaries/MissingValueSymbolSummary_{currentMetaboliteName}.txt")
          )

    },
    MissingValueSymbolSummaries$round1,
    MissingValueSymbolSummaries$round2,
    names(x = MissingValueSymbolSummaries$round1),
    SIMPLIFY = FALSE
      )
  
###### create a VitaminDplinkPhenotypeFile with duplicate quantitations removed to be used in the candidate-gene association analysis
  # all participants must be unique for the candidate gene association analysis
VitaminDplinkPhenotypeFile <-
  VitaminDplinkPhenotypeFileWithDuplicates %>%
  # join the lookup table of all samples from candidate gene sequencing with the metabolite data by 
  # visit code based on the visit codes that are present in the LookupTableOfAllSamples 
  # after replacing visit codes in the LookupTableOfAllSamples with the visit code at the latest metabolite quantitation date for duplicate encounters
  left_join(
    x = LookupTableOfAllSamples %>%
        # replace the visit codes in the LookupTableOfAllSamples with the visit code from the latest metabolite quantitation date
        left_join(
          x = .,
          y = DuplicateMetaboliteQuantitationEncountersStudyCodes %>%
              mutate(
                .data = .,
                Visit.code__laterDate_copy = Visit.code__laterDate
              ),
          by = c("VisitCode" = "Visit.code__laterDate_copy")
        ) %>%
        # deselect unneccesary columns
        select(
          .data = .,
          -Visit.code__laterDate,
          -Participant.ID,
          -Visit.code__EarliestDate
        ),
    y = .,
    by = c("VisitCode" = "VisitCode")
    ) %>%
  # deselect unnecessary columns
  select(
    .data = .,
    -VisitCode
    ) %>%
  # create a family ID and an individual ID column with the VCFfileIDs and select all phenotypes
  select(
    .data = .,
    "FID" = VCFfileID,
    "IID" = VCFfileID,
    everything()
    ) %>%
  # filter out the observation with the extreme outlying 25(OH)D3 concentration value
  # that is 10 standard deviations greater than the median value due to excessive dietary vitamin D supplementation
  filter(
    .data = .,
    primary_25OHD3_ngperml < 170
    )

######## calculate the lower limit of quantitation (LLOQ) for each metabolite

  LowerLimitsOfQuantitation <-
    VitaminDplinkPhenotypeFile %>%
    # obtain the column names from the VitaminDplinkPhenotypeFile
    names(x = .) %>%
    # isolate the columns that contain the string "D3_ngperml" or "D2_ngperml"
    grepl(
      pattern = "(D3_ngperml|D2_ngperml)",
      x = .
      ) %>%
    which(x = .) %>%
    names(x =VitaminDplinkPhenotypeFile)[.] %>%
    # loop through each column
    lapply(
      X = .,
      FUN = function(currentMetaboliteColumn)
      {
          LLOQ <-
            VitaminDplinkPhenotypeFile %>%
            # select the current metabolite
            select(
              .data = .,
              "currentMetabolite" = all_of(x = currentMetaboliteColumn)
              ) %>%
            # change values of -9 to missing (NA)
            mutate(
              .data = .,
              currentMetabolite = if_else(
                                          condition = as.numeric(x = currentMetabolite) == as.numeric(x = -9),
                                          true = as.numeric(x = NA),
                                          false = as.numeric(x = currentMetabolite)
                                            )
              ) %>%
            # pull out the column for the current metabolite
            pull(
              .data = .,
              currentMetabolite
              ) %>%
            # calculate the minimum value: the lower limit of quantitation (LLOQ)
            min(
              .,
              na.rm = TRUE
              )
          
          return(LLOQ)
      }
        ) %>%
    unlist(x = .)
  
  # name the array elements according to metabolite name
  names(x = LowerLimitsOfQuantitation) <-
    VitaminDplinkPhenotypeFile %>%
    # obtain the column names from the VitaminDplinkPhenotypeFile
    names(x = .) %>%
    # isolate the columns that contain the string "D3_ngperml" or "D2_ngperml"
    grepl(
      pattern = "(D3_ngperml|D2_ngperml)",
      x = .
    ) %>%
    which(x = .) %>%
    names(x =VitaminDplinkPhenotypeFile)[.]
  
  # save the LowerLimitsOfQuantitation to the file system
  if(!dir.exists(paths = "./LowerLimitsOfQuantitation/"))
  {
    dir.create(path = "./LowerLimitsOfQuantitation/")
  }
  
  # convert the named array of LowerLimitsOfQuantitation for each metabolite to a tibble and save to the file system
  LowerLimitsOfQuantitation %>%
    enframe(
      x = .,
      name = "Metabolite",
      value = "LowerLimitOfQuantitation"
      ) %>%
    write_tsv(
      x = .,
      file = "./LowerLimitsOfQuantitation/LowerLimitsOfQuantitation.txt"
      )
  
######### identify the number of missing values for each vitamin D metabolite and metabolite ratio 
        # for the participants that are included in the candidate gene association study based on 
        # the participants in the VitaminDplinkPhenotypeFile (N = 468 with the extreme outlying observation removed)

    # identify vitamin D metabolite and metabolite ratio columns
    # by finding column names that contain "ngperml" or "ratio
    MissingMetaboliteValueSummary <-
        VitaminDplinkPhenotypeFile %>%
          names(x = .) %>%
          grepl(
            pattern = "(ratio|ngperml)",
            x = .
            ) %>%
          which(x = .) %>%
          names(x = VitaminDplinkPhenotypeFile)[.] %>%
          VitaminDplinkPhenotypeFile[.] %>%
          # Replace any missing value (indicated with a -9) with a "missing label"
          # If the value is not missing, label as "not missing"
          mutate_all(
            .tbl = .,
            .funs = ~ case_when(
                                .x == as.numeric(x = -9) ~ "missing",
                                TRUE ~ "not missing"
                                )
            ) %>%
          # loop through each metabolite and count the number of missing observations and the non-missing observations 
          mapply(
            FUN = function(currentMetabolite,currentMetaboliteName)
            {
              MissingSummaryToReturn <-
                currentMetabolite %>%
                  # count missing and non-missing observations for the current metabolite
                  table() %>%
                  # convert to a dataframe
                  data.frame() %>%
                  # change the column names
                  select(
                    .data = .,
                    "Status" = .,
                    "Count" = Freq
                    ) %>%
                  # add a column with the name of the metabolite
                  mutate(
                    .data = .,
                    metabolite_or_ratio = currentMetaboliteName
                    )
              
              return(MissingSummaryToReturn)
            },
            .,
            names(x = .),
            SIMPLIFY = FALSE
              ) %>%
          # unname the list 
          unname(obj = .) %>%
          # bind the results by row
          do.call(
            what = "rbind",
            args = .
            ) %>%
          # filter to "not missing" values only
          filter(
            .data = .,
            Status == "not missing"
            ) %>%
          # calculate the number missing by subtracting the non-missing by the total sample size [468 with the extreme 25(OH)D3 outlier removed]
          mutate(
            .data = .,
            Missing_Value_Count = 468 - Count
            ) %>%
          # reorder and select relevant columns
          select(
            .data = .,
            metabolite_or_ratio,
            Missing_Value_Count
            )
    
    # save the missing metabolite value summary to the file system
    if(!dir.exists(paths = "./MissingMetaboliteValueSummary/"))
    {
      dir.create(path = "./MissingMetaboliteValueSummary/")
    }
    
    MissingMetaboliteValueSummary %>%
      write_tsv(
        x = .,
        file = "./MissingMetaboliteValueSummary/MissingMetaboliteValueSummary.txt"
        )
    
      
########### create histograms and box and whisker plots of each vitamin D metabolite and ratio and perform Shapiro-Wilkes tests to check 
          # for deviations from normality on untransformed and log10-transformed data
########### additionally, create plots of pairwise correlations of vitamin D3 and D2 metabolites

# identify vitamin D metabolite and metabolite ratio columns
# by finding column names that contain "ngperml" or "ratio"
VitaminDplinkPhenotypeFile %>%
  names(x = .) %>%
  grepl(
    pattern = "(ratio|ngperml)",
    x = .
    ) %>%
  which(x = .) %>%
  names(x = VitaminDplinkPhenotypeFile)[.] %>%
  # loop through each metabolite or ratio and create a box and whisker plot and un-transformed and log-transformed histograms 
  lapply(
    X = .,
    FUN = function(currentMetabolite)
      {
      
      #### plot values of the current metabolite without any transformations and perform a Shapiro-Wilke's test for normality

        # identify non-missing values to plot
        ValuesToPlot <-
          VitaminDplinkPhenotypeFile %>%
          # select the current metabolite from the VitaminDplinkPhenotypeFile
          select(
            .data = .,
            "value" = all_of(x = currentMetabolite)
          ) %>%
          # change values of -9 to an NA
          mutate(
            .data = .,
            value = if_else(
                            condition = value == as.numeric(x = -9),
                            true = as.numeric(x = NA),
                            false = as.numeric(x = value)
                              )
            ) %>%
          # remove missing values
          na.omit(object = .)
        
        # create a box and whisker plot of the current metabolite values to identify any concerning outliers
        BoxAndWhiskerPlot <-
          ValuesToPlot %>%
          ggplot(
            data = .,
            mapping = aes(
                          x = "",
                          y = value
                            )
              ) +
          geom_boxplot() +
          theme_classic() +
          scale_y_continuous(
            labels = seq(
                         0,
                         max(ValuesToPlot$value,na.rm = TRUE),
                         max(ValuesToPlot$value,na.rm = TRUE)/10
                         ) %>%
                      round(x = .,digits = 4),
            breaks = seq(
                         0,
                         max(ValuesToPlot$value,na.rm = TRUE),
                         max(ValuesToPlot$value,na.rm = TRUE)/10
                         ) %>%
                     round(x = .,digits = 4)
              ) +
          xlab(label = currentMetabolite) +
          ylab(label = "Concentration (ng/mL) or ratio") +
          theme(
            axis.title = element_text(family = "Arial",face = "bold",size = 11),
            axis.text.x = element_text(family = "Arial",face = "bold",size = 11),
            axis.text.y = element_text(family = "Arial",face = "bold",size = 11)
          )
        
        # save the final plot to the file system as a tiff
        if(!dir.exists(paths = "./MetaboliteBoxAndWhiskers/"))
        {
          dir.create(path = "./MetaboliteBoxAndWhiskers/")
        }
        
        tiff(
          filename = glue("./MetaboliteBoxAndWhiskers/{currentMetabolite}_BoxAndWhisker.tiff"),
          width = 5,
          height = 7,
          units = "in",
          compression = "none",
          res = 300
          )

        print(x = BoxAndWhiskerPlot)

        dev.off()
        
        # perform a Shapiro-Wilk's test for normality
        TestForNormalityResult <-
          ValuesToPlot$value %>%
          shapiro.test(x = .)

        # calculate the binwidth required for 30 histogram bins
        HistogramBinWidth <-
          ValuesToPlot %>%
          # calculate value max, min, and the binwidth for 30 bins
          summarise(
            .data = .,
            "MaxValue" = max(value,na.rm = TRUE),
            "MinValue" = min(value,na.rm = TRUE),
            "BinWidth" = (max(value,na.rm = TRUE) - min(value,na.rm = TRUE))/30
          ) %>%
          pull(.data = .,BinWidth)


          PlotToReturn <-
              ValuesToPlot %>%
              # create a histogram
              ggplot(
                data = .,
                mapping = aes(x = value)
              ) +
              geom_histogram(
                binwidth = HistogramBinWidth,
                fill="#69b3a2",
                color="black",
                alpha=0.9
              ) +
              annotate(
                geom = "text",
                x = (max(ValuesToPlot$value) - min(ValuesToPlot$value))/1.2,
                y = 40,
                label = paste0("Non-normality test (p = ",signif(x = TestForNormalityResult$p.value,digits = 3),")"),
                size = 4,
                family = "Arial"
                ) +
              theme_classic() +
              xlab(label = "concentration (ng/mL) or ratio") +
              ylab(label = "sample count") +
              ggtitle(label = currentMetabolite) +
              theme(
                axis.title = element_text(family = "Arial",face = "bold",size = 11),
                axis.text.x = element_text(family = "Arial",face = "bold",size = 11),
                axis.text.y = element_text(family = "Arial",face = "bold",size = 11),
                plot.title = element_text(family = "Arial",face = "bold",size = 11,hjust = 0.5)
              )

          # save the final plot to the file system as a tiff
          if(!dir.exists(paths = "./MetaboliteHistograms/"))
          {
            dir.create(path = "./MetaboliteHistograms/")
          }

          tiff(
            filename = glue("./MetaboliteHistograms/{currentMetabolite}_histogram.tiff"),
            width = 7,
            height = 5,
            units = "in",
            compression = "none",
            res = 300
            )

          print(x = PlotToReturn)

          dev.off()
          
          
          #### plot values of the current metabolite with a log10 transformation and perform a Shapiro-Wilke's test for normality
          
          # identify non-missing values to plot
          ValuesToPlot <-
            VitaminDplinkPhenotypeFile %>%
            # select the current metabolite from the VitaminDplinkPhenotypeFile
            select(
              .data = .,
              "value" = all_of(x = currentMetabolite)
            ) %>%
            # change values of -9 or 0 to an NA, values of 0 cannot be log-transformed
            mutate(
              .data = .,
              value = if_else(
                              condition = value == as.numeric(x = -9) | value == as.numeric(x = 0),
                              true = as.numeric(x = NA),
                              false = as.numeric(x = value)
                                )
              ) %>%
            # remove missing values
            na.omit(object = .) %>%
            # perform a log10 transformation of the metabolite or ratio values
            mutate(
              .data = .,
              value = log10(x = value)
              )

          # perform a Shapiro-Wilk's test for normality on the log10() transformed data
          TestForNormalityResult <-
            ValuesToPlot$value %>%
            shapiro.test(x = .)

          # calculate the binwidth required for 30 histogram bins
          HistogramBinWidth <-
            ValuesToPlot %>%
            # calculate value max, min, and the binwidth for 30 bins
            summarise(
              .data = .,
              "MaxValue" = max(value,na.rm = TRUE),
              "MinValue" = min(value,na.rm = TRUE),
              "BinWidth" = (max(value,na.rm = TRUE) - min(value,na.rm = TRUE))/30
            ) %>%
            pull(.data = .,BinWidth)


            PlotToReturn <-
                ValuesToPlot %>%
                # create a histogram
                ggplot(
                  data = .,
                  mapping = aes(x = value)
                ) +
                geom_histogram(
                  binwidth = HistogramBinWidth,
                  fill="#69b3a2",
                  color="black",
                  alpha=0.9
                ) +
                annotate(
                  geom = "text",
                  x = max(ValuesToPlot$value) - 1,
                  y = 40,
                  label = paste0("Non-normality test (p = ",signif(x = TestForNormalityResult$p.value,digits = 3),")"),
                  size = 4,
                  family = "Arial"
                  ) +
                theme_classic() +
                ggtitle(label = currentMetabolite) +
                xlab(label = "Log(concentration) [ng/mL] or Log(ratio)") +
                ylab(label = "sample count") +
                theme(
                  axis.title = element_text(family = "Arial",face = "bold",size = 11),
                  axis.text.x = element_text(family = "Arial",face = "bold",size = 11),
                  axis.text.y = element_text(family = "Arial",face = "bold",size = 11),
                  plot.title = element_text(family = "Arial",face = "bold",size = 11,hjust = 0.5)
                )

            # save the final plot to the file system as a tiff
            if(!dir.exists(paths = "./MetaboliteHistograms/"))
            {
              dir.create(path = "./MetaboliteHistograms/")
            }

            tiff(
              filename = glue("./MetaboliteHistograms/log10_{currentMetabolite}_histogram.tiff"),
              width = 7,
              height = 5,
              units = "in",
              compression = "none",
              res = 300
              )

            print(x = PlotToReturn)

            dev.off()
    }
  )

# suppress non-problematic warning messages for removing missing values when creating metabolite correlation plots
options(warn = -1)

 # create a scatterplot/correlation matrix of the relevant vitamin D3 metabolites and the vitamin D2 metabolites separately 
 # using the ggpairs() function from GGally package (see https://www.r-graph-gallery.com/199-correlation-matrix-with-ggally.html for more info)
c(
  "D3_ngperml",
  "D2_ngperml"
  ) %>%
  lapply(
    X = .,
    FUN = function(currentMetaboliteSet)
    {
      DataToPlot <-
          # select the metabolite columns that correspond to the currentMetaboliteSet
          VitaminDplinkPhenotypeFile %>%
          names(x = .) %>%
          grepl(
            pattern = currentMetaboliteSet,
            x = .
            ) %>%
          which(x = .) %>%
          names(x = VitaminDplinkPhenotypeFile)[.] %>%
          VitaminDplinkPhenotypeFile[.] %>%
          # change any values with a -9.0 to a numeric NA
          mutate(
            .data = .,
            across(
              .cols = everything(),
              .fns = function(currentColumn)
              {
                ValuesToReturn <-
                    data.frame(
                      "values" = currentColumn
                    ) %>%
                    mutate(
                      .data = .,
                      values = if_else(
                                       condition = values == as.numeric(x = -9),
                                       true = as.numeric(x = NA),
                                       false = as.numeric(x = values)
                                         )
                      ) %>%
                    pull(
                      .data = .,
                      values
                      )

                return(ValuesToReturn)
              }
                )
            )
      
        # remove the "_ngperml" suffix on the column names and the "character_" prefix
        names(x = DataToPlot) <-
          DataToPlot %>%
          names(x = .) %>%
          lapply(
            X = .,
            FUN = function(currentColumnName)
            {
              ValueToReturn <-
                  currentColumnName %>%
                  gsub(
                    pattern = "_ngperml",
                    replacement = "",
                    x = .
                    ) %>%
                  sub(
                    pattern = "(primary_|R_|beta_|active_)",
                    replacement = "",
                    x = .
                    )
              
              return(ValueToReturn)
            }
          ) %>%
          unlist(x = .)
      
      # create the plot with all metabolites
      PlotObject <-
        DataToPlot %>%
          ggpairs(
            data = .,
            xlab = "concentration (ng/mL)",
            ylab = "concentration (ng/mL)"
            ) +
          theme(
            text = element_text(family = "Arial",face = "bold",size = 11)
              ) 
      
      # save the final plot to the file system as a tiff
      if(!dir.exists(paths = "./MetaboliteCorrelations/"))
      {
        dir.create(path = "./MetaboliteCorrelations/")
      }
      
      # save the plot to the file system
      tiff(
        filename = glue("./MetaboliteCorrelations/{currentMetaboliteSet}.tiff"),
        width = 7,
        height = 7,
        units = "in",
        compression = "none",
        res = 300
        )
      print(x = PlotObject)
      dev.off()
    }
      )

# turn warning messages back on
options(warn = 0)

#### perform tests for heteroscedasticity with the D3 and D2 metabolites separately using the 
   # olsrr package (see https://cran.r-project.org/web/packages/olsrr/ for more info)

HeteroscedasticityTestResults <-
# loop through the D3 and D2 metabolites seperately
c("D3_ngperml","D2_ngperml") %>%
  lapply(
    X = .,
    FUN = function(currentMetaboliteSet)
    {
      
      MetabolitesToLoopThrough <-
          # select the metabolite columns that correspond to the currentMetaboliteSet
          VitaminDplinkPhenotypeFile %>%
          names(x = .) %>%
          grepl(
            pattern = currentMetaboliteSet,
            x = .
            ) %>%
          which(x = .) %>%
          names(x = VitaminDplinkPhenotypeFile)[.]
      
      # obtain a metabolite data set to be used in regressions
      MetaboliteData <-
          # select the metabolite columns from the VitaminDplinkPhenotypeFile
          MetabolitesToLoopThrough %>%
          VitaminDplinkPhenotypeFile[.] %>%
          # change any values with a -9.0 to a numeric NA
          mutate(
            .data = .,
            across(
              .cols = everything(),
              .fns = function(currentColumn)
              {
                ValuesToReturn <-
                    data.frame(
                      "values" = currentColumn
                    ) %>%
                    mutate(
                      .data = .,
                      values = if_else(
                                       condition = values == as.numeric(x = -9),
                                       true = as.numeric(x = NA),
                                       false = as.numeric(x = values)
                                         )
                      ) %>%
                    pull(
                      .data = .,
                      values
                      )

                return(ValuesToReturn)
              }
                )
            )
        
        # loop through each metabolite in the current metabolite dataset
        MetabolitesToLoopThrough %>%
          lapply(
            X = .,
            FUN = function(currentOutcomeMetabolite)
            {
                    # perform regressions of all combinations of metabolites and 
                    # test for heteroscedasticity with the F-test for heteroscedasticity from the olsrr package
                    # see https://cran.r-project.org/web/packages/olsrr/vignettes/heteroskedasticity.html for more info on the F-test
                          # The null hypothesis is that the variance is homogenous
                          # The alternative hypothesis is that the variance is not homogenous
                    MetabolitesToLoopThrough %>%
                      lapply(
                        X = .,
                        FUN = function(currentPredictorMetabolite)
                        {
                          
                              # if the currentOutcomeMetabolite is not the same as the currentPredictorMetabolite,
                              # perform a linear regression of the currentOutcomeMetabolite versus the currentPredictorMetabolite
                              # and perform the test for variance homogeneity
                              if(currentOutcomeMetabolite!=currentPredictorMetabolite)
                              {
                                  DataForRegression <-
                                      MetaboliteData %>%
                                      # select columns for the currentOutcomeMetabolite and currentPredictorMetabolite
                                      select(
                                        .data = .,
                                        "currentOutcomeMetabolite" = all_of(x = currentOutcomeMetabolite),
                                        "currentPredictorMetabolite" = all_of(x = currentPredictorMetabolite)
                                        ) %>%
                                        # remove missing values
                                      na.omit(object = .)
                                  
                                  # perform the regression of currentOutcomeMetabolite ~ currentPredictorMetabolite
                                  RegressionModel <-
                                      DataForRegression %>%
                                      lm(
                                        formula = currentOutcomeMetabolite ~ currentPredictorMetabolite,
                                        data = .
                                          )
                                  
                                  # perform the F-test for variance homogeneity of the currentOutcomeMetabolite with the RegressionModel
                                  VarianceHomogeneityTestResult <-
                                    ols_test_f(model = RegressionModel)
                                  
                                  # return a dataframe of one row with the currentOutcomeMetabolite, currentPredictorMetabolite, and 
                                  # the F-test p-value for variance homogeneity
                                  ResultToReturn <- 
                                    data.frame(
                                      "Outcome" = currentOutcomeMetabolite,
                                      "Predictor" = currentPredictorMetabolite,
                                      "Heteroscedastiticy_p.value" = signif(x = VarianceHomogeneityTestResult$p,digits = 3),
                                      "Regression" = paste0(currentOutcomeMetabolite," ~ ",currentPredictorMetabolite)
                                    ) %>%
                                    # add astericks to the p-values for levels of significance
                                    mutate(
                                      .data = .,
                                      Heteroscedastiticy_p.value = Heteroscedastiticy_p.value %>%
                                                                   AddAstericksToPvalues(columnVector = .)
                                      )
                                    
                              } else {
                                  
                                ResultToReturn <- NULL
                              }
                          
                            return(ResultToReturn)
                          
                        }
                          ) %>%
                      # bind the results by row
                      do.call(
                        what = "rbind",
                        args = .
                        )
            }
              ) %>%
          do.call(
            what = "rbind",
            args = .
            )
    }
  ) %>%
  do.call(
    what = "rbind",
    args = .
    ) %>%
  # filter the results to relevant outcome metabolite and predictor metabolite pairs
  filter(
    .data = .,
      (Regression == "primary_25OHD3_ngperml ~ VitaminD3_ngperml") | 
      (Regression == "R_24R25OH2D3_ngperml ~ primary_25OHD3_ngperml") |
      (Regression == "beta_4beta25OH2D3_ngperml ~ primary_25OHD3_ngperml") |
      (Regression == "active_1alpha25OH2D3_ngperml ~ primary_25OHD3_ngperml") |
      (Regression == "primary_25OHD2_ngperml ~ VitaminD2_ngperml") |
      (Regression == "active_1alpha25OH2D2_ngperml ~ primary_25OHD2_ngperml")
    )

# save the HeteroscedasticityTestResults to the file system with the metabolite correlations
HeteroscedasticityTestResults %>%
  write_tsv(
    x = .,
    file = "./MetaboliteCorrelations/HeteroscedasticityTestResults.tsv"
      )

# save a file of all continuous phenotypes to the file system
VitaminDplinkPhenotypeFile %>%
  select(
    .data = .,
    -VitDsufficiency,
    -VitDdeficiency,
    -OverallVitDstatus
  ) %>%
  write_tsv(
    x = .,
    file = "VitaminDplinkPhenotypeFile.txt"
    )

# save a file of binary phenotypes to the file system
VitaminDplinkPhenotypeFile %>%
  select(
    .data = .,
    FID,
    IID,
    VitDsufficiency,
    VitDdeficiency
    ) %>%
  write_tsv(
    x = .,
    file = "VitaminDplinkPhenotypeFileBinary.txt"
  )

DemographicCovariateData <-
# create a covariates files with age, BMI, gender, study date, and study month
ParticipantData %>% 
  # select relevant columns from the participant data
  select(
    .data = .,
    Participant.ID,
    Sex,
    "Blood_Quanta" = Blood.Quanta..1.100..
    ) %>%
  # for any blood quantas that are reported as a fraction, compute the decimal number
  mutate(
    .data = .,
    Blood_Quanta = Blood_Quanta %>%
                   lapply(
                     X = .,
                     FUN = function(currentValue)
                     {
                       print(x = currentValue)
                       
                       # if the value contains a "/."
                       # place a zero in between the "/" and the "." so the denominator is recognized as a number,
                       # then evaluate the fraction
                       if(grepl(pattern = "\\/\\.",x = currentValue)){
                         
                         ValuesToDivide <-
                         currentValue %>%
                           str_replace(
                             string = .,
                             pattern = "\\/\\.",
                             replacement = "/0."
                           ) %>%
                           str_split(
                             string = .,
                             pattern = "(\\/| \\/ )"
                           ) %>%
                           unlist(x = .) %>%
                           as.numeric(x = .) 
                         
                         ValueToReturn <-
                           # divide the numerator by the denominator
                           (ValuesToDivide[1]/ValuesToDivide[2]) %>%
                           as.character(x = .)
                         
                       }
                       # if the current value contains a "/" 
                       # split the string into an array on the "/"
                       # then divide the numerator and denominator
                       else if (grepl(pattern = "\\/",x = currentValue))
                       {
                         
                         ValuesToDivide <-
                         currentValue %>%
                           str_split(
                             string = .,
                             pattern = "(\\/| \\/ )"
                             ) %>%
                           unlist(x = .) %>%
                           as.numeric(x = .) 
                         
                         ValueToReturn <-
                           # divide the numerator by the denominator
                           (ValuesToDivide[1]/ValuesToDivide[2]) %>%
                           as.character(x = .)
                         
                         # if the value is NA, return an NA
                       } else if (is.na(x = currentValue)){
                         
                         ValueToReturn <- NA
                         
                         # otherwise return the current value
                       } else {
                         ValueToReturn <- currentValue
                       }
                       
                       print(x = ValueToReturn)
                  
                       return(ValueToReturn)
                     }
                       ) %>%
                    unlist(x = .)
    ) %>%
  # join the participant data with the relevant columns from the encounters data based on participant ID
  left_join(
    x = .,
    y = EncountersData %>%
        select(
          .data = .,
          Visit.code,
          Participant.ID,
          Study.Date,
          Age.on.Study.Date,
          Height.on.Study.Date,
          Weight.on.Study.Date
        ),
    by = "Participant.ID"
      ) %>%
  # join the resulting file with the lookup table of samples that were sequenced based on visit code to get sample sequencing IDs
  left_join(
    x = LookupTableOfAllSamples,
    y = .,
    by = c("VisitCode" = "Visit.code")
    ) %>%
# remove the duplicate sequencing sample (visit code A-313) that was removed during sequencing quality control
filter(
  .data = .,
  VisitCode != "A-313"
)

# save the demographic covariate data to the file system
write_tsv(
  x = DemographicCovariateData,
  file = "DemographicCovariateData.txt"
    )

# create a PLINK phenotype file for each individual phenotype
# that will be used for rare variant analysis with SKAT  

  # obtain the names of each of the phenotypes
  PhenotypesToLoopThrough <-
    VitaminDplinkPhenotypeFile %>%
    select(
      .data = .,
      -FID,
      -IID,
      -OverallVitDstatus
      ) %>%
    names(x = .)
         
  # select the sample IDs and the phenotype for each phenotype
  # from the current data set and save to the file system,
  # this will create a phenotype file for each individual phenotype
  PhenotypesToLoopThrough %>%
    lapply(
      X = .,
      FUN = function(currentPhenotype)
      {
        # create a directory to save the file to 
        if(!dir.exists(paths = glue("./pheno_files_for_rare_variant_analysis/")))
        {
          dir.create(path = glue("./pheno_files_for_rare_variant_analysis/"))
        }
        
        FileToReturn <-
            VitaminDplinkPhenotypeFile %>%
              select(
                .data = .,
                FID,
                IID,
                all_of(x = currentPhenotype)
                )
        
          FileToReturn %>%
            write_tsv(
              x = .,
              file = glue("./pheno_files_for_rare_variant_analysis/{currentPhenotype}_plink_pheno_for_rareVariantAnalysis.txt")
                )
                
                
          }
        )
  
  print(x = "$$$$$$$$$$$$$$$$$$$ Phenotype files for common and rare variant analysis were computed $$$$$$$$$$$$$$$$$$$$")
