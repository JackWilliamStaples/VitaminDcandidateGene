# this is a script for preparing phenotype and covariate files for linear and logistic regressions

# load necessary packages
source("../scripts/load_R_packages.R")
# load necessary functions
source("../scripts/genotype_phenotype_association/genetic_association_analysis_functions.R")

PhenotypeData <-
# load the vitamin D metabolite data
read.delim(
  file = "./VitaminDplinkPhenotypeFile.txt",
  header = TRUE
  )

PhenotypeData <-
  PhenotypeData %>%
  # change any -9 values in the metabolite columns (indicated with nanogram per mL or ratio) to a numeric NA
  mutate_at(
    .tbl = .,
    .vars = PhenotypeData %>%
            names(x = .) %>%
            grepl(
              pattern = "(ngperml|ratio)",
              x = .
              ) %>%
            which(x = .) %>%
            names(x = PhenotypeData)[.],
    .funs = function(currentMetaboliteColumn)
            {
              ColumnToReturn <-
                  currentMetaboliteColumn %>%
                  as.numeric(x = .) %>%
                  lapply(
                    X = .,
                    FUN = function(currentValue)
                    {
                      if(currentValue==as.numeric(x = -9))
                      {
                        ValueToReturn <- as.numeric(x = NA)
                      } else {
                        ValueToReturn <- currentValue
                      }
                      return(ValueToReturn)
                    }
                    ) %>%
                  unlist(x = .)

              return(ColumnToReturn)
            }
      )

# create a decision table for converting study month to a study season based on
# 4 vitamin D seasons, identified via the cosinor regression in the regression script
StudyMonthToStudySeasonDecisionTable <-
  data.frame(
    "Study_Month" = c(12,1,2,3,4,5,6,7,8,9,10,11),
    "Study_Season" = c(
                      rep_len(x = "C_Dec_Feb_trough",length.out = 3),
                      rep_len(x = "D_Mar_May_ascending",length.out = 3),
                      rep_len(x = "A_June_Aug_peak",length.out = 3),
                      rep_len(x = "B_Sept_Nov_descending",length.out = 3)
                    )
  )

# create a decision table for converting study month to the day of the year
StudyMonthToStudyDayDecisionTable <-
  data.frame(
    "Study_Month" = c(1,2,3,4,5,6,7,8,9,10,11,12),
    # starting with the first day of January as day 1 of the year
    "Cumulative_days_in_year" = c(1,31,28,31,30,31,30,31,31,30,31,30) %>% 
                                cumsum(x = .)
  )

# create a decision table for assigning the metabolite quantitation round (round 1, 2, 3, or 4) based on the visit code
MetaboliteQuantitationBatchDecisionTable <-
  data.frame(
    # Metabolite quantitation rounds: 
      # round 1: visit codes A-001 to A-210
      # round 2: visit codes A-211 to A-332
      # round 3: visit codes A-333 to A-416
      # round 4: visit codes A-417 to A-560
    "first_sample_visit_code" = c(1,211,333,417) %>% 
                                as.numeric(x = .),
    "last_sample_visit_code" = c(210,332,416,560) %>% 
                               as.numeric(x = .),
    "Quantitation_round" = c(
                             "1 (A-001 to A-210)",
                             "2 (A-211 to A-332)",
                             "3 (A-333 to A-416)",
                             "4 (A-417 to A-560)"
                             )
  )

CovariateData <-
# load the demographic covariate data
read.delim(
  file = "./DemographicCovariateData.txt",
  header = TRUE
  ) %>%
  # create a BMI column by dividing the weight by the height squared
  # weight in pounds needs to be converted to kilograms (1 kilogram is 2.20462 pounds)
  # height in inches needs to be converted to meters (1 meter is 39.3701 inches)
  mutate(
    .data = .,
    BMI = (
            (Weight.on.Study.Date*(1/2.20462))
            /
            ((Height.on.Study.Date*(1/39.3701))^2)
            )
    ) %>%
  # create a study month column by finding the numbers before the first backslash in the calendar date
  mutate(
    .data = .,
    StudyMonth = Study.Date %>%
                 lapply(
                   X = .,
                   FUN = function(currentDate)
                   {
                     DateToReturn <-
                     currentDate %>%
                       str_split(
                         string = .,
                         pattern = "/"
                         ) %>%
                       unlist(x = .) %>%
                       head(x = .,1)

                     return(DateToReturn)
                   }
                     ) %>%
                unlist(x = .)
    ) %>%
  # create a day in month column by finding the numbers after the first backslash in the calendar date
  mutate(
    .data = .,
    DayInMonth = Study.Date %>%
                 lapply(
                   X = .,
                   FUN = function(currentDate)
                   {
                     DateAsArray <-
                     currentDate %>%
                       str_split(string = .,pattern = "/") %>%
                       unlist(x = .)

                     DayToReturn <-
                       DateAsArray[2]

                     return(DayToReturn)
                   }
                     ) %>%
                unlist(x = .)
    ) %>%
  # create a day of year by adding the cumulative year day for the current month to the day of the month using the StudyMonthToStudyDayDecisionTable
  mutate(
    .data = .,
    DayOfYear = mapply(
                  FUN = function(currentMonth,currentDayInMonth)
                  {
                    CumulativeMonthDay <-
                    StudyMonthToStudyDayDecisionTable %>%
                      # filter the StudyMonthToStudyDayDecisionTable to the StudyMonth
                      filter(
                        .data = .,
                        Study_Month == currentMonth
                      ) %>%
                      # select the cumulative days of the year corresponding to the month
                      pull(
                        .data = .,
                        Cumulative_days_in_year
                      )

                    # add the cumulative days of the year corresponding to the month and the day in the month
                    DateToReturn <-
                      as.numeric(x = CumulativeMonthDay) + as.numeric(x = currentDayInMonth)

                    return(DateToReturn)
                  },
                  StudyMonth,
                  DayInMonth,
                  SIMPLIFY = FALSE
                    ) %>%
                  unlist(x = .)

    ) %>%
  # create a study season column by filtering the StudyMonthToStudySeasonDecisionTable to the current month
  mutate(
    .data = .,
    StudySeason = StudyMonth %>%
                  lapply(
                    X = .,
                    FUN = function(currentMonth)
                    {

                      SeasonToReturn <-
                      StudyMonthToStudySeasonDecisionTable %>%
                        filter(
                          .data = .,
                          Study_Month == currentMonth
                          ) %>%
                        pull(.data = .,Study_Season)

                      # if there is no season returned, return an NA
                      if(length(x = SeasonToReturn)==0)
                      {
                        SeasonToReturn <- as.character(x = NA)
                      }

                      return(SeasonToReturn)
                    }
                      ) %>%
                    unlist(x = .)

    ) %>%
  # ensure that age and BMI are numerics
  mutate_at(
    .tbl = .,
    .vars = c("Age.on.Study.Date","BMI"),
    .funs = function(currentColumn)
            {
              DataToReturn <-
              as.numeric(x = currentColumn)

              return(DataToReturn)
            }
      ) %>%
  # ensure that gender and study season are categorical
  mutate_at(
    .tbl = .,
    .vars = c("Sex","StudySeason"),
    .funs = function(currentColumn)
            {
              DataToReturn <-
              factor(x = currentColumn)

              return(DataToReturn)
            }
      )

# save the CovariateData that has been prepared to the file system to be used in other scripts
CovariateData %>%
  # save to the file system
  write_tsv(
    x = .,
    file = "./CovariateData_prepared_for_regressions.txt"
  )

PhenotypeAndCovariateData <-
# join the metabolite data with the covariate data
PhenotypeData %>%
  left_join(
    x = .,
    y = CovariateData,
    by = c("FID" = "VCFfileID")
      ) %>%
  # create a visit code column with the "A-" prefix removed
  mutate(
    .data = .,
    VisitCode_noPrefix = VisitCode %>%
                          gsub(
                            pattern = "A-",
                            replacement = "",
                            x = .
                            ) %>%
                          as.numeric(x = .)
    ) %>%
  # create a quantitation round column based on the VisitCode_noPrefix column and the MetaboliteQuantitationBatchDecisionTable
  mutate(
    .data = .,
    Quantitation_round = VisitCode_noPrefix %>%
                         lapply(
                           X = .,
                           FUN = function(currentVisitCode)
                           {
                             QuantitationRoundToReturn <-
                               # filter to the quantitation round based on the current visit code
                                MetaboliteQuantitationBatchDecisionTable %>%
                                    filter(
                                      .data = .,
                                        # if the current visit code is between the first and last sample visit codes
                                        # for the quantitation round, select the quantitation round
                                        (first_sample_visit_code <= currentVisitCode) & 
                                        (last_sample_visit_code >= currentVisitCode)
                                      ) %>%
                                  pull(
                                    .data = .,
                                    Quantitation_round
                                    )
                                 
                              return(QuantitationRoundToReturn) 
                           }
                             ) %>%
                           unlist(x = .) %>%
                           factor(x = .)
                            
    ) %>%
  # deselect the visit code with no prefix now that it is no longer needed
  select(
    .data = .,
    -VisitCode_noPrefix
    )

# save the PhenotypeAndCovariateData to the file system to be used in other scripts
PhenotypeAndCovariateData %>%
  # save to the file system
  write_tsv(
    x = .,
    file = "./PhenotypeAndCovariateData_prepared_for_regressions.txt"
    )

# create a plink format covariate file and save to the file system for rare variant analysis
PhenotypeAndCovariateData %>%
  select(
    .data = .,
    FID,
    IID,
    Age.on.Study.Date,
    BMI,
    Sex,
    StudySeason
    ) %>%
  # change missing values to -9 as required by PLINK
  mutate_at(
    .tbl = .,
    .vars = c(
              "Age.on.Study.Date",
              "BMI",
              "Sex",
              "StudySeason"
              ),
    .funs = function(currentColumn)
    {
        ColumnToReturn <-
          data.frame(
            "ColumnValues" = currentColumn
           ) %>%
          # make sure the column is a character vector
          mutate(
            .data = .,
            ColumnValues = ColumnValues %>% as.character(x = .)
            ) %>%
          # if the column value is missing, convert it to a -9
           mutate(
             .data = .,
             ColumnValues = if_else(
                                    condition = is.na(x = ColumnValues),
                                    true = as.character(x = -9),
                                    false = as.character(x = ColumnValues)
                                    )
             ) %>%
          pull(.data = .,ColumnValues)

        return(ColumnToReturn)

    }
      ) %>%
  # one-hot encode the gender columns
  mutate(
    .data = .,
    Sex_male = Sex,
    Sex_female = Sex
    ) %>%
  mutate(
    .data = .,
    Sex_male = if_else(
                       condition = (Sex_male == "Male") & (Sex_male != as.character(x = -9)),
                       true = as.character(x = 1),
                       false = if_else(
                                       condition = Sex_male == as.character(x = -9),
                                       true = Sex_male,
                                       false = as.character(x = 0)
                                       )
                         ),
    Sex_female = if_else(
                         condition = (Sex_female == "Female") & (Sex_female != as.character(x = -9)),
                         true = as.character(x = 1),
                         false = if_else(
                                         condition = Sex_female == as.character(x = -9),
                                         true = Sex_female,
                                         false = as.character(x = 0)
                                         )
                        )
    ) %>%
  # one-hot encode the 4 seasons which is required for PLINK
      # 1 or 0 for A_June_Aug_peak,
      # 1 or 0 for B_Sept_Nov_descending,
      # 1 or 0 for C_Dec_Feb_trough,
      # 1 or 0 for D_Mar_May_ascending
  mutate(
    .data = .,
    StudySeason_A_June_Aug_peak = StudySeason,
    StudySeason_B_Sept_Nov_descending = StudySeason,
    StudySeason_C_Dec_Feb_trough = StudySeason,
    StudySeason_D_Mar_May_ascending = StudySeason
  ) %>%
  mutate(
    .data = .,
    StudySeason_A_June_Aug_peak = if_else(
                                          condition = (StudySeason_A_June_Aug_peak == "A_June_Aug_peak") & (StudySeason_A_June_Aug_peak!=as.character(x = -9)),
                                          true = as.character(x = 1),
                                          false = if_else(
                                                          condition = StudySeason_A_June_Aug_peak == as.character(x = -9),
                                                          true = StudySeason_A_June_Aug_peak,
                                                          false = as.character(x = 0)
                                                          )
                                            ),
    StudySeason_B_Sept_Nov_descending = if_else(
                                                condition = (StudySeason_B_Sept_Nov_descending == "B_Sept_Nov_descending") & (StudySeason_B_Sept_Nov_descending!=as.character(x = -9)),
                                                true = as.character(x = 1),
                                                false = if_else(
                                                                condition = StudySeason_B_Sept_Nov_descending == as.character(x = -9),
                                                                true = StudySeason_B_Sept_Nov_descending,
                                                                false = as.character(x = 0)
                                                                )
                                              ),
    StudySeason_C_Dec_Feb_trough = if_else(
                                          condition = (StudySeason_C_Dec_Feb_trough == "C_Dec_Feb_trough") & (StudySeason_C_Dec_Feb_trough!=as.character(x = -9)),
                                          true = as.character(x = 1),
                                          false = if_else(
                                                          condition = StudySeason_C_Dec_Feb_trough == as.character(x = -9),
                                                          true = StudySeason_C_Dec_Feb_trough,
                                                          false = as.character(x = 0)
                                                          )
                                        ),
    StudySeason_D_Mar_May_ascending = if_else(
                                              condition = (StudySeason_D_Mar_May_ascending == "D_Mar_May_ascending") & (StudySeason_D_Mar_May_ascending!=as.character(x = -9)),
                                              true = as.character(x = 1),
                                              false = if_else(
                                                              condition = StudySeason_D_Mar_May_ascending == as.character(x = -9),
                                                              true = StudySeason_D_Mar_May_ascending,
                                                              false = as.character(x = 0)
                                                              )
                                            )
    ) %>%
  select(
    .data = .,
    -Sex,
    -StudySeason
    ) %>%
  # save to the file system
  write_tsv(
    x = .,
    file = "./CovariateFile_PLINKformat.txt"
    )

# create an array of haplotype sample IDs from the vitamin D binding protein haplotype sample file
VitaminDbindingProteinHaplotypeSampleIDs <-
  # load the vitamin D binding protein sample file
  "./GC_haplotypes.sample" %>%
  read.table(file = .,header = TRUE) %>%
  # remove the first row that contains a sample ID of "0"
  filter(
    .data = .,
    ID_1 != 0
    ) %>%
  # create columns for haplotype sample IDs
  mutate(
    .data = .,
    haplotypeID1 = paste0("X","_",ID_1,"_","hap1"),
    haplotypeID2 = paste0("X","_",ID_2,"_","hap2")
    ) %>%
  # select the haplotype sample IDs
  select(
    .data = .,
    haplotypeID1,
    haplotypeID2
    ) %>%
  # create an unnamed array of all the haplotype IDs
  unlist(x = .) %>%
  unname(obj = .) %>%
  # sort the array alpha-numerically
  sort(x = .)

# load the vitamin D binding protein haplotype matrix
VitaminDbindingProteinHaplotypeMatrix <-
  "./GC_haplotypes.hap" %>%
  read.table(file = .,header = FALSE)

# update the column names of the vitamin D binding protein haplotype matrix using the haplotype sample IDs
names(x = VitaminDbindingProteinHaplotypeMatrix) <-
  c("chrom","snp","pos","ref","alt") %>%
  c(.,VitaminDbindingProteinHaplotypeSampleIDs)

# filter the haplotype matrix to the common vitamin D binding protein polymorphisms only
VitaminDbindingProteinHaplotypeMatrix <-
  VitaminDbindingProteinHaplotypeMatrix %>%
  # filter the haplotype matrix to the common vitamin D binding protein polymorphisms:
    # rs7041 (4:72618334:A:C) and rs4588 (4:72618323:G:T) variants
  filter(
    .data = .,
    snp == "4:72618334_A_C" | snp == "4:72618323_G_T"
  ) %>%
  # remove unnecessary columns
  select(
    .data = .,
    -chrom,
    -pos,
    -ref,
    -alt
    ) %>%
  # append "snp" to the snp ID so the column can be used to update the rownames
  mutate(
    .data = .,
    snp = paste0("snp_",snp)
    )

# convert the first column to the row names and remove the first column
rownames(x = VitaminDbindingProteinHaplotypeMatrix) <- VitaminDbindingProteinHaplotypeMatrix$snp
VitaminDbindingProteinHaplotypeMatrix$snp <- NULL

VitaminDbindingProteinHaplotypeMatrix <-
  # transpose the rows and columns of the VitaminDbindingProteinHaplotypeMatrix to get a sample x haplotype matrix
  VitaminDbindingProteinHaplotypeMatrix %>%
  as.matrix(x = .) %>%
  t(x = .) %>%
  as.data.frame(x = .) %>%
  # create a column that contains the row names
  mutate(
    .data = .,
    Sample_haplotype_ID = rownames(x = .)
    ) %>%
  # create a column that contains the sample ID also with out the hap1, hap2, or X characters
  mutate(
    .data = .,
    Sample_ID = Sample_haplotype_ID %>%
                gsub(pattern = "X_",replacement = "",x = .) %>%
                gsub(pattern = "(_hap1|_hap2)",replacement = "",x = .)
    ) %>%
  # arrange the snp columns to keep them in the order rs7041 (4:72618334:A:C) followed by rs4588 (4:72618323:G:T)
  select(
    .data = .,
    Sample_ID,
    Sample_haplotype_ID,
    "rs7041_4:72618334_A_C" = `snp_4:72618334_A_C`,
    "rs4588_4:72618323_G_T" = `snp_4:72618323_G_T`
    ) %>%
  # create a column that contains the haplotype, based on the 3 possible combinations of the rs7041 (4:72618334:A:C) and rs4588 (4:72618323:G:T) variants
  mutate(
    .data = .,
    GC_haplotype = if_else(
                          # the 1F allele (rs7041-T/rs4588-C):
                          # this corresponds to the reference allele for rs7041 == 0 and the reference allele for rs4588 == 0
                           condition = (`rs7041_4:72618334_A_C` == 0) & (`rs4588_4:72618323_G_T` == 0),
                           true = "GC1F",
                           false = if_else(
                                 # the 1S allele (rs7041-G/rs4588-C):
                                 # this corresponds to the ALTERNATE allele for rs7041 == 1 and the reference allele for rs4588 == 0
                                          condition = (`rs7041_4:72618334_A_C` == 1) & (`rs4588_4:72618323_G_T` == 0),
                                          true = "GC1S",
                                          false = if_else(
                                            # the 2 allele (rs7041-T/rs4588-A):
                                            # this corresponds to the reference allele for rs7041 == 0 and ALTERNATE allele for rs4588 == 1
                                                          condition = (`rs7041_4:72618334_A_C` == 0) & (`rs4588_4:72618323_G_T` == 1),
                                                          true = "GC2",
                                                          false = "GC_both_alternate"
                                                            )
                                            )
                             )
    ) %>%
  # group by the sample ID
  group_by(
    .data = .,
    Sample_ID
    ) %>%
  # split into a list of dataframes based on sample ID
  group_split(.tbl = .) %>%
  # create columns sample ID, GC_hap1, GC_hap2, and GC_diplotype for each individual sample ID
  lapply(
    X = .,
    FUN = function(currentDataFrame)
    {
      # return a tibble with sample ID, haplotypes, and diplotype
      TibbleToReturn <-
      tibble(
        # use the sample ID as the sample ID
        "Sample_ID" = currentDataFrame$Sample_ID %>% unique(x = .),
        # filter to the sample_haplotype_ID for haplotype 1 and select the haplotype
        "GC_hap1" = currentDataFrame %>%
                    filter(
                      .data = .,
                      Sample_haplotype_ID == paste0("X_",Sample_ID,"_hap1")
                      ) %>%
                    pull(.data = .,GC_haplotype),
        # filter to the sample_haplotype_ID for haplotype 2 and select the haplotype
        "GC_hap2" = currentDataFrame %>%
                    filter(
                      .data = .,
                      Sample_haplotype_ID == paste0("X_",Sample_ID,"_hap2")
                      ) %>%
                    pull(.data = .,GC_haplotype)
      ) %>%
      # combine the GC_hap1 and GC_hap2 columns to obtain the GC diplotype
        # sort the GC_hap1 and GC_hap2 before combining into a diplotype so the genotype permutation is the same for equivalent heterozygotes
      mutate(
        .data = .,
        GC_dip = c(GC_hap1,GC_hap2) %>%
                 sort(x = .) %>%
                 paste0(.,collapse = "/")
        )

      return(TibbleToReturn)

    }
      ) %>%
  # bind the results together by row
  do.call(
    what = "rbind",
    args = .
    )

### save the VitaminDbindingProteinHaplotypeMatrix to the file system to be used in other scripts
VitaminDbindingProteinHaplotypeMatrix %>%
  write_tsv(
    x = .,
    file = "./VitaminDbindingProteinHaplotypeMatrix_preparedForRegressions.txt"
      )

# create an array of metabolites to use as dependent variables in regression
Metabolites <-
PhenotypeAndCovariateData %>%
  names(x = .) %>%
  grepl(pattern = "(ngperml|ratio)",x = .) %>%
  which(x = .) %>%
  names(x = PhenotypeAndCovariateData)[.]

# create the final, clean dataset for performing association tests of metabolite concentration versus diplotype
VitaminDbindingProteinDataForLinearRegression <-
# join the VitaminDbindingProteinHaplotype data with the metabolite phenotype data for each sample to perform regressions against GC diplotype carriers
VitaminDbindingProteinHaplotypeMatrix %>%
  # remove the duplicate copy of the sample ID from the sample ID column
  mutate(
    .data = .,
    Sample_ID = Sample_ID %>%
                str_split(
                  string = .,
                  pattern = "_"
                  ) %>%
                unlist(x = .) %>%
                unique(x = .)
    ) %>%
  # join the VitaminDbindingProteinHaplotypeMatrix with the PhenotypeAndCovariateData by sample sequencing ID
  left_join(
    x = .,
    y = PhenotypeAndCovariateData %>%
        # select the sample ID and metabolite columns from the phenotype and covariate data
        select(
          .data = .,
          FID,
          all_of(x = Metabolites),
          Age.on.Study.Date,
          BMI,
          Sex,
          StudySeason
        ) %>%
        # make sure the sample ID is a character vector
       mutate(
         .data = .,
         FID = FID %>% as.character(x = .)
       ),
     by = c("Sample_ID" = "FID")
    )

### save the VitaminDbindingProteinDataForLinearRegression to the file system to be used in other scripts
VitaminDbindingProteinDataForLinearRegression %>%
  write_tsv(
    x = .,
    file = "./VitaminDbindingProteinHaplotypeMatrix_withPhenotypeData_prepared_for_regressions.txt"
  )

ClumpData <-
# load all of the linear additive association CLUMP files for reducing the independent metabolite common variant associations
  # down to only an independent set of causative SNPs
  (grepl(
    # identify files that do not contain the following patterns
    pattern = "(log|nosex|logistic|qassoc)",
    x = list.files(path = "./regression_results_temp/")
    )==FALSE) %>%
  which(x = .) %>%
  list.files(path = "./regression_results_temp/")[.] %>%
  mclapply(
    X = .,
    FUN = function(currentClumpFileToLoad)
    {
      LoadedFile <-
        # load the current clump file
        glue("./regression_results_temp/{currentClumpFileToLoad}") %>%
        read.table(
          file = .,
          header = TRUE
          ) %>%
        # add a column with the current gene based on the clump file name
        mutate(
          .data = .,
          Gene = currentClumpFileToLoad %>%
                 str_extract(
                   string = .,
                   pattern = "(CYP3A4|DHCR7|CYP2R1|CYP27B1|SULT2A1|UGT1A4|GC|VDR|CASR|CUBN|CYP24A1|RXRG|RXRB|RXRA)"
                     )
          ) %>%
        # add a column with the current metabolite based on the clump file name
        mutate(
          .data = .,
          metabolite = currentClumpFileToLoad %>%
                       str_extract(
                         string = .,
                         pattern = "(active_1alpha25OH2D2_ngperml|active_1alpha25OH2D3_ngperml|beta_4beta25OH2D3_ngperml|Overall_primary_25OHD_ngperml|primary_25OHD2_ngperml|primary_25OHD3_ngperml|R_24R25OH2D3_ngperml|VitaminD2_ngperml|VitaminD3_ngperml|primary_25OHD3_VitaminD3_ratio|active_1alpha25OH2D3_primary_25OHD3_ratio|R_24R25OH2D3_primary_25OHD3_ratio|beta_4beta25OH2D3_primary_25OHD3_ratio|primary_25OHD2_VitaminD2_ratio|active_1alpha25OH2D2_primary_25OHD2_ratio)"
                           )
          )

        return(LoadedFile)
    },mc.cores = detectCores()
      ) %>%
  do.call(
    what = "rbind",
    args = .
    )

#### save the ClumpData to the file system to be used in other scripts
ClumpData %>%
  write_tsv(
    x = .,
    file = "./ClumpData_prepared_for_linear_regressions.txt" 
      )


# load the binary vitamin D status phenotype data
BinaryPhenotypeData <-
  read.delim(
    file = "./VitaminDplinkPhenotypeFileBinary.txt",
    header = TRUE
    ) %>%
  # change -9 values in the vitamin D status columns to a numeric NA
  mutate_at(
    .tbl = .,
    .vars = c("VitDsufficiency","VitDdeficiency"),
    .funs = function(currentPhenotypeColumn)
    {
      ColumnToReturn <-
        currentPhenotypeColumn %>%
        as.numeric(x = .) %>%
        lapply(
          X = .,
          FUN = function(currentValue)
          {
            if(currentValue==as.numeric(x = -9))
            {
              ValueToReturn <- as.numeric(x = NA)
            } else {
              ValueToReturn <- currentValue
            }
            return(ValueToReturn)
          }
          ) %>%
        unlist(x = .)

        return(ColumnToReturn)
    }
      )

#### save the BinaryPhenotypeData to be used in other scripts
BinaryPhenotypeData %>%
  write_tsv(
    x = .,
    file = "./BinaryPhenotypeData_prepared_for_logistic_regressions.txt"
  )

# join the binary phenotype data with the covariate data
BinaryPhenotypeAndCovariateData <-
  BinaryPhenotypeData %>%
  left_join(
    x = .,
    y = CovariateData,
    by = c("FID" = "VCFfileID")
    )

### save the BinaryPhenotypeAndCovariateData to the file system to be used in other scripts
BinaryPhenotypeAndCovariateData %>%
  write_tsv(
    x = .,
    file = "./BinaryPhenotypeAndCovariateData_prepared_for_logistic_regressions.txt"
      )

# load the vitamin D binding protein diplotype data that has been prepared for regressions
VitaminDbindingProteinDataForLogisticRegression <-
  VitaminDbindingProteinHaplotypeMatrix %>%
  # remove the duplicate copy of the sample ID from the sample ID column
  mutate(
    .data = .,
    Sample_ID = Sample_ID %>%
                str_split(
                  string = .,
                  pattern = "_"
                  ) %>%
                unlist(x = .) %>%
                unique(x = .)
    ) %>%
  # join the VitaminDbindingProteinHaplotypeMatrix_preparedForRegressions with the phenotype and covariate data
  left_join(
    x = .,
    y = BinaryPhenotypeAndCovariateData %>%
        # select the sample ID, phenotype, and covariate columns
        select(
          .data = .,
          FID,
          VitDsufficiency,
          VitDdeficiency,
          Age.on.Study.Date,
          BMI,
          Sex,
          StudySeason
         ) %>%
         # make sure the sample ID is a character vector
        mutate(
          .data = .,
          FID = FID %>% as.character(x = .)
        ),
    by = c("Sample_ID" = "FID")
    )

### save the VitaminDbindingProteinDataForLogisticRegression to the file system for use in other scripts
VitaminDbindingProteinDataForLogisticRegression %>%
  write_tsv(
    x = .,
    file = "./VitaminDbindingProteinDataForLogisticRegression.txt"
      )


##### load the causitive independent variants from logistic regression that were identified via P-value based linkage disequilibrium pruning
  ClumpDataLogistic <-
    (grepl(
      # identify files that do not contain the following patterns
      pattern = "(log$|nosex$|means$|adjusted$|qassoc$|additive|logistic$)",
      x = list.files(path = "./regression_results_temp/")
      )==FALSE) %>%
    which(x = .) %>%
    # load all of the logistic association clump files
    list.files(path = "./regression_results_temp/")[.] %>%
    lapply(
      X = .,
      FUN = function(currentClumpFileToLoad)
      {
        LoadedFile <-
          # load the current clump file
          glue("./regression_results_temp/{currentClumpFileToLoad}") %>%
          read.table(
            file = .,
            header = TRUE
            ) %>%
          # add a column with the current gene based on the clump file name
          mutate(
            .data = .,
            Gene = currentClumpFileToLoad %>%
                   str_extract(
                     string = .,
                     pattern = "(CYP3A4|DHCR7|CYP2R1|CYP27B1|SULT2A1|UGT1A4|GC|VDR|CASR|CUBN|CYP24A1|RXRG|RXRB|RXRA)"
                       )
            ) %>%
          # add a column with the current phenotype based on the clump file name
          mutate(
            .data = .,
            phenotype = currentClumpFileToLoad %>%
                         str_extract(
                           string = .,
                           pattern = "(VitDsufficiency|VitDdeficiency)"
                             )
            )

          return(LoadedFile)
      }
        ) %>%
    # bind all of the loaded files by row
    do.call(
      what = "rbind",
      args = .
      )
  
  #### save the ClumpDataLogistic to the file system to be used in other scripts
  ClumpDataLogistic %>%
    write_tsv(
      x = .,
      file = "./ClumpDataLogistic_prepared_for_logistic_regressions.txt"
        )
  
  

