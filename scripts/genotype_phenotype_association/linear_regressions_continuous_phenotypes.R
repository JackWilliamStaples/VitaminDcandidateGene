# This script performs analysis of continuous metabolite concentrations and ratios, the analysis includes:
  # (1) D3 and D2 metabolite phenotype principal components analysis to visualize metabolite quantitation batch effects
  # (2) Cosinor seasonal regression analysis to identify seasonal 25(OH)D3 trends with and without demographic covariates
  # (3) Grand summary statistics and summary statistics stratified by season and gender for demographic factors
  # (4) Summary statistics for clinical vitamin D status
  # (5) Univariate analysis of metabolite concentration and ratio versus individual seasonal and demographic factors
  # (6) Association analysis for significant relationships among seasonal and demographic factors (age, BMI, gender, and season)
  # (7) Tests for confounding relationships between demographic and seasonal factors
  # (8) Tests for confounding of genetic associations by seasonal and demographic factors
  # (9) multiple regressions of significant variants identified in primary analysis with adjustment for all demographic/seasonal covariates
  # (10) multiple regression with all significant causative variants identified in primary analysis to estimate polygenic contribution 
         # to phenotypic variability for each phenotype
  # (11) multiple regression with all significant causative variants identified in primary analysis and seasonal/demographic 
         # covariates together to estimate genetic and environmental contribution to variability of each phenotype

# load necessary packages
source("../scripts/load_R_packages.R")
# load necessary functions
source("../scripts/genotype_phenotype_association/genetic_association_analysis_functions.R")

# load the covariate data prepared for regressions
CovariateData <-
  read.delim(
    file = "./CovariateData_prepared_for_regressions.txt",
    header = TRUE
      )

# load the combined phenotype and covariate data prepared for regressions
PhenotypeAndCovariateData <-
  read.delim(
    file = "./PhenotypeAndCovariateData_prepared_for_regressions.txt",
    header = TRUE
      )

### binarize demographic covariates for performing cosinor regressions,
   # the cosinor package only accepts binary covariates
   # The covariates are binarized as follows:
      # (1) Obesity: BMI greater than or equal to 30 kg/sq. meter = 1, otherwise 0
      # (2) Age greater than or equal to 40 = 1, otherwise 0 based on the same age cut-point in Fohner vitamin D study
      # (3) Gender -- Female = 1, Male = 0
PhenotypeAndCovariateDataForCosinorRegression <-
  PhenotypeAndCovariateData %>%
  mutate(
    .data = .,
    BMI = if_else(
                  # if the BMI is greater than or equal to 30 kg/sq. meter and is NOT missing, code as 1
                  condition = BMI >= 30 & !is.na(x = BMI),
                  true = as.numeric(x = 1),
                  false = if_else(
                                  # if the BMI is less than 30 kg/sq. meter and is NOT missing, code as 0
                                  condition = BMI < 30 & !is.na(x = BMI),
                                  true = as.numeric(x = 0),
                                  # code NA values as NA
                                  false = as.numeric(x = NA)
                                    )
                    ),
    Age.on.Study.Date = if_else(
                                # if the age is greater than or equal to 40 and is NOT missing, code as 1
                                condition = Age.on.Study.Date >= 40 & !is.na(x = Age.on.Study.Date),
                                true = as.numeric(x = 1),
                                false = if_else(
                                                # if the age is less than or equal to 40 and is NOT missing, code as 0
                                                condition = Age.on.Study.Date < 40 & !is.na(x = Age.on.Study.Date),
                                                true = as.numeric(x = 0),
                                                # code NA values as NA
                                                false = as.numeric(x = NA)
                                                  )
                                  ),
    Sex = if_else(
                  # if the gender is female and not missing, code as 0
                  condition = Sex == "Female" & !is.na(x = Sex),
                  true = as.numeric(x = 0),
                  # if the gender is male and not missing, code as 1
                  false = if_else(
                                  condition = Sex == "Male" & !is.na(x = Sex),
                                  true = as.numeric(x = 1),
                                  # code NA values as NA
                                  false = as.numeric(x = NA)
                                    )
                    )
    )

# load the vitamin D binding protein haplotype matrix that has been prepared for regressions
VitaminDbindingProteinHaplotypeMatrix <-
  read.delim(
    file = "./VitaminDbindingProteinHaplotypeMatrix_preparedForRegressions.txt",
    header = TRUE
    )

# load the vitamin D binding protein haplotype matrix with phenotype data added
VitaminDbindingProteinDataForLinearRegression <-
  read.delim(
    file = "./VitaminDbindingProteinHaplotypeMatrix_withPhenotypeData_prepared_for_regressions.txt",
    header = TRUE
    )

# load the dbSNP rsID lookup table
rsIDlookupTable <-
  read.delim(
    file = "./PLINKvariantID_dbSNP_rsID_lookupTable_commonVariants.txt",
    header = TRUE
      )

# load the clump data prepared for linear regressions
ClumpData <-
  read.delim(
    file = "./ClumpData_prepared_for_linear_regressions.txt",
    header = TRUE
    )

################ perform a principal components analysis of the vitamin D metabolite quantitation data
               # Do this for the D3 and D2 metabolites separately grouped by metabolite quantitation round to see if
               # there is any batch effect from the 4 different metabolite quantitation rounds

# loop through the vitamin D3 and vitamin D2 metabolites separately
c(
  "D3_ngperml",
  "D2_ngperml"
  ) %>%
  lapply(
    X = .,
    FUN = function(currentMetaboliteSet)
    {
      
      # obtain the data for the principal components analysis
      DataForPCAplot <-
            PhenotypeAndCovariateData %>%
            # select columns that match the currentMetaboliteSet and the MetaboliteQuantitationRound, and the study season
           select(
             .data = .,
             contains(match = currentMetaboliteSet),
             Quantitation_round,
             StudySeason
             ) %>%
           # change all values of -9 in metabolite columns to an NA in every row/column
          mutate(
            .data = .,
            across(
              .cols = contains(match = currentMetaboliteSet),
              .fns = ~ if_else(
                               condition = .x == as.numeric(x = -9),
                               true = as.numeric(x = NA),
                               false = as.numeric(x = .x)
                              )
                )
            )
      
      # if the current metabolite set is the vitamin D2 metabolite set,
      # remove the 1,25(OH)2D2 metabolite since this metabolite was not quantitated in all quantitation rounds
      if(currentMetaboliteSet == "D2_ngperml")
      {
        DataForPCAplot <-
          DataForPCAplot %>%
          select(
            .data = .,
            -active_1alpha25OH2D2_ngperml
            )
      }

      # create a matrix with missing values removed for the principal components analysis
      MetaboliteValueMatrixForPCA <-
        DataForPCAplot %>%
          na.omit(object = .) %>%
          # deselect the quantitation round and study season
          select(
            .data = .,
            -Quantitation_round,
            -StudySeason
            ) %>%
          as.matrix(x = .)

      # create an array to label the points on the PCA plot with the Quantitation_round
      MetaboliteQuantitationRoundLabels <-
        DataForPCAplot %>%
        na.omit(object = .) %>%
        pull(
          .data = .,
          Quantitation_round
        ) %>%
        as.factor(x = .)
      
      # create an array to label the points on the PCA plot with the StudySeason
      SeasonOfSampleLabels <-
        DataForPCAplot %>%
        na.omit(object = .) %>%
        pull(
          .data = .,
          StudySeason
        ) %>%
        as.factor(x = .)

      # perform the PCA with values centered and scaled
      PCAresults <-
        prcomp(
          x = MetaboliteValueMatrixForPCA,
          center = TRUE,
          scale. = TRUE
          )

      # plot the pca results colored by metabolite quantitation round
      PCAplot <-
        fviz_pca_ind(
          X = PCAresults,
          # color by group
          col.ind = MetaboliteQuantitationRoundLabels,
          # change the legend title
          legend.title = "Quantitation Round",
          # remove individual point labels
          geom = "point"
            ) +
        labs(
          title = currentMetaboliteSet %>%
                  gsub(
                    pattern = "_ngperml",
                    replacement = "",
                    x = .
                    ) %>%
                  paste0(.," metabolite PCA by batch"),
          x = "PC1",
          y = "PC2",
            ) +
        theme_classic() +
        theme(
          axis.title = element_text(family = "Arial",face = "bold",colour = "black",size = 11),
          title = element_text(family = "Arial",face = "bold",colour = "black",size = 11,hjust = 0.5),
          axis.text.x = element_text(family = "Arial",face = "bold",colour = "black",size = 11),
          axis.text.y = element_text(family = "Arial",face = "bold",colour = "black",size = 11),
          legend.text = element_text(family = "Arial",face = "bold",colour = "black",size = 11)
        )
      
      # create a directory for saving the PCA plot if it does not exist
      if(!dir.exists(paths = "./MetabolitePCAanalysis/"))
      {
        dir.create(path = "./MetabolitePCAanalysis/")
      }

      # save the PCA results as a .tiff image
      tiff(
        filename = glue("./MetabolitePCAanalysis/PCA_results_{currentMetaboliteSet}.tiff"),
        width = 7,
        height = 5,
        units = "in",
        compression = "none",
        res = 300
        )

      print(x = PCAplot)

      dev.off()
      
      # plot the pca results colored by season of sampling
      PCAplot <-
        fviz_pca_ind(
          X = PCAresults,
          # color by group
          col.ind = SeasonOfSampleLabels,
          # change the legend title
          legend.title = "Season of sample",
          # remove individual point labels
          geom = "point"
            ) +
        #theme_minimal() +
        labs(
          title = currentMetaboliteSet %>%
                  gsub(
                    pattern = "_ngperml",
                    replacement = "",
                    x = .
                    ) %>%
                  paste0(.," metabolite PCA by season of sample"),
          x = "PC1",
          y = "PC2",
            ) +
        theme_classic() +
        theme(
          axis.title = element_text(family = "Arial",face = "bold",colour = "black",size = 11),
          title = element_text(family = "Arial",face = "bold",colour = "black",size = 11,hjust = 0.5),
          axis.text.x = element_text(family = "Arial",face = "bold",colour = "black",size = 11),
          axis.text.y = element_text(family = "Arial",face = "bold",colour = "black",size = 11),
          legend.text = element_text(family = "Arial",face = "bold",colour = "black",size = 11)
        )
      
      # save the PCA results as a .tiff image
      tiff(
        filename = glue("./MetabolitePCAanalysis/PCA_results_{currentMetaboliteSet}_by_season.tiff"),
        width = 7,
        height = 5,
        units = "in",
        compression = "none",
        res = 300
        )
      
      print(x = PCAplot)

      dev.off()
      
    }
      )

############# cosinor regression analysis of vitamin D and vitamin D metabolite levels 
            # versus time of year to identify seasonal trends


# temporarily suppress non-problematic warning messages for the cosinor analysis boxplot
# there is always a warning message that the seasonal boxplot could not calculate a value for January since there is no data
options(warn = -1) 

# create cosinor regression plots of 25(OH)D3 and 25(OH)D [25(OH)D3 + 25(OH)D2] 
# concentration with data grouped by calendar month
PhenotypeAndCovariateData %>%
  # isolate the column names of the PhenotypeAndCovariateData
  names(x = .) %>%
  # obtain column names that only contain "primary_25OHD3_ngperml" or "Overall_primary_25OHD_ngperml"
  .[grepl(pattern = "(primary_25OHD3_ngperml|Overall_primary_25OHD_ngperml)",x = .)] %>%
  # loop through each of the selected metabolites and create a plot of the cosinor regression
  # line overlayed with the actual data points with a boxplot overlay
  lapply(
    X = .,
    FUN = function(currentMetabolite)
    {
      MetabolitePlot <-
        currentMetabolite %>%
        CreateCosinorModelWithNoCovariatesGroupedByMonth(
          metaboliteOfInterest = .,
          MetaboliteDataToUse = PhenotypeAndCovariateData,
          periodForOneWaveCycle = 12,
          FitToMedian = FALSE
        )

      # create a directory for saving the cosinor plots
      if(!dir.exists(paths = "./seasonal_metabolite_plots/"))
      {
        dir.create(path = "./seasonal_metabolite_plots/")
      }
      
      # save plot as a .tiff
      tiff(
        filename = glue("./seasonal_metabolite_plots/{currentMetabolite}_cosinor_fitToAllDataByMonth_noCovar.tiff"),
        width = 9,
        height = 5,
        compression = "none",
        units = "in",
        res = 300
        )
      
      print(MetabolitePlot)
      dev.off()

    }
  )

# create cosinor regression plots of 25(OH)D3 and 25(OH)D [25(OH)D3 + 25(OH)D2] concentration
# with the cosinor regression fitted to ALL of the data by the exact day of a 365 day year
PhenotypeAndCovariateData %>%
  # isolate the column names of the PhenotypeAndCovariateData
  names(x = .) %>%
  # obtain column names that only contain "primary_25OHD3_ngperml" or "Overall_primary_25OHD_ngperml"
  .[grepl(pattern = "(primary_25OHD3_ngperml|Overall_primary_25OHD_ngperml)",x = .)] %>%
  # loop through each of the selected metabolites and create a plot of the cosinor regression
  # line overlayed with the actual data points overlayed
  lapply(
    X = .,
    FUN = function(currentMetabolite)
    {
      MetabolitePlot <-
        currentMetabolite %>%
        CreateCosinorModelWithNoCovariatesByDayOfYear(
          metaboliteOfInterest = .,
          MetaboliteDataToUse = PhenotypeAndCovariateData
        )

      # create a directory for saving the cosinor plots
      if(!dir.exists(paths = "./seasonal_metabolite_plots/"))
      {
        dir.create(path = "./seasonal_metabolite_plots/")
      }
      
      # save the plot as a .tiff
      tiff(
        filename = glue("./seasonal_metabolite_plots/{currentMetabolite}_cosinor_fitToAllDataByDayOfYear_noCovar.tiff"),
        width = 9,
        height = 5,
        units = "in",
        compression = "none",
        res = 300
        )
      
      print(MetabolitePlot)
      dev.off()

    }
  )


# switch the warning messages back on
options(warn = 0)

# create a directory for saving univariate cosinor regressions if it does not exist
if(!dir.exists(paths = "./univariate_cosinor_regressions/"))
{
  dir.create(path = "./univariate_cosinor_regressions/")
}

# calculate sample sizes for the binarized covariates included in the cosinor regression
  # loop through an array of demographic covariates
  mapply(
    FUN = function(currentCovariate,currentCovariateLabel)
    {
          # select the column for the current covariate with missing observations removed and convert to an array
          BinaryCovariateArray <-
              PhenotypeAndCovariateDataForCosinorRegression %>%
                select(
                  .data = .,
                  all_of(x = currentCovariate)
                   ) %>%
                na.omit(object = .) %>%
                pull(.data = .)
          
          # compute the total number of non-missing observations for the current covariate
          TotalN <-
            BinaryCovariateArray %>%
            length(x = .)
          
          # compute the number of entries with the binary covariate coded as 1, by taking the sum of the array
          CovariateEqualToOneN <-
            BinaryCovariateArray %>%
            sum(.)
          
          # compute the number of entries with the binary covariate coded as 0, by subtracting the total length by the number 
          # of observations coded as 1
          CovariateEqualToZeroN <-
            TotalN - CovariateEqualToOneN
          
          # return the results as a dataframe with a label for describing the covariate
          DataToReturn <-
            data.frame(
              "covariate" = currentCovariateLabel,
              "total_N" = TotalN,
              "covariate_is_1_N" = CovariateEqualToOneN,
              "covariate_is_0_N" = CovariateEqualToZeroN
            )
      
      return(DataToReturn)
    },
    c("Age.on.Study.Date","BMI","Sex"),
    c("Age >= 40 == 1","BMI >= 30 == 1","Gender = Male == 1"),
    SIMPLIFY = FALSE
      ) %>%
    do.call(
      what = "rbind",
      args = .
      ) %>%
    write_tsv(
      x = .,
      file = "./univariate_cosinor_regressions/BinaryCovariateSampleSizes.txt"
    )
  
# evaluate each demographic covariate with the cosinor regression individually for 25(OH)D3 and 25(OH)D [D3 + D2] metabolites
UnivariateCosinorRegressionResults <-
  # loop through the primary 25(OH)D3 and 25(OH)D metabolites
  c(
    "primary_25OHD3_ngperml",
    "Overall_primary_25OHD_ngperml"
    ) %>%
  lapply(
    X = .,
    FUN = function(currentMetabolite)
    {
      DataToReturn <-
      # loop through an array of demographic covariates
      c(
        "Age.on.Study.Date",
        "BMI",
        "Sex"
        ) %>%
        # perform a cosinor regression of metabolite concentration versus each covariate individually
        # with time as the day of the 365 day year
        lapply(
          X = .,
          FUN = function(currentCovariate)
          {
            # create the regression formula
            RegressionFormula <-
              glue("{currentMetabolite} ~ time(DayOfYear) + {currentCovariate} + amp.acro({currentCovariate})") %>%
              as.formula(object = .)

              # obtain the data for the regression with missing values removed
              DataForRegression <-
                PhenotypeAndCovariateDataForCosinorRegression %>%
                select(
                  .data = .,
                  all_of(x = currentMetabolite),
                  all_of(x = currentCovariate),
                  DayOfYear
                  ) %>%
                na.omit(object = .)

              # obtain the regression model object
              RegressionModelObject <-
                DataForRegression %>%
                cosinor.lm(
                  formula = RegressionFormula,
                  period = 365,
                  data = .
                )

              # return the regression summary
              RegressionSummaryToReturn <-
                RegressionModelObject %>%
                summary(object = .)
              
              # return only the transformed cosinor regression summary with amplitude, acrophase, and annual mean
              RegressionSummaryToReturn <-
                 RegressionSummaryToReturn$transformed.table %>%
                  # compute the sum of squared residuals (SSR) and sum of squares total (SST)
                  mutate(
                    .data = .,
                    SSR = (RegressionModelObject$fit$residuals^2) %>% sum(.,na.rm = TRUE),
                    SST = (((DataForRegression %>% pull(.data = .,!!currentMetabolite))-mean(x = DataForRegression %>% pull(.data = .,!!currentMetabolite)))^2) %>%
                          sum(.,na.rm = TRUE)
                    ) %>%
                  # compute the R-squared = 1 - (SSR/SST) and round to 2 decimals
                  mutate(
                    .data = .,
                    R.squared = (1 - (SSR/SST)) %>% round(x = .,digits = 3)
                    ) %>%
                  # deselect the SSR and SST columns now that R-squared has been calculated
                  select(
                    .data = .,
                    -SSR,
                    -SST
                    ) %>%
                  # add a column with a label for the metabolite
                  mutate(
                    .data = .,
                    Metabolite = currentMetabolite
                    ) %>%
                  # add a column with a label for the covariate
                  mutate(
                    .data = .,
                    Covariate = currentCovariate
                    ) %>%
                  # add a column with the sample size included in the regression
                  mutate(
                    .data = .,
                    N = DataForRegression %>% nrow(x = .)
                    ) %>%
                  # create a column of the rownames
                  mutate(
                    .data = .,
                    RegressionTerm = rownames(x = RegressionSummaryToReturn$transformed.table)
                  ) %>%
                 # rearrange the columns
                  select(
                    .data = .,
                    RegressionTerm,
                    everything()
                  ) %>%
                # combine the estimate and standard error columns and round the values
                # and round the upper and lower confidence intervals
                mutate(
                  .data = .,
                  estimate_standard.error = paste(
                                                  signif(x = estimate,digits = 3),
                                                  "±",
                                                  signif(x = standard.error,digits = 3)
                                                  ),
                  upper.CI = upper.CI %>% signif(x = .,digits = 3),
                  lower.CI = lower.CI %>% signif(x = .,digits = 3)
                  )

              # remove the rownames, now that they have been converted to a column
              rownames(x = RegressionSummaryToReturn) <- NULL
              
            return(RegressionSummaryToReturn)
          }
          ) %>%
        do.call(
          what = "rbind",
          args = .
          )

      return(DataToReturn)
    }
      ) %>%
  do.call(
    what = "rbind",
    args = .
    ) %>%
  # round and add astericks to P-values
  mutate(
    .data = .,
    p.value = p.value %>%
              signif(x = .,digits = 3) %>%
              AddAstericksToPvalues(columnVector = .)
    ) %>%
  # group by the metabolite
  group_by(
    .data = .,
    Metabolite
    ) %>%
  # loop through the list of regression summaries by metabolite
  group_split(.tbl = .) %>%
  lapply(
    X = .,
    FUN = function(currentDataSet)
    {
      # identify the current metabolite
      currentMetabolite <-
        unique(currentDataSet$Metabolite)

      currentDataSet %>%
      tibble() %>%
      # save to the file system and label with the metabolite
      write_tsv(
        x = .,
        file = glue("./univariate_cosinor_regressions/univariate_cosinor_regression_{currentMetabolite}.tsv")
        )
    }
  )

########### summary statistics, stratified by season of sampling and gender:
          # code for finding grand average BMI, average age, gender count, and participant blood quanta 
          # code for finding average BMI, age, and count of each gender in each season of sampling,
          # code for finding average age, average BMI, average primary 25(OH)D3 and 25(OH)D [D3 + D2], and count of each gender

# create a directory for saving the demographic summary statistics
if(!dir.exists(paths = "./demographic_and_seasonal_summary_stats/"))
{
  dir.create(path = "./demographic_and_seasonal_summary_stats/")
}

# calculate grand summary statistics
SummaryStatisticsGrand <-
  PhenotypeAndCovariateData %>%
  # calculate the mean and standard deviation for age, blood quanta, BMI, 25OHD, and 25OHD3 concentration
  # calculate the counts of each gender also
  SummarizeDemographicData(DataFrameToUse = .) %>%
  # create a study season column and label it as grand so this table can be combined with the SummaryStatisticsBySeason table below
  mutate(
    .data = .,
    StudySeason = "grand"
    )

print(x = "######################## grand summary statistics ##################################")
print(x = SummaryStatisticsGrand)

# save the grand summary statistics to the file system
SummaryStatisticsGrand %>%
  write_tsv(
    x = .,
    file = "./demographic_and_seasonal_summary_stats/SummaryStatisticsGrand.tsv"
  )

# calculate summary statistics grouped by season
SummaryStatisticsBySeason <-
PhenotypeAndCovariateData %>%
  # group by the study season
  group_by(
    .data = .,
    StudySeason
    ) %>%
  # calculate the mean and standard deviation for age, blood quanta, BMI, 25OHD, and 25OHD3 concentration
  # calculate the counts of each gender also
  SummarizeDemographicData(DataFrameToUse = .) %>%
  # combine the summary statistics by season with the grand summary statistics by row
  list(
    .,
    SummaryStatisticsGrand
    ) %>%
  do.call(
    what = "rbind",
    args = .
      )
  

# save the seasonal demographic summary stats
SummaryStatisticsBySeason %>%
write_tsv(
  x = .,
  file = "./demographic_and_seasonal_summary_stats/SummaryStatisticsBySeason.tsv"
    )

# compute average age, BMI, overall 25(OH)D, and primary_25(OH)D3 concentration stratified by gender
SummaryStatisticsByGender <-
# select the relevant demographic data and 25(OH)D3 concentration
PhenotypeAndCovariateData %>%
  # group by gender 
  group_by(
    .data = .,
    Sex
    ) %>%
  # calculate the mean and standard deviation for age, blood quanta, BMI, 25OHD, and 25OHD3 concentration
  # calculate the counts of each gender also
  SummarizeDemographicData(DataFrameToUse = .)

# save the summary stats by gender to the file system
SummaryStatisticsByGender %>%
  write_tsv(
    x = .,
    file = "./demographic_and_seasonal_summary_stats/SummaryStatisticsByGender.tsv"
      )

# calculate percentage of participants that have:
    # (1) 25OHD levels below sufficiency (≤ 20 ng/mL)
    # (2) insufficient 25OHD levels (between 12 and 20 ng/mL)
    # (3) 25OHD deficiency (below 12 ng/mL)
VitaminDStatuses <-
PhenotypeAndCovariateData %>%
  select(
    .data = .,
    Overall_primary_25OHD_ngperml
    ) %>%
  na.omit(object = .) %>%
  # create binarized columns for BelowSufficiency, Insufficiency, and Deficiency
  mutate(
    .data = .,
    Insufficiency_percent = ((Overall_primary_25OHD_ngperml >= 12) & (Overall_primary_25OHD_ngperml <= 20)),
    Deficiency_percent = Overall_primary_25OHD_ngperml < 12
    ) %>%
  # deselect the concentration column now that it is no longer needed
  select(
    .data = .,
    -Overall_primary_25OHD_ngperml
    ) %>%
  # calculate the percentage for each category
  lapply(
    X = .,
    FUN = function(currentColumn)
    {
      SampleSize <-
        # find the sample size based on observations without missing data
        PhenotypeAndCovariateData %>%
        select(
          .data = .,
          Overall_primary_25OHD_ngperml
        ) %>%
        na.omit(object = .) %>%
        nrow(x = .)

      PercentToReturn <-
        round(x = (sum(currentColumn,na.rm = TRUE)/SampleSize)*100,digits = 2)

      return(PercentToReturn)
    }
  ) %>%
  do.call(
    what = "cbind",
    args = .
    ) %>%
  data.frame() %>%
  # calculate the percentage with vitamin D status below sufficiency by adding the Insufficiency_percent and Deficiency_percent
  mutate(
    .data = .,
    BelowSufficiency_percent = round(x = Insufficiency_percent + Deficiency_percent,digits = 2)
    ) %>%
  # move the BelowSufficiency_percent column to the front
  select(
    .data = .,
    BelowSufficiency_percent,
    everything()
  ) %>%
  # create a column with the group label as "All"
  mutate(
    .data = .,
    Group = "All"
    ) %>%
  # create a column with the total sample size
  mutate(
    .data = .,
    total_N = PhenotypeAndCovariateData %>%
              select(
                .data = .,
                Overall_primary_25OHD_ngperml
              ) %>%
              na.omit(object = .) %>%
              nrow(x = .)
    )

# calculate percentage of participants grouped by gender that have:
    # (1) 25OHD levels below sufficiency (≤ 20 ng/mL)
    # (2) insufficient 25OHD levels (between 12 and 20 ng/mL)
    # (3) 25OHD deficiency (below 12 ng/mL)
VitaminDStatusesByGender <-
  PhenotypeAndCovariateData %>%
    CalculateVitaminDstatusByGroup(
      DataFrameToUse = .,
      GroupingVariable = "Sex"
      )

# calculate percentage of participants grouped by season that have:
    # (1) 25OHD levels below sufficiency (≤ 20 ng/mL)
    # (2) insufficient 25OHD levels (between 12 and 20 ng/mL)
    # (3) 25OHD deficiency (below 12 ng/mL)
VitaminDStatusesBySeason <-
  PhenotypeAndCovariateData %>%
    CalculateVitaminDstatusByGroup(
      DataFrameToUse = .,
      GroupingVariable = "StudySeason"
      )

# bind vitamin D statuses for all participants, each gender, and each season together by row
VitaminDstatusPercentages <-
  list(
    VitaminDStatuses,
    VitaminDStatusesByGender,
    VitaminDStatusesBySeason
   ) %>%
  do.call(
    what = "rbind",
    args = .
    ) %>%
  data.frame() %>%
  # rearrange the columns
  select(
    .data = .,
    Group,
    total_N,
    BelowSufficiency_percent,
    Insufficiency_percent,
    Deficiency_percent
    )

  # save the results to the file system
  VitaminDstatusPercentages %>%
    write_tsv(
      x = .,
      file = "./demographic_and_seasonal_summary_stats/VitaminDstatusPercentages.tsv"
        )

print(x = "################# Vitamin D status percentages #####################")
print(x = VitaminDstatusPercentages)

############## univariate linear regressions for demographic and seasonal factors

# create an array of metabolites to use as independent variables in regression
Metabolites <-
PhenotypeAndCovariateData %>%
  names(x = .) %>%
  grepl(pattern = "(ngperml|ratio)",x = .) %>%
  which(x = .) %>%
  names(x = PhenotypeAndCovariateData)[.]

# create an array of covariates to use in univariate linear regression
Covariates <- 
  c("BMI","StudySeason","Sex","Age.on.Study.Date") %>% 
  sort(x = .)

UnivariateDemographicRegressionResults <-
    # loop through each metabolite and perform a regression with each covariate independently included as a covariate
    Metabolites %>%
    mclapply(
      X = .,
      FUN = function(currentMetabolite)
      {
            ResultsToReturn <-
              # loop through each of the covariates
             Covariates %>%
              lapply(
                X = .,
                FUN = function(currentCovariate)
                {
                  
                  # select columns for the current metabolite and current covariate
                  DataForRegression <-
                    PhenotypeAndCovariateData %>%
                    select(
                      .data = .,
                      "X" = all_of(x = currentCovariate),
                      "Y" = all_of(x = currentMetabolite)
                    ) %>%
                    # remove observations with missing values
                    na.omit(object = .)
                  
                  # create a robust standard errors regression model with the current metabolite and current covariate
                  # with the metabolite as the dependent variable
                  RegressionModel <-
                    DataForRegression %>%
                      lm_robust(
                        formula = Y ~ X,
                        data = .
                        )
                  
                  # perform an F-test for variance homogeneity to see if there is any heteroscedasticity with
                  # the lm() function; the lm_robust() function cannot be used in the ols_test_f() function
                  VarianceHomogeneityTestResult <-
                    DataForRegression %>%
                    lm(
                      formula = Y ~ X,
                      data = .
                    ) %>%
                    ols_test_f(
                      model = .
                      )
                  
                  ResultsToReturn <-
                    # create a tidy summary of the linear robust standard errors regression
                    TidyLinearRobustStandardErrorsRegressionResults(
                      DataForRegression = DataForRegression,
                      RegressionModel = RegressionModel
                      ) %>%
                      # add a column with the label for the current covariate
                    mutate(
                      .data = .,
                      Covariate_Included = rep_len(x = currentCovariate,length.out = nrow(x = .))
                      ) %>%
                      # add a column with the label for the current metabolite
                    mutate(
                      .data = .,
                      Metabolite_Included = currentMetabolite
                      ) %>%
                    mutate(
                        .data = .,
                        # the P-value from the test for Heteroscedastiticy with astericks for levels of significance
                        Heteroscedastiticy_p.value = VarianceHomogeneityTestResult$p %>%
                                                     signif(x = .,digits = 3) %>%
                                                     AddAstericksToPvalues(columnVector = .)
                        ) %>%
                        # add a column for the regression that was performed
                        mutate(
                          .data = .,
                          Regression = paste(
                                              currentMetabolite,
                                              "~",
                                              currentCovariate
                                            )
                        ) %>%
                        # add a column named "Confounder" and fill every entry with "none",
                        # this will make it easier to combine the UnivariateDemographicRegressionResults table with the 
                        # DemographicAndSeasonalConfoundingRegressionResults computed below to assess if demographic/seasonal covariates 
                        # are confounding vitamin D concentration associations
                        mutate(
                          .data = .,
                          Confounder = "none"
                          ) %>%
                        # add a regression type column and label it as univariate
                        mutate(
                          .data = .,
                          Regression_type = "univariate"
                          )
    
                  return(ResultsToReturn)
                }
                  ) %>%
              # bind all of the results together by row
              do.call(
                what = "rbind",
                args = .
                  )
    
            return(ResultsToReturn)
      },mc.cores = detectCores()
        ) %>%
      # bind all regression results for all metabolites and all covariates together by row
      do.call(
        what = "rbind",
        args = .
          ) %>%
      # filter out the intercept term from the regressions 
      filter(
        .data = .,
        term != "(Intercept)"
        )

# save the univariate demographic regression results to the file system
if(!dir.exists(paths = "./univariate_demographic_metabolite_regressions/"))
{
  dir.create(path = "./univariate_demographic_metabolite_regressions/")
}

UnivariateDemographicRegressionResults %>%
write_tsv(
  x = .,
  file = "./univariate_demographic_metabolite_regressions/UnivariateDemographicRegressionResults.txt"
  )

############ code for making plots of vitamin D concentration stratified by age, BMI, and gender

  # create a decision table for labeling age by decade
  AgeLabelDecisionTable <-
    data.frame(
      "AgeCutoff" = c(19,29,39,49,59,69,79,89),
      "DecadeLabel" = c("18-19","20-29","30-39","40-49","50-59","60-69","70-79","80-89")
    )
  
  # create a decision table for labeling BMI in units of 10
  BMIlabelDecisionTable <-
    data.frame(
      "BMIcutoff" = c(19,29,39,49,59,69),
      "BMIlabel" = c("10-19","20-29","30-39","40-49","50-59","60-69")
    )
  
# temporarily suppress non-problematic warning messages for removing missing values from the boxplots below
options(warn = -1)   
  
 ### serum 25OHD versus age summary stats
 # create a plot of serum 25(OH)D versus age by decade
DataFor25OHDversusAgePlot <-
 # select the primary vitamin D metabolite and the age on study date columns
 PhenotypeAndCovariateData %>%
   select(
     .data = .,
     Overall_primary_25OHD_ngperml,
     Age.on.Study.Date
     ) %>%
   na.omit(object = .) %>%
  # create a column with vitamin D status categories
    mutate(
      .data = .,
      primaryDstatus = if_else(
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
                                                                false = "Sufficient (> 20 & ≤ 50 ng/mL)"
                                                              )
                                                  )
                                    )
    ) %>%
   # create age by decade categories 
   mutate(
     .data = .,
     AgeCategories = Age.on.Study.Date %>%
                     lapply(
                       X = .,
                       FUN = function(currentAge)
                       {
                         
                           DecadeToReturn <-
                             AgeLabelDecisionTable %>%
                             # filter the AgeLabelDecisionTable to the AgeCutoff that
                             # is greater than or equal to the current age
                             filter(
                               .data = .,
                               AgeCutoff >= currentAge
                             ) %>%
                             # select the first row of the table
                             head(x = .,1) %>%
                             # pull out the decade label for that AgeCutoff
                             pull(
                               .data = .,
                               DecadeLabel
                             )
                           
                         return(DecadeToReturn)
                       }
                     ) %>%
                    unlist(x = .)
     ) %>%
   # convert the age by decade categories as factors
   mutate(
     .data = .,
     AgeCategories = AgeCategories %>% as.factor(x = .)
     )

Serum25OHDversusAgeBoxplot <-
  DataFor25OHDversusAgePlot %>%
   # create the boxplot
   ggplot(
     data = .,
     mapping = aes(
                   x = factor(x = AgeCategories),
                   y = Overall_primary_25OHD_ngperml
                      )
     ) +
   geom_boxplot(
     # make the boxes in the boxplot transparent
     alpha = 0.1,
     # remove the outliers on the boxplot since they are already present in the data
     outlier.shape = NA
     ) +
   geom_point(
     aes(
         color = factor(
                   x = primaryDstatus,
                   levels =
                           c(
                            "High (> 50 ng/mL)",
                            "Sufficient (> 20 & ≤ 50 ng/mL)",
                            "Insufficient (between 12 & 20 ng/mL)",
                            "Deficient (< 12 ng/mL)"
                           )
                     )
         ),
     position = position_jitter(width = 0.2,height = 0.1),
     size = 1.5
     ) +
   theme_classic() +
   labs(
     y = "Concentration (ng/mL)",
     x = "Age (years)"
     ) +
   scale_color_discrete(name = "Clinical Status") +
   theme(
     axis.title = element_text(family = "Arial",face = "bold",size = 11),
     axis.text.x = element_text(family = "Arial",face = "bold",size = 11),
     axis.text.y = element_text(family = "Arial",face = "bold",size = 11),
     legend.text = element_text(family = "Arial",face = "bold",size = 8),
     legend.title = element_text(family = "Arial",face = "bold",size = 8,hjust = 0.5),
     legend.position = c(.99, .99),
     legend.justification = c("right", "top"),
     legend.box.just = "right",
     legend.margin = margin(6, 6, 6, 6),
     legend.box.background = element_rect(color="black", size=1)
   ) +
   scale_x_discrete(
     # annotate the x axis with the age categories as well as the sample size for each category
     labels = DataFor25OHDversusAgePlot %>%
              group_by(
                .data = .,
                AgeCategories
                ) %>%
              summarise(
                .data = .,
                N = n()
                ) %>%
              mutate(
                .data = .,
                label = paste0(AgeCategories,"\n","(N = ",N,")")
                ) %>%
              pull(.data = .,label)
   ) +
   scale_y_continuous(
     breaks = c(0,12,20,30,40,50,60,70,80,90,100),
     labels = c(0,12,20,30,40,50,60,70,80,90,100),
     limits = c(0,100)
   ) # +
  # # label the plot with the beta-coefficient and p-value from the univariate regression
  # annotate(
  #   geom = "text",
  #   x = "20-29",
  #   y = 95,
  #   # select p-value and beta-coefficient labels for Overall_primary_25OHD_ngperml and Age.on.Study.Date
  #   label = UnivariateDemographicRegressionResults %>%
  #           filter(
  #             .data = .,
  #             (Metabolite_Included == "Overall_primary_25OHD_ngperml") & 
  #             (Covariate_Included == "Age.on.Study.Date") & 
  #             (term != "(Intercept)")
  #             ) %>%
  #           pull(
  #             .data = .,
  #             p.value
  #             ) %>%
  #           paste0(
  #             "P = ",
  #             .,
  #             "\n",
  #             "β = ",
  #             UnivariateDemographicRegressionResults %>%
  #               filter(
  #                 .data = .,
  #                 (Metabolite_Included == "Overall_primary_25OHD_ngperml") & 
  #                 (Covariate_Included == "Age.on.Study.Date") & 
  #                 (term != "(Intercept)")
  #               ) %>%
  #               pull(
  #                 .data = .,
  #                 estimate_std.error
  #               )
  #             ),
  #   fontface = "bold",
  #   size = 4,
  #   family = "Arial"
  #   ) 

 ### serum 25(OH)D versus BMI summary stats
 # select the primary vitamin D metabolite and the BMI on study date columns
Serum25OHDversusBMIBoxplotData <-
   PhenotypeAndCovariateData %>%
   select(
     .data = .,
     Overall_primary_25OHD_ngperml,
     BMI
   ) %>%
   na.omit(object = .) %>%
   # create BMI categories in increments of 10
   mutate(
     .data = .,
     BMIcategories = BMI %>%
                     lapply(
                            X = .,
                            FUN = function(currentBMIvalue)
                            {
                              
                              BMIlabelToReturn <-
                              BMIlabelDecisionTable %>%
                                # filter the BMIlabelDecisionTable to the BMI cutoff that is greater 
                                # than or equal to the currentBMIvalue
                                filter(
                                  .data = .,
                                  BMIcutoff >= currentBMIvalue
                                  ) %>%
                                # select the first row
                                head(x = .,1) %>%
                                # pull out the BMI label
                                pull(.data = .,BMIlabel)
                              
                              return(BMIlabelToReturn)
                            }
                              ) %>%
                            unlist(x = .)
   ) %>%
   # make sure BMI categories are factors
   mutate(
     .data = .,
     BMIcategories = BMIcategories %>% as.factor(x = .)
   ) %>%
  # create a column with vitamin D status categories
    mutate(
      .data = .,
      primaryDstatus = if_else(
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
                                                                false = "Sufficient (> 20 & ≤ 50 ng/mL)"
                                                              )
                                                  )
                                    )
    )

Serum25OHDversusBMIBoxplot <-
  Serum25OHDversusBMIBoxplotData %>%
   # create the boxplot
   ggplot(
     data = .,
     mapping = aes(
       y = Overall_primary_25OHD_ngperml,
       x = factor(x = BMIcategories)
     )
   ) +
  geom_boxplot(
    # remove the outliers on the boxplot since they are already present in the data
    outlier.shape = NA
    ) +
  geom_point(
    aes(
        color = factor(
                  x = primaryDstatus,
                  levels =
                          c(
                           "High (> 50 ng/mL)",
                           "Sufficient (> 20 & ≤ 50 ng/mL)",
                           "Insufficient (between 12 & 20 ng/mL)",
                           "Deficient (< 12 ng/mL)"
                          )
                    )
        ),
    position = position_jitter(width = 0.2,height = 0.1),
    size = 1.5
    ) +
  theme_classic() +
   labs(
     y = "Concentration (ng/mL)",
     x = bquote(bold('BMI'~(kg/m^2)))
     ) +
  scale_color_discrete(name = "Clinical Status") +
   theme(
     axis.title = element_text(family = "Arial",face = "bold",size = 11),
     axis.text.x = element_text(family = "Arial",face = "bold",size = 11),
     axis.text.y = element_text(family = "Arial",face = "bold",size = 11),
     legend.text = element_text(family = "Arial",face = "bold",size = 8),
     legend.title = element_text(family = "Arial",face = "bold",size = 8,hjust = 0.5),
     legend.position = c(.99, .99),
     legend.justification = c("right", "top"),
     legend.box.just = "right",
     legend.margin = margin(6, 6, 6, 6),
     legend.box.background = element_rect(color="black", size=1)
   ) +
   scale_x_discrete(
     labels = Serum25OHDversusBMIBoxplotData %>%
              group_by(
                .data = .,
                BMIcategories
                ) %>%
              summarize(
                .data = .,
                N = n()
                ) %>%
              mutate(
                .data = .,
                label = paste0(BMIcategories,"\n","(N = ",N,")")
                ) %>%
              pull(
                .data = .,
                label
                )
   ) +
   scale_y_continuous(
     breaks = c(0,12,20,30,40,50,60,70,80,90,100),
     labels = c(0,12,20,30,40,50,60,70,80,90,100),
     limits = c(0,100)
   ) # +
  # # label the plot with the beta-coefficient and p-value from the univariate regression
  # annotate(
  #   geom = "text",
  #   x = "10-19",
  #   y = 95,
  #   # select p-value and beta-coefficient labels for Overall_primary_25OHD_ngperml and BMI
  #   label = UnivariateDemographicRegressionResults %>%
  #           filter(
  #             .data = .,
  #             (Metabolite_Included == "Overall_primary_25OHD_ngperml") & 
  #             (Covariate_Included == "BMI") & 
  #             (term != "(Intercept)")
  #             ) %>%
  #           pull(
  #             .data = .,
  #             p.value
  #             ) %>%
  #           paste0(
  #             "P = ",
  #             .,
  #             "\n",
  #             "β = ",
  #             UnivariateDemographicRegressionResults %>%
  #               filter(
  #                 .data = .,
  #                 (Metabolite_Included == "Overall_primary_25OHD_ngperml") & 
  #                 (Covariate_Included == "BMI") & 
  #                 (term != "(Intercept)")
  #               ) %>%
  #               pull(
  #                 .data = .,
  #                 estimate_std.error
  #               )
  #             ),
  #   fontface = "bold",
  #   size = 4,
  #   family = "Arial"
  #   )
 
 ### serum 25(OH)D versus gender summary stats
 # select gender and 25(OH)D columns
 Serum25OHDgroupedByGenderBoxPlot <-
   PhenotypeAndCovariateData %>%
   select(
     .data = .,
     Sex,
     Overall_primary_25OHD_ngperml
   ) %>%
   na.omit(object = .) %>%
   # convert gender to a factor
   mutate(
     .data = .,
     Sex = Sex %>% as.factor(x = .)
   ) %>%
   # create a column with vitamin D status categories
     mutate(
       .data = .,
       primaryDstatus = if_else(
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
                                                                 false = "Sufficient (> 20 & ≤ 50 ng/mL)"
                                                               )
                                                   )
                                     )
     ) %>%
   # create a boxplot of 25(OH)D grouped by gender
   # create the boxplot
   ggplot(
     data = .,
     mapping = aes(
       y = Overall_primary_25OHD_ngperml,
       x = Sex
     )
   ) +
   geom_boxplot(
     # remove the outliers on the boxplot since they are already present in the data
     outlier.shape = NA
     ) +
   geom_point(
     aes(
         color = factor(
                   x = primaryDstatus,
                   levels =
                           c(
                            "High (> 50 ng/mL)",
                            "Sufficient (> 20 & ≤ 50 ng/mL)",
                            "Insufficient (between 12 & 20 ng/mL)",
                            "Deficient (< 12 ng/mL)"
                           )
                     )
         ),
     position = position_jitter(width = 0.2,height = 0.1),
     size = 1.5
     ) +
   theme_classic() +
   scale_color_discrete(name = "Clinical Status") +
   labs(y = "Concentration (ng/mL)",x = "Gender") +
   theme(
     axis.title = element_text(family = "Arial",face = "bold",size = 11),
     axis.text.x = element_text(family = "Arial",face = "bold",size = 11),
     axis.text.y = element_text(family = "Arial",face = "bold",size = 11),
     legend.text = element_text(family = "Arial",face = "bold",size = 8),
     legend.title = element_text(family = "Arial",face = "bold",size = 8,hjust = 0.5),
     legend.position = c(.95, .99),
     legend.justification = c("right", "top"),
     legend.box.just = "right",
     legend.margin = margin(6, 6, 6, 6),
     legend.box.background = element_rect(color="black", size=1)
   ) +
   scale_y_continuous(
     breaks = c(0,12,20,30,40,50,60,70,80,90,100),
     labels = c(0,12,20,30,40,50,60,70,80,90,100),
     limits = c(0,120)
   ) +
   scale_x_discrete(
     # add sample sizes to the axis labels
     labels = PhenotypeAndCovariateData %>%
              select(
                .data = .,
                Sex
                ) %>%
              group_by(
                .data = .,
                Sex
                ) %>%
              summarize(
                .data = .,
                N = n()
                ) %>%
              mutate(
                .data = .,
                label = paste0(Sex,"\n","(N = ",N,")")
                ) %>%
              pull(.data = .,label)
   ) # +
   # # label the plot with the beta-coefficient and p-value from the univariate regression
   # annotate(
   #   geom = "text",
   #   x = "Male",
   #   y = 70,
   #   # select p-value and beta-coefficient labels for Overall_primary_25OHD_ngperml and gender
   #   label = UnivariateDemographicRegressionResults %>%
   #           filter(
   #             .data = .,
   #             (Metabolite_Included == "Overall_primary_25OHD_ngperml") & 
   #             (Covariate_Included == "Sex") & 
   #             (term != "(Intercept)")
   #             ) %>%
   #           pull(
   #             .data = .,
   #             p.value
   #             ) %>%
   #           paste0(
   #             "P = ",
   #             .,
   #             "\n",
   #             "β = ",
   #             UnivariateDemographicRegressionResults %>%
   #               filter(
   #                 .data = .,
   #                 (Metabolite_Included == "Overall_primary_25OHD_ngperml") & 
   #                 (Covariate_Included == "Sex") & 
   #                 (term != "(Intercept)")
   #               ) %>%
   #               pull(
   #                 .data = .,
   #                 estimate_std.error
   #               )
   #             ),
   #   fontface = "bold",
   #   size = 4,
   #   family = "Arial"
   #   )
   
 
 # aggregate all of the boxplots as a list and save each one to the file system
 BoxplotsToSave <-
   list(
     "Serum25OHDversusAgeBoxplot" = Serum25OHDversusAgeBoxplot,
     "Serum25OHDversusBMIBoxplot" = Serum25OHDversusBMIBoxplot,
     "Serum25OHDgroupedByGenderBoxPlot" = Serum25OHDgroupedByGenderBoxPlot
     ) 
 
 mapply(
   FUN = function(currentPlot,nameOfCurrentPlot)
   {
     # create a directory for saving everything
     if(!dir.exists(paths = "./Serum25OHDDemographicBoxplots/"))
     {
       dir.create(path = "./Serum25OHDDemographicBoxplots/")
     }
     
     tiff(
       filename = glue("./Serum25OHDDemographicBoxplots/{nameOfCurrentPlot}.tiff"),
       width = 7,
       height = 5,
       units = "in",
       compression = "none",
       res = 300
       )
     
     print(x = currentPlot)
     dev.off()
     
    },
   BoxplotsToSave,
   names(x = BoxplotsToSave),
   SIMPLIFY = FALSE
   )
 
# turn warning messages back on
options(warn = 0)
 
##### test for pairwise associations of demographic and seasonal covariates to identify potential confounding relationships

  # A confounding relationship is when the following conditions are satisfied:
  # (1) covariates A and B are independently associated with 25(OH)D3 concentration
  # (2) covariates A and B are associated with each other
  # (3) when comparing linear regressions of outcome ~ covariate A and outcome ~ covariate A + covariate B,
        # the effect size of outcome ~ covariate A changes by 10% or more with the addition of covariate B

  # the tests that are performed here are:
    # (1) regression of BMI versus age, gender, and season independently
    # (2) regression of Age versus gender and season independently
    # (3) chi-squared test for independence of gender count and season of sampling
 
AgeAndBMIVersusOtherCovariatesRegressionResults <-
  # loop through the two possible outcomes (age and BMI), the dependent variables
  c("Age.on.Study.Date","BMI") %>%
  lapply(
    X = .,
    FUN = function(currentOutcome)
    {
      RegressionResultsForCurrentOutcome <-
      # loop through the possible covariates that could be associated with the outcome, the independent variables
      c(
        "Age.on.Study.Date",
        "BMI",
        "Sex",
        "StudySeason"
        ) %>%
      lapply(
        X = .,
        FUN = function(currentCovariate)
        {

          DataToIncludeInRegression <-
          PhenotypeAndCovariateData %>%
            # select the outcome column and the column of the current covariate
            select(
              .data = .,
              "independentVariable" = all_of(x = currentCovariate),
              "dependentVariable" = all_of(x = currentOutcome)
              ) %>%
            # remove observations with missing values
            na.omit(object = .)
          
          # obtain the robust standard errors regression model
          RegressionModel <-
            DataToIncludeInRegression %>%
              # perform the regression of the currentOutcome (dependent variable) versus the currentCovariate (independent variable)
              lm_robust(
                formula = dependentVariable ~ independentVariable,
                data = .
                ) 
          
          # test for variance heteroscedasticity using the lm() function. 
          # the lm_robust() function cannot be used with the ols_test_f() function
          VarianceHomogeneityTestResult <-
            DataToIncludeInRegression %>%
            lm(
              formula = dependentVariable ~ independentVariable,
               data = .
               ) %>%
            ols_test_f(
              model = .
              )
            
          # obtain the regression summary
            RegressionSummary <-
              RegressionModel %>%
              # summarize the regression
              summary(object = .)
            
          # create a tidy table of regression results to return with
          # beta-coefficients, standard errors, p-values, unadjusted and adjusted R-squared, and sample size
          RegressionResultsToReturn <-
            TidyLinearRobustStandardErrorsRegressionResults(
              DataForRegression = DataToIncludeInRegression,
              RegressionModel = RegressionModel
                ) %>%
            # remove the intercept term
            filter(
              .data = .,
              term != "(Intercept)"
              ) %>%
            # create a column with a label for the current outcome
            # create a column with a label for the current covariate,
            mutate(
              .data = .,
              Outcome = currentOutcome,
              Covariate = currentCovariate,
              # the P-value from the test for Heteroscedastiticy with astericks added for level of significance
              Heteroscedastiticy_p.value = VarianceHomogeneityTestResult$p %>% 
                                            signif(x = .,digits = 3) %>%
                                            AddAstericksToPvalues(columnVector = .)
            ) %>%
            # remove rows where the outcome and the covariate are equivalent
            filter(
              .data = .,
              Outcome != Covariate
              ) 
          
          return(RegressionResultsToReturn)
        }
          ) %>%
        do.call(
          what = "rbind",
          args = .
          )

      return(RegressionResultsForCurrentOutcome)
    }
    ) %>%
  do.call(
    what = "rbind",
    args = .
    ) %>%
  # select useful columns, rearrange, and rename them so they are easier to read
  select(
    .data = .,
    Outcome,
    Covariate,
    term,
    "beta_std.error" = estimate_std.error,
    p.value,
    F_stat,
    F_stat_p_value,
    R2_adj,
    N,
    Heteroscedastiticy_p.value
    ) 

# the final test that is needed to analyze for confounding is differences in gender count in each season of sampling
# this is accomplished by performing a chi-squared test of independence between gender count and season
DataForChiSquaredTest <-
  # select the gender count columns from the SummaryStatisticsBySeason table
  SummaryStatisticsBySeason %>%
  select(
    .data = .,
    Females_n,
    Males_n
    ) %>%
  # make sure the table is a matrix
  as.matrix(x = .)

  # name the rows of the matrix with the corresponding season,
  # a matrix format is required for the chi-squared test
rownames(x = DataForChiSquaredTest) <-
  SummaryStatisticsBySeason$StudySeason %>%
  as.character(x = .)

# perform the chi-squared test of independence of gender count and season
ChiSquareTestResultsGenderVersusSeason <-
  DataForChiSquaredTest %>%
  chisq.test(x = .) %>%
  broom::tidy(x = .) %>%
  # round the p.value and add astericks for significance
  mutate(
    .data = .,
    p.value = p.value %>% 
              signif(x = .,digits = 3) %>% 
              AddAstericksToPvalues(columnVector = .)
    ) %>%
  select(
    .data = .,
    "chiSq" = statistic,
    p.value,
    "df" = parameter,
    method
    )

# create a directory for saving the pairwise associations of demographic and seasonal covariates
if(!dir.exists(paths = "./pairwiseDemographicAndSeasonalCovariateAssociations/"))
{
  dir.create(path = "./pairwiseDemographicAndSeasonalCovariateAssociations/")
}

# save the AgeAndBMIVersusOtherCovariatesRegressionResults and ChiSquareTestResultsGenderVersusSeason to the file system
AgeAndBMIVersusOtherCovariatesRegressionResults %>%
write_tsv(
  x = .,
  file = "./pairwiseDemographicAndSeasonalCovariateAssociations/AgeAndBMIVersusOtherCovariatesRegressionResults.txt"
  )

ChiSquareTestResultsGenderVersusSeason %>%
write_tsv(
  x = .,
  file = "./pairwiseDemographicAndSeasonalCovariateAssociations/ChiSquareTestResultsGenderVersusSeason.txt"
  )

# Final tests for confounding associations with 25(OH)D3 concentration as well as all other metabolites:
  # Example for 25(OH)D3:
      # since age was associated with 25(OH)D3 and age was associated with BMI,
      # need to see if BMI confounds the association of age with 25(OH)D3
          # Thus, compare [25(OH)D3 ~ age] to [25(OH)D3 ~ age + BMI]

DemographicAndSeasonalConfoundingRegressionResults <-
# loop through the all of the potential covariates
c(
  "Age.on.Study.Date",
  "BMI",
  "Sex",
  "StudySeason"
  ) %>%
  lapply(
    X = .,
    FUN = function(currentCovariateOfInterest)
    {
      # loop through the candidate confounders
      c(
        "Age.on.Study.Date",
        "BMI",
        "Sex",
        "StudySeason"
        ) %>%
        lapply(
          X = .,
          FUN = function(currentPotentialConfounder)
          {
            # loop through all of the metabolites 
            Metabolites %>%
              lapply(
                X = .,
                FUN = function(currentMetabolite)
                {
                  # perform the regression only in cases when the currentCovariateOfInterest is not the same as the currentPotentialConfounder
                  # because a variable cannot confound itself
                  if(currentCovariateOfInterest!=currentPotentialConfounder)
                  {
                    # select columns corresponding to the current metabolite, the currentCovariateOfInterest, and the currentPotentialConfounder
                    # this will be the data for the regression of the metabolite versus the covariate plus the confounder
                    DataToIncludeInRegressions <-
                      PhenotypeAndCovariateData %>%
                      select(
                        .data = .,
                        "Metabolite" = all_of(x = currentMetabolite),
                        "Covariate" = all_of(x = currentCovariateOfInterest),
                        "Confounder" = all_of(x = currentPotentialConfounder)
                      ) %>%
                      # remove observations with missing values
                      na.omit(object = .)
                    
                    # perform the robust standard errors regression of the metabolite versus Covariate + Confounder
                    RegressionModel <-
                      DataToIncludeInRegressions %>%
                      lm_robust(
                        formula = Metabolite ~ Covariate + Confounder,
                        data = .
                      )
                    
                    # tidy the BivariateRegressionSummary and add useful test statistic columns
                    TidyRegressionSummaries <-
                            TidyLinearRobustStandardErrorsRegressionResults(
                              DataForRegression = DataToIncludeInRegressions,
                              RegressionModel = RegressionModel
                              ) %>%
                            # remove the intercept term 
                            filter(
                              .data = .,
                              term != "(Intercept)"
                              ) %>%
                            # label the confounder with the actual confounder that was included
                            mutate(
                              .data = .,
                              term = if_else(
                                             condition = term == "Confounder",
                                             true = currentPotentialConfounder,
                                             false = term
                                               )
                              ) %>%
                          # label the covariate with the actual covariate that was included
                          mutate(
                            .data = .,
                            term = if_else(
                                          condition = term == "Covariate",
                                          true = currentCovariateOfInterest,
                                          false = term
                                            )
                            ) %>%
                        # add a column for the regression that was performed
                        mutate(
                          .data = .,
                          Regression = paste(
                                            currentMetabolite,
                                            "~",
                                            currentCovariateOfInterest,
                                            "+",
                                            currentPotentialConfounder
                                            )
                            ) %>%
                        # add a column with the metabolite included
                        mutate(
                          .data = .,
                          Metabolite_Included = currentMetabolite
                          ) %>%
                        # add a column with the covariate included
                        mutate(
                          .data = .,
                          Covariate_Included = currentCovariateOfInterest
                          ) %>%
                        # add a column with the potential confounder included
                        mutate(
                          .data = .,
                          Confounder = currentPotentialConfounder
                          ) %>%
                        # add a regression type column and label it as bivariate
                        mutate(
                          .data = .,
                          Regression_type = "bivariate"
                        )
                    
                    return(TidyRegressionSummaries)

                  }
                  
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
          )
    }
    ) %>%
  do.call(
    what = "rbind",
    args = .
    )

# create a directory for saving the DemographicAndSeasonalConfoundingRegressionResults
if(!dir.exists(paths = "./DemographicAndSeasonalConfoundingRegressionResults/"))
{
  dir.create(path = "./DemographicAndSeasonalConfoundingRegressionResults/")
}

DemographicAndSeasonalConfoundingRegressionResults %>%
  write_tsv(
    x = .,
    file = "./DemographicAndSeasonalConfoundingRegressionResults/DemographicAndSeasonalConfoundingRegressionResults.txt"
      )

# bind the UnivariateDemographicRegressionResults and the DemographicAndSeasonalConfoundingRegressionResults together by row
UnivariateAndBivariateMetaboliteRegressionAnalysisForPotentialConfounding <-
  list(
    UnivariateDemographicRegressionResults %>%
      select(
        .data = .,
        -Heteroscedastiticy_p.value
        ),
    DemographicAndSeasonalConfoundingRegressionResults
    ) %>%
    do.call(
      what = "rbind",
      args = .
      ) %>%
    data.frame() %>%
    # group by the dependent metabolite variable and the independent covariate variable
    group_by(
      .data = .,
      Metabolite_Included,
      Covariate_Included
      ) %>%
    # split the dataframe into a list of dataframes based on the metabolite variable and the independent covariate variable
    group_split(
      .tbl = .
      ) %>%
    # bind the list of dataframes back into a dataframe now that is has been sorted by metabolite variable and the independent covariate variable
    do.call(
      what = "rbind",
      args = .
      ) %>%
    data.frame()

  # save the final results to the file system
  UnivariateAndBivariateMetaboliteRegressionAnalysisForPotentialConfounding %>%
    write_tsv(
      x = .,
      file = "./DemographicAndSeasonalConfoundingRegressionResults/UnivariateAndBivariateMetaboliteRegressionAnalysisForPotentialConfounding.txt"
    )
  
  # identify the causative associated variants by filtering the clump data to variants with a bonferroni adjusted P-value of 0.05 or less
  CausativeAssociatedVariants <-
    ClumpData %>%
    # group the results by metabolite to create a column with bonferroni adjusted P-values
    group_by(
      .data = .,
      metabolite
      ) %>%
    # split into a list of dataframes by metabolite
    group_split(.tbl = .) %>%
    # loop through each dataframe and compute bonferroni adjusted P-values
    mclapply(
      X = .,
      FUN = function(currentMetaboliteAssociations)
      {
        
        DataToReturn <-
        currentMetaboliteAssociations %>%
          # the bonferroni adjusted P-values are computed by multiplying the unadjusted P-value 
          # by the number of independent causative variants that were tested
          mutate(
            .data = .,
            P_bonferroni = nrow(x = currentMetaboliteAssociations)*P
          ) %>%
          # add a column with the number of snps used for the bonferroni adjustment
          mutate(
            .data = .,
            Bonferroni_multiplier = nrow(x = currentMetaboliteAssociations)
            )

        return(DataToReturn)
      },mc.cores = detectCores()
      ) %>%
    # bind the list of dataframes back together by row
    do.call(
      what = "rbind",
      args = .
      ) %>%
    # filter to results with a Bonferroni adjusted P-value of less than or equal to 0.05
    filter(
      .data = .,
      P_bonferroni <= 0.05
      ) %>%
    # select relevant columns
    select(
      .data = .,
      CHR,
      SNP,
      Gene,
      metabolite,
      P_bonferroni,
      Bonferroni_multiplier
    ) %>%
    # join the clump results with the rsIDlookupTable
    left_join(
      x = .,
      y = rsIDlookupTable,
      by = "SNP"
        ) %>%
    # group by gene and metabolite
    group_by(
      .data = .,
      Gene,
      metabolite
      ) %>%
    # split into a list of dataframes by gene and metabolite
    group_split(.tbl = .) %>%
    mclapply(
      X = .,
      FUN = function(currentDataFrame)
      {

        DataToReturn <-
        currentDataFrame %>%
        # join the causative associated variants with the linear regression data with annotation data from
        # the original additive association with PLINK to obtain beta-coefficients, standard errors, and sample sizes
        left_join(
          x = .,
          y = "./linear_regression_results_file_and_annotation_common_variants.txt" %>%
              # load the regression results file and annotation data from the original additive association
              read.delim(
                file = .,
                header = TRUE
              ) %>%
              # add a space before the intron and exon number so the format is not changed to a date when opened in microsoft excel
              mutate(
                .data = .,
                intron_number_VEP = intron_number_VEP %>%
                                    paste0(" ",.),
                exon_number_VEP = exon_number_VEP %>%
                                  paste0(" ",.)
                ) %>%
              # select the SNP, beta coefficient, standard error, sample size,  metabolite, gene, and unadjusted P-value from the PLINK association
              select(
                .data = .,
                SNP,
                BETA,
                SE,
                "Pvalue_unadjust_plink" = Pvalue,
                NMISS,
                Metabolite,
                "Variant ID" = existing_variant_VEP,
                "Star Allele" = starAllele,
                "cDNA" = HGVSc,
                "Exon" = exon_number_VEP,
                "Intron" = intron_number_VEP,
                "Gene" = Gene_VEP,
                "Consequence" = Consequence_VEP,
                "Amino Acid" = HGVSp,
                "CADD (Phred scale)" = CADD_PHRED_VEP,
                "ADME Framework Score" = ADME_optimized_prediction,
                "PharmGKB Evidence" = pharmGKB_EvidenceLevel,
                "AAF**" = CSKT_AF,
                "ALL" = X1000g2015aug_all,
                "AFR" = X1000g2015aug_afr,
                "AMR" = X1000g2015aug_amr,
                "EAS" = X1000g2015aug_eas,
                "EUR" = X1000g2015aug_eur,
                "SAS" = X1000g2015aug_sas,
                "ClinVar Descript." = CLNDN,
                "ClinVar Signif." = CLNSIG,
                "Pubmed ID" = Pubmed_ID
                ) %>%
              # filter to the current gene and metabolite that is in the current data frame
              filter(
                .data = .,
                (Metabolite == unique(x = currentDataFrame$metabolite)) &
                  (Gene == unique(x = currentDataFrame$Gene))
                ) %>%
               # deselect the metabolite column and gene column so there aren't duplicate columns in the final table
              select(
                .data = .,
                -Metabolite,
                -Gene
              ),
          by = "SNP"
          )

        return(DataToReturn)
      },mc.cores = detectCores()
        ) %>%
    # bind the results back together by row
    do.call(
      what = "rbind",
      args = .
      ) %>%
    # round the Bonferroni adjusted P-values and the undajusted P-values and add astericks with levels of significance
    mutate(
      .data = .,
      P_bonferroni = P_bonferroni %>% signif(x = .,digits = 3) %>% AddAstericksToPvalues(columnVector = .),
      Pvalue_unadjust_plink = Pvalue_unadjust_plink %>% signif(x = .,digits = 3) %>% AddAstericksToPvalues(columnVector = .),
      # paste the beta coefficient and standard error columns together
      BETA_SE = paste0(BETA,"±",SE)
      ) 

  # save the causative associated variants to the file system
  if(!dir.exists(paths = "./CausativeVariantAssociations/"))
  {
    dir.create(path = "./CausativeVariantAssociations/")
  }

  CausativeAssociatedVariants %>%
    write_tsv(
      x = .,
      "./CausativeVariantAssociations/CausativeIndependentVariants.tsv"
      )
  
# Perform univariate, robust standard errors linear regressions of phenotype versus genotype for significant variant/metabolite pairs only.
# These regression results need to be compared to regressions of phenotype ~ genotype + each individual demographic and seasonal covariate to see if 
# there is any confounding going on
UnivarateGenotypeRegressionsWithRobustStandardErrors <-
      list.files(path = "./genotype_matrices_with_metabolite_data/") %>%
      lapply(
        X = .,
        FUN = function(currentFileName)
        {

          GentoypeMatrixWithMetaboliteData <-
            # load the genotype matrix with metabolite data added
            read.delim(
              file = glue("./genotype_matrices_with_metabolite_data/{currentFileName}"),
              header = TRUE
            )

         # identify the variants that are in the genotype matrix
          VariantsInGenotypeMatrix <-
            GentoypeMatrixWithMetaboliteData %>%
            names(x = .) %>%
            grepl(pattern = "X",x = .) %>%
            which(x = .) %>%
            names(x = GentoypeMatrixWithMetaboliteData)[.]
          
          # identify the variants in the genotype matrix that have significant associations with a metabolite in primary analysis
          VariantsToLoopThrough <-
            (VariantsInGenotypeMatrix %in% CausativeAssociatedVariants$VariantID_GenotypeMatrixFormat) %>%
            which(x = .) %>%
            VariantsInGenotypeMatrix[.]

          # identify the current gene based on the file name
          currentGene <-
            currentFileName %>%
            str_extract(
              string = .,
              pattern = "(CYP3A4|DHCR7|CYP2R1|CYP27B1|SULT2A1|UGT1A4|GC|VDR|CASR|CUBN|CYP24A1|RXRG|RXRB|RXRA)"
                )

          print(x = currentGene)

          RegressionResults <-
           # perform a univariate linear regression with robust standard errors of phenotype ~ genotype
           # with each significant variant in the genotype matrix with all appropriate metabolites
          VariantsToLoopThrough %>%
            mclapply(
              X = .,
              FUN = function(currentVariantToIncludeInRegression)
              {
                print(x = currentVariantToIncludeInRegression)
                
                # identify the metabolites to loop through based on the currentVariantToIncludeInRegression
                MetabolitesToLoopThrough <-
                  CausativeAssociatedVariants %>%
                    filter(
                      .data = .,
                      VariantID_GenotypeMatrixFormat == currentVariantToIncludeInRegression
                      ) %>%
                    pull(
                      .data = .,
                      metabolite
                      )

                MetabolitesToLoopThrough %>%
                  lapply(
                    X = .,
                    FUN = function(currentMetabolite)
                    {

                        print(x = currentMetabolite)


                      DataToUseInRegression <-
                        # select the columns of the metabolite and the current variant
                        GentoypeMatrixWithMetaboliteData %>%
                        select(
                          .data = .,
                          "Phenotype" = all_of(currentMetabolite),
                          "Genotype" = all_of(currentVariantToIncludeInRegression)
                        ) %>%
                        # ensure the genotype is continuous
                        mutate(
                          .data = .,
                          Genotype = Genotype %>% as.numeric(x = .)
                        ) %>%
                        # ensure the phenotype is a continuous numeric
                        mutate(
                          .data = .,
                          Phenotype = Phenotype %>% as.numeric(x = .)
                        ) %>%
                        # remove any missing values
                        na.omit(object = .)

                      # if there is more than one genotype value, perform the regression (e.g., 0 and 1)
                      # if there is only one genotype value (e.g., 0 only),
                        # return a NULL
                      if(
                        length(x = unique(x = as.numeric(x = DataToUseInRegression$Genotype)))>1
                        )
                      {

                        # obtain the robust standard errors regression model
                        RegressionModel <-
                            DataToUseInRegression %>%
                            # define the robust standard errors regression model
                            lm_robust(
                              formula = Phenotype ~ Genotype,
                              data = .
                            )

                         # perform an F-test to see if there is variance homogeneity in the regression of Phenotype ~ Genotype
                         # The lm() function must be used for this, not the lm_robust() function
                        # test for variance heteroscedasticity
                        VarianceHomogeneityTestResult <-
                          ols_test_f(
                                     model = lm(
                                                formula = Phenotype ~ Genotype,
                                                data = DataToUseInRegression
                                                  )
                                       )

                        RegressionResultsToReturn <-
                          # tidy the model output
                          TidyLinearRobustStandardErrorsRegressionResults(
                            DataForRegression = DataToUseInRegression,
                            RegressionModel = RegressionModel
                            ) %>%
                          # create a column label for the current gene
                          mutate(
                            .data = .,
                            Gene = rep_len(x = currentGene,length.out = nrow(x = .))
                          ) %>%
                          # create a column with a label for the metabolite
                          mutate(
                            .data = .,
                            metabolite = rep_len(x = currentMetabolite,length.out = nrow(x = .))
                          ) %>%
                          # create a column with a label for the variant
                          mutate(
                            .data = .,
                            variant = rep_len(x = currentVariantToIncludeInRegression,length.out = nrow(x = .))
                          ) %>%
                          # create a column with a label for the covariate (none in this case)
                          mutate(
                            .data = .,
                            Covariate = "none"
                            ) %>%
                          # create a column with a label for the regression
                          mutate(
                            .data = .,
                            Regression = paste0(
                                                currentMetabolite,
                                                " ~ ",
                                                currentVariantToIncludeInRegression
                                                )
                            ) %>%
                          # create a column with a label for the regression type
                          mutate(
                            .data = .,
                            Regression_type = "univariate"
                            ) %>%
                        mutate(
                          .data = .,
                          # add the p-value for the F Test for Heteroskedasticity of Phenotype ~ Genotype
                          Heteroscedastiticy_genotype_p.value = VarianceHomogeneityTestResult$p %>%
                                                                signif(x = .,digits = 3) %>%
                                                                AddAstericksToPvalues(columnVector = .)
                          )

                      } else {
                        RegressionResultsToReturn <- NULL
                      }

                      return(RegressionResultsToReturn)

                    }
                      ) %>%
                  do.call(
                    what = "rbind",
                    args = .
                    )


              },mc.cores = detectCores()
                ) %>%
            do.call(
              what = "rbind",
              args = .
              )

          return(RegressionResults)
        }
          ) %>%
        do.call(
          what = "rbind",
          args = .
          ) %>%
      filter(
        .data = .,
        term != "(Intercept)"
        )

# Perform robust standard errors linear regressions of phenotype versus genotype for significant variant/metabolite pairs only plus 
# each individual demographic and seasonal covariate to see if the relationship between phenotype and genotype is being confounded

BivariateGenotypeRegressions <-
  list.files(path = "./genotype_matrices_with_metabolite_data/") %>%
  lapply(
    X = .,
    FUN = function(currentFileName)
    {
  
      GentoypeMatrixWithMetaboliteAndCovariateData <-
        read.delim(
          file = glue("./genotype_matrices_with_metabolite_data/{currentFileName}"),
          header = TRUE
        ) %>%
        # upon loading the genotype matrix join it with the covariate data by sequencing ID
        left_join(
          x = .,
          y = CovariateData,
          by = c("sampleID" = "VCFfileID")
            )
  
     # identify the variants that are in the genotype matrix
      VariantsInGenotypeMatrix <-
        GentoypeMatrixWithMetaboliteAndCovariateData %>%
        names(x = .) %>%
        grepl(pattern = "X",x = .) %>%
        which(x = .) %>%
        names(x = GentoypeMatrixWithMetaboliteAndCovariateData)[.]
      
      # identify the variants in the genotype matrix that have significant associations with a metabolite in primary analysis
      VariantsToLoopThrough <-
        (VariantsInGenotypeMatrix %in% CausativeAssociatedVariants$VariantID_GenotypeMatrixFormat) %>%
        which(x = .) %>%
        VariantsInGenotypeMatrix[.]
  
      # identify the current gene based on the file name
      currentGene <-
        currentFileName %>%
        str_extract(
          string = .,
          pattern = "(CYP3A4|DHCR7|CYP2R1|CYP27B1|SULT2A1|UGT1A4|GC|VDR|CASR|CUBN|CYP24A1|RXRG|RXRB|RXRA)"
            )
  
      print(x = currentGene)
  
      RegressionResults <-
      # perform a multiple linear regression with robust standard errors with each 
      # significant variant in the genotype matrix with all appropriate metabolites
      VariantsToLoopThrough %>%
        mclapply(
          X = .,
          FUN = function(currentVariantToIncludeInRegression)
          {
            print(x = currentVariantToIncludeInRegression)
            
            # identify the metabolites to loop through based on the currentVariantToIncludeInRegression
            MetabolitesToLoopThrough <-
            CausativeAssociatedVariants %>%
              filter(
                .data = .,
                VariantID_GenotypeMatrixFormat == currentVariantToIncludeInRegression
                ) %>%
              pull(
                .data = .,
                metabolite
                )
  
            MetabolitesToLoopThrough %>%
              lapply(
                X = .,
                FUN = function(currentMetabolite)
                {
                  
                    print(x = currentMetabolite)
                    
                    
                      # loop through each individual seasonal and demographic covariate and include
                      # them in the metabolite regression with genotype
                      c("Age.on.Study.Date","StudySeason","BMI","Sex") %>%
                      sort(x = .) %>%
                        lapply(
                          X = .,
                          FUN = function(currentCovariate)
                          {
                            
                              DataToUseInRegression <-
                                # select the columns of the metabolite, the current variant, and the current covariate
                                GentoypeMatrixWithMetaboliteAndCovariateData %>%
                                select(
                                  .data = .,
                                  "Phenotype" = all_of(currentMetabolite),
                                  "Genotype" = all_of(currentVariantToIncludeInRegression),
                                  "Covariate" = all_of(currentCovariate)
                                ) %>%
                                # ensure the genotype is continuous
                                mutate(
                                  .data = .,
                                  Genotype = Genotype %>% as.numeric(x = .)
                                ) %>%
                                # ensure the phenotype is a continuous numeric
                                mutate(
                                  .data = .,
                                  Phenotype = Phenotype %>% as.numeric(x = .)
                                ) %>%
                                # remove any missing values
                                na.omit(object = .)
                              
                              # if the current covariate is StudySeason or Sex, make sure it is categorical
                              if(currentCovariate == "StudySeason" | currentCovariate == "Sex")
                              {
                                DataToUseInRegression <-
                                  DataToUseInRegression %>%
                                  mutate(
                                    .data = .,
                                    Covariate = Covariate %>%
                                                as.factor(x = .)
                                    )
                              }

                              # if the current covariate is Age.on.Study.Date or BMI, make sure it is numeric
                              if(currentCovariate == "Age.on.Study.Date" | currentCovariate == "BMI")
                              {
                                DataToUseInRegression <-
                                  DataToUseInRegression %>%
                                  mutate(
                                    .data = .,
                                    Covariate = Covariate %>%
                                                as.numeric(x = .)
                                  )
                              }
                              
                              # if there is more than one genotype value (e.g., 0 and 1) and the covariate is continuous,
                                 # perform the regression
                              # if there is more than one genotype value (e.g., 0 and 1) and the categorical covariate has more than one group,
                                 # perform the regression
                              # otherwise, return a NULL because the regression can't be performed
                              if(
                                  (
                                    (length(x = unique(x = as.numeric(x = DataToUseInRegression$Genotype)))>1) &
                                    (is.numeric(x = DataToUseInRegression$Covariate))
                                  ) |
                                  (
                                    (length(x = unique(x = as.numeric(x = DataToUseInRegression$Genotype)))>1) &
                                    (length(x = unique(x = as.factor(x = DataToUseInRegression$Covariate)))>1)
                                  )
                                )
                              {

                                # obtain the robust standard errors regression model
                                RegressionModel <-
                                    DataToUseInRegression %>%
                                    # define the robust standard errors multiple regression model 
                                    lm_robust(
                                      formula = Phenotype ~ Genotype + Covariate,
                                      data = .
                                    )

                                RegressionResultsToReturn <-
                                  # tidy the model output 
                                  TidyLinearRobustStandardErrorsRegressionResults(
                                    DataForRegression = DataToUseInRegression,
                                    RegressionModel = RegressionModel
                                    ) %>%
                                  # create a column label for the current gene
                                  mutate(
                                    .data = .,
                                    Gene = rep_len(x = currentGene,length.out = nrow(x = .))
                                  ) %>%
                                  # create a column with a label for the metabolite
                                  mutate(
                                    .data = .,
                                    metabolite = rep_len(x = currentMetabolite,length.out = nrow(x = .))
                                  ) %>%
                                  # create a column with a label for the variant
                                  mutate(
                                    .data = .,
                                    variant = rep_len(x = currentVariantToIncludeInRegression,length.out = nrow(x = .))
                                  ) %>%
                                  # create a column with a label for the covariate
                                  mutate(
                                    .data = .,
                                    Covariate = currentCovariate
                                    ) %>%
                                  # create a column with the regression performed
                                  mutate(
                                    .data = .,
                                    Regression = paste0(
                                                        currentMetabolite,
                                                        " ~ ",
                                                        currentVariantToIncludeInRegression,
                                                        " + ",
                                                        currentCovariate
                                                        )
                                    ) %>%
                                  # create a column with a label for the regression type
                                  mutate(
                                    .data = .,
                                    Regression_type = "bivariate"
                                    ) 

                              } else {
                                RegressionResultsToReturn <- NULL
                              }

                              return(RegressionResultsToReturn)
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
                )
  
  
          },mc.cores = detectCores()
            ) %>%
        do.call(
          what = "rbind",
          args = .
          )
  
      return(RegressionResults)
    }
      ) %>%
    do.call(
      what = "rbind",
      args = .
      ) %>%
    filter(
      .data = .,
      term != "(Intercept)"
      )

# bind the UnivariateGenotypeRegressions and the BivariateGenotypeRegressions together by row
UnivariateAndBivariateRegressionAnalysisForPotentialGenotypeConfoundingWithRobustStandardErrors <-
  list(
    UnivarateGenotypeRegressionsWithRobustStandardErrors %>%
      select(
        .data = .,
        -Heteroscedastiticy_genotype_p.value
        ),
    BivariateGenotypeRegressions
    ) %>%
    do.call(
      what = "rbind",
      args = .
      ) %>%
    data.frame() %>%
    # group by the metabolite and the variant
    group_by(
      .data = .,
      metabolite,
      variant
      ) %>%
    # split the dataframe into a list of dataframes based on the metabolite and the variant
    group_split(
      .tbl = .
      ) %>%
    # bind the list of dataframes back into a dataframe now that is has been sorted by metabolite and variant
    do.call(
      what = "rbind",
      args = .
      ) %>%
    data.frame()

# create a file for saving genotype confounding regression results if it does not exist
if(!dir.exists(paths = "./GenotypeConfoundingRegressionAnalysis/"))
{
  dir.create(path = "./GenotypeConfoundingRegressionAnalysis/")
}

# save the final results to the file system
UnivariateAndBivariateRegressionAnalysisForPotentialGenotypeConfoundingWithRobustStandardErrors %>%
  write_tsv(
    x = .,
    file = "./GenotypeConfoundingRegressionAnalysis/UnivariateAndBivariateRegressionAnalysisForPotentialGenotypeConfoundingWithRobustStandardErrors.txt"
  )

###### combine the univariate genotype regression summary statistics from the univariate genotype regressions with ordinary least squares (OLS) and robust standard errors (RSE)
     # the robust standard errors regressions adjust for heteroscedasticity in the regression

RegressionResultsFile_OLS_RSE <-
  # join the CausativeAssociatedVariants ordinary least squares regression summary statistics with the robust standard errors summary statistics
  # for the same variant/metabolite pairs
  CausativeAssociatedVariants %>%
  # select the useful columns only and add "_OLS" to the column names to indicate ordinary least squares (OLS)
  select(
    .data = .,
    existing_variant_VEP,
    CHR,
    SNP,
    "NMISS_OLS" = NMISS,
    "BETA_OLS" = BETA,
    "SE_OLS" = SE,
    "Pvalue_unadjust_plink_OLS" = Pvalue_unadjust_plink,
    "Bonferroni_multiplier_OLS" = Bonferroni_multiplier,
    "P_bonferroni_OLS" = P_bonferroni,
    Gene,
    metabolite,
    VariantID_GenotypeMatrixFormat,
    avsnp138,
    avsnp142,
    `Star Allele`:`Pubmed ID`
    ) %>%
    # group the ordinary least squares regression summary statistics by the metabolite
    group_by(
      .data = .,
      metabolite
      ) %>%
    # split into a list of dataframes based on the metabolite groupings
    group_split(.tbl = .) %>%
    # loop through each dataframe and join with the regression results file with robust standard errors
    # for each metabolite based on the variant ID
    lapply(
      X = .,
      FUN = function(currentDataFrame)
      {
        
        DataToReturn <-
            currentDataFrame %>%
              # join the currentDataFrame with the univariate robust standard errors genotype regressions based on the
              # variant IDs and the current metabolite
              left_join(
                x = .,
                y = UnivarateGenotypeRegressionsWithRobustStandardErrors %>%
                    # filter to univariate regressions of phenotype versus genotype only
                    filter(
                      .data = .,
                      Regression_type == "univariate"
                      ) %>%
                    # select useful columns and add "_RSE" to the column names to indicate robust standard errors (RSE)
                    select(
                      .data = .,
                      variant,
                      "N_RSE" = N,
                      "estimate_std.error_RSE" = estimate_std.error,
                      Heteroscedastiticy_genotype_p.value,
                      "F_stat_p_value_unadjusted_RSE" = F_stat_p_value,
                      "metabolite_RSE" = metabolite
                    ) %>%
                    # filter to the current metabolite in the currentDataFrame
                    filter(
                      .data = .,
                      metabolite_RSE == unique(x = currentDataFrame$metabolite)
                      ),
                by = c("VariantID_GenotypeMatrixFormat" = "variant")
                  ) %>%
              # remove the astericks from the unadjusted robust standard errors F-statistic P-value
              # and convert the column to a numeric
              mutate(
                .data = .,
                F_stat_p_value_unadjusted_RSE = F_stat_p_value_unadjusted_RSE %>%
                                                gsub(
                                                  pattern = "\\*",
                                                  replacement = "",
                                                  x = .
                                                  ) %>%
                                                 as.numeric(x = .)
                ) %>%
              # multiply the unadjusted robust standard errors F-statistic P-value by the Bonferroni multiplier
              # to adjust the P-value for multiple comparisons
              mutate(
                .data = .,
                F_stat_p_value_adjusted_RSE = F_stat_p_value_unadjusted_RSE*Bonferroni_multiplier_OLS
                ) %>%
              # round P-values and add astericks to P-values to indicate level of significance
              mutate(
                .data = .,
                F_stat_p_value_unadjusted_RSE = F_stat_p_value_unadjusted_RSE %>%
                                                signif(x = .,digits = 3) %>%
                                                AddAstericksToPvalues(columnVector = .),
                F_stat_p_value_adjusted_RSE = F_stat_p_value_adjusted_RSE %>%
                                              signif(x = .,digits = 3) %>%
                                              AddAstericksToPvalues(columnVector = .)
                )

        return(DataToReturn)
      }
        ) %>%
    do.call(
      what = "rbind",
      args = .
      ) %>%
  # reorganize the columns
  select(
    .data = .,
    existing_variant_VEP,
    Gene,
    CHR,
    SNP,
    VariantID_GenotypeMatrixFormat,
    "NonMiss_OLS" = NMISS_OLS,
    BETA_OLS,
    SE_OLS,
    Pvalue_unadjust_plink_OLS,
    Bonferroni_multiplier_OLS,
    P_bonferroni_OLS,
    Heteroscedastiticy_genotype_p.value,
    N_RSE,
    estimate_std.error_RSE,
    F_stat_p_value_unadjusted_RSE,
    F_stat_p_value_adjusted_RSE,
    metabolite,
    everything()
  )

# save the regression summary statistics with ordinary least squares and robust standard errors to the file system
RegressionResultsFile_OLS_RSE %>%
  write_tsv(
    x = .,
    "./CausativeVariantAssociations/RegressionResultsFile_OLS_and_RSE.tsv"
    )

############## multiple linear regressions for all demographic factors with each individual variant

# create a directory for saving multiple regression results
if(!dir.exists(paths = "./multiple_regression_results_with_genotype/"))
{
  dir.create(path = "./multiple_regression_results_with_genotype/")
}

# perform multiple regression with robust standard errors with all demographic factors along with genotype for significant variant/metabolite pairs only
# load the genotype matrices with metabolite data already previously computed from performing additive allele association
# list the names of the genotype matrix files and load each one individually
multipleRegressionResults <-
    list.files(path = "./genotype_matrices_with_metabolite_data/") %>%
    lapply(
      X = .,
      FUN = function(currentFileName)
      {

        GentoypeMatrixWithMetaboliteAndCovariateData <-
            read.delim(
              file = glue("./genotype_matrices_with_metabolite_data/{currentFileName}"),
              header = TRUE
            ) %>%
            # upon loading the genotype matrix join it with the covariate data by sequencing ID
            left_join(
              x = .,
              y = CovariateData,
              by = c("sampleID" = "VCFfileID")
                )

       # identify the variants that are in the genotype matrix
        VariantsInGenotypeMatrix <-
          GentoypeMatrixWithMetaboliteAndCovariateData %>%
          names(x = .) %>%
          grepl(pattern = "X",x = .) %>%
          which(x = .) %>%
          names(x = GentoypeMatrixWithMetaboliteAndCovariateData)[.]
        
        # identify the variants in the genotype matrix that have significant associations with a metabolite in primary analysis
        VariantsToLoopThrough <-
          (VariantsInGenotypeMatrix %in% CausativeAssociatedVariants$VariantID_GenotypeMatrixFormat) %>%
          which(x = .) %>%
          VariantsInGenotypeMatrix[.]

        # identify the current gene based on the file name
        currentGene <-
          currentFileName %>%
          str_extract(
            string = .,
            pattern = "(CYP3A4|DHCR7|CYP2R1|CYP27B1|SULT2A1|UGT1A4|GC|VDR|CASR|CUBN|CYP24A1|RXRG|RXRB|RXRA)"
              )

        print(x = currentGene)

        RegressionResults <-
         # perform a multiple linear regression with each significant variant in the genotype matrix with all appropriate metabolites
        VariantsToLoopThrough %>%
          mclapply(
            X = .,
            FUN = function(currentVariantToIncludeInRegression)
            {
              print(x = currentVariantToIncludeInRegression)
              
              # identify the metabolites to loop through based on the currentVariantToIncludeInRegression
              MetabolitesToLoopThrough <-
                  CausativeAssociatedVariants %>%
                    filter(
                      .data = .,
                      VariantID_GenotypeMatrixFormat == currentVariantToIncludeInRegression
                      ) %>%
                    pull(
                      .data = .,
                      metabolite
                      )

              MetabolitesToLoopThrough %>%
                lapply(
                  X = .,
                  FUN = function(currentMetabolite)
                  {
                      print(x = currentMetabolite)
                    
                    DataToUseInRegression <-
                      # select the columns of the metabolite,
                      # the current variant, and demographic factors: season, age, BMI, and gender
                      GentoypeMatrixWithMetaboliteAndCovariateData %>%
                      select(
                        .data = .,
                        "Phenotype" = all_of(currentMetabolite),
                        "Genotype" = all_of(currentVariantToIncludeInRegression),
                        "Age" = Age.on.Study.Date,
                        StudySeason,
                        BMI,
                        Sex
                      ) %>%
                      # ensure the genotype is continuous
                      mutate(
                        .data = .,
                        Genotype = Genotype %>% as.numeric(x = .)
                      ) %>%
                      # ensure the phenotype is a continuous numeric
                      mutate(
                        .data = .,
                        Phenotype = Phenotype %>% as.numeric(x = .)
                      ) %>%
                      # remove any missing values
                      na.omit(object = .)

                    # if there is more than one genotype value, perform the regression (e.g., 0 and 1)
                    # if there is only one genotype value (e.g., 0 only),
                      # return a NULL
                    # if there is only one study season as well, return a NULL,
                    # there must be more than one category to perform a regression with a categorical variable,
                    # if there is only one gender, return a NULL
                    if(
                      length(x = unique(x = as.numeric(x = DataToUseInRegression$Genotype)))>1 &
                      length(x = unique(x = as.factor(x = DataToUseInRegression$StudySeason)))>1 &
                      length(x = unique(x = as.factor(x = DataToUseInRegression$Sex)))>1
                      )
                    {

                      # obtain the robust standard error regression model
                      RegressionModel <-
                          DataToUseInRegression %>%
                          # define the multiple robust standard error regression model
                          lm_robust(
                            formula = Phenotype ~ Genotype + BMI + StudySeason + Age + Sex,
                            data = .
                          )

                      RegressionResultsToReturn <-
                        # tidy the model output 
                        TidyLinearRobustStandardErrorsRegressionResults(
                          DataForRegression = DataToUseInRegression,
                          RegressionModel = RegressionModel
                          ) %>%
                        # create a column label for the current gene
                        mutate(
                          .data = .,
                          Gene = rep_len(x = currentGene,length.out = nrow(x = .))
                        ) %>%
                        # create a column with a label for the metabolite
                        mutate(
                          .data = .,
                          metabolite = rep_len(x = currentMetabolite,length.out = nrow(x = .))
                        ) %>%
                        # create a column with a label for the variant
                        mutate(
                          .data = .,
                          variant = rep_len(x = currentVariantToIncludeInRegression,length.out = nrow(x = .))
                        ) 

                    } else {
                      RegressionResultsToReturn <- NULL
                    }
                    

                    return(RegressionResultsToReturn)

                  }
                    ) %>%
                do.call(
                  what = "rbind",
                  args = .
                  )


            },mc.cores = detectCores()
              ) %>%
          do.call(
            what = "rbind",
            args = .
            )

        return(RegressionResults)
      }
        ) %>%
      do.call(
        what = "rbind",
        args = .
        )

multipleRegressionResults %>%
write_tsv(
  x = .,
  file = "./multiple_regression_results_with_genotype/fully_adjusted_regression_with_genotype_results_individual_variants.txt"
    )

print(x = "multiple regression complete !")

############# identify common vitamin D binding protein haplotypes/diplotypes and include the haplotype/diplotypes in univariate metabolite regression #################

# create frequency plots of GC haplotypes and GC diplotypes
HaplotypePlot <-
  VitaminDbindingProteinHaplotypeMatrix %>%
  # select the two haplotype columns
  select(
    .data = .,
    GC_hap1,
    GC_hap2
    ) %>%
  # put the haplotype columns into one array
  unlist(x = .) %>%
  # remove the column names
  unname(obj = .) %>%
  # make the array a single column in a tibble
  tibble("haplotypes" = .) %>%
  # group by the haplotypes that are in the column
  group_by(
    .data = .,
    haplotypes
    ) %>%
  # count each of the haplotypes
  summarise(
    .data = .,
    haplotypeCount = n()
    ) %>%
  # ungroup the haplotype column
  ungroup(x = .) %>%
  # create a gene column just as a place holder for the bar in the graph
  mutate(
    .data = .,
    gene = "VDBP"
    ) %>%
  # create a bar chart filled by frequency of the haplotype
  ggplot(
    data = .,
    aes(
      x = gene,
      y = haplotypeCount,
      fill = haplotypes
      )
    ) +
  geom_col(position = "fill") +
  theme_classic() +
  ylab(label = "Haplotype frequency") +
  labs(fill = "GC haplotype") +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(family = "Arial",face = "bold",size = 11),
    plot.title = element_text(family = "Arial",face = "bold",size = 11,hjust = 0.5),
    axis.title.y = element_text(family = "Arial",face = "bold",size = 11,margin = margin(r = 15)),
    legend.text = element_text(family = "Arial",face = "bold",size = 11),
    legend.title = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.ticks.length.x = unit(x = 3,units = "pt"),
    plot.margin=unit(c(1,1,1,1),"cm")
      ) +
  scale_y_continuous(expand = c(0,0))

# create a haplotype frequency table as well using the data for the haplotype frequency plot
HaplotypeFrequencyTable <-
  VitaminDbindingProteinHaplotypeMatrix %>%
    # select the two haplotype columns
    select(
      .data = .,
      GC_hap1,
      GC_hap2
    ) %>%
    # put the haplotype columns into one array
    unlist(x = .) %>%
    # remove the column names
    unname(obj = .) %>%
    # make the array a single column in a tibble
    tibble("haplotypes" = .) %>%
    # group by the haplotypes that are in the column
    group_by(
      .data = .,
      haplotypes
    ) %>%
    # count each of the haplotypes
    summarise(
      .data = .,
      haplotypeCount = n()
    ) %>%
    # ungroup the table
    ungroup(x = .) %>%
    # compute the frequency of the haplotypes
    mutate(
      .data = .,
      haplotypeFrequency = haplotypeCount/sum(haplotypeCount,na.rm = TRUE)
      )

# create a GC diplotype plot
DiplotypePlot <-
  VitaminDbindingProteinHaplotypeMatrix %>%
  # select the diplotype column
  select(
    .data = .,
    GC_dip
    ) %>%
  # group by the diplotype
  group_by(
    .data = .,
    GC_dip
    ) %>%
  # count the number of each diplotype
  summarize(
    .data = .,
    DiplotypeCount = n()
    ) %>%
  # ungroup the dataframe
  ungroup(x = .) %>%
  # create a gene column as a place holder for the bar in the bar chart
  mutate(
    .data = .,
    gene = "VDBP"
    ) %>%
  # create a bar chart filled by frequency of the diplotype
  ggplot(
    data = .,
    aes(
      x = gene,
      y = DiplotypeCount,
      fill = GC_dip
      )
    ) +
  geom_col(position = "fill") +
  theme_classic() +
  ylab(label = "Diplotype frequency") +
  labs(fill = "GC diplotype") +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(family = "Arial",face = "bold",size = 11),
    plot.title = element_text(family = "Arial",face = "bold",size = 11,hjust = 0.5),
    axis.title.y = element_text(family = "Arial",face = "bold",size = 11,margin = margin(r = 15)),
    legend.text = element_text(family = "Arial",face = "bold",size = 11),
    legend.title = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.ticks.length.x = unit(x = 3,units = "pt"),
    plot.margin=unit(c(1,1,1,1),"cm")
      ) +
  scale_y_continuous(expand = c(0,0))

# create a diplotype frequency table as well using the data for the diplotype frequency plot
DiplotypeFrequencyTable <-
    VitaminDbindingProteinHaplotypeMatrix %>%
    # select the diplotype column
    select(
      .data = .,
      GC_dip
      ) %>%
    # group by the diplotype
    group_by(
      .data = .,
      GC_dip
      ) %>%
    # count the number of each diplotype
    summarize(
      .data = .,
      DiplotypeCount = n()
      ) %>%
    # ungroup the dataframe
    ungroup(x = .) %>%
    # compute diplotype frequencies
    mutate(
      .data = .,
      DiplotypeFrequency = DiplotypeCount/sum(DiplotypeCount,na.rm = TRUE)
      )

# create a directory for saving vitamin D binding protein haplotype analysis results
if(!dir.exists(paths = "./VDBP_haplotype_analysis/"))
{
  dir.create(path = "./VDBP_haplotype_analysis/")
}

# save the haplotype and diplotype frequency plots to the VDBP_haplotype_analysis directory
tiff(
  filename = "./VDBP_haplotype_analysis/HaplotypePlot.tiff",
  width = 7,
  height = 5,
  units = "in",
  compression = "none",
  res = 300
  )
HaplotypePlot
dev.off()

tiff(
  filename = "./VDBP_haplotype_analysis/DiplotypePlot.tiff",
  width = 7,
  height = 5,
  units = "in",
  compression = "none",
  res = 300
  )
DiplotypePlot
dev.off()

# save the haplotype and diplotype frequency tables to the VDBP_haplotype_analysis directory
mapply(
  FUN = function(currentTable,currentTableName)
  {
    currentTable %>%
      write_tsv(
        x = .,
        file = glue("./VDBP_haplotype_analysis/{currentTableName}.txt")
          )
  },
  list(HaplotypeFrequencyTable,DiplotypeFrequencyTable),
  c("HaplotypeFrequencyTable","DiplotypeFrequencyTable")
  )

# perform a robust standard errors regression to determine if there is any difference
# in metabolite concentration based on vitamin D binding protein diplotype groups
VDBPregressionTestResults <-
  Metabolites %>%
    lapply(
      X = .,
      FUN = function(currentMetabolite)
      {
        # select the metabolite, vitamin D binding protein diplotype, and covariate columns
        DataForRegression <-
            VitaminDbindingProteinDataForLinearRegression %>%
              select(
                .data = .,
                "Metabolite" = all_of(x = currentMetabolite),
                GC_dip,
                Age.on.Study.Date,
                BMI,
                Sex,
                StudySeason
              ) %>%
              na.omit(object = .) %>%
              mutate(
                .data = .,
                GC_dip = GC_dip %>% factor(x = .),
                Age.on.Study.Date = Age.on.Study.Date %>% as.numeric(x = .),
                BMI = BMI %>% as.numeric(x = .),
                Sex = Sex %>% factor(x = .),
                StudySeason = StudySeason %>% factor(x = .)
              )

        # perform the robust standard errors regression without covariates
        VDBPregressionResultsUnadjusted <-
          DataForRegression %>%
          lm_robust(
              Metabolite ~ GC_dip,
              data = .
              )

        # perform the robust standard errors regression with all covariates
        VDBPregressionResultsFullyAdjusted <-
          DataForRegression %>%
            lm_robust(
              Metabolite ~ GC_dip + Age.on.Study.Date + BMI + Sex + StudySeason,
              data = .
            )

        # test for heteroscedasticity of Metabolite concentration ~ GC diplotype
        VarianceHomogeneityTestResult <-
          DataForRegression %>%
          lm(
            formula = Metabolite ~ GC_dip,
            data = .
          ) %>%
          ols_test_f(model = .)

        ResultsToReturnUnadjusted <-
          TidyLinearRobustStandardErrorsRegressionResults(
            DataForRegression = DataForRegression,
            RegressionModel = VDBPregressionResultsUnadjusted
            ) %>%
          mutate(
            .data = .,
            # add the p-value from the variance homogeneity test
            "Heteroscedasticity_p.value" = VarianceHomogeneityTestResult$p %>%
                                           signif(x = .,digits = 3) %>%
                                           AddAstericksToPvalues(columnVector = .),
            # add the current metabolite
            "Metabolite" = currentMetabolite,
            # add a regression type label
            RegressionType = "unadjusted"
            )

        ResultsToReturnFullyAdjusted <-
          TidyLinearRobustStandardErrorsRegressionResults(
            DataForRegression = DataForRegression,
            RegressionModel = VDBPregressionResultsFullyAdjusted
              ) %>%
          mutate(
            .data = .,
            # add an NA p-value for the variance homogeneity test
            Heteroscedasticity_p.value = as.character(x = NA),
            # add the current metabolite
            "Metabolite" = currentMetabolite,
            RegressionType = "adjusted"
            )

        # create a dataframe of the unadjusted and adjusted regression results
        ResultsToReturn <-
          list(
            ResultsToReturnUnadjusted,
            ResultsToReturnFullyAdjusted
          ) %>%
          do.call(
            what = "rbind",
            args = .
            ) %>%
          data.frame()

        # temporarily suppress non-problematic warning messages for plotting
        options(warn = -1)

        VitaminDBindingProteinRegressionplot <-
            # create a plot of metabolite concentration stratified by diplotype
            VitaminDbindingProteinDataForLinearRegression %>%
              # group by the six possible diplotypes
              group_by(
                .data = .,
                GC_dip
                ) %>%
              # calculate the mean and standard deviation of the currentMetabolite and the diplotype count
              summarise(
                .data = .,
                "mean_conc" = mean(
                                  x = eval(expr = parse(text = currentMetabolite)),
                                  na.rm = TRUE
                                  ),
                "stdev" = sd(
                            x = eval(expr = parse(text = currentMetabolite)),
                            na.rm = TRUE
                            ),
                diplotypeCount = n()
                ) %>%
              # ungroup by the diplotype count
              ungroup(x = .) %>%
              # compute the standard error using the standard deviation and square root of the diplotype count
              mutate(
                .data = .,
                sterror = stdev/sqrt(x = diplotypeCount)
                ) %>%
              ggplot(
                data = .
              ) +
              geom_bar(
                mapping = aes(
                              x = GC_dip,
                              y = mean_conc
                              ),
                              stat = "identity",
                              fill="skyblue",
                              alpha = 0.7
                  ) +
              geom_errorbar(
                mapping = aes(
                  x=GC_dip,
                  ymin=mean_conc-sterror,
                  ymax=mean_conc+sterror
                  ),
                width=0.4,
                alpha=0.9,
                size=1.3,
                colour="orange"
                  ) +
              geom_point(
                data = VitaminDbindingProteinDataForLinearRegression %>%
                       select(
                         .data = .,
                         GC_dip,
                         !!currentMetabolite
                         ) %>%
                        na.omit(object = .),
                mapping = aes(
                              x = GC_dip,
                              y = eval(expr = parse(text = currentMetabolite))
                              ),
                position = "jitter"
                  )

          # identify the maximum metabolite value to decide how the y-axis limits should be formatted
          MaxMetaboliteValue <-
            VitaminDbindingProteinDataForLinearRegression %>%
            select(
              .data = .,
              GC_dip,
              !!currentMetabolite
            ) %>%
            na.omit(object = .) %>%
            pull(
              .data = .,
              eval(expr = parse(text = currentMetabolite))
              ) %>%
            max()

        VitaminDBindingProteinRegressionplot <-
          VitaminDBindingProteinRegressionplot +
          theme_classic() +
          xlab(label = "GC Diplotype") +
          ylab(label = "Concentration (ng/mL)") +
          ggtitle(label = paste0(currentMetabolite)) +
          theme(
            plot.title = element_text(family = "Arial",face = "bold",size = 11,hjust = 0.5),
            axis.title = element_text(
                                      face = "bold",
                                      size = 11,
                                      family = "Arial"
                                      ),
            axis.text.x = element_text(
                                      face = "bold",
                                      size = 11,
                                      family = "Arial",
                                      angle = 45,
                                      hjust = 1
                                      ),
            axis.text.y = element_text(
                                      face = "bold",
                                      size = 11,
                                      family = "Arial",
                                      margin = margin(t = 0,r = 20,b = 0,l = 0)
                                      )
            ) +
          # annotate the plot with the P-value from the regression without covariates performed above
          annotate(
            geom = "text",
            x = "GC1F/GC2",
            y = MaxMetaboliteValue/1.5,
            label = paste0("F-stat (P = ",ResultsToReturnUnadjusted$F_stat_p_value,")"),
            fontface = "bold",
            size = 4
          )

        tiff(
          filename = paste0("./VDBP_haplotype_analysis/","VitaminDBindingProteinRegressionplot_",currentMetabolite,".tiff"),
          width = 7,
          height = 5,
          units = "in",
          compression = "none",
          res = 300
          )

        print(x = VitaminDBindingProteinRegressionplot)
        dev.off()

        # turn global warning messages back on
        options(warn = 0)


        return(ResultsToReturn)
      }
        ) %>%
     do.call(
       what = "rbind",
       args = .
       )

VDBPregressionTestResults %>%
  write_tsv(
    x = .,
    file = "./VDBP_haplotype_analysis/VDBPregressionTestResults.tsv"
    )

print(x = "###############################################################################################")
print(x = "################                                                          #####################")
print(x = "################ vitamin D binding protein haplotype analysis is complete!#####################")
print(x = "################                                                          #####################")
print(x = "###############################################################################################")

############## perform a multiple linear regression with all independent causative variants (identified via clumping and multiple testing correction) 
             # that were significantly associated with metabolite concentration in robust standard errors linear regression
             # for each metabolite to estimate polygenic contribution to phenotypic variability

# create a dataset of independent/causative significantly associated variants from robust standard errors regression
# with seasonal and demographic covariates for all metabolites
PolygenicmultipleRegressionDataSet <-
  RegressionResultsFile_OLS_RSE %>%
  # filter to results that were significant (there is an astericks on the adjusted P-value) in robust standard errors regression
  filter(
    .data = .,
    grepl(
      pattern = "\\*",
      x = F_stat_p_value_adjusted_RSE
      )
    ) %>%
    # group the causative associated variants by metabolite
    group_by(
      .data = .,
      metabolite
      ) %>%
    # split into a list of dataframes by metabolite
    group_split(
      .tbl = .
      ) %>%
    # loop through each dataframe
    mclapply(
      X = .,
      FUN = function(currentDataSet)
      {
        # identify the genes with significant variants for the current metabolite
        GenesWithSignificantVariants <- currentDataSet$Gene %>% unique(x = .)
        # identify the current set of variants with the genotype matrix variant ID format
        currentSetOfVariants <- currentDataSet$VariantID_GenotypeMatrixFormat %>% unique(x = .)
        # identify the current metabolite
        currentMetabolite <- currentDataSet$metabolite %>% unique(x = .)
        
        # disable the non-problematic warning message that duplicate columns are introduced when joining dataframes below
        options(warn = -1)
        
          # obtain all significant variant genotypes for the current metabolite and metabolite values for the current metabolite
          DataToReturn <-
            GenesWithSignificantVariants %>%
            lapply(
              X = .,
              FUN = function(currentGene)
              {
                
                  # load the genotype matrix corresponding to the current gene
                  GenotypeMatrix <-
                    glue("./genotype_matrices_with_metabolite_data/{currentGene}_geno_matrix.raw") %>%
                    read.delim(
                      file = .,
                      header = TRUE
                      )

                  # identify the column names of the GenotypeMatrix
                  GenotypeMatrixColumnNames <- names(x = GenotypeMatrix)

                  # identify the column names in the genotype matrix that are in the currentSetOfVariants
                  VariantColumnsToSelect <-
                    (GenotypeMatrixColumnNames %in% currentSetOfVariants) %>%
                      which(x = .) %>%
                      GenotypeMatrixColumnNames[.]

                  DataToReturn <-
                    GenotypeMatrix %>%
                    # select the samples IDs and the current set of significant variants for the current metabolite
                    # and the current metabolite
                    select(
                      .data = .,
                      sampleID,
                      all_of(x = VariantColumnsToSelect),
                      all_of(x = currentMetabolite)
                    )

                  return(DataToReturn)
              }
                ) %>%
            # join all of the genotype data together by sample ID,
            # this will give a genotype matrix of all relevant variants for the current metabolite
            Reduce(
              f = JoinDataFrames,
              x = .
              )

          # turn warning messages back on
          options(warn = 0)
          
            # identify the phenotype column to select, the phenotype column to select will either have the label of the currentMetabolite
            # or the phenotype column to select will have the label of the currentMetabolite with a ".x" suffix in the case when there
            # are duplicated columns present after joining the genotype matrices above
            PhenotypeColumns <-
              names(x = DataToReturn) %>%
              grepl(
                pattern = currentMetabolite,
                x = .
                  ) %>%
              which(x = .) %>%
              names(x = DataToReturn)[.]

            # if any of the phenotype columns contains the suffix ".x",
            # select the phenotype column with the .x suffix
            if( any(PhenotypeColumns == paste0(currentMetabolite,".x")) )
            {
              PhenotypeColumnToSelect <-
                paste0(currentMetabolite,".x")
            # otherwise, select the currentMetabolite without the .x suffix
            } else {
              PhenotypeColumnToSelect <- currentMetabolite
            }

            DataToReturn <-
            DataToReturn %>%
              # select the sample ID, genotypes, and the current metabolite
              select(
                .data = .,
                sampleID,
                starts_with(match = "X"),
                !!currentMetabolite := !!PhenotypeColumnToSelect
                ) %>%
            # create a column with a label for the current metabolite
            mutate(
              .data = .,
              Metabolite = currentMetabolite
              ) %>%
            # join the genotypes and phenotype data with the relevant covariate data by sample ID
            left_join(
              x = .,
              y = CovariateData %>%
                  select(
                    .data = .,
                    VCFfileID,
                    Sex,
                    Age.on.Study.Date,
                    BMI,
                    StudySeason,
                    DayOfYear,
                    StudyMonth
                  ),
              by = c("sampleID" = "VCFfileID")
              )

        return(DataToReturn)
      },mc.cores = detectCores()
  ) 

# perform a multiple linear regression with robust standard errors of metabolite concentration versus significant variants for each metabolite 
# with sets of variants that were independently associated with metabolite concentration to obtain
# an estimate of the polygenic contribution (i.e., heritability) to metabolite variability for each individual metabolite

GeneticContributionRegressionSummary <-
  # loop through each dataset in the PolygenicmultipleRegressionDataSet list
  PolygenicmultipleRegressionDataSet %>%
    mclapply(
      X = .,
      FUN = function(currentDataSet)
      {
        
        # first specify the linear regression model formula with genotypes only
        RegressionModel <-
          currentDataSet %>%
          names(x = .) %>%
          grepl(pattern = "X",x = .) %>%
          which(x = .) %>%
          names(x = currentDataSet)[.] %>%
          paste(
            .,
            collapse = " + "
            ) %>%
          paste(
            unique(x = currentDataSet$Metabolite),
            .,
            sep = " ~ "
            )
        
        GeneticContributionModelFormula <-
          RegressionModel %>%
          as.formula(object = .)
        
        # remove any observations with missing values for the regressions
        DataForRegression <-
          currentDataSet %>%
          na.omit(object = .)
        
        # perform the multiple regression with robust standard errors for the genetic contribution only
        GeneticContributionModel <-
          DataForRegression %>%
          lm_robust(
            formula = GeneticContributionModelFormula,
            data = .
          )
        
        TidyGeneticContributionRegressionSummary <-
        # tidy the regression model output
          TidyLinearRobustStandardErrorsRegressionResults(
            DataForRegression = DataForRegression,
            RegressionModel = GeneticContributionModel
            ) %>%
          # add a column with the regression formula and the metabolite included
          mutate(
            .data = .,
            Regression_multiple_genotype = RegressionModel,
            Metabolite = unique(x = currentDataSet$Metabolite)
          ) 
        
      },mc.cores = detectCores()
        )

# identify each metabolite in the GeneticContributionRegressionSummary list
GeneticContributionRegressionSummaryNames <-
GeneticContributionRegressionSummary %>%
  mclapply(
    X = .,
    FUN = function(currentDataSet)
    {
      ValueToReturn <-
      currentDataSet$Metabolite %>%
        unique(x = .)

      return(ValueToReturn)
    },mc.cores = detectCores()
      )
# name each element of the GeneticContributionRegressionSummary list based on the metabolite
names(x = GeneticContributionRegressionSummary) <- GeneticContributionRegressionSummaryNames

# save the regression of the total genetic contribution to the file system for each metabolite
mapply(
  FUN = function(currentDataSet,currentDataSetName)
  {
    currentDataSet %>%
    write_tsv(
      x = .,
      file = glue("./multiple_regression_results_with_genotype/TotalGeneticContributionRegression_{currentDataSetName}.txt")
    )
  },
  GeneticContributionRegressionSummary,
  GeneticContributionRegressionSummaryNames,
  SIMPLIFY = FALSE
  )

# perform the multiple regression using the independent causative variants again for each metabolite
# with the additional seasonal and demographic factors to control for any potential confounding of the genetic associations
GeneticContributionWithCovariatesRegressionSummary <-
PolygenicmultipleRegressionDataSet %>%
  mclapply(
    X = .,
    FUN = function(currentDataSet)
    {
        # identify the current metabolite
        currentMetabolite <-
          currentDataSet$Metabolite %>%
          unique(x = .)
      
        # first specify the linear regression model formula with genotypes and seasonal/demographic factors
        RegressionFormula <-
          currentDataSet %>%
          # remove the sample ID, phenotype columns, and other unnecessary columns
          select(
            .data = .,
            -sampleID,
            -Metabolite,
            -!!currentMetabolite,
            -DayOfYear,
            -StudyMonth
            ) %>%
          # obtain an array of genotype and seasonal/demographic variables
          names(x = .) %>%
          paste(
            .,
            collapse = " + "
          ) %>%
          paste(
            currentMetabolite,
            .,
            sep = " ~ "
          )
          
        GeneticContributionWithCovariatesModelFormula <-
          RegressionFormula %>%
          as.formula(object = .)
        
        # remove any observations with missing values for the regression dataset
        DataForRegression <-
          currentDataSet %>%
            na.omit(object = .) 
          
        # perform the multiple linear regression with robust standard errors
        GeneticContributionWithCovariatesModel <-
          DataForRegression %>%
          lm_robust(
            formula = GeneticContributionWithCovariatesModelFormula,
            data = .
            )
        
        TidyGeneticContributionWithCovariatesRegressionSummary <-
          # tidy the summary 
          TidyLinearRobustStandardErrorsRegressionResults(
            DataForRegression = DataForRegression,
            RegressionModel = GeneticContributionWithCovariatesModel
            ) %>%
          # add columns with the regression formula and metabolite
          mutate(
            .data = .,
            Metabolite = unique(x = currentDataSet$Metabolite),
            Regression = RegressionFormula
          ) %>%
          # add a column with rs numbers and the gene corresponding to each variant ID
          left_join(
            x = .,
            y = CausativeAssociatedVariants %>%
                filter(
                  .data = .,
                  metabolite == unique(x = currentDataSet$Metabolite)
                  ) %>%
                select(
                  .data = .,
                  existing_variant_VEP,
                  VariantID_GenotypeMatrixFormat
                ),
            by = c("term" = "VariantID_GenotypeMatrixFormat")
          ) 
        
        return(TidyGeneticContributionWithCovariatesRegressionSummary)
        
    },mc.cores = detectCores()
      )

# identify each metabolite in the GeneticContributionWithCovariatesRegressionSummary list
GeneticContributionWithCovariatesRegressionSummaryNames <-
  GeneticContributionWithCovariatesRegressionSummary %>%
  lapply(
    X = .,
    FUN = function(currentDataSet)
    {
      ValueToReturn <-
      currentDataSet$Metabolite %>%
        unique(x = .)

      return(ValueToReturn)
    }
      )
# name each element of the GeneticContributionWithCovariatesRegressionSummary list based on the metabolite
names(x = GeneticContributionWithCovariatesRegressionSummary) <- GeneticContributionWithCovariatesRegressionSummaryNames

# save the regression of the total genetic contribution with covariates added to the file system for each metabolite
mapply(
  FUN = function(currentDataSet,currentDataSetName)
  {
    currentDataSet %>%
    write_tsv(
      x = .,
      file = glue("./multiple_regression_results_with_genotype/GeneticContributionWithDemographicAndSeasonalCovariates_{currentDataSetName}.txt")
    )
  },
  GeneticContributionWithCovariatesRegressionSummary,
  GeneticContributionWithCovariatesRegressionSummaryNames,
  SIMPLIFY = FALSE
  )
