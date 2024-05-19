# This script performs analysis of binary vitamin D status phenotypes [sufficiency and deficiency], the analysis includes:
  # (1) Univariate and multiple logistic regressions of demographic/seasonal covariates with binary vitamin D status phenotypes
  # (2) multiple logistic regressions with individual genetic variants included with adjustment for all seasonal and demographic covariates
  # (3) Logistic regression of binary phenotypes versus all significant causative variants with and without seasonal/demographic covariates included
  # (4) Logistic regression of binary phenotypes versus common vitamin D binding protein diplotypes

# load necessary R packages
source("../scripts/load_R_packages.R")
# load necessary functions
source("../scripts/genotype_phenotype_association/genetic_association_analysis_functions.R")

# load the binary vitamin D status phenotype data
BinaryPhenotypeData <-
  read.delim(
             file = "./BinaryPhenotypeData_prepared_for_logistic_regressions.txt",
             header = TRUE
             )

# load the covariate file that was prepared in the linear_regressions_continuous_phenotypes.R script
CovariateData <-
  read.delim(
    file = "./CovariateData_prepared_for_regressions.txt",
    header = TRUE
    )

# load the binary phenotype data with covariate data added
BinaryPhenotypeAndCovariateData <-
  read.delim(
    file = "./BinaryPhenotypeAndCovariateData_prepared_for_logistic_regressions.txt",
    header = TRUE
    )
  
# load the vitamin D binding protein data with logistic regression phenotypes added
VitaminDbindingProteinDataForLogisticRegression <-
  read.delim(
    file = "./VitaminDbindingProteinDataForLogisticRegression.txt",
    header = TRUE
    )

# load the dbSNP rsID lookup table
rsIDlookupTable <-
  read.delim(
    file = "./PLINKvariantID_dbSNP_rsID_lookupTable_commonVariants.txt",
    header = TRUE
      )

# load the causitive independent variants from logistic regression that were identified via P-value based linkage disequilibrium pruning
ClumpData <-
  read.delim(
    file = "./ClumpDataLogistic_prepared_for_logistic_regressions.txt",
    header = TRUE
    )

# create a file for saving all vitamin D status logistic regressions if it does not exist
if(!dir.exists(paths = "./logistic_regressions_with_covariates/"))
{
  dir.create(path = "./logistic_regressions_with_covariates/")
}

########## perform logistic regressions of binary phenotypes versus each individual demographic and seasonal covariate
  
UnivariateDemographicAndSeasonalLogisticRegressions <-
  # loop through the binary phenotypes
  c(
    "VitDsufficiency",
    "VitDdeficiency"
    ) %>%
    lapply(
      X = .,
      FUN = function(currentPhenotype)
      {
         ResultsToReturn <-
            # loop through the covariates
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
                        # select columns from the BinaryPhenotypeAndCovariateData 
                        # for the currentPhenotype and currentCovariate
                        DataForRegression <-
                          BinaryPhenotypeAndCovariateData %>%
                          select(
                            .data = .,
                            "Phenotype" = all_of(x = currentPhenotype),
                            "Covariate" = all_of(x = currentCovariate)
                            ) %>%
                          # remove missing observations 
                          na.omit(object = .) 
                        
                        # set a random number seed in case one is required
                        # for reproducibility of the logistic regression
                        set.seed(seed = 777)
                        
                        # define the logistic regression model
                        RegressionModel <-
                          DataForRegression %>%
                          glm(
                            formula = Phenotype ~ Covariate,
                            family = "binomial",
                            data = .
                            )
                        
                        # If the current covariate is a categorical covariate,
                        # test for the overall significance of the covariate with the wald.test function from the aod package
                        # For more information on this, see: https://stats.idre.ucla.edu/mplus/dae/logit-regression/
                        if(
                          (currentCovariate == "Sex") |
                          (currentCovariate == "StudySeason")
                          )
                        {
                            # identify the column indices of the variance-covariance matrix
                            # columns that correspond to the categorical covariate
                            IndicesOfTermsToSelectFromVarianceCovarianceMatrix <-
                              vcov(object = RegressionModel) %>%
                              colnames(x = .) %>%
                              grepl(pattern = "Covariate",x = .) %>%
                              which(x = .)

                            # Compute the Wald test with the model coefficients from the logistic regression model,
                            # the variance-covariance matrix of the regression model,
                            # and the categorical covariate terms from the variance-covariance matrix
                            WaldTestForCategoricalCovariateSignificanceResult <-
                              wald.test(
                                b = coef(object = RegressionModel),
                                Sigma = vcov(object = RegressionModel),
                                Terms = IndicesOfTermsToSelectFromVarianceCovarianceMatrix
                                  )

                            # extract the P-value from the overall significance test for the categorical covariate
                            WaldTestForCategoricalCovariateSignificanceResultPvalue <-
                              WaldTestForCategoricalCovariateSignificanceResult$result$chi2["P"] %>%
                              unname(obj = .)

                        }
                        
                        # tidy the regression summary and return the results
                        ResultsToReturn <-
                          TidyLogisticRegressionResults(
                            DataForRegression = DataForRegression,
                            RegressionModel = RegressionModel
                            ) %>%
                           mutate(
                            .data = .,
                            # add columns with the phenotype and covariate
                            Phenotype = currentPhenotype,
                            Covariate = currentCovariate,
                            # add a column for the confounder ("none" in this case)
                            Confounder = "none",
                            # add a column for the overall significance of the categorical confounder ("not_available" in this case since there isn't one)
                            CategoricalConfounderSignificance = "not_available",
                            # add a column for the regression
                            Regression = paste0(currentPhenotype," ~ ",currentCovariate),
                            # add a column for the regression type
                            RegressionType = "univariate"
                            )
                        
                        # If the current covariate is a categorical covariate,
                        # add the P-value from the wald test for overall significance for the categorical covariate,
                        # otherwise, label the covariate as non-categorical
                        if(
                          (currentCovariate == "Sex") |
                          (currentCovariate == "StudySeason")
                        )
                        {

                            ResultsToReturn <-
                                ResultsToReturn %>%
                                mutate(
                                  .data = .,
                                  CategoricalCovariateSignificance = WaldTestForCategoricalCovariateSignificanceResultPvalue %>%
                                                                     signif(x = .,digits = 3) %>%
                                                                     AddAstericksToPvalues(columnVector = .)
                                       )

                        } else {

                            ResultsToReturn <-
                              ResultsToReturn %>%
                              mutate(
                                .data = .,
                                CategoricalCovariateSignificance = "covariate_is_not_categorical"
                                     )

                        }

                        return(ResultsToReturn)
                        
                  }
                    ) %>%
                  do.call(
                    what = "rbind",
                    args = .
                    )
         
            return(ResultsToReturn)
      }
        ) %>%
    do.call(
      what = "rbind",
      args = .
      )

########## perform logistic regressions of binary phenotypes versus each possible pair of 
########## demographic and seasonal covariates to identify any possible confounding relationships


BivariateDemographicAndSeasonalLogisticRegressions <-
    # loop through the binary phenotypes
   c(
     "VitDsufficiency",
     "VitDdeficiency"
     ) %>%
     lapply(
       X = .,
       FUN = function(currentPhenotype)
       {
         ResultsToReturn <-
            # loop through the seasonal and demographic covariates
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
                      ResultsToReturn <-
                        # loop through the potential confounding variables (the same seasonal and demographic factors)
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
                                
                                  # if the currentCovariate is not the same as the currentPotentialConfounder,
                                  # perform a regression of phenotype ~ currentCovariate + currentPotentialConfounder
                                  if(currentCovariate!=currentPotentialConfounder)
                                  {
                                      
                                      DataForRegression <-
                                        # select the columns for the phenotype, currentCovariate, and currentPotentialConfounder
                                        BinaryPhenotypeAndCovariateData %>%
                                        select(
                                          .data = .,
                                          "phenotype" = all_of(x = currentPhenotype),
                                          "covariate" = all_of(x = currentCovariate),
                                          "confounder" = all_of(x = currentPotentialConfounder)
                                          ) %>%
                                          # remove missing values
                                         na.omit(object = .)
                                      
                                        # set a random number seed in case one is required
                                        # for reproducibility of the logistic regression
                                        set.seed(seed = 777)
                                      
                                        # define the regression model
                                        RegressionModel <-
                                          DataForRegression %>%
                                          glm(
                                            formula = phenotype ~ covariate + confounder,
                                            family = "binomial",
                                            data = .
                                            )
                                        
                                        # If the current covariate is a categorical covariate,
                                        # test for the overall significance of the covariate with the wald.test function from the aod package
                                            # For more information on this, see: https://stats.idre.ucla.edu/mplus/dae/logit-regression/
                                        if(
                                          (currentCovariate == "Sex") |
                                          (currentCovariate == "StudySeason")
                                          )
                                        {
                                            # identify the column indices of the variance-covariance matrix
                                            # columns that correspond to the categorical covariate
                                            IndicesOfTermsToSelectFromVarianceCovarianceMatrix <-
                                              vcov(object = RegressionModel) %>%
                                              colnames(x = .) %>%
                                              grepl(pattern = "covariate",x = .) %>%
                                              which(x = .)
                                          
                                            # Compute the Wald test with the model coefficients from the logistic regression model,
                                            # the variance-covariance matrix of the regression model,
                                            # and the categorical covariate terms from the variance-covariance matrix
                                            WaldTestForCategoricalCovariateSignificanceResult <-
                                              wald.test(
                                                b = coef(object = RegressionModel),
                                                Sigma = vcov(object = RegressionModel),
                                                Terms = IndicesOfTermsToSelectFromVarianceCovarianceMatrix
                                                  )

                                            # extract the P-value from the overall significance test for the categorical covariate
                                            WaldTestForCategoricalCovariateSignificanceResultPvalue <-
                                              WaldTestForCategoricalCovariateSignificanceResult$result$chi2["P"] %>%
                                              unname(obj = .)

                                        }

                                        # If the current confounder is a categorical variable,
                                        # test for the overall significance of the confounder with the wald.test function from the aod package
                                            # For more information on this, see: https://stats.idre.ucla.edu/mplus/dae/logit-regression/
                                        if(
                                          (currentPotentialConfounder == "Sex") |
                                          (currentPotentialConfounder == "StudySeason")
                                          )
                                        {
                                            # identify the column indices of the variance-covariance matrix
                                            # columns that correspond to the categorical confounder
                                            IndicesOfTermsToSelectFromVarianceCovarianceMatrix <-
                                              vcov(object = RegressionModel) %>%
                                              colnames(x = .) %>%
                                              grepl(pattern = "confounder",x = .) %>%
                                              which(x = .)
                                          
                                            # Compute the Wald test with the model coefficients from the logistic regression model,
                                            # the variance-covariance matrix of the regression model,
                                            # and the categorical confounder terms from the variance-covariance matrix
                                            WaldTestForCategoricalConfounderSignificanceResult <-
                                              wald.test(
                                                b = coef(object = RegressionModel),
                                                Sigma = vcov(object = RegressionModel),
                                                Terms = IndicesOfTermsToSelectFromVarianceCovarianceMatrix
                                                  )

                                            # extract the P-value from the overall significance test for the categorical confounder
                                            WaldTestForCategoricalConfounderSignificanceResultPvalue <-
                                              WaldTestForCategoricalConfounderSignificanceResult$result$chi2["P"] %>%
                                              unname(obj = .)

                                        }
                                        
                                        ResultsToReturn <-
                                          # tidy the regression summary and return the results 
                                          TidyLogisticRegressionResults(
                                            DataForRegression = DataForRegression,
                                            RegressionModel = RegressionModel
                                            ) %>%
                                            # add columns for the phenotype, covariate, and confounder
                                           mutate(
                                             .data = .,
                                             Phenotype = currentPhenotype,
                                             Covariate = currentCovariate,
                                             Confounder = currentPotentialConfounder,
                                             # add a column for the regression and regression type
                                             Regression = paste0(currentPhenotype," ~ ",currentCovariate," + ",currentPotentialConfounder),
                                             RegressionType = "bivariate"
                                             ) 
                                        
                                        # If the current covariate is categorical,
                                        # add the Wald-test results for overall significance for the categorical term to the results to return
                                        if(
                                          (currentCovariate == "Sex") |
                                          (currentCovariate == "StudySeason") 
                                        ) 
                                        {
                                          ResultsToReturn <-
                                            ResultsToReturn %>%
                                            mutate(
                                              .data = .,
                                              CategoricalCovariateSignificance = WaldTestForCategoricalCovariateSignificanceResultPvalue %>%
                                                                                 signif(x = .,digits = 3) %>%
                                                                                AddAstericksToPvalues(columnVector = .)
                                              )
                                        } else {
                                          ResultsToReturn <-
                                            ResultsToReturn %>%
                                            mutate(
                                              .data = .,
                                              CategoricalCovariateSignificance = "covariate_is_not_categorical"
                                                   )
                                        }
                                        
                                      # If the current potential confounder is categorical,
                                      # add the Wald-test results for overall significance for the categorical term to the results to return
                                      if(
                                        (currentPotentialConfounder == "Sex") |
                                        (currentPotentialConfounder == "StudySeason") 
                                      ) 
                                      {
                                        
                                        ResultsToReturn <-
                                          ResultsToReturn %>%
                                          mutate(
                                            .data = .,
                                            CategoricalConfounderSignificance = WaldTestForCategoricalConfounderSignificanceResultPvalue %>%
                                                                                signif(x = .,digits = 3) %>%
                                                                                AddAstericksToPvalues(columnVector = .)
                                            )
                                        
                                      } else {
                                        
                                        ResultsToReturn <-
                                          ResultsToReturn %>%
                                          mutate(
                                            .data = .,
                                            CategoricalConfounderSignificance = "confounder_is_not_categorical"
                                                 )
                                        
                                      }
                                        
                                    } else {
                                    
                                        ResultsToReturn <- NULL
                                    } 
                            
                                    return(ResultsToReturn)
                                  
                            
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

##### bind the univariate and bivariate DemographicAndSeasonalLogisticRegressions together to easily compare odds ratios 
      # for covariates in the presence and absence of a confounder. 
  UnivariateAndBivariateDemographicAndSeasonalLogisticRegressions <-
      list(
        UnivariateDemographicAndSeasonalLogisticRegressions,
        BivariateDemographicAndSeasonalLogisticRegressions
      ) %>%
      # bind the dataframes by row
      do.call(
        what = "rbind",
        args = .
        ) %>%
      data.frame() %>%
      # group by the phenotype, then by the covariate
      group_by(
        .data = .,
        Phenotype,
        Covariate
        ) %>%
      # split into a list of dataframes based on the phenotype and covariate
      group_split(.tbl = .) %>%
      # bind the list of dataframes back together to keep the list sorting
      do.call(
        what = "rbind",
        args = .
        ) %>%
      data.frame() 
  
##### identify the causative associated variants by filtering the clump data to variants with a bonferroni adjusted P-value of 0.05 or less
  CausativeAssociatedVariants <-
      ClumpData %>%
      # group the results by phenotype to create a column with bonferroni adjusted P-values
      group_by(
        .data = .,
        phenotype
        ) %>%
      # split into a list of dataframes by phenotype
      group_split(.tbl = .) %>%
      # loop through each dataframe and compute bonferroni adjusted P-values
      lapply(
        X = .,
        FUN = function(currentPhenotypeAssociations)
        {
          DataToReturn <-
            currentPhenotypeAssociations %>%
            # the bonferroni adjusted P-values are computed by multiplying the unadjusted 
            # P-value by the number of independent causative variants that were tested
            mutate(
              .data = .,
              P_bonferroni = nrow(x = currentPhenotypeAssociations)*P
            ) %>%
            # add a column with the number of snps used for the bonferroni adjustment
            mutate(
              .data = .,
              Bonferroni_multiplier = nrow(x = currentPhenotypeAssociations)
              )
  
          return(DataToReturn)
        }
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
        phenotype,
        P_bonferroni,
        Bonferroni_multiplier
      ) %>%
      # join the clump results with the rsIDlookupTable
      left_join(
        x = .,
        y = rsIDlookupTable,
        by = "SNP"
          ) %>%
      # group by gene and phenotype
      group_by(
        .data = .,
        Gene,
        phenotype
        ) %>%
      # split into a list of dataframes by gene and phenotype
      group_split(.tbl = .) %>%
      lapply(
        X = .,
        FUN = function(currentDataFrame)
        {
          
          # if the gene label in the current dataframe is UGT1A4, change it to UGT1A1
          # so the join in the code below is completed correctly
          currentDataFrame <-
            currentDataFrame %>%
            mutate(
              .data = .,
              Gene = if_else(
                             condition = Gene == "UGT1A4",
                             true = "UGT1A1",
                             false = Gene
                               )
              )
          
          DataToReturn <-
              currentDataFrame %>%
              # join the causative associated variants with the logistic regression data from
              # the original logistic association with PLINK to obtain the odds ratio, statistic, unadjusted p-value, and sample sizes
              left_join(
                x = .,
                y = "./logistic_regression_results_file_and_annotation_common_variants.txt" %>%
                    # load the logistic regression results file and annotation data from the original logistic association
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
                    # select the SNP, odds-ratio, std. errror, confidence interval, stat, sample size, phenotype, gene, and unadjusted P-value from the PLINK association
                    select(
                      .data = .,
                      SNP,
                      OR,
                      SE,
                      L95,
                      U95,
                      STAT,
                      "Pvalue_unadjust_plink" = Pvalue,
                      NMISS,
                      Phenotype,
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
                    # filter to the current gene and phenotype that is in the current data frame
                    filter(
                      .data = .,
                      (Phenotype == unique(x = currentDataFrame$phenotype)) #&
                        #(Gene == unique(x = currentDataFrame$Gene))
                      ) %>%
                     # deselect the phenotype column and gene column so there aren't duplicate columns in the final table
                    select(
                      .data = .,
                      -Phenotype,
                      -Gene
                    ),
                by = "SNP"
                )

          return(DataToReturn)
        }
          ) %>%
      # bind the results back together by row
      do.call(
        what = "rbind",
        args = .
        ) %>%
      # round the Bonferroni adjusted P-values and the unadjusted P-values and add astericks with levels of significance
      mutate(
        .data = .,
        P_bonferroni = P_bonferroni %>% 
                       signif(x = .,digits = 3) %>% 
                       AddAstericksToPvalues(columnVector = .),
        Pvalue_unadjust_plink = Pvalue_unadjust_plink %>% 
                                signif(x = .,digits = 3) %>% 
                                AddAstericksToPvalues(columnVector = .),
        # paste the odds ratio and confidence interval columns together
        OR_CI = paste0(OR," (",L95,"-",U95,")")
        )
  
##### perform univariate logistic regressions of phenotype versus each individual, significant, causative significant variant
  
  UnivariateGenotypeLogisticRegressions <-
    # loop through each of the genotype matrices with binary phenotypes added
    list.files(path = "./genotype_matrices_with_binary_phenotypes/")  %>%
    lapply(
      X = .,
      FUN = function(currentFileName)
      { 
          GenotypeMatrixWithPhenotypeData <-
            # load the genotype matrix with binary phenotype data
            read.delim(
              file = glue("./genotype_matrices_with_binary_phenotypes/{currentFileName}"),
              header = TRUE
              )

          # identify the variants that are in the genotype matrix based on the "X" prefix
          VariantsInGenotypeMatrix <-
            GenotypeMatrixWithPhenotypeData %>%
            names(x = .) %>%
            grepl(pattern = "X",x = .) %>%
            which(x = .) %>%
            names(x = GenotypeMatrixWithPhenotypeData)[.]

          # identify the variants in the genotype matrix that have significant associations with a binary phenotype in primary analysis
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
          
          # perform a univariate logistic regression of phenotype ~ genotype if there are any significant variants
          # with each significant variant in the genotype matrix with all appropriate phenotypes
          if(length(x = VariantsToLoopThrough)>0)
          {
            
            ResultsToReturn <-
                VariantsToLoopThrough %>%
                   lapply(
                     X = .,
                     FUN = function(currentVariantToIncludeInRegression)
                     {
                       
                          # identify the phenotypes to loop through based on the currentVariantToIncludeInRegression
                          PhenotypesToLoopThrough <-
                            CausativeAssociatedVariants %>%
                              filter(
                                .data = .,
                                VariantID_GenotypeMatrixFormat == currentVariantToIncludeInRegression
                                ) %>%
                              pull(
                                .data = .,
                                phenotype
                                )
                          
                        ResultsToReturn <-
                            PhenotypesToLoopThrough %>%
                              lapply(
                                X = .,
                                FUN = function(currentPhenotype)
                                {
                                      # select the columns of the current phenotype and current variant genotype
                                      DataForRegression <-
                                        GenotypeMatrixWithPhenotypeData %>%
                                           select(
                                             .data = .,
                                             "Phenotype" = all_of(x = currentPhenotype),
                                             "Genotype" = all_of(x = currentVariantToIncludeInRegression)
                                             ) %>%
                                            # ensure genotype is continuous and phenotype is binary
                                            mutate(
                                              .data = .,
                                              Genotype = Genotype %>% as.numeric(x = .),
                                              Phenotype = Phenotype %>% factor(x = .)
                                              ) %>%
                                            # remove missing values
                                            na.omit(object = .)
                                      
                                      
                                      # if there is more than one genotype value (e.g., 0 and 1), perform the logistic regression
                                      if(
                                         length(x = unique(x = as.numeric(x = DataForRegression$Genotype)))>1
                                        )
                                      {
                                         
                                          # set a random number seed for logistic regression
                                          set.seed(777)
                                        
                                          # define the regression model
                                          RegressionModel <-
                                            DataForRegression %>%
                                            glm(
                                              formula = Phenotype ~ Genotype,
                                              family = "binomial",
                                              data = .
                                              )
                                          
                                          ResultsToReturn <-
                                            # tidy the logistic regression model results
                                            TidyLogisticRegressionResults(
                                              DataForRegression = DataForRegression,
                                              RegressionModel = RegressionModel
                                              ) %>%
                                            # create a column for the current gene
                                            mutate(
                                              .data = .,
                                              Gene = currentGene,
                                              # the current phenotype
                                              Phenotype = currentPhenotype,
                                              # the current variant
                                              Variant = currentVariantToIncludeInRegression,
                                              # the current covariate (none in this case)
                                              Covariate = "none",
                                              # the overall significance of the categorical covariate in logistic regression ("not_available" because there isn't one)
                                              CategoricalCovariateSignificance = "not_available",
                                              # a label for the regression
                                              Regression = paste0(currentPhenotype," ~ ",currentVariantToIncludeInRegression),
                                              Regression_type = "univariate"
                                              )
                                          
                                      } else {
                                        
                                         ResultsToReturn <- NULL
                                      }
                                      
                                      return(ResultsToReturn)
                                }
                                  ) %>%
                                do.call(
                                  what = "rbind",
                                  args = .
                                  )
                        
                        
                                return(ResultsToReturn)
                     }
                       ) %>%
                      do.call(
                        what = "rbind",
                        args = .
                        )
              
            
          } else {
            
            ResultsToReturn <- NULL
          }
          
          return(ResultsToReturn)
      }
        ) %>%
     do.call(
       what = "rbind",
       args = .
       )
  
 ###### perform logistic regression of phenotype versus genotype for significant genotype/phenotype pairs plus each
        # individual demographic and seasonal covariate to see if the relationship between phenotype and genotype is being confounded
    BivariateGenotypeLogisticRegressions <-
      # loop through each of the genotype matrices with binary phenotypes added
      list.files(path = "./genotype_matrices_with_binary_phenotypes/")  %>%
        lapply(
          X = .,
          FUN = function(currentFileName)
          {
              GenotypeMatrixWithPhenotypeAndCovariateData <-
                # load the genotype matrix with binary phenotype data
                read.delim(
                  file = glue("./genotype_matrices_with_binary_phenotypes/{currentFileName}"),
                  header = TRUE
                  ) %>%
                # join the genotype matrix with covariate data by sequencing sample ID
                left_join(
                  x = .,
                  y = CovariateData,
                  by = c("sampleID" = "VCFfileID")
                  )

              # identify the variants that are in the genotype matrix based on the "X" prefix
              VariantsInGenotypeMatrix <-
                GenotypeMatrixWithPhenotypeAndCovariateData %>%
                names(x = .) %>%
                grepl(pattern = "X",x = .) %>%
                which(x = .) %>%
                names(x = GenotypeMatrixWithPhenotypeAndCovariateData)[.]

              # identify the variants in the genotype matrix that have significant associations with a binary phenotype in primary analysis
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

              # perform a logistic regression of phenotype ~ genotype + each individual covariate
              # if there are any significant variants in the genotype matrix
              if(length(x = VariantsToLoopThrough)>0)
              {
                 ResultsToReturn <-
                    # loop through each variant
                    VariantsToLoopThrough %>%
                      lapply(
                        X = .,
                        FUN = function(currentVariantToIncludeInRegression)
                        {
                            # identify the phenotypes to loop through based on the currentVariantToIncludeInRegression
                            PhenotypesToLoopThrough <-
                              CausativeAssociatedVariants %>%
                                filter(
                                  .data = .,
                                  VariantID_GenotypeMatrixFormat == currentVariantToIncludeInRegression
                                  ) %>%
                                pull(
                                  .data = .,
                                  phenotype
                                  )

                            ResultsToReturn <-
                                # loop through each possible phenotype
                                PhenotypesToLoopThrough %>%
                                  lapply(
                                    X = .,
                                    FUN = function(currentPhenotype)
                                    {
                                        ResultsToReturn <-
                                            # loop through each possible covariate
                                            c("Age.on.Study.Date","StudySeason","BMI","Sex") %>%
                                            sort(x = .) %>%
                                            lapply(
                                              X = .,
                                              FUN = function(currentCovariate)
                                              {
                                                      
                                                    # select columns for the phenotype, genotype, and covariate
                                                    DataForRegression <-
                                                      GenotypeMatrixWithPhenotypeAndCovariateData %>%
                                                      select(
                                                        .data = .,
                                                        "Phenotype" = all_of(x = currentPhenotype),
                                                        "Genotype" = all_of(x = currentVariantToIncludeInRegression),
                                                        "Covariate" = all_of(x = currentCovariate)
                                                        ) %>%
                                                       # make sure genotype is continuous and phenotype is categorical
                                                      mutate(
                                                        .data = .,
                                                        Genotype = Genotype %>% as.numeric(x = .),
                                                        Phenotype = Phenotype %>% factor(x = .)
                                                        ) %>%
                                                        # remove missing values
                                                      na.omit(object = .)
                                                    
                                                    # if the current covariate is StudySeason or Sex, make sure it is categorical
                                                    if(currentCovariate == "StudySeason" | currentCovariate == "Sex")
                                                    {
                                                      DataForRegression <-
                                                        DataForRegression %>%
                                                        mutate(
                                                          .data = .,
                                                          Covariate = Covariate %>%
                                                                      as.factor(x = .)
                                                          )
                                                    }

                                                    # if the current covariate is Age.on.Study.Date or BMI, make sure it is numeric
                                                    if(currentCovariate == "Age.on.Study.Date" | currentCovariate == "BMI")
                                                    {
                                                      DataForRegression <-
                                                        DataForRegression %>%
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
                                                          (length(x = unique(x = as.numeric(x = DataForRegression$Genotype)))>1) &
                                                          (is.numeric(x = DataForRegression$Covariate))
                                                        ) |
                                                        (
                                                          (length(x = unique(x = as.numeric(x = DataForRegression$Genotype)))>1) &
                                                          (length(x = unique(x = as.factor(x = DataForRegression$Covariate)))>1)
                                                        )
                                                      )
                                                    {
                                                            # set a random number seed for reproducibility 
                                                            set.seed(777)
                                                      
                                                            # define the regression model
                                                            RegressionModel <-
                                                              DataForRegression %>%
                                                              glm(
                                                                formula = Phenotype ~ Genotype + Covariate,
                                                                family = "binomial",
                                                                data = .
                                                                ) 
                                                            
                                                            # If the current covariate is a categorical covariate,
                                                            # test for the overall significance of the covariate with the wald.test function from the aod package
                                                                # For more information on this, see: https://stats.idre.ucla.edu/mplus/dae/logit-regression/
                                                            if(
                                                              (currentCovariate == "Sex") |
                                                              (currentCovariate == "StudySeason")
                                                              )
                                                            {
                                                              
                                                                # identify the column indices of the variance-covariance matrix
                                                                # columns that correspond to the categorical covariate
                                                                IndicesOfTermsToSelectFromVarianceCovarianceMatrix <-
                                                                  vcov(object = RegressionModel) %>%
                                                                  colnames(x = .) %>%
                                                                  grepl(pattern = "Covariate",x = .) %>%
                                                                  which(x = .)

                                                                # Compute the Wald test with the model coefficients from the logistic regression model,
                                                                # the variance-covariance matrix of the regression model,
                                                                # and the categorical covariate terms from the variance-covariance matrix
                                                                WaldTestForCategoricalCovariateSignificanceResult <-
                                                                  wald.test(
                                                                    b = coef(object = RegressionModel),
                                                                    Sigma = vcov(object = RegressionModel),
                                                                    Terms = IndicesOfTermsToSelectFromVarianceCovarianceMatrix
                                                                      )

                                                                # extract the P-value from the overall significance test for the categorical covariate
                                                                WaldTestForCategoricalCovariateSignificanceResultPvalue <-
                                                                  WaldTestForCategoricalCovariateSignificanceResult$result$chi2["P"] %>%
                                                                  unname(obj = .)

                                                            }
                                                              
                                                             # tidy the regression model results and return them
                                                             ResultsToReturn <-
                                                               TidyLogisticRegressionResults(
                                                                 DataForRegression = DataForRegression,
                                                                 RegressionModel = RegressionModel
                                                                 ) %>%
                                                               mutate(
                                                                 .data = .,
                                                                 # label the current gene
                                                                 Gene = currentGene,
                                                                 # the current phenotype
                                                                 Phenotype = currentPhenotype,
                                                                 # the current variant
                                                                 Variant = currentVariantToIncludeInRegression,
                                                                 # the covariate
                                                                 Covariate = currentCovariate,
                                                                 # the regression model
                                                                 Regression = paste0(
                                                                                     currentPhenotype,
                                                                                     " ~ ",
                                                                                     currentVariantToIncludeInRegression,
                                                                                     " + ",
                                                                                     currentCovariate
                                                                                     ),
                                                                 # the regression type
                                                                 Regression_type = "bivariate"
                                                                 )
                                                             
                                                             # If the current covariate is a categorical covariate,
                                                             # add the P-value from the wald test for overall significance for the categorical covariate,
                                                             # otherwise, label the covariate as non-categorical
                                                             if(
                                                               (currentCovariate == "Sex") |
                                                               (currentCovariate == "StudySeason")
                                                             )
                                                             {

                                                                 ResultsToReturn <-
                                                                     ResultsToReturn %>%
                                                                     mutate(
                                                                       .data = .,
                                                                       CategoricalCovariateSignificance = WaldTestForCategoricalCovariateSignificanceResultPvalue %>%
                                                                                                          signif(x = .,digits = 3) %>%
                                                                                                          AddAstericksToPvalues(columnVector = .)
                                                                            )

                                                             } else {

                                                                 ResultsToReturn <-
                                                                   ResultsToReturn %>%
                                                                   mutate(
                                                                     .data = .,
                                                                     CategoricalCovariateSignificance = "covariate_is_not_categorical"
                                                                          )

                                                             }
                                                             
                                                    } else {
                                                      
                                                        ResultsToReturn <- NULL
                                                    }
                                                    
                                                    return(ResultsToReturn)
                                                    
                                              }
                                                ) %>%
                                              do.call(
                                                what = "rbind",
                                                args = .
                                                )
                                          
                                        return(ResultsToReturn)
                                    }
                                    ) %>%
                                   do.call(
                                     what = "rbind",
                                     args = .
                                     )
                            
                            return(ResultsToReturn)
                        }
                          ) %>%
                        do.call(
                          what = "rbind",
                          args = .
                          )
                 
              } else {
                
                ResultsToReturn <- NULL
              }
              
           return(ResultsToReturn)
              
          }
          ) %>%
      do.call(
        what = "rbind",
        args = .
        )
    
    
##### bind the results of the UnivariateGenotypeLogisticRegressions and BivariateGenotypeLogisticRegressions by row
UnivariateAndBivariateGenotypeLogisticRegressions <-
    list(
      UnivariateGenotypeLogisticRegressions,
      BivariateGenotypeLogisticRegressions
    ) %>%
    # bind the results by row
    do.call(
      what = "rbind",
      args = .
      ) %>%
    data.frame() %>%
    # group by the phenotype and the variant
    group_by(
      .data = .,
      Phenotype,
      Variant
      ) %>%
    # split into a list of dataframes by phenotype and variant
   group_split(.tbl = .) %>%
   do.call(
     what = "rbind",
     args = .
     ) %>%
   data.frame()

##### perform multiple regressions of significant causative variants with all seasonal and demographic covariates included
  
multipleGenotypeRegressions <-
  # loop through each of the genotype matrices with binary phenotypes added
  list.files(path = "./genotype_matrices_with_binary_phenotypes/")  %>%
    lapply(
      X = .,
      FUN = function(currentFileName)
      {
          GenotypeMatrixWithPhenotypeAndCovariateData <-
            # load the genotype matrix with binary phenotype data
            read.delim(
              file = glue("./genotype_matrices_with_binary_phenotypes/{currentFileName}"),
              header = TRUE
              ) %>%
            # join the genotype matrix with covariate data by sequencing sample ID
            left_join(
              x = .,
              y = CovariateData,
              by = c("sampleID" = "VCFfileID")
              )

          # identify the variants that are in the genotype matrix based on the "X" prefix
          VariantsInGenotypeMatrix <-
            GenotypeMatrixWithPhenotypeAndCovariateData %>%
            names(x = .) %>%
            grepl(pattern = "X",x = .) %>%
            which(x = .) %>%
            names(x = GenotypeMatrixWithPhenotypeAndCovariateData)[.]

          # identify the variants in the genotype matrix that have significant associations with a binary phenotype in primary analysis
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

          # perform a logistic regression of phenotype ~ genotype + all seasonal and demographic covariates
          # if there are any significant variants in the genotype matrix
          if(length(x = VariantsToLoopThrough)>0)
          {
             ResultsToReturn <-
                # loop through each variant
                VariantsToLoopThrough %>%
                  lapply(
                    X = .,
                    FUN = function(currentVariantToIncludeInRegression)
                    {
                        # identify the phenotypes to loop through based on the currentVariantToIncludeInRegression
                        PhenotypesToLoopThrough <-
                          CausativeAssociatedVariants %>%
                            filter(
                              .data = .,
                              VariantID_GenotypeMatrixFormat == currentVariantToIncludeInRegression
                              ) %>%
                            pull(
                              .data = .,
                              phenotype
                              )

                        ResultsToReturn <-
                            # loop through each possible phenotype
                            PhenotypesToLoopThrough %>%
                              lapply(
                                X = .,
                                FUN = function(currentPhenotype)
                                {
                                    # select columns for the phenotype, genotype, and all seasonal/demographic covariates
                                    DataForRegression <-
                                      GenotypeMatrixWithPhenotypeAndCovariateData %>%
                                      select(
                                        .data = .,
                                        "Phenotype" = all_of(x = currentPhenotype),
                                        "Genotype" = all_of(x = currentVariantToIncludeInRegression),
                                        "Age" = Age.on.Study.Date,
                                        StudySeason,
                                        BMI,
                                        Sex
                                        ) %>%
                                       # make sure genotype is continuous and phenotype is categorical
                                       mutate(
                                         .data = .,
                                         Genotype = Genotype %>% as.numeric(x = .),
                                         Phenotype = Phenotype %>% factor(x = .)
                                         ) %>%
                                       # remove missing observations
                                      na.omit(object = .)
                                    
                                    # if there is more than one genotype value, perform the regression (e.g., 0 and 1)
                                    # if there is only one genotype value (e.g., 0 only),
                                      # return a NULL
                                    # if there is only one study season as well, return a NULL,
                                    # there must be more than one category to perform a regression with a categorical variable,
                                    # if there is only one gender, return a NULL
                                    if(
                                      length(x = unique(x = as.numeric(x = DataForRegression$Genotype)))>1 &
                                      length(x = unique(x = as.factor(x = DataForRegression$StudySeason)))>1 &
                                      length(x = unique(x = as.factor(x = DataForRegression$Sex)))>1
                                      )
                                    {
                                            # set a random number seed for reproducibility
                                            set.seed(777)
                                            
                                            # obtain the logistic regression model
                                            RegressionModel <-
                                              DataForRegression %>%
                                              glm(
                                                formula = Phenotype ~ Genotype + BMI + StudySeason + Age + Sex,
                                                family = "binomial",
                                                data = .
                                                )
                                            
                                            # calculate P-values for overall term significance for categorical covariates: sex and study season
                                            # with a Wald test
                                            OverallCategoricalCovariateSignificance <-
                                                c(
                                                  "Sex",
                                                  "StudySeason"
                                                  ) %>%
                                                lapply(
                                                  X = .,
                                                  FUN = function(currentCategoricalCovariate)
                                                  {
                                                    # identify the column indices of the variance-covariance matrix
                                                    # columns that correspond to the categorical covariate
                                                    IndicesOfTermsToSelectFromVarianceCovarianceMatrix <-
                                                      vcov(object = RegressionModel) %>%
                                                      colnames(x = .) %>%
                                                      grepl(
                                                        pattern = currentCategoricalCovariate,
                                                        x = .
                                                        ) %>%
                                                      which(x = .)

                                                    # Compute the Wald test with the model coefficients from the logistic regression model,
                                                    # the variance-covariance matrix of the regression model,
                                                    # and the categorical covariate terms from the variance-covariance matrix
                                                    WaldTestForCategoricalCovariateSignificanceResult <-
                                                      wald.test(
                                                        b = coef(object = RegressionModel),
                                                        Sigma = vcov(object = RegressionModel),
                                                        Terms = IndicesOfTermsToSelectFromVarianceCovarianceMatrix
                                                          )

                                                    # extract the P-value from the overall significance test for the categorical covariate
                                                    WaldTestForCategoricalCovariateSignificanceResultPvalue <-
                                                      WaldTestForCategoricalCovariateSignificanceResult$result$chi2["P"] %>%
                                                      unname(obj = .)

                                                    # return the P-value
                                                    return(WaldTestForCategoricalCovariateSignificanceResultPvalue)

                                                  }
                                                    ) %>%
                                              unlist(x = .)

                                            # update the names of the elements of the OverallCategoricalCovariateSignificance array
                                            names(x = OverallCategoricalCovariateSignificance) <- c("Sex","StudySeason")
                                            
                                            # tidy the logistic regression model results and return them 
                                            ResultsToReturn <-
                                              TidyLogisticRegressionResults(
                                                DataForRegression = DataForRegression,
                                                RegressionModel = RegressionModel
                                                ) %>%
                                              # label the current gene
                                              mutate(
                                                .data = .,
                                                Gene = currentGene,
                                                # the current phenotype
                                                Phenotype = currentPhenotype,
                                                # the current variant
                                                Variant = currentVariantToIncludeInRegression,
                                                # overall significance for categorical regression terms (sex and study season)
                                                CategoricalCovariateSignificance_sex = OverallCategoricalCovariateSignificance["Sex"] %>%
                                                                                       signif(x = .,digits = 3) %>%
                                                                                       AddAstericksToPvalues(columnVector = .),
                                                CategoricalCovariateSignificance_studySeason = OverallCategoricalCovariateSignificance["StudySeason"] %>%
                                                                                               signif(x = .,digits = 3) %>%
                                                                                               AddAstericksToPvalues(columnVector = .)
                                                ) 
                                              
                                    } else {
                                              ResultsToReturn <- NULL
                                    }
                                              return(ResultsToReturn)
                                  
                                }
                              ) %>%
                              do.call(
                                what = "rbind",
                                args = .
                                )
                        
                        return(ResultsToReturn)
                      }
                    ) %>%
                  do.call(
                    what = "rbind",
                    args = .
                    )
             
            } else {
              
              ResultsToReturn <- NULL
              
            }
          
        return(ResultsToReturn)
      }
    ) %>%
  do.call(
    what = "rbind",
    args = .
    )
  
##### create a data set for performing polygenic multiple regressions with all significant causative variants (identified via clumping)
      # that were associated with each individual binary phenotype

# create a dataset of independent/causative associated variants with seasonal/demographic covariates for all phenotypes
PolygenicmultipleRegressionDataSet <-
  CausativeAssociatedVariants %>%
    # group the causative associated variants by phenotype
    group_by(
      .data = .,
      phenotype
      ) %>%
    # split into a list of dataframes by phenotype
    group_split(
      .tbl = .
      ) %>%
  lapply(
    X = .,
    FUN = function(currentDataSet)
    {
      # identify the genes with significant variants for the current phenotype
      GenesWithSignificantVariants <- currentDataSet$Gene %>% unique(x = .)
      # identify the current set of variants with the genotype matrix variant ID format
      currentSetOfVariants <- currentDataSet$VariantID_GenotypeMatrixFormat %>% unique(x = .)
      # identify the current phenotype
      currentPhenotype <- currentDataSet$phenotype %>% unique(x = .)
      
      # disable the non-problematic warning message that duplicate columns are introduced when joining dataframes below
      options(warn = -1)
      
        # obtain all significant genotypes for the current phenotype and phenotype values for the current phenotype
        DataToReturn <-
          GenesWithSignificantVariants %>%
          lapply(
            X = .,
            FUN = function(currentGene)
            {
              
                # if the current gene is ugt1a1, change it to ugt1a4 for compatibility with the ugt1a4 genotype matrix file name
                if(currentGene == "UGT1A1")
                {
                  currentGene <- "UGT1A4"
                }
              
                # load the genotype matrix corresponding to the current gene
                GenotypeMatrix <-
                  glue("./genotype_matrices_with_binary_phenotypes/{currentGene}_geno_matrix.raw") %>%
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
                  # select the samples IDs and the current set of significant variants for the current phenotype
                  # and the current phenotype
                  select(
                    .data = .,
                    sampleID,
                    all_of(x = VariantColumnsToSelect),
                    all_of(x = currentPhenotype)
                  )

                return(DataToReturn)
            }
              ) %>%
          # join all of the genotype data together by sample ID,
          # this will give a genotype matrix of all relevant variants for the current phenotype
          Reduce(
            f = JoinDataFrames,
            x = .
            )

        # turn warning messages back on
        options(warn = 0)

          # identify the phenotype column to select, the phenotype column to select will either have the label of the currentPhenotype
          # or the phenotype column to select will have the label of the currentPhenotype with a ".x" suffix in the case when there
          # are duplicated columns present after joining the genotype matrices above
          PhenotypeColumns <-
            names(x = DataToReturn) %>%
            grepl(
              pattern = currentPhenotype,
              x = .
                ) %>%
            which(x = .) %>%
            names(x = DataToReturn)[.]

          # if any of the phenotype columns contains the suffix ".x",
          # select the phenotype column with the .x suffix
          if( any(PhenotypeColumns == paste0(currentPhenotype,".x")) )
          {
            PhenotypeColumnToSelect <-
              paste0(currentPhenotype,".x")
          # otherwise, select the currentPhenotype without the .x suffix
          } else {
            PhenotypeColumnToSelect <- currentPhenotype
          }
          
          DataToReturn <-
              DataToReturn %>%
                # select the sample ID, genotypes, and the current phenotype
                select(
                  .data = .,
                  sampleID,
                  starts_with(match = "X"),
                  !!currentPhenotype := !!PhenotypeColumnToSelect
                  ) %>%
              # create a column with a label for the current phenotype
              mutate(
                .data = .,
                Phenotype = currentPhenotype
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
      
    }
      )

#### perform multiple logistic regressions of phenotype ~ set of all signficantly associated variants for each phenotype with and without covariates

   PolygenicContributionRegressionSummaryWithAndWithoutCovarariates <-
     # loop through each dataset in the PolygenicmultipleRegressionDataSet list
     PolygenicmultipleRegressionDataSet %>%
     lapply(
       X = .,
       FUN = function(currentDataSet)
       {
         
              # remove missing observations from the dataset
              DataForRegression <-
                currentDataSet %>%
                na.omit(object = .) %>%
                # make sure the phenotype is categorical 
                mutate(
                  .data = .,
                  "Phenotype_value" = factor(x = eval(parse(text = currentDataSet$Phenotype)))
                  ) %>%
                # make sure age and BMI are continuous
                # make sure sex and study season are categorical
                mutate(
                  .data = .,
                  Age.on.Study.Date = Age.on.Study.Date %>% as.numeric(x = .),
                  BMI = BMI %>% as.numeric(x = .),
                  Sex = Sex %>% factor(x = .),
                  StudySeason = StudySeason %>% factor(x = .)
                  )
              
              # define the regression model formula with genotypes only
              RegressionModelFormulaGenotypesOnly <-
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
                  "Phenotype_value",
                  .,
                  sep = " ~ "
                  )
              
              # set a random number seed for reproducibility
              set.seed(777)
              
              # perform the multiple logistic regression with genotypes
              RegressionModelGenotypesOnly <-
                  glm(
                    formula = as.formula(object = RegressionModelFormulaGenotypesOnly),
                    family = "binomial",
                    data = DataForRegression
                    )
                
              # tidy the logistic regression results and return them 
              ResultsToReturnGenotypesOnly <-
                TidyLogisticRegressionResults(
                  DataForRegression = DataForRegression,
                  RegressionModel = RegressionModelGenotypesOnly
                  ) %>%
                # add a column with the regression model formula and the phenotype included
                mutate(
                  .data = .,
                  Regression = RegressionModelFormulaGenotypesOnly,
                  Phenotype = unique(x = currentDataSet$Phenotype)
                  ) %>%
                # add a column with a description of the regression
                mutate(
                  .data = .,
                  Regression_descript = "Genotypes only"
                  ) %>%
                mutate(
                  .data = .,
                  # overall significance for categorical regression terms (sex and study season): "not_available" since they weren't included in this regression
                  CategoricalCovariateSignificance_sex = "not_available",
                  CategoricalCovariateSignificance_studySeason = "not_available"
                )
              
              # define the regression model formula with genotypes and all other covariates
              RegressionModelFormulaFullyAdjustedModel <-
                currentDataSet %>%
                names(x = .) %>%
                grepl(pattern = "X",x = .) %>%
                which(x = .) %>%
                names(x = currentDataSet)[.] %>%
                c("Age.on.Study.Date","BMI","Sex","StudySeason") %>%
                paste(
                  .,
                  collapse = " + "
                  ) %>%
                paste(
                  "Phenotype_value",
                  .,
                  sep = " ~ "
                  )
              
              # perform the multiple logistic regression with genotypes and all covariates
              RegressionModelFullyAdjustedModel <-
                  glm(
                    formula = as.formula(object = RegressionModelFormulaFullyAdjustedModel),
                    family = "binomial",
                    data = DataForRegression
                    )
              
              # calculate P-values for overall term significance for categorical covariates: sex and study season
              # with a Wald test
              OverallCategoricalCovariateSignificance <-
                  c(
                    "Sex",
                    "StudySeason"
                    ) %>%
                  lapply(
                    X = .,
                    FUN = function(currentCategoricalCovariate)
                    {
                      # identify the column indices of the variance-covariance matrix
                      # columns that correspond to the categorical covariate
                      IndicesOfTermsToSelectFromVarianceCovarianceMatrix <-
                        vcov(object = RegressionModelFullyAdjustedModel) %>%
                        colnames(x = .) %>%
                        grepl(
                          pattern = currentCategoricalCovariate,
                          x = .
                          ) %>%
                        which(x = .)

                      # Compute the Wald test with the model coefficients from the logistic regression model,
                      # the variance-covariance matrix of the regression model,
                      # and the categorical covariate terms from the variance-covariance matrix
                      WaldTestForCategoricalCovariateSignificanceResult <-
                        wald.test(
                          b = coef(object = RegressionModelFullyAdjustedModel),
                          Sigma = vcov(object = RegressionModelFullyAdjustedModel),
                          Terms = IndicesOfTermsToSelectFromVarianceCovarianceMatrix
                            )

                      # extract the P-value from the overall significance test for the categorical covariate
                      WaldTestForCategoricalCovariateSignificanceResultPvalue <-
                        WaldTestForCategoricalCovariateSignificanceResult$result$chi2["P"] %>%
                        unname(obj = .)

                      # return the P-value
                      return(WaldTestForCategoricalCovariateSignificanceResultPvalue)

                    }
                      ) %>%
                unlist(x = .)

              # update the names of the elements of the OverallCategoricalCovariateSignificance array
              names(x = OverallCategoricalCovariateSignificance) <- c("Sex","StudySeason")
              
              # tidy the logistic regression results and return them
              ResultsToReturnFullyAdjustedModel <-
                TidyLogisticRegressionResults(
                  DataForRegression = DataForRegression,
                  RegressionModel = RegressionModelFullyAdjustedModel
                  ) %>%
                # add a column with the regression model formula and phenotype included
                mutate(
                  .data = .,
                  Regression = RegressionModelFormulaFullyAdjustedModel,
                  Phenotype = unique(x = currentDataSet$Phenotype)
                  ) %>%
                # add a column with a description of the regression
                mutate(
                  .data = .,
                  Regression_descript = "Genotypes and covariates"
                ) %>%
                mutate(
                  .data = .,
                  # overall significance for categorical regression terms (sex and study season)
                  CategoricalCovariateSignificance_sex = OverallCategoricalCovariateSignificance["Sex"] %>%
                                                         signif(x = .,digits = 3) %>%
                                                         AddAstericksToPvalues(columnVector = .),
                  CategoricalCovariateSignificance_studySeason = OverallCategoricalCovariateSignificance["StudySeason"] %>%
                                                                 signif(x = .,digits = 3) %>%
                                                                 AddAstericksToPvalues(columnVector = .)
                )
              
              
              ResultsToReturn <-
                # bind the polygenic regressions with and without covariates together by row
                # and return the results
                list(
                  ResultsToReturnGenotypesOnly,
                  ResultsToReturnFullyAdjustedModel
                ) %>%
                do.call(
                  what = "rbind",
                  args = .
                  ) %>%
                data.frame() 
              
              return(ResultsToReturn)
       }
         )
   
   # name the polygenic multiple regression results according to phenotype
  names(x = PolygenicContributionRegressionSummaryWithAndWithoutCovarariates) <-
                               PolygenicmultipleRegressionDataSet %>%
                                 lapply(
                                   X = .,
                                   FUN = function(currentDataSet)
                                   {
                                     # isolate the phenotype from the current data set
                                     Phenotype <-
                                       unique(x = currentDataSet$Phenotype) %>%
                                       # paste the phenotype with the string "PolygenicRegressionResults"
                                       paste0(
                                         "PolygenicRegressionResults_",
                                         .
                                       )

                                     return(Phenotype)
                                   }
                                     ) %>%
                                  unlist(x = .)
  
##### test common vitamin D binding protein diplotypes for associations with binary vitamin D status phenotypes with logistic regression
  VDBPregressionTestResults <-
      # loop through the possible phenotypes
      c(
        "VitDsufficiency",
        "VitDdeficiency"
        ) %>%
        lapply(
          X = .,
          FUN = function(currentPhenotype)
          {
                # select the phenotype, vitamin D binding protein diplotype, and covariate columns
                DataForRegression <-
                  VitaminDbindingProteinDataForLogisticRegression %>%
                  select(
                    .data = .,
                    "Phenotype" = all_of(x = currentPhenotype),
                    GC_dip,
                    Age.on.Study.Date,
                    BMI,
                    Sex,
                    StudySeason
                    ) %>%
                    # remove missing observations
                   na.omit(object = .)
                
                # set a random number seed for reproducibility
                set.seed(777)
                
                # define the regression model with diplotype only
                RegressionModelUnadjusted <-
                  DataForRegression %>%
                  glm(
                    formula = Phenotype ~ GC_dip,
                    family = "binomial",
                    data = .
                    )
                
                # calculate a P-value for overall significance of the GC diplotype term with a Wald test
                # identify the column indices of the variance-covariance matrix that correspond to the diplotype
                IndicesOfTermsToSelectFromVarianceCovarianceMatrix <-
                  vcov(object = RegressionModelUnadjusted) %>%
                  colnames(x = .) %>%
                  grepl(
                    pattern = "GC_dip",
                    x = .
                    ) %>%
                  which(x = .)
    
                # Compute the Wald test with the model coefficients from the logistic regression model,
                # the variance-covariance matrix of the regression model,
                # and the diplotype terms from the variance-covariance matrix
                WaldTestForCategoricalDiplotypeSignificanceResult <-
                  wald.test(
                    b = coef(object = RegressionModelUnadjusted),
                    Sigma = vcov(object = RegressionModelUnadjusted),
                    Terms = IndicesOfTermsToSelectFromVarianceCovarianceMatrix
                      )
    
                # extract the P-value from the overall significance test for the categorical diplotype
                WaldTestForCategoricalDiplotypeSignificanceResultPvalue <-
                  WaldTestForCategoricalDiplotypeSignificanceResult$result$chi2["P"] %>%
                  unname(obj = .)
                
                # tidy the unadjusted regression results
                ResultsToReturnUnadjusted <-
                  TidyLogisticRegressionResults(
                    DataForRegression = DataForRegression,
                    RegressionModel = RegressionModelUnadjusted
                    ) %>%
                  # add a label for the regression type
                   mutate(
                     .data = .,
                     RegressionType = "unadjusted"
                     ) %>%
                  mutate(
                    .data = .,
                    # overall significance for categorical GC diplotype
                    OverallDiplotypeSignificance = WaldTestForCategoricalDiplotypeSignificanceResultPvalue %>%
                                                   signif(x = .,digits = 3) %>%
                                                   AddAstericksToPvalues(columnVector = .),
                    # overall significance for categorical regression terms (sex and study season): "not_available" since they weren't included in this regression
                    CategoricalCovariateSignificance_sex = "not_available",
                    CategoricalCovariateSignificance_studySeason = "not_available"
                  )
                
                # define the regression model with diplotype and all covariates
                RegressionModelFullyAdjusted <-
                  DataForRegression %>%
                  glm(
                    formula = Phenotype ~ GC_dip + Age.on.Study.Date + BMI + Sex + StudySeason,
                    family = "binomial",
                    data = .
                    )
                
                # calculate a P-value for overall significance of the GC diplotype term with a Wald test
                # identify the column indices of the variance-covariance matrix that correspond to the diplotype
                IndicesOfTermsToSelectFromVarianceCovarianceMatrix <-
                  vcov(object = RegressionModelFullyAdjusted) %>%
                  colnames(x = .) %>%
                  grepl(
                    pattern = "GC_dip",
                    x = .
                    ) %>%
                  which(x = .)
    
                # Compute the Wald test with the model coefficients from the logistic regression model,
                # the variance-covariance matrix of the regression model,
                # and the diplotype terms from the variance-covariance matrix
                WaldTestForCategoricalDiplotypeSignificanceResult <-
                  wald.test(
                    b = coef(object = RegressionModelFullyAdjusted),
                    Sigma = vcov(object = RegressionModelFullyAdjusted),
                    Terms = IndicesOfTermsToSelectFromVarianceCovarianceMatrix
                      )
    
                # extract the P-value from the overall significance test for the categorical diplotype
                WaldTestForCategoricalDiplotypeSignificanceResultPvalue <-
                  WaldTestForCategoricalDiplotypeSignificanceResult$result$chi2["P"] %>%
                  unname(obj = .)
                
                # calculate P-values for overall term significance for categorical covariates: sex and study season
                # with a Wald test
                OverallCategoricalCovariateSignificance <-
                    c(
                      "Sex",
                      "StudySeason"
                      ) %>%
                    lapply(
                      X = .,
                      FUN = function(currentCategoricalCovariate)
                      {
                        # identify the column indices of the variance-covariance matrix
                        # columns that correspond to the categorical covariate
                        IndicesOfTermsToSelectFromVarianceCovarianceMatrix <-
                          vcov(object = RegressionModelFullyAdjusted) %>%
                          colnames(x = .) %>%
                          grepl(
                            pattern = currentCategoricalCovariate,
                            x = .
                            ) %>%
                          which(x = .)
    
                        # Compute the Wald test with the model coefficients from the logistic regression model,
                        # the variance-covariance matrix of the regression model,
                        # and the categorical covariate terms from the variance-covariance matrix
                        WaldTestForCategoricalCovariateSignificanceResult <-
                          wald.test(
                            b = coef(object = RegressionModelFullyAdjusted),
                            Sigma = vcov(object = RegressionModelFullyAdjusted),
                            Terms = IndicesOfTermsToSelectFromVarianceCovarianceMatrix
                              )
    
                        # extract the P-value from the overall significance test for the categorical covariate
                        WaldTestForCategoricalCovariateSignificanceResultPvalue <-
                          WaldTestForCategoricalCovariateSignificanceResult$result$chi2["P"] %>%
                          unname(obj = .)
    
                        # return the P-value
                        return(WaldTestForCategoricalCovariateSignificanceResultPvalue)
    
                      }
                        ) %>%
                  unlist(x = .)
                
                # update the names of the elements of the OverallCategoricalCovariateSignificance array
                names(x = OverallCategoricalCovariateSignificance) <- c("Sex","StudySeason")
                
                # tidy the fully adjusted regression results
                ResultsToReturnFullyAdjusted <-
                  TidyLogisticRegressionResults(
                    DataForRegression = DataForRegression,
                    RegressionModel = RegressionModelFullyAdjusted
                    ) %>%
                 # add a label for the regression type
                  mutate(
                    .data = .,
                    RegressionType = "adjusted"
                    ) %>%
                  mutate(
                    .data = .,
                    # overall significance for categorical GC diplotype
                    OverallDiplotypeSignificance = WaldTestForCategoricalDiplotypeSignificanceResultPvalue %>%
                                                   signif(x = .,digits = 3) %>%
                                                   AddAstericksToPvalues(columnVector = .),
                    # overall significance for categorical regression terms (sex and study season)
                    CategoricalCovariateSignificance_sex = OverallCategoricalCovariateSignificance["Sex"] %>%
                                                           signif(x = .,digits = 3) %>%
                                                           AddAstericksToPvalues(columnVector = .),
                    CategoricalCovariateSignificance_studySeason = OverallCategoricalCovariateSignificance["StudySeason"] %>%
                                                                   signif(x = .,digits = 3) %>%
                                                                   AddAstericksToPvalues(columnVector = .)
                    )
                
                # bind the unadjusted and fully adjusted regression results by row
                ResultsToReturn <- 
                    list(
                      ResultsToReturnUnadjusted,
                      ResultsToReturnFullyAdjusted
                    ) %>%
                    do.call(
                      what = "rbind",
                      args = .
                      ) %>%
                    data.frame() %>%
                    # create a column for the current phenotype
                    mutate(
                      .data = .,
                      Phenotype = currentPhenotype
                      ) 
                
                return(ResultsToReturn)
          }
            )
  
   # name the vitamin D binding protein regression results according to phenotype
  names(x = VDBPregressionTestResults) <-
                    VDBPregressionTestResults %>%
                     lapply(
                       X = .,
                       FUN = function(currentDataSet)
                       {
                         # isolate the phenotype from the current data set
                         Phenotype <-
                           unique(x = currentDataSet$Phenotype) %>%
                           # paste the phenotype with the string "VDBPregressionTestResults_"
                           paste0(
                             "VDBPregressionTestResults_",
                             .
                           )

                         return(Phenotype)
                       }
                         ) %>%
                      unlist(x = .)
  
  
#### save all of the logistic regression results to the file system
  resultsToSave <-
    list(
      "UnivariateAndBivariateDemographicAndSeasonalLogisticRegressions" = UnivariateAndBivariateDemographicAndSeasonalLogisticRegressions,
      "CausativeAssociatedVariants" = CausativeAssociatedVariants,
      "UnivariateAndBivariateGenotypeLogisticRegressions" = UnivariateAndBivariateGenotypeLogisticRegressions,
      "multipleGenotypeRegressions" = multipleGenotypeRegressions
      ) %>%
    c(
      .,
      PolygenicContributionRegressionSummaryWithAndWithoutCovarariates,
      VDBPregressionTestResults
     ) 
  
  mapply(
    FUN = function(currentResultToSave,fileName)
    {
      currentResultToSave %>%
        write_tsv(
          x = .,
          file = glue("./logistic_regressions_with_covariates/{fileName}.txt")
            )
    },
    resultsToSave,
    names(x = resultsToSave)
      )

