#### this is a script for storing all functions that are used in the vitamin D regulatory pathway genetic association analysis

##################################################################################################################################
  ### define a function for adding levels of significance to P-value columns with astericks
    # this function is designed to be used with tidyverse mutate on an input column vector
  AddAstericksToPvalues <-
    function(
      columnVector = NULL # a column vector from a dataframe
      )
    {
        VectorToReturn <-
          columnVector %>%
          # convert to a numeric so the if-else statements below can be evaluated
          as.numeric(x = .) %>%
          lapply(
            X = .,
            FUN = function(currentPvalue)
            {
                # if the current P-value is not NA, add astericks for level of significance,
                # otherwise, return an NA
                if(!is.na(x = currentPvalue))
                {
                     # if the P-value is less than 0.001, add 3 astericks
                      if (currentPvalue < 0.001) {
                        ValueToReturn <-
                          paste0(currentPvalue,"***")
                      } # if the P-value is less than 0.01, add 2 astericks
                      else if (currentPvalue < 0.01) {
                        ValueToReturn <-
                          paste0(currentPvalue,"**")
                      } # if the P-value is less than 0.05, add 1 astericks
                      else if (currentPvalue < 0.05)
                    {
                      ValueToReturn <-
                        paste0(currentPvalue,"*")
                    } else {
                      ValueToReturn <- as.character(x = currentPvalue)
                    }

                } else {
                  ValueToReturn <- as.character(x = NA)
                }

                 return(ValueToReturn)
          }
          ) %>%
          unlist(x = .)

        return(VectorToReturn)
    }


##################################################################################################################################
  ### define a function for calculating percentage of participants within a group of participants with 
  # vitamin D levels below sufficiency, levels of insufficiency, and deficiency
  # example groups could be gender or seasonal groupings

  # calculate percentage of participants grouped by gender that have:
      # (1) 25OHD levels below sufficiency (≤ 20 ng/mL)
      # (2) insufficient 25OHD levels (between 12 and 20 ng/mL)
      # (3) 25OHD deficiency (below 12 ng/mL)
    CalculateVitaminDstatusByGroup <-
      function(
        DataFrameToUse = NULL, # an arbitrary dataframe 
        GroupingVariable = NULL # a string with the column name to group the results by from the arbitrary dataframe 
        )
      {
        
        ResultsToReturn <-
          DataFrameToUse %>%
            # select the overall vitamin D concentration column and the GroupingVariable
            select(
              .data = .,
              !!GroupingVariable,
              Overall_primary_25OHD_ngperml
              ) %>%
            # remove missing observations
            na.omit(object = .) %>%
            # create binarized columns for BelowSufficiency, Insufficiency, and Deficiency
            mutate(
              .data = .,
              Insufficiency_percent = ((Overall_primary_25OHD_ngperml >= 12) & (Overall_primary_25OHD_ngperml <= 20)),
              Deficiency_percent = Overall_primary_25OHD_ngperml < 12
              ) %>%
            # group by GroupingVariable
            group_by(
              .data = .,
              !!rlang::sym(x = GroupingVariable)
              ) %>%
            # deselect the concentration column now that it is no longer needed
            select(
              .data = .,
              -Overall_primary_25OHD_ngperml
              ) %>%
            group_split(.tbl = .) %>%
            # calculate the percentage for each category
            lapply(
              X = .,
              FUN = function(currentDataFrame)
              {
                
                # ungroup the current data frame
                currentDataFrame <-
                  currentDataFrame %>%
                  ungroup(x = .)
                
                DataToReturn <-
                  currentDataFrame %>%
                  # deselect the GroupingVariable
                  select(
                    .data = .,
                    -!!GroupingVariable
                    ) %>%
                    # loop through each boolean vitamin D status column and calculate percentages
                    lapply(
                      X = .,
                      FUN = function(currentColumn)
                      {
                        
                        SampleSize <-
                          # find the sample size based on observations without missing data
                          DataFrameToUse %>%
                          select(
                            .data = .,
                            !!GroupingVariable,
                            Overall_primary_25OHD_ngperml
                          ) %>%
                          na.omit(object = .) %>%
                          # filter to observations that match the GroupingVariable of the currentDataFrame
                          filter(
                            .data = .,
                            !!rlang::sym(x = GroupingVariable) == eval(parse(text = paste0("unique(x = currentDataFrame$",GroupingVariable,")"))) 
                            ) %>%
                          nrow(x = .)

                        PercentToReturn <-
                          round(x = (sum(currentColumn,na.rm = TRUE)/SampleSize)*100,digits = 1)

                        return(PercentToReturn)
                      }
                      ) %>%
                  # bind results together by column
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
                  # create a column with a label for the current GroupingVariable
                  mutate(
                    .data = .,
                    Group = eval(parse(text = paste0("unique(x = currentDataFrame$",GroupingVariable,")"))) 
                    ) %>%
                  # create a column with the sample size
                  mutate(
                    .data = .,
                    total_N = DataFrameToUse %>%
                              select(
                                .data = .,
                                !!GroupingVariable,
                                Overall_primary_25OHD_ngperml
                              ) %>%
                              na.omit(object = .) %>%
                              # filter to observations that match the GroupingVariable of the currentDataFrame
                              filter(
                                .data = .,
                                !!rlang::sym(x = GroupingVariable) == eval(parse(text = paste0("unique(x = currentDataFrame$",GroupingVariable,")"))) 
                              ) %>%
                              nrow(x = .)
                    )

                return(DataToReturn)

              }
            ) %>%
          do.call(
            what = "rbind",
            args = .
            )
          
      }

##################################################################################################################################
    ### define a function for converting the column of a dataframe that is in picograms per mL to something in nanograms per mL
    # unit conversion: 1 pg * (1 g/1,000,000,000,000 pg) * (1,000,000,000 ng/1 g)
    ConvertPicoGramsToNanoGrams <-
      function(
        ColumnToUse = NULL # this is a numeric vector in units of picograms per mL
      )
      {
        
        ArrayToReturn <-
          # convert the array to a single column tibble
          ColumnToUse %>%
          as.numeric(x = .) %>%
          tibble(ColumnToUse = .) %>%
          # convert the units of the single column from picogram/mL to nanogram/mL,
          # if the number is not an NA
          mutate(
            .data = .,
            ColumnToUse = if_else(
                                  condition = is.na(x = ColumnToUse),
                                  true = ColumnToUse,
                                  false = ColumnToUse*(1/1000000000000)*(1000000000)
                                )
          ) %>%
          pull(
            .data = .,
            ColumnToUse
          )
        
        return(ArrayToReturn)
        
      }
    
######################################################################################################################################
    
    ##### this is a function for creating files of snps that will be used to filter VCF files for locusZoom plots
          # this will create a list of snps for each individual gene and phenotype 
    CreateLocusZoomSNPlists <-
      function(
        DataFrameToUse = NULL, # this takes in a dataframe of either linear or logistic regression results computed with PLINK for multiple genes and phenotypes
        PhenotypeColumnName = NULL, # a string with the name of the phenotype column in the DataFrameToUse
        GeneColumnName = NULL # a string with the name of the gene column in the DataFrameToUse
      )
      {
          # save a list of variant IDs to filter the VCF file used for calculating CSKT linkage disequilibrium
          # that are the same exact variants that are plottted for the LocusZoom association results
          # for multi-allelic sites, the variant is selected based on the lowest P-value for the multi-allelic site
          # do this for each phenotype, for each gene
          DataFrameToUse %>%
            # group the data frame by gene and phenotype
            group_by(
              .data = .,
              eval(parse(text = PhenotypeColumnName)),
              eval(parse(text = GeneColumnName))
            ) %>%
            # split into a list of dataframes based on the groupings
            group_split(
              .tbl = .,
              .keep = TRUE
            ) %>%
            lapply(
              X = .,
              FUN = function(currentDataSet)
              {
                FileToSave <-
                  currentDataSet %>%
                  # ensure that there are no groupings that remain in the current data set
                  ungroup(x = .) %>%
                  # select the SNP and Pvalue columns
                  select(
                    .data = .,
                    SNP,
                    Pvalue
                  ) %>%
                  # change the format of the SNP  from chr:pos:ref:alt to chr:pos as required for locusZoom in a new column named markerName
                  # remove the ref and alt allele suffix from the SNP
                  mutate(
                    .data = .,
                    MarkerName = SNP %>%
                      # replace the ":ref_allele:alt_allele" suffix with an empty string
                      # the alt allele can be a string of capital letters for SNPs and insertions or
                      # an astericks (*) (punctuation) for deletions
                      gsub(
                        pattern = "[[:punct:]]+[[:upper:]]+[[:punct:]]+([[:upper:]]+|[[:punct:]]+)",
                        replacement = "",
                        x = .
                      )
                  ) %>%
                  # remove any rows with an NA P-value
                  na.omit(object = .)
  
                # if there are any SNPs remaining after removing NA rows, save the results,
                # otherwise do not save any file
                if(nrow(x = FileToSave)>0)
                {
                  FileToSave <-
                    FileToSave %>%
                    # filter to the lowest P-value for the marker name for multi-allelic sites that have the same chr:position marker name
                    group_by(
                      .data = .,
                      MarkerName
                    ) %>%
                    filter(
                      .data = .,
                      Pvalue == min(Pvalue,na.rm = TRUE)
                    ) %>%
                    # ungroup by markerName now that all variants are unique
                    ungroup(x = .) %>%
                    # select the SNP column
                    select(
                      .data = .,
                      SNP
                    ) %>%
                    # make sure all SNPs are unique
                    unique(x = .)
  
                  # obtain the label of the current phenotype and gene
                  currentPhenotype <- eval(parse(text = paste0("currentDataSet$",PhenotypeColumnName))) %>% unique(x = .)
                  currentGene <- eval(parse(text = paste0("currentDataSet$",GeneColumnName))) %>% unique(x = .)
  
                  # make a directory for saving the results
                  if(!dir.exists(paths = glue("./Results_for_Calculating_LD_for_LocusZoom/")))
                  {
                    dir.create(path = glue("./Results_for_Calculating_LD_for_LocusZoom/"))
                  }
  
                  # save the results to the directory with the file labeled with the phenotype and the gene
                  FileToSave %>%
                    write_tsv(
                      x = .,
                      file = glue(
                        "./Results_for_Calculating_LD_for_LocusZoom/SNPsForLDcalculation_{currentPhenotype}_{currentGene}.txt"),
                      col_names = FALSE
                    )
                  # print a message saying the results were saved successfully
                  print(x = glue("SNPs for calculating linkage disequilibrium for {currentPhenotype} {currentGene} were saved"))
  
                } else {
  
                  # obtain the label of the current phenotype and gene
                  currentPhenotype <- eval(parse(text = paste0("currentDataSet$",PhenotypeColumnName))) %>% unique(x = .)
                  currentGene <- eval(parse(text = paste0("currentDataSet$",GeneColumnName))) %>% unique(x = .)
  
                  print(x = glue("There are no SNPs do not have NA P-values for {currentPhenotype} {currentGene}"))
                  print(x = "No file was saved")
                }
  
  
              }
            )
      }
      
    
######################################################################################################################################
    
    #### this is a function for creating files with snps and P-values for LocusZoom for individual genes and phenotypes contained within the input DataFrameToUse
    CreateMarkerNameAndPvalueFilesForLocusZoom <-
      function(
        DataFrameToUse = NULL, # this takes in a dataframe of either linear or logistic regression results computed with PLINK for variants across multiple genes and phenotypes
        PhenotypeColumnName = NULL, # a string with the name of the phenotype column in the DataFrameToUse
        GeneColumnName = NULL # a string with the name of the gene column in the DataFrameToUse
      )
      {
        
          # save a text file with the variant name (MarkerName) and p values (P-value) for LocusZoom
          # for each phenotype and for each gene
          DataFrameToUse %>%
            # group the data frame by gene and phenotype
            group_by(
              .data = .,
              eval(parse(text = PhenotypeColumnName)),
              eval(parse(text = GeneColumnName))
            ) %>%
            # split into a list of dataframes based on the groupings
            group_split(
              .tbl = .,
              .keep = TRUE
            ) %>%
            lapply(
              X = .,
              FUN = function(currentDataSet)
              {
                FileToSave <-
                  currentDataSet %>%
                  # ensure that there are no groupings that remain in the current data set
                  ungroup(x = .) %>%
                  # select the marker name and P-value columns
                  select(
                    .data = .,
                    SNP,
                    Pvalue
                  )  %>%
                  # remove any rows with an NA P-value
                  na.omit(object = .) %>%
                  mutate(
                    .data = .,
                      # replace the ":ref_allele:alt_allele" suffix with an empty string
                      # the alt allele can be a string of capital letters for SNPs and insertions or
                      # an astericks (*) (punctuation) for deletions
                    MarkerName = SNP %>%
                      gsub(
                        pattern = "[[:punct:]]+[[:upper:]]+[[:punct:]]+([[:upper:]]+|[[:punct:]]+)",
                        replacement = "",
                        x = .
                        )
                    )
  
                # if there are any SNPs remaining after removing NA rows, save the results,
                # otherwise do not save any file
                if(nrow(x = FileToSave)>0)
                {
  
                  FileToSave <-
                    FileToSave %>%
                    # filter to the lowest P-value for the marker name for multi-allelic sites that have the same chr:position marker name
                    group_by(
                      .data = .,
                      MarkerName
                    ) %>%
                    filter(
                      .data = .,
                      Pvalue == min(Pvalue,na.rm = TRUE)
                    ) %>%
                    # ungroup by markerName now that all variants are unique
                    ungroup(x = .) %>%
                    # drop the MarkerName column
                    select(
                      .data = .,
                      -MarkerName
                      ) %>%
                    # change the name of P to P-value as required for locusZoom
                    # rename the SNP column as MarkerName as required for locusZoom
                    select(
                      .data = .,
                      "MarkerName" = SNP,
                      "P-value" = Pvalue
                    )
                  
                  # obtain the label of the current phenotype and gene
                  currentPhenotype <- eval(parse(text = paste0("currentDataSet$",PhenotypeColumnName))) %>% unique(x = .)
                  currentGene <- eval(parse(text = paste0("currentDataSet$",GeneColumnName))) %>% unique(x = .)
  
                  # make a directory for saving the results if it does not already exist
                  if(!dir.exists(paths = glue("./Results_for_LocusZoom/")))
                  {
                    dir.create(path = glue("./Results_for_LocusZoom/"))
                  }
  
  
  
                  # save the results to the directory with the file labeled with the phenotype and the gene
                  FileToSave %>%
                    write_tsv(
                      x = .,
                      file = glue("./Results_for_LocusZoom/readyToPlot_{currentPhenotype}_{currentGene}_Results_for_LocusZoom.txt")
                    )
                  # print a message saying the results were saved successfully
                  print(x = glue("LocusZoom results for {currentPhenotype} {currentGene} were saved"))
  
                  # create a file with the SNP with the lowest P-value
                  # this SNP will be used as the reference SNP for calculating linkage disequilibrium (LD for all SNPs will be computed with reference to this one)
                  ReferenceSNP <-
                    currentDataSet %>%
                    # ensure that there are no groupings that remain in the current data set
                    ungroup(x = .) %>%
                    mutate(
                      .data = .,
                        # replace the ":ref_allele:alt_allele" suffix with an empty string
                        # the alt allele can be a string of capital letters for SNPs and insertions or
                        # an astericks (*) (punctuation) for deletions
                      chr_pos = SNP %>%
                        gsub(
                          pattern = "[[:punct:]]+[[:upper:]]+[[:punct:]]+([[:upper:]]+|[[:punct:]]+)",
                          replacement = "",
                          x = .
                          )
                      ) %>%
                    # select the marker name and P-value columns
                    select(
                      .data = .,
                      chr_pos,
                      "MarkerName" = SNP,
                      "P-value" = Pvalue
                    ) %>%
                    # remove any variants with NA p-values
                    na.omit(object = .) %>%
                    # filter to variants with the lowest P-value for multiallelic sites
                    group_by(
                      .data = .,
                      chr_pos
                      ) %>%
                    filter(
                      .data = .,
                      `P-value` == min(`P-value`,na.rm = TRUE)
                    ) %>%
                    # ungroup the table
                    ungroup(x = .) %>%
                    # now filter to the variant with the lowest P-value for ALL variants
                    filter(
                      .data = .,
                      `P-value` == min(`P-value`,na.rm = TRUE)
                    ) %>%
                    select(
                      .data = .,
                      MarkerName
                      ) %>%
                    unique(x = .)
  
                  # if it is the case that there is more than one SNP that is tied for the lowest P-value,
                  # arbitrarily choose the first SNP to be the reference SNP for the linkage disequilibrium calculation
                  # the other variants that all have the same P-values as the reference SNP are likely in tight linkage disequilibrium with the reference
                  if(nrow(x = ReferenceSNP) > 1)
                  {
                    ReferenceSNP <-
                      ReferenceSNP[1,]
                  }
  
                  # create a directory for storing the reference SNP
                  if(!dir.exists(paths = glue("./ReferenceSNPs_for_LinkDis/")))
                  {
                    dir.create(path = glue("./ReferenceSNPs_for_LinkDis/"))
                  }
  
                  # save the results to the directory with the file labeled with the phenotype and the gene
                  ReferenceSNP %>%
                    write_tsv(
                      x = .,
                      file = glue("./ReferenceSNPs_for_LinkDis/referenceSNP_forLDcalc_{currentPhenotype}_{currentGene}.txt"),
                      col_names = FALSE
                        )
  
                } else {
  
                  # obtain the label of the current phenotype and gene
                  currentPhenotype <- eval(parse(text = paste0("currentDataSet$",PhenotypeColumnName))) %>% unique(x = .)
                  currentGene <- eval(parse(text = paste0("currentDataSet$",GeneColumnName))) %>% unique(x = .)
  
                  print(x = glue("There are no SNPs that do not have NA P-values for {currentPhenotype} {currentGene}"))
                  print(x = "No file was saved")
                }
  
  
              }
            )
      }

###################################################################################################################################### 
    
    #### create a function that converts a PLINK variant ID column that has the format chrom:pos:ref:alt to the format 
    # of a PLINK genotype matrix (Xchrom.pos.ref.alt_alt)
    ConvertPLINKvariantIDtoGenotypeMatrixFormat <-
      function(
        DataFrameToUse = NULL # an arbitrary dataframe with a SNP column with a variant ID in the format chrom:pos:ref:alt
      ) 
      {
          # create a column with the PLINK variant ID converted to the format that the genotype matrices are in
          # the PLINK variant ID column has the format chrom:pos:ref:alt
          # the genotype matrices have the format Xchrom.pos.ref.alt_alt
          DataToReturn <-
              DataFrameToUse %>%
              mutate(
                .data = .,
                VariantID_GenotypeMatrixFormat = SNP %>%
                                                 lapply(
                                                   X = .,
                                                   FUN = function(currentVariant)
                                                   {
                                                     # obtain the prefix string, which is an "X" concatenated with
                                                     # the current PLINK variant ID with colons replaced with periods
                                                     PrefixString <-
                                                     currentVariant %>%
                                                       # replace all colons with a period
                                                       gsub(
                                                         pattern = ":",
                                                         replacement = ".",
                                                         x = .
                                                       ) %>%
                                                       # append an X to the front of the string
                                                       paste0("X",.)

                                                     # split the current PLINK variant ID into an array based on the colon
                                                     # find the index of the last element in the resulting array, which is the index for the alternate allele
                                                     LengthOfCurrentVariantStringAfterSplittingOnColon <-
                                                     # split the current variant string into an array based on the colon
                                                       currentVariant %>%
                                                       str_split(string = .,pattern = ":") %>%
                                                       unlist(x = .) %>%
                                                      # obtain the length of the array
                                                       length(x = .)

                                                     # obtain the alternate allele using the index of the alternate allele
                                                     # from the current variant ID split on the colon
                                                     AlternateAllele <-
                                                       currentVariant %>%
                                                       str_split(string = .,pattern = ":") %>%
                                                       unlist(x = .) %>%
                                                       .[LengthOfCurrentVariantStringAfterSplittingOnColon]

                                                     # paste the prefix string and alternate allele together separated by an underscore
                                                     StringToReturn <-
                                                       paste0(PrefixString,"_",AlternateAllele)

                                                     # if the final string to return contains an asterisk (the symbol for a deletion),
                                                     # replace all astericks with a period
                                                     if(grepl(pattern = "\\*",x = StringToReturn)==TRUE)
                                                     {
                                                       StringToReturn <-
                                                         StringToReturn %>%
                                                         gsub(pattern = "\\*",replacement = ".",x = .)
                                                     }

                                                     return (StringToReturn)

                                                   }
                                                     ) %>%
                                                unlist(x = .)
          
                )
        
        
                return(DataToReturn)
      }
    
################################################################################################################################## 
   
    ### define a function for cleaning the columns that have metabolite concentrations
    # By cleaning, I mean.......
      # change all of the BLQs to NA,
      # the NDs to NA
      # anything with a less than sign to NA as well
      # anything with a N/A to NA
      
    DataCleaningFunctionToApplyToAllMetaboliteColumns <-
      function(columnName = NULL)
      {
        ValuesToReturn <-
            columnName %>%
              # convert to a character vector
              as.character(x = .) %>%
              # convert the column name array to a tibble of a single column
              tibble(columnName = .) %>%
              # clean the single column based on the if_else statements below
              mutate(
                .data = .,
                columnName =
                  if_else(
                    # replace any NDs with a character NA value
                    condition = columnName == "ND",
                    true = as.character(x = NA),
                    false =
                      if_else(
                        # replace any BLQs with a character NA value
                        condition = columnName == "BLQ",
                        true = as.character(x = NA),
                        false =
                          if_else(
                            # replace any blqs with a character NA value
                            condition = columnName == "blq",
                            true = as.character(x = NA),
                            false =
                              if_else(
                                # replace any less than symbol followed by a space with a character NA value
                                condition = startsWith(x = columnName,prefix = "< "),
                                true = as.character(x = NA),
                                false =
                                  if_else(
                                    # replace any less than symbol without a space with a character NA value
                                    condition = startsWith(x = columnName,prefix = "<"),
                                    true = as.character(x = NA),
                                    false =
                                      if_else(
                                        # replace any "N/A" with an "NA"
                                        condition = columnName == "N/A",
                                        true = as.character(x = NA),
                                        false = columnName
                                      )
                                  )
                              )
                          )
                      )
                  )
              ) %>%
              # convert the column name column to an array
              pull(
                .data = .,
                columnName
              ) %>%
              # convert the array to a numeric
              as.numeric(x = .) 
        
        return(ValuesToReturn)
      }
    
##################################################################################################################################
    
    
    ### create a function for joining two dataframes based on a join column with the merge function
      # this function is used with the Reduce function to join a list of dataframes
    JoinDataFrames <- 
      function(
        df1 = NULL, # arbitrary dataframe
        df2 = NULL # arbitrary dataframe
          )
        { 
           DataFrameToReturn <- 
             merge(
               df1, 
               df2, 
               by = "sampleID"
               ) 
           
           return(DataFrameToReturn)
           
        }
    
    
##############################################################################################################################
    
    # This is modified function code from Michael Sach's original function [ggplot.cosinor.lm()] for plotting of the original data points with the fitted cosinor model.
    # The code for this function is from Michael Sach's cosinor package, I have just made slight modifications to it.
    # The user must supply the cosinor model object from fitting the model to their data and an array of time points from their original data set.
    # The user must also supply a number between 0 and 1 for the scatter plot point transparency -- a value of 1 is max transparency
    CosinorPlottingFunction_ModifiedCodeFromMichaelSachs <-
      function (
        object, 
        x_str = NULL,
        originalTimePointArray = NULL,
        scatterPlotPointTransparency=NULL
        )
      {
        # modification by Jack Staples: changed timeax from 0 to the period in increments of 1 instead of 200 uniform increments
        timeax <- seq(0, object$period) #, length.out = 200)
        covars <- grep("(rrr|sss)", attr(object$fit$terms, "term.labels"),
                       invert = TRUE, value = TRUE)
        newdata <- data.frame(time = timeax, rrr = cos(2 * pi * timeax/object$period),
                              sss = sin(2 * pi * timeax/object$period))
        for (j in covars) {
          newdata[, j] <- 0
        }
        if (!is.null(x_str)) {
          for (d in x_str) {
            tdat <- newdata
            tdat[, d] <- 1
            newdata <- rbind(newdata, tdat)
          }
          newdata$levels <- ""
          for (d in x_str) {
            newdata$levels <- paste(newdata$levels, paste(d,
                                                          "=", newdata[, d]))
          }
        }

        newdata$Y.hat <- predict(object$fit, newdata = newdata)
        
        ### JS code:
          # identify the time of the peak (max) and trough (min) cosinor model predictions
          PeakPredictionTime <<-
            # filter the newdata with cosinor model predicted values to the maximum Y.hat value
            newdata %>%
            filter(
              .data = ., 
              Y.hat == max(Y.hat,na.rm = TRUE)  
            ) %>%
            # pull out the time
            pull(
              .data = .,
              time
              )
        
          TroughPredictionTime <<-
            # filter the newdata with cosinor model predicted values to the minimum Y.hat value
            newdata %>%
              filter(
                .data = .,
                Y.hat == min(Y.hat,na.rm = TRUE)
              ) %>%
            # pull out the time
            pull(
              .data = .,
              time
              )
            
            # create a dataframe of the original data points to plot
            originalDataPointsToPlot <-
              data.frame(
                # the actual time values from the data
                "time" = originalTimePointArray,
                # This is the actual metabolite data values, not the predicted metabolite data values,
                # they are just given the Y.hat label here for 'predicted' so they can be plotted correctly 
                # with the predicted metabolite values from the cosinor curve in the plotting code below
                "Y.hat" = object$fit$model$Y
              )
            
        ###

        if (missing(x_str) || is.null(x_str)) {

          ggplot(newdata, aes_string(x = "time", y = "Y.hat")) +
            geom_line(colour = "darkred",size = 1) +
            # JS code for adding data points as a scatter plot
            geom_point(
              data = originalDataPointsToPlot,
              colour = "black",
              fill = "#A4A4A4",
              shape = 21,
              size = 3.5,
              stroke = 1,
              alpha = scatterPlotPointTransparency
            )
        }
        else {

          ggplot(newdata, aes_string(x = "time", y = "Y.hat", col = "levels")) +
            geom_line() +
            # JS code for adding data points as a scatter plot
            geom_point(
              data = originalDataPointsToPlot,
              colour = "black",
              fill = "#A4A4A4",
              shape = 21,
              size = 3.5,
              stroke = 1,
              alpha = scatterPlotPointTransparency
            )
        }


      }
    
##################################################################################################################################    

    # Define a function where a user can select a metabolite of interest from 
    # the metabolite data and a cosinor model will be fit with no covariates included.
    # The cosinor model can be fit to all the data grouped by month OR 
    # the median of each month if the FitToMedian parameter is TRUE.
    CreateCosinorModelWithNoCovariatesGroupedByMonth <-
      function(
        metaboliteOfInterest = NULL, # a string with the name of the metabolite of interest from the metabolite data (example: "primary_25OHD3_ngperml")
        MetaboliteDataToUse = NULL, # metabolite data dataframe
        periodForOneWaveCycle = 12, # period for one sinusoidal wave cycle (12 months)
        FitToMedian = FALSE # option to fit the cosinor model to the median of each month instead of all of the data grouped by month
      )
      {

        # if the FitToMedian parameter is set to TRUE, fit the cosinor model to the median for each individual month,
        # if the FitToMedian parameter is set to FALSE, fit the cosinor model to ALL of the data, grouped by month
       if(FitToMedian == TRUE)
       {

         DataToFitCosinorModelTo <-
           # select the metabolite of interest and time columns, time is the month in this case
           MetaboliteDataToUse %>%
           select(
             .data = .,
             "Y" = all_of(x = metaboliteOfInterest),
             "time" = StudyMonth
           ) %>%
           # ensure that the Y and time are numerics
           mutate(
             .data = .,
             Y = Y %>% as.numeric(x = .),
             time = time %>% as.numeric(x = .)
           ) %>%
           # remove any rows that contain NA values
           na.omit(object = .) %>%
           # calculate the median metabolite concentration for each month
           group_by(
             .data = .,
             time
           ) %>%
           mutate(
             .data = .,
             MedianConcentration = median(x = Y)
           ) %>%
           # ungroup the table after calculating the median
           ungroup(x = .) %>%
           # drop the column with the metabolite values
           select(
             .data = .,
             -Y
           ) %>%
           # rename the median concentration column as Y
           select(
             .data = .,
             "Y" = MedianConcentration,
             time
           ) %>%
           # obtain unique rows only
           distinct(.data = .) %>%
           # convert to a dataframe, the cosinor.lm function will not accept a tibble
           as.data.frame(x = .)

       } else {

         DataToFitCosinorModelTo <-
           # select the metabolite of interest and time columns, time is the month in this case
           MetaboliteDataToUse %>%
           select(
             .data = .,
             "Y" = all_of(x = metaboliteOfInterest),
             "time" = StudyMonth
           ) %>%
           # ensure that the Y and time are numerics
           mutate(
             .data = .,
             Y = Y %>% as.numeric(x = .),
             time = time %>% as.numeric(x = .)
           ) %>%
           # remove any rows that contain NA values
           na.omit(object = .) %>%
           # convert to a dataframe, the cosinor.lm function will not accept a tibble
           as.data.frame(x = .)
       }

        # fit the cosinor model with the data using the user defined period
        CosinorModel <-
          DataToFitCosinorModelTo %>%
          cosinor.lm(
            formula = Y ~ time(time),
            data = .,
            period = periodForOneWaveCycle
          )

        # calculate the sum of squared residuals (SSR) and the sum of squares total (SST)
        SSR <- ((CosinorModel$fit$residuals)^2) %>% sum(.,na.rm = TRUE)
        SST <- ((DataToFitCosinorModelTo$Y - mean(x = DataToFitCosinorModelTo$Y,na.rm = TRUE))^2) %>% sum(.,na.rm = TRUE)

        # calculate the coefficient of determination (R-squared) as 1-(SSR/SST)
        R.squared <- (1-(SSR/SST)) %>% round(x = .,digits = 4)
        
        # obtain the transformed model coeffiencts from the cosinor regression summary,
        # which includes the annual mean concentration estimate, the amplitude, and acrophase
        CosinorModelSummary <-
          CosinorModel %>%
          summary(object = .) %>%
          .$transformed.table %>%
          # convert rownames to a column
          mutate(
            .data = .,
            Term = rownames(x = .)
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

        # identify the cosinor model predicted annual mean concentration
        AnnualMeanConcentration <-
          CosinorModelSummary %>%
          filter(
            .data = .,
            Term == "(Intercept)"
            ) %>%
          pull(
            .data = .,
            estimate_standard.error
            )

        # identify the cosinor model predicted amplitude
        Amplitude <-
          CosinorModelSummary %>%
          filter(
            .data = .,
            Term == "amp"
          ) %>%
          pull(
            .data = .,
            estimate_standard.error
          )

        # identify the mean minimum and maximum values from the cosinor model fitted values
        MeanMinimumConc <-
          CosinorModel$fit$fitted.values %>% min()

        MeanMaximumConc <-
          CosinorModel$fit$fitted.values %>% max()
        
        # if the metabolite being plotted is the "primary_25OHD3_ngperml" or the "Overall_primary_25OHD_ngperml" metabolite,
        # conditionally color the points of the graph based on vitamin D toxicity, sufficiency, insufficiency, and deficiency
        if(metaboliteOfInterest=="primary_25OHD3_ngperml" | metaboliteOfInterest=="Overall_primary_25OHD_ngperml")
        {
          VitaminDstatusCategories <-
            MetaboliteDataToUse %>%
            select(
              .data = .,
              "Y" = all_of(x = metaboliteOfInterest),
              "time" = StudyMonth
            ) %>%
            # ensure that the Y and time are numerics
            mutate(
              .data = .,
              Y = Y %>% as.numeric(x = .),
              time = time %>% as.numeric(x = .)
            ) %>%
            # remove any rows that contain NA values
            na.omit(object = .) %>%
            # create a column with vitamin D status categories
            mutate(
              .data = .,
                        primaryDstatus = if_else(
                                                  # if the value is less than 12, code as deficient
                                                  condition = Y < 12,
                                                  true = "Deficient (< 12 ng/mL)",
                                                  false = if_else(
                                                                  # if the value is between 12 and 20, code as insufficient
                                                                  condition = Y >= 12 & Y <= 20,
                                                                  true = "Insufficient (between 12 & 20 ng/mL)",
                                                                  false = if_else(
                                                                                  # if the value is greater than 50, code as high
                                                                                  condition = Y > 50,
                                                                                  true = "High (> 50 ng/mL)",
                                                                                  # otherwise code as sufficient
                                                                                  false = "Sufficient (> 20 & ≤ 50 ng/mL)"
                                                                                )
                                                                    )
                                                      )
            ) %>%
            # pull out the primaryDstatus to an array
            pull(.data = .,primaryDstatus)
        }
        
        # plot the cosinor model fit with the actual data points overlayed
        # I have modified Michael Sach's original plotting function code in the cosinor package to allow for this
        PlotToReturn <-
          CosinorModel %>%
          CosinorPlottingFunction_ModifiedCodeFromMichaelSachs(
            object = .,
            x_str = NULL,
            originalTimePointArray = DataToFitCosinorModelTo$time,
            scatterPlotPointTransparency = 0
          ) +
          # add a layer of points for the original data
          geom_point(
            mapping = aes(
              x = time,
              y = Y,
              # if the metabolite of interest is the primary D metabolite,
              # color the graphed points by clinical D status defined above,
              # otherwise do not add color groups to the points
              color = if(metaboliteOfInterest=="primary_25OHD3_ngperml" | metaboliteOfInterest=="Overall_primary_25OHD_ngperml")
              {
                factor(
                  x = VitaminDstatusCategories,
                  levels = c(
                    "High (> 50 ng/mL)",
                    "Sufficient (> 20 & ≤ 50 ng/mL)",
                    "Insufficient (between 12 & 20 ng/mL)",
                    "Deficient (< 12 ng/mL)"
                  )
                )
              } else {
                NULL
              }
            ),
            data = MetaboliteDataToUse %>%
                    select(
                      .data = .,
                      "Y" = all_of(x = metaboliteOfInterest),
                      "time" = StudyMonth
                    ) %>%
                  # ensure that the Y and time are numerics
                  mutate(
                    .data = .,
                    Y = Y %>% as.numeric(x = .),
                    time = time %>% as.numeric(x = .)
                  ) %>%
                  # remove any rows that contain NA values
                  na.omit(object = .),
            position = position_jitter(width = 0.1),
            size = 1.5
          ) +
          # add a boxplot layer for each month
          # this will generate a non-problematic warning for each metabolite that 1 row was removed because
          # there is no January data currently
          geom_boxplot(
            aes(
              y = Y,
              group = time
              ),
            data = MetaboliteDataToUse %>%
                    select(
                      .data = .,
                      "Y" = all_of(x = metaboliteOfInterest),
                      "time" = StudyMonth
                    ) %>%
                    # ensure that the Y and time are numerics
                    mutate(
                      .data = .,
                      Y = Y %>% as.numeric(x = .),
                      time = time %>% as.numeric(x = .)
                    ) %>%
                    # remove any rows that contain NA values
                    na.omit(object = .),
            # make the boxes in the boxplot transparent
            alpha = 0.1,
            # remove the outliers on the boxplot since they are already present in the data
            outlier.shape = NA
          ) +
          annotate(
            geom = "text",
            family = "Arial",
            x = 2,
            y = 98,
            label = paste(
                          "R^2 =",
                          round(
                                x = R.squared,
                                digits = 3
                                )
                          ),
            fontface = "bold",
            size = 4
          ) +
          theme_classic() +
          scale_y_continuous(
            breaks = c(0,12,20,30,40,50,60,70,80,90,100),
            labels = c(0,12,20,30,40,50,60,70,80,90,100),
            limits = c(0,100)
          ) +
          labs(
               y = "Concentration (ng/mL)",
               x = "Calendar month"
               ) +
          theme(
            axis.title = element_text(family = "Arial",face = "bold",size = 11),
            axis.text.x = element_text(family = "Arial",face = "bold",size = 11),
            axis.text.y = element_text(family = "Arial",face = "bold",size = 11),
            legend.text = element_text(family = "Arial",face = "bold",size = 8),
            legend.title = element_text(family = "Arial",face = "bold",size = 8,hjust = 0.5),
            legend.position = c(.90, .90),
            legend.justification = c("right", "top"),
            legend.box.just = "right",
            legend.margin = margin(6, 6, 6, 6),
            legend.box.background = element_rect(color="black", size=1)
          ) +
          scale_x_continuous(
            breaks = c(1,2,3,4,5,6,7,8,9,10,11,12),
            # label each month with the name of the month and the number of observations for each month
            labels = MetaboliteDataToUse %>%
                      select(
                        .data = .,
                        "Y" = all_of(x = metaboliteOfInterest),
                        "time" = StudyMonth
                      ) %>%
                      # remove any rows that contain NA values
                      na.omit(object = .) %>%
                      mutate(
                        .data = .,
                        time = time %>% as.numeric(x = .)
                        ) %>%
                      select(.data = .,time) %>%
                      # arrange the time by month in numeric order
                      arrange(
                        .data = .,
                        time
                        ) %>%
                      # count the number of observations in each month
                      group_by(
                        .data = .,
                        time
                        ) %>%
                      summarize(
                        .data = .,
                        N = n()
                        ) %>%
                        # add a row for january with a count of 0
                        rbind(
                          data.frame(
                            "time" = c(1),
                            "N" = c(0)
                            ),
                          .
                        ) %>%
                        # add a column with the month abbreviations
                        mutate(
                          .data = .,
                          Month = c("JAN","FEB","MAR","APR","MAY","JUN","JUL","AUG","SEPT","OCT","NOV","DEC")
                          ) %>%
                        # add a column with the label which is the month and the sample count separated by a new line
                        mutate(
                          .data = .,
                          label = paste0(Month,"\n","(N = ",N,")")
                          ) %>%
                        pull(.data = .,label),
                        limits = c(1,13)
            )
        
        # if the metabolite of interest is primary_25OHD3_ngperml or overall 25OHD, update the legend title
        if(metaboliteOfInterest=="primary_25OHD3_ngperml" | metaboliteOfInterest=="Overall_primary_25OHD_ngperml")
        {
          PlotToReturn <-
          PlotToReturn +
          scale_color_discrete(name = "Clinical Status")
        }
        
        # print a message about the R^2 value for the cosinor regression model and the peak and trough prediction time
        print("***************************************************************************************************************************************************")
        print(glue("****For cosinor model with no covariates with data modeled by month, coefficient of determination (R^2) for {metaboliteOfInterest} is {R.squared} ******"))
        print(paste("Model predicted mean annual concentration is:",AnnualMeanConcentration))
        print(paste("Model predicted amplitude is:",Amplitude))
        print(paste("Model predicted mean maximum concentration is:",MeanMaximumConc))
        print(paste("Peak prediction month is:",PeakPredictionTime))
        print(paste("Model predicted mean minimum concentration is:",MeanMinimumConc))
        print(paste("Trough prediction month is:",TroughPredictionTime))
        print("***************************************************************************************************************************************************")
        PlotToReturn
        return(PlotToReturn)
      }
    
##################################################################################################################################
    
    
    ### define a function for trimming leading and trailing white space from all entries in an arbitrary dataframe
    ### as well as converting any empty strings into missing values
    RemoveWhiteSpaceAndChangeEmptyStringsToNA <-
      function(
        DataFrameToUse = NULL
      )
      {
        DataFrameToReturn <-
          DataFrameToUse %>%
        # remove any starting or trailing whitespace from each column in the table,
        # if the entry in the dataframe is an empty string, convert it to an NA
          lapply(
            X = .,
            FUN = function(currentColumnVector)
            {
              currentColumnVector %>%
                # loop through each entry in the currentColumnVector
                lapply(
                  X = .,
                  FUN = function(currentElement)
                  {
                    # trim the whitespace from the front and back of the string
                    StringToReturn <-
                      currentElement %>%
                      str_trim(
                        string = .,
                        side = "both"
                        )
      
                    # if the current element is an empty string or an NA, make it an NA
                    if((StringToReturn=="") | is.na(x = StringToReturn))
                    {
                      StringToReturn  <- as.character(x = NA)
                    }
      
                    return(StringToReturn)
                  }
                    ) %>%
                unlist(x = .)
            }
              ) %>%
          # convert the final list back to a dataframe
          as.data.frame(x = .)
        
        return(DataFrameToReturn)
      }
    
##################################################################################################################################
    
    # define a function where a user can select a metabolite of interest from the metabolite data and a cosinor model 
    # will be fit with no covariates included by the day of the year for a 365 day year
    CreateCosinorModelWithNoCovariatesByDayOfYear <-
      function(
        metaboliteOfInterest = NULL, # a string with the name of the metabolite of interest from the metabolite data (example: "primary_25OHD3_ngperml")
        MetaboliteDataToUse = NULL # metabolite data dataframe
        )
      {

        DataToFitCosinorModelTo <-
        # select the metabolite of interest and time columns, time is the day of the year in this case
        MetaboliteDataToUse %>%
          select(
            .data = .,
            "Y" = all_of(x = metaboliteOfInterest),
            "time" = DayOfYear
            ) %>%
          # ensure that the Y and time are numerics
          mutate(
            .data = .,
            Y = Y %>% as.numeric(x = .),
            time = time %>% as.numeric(x = .)
            ) %>%
          # remove any rows that contain NA values
          na.omit(object = .)

        # fit the cosinor model with the data using a period of 365 days
        CosinorModel <-
          DataToFitCosinorModelTo %>%
          cosinor.lm(
            formula = Y ~ time(time),
            data = .,
            period = 365
            )
        
        # calculate the sum of squared residuals (SSR) and the sum of squares total (SST)
        SSR <- ((CosinorModel$fit$residuals)^2) %>% sum(.,na.rm = TRUE)
        SST <- ((DataToFitCosinorModelTo$Y - mean(x = DataToFitCosinorModelTo$Y,na.rm = TRUE))^2) %>% sum(.,na.rm = TRUE)

        # calculate the coefficient of determination (R-squared) as 1-(SSR/SST)
        R.squared <- (1-(SSR/SST)) %>% round(x = .,digits = 3)

        # obtain the transformed model coeffiencts from the cosinor regression summary,
        # which includes the annual mean concentration estimate, the amplitude, and acrophase
        CosinorModelSummary <-
          CosinorModel %>%
          summary(object = .) %>%
          .$transformed.table %>%
          # convert rownames to a column
          mutate(
            .data = .,
            Term = rownames(x = .)
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

        # identify the cosinor model predicted annual mean concentration
        AnnualMeanConcentration <-
          CosinorModelSummary %>%
          filter(
            .data = .,
            Term == "(Intercept)"
            ) %>%
          pull(
            .data = .,
            estimate_standard.error
            )

        # identify the cosinor model predicted amplitude
        Amplitude <-
          CosinorModelSummary %>%
          filter(
            .data = .,
            Term == "amp"
          ) %>%
          pull(
            .data = .,
            estimate_standard.error
          )

        # identify the mean minimum and maximum values from the cosinor model fitted values
        MeanMinimumConc <-
          CosinorModel$fit$fitted.values %>% min()

        MeanMaximumConc <-
          CosinorModel$fit$fitted.values %>% max()
        
        # if the metabolite being plotted is the "primary_25OHD3_ngperml" or "Overall_primary_25OHD_ngperml" metabolite,
        # conditionally color the points of the graph based on vitamin D toxicity, sufficiency, insufficiency, and deficiency
        if(metaboliteOfInterest=="primary_25OHD3_ngperml" | metaboliteOfInterest=="Overall_primary_25OHD_ngperml")
        {
          VitaminDstatusCategories <-
            MetaboliteDataToUse %>%
            select(
              .data = .,
              "Y" = all_of(x = metaboliteOfInterest),
              "time" = DayOfYear
            ) %>%
            # ensure that the Y and time are numerics
            mutate(
              .data = .,
              Y = Y %>% as.numeric(x = .),
              time = time %>% as.numeric(x = .)
            ) %>%
            # remove any rows that contain NA values
            na.omit(object = .) %>%
            # create a column with vitamin D status categories
            mutate(
              .data = .,
              primaryDstatus = if_else(
                                        # if the value is less than 12, code as deficient
                                        condition = Y < 12,
                                        true = "Deficient (< 12 ng/mL)",
                                        false = if_else(
                                                        # if the value is between 12 and 20, code as insufficient
                                                        condition = Y >= 12 & Y <= 20,
                                                        true = "Insufficient (between 12 & 20 ng/mL)",
                                                        false = if_else(
                                                                        # if the value is greater than 50, code as high
                                                                        condition = Y > 50,
                                                                        true = "High (> 50 ng/mL)",
                                                                        # otherwise code as sufficient
                                                                        false = "Sufficient (> 20 & ≤ 50 ng/mL)"
                                                                      )
                                                          )
                                            )
            ) %>%
            # pull out the primaryDstatus to an array
            pull(.data = .,primaryDstatus)
        }


         # Plot the cosinor model overlayed with the original data to see how well the model fits the original data
         # I have modified Michael Sach's original plotting function code in the cosinor package to allow for this
         PlotToReturn <-
           CosinorModel %>%
           CosinorPlottingFunction_ModifiedCodeFromMichaelSachs(
             object = .,
             x_str = NULL,
             originalTimePointArray = DataToFitCosinorModelTo$time,
             scatterPlotPointTransparency = 0
             ) +
           # add a layer of points for the original data
           geom_point(
             mapping = aes(
               x = time,
               y = Y,
               # if the metabolite of interest is the primary D metabolite,
               # color the graphed points by clinical D status defined above
               # otherwise do not add color groups to the points
               color = if(metaboliteOfInterest=="primary_25OHD3_ngperml" | metaboliteOfInterest=="Overall_primary_25OHD_ngperml")
               {
                 factor(
                   x = VitaminDstatusCategories,
                   levels = c(
                     "High (> 50 ng/mL)",
                     "Sufficient (> 20 & ≤ 50 ng/mL)",
                     "Insufficient (between 12 & 20 ng/mL)",
                     "Deficient (< 12 ng/mL)"
                   )
                 )
               } else {
                 NULL
               }
             ),
             data = MetaboliteDataToUse %>%
               select(
                 .data = .,
                 "Y" = all_of(x = metaboliteOfInterest),
                 "time" = DayOfYear
               ) %>%
               # ensure that the Y and time are numerics
               mutate(
                 .data = .,
                 Y = Y %>% as.numeric(x = .),
                 time = time %>% as.numeric(x = .)
               ) %>%
               # remove any rows that contain NA values
              na.omit(object = .),
             size = 1.5
           ) +
           annotate(
             geom = "text",
             x = 5,
             y = 98,
             label = paste(
                           "R^2 =",
                           round(x = R.squared,digits = 3)
                           ),
             fontface = "bold",
             size = 4,
             family = "Arial"
           ) +
           theme_classic() +
           labs(
             y = "Concentration (ng/mL)",
             x = "Calendar date (month/day)"
             ) +
           theme(
             axis.title = element_text(family = "Arial",face = "bold",size = 11),
             axis.text.x = element_text(family = "Arial",face = "bold",size = 11),
             axis.text.y = element_text(family = "Arial",face = "bold",size = 11),
             legend.text = element_text(family = "Arial",face = "bold",size = 8),
             legend.title = element_text(family = "Arial",face = "bold",size = 8,hjust = 0.5),
             legend.position = c(.90, .90),
             legend.justification = c("right", "top"),
             legend.box.just = "right",
             legend.margin = margin(6, 6, 6, 6),
             legend.box.background = element_rect(color="black", size=1)
             ) +
           scale_x_continuous(
             breaks = c(1,31,28,31,30,31,30,31,31,30,31,30) %>% cumsum(x = .) %>% c(.,365),
             labels = c("1/1","2/1","3/1","4/1","5/1","6/1","7/1","8/1","9/1","10/1","11/1","12/1","12/31"),
             limits = c(0,365)
               ) +
           scale_y_continuous(
             breaks = c(0,12,20,30,40,50,60,70,80,90,100),
             labels = c(0,12,20,30,40,50,60,70,80,90,100),
             limits = c(0,100)
               )

         # if the metabolite of interest is primary_25OHD3_ngperml, update the legend title
         if(metaboliteOfInterest=="primary_25OHD3_ngperml" | metaboliteOfInterest=="Overall_primary_25OHD_ngperml")
         {
           PlotToReturn <-
           PlotToReturn +
           scale_color_discrete(name = "Clinical Status")
         }
         
         # print a message about the R^2 value for the cosinor regression model and the peak and trough prediction time
         print("***************************************************************************************************************************************************")
         print(glue("****For cosinor model with no covariates with data modeled by day of 365 day year, coefficient of determination (R^2) for {metaboliteOfInterest} is {R.squared} ******"))
         print(paste("Model predicted mean annual concentration is:",AnnualMeanConcentration))
         print(paste("Model predicted amplitude is:",Amplitude))
         print(paste("Model predicted mean maximum concentration is:",MeanMaximumConc))
         print(paste("Peak prediction day is:",PeakPredictionTime))
         print(paste("Model predicted mean minimum concentration is:",MeanMinimumConc))
         print(paste("Trough prediction day is:",TroughPredictionTime))
         print("***************************************************************************************************************************************************")
         PlotToReturn
         
         return(PlotToReturn)
      }
    
##################################################################################################################################
    
##################################################################################################################################
    
    ### create a function for summarizing demographic data and any metabolite data in a table
    # by computing means and standard deviations of continuous data and counts of categorical data
    SummarizeDemographicData <-
      function(
        DataFrameToUse = NULL # a dataframe with column names: Age.on.Study.Date, Blood_Quanta, BMI, Sex, Overall_primary_25OHD_ngperml, and primary_25OHD3_ngperml
      ) {
          
          DataFrameToReturn <-
            DataFrameToUse %>%
            # select relevant seasonal and demographic data as well as vitamin D and vitamin D metabolite concentrations
            select(
              .data = .,
              Age.on.Study.Date,
              Blood_Quanta,
              BMI,
              Sex,
              Overall_primary_25OHD_ngperml,
              primary_25OHD3_ngperml
            ) %>%
            # convert continuous data to numeric arrays
            mutate_at(
              .tbl = .,
              .vars = c(
                        "Age.on.Study.Date",
                        "Blood_Quanta",
                        "BMI",
                        "Overall_primary_25OHD_ngperml",
                        "primary_25OHD3_ngperml"
                        ),
              .funs = function(currentColumn)
                      {
                        DataToReturn <-
                          currentColumn %>%
                          as.numeric(x = .)

                        return(DataToReturn)
                      }
              ) %>%
            # compute sample counts, standard deviations and means
            summarise(
              .data = .,
              "N" = n(),
              "Age_mean" = Age.on.Study.Date %>% mean(x = .,na.rm = TRUE) %>% round(x = .,digits = 0),
              "Age_stdev" = Age.on.Study.Date %>% sd(x = .,na.rm = TRUE) %>% round(x = .,digits = 0),
              "Blood_Quanta_mean" = Blood_Quanta %>% mean(x = .,na.rm = TRUE) %>% round(x = .,digits = 3),
              "Blood_Quanta_max" = Blood_Quanta %>% max(.,na.rm = TRUE),
              "Blood_Quanta_min" = Blood_Quanta %>% min(.,na.rm = TRUE),
              "Blood_Quanta_stdev" = Blood_Quanta %>% sd(x = .,na.rm = TRUE) %>% round(x = .,digits = 3),
              "BMI_mean" = BMI %>% mean(x = .,na.rm = TRUE) %>% round(x = .,digits = 0),
              "BMI_stdev" = BMI %>% sd(x = .,na.rm = TRUE) %>% round(x = .,digits = 0),
              "Overall_primary_25OHD_mean" = Overall_primary_25OHD_ngperml %>% mean(x = .,na.rm = TRUE) %>% round(x = .,digits = 1),
              "Overall_primary_25OHD_stdev" = Overall_primary_25OHD_ngperml %>% sd(x = .,na.rm = TRUE) %>% round(x = .,digits = 1),
              "primary_25OHD3_mean" = mean(x = primary_25OHD3_ngperml,na.rm = TRUE) %>% round(x = .,digits = 1),
              "primary_25OHD3_stdev" = sd(x = primary_25OHD3_ngperml,na.rm = TRUE) %>% round(x = .,digits = 1),
              "Males_n" = length(x = which(x = Sex == "Male")),
              "Females_n" = length(x = which(x = Sex == "Female"))
              ) %>%
            # paste the mean and standard deviation columns together
            mutate(
              .data = .,
              "Age_y" = paste(Age_mean,"±",Age_stdev),
              "BloodQuant" = paste(Blood_Quanta_mean,"±",Blood_Quanta_stdev),
              "BMI_kg/m2" = paste(BMI_mean,"±",BMI_stdev),
              "Overall_25OHD_ng/ml" = paste(Overall_primary_25OHD_mean,"±",Overall_primary_25OHD_stdev),
              "Serum_25OHD3_ngperml" = paste(primary_25OHD3_mean,"±",primary_25OHD3_stdev)
              ) %>%
            # deselect unnecessary columns
            select(
              .data = .,
              -Age_mean,
              -Age_stdev,
              -Blood_Quanta_mean,
              -Blood_Quanta_stdev,
              -BMI_mean,
              -BMI_stdev,
              -Overall_primary_25OHD_mean,
              -Overall_primary_25OHD_stdev,
              -primary_25OHD3_mean,
              -primary_25OHD3_stdev
              )
          
          return(DataFrameToReturn)
      }
    

#########################################################################################################################################################################
    
    ### create a function for creating a tidy summary of the linear robust standard errors regression results
    # that are output from the R lm_robust() function from the estimatrr package

    TidyLinearRobustStandardErrorsRegressionResults <-
      function(
        DataForRegression = NULL, # the dataframe used for the linear regression with observations with missing values removed
        RegressionModel = NULL # the lm_robust() regression model object
      )
      {
          # summarize the regression model
          ModelSummary <-
            RegressionModel %>%
            summary(object = .)

          ResultsToReturn <-
            RegressionModel %>%
            # tidy the regression model with the etimatr tidy() function
            estimatr::tidy(x = .) %>%
            mutate(
              .data = .,
              # add the sample size
              N = DataForRegression %>% na.omit(object = .) %>% nrow(x = .),
              # add the unadjusted and adjusted R-squared
              R2 = ModelSummary$r.squared %>% round(x = .,digits = 3),
              R2_adj = ModelSummary$adj.r.squared %>% round(x = .,digits = 3),
              # add a column with the estimate and standard error pasted together,
              estimate_std.error = paste(
                                        signif(x = estimate,digits = 3),
                                        "±",
                                        signif(x = std.error,digits = 3)
                                        ),
              # round the p.value and add astericks for level of significance
              p.value = p.value %>% signif(x = .,digits = 3) %>% AddAstericksToPvalues(columnVector = .),
              # add the F-statistic
              F_stat = ModelSummary$fstatistic[["value"]] %>% signif(x = .,digits = 3),
              # add the F-statistic p-value by using the pf() function and add astericks for level of significance
              F_stat_p_value = pf(
                                  q = ModelSummary$fstatistic["value"],
                                  df1 = ModelSummary$fstatistic["numdf"],
                                  df2 = ModelSummary$fstatistic["dendf"],
                                  lower.tail = FALSE
                                  ) %>%
                                signif(
                                  x = .,
                                  digits = 3
                                  ) %>%
                                AddAstericksToPvalues(columnVector = .)
            )


          return(ResultsToReturn)
      }
    
#########################################################################################################################################################################

    ### create a function for creating a tidy summary of the logistic regression results that are output from the R glm(family = "binomial") function

    TidyLogisticRegressionResults <-
      function(
        DataForRegression = NULL, # the dataframe used for the logistic regression with observations with missing values removed
        RegressionModel = NULL # the glm() regression model object
      )
      {
          # summarize the regression
          RegressionSummary <-
            RegressionModel %>%
            summary(object = .)

          ResultsToReturn <-
              # convert the regression coefficients to a dataframe
              RegressionSummary$coefficients %>%
              data.frame(
                .,
                row.names = NULL
                ) %>%
              # convert the rownames to an actual column
              mutate(
                .data = .,
                term = rownames(x = RegressionSummary$coefficients)
                ) %>%
              # rearrange and rename columns
              select(
                .data = .,
                term,
                "p.value" = Pr...z..,
                everything()
                ) %>%
              mutate(
                .data = .,
                # add a column with the sample size
                N = DataForRegression %>% na.omit(object = .) %>% nrow(x = ),
                # add a column with the null deviance
                Null.deviance = RegressionSummary$null.deviance,
                # add a column with the residual deviance
                Deviance = RegressionSummary$deviance,
                # add a column with the Akaike Information Criterion (AIC)
                AIC = RegressionSummary$aic,
                # add a column with the estimate and standard error pasted together
                estimate_std.error = paste0(signif(x = Estimate,digits = 3)," ± ",signif(x = Std..Error,digits = 3)),
                # round the p-value and add astericks with levels of significance
                p.value = p.value %>% signif(x = .,digits = 3) %>% AddAstericksToPvalues(columnVector = .),
                # Compute the P-value for significance of the overall model by comparing the model with predictors to the null regression model
                # with a distributed Chi-squared test with degrees of freedom = degrees of freedom of model with predictors - degrees of freedom on NULL model.
                # See: https://stats.idre.ucla.edu/mplus/dae/logit-regression/ for more details on this.
                ChiSq_p_value = with(
                                    RegressionModel,
                                    pchisq(
                                      null.deviance - deviance,
                                      df = df.null - df.residual,
                                      lower.tail = FALSE
                                      )
                                    ) %>%
                                  signif(x = .,digits = 3) %>%
                                  AddAstericksToPvalues(columnVector = .)
               ) %>%
              # calculate the odds-ratio and 95% confidence interval by exponentiating the model coefficients and the confidence interval
                  # bind the odds-ratio and 95% confidence intervals with the regression summary by column
              cbind(
                  .,
                  cbind(
                    "OR" = coef(object = RegressionModel),
                    confint(object = RegressionModel)
                    ) %>%
                  exp(x = .) %>%
                  data.frame(
                    .,
                    row.names = NULL
                    )
                ) %>%
              # rename the 2.5 % and 97.5 % confidence interval columns
              dplyr::rename(
                .data = .,
                CI_2.5 = X2.5..,
                CI_97.5 = X97.5..
                ) %>%
              # add a column with the odds-ratio and confidence intervals pasted together
              mutate(
                .data = .,
                OR_CI = paste0(signif(x = OR,digits = 3)," (",signif(x = CI_2.5,digits = 3),"-",signif(x = CI_97.5,digits = 3),")")
                ) %>%
              # filter out the intercept term
             filter(
               .data = .,
               term != "(Intercept)"
               )

          return(ResultsToReturn)
      }
