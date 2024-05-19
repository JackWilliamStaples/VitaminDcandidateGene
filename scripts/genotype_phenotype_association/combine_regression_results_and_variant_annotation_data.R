# this is a script for combining PLINK regression association results with the variant annotation data

# load necessary packages
source("../scripts/load_R_packages.R")
# load necessary functions
source("../scripts/genotype_phenotype_association/genetic_association_analysis_functions.R")

# list the files in the genotype and phenotype association regression results directory
FilesInDirectory <-
list.files(
  path = "./regression_results_temp/"
  ) 

#### load the additive linear allele association model results for each gene from the regression output directory
FilesToLoadLinear <-
FilesInDirectory %>%
  # find indexes of the files containing the qassoc suffix
  grepl(
    pattern = "qassoc$",
    x = .
      ) %>%
  which(x = .) %>%
  # filter the files in the directory to the qassoc indexes
  FilesInDirectory[.]

### load the logistic regression association model results for each gene
FilesToLoadLogistic <-
  FilesInDirectory %>%
  grepl(
    pattern = "logistic$",
    x = .
    ) %>%
  which(x = .) %>%
  # filter to files with the logistic index
  FilesInDirectory[.]

# combine all of the linear metabolite regression results for each gene/metabolite
MetaboliteRegressionResultsLinear <-
# load all of the files into memory
FilesToLoadLinear %>%
  lapply(
    X = .,
    FUN = function(currentFileToLoad)
    {
      # load the file with no header
      LoadedFile <-
        read.table(
          file = paste0("./regression_results_temp/",currentFileToLoad),
          header = FALSE
          )

      # change the column names
      names(x = LoadedFile) <- c("CHR","SNP","BP","NMISS","BETA","SE","R2","Tstat","Pvalue")
      
      # remove the first row of the file now that the column names have been updated
      # the first row is the name of the columns with the header=FALSE option used in read.table
      # the file did not load in correctly when the header=TRUE argument was used even though there was a header
      LoadedFile <-
        LoadedFile[-1,]

      LoadedFile <-
      # create a column with the label for the gene and the metabolite
      LoadedFile %>%
        # isolate the gene name from the possible gene options
        mutate(
          .data = .,
          Gene = currentFileToLoad %>%
                  str_extract(
                    pattern = "(CYP3A4|DHCR7|CYP2R1|CYP27B1|SULT2A1|UGT1A4|GC|VDR|CASR|CUBN|CYP24A1|RXRG|RXRB|RXRA)",
                    string = .
                      ) %>%
                  # repeat the gene label for the number of rows of the table
                  rep_len(
                    x = .,
                    length.out = nrow(x = LoadedFile)
                    )
                ) %>%
        # isolate the metabolite
        mutate(
          .data = .,
          Metabolite = currentFileToLoad %>%
                        str_extract(
                          string = .,
                          pattern = "(active_1alpha25OH2D2_ngperml|active_1alpha25OH2D3_ngperml|beta_4beta25OH2D3_ngperml|Overall_primary_25OHD_ngperml|primary_25OHD2_ngperml|primary_25OHD3_ngperml|R_24R25OH2D3_ngperml|VitaminD2_ngperml|VitaminD3_ngperml|primary_25OHD3_VitaminD3_ratio|active_1alpha25OH2D3_primary_25OHD3_ratio|R_24R25OH2D3_primary_25OHD3_ratio|beta_4beta25OH2D3_primary_25OHD3_ratio|primary_25OHD2_VitaminD2_ratio|active_1alpha25OH2D2_primary_25OHD2_ratio)"
                            ) %>%
                        rep_len(
                          x = .,
                          length.out = nrow(x = LoadedFile)
                            )
          )



      # return the loaded file
        return(LoadedFile)

    }
      ) %>%
  # bind all of the files together by row
  do.call(
    what = "rbind",
    args = .
    ) %>%
  # arrange the rows by chromosome in ascending order followed by position in ascending order
  # followed by metabolite
  # followed by P-value in ascending order
  arrange(
    .data = .,
    CHR,
    BP,
    SNP,
    Metabolite,
    Pvalue
  )

# combine all of the logistic regression results for each gene/phenotype
MetaboliteRegressionResultsLogistic <-
# load all of the files into memory
FilesToLoadLogistic %>%
  lapply(
    X = .,
    FUN = function(currentFileToLoad)
    {
      
      # load the file with no header
      LoadedFile <-
        read.table(
          file = paste0("./regression_results_temp/",currentFileToLoad),
          header = FALSE
          )
      
      # change the column names
      names(x = LoadedFile) <- c("CHR","SNP","BP","A1","TEST","NMISS","OR","SE","L95","U95","STAT","Pvalue")

      # remove the first row of the file now that the column names have been updated
      # the first row is the name of the columns with the header=FALSE option used in read.table
      # the file did not load in correctly when the header=TRUE argument was used even though there was a header
      LoadedFile <-
        LoadedFile[-1,]
      
      LoadedFile <-
      # create a column with the label for the the gene and the phenotype
      LoadedFile %>%
        # isolate the gene name from the possible gene options
        mutate(
          .data = .,
          Gene = currentFileToLoad %>%
                  str_extract(
                    pattern = "(CYP3A4|DHCR7|CYP2R1|CYP27B1|SULT2A1|UGT1A4|GC|VDR|CASR|CUBN|CYP24A1|RXRG|RXRB|RXRA)",
                    string = .
                      ) %>%
                  # repeat the gene label for the number of rows of the table
                  rep_len(
                    x = .,
                    length.out = nrow(x = LoadedFile)
                    )
                ) %>%
        # isolate the vitamin D status phenotype
        mutate(
          .data = .,
          Phenotype = currentFileToLoad %>%
                       str_extract(
                        pattern = "(VitDsufficiency|VitDdeficiency)",
                        string = .
                        )
          )



      # return the loaded file
        return(LoadedFile)

    }
      ) %>%
  # bind all of the files together by row
  do.call(
    what = "rbind",
    args = .
    ) %>%
  # arrange the rows by chromosome in ascending order followed by position in ascending order
  # followed by P-value in ascending order
  arrange(
    .data = .,
    CHR,
    BP,
    SNP,
    Pvalue
  )

VariantAnnotationData <-
# load in the ANNOVAR/pharmGKB/VEP annotation data
"../annotation_output/allVariants_annovar_pharmGKB_and_VEP_data.tsv" %>%
read.delim(
  file = .,
  header = TRUE
    ) %>%
  # select relevant columns from the loaded file
  select(
    .data = .,
    Chr:VariantID_fromInputFile,
    LRT_ADME:starAllele
    ) 

# join the ANNOVAR annotation data with the linear and logistic regression results data and save to the file system
mapply(
  FUN = function(currentRegressionDataSet,currentRegressionDataSetName)
  {
    currentRegressionDataSet %>%
      # create a duplicate variant ID column for joining in the metabolite regression results
      mutate(
        .data = .,
        variant_ID_join_column = SNP
        ) %>%
      # join the metabolite regression results with the ANNOVAR data based on the variant IDs present in the regression data set
      left_join(
        x = .,
        y = VariantAnnotationData %>%
            mutate(
              .data = .,
              variant_ID_join_column = VariantID_fromInputFile
            ),
        by = "variant_ID_join_column"
          ) %>%
      # arrange the rows by chromosome in ascending order followed by position in ascending order
      # followed by variant ID
      # followed by P-value in ascending order
      arrange(
        .data = .,
        CHR,
        BP,
        SNP,
        Pvalue
      ) %>%
      # create a column with the PLINK variant ID converted to the format that the genotype matrices are in
      # the PLINK variant ID column has the format chrom:pos:ref:alt
      # the genotype matrices have the format Xchrom.pos.ref.alt_alt
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
        ) %>%
      write_tsv(
        x = .,
        file = paste0(currentRegressionDataSetName,"_regression_results_file_and_annotation_common_variants.txt")
        )
    
  },
  list(MetaboliteRegressionResultsLinear,MetaboliteRegressionResultsLogistic),
  c("linear","logistic"),
  SIMPLIFY = FALSE
)
  



print(x = "*************************************************************************************")

print(x = "!!!!!!!!REGRESSION RESULTS WITH ANNOVAR/pharmGKB/VEP ANNOTATION SAVED TO THE FILE SYSTEM!!!!!!!!!")
  