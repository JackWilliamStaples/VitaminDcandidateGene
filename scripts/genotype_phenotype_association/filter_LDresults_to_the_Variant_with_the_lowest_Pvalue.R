# load necessary packages
source("../scripts/load_R_packages.R")

FilesInDirectory <-
  # list the files in the locus zoom temporary directory
  list.files(
    path = "./locusZoom_temp/"
  ) 

# identify the linkage disequilibrium files that need to be loaded based on the "linkage_disequilibrium.txt" string
FilesToLoad <-  
  FilesInDirectory %>%
  # identify files that end with the pattern "linkage_disequilibrium.txt"
  grepl(
    pattern = "linkage_disequilibrium.txt",
    x = .
  ) %>%
  which(x = .) %>%
  FilesInDirectory[.]

# load the linkage disequilibrium data that was computed in the previous step
LinkageDisequlibriumData <-
  FilesToLoad %>%
  lapply(
    X = .,
    FUN = function(currentFile)
    {
      
      LoadedFile <-
        glue("./locusZoom_temp/{currentFile}") %>%
        read.table(file = .,header = TRUE) %>%
        # create a column for the Phenotype
        mutate(
          .data = .,
          Phenotype = 
            str_extract(
              string = currentFile,
              pattern = "(active_1alpha25OH2D2_ngperml|active_1alpha25OH2D3_ngperml|beta_4beta25OH2D3_ngperml|Overall_primary_25OHD_ngperml|primary_25OHD2_ngperml|primary_25OHD3_ngperml|R_24R25OH2D3_ngperml|VitaminD2_ngperml|VitaminD3_ngperml|primary_25OHD3_VitaminD3_ratio|active_1alpha25OH2D3_primary_25OHD3_ratio|R_24R25OH2D3_primary_25OHD3_ratio|beta_4beta25OH2D3_primary_25OHD3_ratio|primary_25OHD2_VitaminD2_ratio|active_1alpha25OH2D2_primary_25OHD2_ratio|VitDsufficiency|VitDdeficiency)"
            )
        ) %>%
        # create a column for the gene
        mutate(
          .data = .,
          Gene = str_extract(
            string = currentFile,
            pattern = "(CYP3A4|DHCR7|CYP2R1|CYP27B1|SULT2A1|UGT1A4|GC|VDR|CASR|CUBN|CYP24A1|RXRG|RXRB|RXRA)"
          )
        )
      
      return(LoadedFile)
    }
  ) %>%
  # combine all the results into one dataframe
  do.call(
    what = "rbind",
    args = .
  )

# load the variants and their candidate-gene association p values that will be plotted with LocusZoom
FilesOfVariantsAndPvalues <-
  list.files(path = "./Results_for_LocusZoom/")

# load the locusZoom data that will be plotted for each gene and phenotype
VariantsAndPvalues <-
  FilesOfVariantsAndPvalues %>%
  mclapply(
    X = .,
    FUN = function(currentFile)
    {
        FileToReturn <-
        # load the current file
              glue("./Results_for_LocusZoom/{currentFile}") %>%
              read.table(
                file = .,
                header = TRUE
              ) %>%
              # add a column for the Phenotype
              mutate(
                .data = .,
                Phenotype = 
                  str_extract(
                    string = currentFile,
                    pattern = "(active_1alpha25OH2D2_ngperml|active_1alpha25OH2D3_ngperml|beta_4beta25OH2D3_ngperml|Overall_primary_25OHD_ngperml|primary_25OHD2_ngperml|primary_25OHD3_ngperml|R_24R25OH2D3_ngperml|VitaminD2_ngperml|VitaminD3_ngperml|primary_25OHD3_VitaminD3_ratio|active_1alpha25OH2D3_primary_25OHD3_ratio|R_24R25OH2D3_primary_25OHD3_ratio|beta_4beta25OH2D3_primary_25OHD3_ratio|primary_25OHD2_VitaminD2_ratio|active_1alpha25OH2D2_primary_25OHD2_ratio|VitDsufficiency|VitDdeficiency)"
                  )
              ) %>%
              # add a column for the gene
              mutate(
                .data = .,
                Gene = 
                  str_extract(
                    string = currentFile,
                    pattern = "(CYP3A4|DHCR7|CYP2R1|CYP27B1|SULT2A1|UGT1A4|GC|VDR|CASR|CUBN|CYP24A1|RXRG|RXRB|RXRA)"
                  )
              ) %>%
              # add a column with the nucleotide position extracted from the variant ID
              mutate(
                .data = .,
                NucleotidePostion = MarkerName %>%
                  # isolate the nucleotide position from the variant ID
                  # by first removing the chromosome number, then by removing the nucleotides
                  sub(
                    pattern = "[[:digit:]]+[[:punct:]]+",
                    replacement = "",
                    x = .
                  ) %>%
                  sub(
                    pattern = "[[:punct:]]+[[:alpha:]]+[[:punct:]]+([[:alpha:]]+|[[:punct:]]+)",
                    replacement = "",
                    x = .
                    )
                ) 
      
      
      
      return(FileToReturn)
    },
    mc.cores = detectCores()
  )

# name the list of dataframes of variants and P-values with the original file name they came from
names(x = VariantsAndPvalues) <- FilesOfVariantsAndPvalues

# save tables to the file system with R^2 linkage disequilibrium values for all SNPs that are computed with reference to the 
# reference SNP with the minimum P-value for the gene for cases where PLINK was able to calculate linkage disequilibrium 
# with reference to the lowest possible P-value SNP for the gene for each gene/phenotype pair
 VariantsAndPvalues %>%
    unname(obj = .) %>%
    lapply(
      X = .,
      FUN = function(currentDataSetOfVariantsAndPvalues)
      {
        # join the current data set of variants and P-values with the linkage disequilibrium data to obtain linkage disequilibrium (R^2) values
        # for the variants to be plotted with reference to the SNP with the lowest P-value for the gene
        LDresultsWithReferenceToMostSignificantVariant <-
            currentDataSetOfVariantsAndPvalues %>%
              left_join(
                x = .,
                y = LinkageDisequlibriumData %>%
                    # filter the LinkageDisequlibriumData to the current Phenotype and gene in the currentDataSetOfVariantsAndPvalues
                    filter(
                      .data = .,
                      Phenotype == unique(x = currentDataSetOfVariantsAndPvalues$Phenotype) & 
                      Gene == unique(x = currentDataSetOfVariantsAndPvalues$Gene)
                      ),
                by = c("MarkerName" = "SNP_A")
              ) %>%
              # remove any results where linkage disequilibrium could not be computed (i.e., LD R^2 is NA)
              na.omit(object = .) %>%
              # remove any rows where the MarkerName is the same as SNP_B (linkage disequilibrium between the SNP and itself is not necessary)
              filter(
                .data = .,
                MarkerName != SNP_B
              ) %>%
              # group by the Phenotype, then by the gene
              group_by(
                .data = .,
                Phenotype.x,
                Gene.x
              ) %>%
              # filter the results for each Phenotype/gene to the lowest P-value variant only
              # resulting in linkage disequilibrium between all variants and the reference variant with the lowest P-value
              filter(
                .data = .,
                P.value == min(P.value,na.rm = TRUE)
              ) %>%
              # split the dataframe into a list of dataframes based on the Phenotype/gene grouping
              group_split(
                .tbl = .,
                .keep = TRUE
              ) %>%
              # loop through what is now a list of dataframes
              lapply(
                X = .,
                FUN = function(currentDataFrame)
                {
                  # save the name of the current Phenotype and gene
                  Phenotype <- currentDataFrame$Phenotype.x %>% unique(x = .)
                  Gene <- currentDataFrame$Gene.x %>% unique(x = .)
        
        
                  # update the columns to match the format required by LocusZoom
                  currentDataFrame %>%
                    # ensure that the dataframe is no longer grouped
                    ungroup(x = .) %>%
                    # filter to rows that match the current Phenotype and current gene
                    filter(
                      .data = .,
                      (Phenotype.x == Phenotype.y) & (Gene.x == Gene.y)
                    ) %>%
                    # create a dprime column filled with missing values, only r2 is required for locusZoom linkage disequilibrium
                    mutate(
                      .data = .,
                      # fill the column with NAs,
                      # column length is the same as any other column in the dataframe
                      dprime = rep_len(
                        x = as.numeric(x = NA),
                        length.out = length(x = R2)
                      )
                    ) %>%
                    select(
                      .data = .,
                      # snp1 is Any SNP in your plotting region.
                      "snp1" = SNP_B,
                      # "snp2 should always be the reference SNP in the region"
                      "snp2" = MarkerName,
                      # D' between snp2 (reference SNP) and snp1.(this can be missing if not known)
                      dprime,
                      # r2 between snp2 (reference SNP) and snp1.
                      "rsquare" = R2
                    )  %>%
                    # create a column for the Phenotype and gene
                    mutate(
                      .data = .,
                      Phenotype = Phenotype
                      ) %>%
                    mutate(
                      .data = .,
                      Gene = Gene
                      )
                }
              )
        
        # save the final file with linkage disequilibrium between the reference SNP with the lowest P-value and all SNPs for the current gene/phenotype
        LDresultsWithReferenceToMostSignificantVariant %>%
        lapply(
          X = .,
          FUN = function(currentSetOfLDresults)
          {
                # obtain the current Phenotype and gene for saving the final output file
                Phenotype <- currentSetOfLDresults$Phenotype %>% unique(x = .)
                Gene <- currentSetOfLDresults$Gene %>% unique(x = .)
                
                FileToSave <-
                  currentSetOfLDresults %>%
                  # select the relevant columns
                  select(
                    .data = .,
                    snp1,
                    snp2,
                    dprime,
                    rsquare
                    ) %>%
                  # remove the nucleotides appended to snp1 and snp2 (locusZoom format is chr:position without nucleotides)
                  # the final character can be a nucleotide or a "*" symbol for deletions
                  mutate(
                    .data = .,
                    snp1 = snp1 %>%
                           sub(
                             pattern = "[[:punct:]]+[[:alpha:]]+[[:punct:]]+([[:alpha:]]+|[[:punct:]]+)",
                             replacement = "",
                             x = .
                             ),
                    snp2 = snp2 %>%
                            sub(
                              pattern = "[[:punct:]]+[[:alpha:]]+[[:punct:]]+([[:alpha:]]+|[[:punct:]]+)",
                              replacement = "",
                              x = .
                            )
                    ) %>%
                  # ensure the rows are unique
                  unique(x = .)
                
                  # if the data frame being saved has at least one row, save it to the file system,
                  # otherwise save nothing
                  if(nrow(x = FileToSave)>0)
                  {
                    print(x = glue("linkageDisequilibriumAndPvalues_forLocusZoom_{Phenotype}_{Gene}.txt is about to be saved"))
                  
                    
                    # save whitespace delimited results to the file system labeled with the Phenotype and gene
                    FileToSave %>%
                      write.table(
                        x = .,
                        file = glue("./locusZoom_temp/linkageDisequilibriumAndPvalues_forLocusZoom_{Phenotype}_{Gene}.txt"),
                        col.names = TRUE,
                        quote = FALSE,
                        row.names = FALSE
                      )
                    
                      print(x = glue("linkageDisequilibriumAndPvalues_forLocusZoom_{Phenotype}_{Gene}.txt saved to file system"))
                  } else {
                    
                      print(x = "~~~~~Nothing to be saved~~~~")
        
                  } 
            
           }
         )
      }
    )
    
    
    # save the variants and P-values files back to their original directory with their original file name only this time with the
    # nucleotides removed from the marker name column
    mapply(
      FUN = function(currentFile,currentFileName)
        {
              currentFile %>%
              select(
                .data = .,
                MarkerName,
                "P-value" = P.value
                ) %>%
              mutate(
                .data = .,
                MarkerName = MarkerName %>%
                              sub(
                                pattern = "[[:punct:]]+[[:alpha:]]+[[:punct:]]+([[:alpha:]]+|[[:punct:]]+)",
                                replacement = "",
                                x = .
                              )
                ) %>% 
              write_tsv(
                x = .,
                file = glue("./Results_for_LocusZoom/{currentFileName}")
              )
        },
      VariantsAndPvalues,
      names(x = VariantsAndPvalues),
      SIMPLIFY = FALSE
        )

