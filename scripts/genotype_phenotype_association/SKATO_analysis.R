# load necessary packages
source("../scripts/load_R_packages.R")
# load genetic association analysis functions
source(file = "../scripts/genotype_phenotype_association/genetic_association_analysis_functions.R")

# set a random number seed for reproducibility
set.seed(seed = 17)

# read in the targeted capture regions for each gene from the original sequencing
Gene_list_regions <-
  read.table(
    file = "../invoice_gene_list_BED_and_coverage_summaries_batch_1/gene_list.bed.txt",
    header = FALSE
  )

# name the columns of the gene list region file
names(x = Gene_list_regions) <- c("Chrom","Start","Stop","Gene")

# ensure that the start and stop base pair positions in the Gene_list_regions table are numerics
Gene_list_regions <-
  Gene_list_regions %>%
  mutate(
    .data = .,
    Start = Start %>% as.numeric(x = .),
    Stop = Stop %>% as.numeric(x = .)
  )

# load the file of variant IDs and allele frequencies that were computed with PLINK
read.table(
  file = "./plink.frq",
  header = TRUE
  ) %>%
  # select columns with the variant and allele frequencies only
  select(
    .data = .,
    SNP,
    MAF
    ) %>%
  # For MAFs that are greater than 0.5, subtract the frequency from 1.
  # The MAF column in the plink.frq files does NOT actually refer to the minor allele frequency but instead refers to 
  # the alternate allele frequency, which is not always the minor (less frequent) allele. Thus, frequencies above 0.5 need to be subtracted from 1
  mutate(
    .data = .,
    MAF = if_else(
                  condition = MAF > 0.5,
                  true = 1 - MAF,
                  false = MAF
                    )
    ) %>%
  # compute SNP-weights with the minor allele frequency using the Madsen and Browning (2009) method: 1/sqrt(p*(1-p))
  # where p is the minor allele frequency (MAF)
  mutate(
    .data = .,
    weight_value = 1/sqrt(MAF*(1-MAF)) 
    ) %>%
  # deselect the minor allele frequency column now that it is no longer needed
  select(
    .data = .,
    -MAF
    ) %>%
  # save the file of SNP IDs and weights to the file system without a header
  write_tsv(
    x = .,
    file = "./SKAT_snp_weights.txt",
    col_names = FALSE
      )

# read in the snp weight file for the SKAT association tests
obj.SNPWeight <-
  Read_SNP_WeightFile(
    FileName = "./SKAT_snp_weights.txt"
    )

# load all of the PLINK rare variant .bim files
RareVariantBIMfilesToLoad <-
  list.files(path = "./rareVariantAnalysisIntermediates/") %>%
  grepl(pattern = ".bim",x = .) %>%
  which(x = .) %>%
  list.files(path = "./rareVariantAnalysisIntermediates/")[.] %>%
  mclapply(
    X = .,
    FUN = function(currentFileToLoad)
    {
          # identify the phenotype based on the file name
          currentPhenotype <-
            currentFileToLoad %>%
            str_extract(
              string = .,
              pattern = "(active_1alpha25OH2D2_ngperml|active_1alpha25OH2D3_ngperml|beta_4beta25OH2D3_ngperml|Overall_primary_25OHD_ngperml|primary_25OHD2_ngperml|primary_25OHD3_ngperml|R_24R25OH2D3_ngperml|VitaminD2_ngperml|VitaminD3_ngperml|primary_25OHD3_VitaminD3_ratio|active_1alpha25OH2D3_primary_25OHD3_ratio|R_24R25OH2D3_primary_25OHD3_ratio|beta_4beta25OH2D3_primary_25OHD3_ratio|primary_25OHD2_VitaminD2_ratio|active_1alpha25OH2D2_primary_25OHD2_ratio|VitDsufficiency|VitDdeficiency)" 
            )
          
          # read in the bim file
          glue("./rareVariantAnalysisIntermediates/{currentFileToLoad}") %>%
          read.table(file = .,header = FALSE) %>%
          # select the SNPid column and the base position column
          select(
            .data = .,
            "SNP_ID" = V2,
            "BP" = V4
          ) %>%
          # ensure that the base position column is a numeric
          mutate(
            .data = .,
            BP = BP %>% as.numeric(x = .)
          ) %>%
          # create a gene column using the Gene_list_regions file above
          mutate(
            .data = .,
            # filter the rows of the Gene_list_regions file to positions that contain the current base pair position
            SetID_Gene = BP %>%
              lapply(
                X = .,
                FUN = function(currentBasePosition)
                {
                  GeneToReturn <-
                    Gene_list_regions %>%
                    filter(
                      .data = .,
                      Start <= currentBasePosition & Stop >= currentBasePosition
                    ) %>%
                    pull(
                      .data = .,
                      Gene
                    )
                  # return the specific gene that contains the SNP
                  return(GeneToReturn)
                }
              ) %>%
              # convert from an array to a list
              unlist(x = .)
          ) %>%
          # deselect the base position column and change the order of the setID and SNPid columns
          select(
            .data = .,
            SetID_Gene,
            SNP_ID,
            -BP
          ) %>%
          # save the file to the file system with no header
          write.table(
            x = .,
            file = glue("./rareVariantAnalysisIntermediates/rare_variants_for_SKAT_analysis_{currentPhenotype}.SetID"),
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE
          )
          
          
          
    },
    mc.cores = detectCores()
    )

print(x = "$$$$$$$$$$$$$ SetID files for O-SKAT analysis were saved $$$$$$$$$$$")

# disable non-problematic warning messages for when observations with missing phenotypes or covariates are removed from the SKAT analysis
options(warn = -1)

SKATresultsForAllPhenotypes <-
# run the SKAT model for each phenotype
      c(
        "active_1alpha25OH2D2_ngperml",
        "active_1alpha25OH2D3_ngperml",
        "beta_4beta25OH2D3_ngperml",
        "Overall_primary_25OHD_ngperml",
        "primary_25OHD2_ngperml",
        "primary_25OHD3_ngperml",
        "R_24R25OH2D3_ngperml",
        "VitaminD2_ngperml",
        "VitaminD3_ngperml",
        "primary_25OHD3_VitaminD3_ratio",
        "active_1alpha25OH2D3_primary_25OHD3_ratio",
        "R_24R25OH2D3_primary_25OHD3_ratio",
        "beta_4beta25OH2D3_primary_25OHD3_ratio",
        "primary_25OHD2_VitaminD2_ratio",
        "active_1alpha25OH2D2_primary_25OHD2_ratio",
        "VitDsufficiency",
        "VitDdeficiency"
        ) %>%
        lapply(
          X = .,
          FUN = function(currentPhenotype)
          {
                # define the PLINK file names that will be loaded into SKAT and the names of the output SNP set data files (SSD) files
                Bed_file_to_load <- glue("./rareVariantAnalysisIntermediates/rare_variants_for_SKAT_analysis_{currentPhenotype}.bed")
                Bim_file_to_load <- glue("./rareVariantAnalysisIntermediates/rare_variants_for_SKAT_analysis_{currentPhenotype}.bim")
                Fam_file_to_load <- glue("./rareVariantAnalysisIntermediates/rare_variants_for_SKAT_analysis_{currentPhenotype}.fam")
                SetID_file_to_load <- glue("./rareVariantAnalysisIntermediates/rare_variants_for_SKAT_analysis_{currentPhenotype}.SetID")
                SSD_file_to_make <- glue("./rareVariantAnalysisIntermediates/rare_variants_for_SKAT_analysis_{currentPhenotype}.SSD")
                SSDinfo_file_to_make <- glue("./rareVariantAnalysisIntermediates/rare_variants_for_SKAT_analysis_{currentPhenotype}.SSD.info")
                      
                # if the .SetID file was successfully created, perform the SKAT association, otherwise print a message that the SKAT association could not be performed
                if(file.exists(glue("./rareVariantAnalysisIntermediates/rare_variants_for_SKAT_analysis_{currentPhenotype}.SetID")))
                { 
                      # create SKAT snp data set files (SSD) for each phenotype
                      Generate_SSD_SetID(
                        File.Bed = Bed_file_to_load,
                        File.Bim = Bim_file_to_load,
                        File.Fam = Fam_file_to_load,
                        File.SetID = SetID_file_to_load,
                        File.SSD = SSD_file_to_make,
                        File.Info = SSDinfo_file_to_make
                          )
                      
                      # if the phenotype is binary, set the binary phenotype option to TRUE
                      if(currentPhenotype == "VitDsufficiency" | currentPhenotype == "VitDdeficiency")
                      {
                          BinaryPhenotypeOption <- TRUE
                      } else {
                          BinaryPhenotypeOption <- FALSE
                      }
                      
                      # load the covariate file and the phenotype file of continuous or binary phenotypes
                      FAM_Cov <-
                        Read_Plink_FAM_Cov(
                                           Filename = Fam_file_to_load,
                                           File_Cov = "./CovariateFile_PLINKformat.txt",
                                           Is.binary = BinaryPhenotypeOption,
                                           cov_header = TRUE
                                          )
                      
                      # create arrays of covariates
                      Age <- FAM_Cov$Age.on.Study.Date
                      BMI <- FAM_Cov$BMI
                      Sex_male <- FAM_Cov$Sex_male
                      Sex_female <- FAM_Cov$Sex_female
                      StudySeason_A_June_Aug_peak <- FAM_Cov$StudySeason_A_June_Aug_peak
                      StudySeason_B_Sept_Nov_descending <- FAM_Cov$StudySeason_B_Sept_Nov_descending
                      StudySeason_C_Dec_Feb_trough <- FAM_Cov$StudySeason_C_Dec_Feb_trough
                      StudySeason_D_Mar_May_ascending <- FAM_Cov$StudySeason_D_Mar_May_ascending
                      
                      # create an array of the phenotype values from the loaded FAM file
                      PhenotypeValues <-
                        FAM_Cov$Phenotype
                      
                        # open the SSD file
                        Rare_VariantsSSD.INFO <-
                          Open_SSD(
                            File.SSD = SSD_file_to_make,
                            File.Info = SSDinfo_file_to_make
                          )
                        
                        # if the phenotype is binary, define the NULL SKAT model with a binary outcome
                        if(currentPhenotype == "VitDsufficiency" | currentPhenotype == "VitDdeficiency")
                        {
                          # define the NULL SKAT model on the dataset with all one-hot encoded seasonal and demographic covariates included
                          # the outcome is dichotomous in this case
                          NullSKATmodelObject <-
                            SKAT_Null_Model(
                              formula = PhenotypeValues ~ Age + BMI + Sex_male + Sex_female + StudySeason_A_June_Aug_peak + StudySeason_B_Sept_Nov_descending + StudySeason_C_Dec_Feb_trough + StudySeason_D_Mar_May_ascending,
                              out_type = "D"
                            )
                          
                        } else {
                          
                            # define the NULL SKAT model on the dataset with all one-hot encoded seasonal and demographic covariates included
                            # the outcome is continuous in this case
                            NullSKATmodelObject <-
                              SKAT_Null_Model(
                                formula = PhenotypeValues ~ Age + BMI + Sex_male + Sex_female + StudySeason_A_June_Aug_peak + StudySeason_B_Sept_Nov_descending + StudySeason_C_Dec_Feb_trough + StudySeason_D_Mar_May_ascending,
                                out_type = "C"
                              )
                        }
          
                        # Run the SKATO (optimal combination of burden and SKAT test) association on the NULL model defined above
                        # using the Madsen and Browning weights for minor allele frequency.
                        # The SKAT-O model will run for each individual gene defined in the the SETid file
                        SKAT_association_results <-
                          SKAT.SSD.All(
                            SSD.INFO = Rare_VariantsSSD.INFO,
                            obj = NullSKATmodelObject,
                            method = "optimal.adj",
                            kernel = "linear.weighted",
                            impute.method = "fixed",
                            obj.SNPWeight = obj.SNPWeight
                          )
                        
                        # save the SKAT association results table to the file system after adding a column of bonferroni adjusted P-values
                        SKAT_association_results$results <-
                          SKAT_association_results$results %>%
                          # create a column with Bonferroni-adjusted P-values for the number of SNP sets (i.e., number of genes) that were tested
                          mutate(
                            .data = .,
                            P.value.Bonferroni = nrow(x = SKAT_association_results$results)*P.value
                          ) %>%
                          # round the P.value columns and the P.value.Bonferroni columns to 3 places after the decimal
                          mutate(
                            .data = .,
                            P.value = P.value %>% signif(x = .,digits = 3),
                            P.value.Bonferroni = P.value.Bonferroni %>% signif(x = .,digits = 3),
                            # calculate the -log10() of the Bonferroni adjusted P-value
                            NegativeLogP.value.Bonferroni = P.value.Bonferroni %>% log10(x = .)*(-1)
                          ) %>%
                          # create a column with an arbitrary index for coloring points by gene
                          mutate(
                            .data = .,
                            ColorIndex = seq_len(length.out = nrow(x = .)) %>% as.factor(x = .)
                          ) %>%
                          # create a column with a label for the current Phenotype
                          mutate(
                            .data = .,
                            Phenotype = currentPhenotype
                            ) %>%
                          # arrange the genes alphabetically
                          arrange(
                            .data = .,
                            SetID
                            )
                        
                        # create a directory for saving the SKAT results
                        if(!dir.exists(paths = "./SKATO_results/"))
                        {
                          dir.create(path = "./SKATO_results/")
                        }
          
                        # save the SKAT results to the SKAT results directory
                        SKAT_association_results$results %>%
                          write.table(
                            x = .,
                            file = glue("./SKATO_results/{currentPhenotype}.txt"),
                            row.names = FALSE,
                            col.names = TRUE,
                            quote = FALSE
                          )
                        
                        # plot the bonferroni adjusted P-values of the SKAT association results and save to the file system
                        SKAT_association_results_plot <-
                          SKAT_association_results$results %>%
                          # plot the results using a dot chart as shown in this article: http://www.sthda.com/english/articles/24-ggpubr-publication-ready-plots/
                          ggdotchart(
                            data = .,
                            x = "SetID",
                            y = "NegativeLogP.value.Bonferroni",
                            color = "ColorIndex",                         # Color by groups
                            rotate = TRUE,                                # Rotate vertically
                            dot.size = 4,                                 # Large dot size
                            ggtheme = theme_pubr()                        # ggplot2 theme
                          )+
                          geom_hline(
                            yintercept = -log10(x = 0.05),
                            linetype="solid",
                            color="green",
                            size=1.5
                            ) +
                          theme_cleveland() + # Add dashed grids
                          scale_fill_viridis(discrete = TRUE) +
                          rremove(object = "legend") +
                          theme(
                            axis.title.x = element_text(
                              margin=margin(5,0,0,0),
                              family = "Arial",
                              face = "bold",
                              size = 11
                            ),
                            axis.title.y = element_text(
                                                        margin=margin(0,5,0,0),
                                                        family = "Arial",
                                                        face = "bold",
                                                        size = 11
                                                      )
                          )
                        # add additional plot formatting and change names of titles
                        SKAT_association_results_plot <-
                          ggpar(
                            p = SKAT_association_results_plot,
                            xlab = "Gene",
                            font.x = "bold",
                            font.xtickslab = "bold",
                            font.family = "Arial",
                            ylab = "-log10(p-value)",
                            font.y = "bold",
                            font.ytickslab = "bold.italic",
                            yticks.by = 1
                          )
          
                        # save the final plot to the file system as a tiff
                        tiff(
                          filename = glue("./SKATO_results/{currentPhenotype}.tiff"),
                          width = 9,
                          height = 5,
                          units = "in",
                          compression = "none",
                          res = 300
                          )
                        
                        print(x = SKAT_association_results_plot)
                        dev.off()
                        
                        return(SKAT_association_results$results)
                        
            } else {
                  
                        return(NULL)
              
                        print(x = "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
                        print(x = "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
                        
                        print(x = glue("The ./rareVariantAnalysisIntermediates/rare_variants_for_SKAT_analysis_{currentPhenotype}.SetID does not exist"))
                        
                        print(x = "SKAT analysis could not be performed")
                        
                        print(x = "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
                        print(x = "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
              
              
            }
        
    }
      ) %>%
  # bind the SKAT association results together by row
  do.call(
    what = "rbind",
    args = .
    ) %>%
  # group by the Phenotype and split into a list of dataframes
  group_by(
    .data = .,
    Phenotype
    ) %>%
  group_split(.tbl = .) %>%
  # Loop through the list of dataframes and append the current Phenotype to all column names.
  # This prevents duplicate columns from appearing when joining all SKAT results tables together by gene
  lapply(
    X = .,
    FUN = function(currentDataFrame)
    {
      # extract the names of the current dataframe
      dataframeNames <-
        currentDataFrame %>%
        names(x = .)
      
      # add the Phenotype to the column name for all columns except for the SetID column
      ModifiedColumnNames <-
            dataframeNames %>%
              lapply(
                X = .,
                FUN = function(currentColumnName)
                {
                    if(currentColumnName!="SetID")
                    {
                      ValueToReturn <- paste0(currentColumnName,"_",unique(x = currentDataFrame$Phenotype))
                    } else {
                      ValueToReturn <- currentColumnName
                    }
                    return(ValueToReturn)
                }
                  ) %>%
               unlist(x = .)
      
       # rename the columns of the current dataframe with the dataframeNames
       names(x = currentDataFrame) <- ModifiedColumnNames
       
       return(currentDataFrame)
    }
      ) %>%
  # join the list of dataframes together by gene (SetID)
  purrr::reduce(
    .x = .,
    .f = inner_join,
    by = "SetID"
    ) %>%
  # select the gene column and the columns with bonferroni adjusted p-values
  select(
    .data = .,
    "Gene" = SetID,
    starts_with(match = "P.value")
    ) %>%
  # add astericks to P-value columns to indicate levels of significance
  mutate_at(
    .tbl = .,
    .vars = vars(starts_with(match = "P.value.Bonferroni")),
    .funs = ~ AddAstericksToPvalues(columnVector = .x)
      ) %>%
  # arrange the results alphabetically by gene
  arrange(
    .data = .,
    Gene
    )

options(warn = 0)

# save the SKATresultsForAllPhenotypes to the file system
SKATresultsForAllPhenotypes %>%
  write_tsv(
    x = .,
    file = "./SKATO_results/SKATresultsForAllPhenotypes.tsv"
      )




