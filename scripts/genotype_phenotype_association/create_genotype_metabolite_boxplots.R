
# load necessary packages
source("../scripts/load_R_packages.R")
# load functions 
source("../scripts/genotype_phenotype_association/genetic_association_analysis_functions.R")

PLINKvariantID_dbSNP_rsID_lookupTable <-
# load the the lookup table of PLINK format variant IDs and dbSNP rsIDs for common variants
read.delim(
  file = "./PLINKvariantID_dbSNP_rsID_lookupTable_commonVariants.txt",
  header = TRUE
    )

FilesInDirectory <-
# list the files in the genotype matrices directory
list.files(path = "./genotype_matrices_with_metabolite_data/")

# load each genotype matrix in the directory
GenotypeMatricesWithMetaboliteDataAdded <-
  FilesInDirectory %>%
  mclapply(
    X = .,
    FUN = function(currentFileToLoad)
    {
      LoadedFile <-
        # load the current file
        read.delim(
          file = paste0("./genotype_matrices_with_metabolite_data/",currentFileToLoad),
          header = TRUE
        )
      
      # only select the columns that contain a significant variant as well as metabolite columns,
      # based on the PLINKvariantID_dbSNP_rsID_lookupTable, which only contains significant variants
      
      # identify column names that contain genotypes, indicated with an "X"
      ColumnNames <-
      LoadedFile %>%
        names(x = .)
      
      VariantColumns <-
        ColumnNames %>%
        grepl(
          pattern = "X",
          x = .
          ) %>%
        which(x = .) %>%
        ColumnNames[.]
      
      # identify variant columns that are in the PLINKvariantID_dbSNP_rsID_lookupTable
      VariantColumnsToSelect <-
        (VariantColumns %in% PLINKvariantID_dbSNP_rsID_lookupTable$VariantID_GenotypeMatrixFormat) %>%
        which(x = .) %>%
        VariantColumns[.]
      
      # create an array of columns to select and select those columns from the LoadedFile and return the LoadedFile
      ColumnsToSelect <-
        # combine the non-variant columns and the VariantColumnsToSelect
      ColumnNames %>%
        grepl(
          pattern = "X",
          x = .
          ) %>%
        which(x = {.} == FALSE) %>%
        ColumnNames[.] %>%
        c(
          VariantColumnsToSelect,
          .
          )
      
      FileToReturn <-
        ColumnsToSelect %>%
        LoadedFile[.]
      
      return(FileToReturn)
    },mc.cores = detectCores()
  )

# name each genotype matrix with metabolite data added according to its file names
names(x = GenotypeMatricesWithMetaboliteDataAdded) <- FilesInDirectory

# load the ordinary least squares and robust standard errors genotype regression summary statistics for annotating genotype boxplots
RegressionResultsFile_OLS_RSE <-
  read.delim(
    file = "./CausativeVariantAssociations/RegressionResultsFile_OLS_and_RSE.tsv",
    header = TRUE
    ) %>%
  # filter to the results that were significant (astericks on P-value) in robust standard errors regression after Bonferroni P-value adjustment
  filter(
    .data = .,
    grepl(
      pattern = "\\*",
      x = F_stat_p_value_adjusted_RSE
      )
    )

# temporarily suppress non-problematic warning messages for cases where missing metabolite values are removed in the boxplots
options(warn = -1) 

# make boxplots of metabolite levels for each variant with genotypes as 0, 1, or 2 copies of the alternate allele for each significant variant/metabolite pair
# Note: these plots take up a ton of space in R studio memory so only assign to an object when needed
mapply(
  FUN = function(currentGenoMatrix,nameOfcurrentGenoMatrix)
  {
    print(x = nameOfcurrentGenoMatrix)

    # obtain an array of all of the variant columns which all start with an "X" followed by the variant ID
    VariantsToPlot <-
      names(x = currentGenoMatrix) %>%
      grepl(pattern = "X",x = .) %>%
      names(x = currentGenoMatrix)[.]

    # loop through each variant and each metabolite and create a boxplot of the metabolite data grouped by genotype (0,1,2)
    VariantsToPlot %>%
      mclapply(
        X = .,
        FUN = function(currentVariant)
        {
          print(currentVariant)
          
          # identify the metabolites that need to be plotted for the current variant based on the significant robust standard errors 
          # regression results in the RegressionResultsFile_OLS_RSE 
          MetabolitesToPlot <-
            RegressionResultsFile_OLS_RSE %>% 
            filter(
              .data = .,
              VariantID_GenotypeMatrixFormat == currentVariant
              ) %>%
            pull(
              .data = .,
              metabolite
              )
          
          # loop through each metabolite for the current variant
          MetabolitesToPlot %>%
            lapply(
              X = .,
              FUN = function(currentMetabolite)
              {
                print(currentMetabolite)

                
                # obtain a variant ID column to label the plot with
                VariantID <-
                # filter the PLINKvariantID_dbSNP_rsID_lookupTable to the currentVariant
                  # based on the VariantID_GenotypeMatrixFormat column
                PLINKvariantID_dbSNP_rsID_lookupTable %>%
                  filter(
                    .data = .,
                    VariantID_GenotypeMatrixFormat == currentVariant
                    ) %>%
                  # create a variant ID column to label the plot with
                  # if the variant has a dbSNP rsID, use the rsID as the label
                  # if there is no rsID, use the SNP column for the rsID
                  mutate(
                    .data = .,
                    VariantIDtoLabelWith = if_else(
                                                   # if there is NO rsID available for the variant, use the SNP column
                                                   condition = (is.na(x = avsnp138)==TRUE) & (is.na(x = avsnp142)==TRUE) & (is.na(x = existing_variant_VEP)==TRUE),
                                                   true = SNP,
                                                   # if there is an rsID available for the variant from all dbSNP columns use the rs number from the VEP annotation data
                                                   false = if_else(
                                                                    condition = (is.na(x = avsnp138)==FALSE) & (is.na(x = avsnp142)==FALSE) & (is.na(x = existing_variant_VEP)==FALSE),
                                                                    true = existing_variant_VEP,
                                                                    # otherwise, use the existing_variant_VEP column, VEP gives a better rs number annotation compared to annovar
                                                                    false = existing_variant_VEP
                                                                      )
                                                )
                    ) %>%
                  pull(
                    .data = .,
                    VariantIDtoLabelWith
                    )

                print(paste0("Variant ID is........",VariantID))


                # identify the reference allele for the current variant
                ReferenceAllele <-
                  PLINKvariantID_dbSNP_rsID_lookupTable %>%
                  filter(
                    .data = .,
                    VariantID_GenotypeMatrixFormat == currentVariant
                  ) %>%
                  pull(.data = .,SNP) %>%
                  str_split(string = .,pattern = ":") %>%
                  unlist(x = .) %>%
                  .[length(x = .)-1]

                # identify the alternate allele for the current variant
                AlternateAllele <-
                  PLINKvariantID_dbSNP_rsID_lookupTable %>%
                  filter(
                    .data = .,
                    VariantID_GenotypeMatrixFormat == currentVariant
                  ) %>%
                  pull(.data = .,SNP) %>%
                  str_split(string = .,pattern = ":") %>%
                  unlist(x = .) %>%
                  tail(x = .,n = 1)
                
                # filter the regression results file to the current variant and metabolite
                SNPforFilteringRegressionResults <-
                PLINKvariantID_dbSNP_rsID_lookupTable %>%
                  filter(
                    .data = .,
                    VariantID_GenotypeMatrixFormat == currentVariant
                  ) %>%
                  pull(.data = .,SNP)
                
                SummaryStatisticsForLabelingPlot <-
                  RegressionResultsFile_OLS_RSE %>%
                    filter(
                      .data = .,
                      SNP == SNPforFilteringRegressionResults & metabolite == currentMetabolite
                      ) %>%
                    # select the beta-coefficient, standard error, N, and Bonferroni adjusted P-value
                    # from the robust standard errors (RSE) regression 
                    select(
                      .data = .,
                      estimate_std.error_RSE,
                      F_stat_p_value_adjusted_RSE,
                      N_RSE
                      ) %>%
                  # create labels with the beta, SE, Bonferroni adjusted P-value, and N
                  mutate(
                    .data = .,
                    BETA_label = paste0("b = ",estimate_std.error_RSE),
                    Pvalue_label = paste0("P = ",F_stat_p_value_adjusted_RSE),
                    N_label = paste0("N = ",N_RSE)
                    ) %>%
                  # deselect columns that are no longer needed
                  select(
                    .data = .,
                    -estimate_std.error_RSE,
                    -F_stat_p_value_adjusted_RSE,
                    -N_RSE
                    )

                # create a table of genotype codes (0,1,2) and their corresponding nucleotide genotypes
                GenotypeLabelTable <-
                  data.frame(
                             "Genotype_Numeric" = c(0,1,2),
                             "Genotype_alleles" = c(
                                                    paste0(ReferenceAllele,"/",ReferenceAllele),
                                                    paste0(ReferenceAllele,"/",AlternateAllele),
                                                    paste0(AlternateAllele,"/",AlternateAllele)
                                                    )
                            )
                
                DataToPlot <-
                # select the columns from the genotype matrix joined with the metabolite data for the current variant and metabolite
                currentGenoMatrix %>%
                  select(
                    .data = .,
                    "Variant" = all_of(currentVariant),
                    "Metabolite" = all_of(currentMetabolite)
                    ) %>%
                  # remove observations with missing values
                  na.omit(
                    object = .
                  )
                
                PlotToReturn <-
                DataToPlot %>%
                    ggplot(
                      data = .,
                      mapping = aes(
                                    x = factor(x = Variant),
                                    y = Metabolite
                                   )
                      ) +
                    geom_boxplot(
                      lwd = 1.5,
                      fatten = 1,
                      outlier.shape = NA
                      ) +
                  scale_y_log10(
                    limits = c(0.0001,180),
                    breaks = c(0.0001,0.001,0.01,0.1,1,5,15,30,75,160),
                    labels = c(0.0001,0.001,0.01,0.1,1,5,15,30,75,160)
                    ) +
                    geom_point(
                      shape = 21,
                      colour = "black",
                      fill = "firebrick4",
                      size = 3,
                      stroke = 2/3,
                      alpha = 0.5,
                      position = position_jitter()
                      ) +
                    theme_classic() +
                    # convert the genotypes of 0/1/2 on the x-axis to genotypes of ref,ref/ref,alt/alt,alt
                    # using the genotype label table define above
                    scale_x_discrete(
                      labels =  GenotypeLabelTable %>%
                                filter(
                                  .data = .,
                                  Genotype_Numeric %in% sort(x = unique(x = DataToPlot$Variant))
                                  ) %>%
                                 pull(.data = .,Genotype_alleles)
                        ) +
                    xlab(label = VariantID) +
                    ylab(label = currentMetabolite) +
                    theme(
                      axis.title = element_text(
                                                face = "bold",
                                                size = 44#,
                                                #family = "Arial"
                                                ),
                      axis.text.x = element_text(
                                                face = "bold",
                                                size = 44,
                                                #family = "Arial",
                                                angle = 45,
                                                hjust = 1
                                                ),
                      axis.text.y = element_text(
                                                face = "bold",
                                                size = 44#,
                                                #family = "Arial"
                                                )
                      ) +
                  annotate(
                    geom = "text",
                    x = 2,
                    y = 0.1,
                    label = paste0(
                                    SummaryStatisticsForLabelingPlot$Pvalue_label,
                                    "\n",
                                    SummaryStatisticsForLabelingPlot$BETA_label,
                                    "\n",
                                    SummaryStatisticsForLabelingPlot$N_label
                                   ),
                    fontface = "bold", 
                    size = 15
                  )
                
                # extract the name of the gene from the nameOfcurrentGenoMatrix
                currentGene <-
                  nameOfcurrentGenoMatrix %>%
                  str_extract(string = .,pattern = "(CYP3A4|DHCR7|CYP2R1|CYP27B1|SULT2A1|UGT1A4|GC|VDR|CASR|CUBN|CYP24A1|RXRG|RXRB|RXRA)")

                # create a directory labeled with the gene and the metabolite, if it doesn't already exist
                if(!dir.exists(paths = glue("./{currentGene}_{currentMetabolite}_boxplots/")))
                {
                  dir.create(path = glue("./{currentGene}_{currentMetabolite}_boxplots/"))
                }

                
                # save the plot to the directory labeled with the current variant ID
                png(
                  filename = glue("./{currentGene}_{currentMetabolite}_boxplots/{currentVariant}_plot.png"),
                  width = 1000,
                  height = 1000
                    )
                print(x = PlotToReturn)
                dev.off()
                
                return(PlotToReturn)

              }
              )
        },mc.cores = detectCores()
        )

  },
  GenotypeMatricesWithMetaboliteDataAdded,
  names(x = GenotypeMatricesWithMetaboliteDataAdded),
  SIMPLIFY = FALSE
)

# switch the warning messages back on now that the plots have been created
options(warn = 0)
