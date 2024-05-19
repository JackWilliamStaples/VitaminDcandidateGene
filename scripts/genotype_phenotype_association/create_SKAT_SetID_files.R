# load necessary packages
source("./scripts/load_R_packages.R")


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

# load all of the PLINK rare variant .bim files
RareVariantBIMfilesToLoad <-
list.files(path = "./rareVariantAnalysisIntermediates/") %>%
  grepl(pattern = ".bim",x = .) %>%
  which(x = .) %>%
  list.files(path = "./rareVariantAnalysisIntermediates/")[.] %>%
  lapply(
    X = .,
    FUN = function(currentFileToLoad)
    {
          # identify the batch based on the file name
          currentBatch <-
            currentFileToLoad %>%
            str_extract(
              string = .,
              pattern = "(batch1|batch2|merge)"
                )
          
          # identify the metabolite based on the file name
          currentMetabolite <-
            currentFileToLoad %>%
            str_extract(
              string = .,
              pattern = "(active_1alpha25OH2D2_ngperml|active_1alpha25OH2D3_ngperml|beta_4beta25OH2D3_ngperml|Overall_primary_25OHD_ngperml|primary_25OHD2_active_1alpha25OH2D2_ratio|primary_25OHD2_ngperml|primary_25OHD3_active_1alpha25OH2D3_ratio|beta_4beta25OH2D3_primary_25OHD3_ratio|primary_25OHD3_ngperml|primary_25OHD3_R_24R25OH2D3_ratio,R_24R25OH2D3_ngperml|VitaminD2_ngperml|VitaminD2_primary_25OHD2_ratio|VitaminD3_ngperml|VitaminD3_primary_25OHD3_ratio)"
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
            file = glue("./rareVariantAnalysisIntermediates/rare_variants_for_SKAT_analysis_{currentBatch}_{currentMetabolite}.SetID"),
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE
          )
          
          
          
    }
    )


print(x = "$$$$$$$$$$$$$ SetID files for O-SKAT analysis were saved $$$$$$$$$$$")



# create SKAT snp data set files (SSD) for each batch/metabolite
c("batch1","batch2","merge") %>%
  lapply(
    X = .,
    FUN = function(currentBatch)
    {
      c(
        "active_1alpha25OH2D2_ngperml",
        "active_1alpha25OH2D3_ngperml",
        "beta_4beta25OH2D3_ngperml",
        "Overall_primary_25OHD_ngperml",
        "primary_25OHD2_active_1alpha25OH2D2_ratio",
        "primary_25OHD2_ngperml",
        "primary_25OHD3_active_1alpha25OH2D3_ratio",
        "beta_4beta25OH2D3_primary_25OHD3_ratio",
        "primary_25OHD3_ngperml",
        "primary_25OHD3_R_24R25OH2D3_ratio",
        "R_24R25OH2D3_ngperml",
        "VitaminD2_ngperml",
        "VitaminD2_primary_25OHD2_ratio",
        "VitaminD3_ngperml",
        "VitaminD3_primary_25OHD3_ratio"
        ) %>%
        lapply(
          X = .,
          FUN = function(currentMetabolite)
          {
            # define the PLINK file names that will be loaded into SKAT and the names of the output SNP set data files (SSD) files
            Bed_file_to_load <- glue("./rareVariantAnalysisIntermediates/rare_variants_for_SKAT_analysis_{currentBatch}_{currentMetabolite}.bed")
            Bim_file_to_load <- glue("./rareVariantAnalysisIntermediates/rare_variants_for_SKAT_analysis_{currentBatch}_{currentMetabolite}.bim")
            Fam_file_to_load <- glue("./rareVariantAnalysisIntermediates/rare_variants_for_SKAT_analysis_{currentBatch}_{currentMetabolite}.fam")
            SetID_file_to_load <- glue("./rareVariantAnalysisIntermediates/rare_variants_for_SKAT_analysis_{currentBatch}_{currentMetabolite}.SetID")
            SSD_file_to_make <- glue("./rareVariantAnalysisIntermediates/rare_variants_for_SKAT_analysis_{currentBatch}_{currentMetabolite}.SSD")
            SSDinfo_file_to_make <- glue("./rareVariantAnalysisIntermediates/rare_variants_for_SKAT_analysis_{currentBatch}_{currentMetabolite}.SSD.info")
            
            # create SKAT snp data set files (SSD) for each batch/metabolite
            Generate_SSD_SetID(
              File.Bed = Bed_file_to_load,
              File.Bim = Bim_file_to_load,
              File.Fam = Fam_file_to_load,
              File.SetID = SetID_file_to_load,
              File.SSD = SSD_file_to_make,
              File.Info = SSDinfo_file_to_make
                )
            
            # read in the PLINK .fam file with continuous phenotypes (Is.binary = FALSE)
            Rare_VariantsFAMfile <- 
              Read_Plink_FAM(
                Filename = Fam_file_to_load,
                Is.binary = FALSE
              )
            
            # create an array of the phenotype values from the loaded FAM file
            PhenotypeValues <-
              Rare_VariantsFAMfile %>%
              select(
                .data = .,
                Phenotype
              ) %>%
              pull(.data = .)
            
            # open the SSD file
            Rare_VariantsSSD.INFO <-
              Open_SSD(
                File.SSD = SSD_file_to_make,
                File.Info = SSDinfo_file_to_make
              )
            
            
            # define the NULL SKAT model on the dataset with no covariates included
            # the outcome is continuous in this case
            NullSKATmodelObject <-
              SKAT_Null_Model(
                formula = PhenotypeValues ~ 1,
                out_type = "C"
              )
            
            # run the SKAT association on the NULL model defined above
            # the SKAT model will run for each individual gene defined in the the SETid file
            SKAT_association_results <-
              SKAT.SSD.All(
                SSD.INFO = Rare_VariantsSSD.INFO,
                obj = NullSKATmodelObject 
              )
            
            # save the SKAT association results table to the file system after adding a column of bonferroni adjusted P-values
            SKAT_association_results$results <-
              SKAT_association_results$results %>% 
              # create a column with Bonferroni-adjusted P-values for the number of SNP sets (i.e., number of genes) that were tested
              mutate(
                .data = .,
                P.value.Bonferroni = nrow(x = SKAT_association_results$results)*P.value
              ) %>%
              # round the P.value columns and the P.value.Bonferroni columns to 5 places after the decimal
              mutate(
                .data = .,
                P.value = P.value %>% round(x = .,digits = 5),
                P.value.Bonferroni = P.value %>% round(x = .,digits = 5)
              ) %>%
              # create a column with an arbitrary index for coloring points by gene
              mutate(
                .data = .,
                ColorIndex = seq_len(length.out = nrow(x = .)) %>% as.factor(x = .)
              )
            
            # create a directory for saving the SKAT results
            if(!dir.exists(paths = "./OSKAT_results/"))
            {
              dir.create(path = "./OSKAT_results/")
            }
            
            # save the SKAT results to the SKAT results directory
            SKAT_association_results$results %>%
              write.table(
                x = .,
                file = glue("./OSKAT_results/{currentBatch}_{currentMetabolite}.txt"),
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
                y = "P.value.Bonferroni",
                color = "ColorIndex",                                # Color by groups
                #palette = "viridis", # Custom color palette
                sorting = "descending",                       # Sort value in descending order
                rotate = TRUE,                                # Rotate vertically
                dot.size = 4,                                 # Large dot size
                #y.text.col = TRUE,                            # Color y text by groups
                ggtheme = theme_pubr()                        # ggplot2 theme
              )+
              geom_hline(yintercept = 0.05,linetype="solid",color="red",size=1.5) +
              theme_cleveland() + # Add dashed grids
              scale_fill_viridis(discrete = TRUE) +
              rremove(object = "legend") +
              theme(
                plot.title = element_text(hjust = 0.5),
                axis.title.x = element_text(
                  #family = "sans", size = 15, 
                  margin=margin(5,0,0,0)
                ), 
                axis.title.y = element_text(
                  #family = "sans", size = 15, 
                  margin=margin(0,5,0,0)
                )
              )
            # add additional plot formatting and change names of titles
            SKAT_association_results_plot <-
              ggpar(
                p = SKAT_association_results_plot,
                font.main = "bold",
                main = glue("{currentBatch} {currentMetabolite}"),
                xlab = "Gene",
                font.x = "bold",
                ylab = "P-value (Bonf. adjust for number of genes)",
                font.y = "bold",
                font.ytickslab = "italic",
                yticks.by = 0.05
              )
          
            # save the final plot to the file system
            png(
              filename = glue("./OSKAT_results/{currentBatch}_{currentMetabolite}.png"),
              width = 10,
              height = 10,
              units = 'in',
              res = 300
              )
            print(x = SKAT_association_results_plot)
            dev.off()
            
          }
            )
        
    }
      )

