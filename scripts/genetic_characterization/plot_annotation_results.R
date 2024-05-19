# this is a script for plotting annotation results from ANNOVAR, pharmGKB, and VEP annotation of variant calling data

# load necessary packages
source("../scripts/load_R_packages.R")
source("../scripts/genetic_characterization/genetic_characterization_functions.R")

# load the ANNOVAR/pharmGKB/VEP annotation files for all variants and study participants, all singleton variants, and all INDELs
# define the relative paths to the files to load for the the annotated variant calling data
FileNamesToLoad <- 
  c(
    # all variant calling data
    "../annotation_output/allVariants_annovar_pharmGKB_and_VEP_data.tsv",
    # singleton variants
    "../annotation_output/singletons_annovar_pharmGKB_and_VEP_data.tsv",
    # INDEL variants
    "../annotation_output/indels_annovar_pharmGKB_and_VEP_data.tsv"
  )

AnnotationDataForPlots <-
FileNamesToLoad %>%
  lapply(
    X = .,
    FUN = function(currentFileToLoad)
    {
      LoadedFile <-
      currentFileToLoad %>%
      read.delim(file = .,header = TRUE)
      
      return(LoadedFile)
    }
      )

# name the list element that contains the loaded annotation dataset
names(x = AnnotationDataForPlots) <- c("allVariants","singletons","indels")

# create a plot of exonic variants and their exonic function class grouped by gene 
# also, create a table with the total number of exonic variants, exonic singletons, and exonic INDELs
  mapply(
    FUN = function(currentDataSet,currentNameOfDataSet)
    {
      ExonicVariants <-
      currentDataSet %>%
        # filter to exonic variants based on the Consequence_VEP column
        FilterToExonicVariantsFromVEPannotation(
          InputDataFrame = .,
          NameOfVEPconsequenceColumn = "Consequence_VEP"
          ) %>%
        # replace any occurrance of "UGT1A1" in the gene column with "UGT1A" since UGT1A1 was not truly sequenced
        # UGT1A4 was actually sequenced but the UGT1A genes (UGT1A1,3-10) have overlapping transcripts so the whole region is annotated as UGT1A1
        UpdateUGT1AfamilyGeneNames(
          InputDataFrame = .,
          VEPgeneColumnName = "Gene_VEP"
            ) %>%
        # remove underscores from the VEP consequence labels and put a space after any commas
        mutate(
          .data = .,
          Consequence_VEP = Consequence_VEP %>%
                            lapply(
                              X = .,
                              FUN = function(currentLabel)
                              {
                                LabelToReturn <-
                                currentLabel %>%
                                  gsub(pattern = "_",replacement = " ",x = .) %>%
                                  gsub(pattern = ",",replacement = ", ",x = .)
                                
                                return(LabelToReturn)
                              } 
                                ) %>%
                              unlist(x = .)
          )
      
      # print a message with number of exonic variants present in the current data set
      print(x = "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$")
      print(x = paste0("There are ",nrow(x = ExonicVariants)," exonic variants in the ",currentNameOfDataSet," dataset !"))
      print(x = "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$")
      
      DataToPlot <-
        ExonicVariants %>%
        # select relevant gene columns and the VEP consequence columns
        select(
          .data = .,
          Gene_VEP,
          Consequence_VEP
          ) %>%
        # group the results by gene, then by exonic function
        group_by(
          .data = .,
          Gene_VEP,
          Consequence_VEP
        ) %>%
        # count the number of annotations grouped by gene and exonic function
        summarise(
          .data = .,
          AnnotationCount = n()
        ) %>%
        # filter the annotation count to counts greater than 0
        filter(
          .data = .,
          AnnotationCount > 0
        )
      
      # print a message with the number of exonic variants grouped by gene and by exonic function class for the current dataset
      print(x = "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$")
      print(x = paste0("Exonic variant count summary for ",currentNameOfDataSet,"........"))
      print(x = DataToPlot)
      print(x = "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$")
      
      PlotToReturn <-
      # create a stacked bar chart of the results grouped by gene and colored by exonic function class
        # with ggplot2 from https://www.r-graph-gallery.com/48-grouped-barplot-with-ggplot2.html
      DataToPlot %>%
        ggplot(
          data = .,
          mapping =
            aes(
              # group the bars by gene
              x = Gene_VEP,
              # height of the bars is the annotation count
              y = AnnotationCount,
              # fill the bars by exonic functional consequence
              fill = Consequence_VEP
            )
        ) +
        # create a stacked bar chart with the height of the bar as the frequency (identity)
        geom_bar(
          position = "stack",
          stat = "identity",
          color = "black"
        ) +
        scale_fill_viridis(discrete = TRUE) +
        scale_y_continuous(
          limits = c(0,90),
          breaks = c(0,5,10,20,30,40,50,60,70,80,90),
          labels = c(0,5,10,20,30,40,50,60,70,80,90)
            ) +
        theme_classic() +
        xlab(label = "Gene Names") +
        ylab(label = "Number of Variants") +
        labs(fill = "Exonic Function") +
        theme(
          axis.text.x =
                        element_text(
                            family = "Arial",
                            angle = 45,
                            face = "bold.italic",
                            size =8,
                            hjust = 1
                                ),
          axis.text.y = element_text(
                                     family = "Arial",
                                     size = 11
                                     ),
          plot.title = element_text(family = "Arial",face = "bold",size = 11,hjust = 0.5),
          axis.title = element_text(family = "Arial",face = "bold",size = 11),
          legend.text = element_text(family = "Arial",face = "bold",size = 8),
          legend.title = element_text(family = "Arial",face = "bold",size = 8,hjust = 0.5),
          legend.position = c(.99, .99),
          legend.justification = c("right", "top"),
          legend.box.just = "right",
          legend.margin = margin(6, 6, 6, 6),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.ticks.x = element_line(colour = "black"),
          axis.ticks.length.x = unit(x = 3,units = "pt")
        )
      
      # create a directory for storing annotation plots if it does not already exist
      if(!dir.exists(paths = "./plots/"))
      {
          dir.create(path = "./plots/")

          # create a subdirectory for exonic annotation plots if it does not exist
          if(!dir.exists(paths = "./plots/exonic_variant_summary/"))
          {
            dir.create(path = "./plots/exonic_variant_summary/")
          }

      }
      
      # save a table with the number of exonic variants grouped by gene and by exonic function class for the current dataset
      DataToPlot %>%
        write_tsv(
          x = .,
          file = glue("./plots/exonic_variant_summary/exonic_variant_count_by_gene_and_function_{currentNameOfDataSet}.txt")
        )
      

      # save the results as a .tiff file labeled by the file name
      glue("./plots/exonic_variant_summary/{currentNameOfDataSet}_exonic_summary_by_gene.tiff") %>%
      tiff(
        filename = .,
        compression = "none",
        res = 300,
        units = "in",
        width = 7,
        height = 5
        )
      
      print(x = PlotToReturn)
      dev.off()

      print(x = "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^")
      print(x = glue("{currentNameOfDataSet} plot of exonic variant annotation summary for each gene was created"))
      
      # return a table with the number of exonic variants in the current data set
      DataToReturn <-
        tibble(
          "NumberOfExonicVariants" = nrow(x = ExonicVariants),
          "DataSet" = currentNameOfDataSet
        )
      
      return(DataToReturn)

    },
    AnnotationDataForPlots,
    names(x = AnnotationDataForPlots),
    SIMPLIFY = FALSE
      ) %>%
    do.call(
      what = "rbind",
      args = .
      ) %>%
    # save a table with the counts of the exonic variants for all variants, singletons, and indels to the file system
    write_tsv(
      x = .,
      file = "./plots/exonic_variant_summary/ExonicVariantCounts.txt"
        )

