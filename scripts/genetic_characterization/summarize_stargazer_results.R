# this is a script for making plots of stargazer haplotype and phenotype calls on the vitamin D variant calling data

# load necessary packages
source("../scripts/load_R_packages.R")

StarAlleleCallingData <-
# load all of the stargazer star allele calling results
c(
  "cyp2r1",
  "ugt1a4"
  ) %>%
  lapply(
    X = .,
    FUN = function(currentGene)
    {
      DataToReturn <-
      glue("{currentGene}_results.stargazer-genotype.txt") %>%
        read.delim(
          file = .,
          header = TRUE
            )
      
      return(DataToReturn)
    }
      )
# name the list elements of the star allele calling data according to the gene
names(x = StarAlleleCallingData) <- c(
                                      "cyp2r1",
                                      "ugt1a4"
                                      )

# create bar charts of haplotype calls and phenotype calls for each gene and create tables
  mapply(
    FUN = function(currentDataSet,currentGene)
    {
      ### create a plot of haplotype frequency
      DataForHaplotypePlot <-
        currentDataSet %>%
        # select the two haplotype columns
        select(
          .data = .,
          hap1_main,
          hap2_main
          ) %>%
        # put the haplotype columns into one array
        unlist(x = .) %>%
        # remove the column names
        unname(obj = .) %>%
        # make the array a single column in a tibble
        tibble("haplotypes" = .) %>%
        # group by the star alleles that are in the column
        group_by(
          .data = .,
          haplotypes
          ) %>%
        # count each of the star alleles
        summarise(
          .data = .,
          haplotypeCount = n()
          ) %>%
        # ungroup the star allele column
        ungroup(x = .) %>%
        # create a column with a label for the current gene
        mutate(
          .data = .,
          gene = rep_len(
                         x = toupper(x = currentGene),
                         length.out = nrow(x = .)
                         )
          ) %>%
        # create a column with haplotype frequencies
        mutate(
          .data = .,
          haplotypeFrequency = round(x = haplotypeCount/sum(haplotypeCount),digits = 4)
          ) %>%
        # create a column with the total haplotypes (should be 2*the number of individuals in the study)
        mutate(
          .data = .,
          TotalHaplotypes = sum(haplotypeCount)
          ) 
      
      HaplotypePlot <-
        DataForHaplotypePlot %>%
        # create a bar chart filled by frequency of the star allele
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
        ylab(label = "Allele Frequency") +
        labs(fill = "Star Allele") +
        ggtitle(
          label = paste0(
                        toupper(x = currentGene),
                        " ",
                        "Haplotypes"
                        )
          ) +
        theme(
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y = element_text(family = "Arial",face = "bold",size = 30),
          plot.title = element_text(family = "Arial",face = "bold",size = 40,hjust = 0.5),
          axis.title.y = element_text(family = "Arial",face = "bold",size = 40,margin = margin(r = 15)),
          legend.text = element_text(family = "Arial",face = "bold",size = 40),
          legend.title = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.ticks.length.x = unit(x = 3,units = "pt"),
          plot.margin=unit(c(1,1,1,1),"cm")
            ) +
        scale_y_continuous(expand = c(0,0)) +
        scale_fill_manual(values = c("#D55E00", "#56B4E9", "#009E73","#F0E442"))
      
      
      ### create a plot of phenotype frequency based on the star allele diplotype
      DataForPhenotypePlot <-
        currentDataSet %>%
        # select the phenotype column
        select(
          .data = .,
          phenotype
          ) %>%
        # remove the _metabolizer suffix from the phenotype string
        mutate(
          .data = .,
          phenotype = phenotype %>%
                      gsub(
                        pattern = "_metabolizer",
                        replacement = "",
                        x = .
                          )
          ) %>%
        # group by the phenotype
        group_by(
          .data = .,
          phenotype
          ) %>%
        # count each phenotype
        summarise(
          .data = .,
          "count" = n()
          ) %>%
        # ungroup by phenotype
        ungroup(x = .) %>%
        # create a column with a label for the current gene
        mutate(
          .data = .,
          gene = rep_len(
                         x = toupper(x = currentGene),
                         length.out = nrow(x = .)
                         )
          ) %>%
        # compute the frequency of each phenotype and the total sample size
        mutate(
          .data = .,
          Frequency = round(x = count/sum(count),digits = 4),
          N = sum(count)
          )
      
      PhenotypePlot <-
        DataForPhenotypePlot %>%
      # create a bar chart filled by frequency of the phenotype call
      ggplot(
        data = .,
        aes(
          x = gene,
          y = count,
          fill = phenotype,
          )
        ) +
      geom_col(position = "fill") +
      theme_classic() +
      ylab(label = "Phenotype Frequency") +
      labs(fill = "Metabolizer") +
      ggtitle(
        label = paste0(
                      toupper(x = currentGene),
                      " ",
                      "Phenotypes"
                      )
        ) +
      theme(
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(family = "Arial",face = "bold",size = 30),
        plot.title = element_text(family = "Arial",face = "bold",size = 40,hjust = 0.5),
        axis.title.y = element_text(family = "Arial",face = "bold",size = 40,margin = margin(r = 15)),
        legend.text = element_text(family = "Arial",face = "bold",size = 40),
        legend.title = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.ticks.length.x = unit(x = 3,units = "pt"),
        plot.margin=unit(c(1,1,1,1),"cm")
          ) +
      scale_y_continuous(expand = c(0,0)) +
      scale_fill_manual(
        breaks = c("rapid","normal","intermediate","unknown"),
        values = c("#E69F00", "#56B4E9", "#009E73","#F0E442")
        )
      
      # save the haplotype frequency and phenotype frequency plots to the file system
      glue("./{currentGene}_haplotypes.png") %>%
        png(
          filename = .,
          width = 1000,
          height = 1000
          )
      print(x = HaplotypePlot)
      dev.off()
      
      glue("./{currentGene}_phenotypes.png") %>%
        png(
          filename = .,
          width = 1000,
          height = 1000
          )
      print(x = PhenotypePlot)
      dev.off()
      
      # create a haplotype table to save to the file system
      HaplotypeTableToSave <-
        DataForHaplotypePlot %>%
        select(
          .data = .,
          "Gene" = gene,
          "Haplotype" = haplotypes,
          "Haplotype Frequency (N = 938)" = haplotypeFrequency
        )
      
      # save the haplotype frequency table to the file system
      HaplotypeTableToSave %>%
        write_tsv(
          x = .,
          file = glue("./{currentGene}_HaplotypeFrequencyTable.tsv")
            )
      
      # create a diplotype/phenotype table to save to the file system
      DiplotypeTableToSave <-
        currentDataSet %>%
        # create a diplotype column by pasting together the main haplotypes
        mutate(
          .data = .,
          Diplotype = mapply(
                            FUN = function(hap1,hap2)
                            {
                              DiplotypeToReturn <-
                              c(
                                hap1,
                                hap2
                               ) %>%
                              # sort before pasting so equivalent diplotype permutations
                              # are combined into one
                              sort(x = .) %>%
                              paste0(.,collapse = "|")
                              
                              return(DiplotypeToReturn)
                            },
                            hap1_main,
                            hap2_main,
                            SIMPLIFY = FALSE
                            ) %>%
                            unlist(x = .)
        ) %>%
        # select the diplotype and phenotype columns
        select(
          .data = .,
          Diplotype,
          phenotype
          ) %>%
        # group by diplotype and phenotype
        group_by(
          .data = .,
          Diplotype,
          phenotype
          ) %>%
        # tally the diplotype and phenotype combination
        summarise(
          .data = .,
          "Count" = n()
          ) %>%
        # ungroup the table
        ungroup(x = .) %>%
        # create a column with a label for the current gene
        mutate(
          .data = .,
          Gene = rep_len(
                        x = toupper(x = currentGene),
                        length.out = nrow(x = .)
                        )
        ) %>%
        # compute the diplotype/phenotype combination frequency
        mutate(
          .data = .,
          Frequency = round(x = Count/sum(Count),digits = 4)
          ) %>%
        # rearrange and relabel columns
        select(
          .data = .,
          Gene,
          Diplotype,
          "Phenotype" = phenotype,
          "Frequency (N = 469)" = Frequency
        ) %>%
        # remove underscores in the phenotype column and replace with spaces
        mutate(
          .data = .,
          Phenotype = Phenotype %>%
                      gsub(
                        pattern = "_",
                        replacement = " ",
                        x = .
                        )
          )
      
      # save the diplotype frequency table to the file system
      DiplotypeTableToSave %>%
        write_tsv(
          x = .,
          file = glue("./{currentGene}_DiplotypeFrequencyTable.tsv")
          )
      
      print(x = "~~~~~~~~~~~~~ StarGazer summary ~~~~~~~~~~~~~~~~~")
   
      print(x = HaplotypeTableToSave)
      print(x = DiplotypeTableToSave)
      
      print(x = "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
     
    },
    StarAlleleCallingData,
    names(x = StarAlleleCallingData)
      )
  
  print(x = "###################################################")
  print(x = "##        stargazer plots were created !         ##")
  print(x = "###################################################")