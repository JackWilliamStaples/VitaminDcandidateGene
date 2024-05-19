# this is a script for assigning CYP3A4 star allele haplotypes using:
  # (1) CYP3A4 star allele data curated by pharmVar for GRCh37 and 
  # (2) phased sample x haplotype matrices of CYP3A4 variant calling data created with BCFtools

# load necessary packages
source("../scripts/load_R_packages.R")
# load genetic characterization functions
source("../scripts/genetic_characterization/genetic_characterization_functions.R")

RawCYP3A4DataFromPharmVar <-
# load the cyp3A4 star allele definition table from PharmVar
read.delim(
  file = "../cyp3a4_pharmVar_starAlleleData/CYP3A4-4.3/GRCh37/CYP3A4.NC_000007.13.haplotypes.tsv",
  header = TRUE,
  skip = 3
    )

# add column names to the PharmVar data
names(x = RawCYP3A4DataFromPharmVar) <- c("Haplotype.Name","Gene","rsID","sequence","start","end","ref","alt","type")

# convert the data from pharmVAR to a one-hot encoded star allele X snpID format for assigning haplotypes
StarAlleleSNPmatrix <-
  RawCYP3A4DataFromPharmVar %>%
  select(
    .data = .,
    Haplotype.Name,
    rsID
    ) %>%
  mutate(
    .data = .,
    Haplotype.Name = Haplotype.Name %>% as.character(x = .),
    rsID = rsID %>% as.factor(x = .)
    ) %>%
  as.data.table(x = .) %>%
  one_hot(dt = .) %>%
  as_tibble(x = .) %>%
  # add two rows to the wildtype CYP3A4 star alleles, denoted as REFERENCE in the pharmVAR file with a 0 for every rsID
  rbind(
    c("CYP3A4*1",rep_len(x = 0,length.out = nrow(x = .)-1)),
    c("CYP3A4*1.001",rep_len(x = 0,length.out = nrow(x = .)-1)),
    .
  )

# removed the "rsID_" that is appended to every rs number
names(x = StarAlleleSNPmatrix) <-
  StarAlleleSNPmatrix %>%
  names(x = .) %>%
  lapply(
    X = .,
    FUN = function(currentString)
    {
      StringToReturn <-
       currentString %>%
        str_replace(
          string = .,
          pattern = "rsID_",
          replacement = ""
          )
      
      return(StringToReturn)
    }
    ) %>%
  unlist(x = .)

CSKTvariantsWithStarAlleles <-
# load in a table of variant IDs and rs numbers from the CSKT data
read.delim(
  file = "../annotation_output/allVariants_annovar_pharmGKB_and_VEP_data.tsv",
  header = TRUE
  ) %>%
  select(
    .data = .,
    snp138,
    avsnp138,
    avsnp142,
    existing_variant_VEP,
    "variantID" = VariantID_fromInputFile
    ) %>%
  # remove variant IDs that are not rs-numbers from the existing_variant_VEP column
  RemoveAlternateVariantIDsAndReplaceMissingRSnumbers(
    InputDataFrame = .,
    rsIDcolumnName = "existing_variant_VEP"
    ) %>%
  # remove variants without rs numbers
  filter(
    .data = .,
    snp138 != "." & avsnp138 != "." & avsnp142 != "." & existing_variant_VEP != "*rsNA"
    ) %>%
  # filter to variants that are part of existing star alleles
  filter(
    .data = .,
    (snp138 %in% (StarAlleleSNPmatrix %>% names(x = .) %>% .[2:length(x = .)])) |
    (avsnp138 %in% (StarAlleleSNPmatrix %>% names(x = .) %>% .[2:length(x = .)])) |
    (avsnp142 %in% (StarAlleleSNPmatrix %>% names(x = .) %>% .[2:length(x = .)])) |
    (existing_variant_VEP %in% (StarAlleleSNPmatrix %>% names(x = .) %>% .[2:length(x = .)]))
    ) %>%
  # convert the variant ID to the same format as in the CYP3A4hapSampleMatrix (chr:position_ref_alt)
  mutate(
    .data = .,
    variantID = variantID %>%
                lapply(
                  X = .,
                  FUN = function(currentVariantID)
                  {
                    SplitVariantID <-
                    currentVariantID %>%
                      str_split(string = .,pattern = ":") %>%
                      unlist(x = .)

                    names(x = SplitVariantID) <-
                      c("chr","position","ref","alt")

                   VariantToReturn <-
                     paste0(
                       SplitVariantID["chr"],
                       ":",
                       SplitVariantID["position"],
                       "_",
                       SplitVariantID["ref"],
                       "_",
                       SplitVariantID["alt"]
                       )

                   return(VariantToReturn)

                  }
                    ) %>%
                  unlist(x = .)
    )

# load in the sample file that corresponds to the CYP3A4 haplotype matrix
# and create an array of sample haplotype IDs
SampleIDsForHaplotypeMatrix <-
read.table(
  file = "./CYP3A4_haplotypes.sample",
  header = TRUE
    ) %>%
  filter(
    .data = .,
    ID_1 != "0"
    ) %>%
  mutate(
    .data = .,
    hap1_label = ID_1 %>% paste0(.,"_hap1"),
    hap2_label = ID_2 %>% paste0(.,"_hap2")
    ) %>%
  select(
    .data = .,
    -ID_1,
    -ID_2,
    -missing
  ) %>%
  unlist(x = .) %>%
  unname(obj = .) %>%
  sort(x = .)

# load in the phased CYP3A4 sample x haplotype matrix for the vitamin D study participants
CYP3A4hapSampleMatrix <-
read.table(
  file = "./CYP3A4_haplotypes.hap",
  header = FALSE
    ) %>%
  # filter to the variants that actually have star alleles
  filter(
    .data = .,
    V2 %in% CSKTvariantsWithStarAlleles$variantID
    )

# name the columns of the CYP3A4hapSampleMatrix with the variant information and sample IDs
names(x = CYP3A4hapSampleMatrix) <- c(
                                      c("chr","snpID","position","ref","alt"),
                                      SampleIDsForHaplotypeMatrix
                                      )

# add rs numbers based on the variant positions and alleles by joining the CSKTvariantsWithStarAlleles with the CYP3A4hapSampleMatrix
CYP3A4hapSampleMatrix <-
  CYP3A4hapSampleMatrix %>%
  mutate(
    .data = .,
    rsID = CSKTvariantsWithStarAlleles %>%
            select(
              .data = .,
              variantID,
              avsnp142
              ) %>%
            left_join(
              x = .,
              y = CYP3A4hapSampleMatrix,
              by = c("variantID" = "snpID")
                ) %>%
             select(
               .data = .,
               avsnp142
               ) %>%
             unlist(x = .) %>%
             unname(obj = .)
    )

# create a function for filling in sample haplotypes with star alleles from a star allele column
# based on whether there is a 1 or a 0 for a star allele haplotype
AssignStarAllele <-
  function(
    possibleStarAlleles = c("CYP3A4*3","CYP3A4*1G","CYP3A4*22","CYP3A4*10"),
    currentSampleHaplotypeVector = NULL
    )
  {
    StarAlleleVectorToReturn <-
    data.frame(
      "StarAllele" = possibleStarAlleles,
      "SNPs" = currentSampleHaplotypeVector
      ) %>%
      mutate(
        .data = .,
        SNPs = SNPs %>% as.numeric(x = .)
        ) %>%
      mutate(
        .data = .,
        SNPs = if_else(
                        condition = SNPs == as.numeric(x = 1),
                        true = StarAllele,
                        false = "0"
                          )
        ) %>%
      pull(
        .data = .,
        SNPs
        )
    
    return(StarAlleleVectorToReturn)
  
  }

# add star alleles by joining the RawCYP3A4DataFromPharmVar with the CYP3A4hapSampleMatrix by rsID
CYP3A4hapSampleMatrix <-
  CYP3A4hapSampleMatrix %>%
  mutate(
    .data = .,
    starAllele = RawCYP3A4DataFromPharmVar %>%
                  select(
                    .data = .,
                    Haplotype.Name,
                    rsID
                    ) %>%
                  left_join(
                    x = CYP3A4hapSampleMatrix,
                    y = .,
                    by = "rsID"
                    ) %>%
                  select(
                    .data = .,
                    Haplotype.Name
                    ) %>%
                  filter(
                    .data = .,
                    !grepl(
                      x = Haplotype.Name,
                      pattern = ".00"
                    )
                    ) %>%
                  unlist(x = .) %>%
                  unname(obj = .)
    ) %>%
  # now that the star alleles have been added to the CYP3A4hapSampleMatrix, 
  # select only the star allele column and all sample haplotyple columns by deselecting the other columns
  select(
    .data = .,
    -chr,
    -snpID,
    -rsID,
    -position,
    -ref,
    -alt
    ) %>%
  select(
    .data = .,
    starAllele,
    everything()
    ) %>%
  mutate_at(
    .tbl = .,
    .vars = vars(matches(match = "(hap1|hap2)")),
    .funs = ~ AssignStarAllele(currentSampleHaplotypeVector = .x)
      ) %>%
# then transpose the result matrix to get a haplotype sample X star allele dataframe
  t(x = .) %>%
  as.data.frame(x = .)

ColumnNames <-
  CYP3A4hapSampleMatrix[1,] %>%
  unlist(x = .) %>%
  unname(obj = .)

CYP3A4hapSampleMatrix <-
# make the rownames into an actual row in the dataframe
CYP3A4hapSampleMatrix %>%
  mutate(
    .data = .,
    SampleID = rownames(x = .)
    ) %>%
  # remove the first row of the dataframe now that it is no longer needed
  filter(
    .data = .,
    SampleID != "starAllele"
  ) %>%
  tibble()

# replace the column names with the contents of the first row of the dataframe
names(x = CYP3A4hapSampleMatrix) <- c(ColumnNames,"SampleID")
  
### assign the star allele haplotype for each chromosome for each participant
CYP3A4hapSampleMatrix <-
  CYP3A4hapSampleMatrix %>%
  # assign the final star allele haplotype for each participant
  mutate(
    .data = .,
    Haplotype = mapply(
                       FUN = function(
                                      star1Gvector,
                                      star3vector,
                                      star10vector,
                                      star22vector
                                      )
                       {
                         # combine each star allele column into a single array
                         currentHaplotype <-
                         c(
                           star1Gvector,
                           star3vector,
                           star10vector,
                           star22vector
                          )
                        
                        # if *1G, *3, *10, and *22 are not present (i.e., all elements in currentHaplotype are "0") 
                        # in the current haplotype, set as the wildtype *1 haplotype 
                        if(all(currentHaplotype=="0")==TRUE)
                        {
                          haplotypeToReturn <- "*1"
                        } else {
                         
                          # identify the non-zero star alleles in the current haplotype,
                          # remove the "CYP3A4" string from the star allele and return the final result
                          haplotypeToReturn <-
                            (currentHaplotype == "0") %>%
                            # identify elements that are NOT zero
                            which(x = !.) %>%
                            unlist(x = .) %>%
                            # subset the currentHaplotype with the elements that are NOT zero
                            currentHaplotype[.] %>%
                            # remove the "CYP3A4" string in the haplotype
                            lapply(
                              X = .,
                              FUN = function(currentStarAllele)
                              {
                                DataToReturn <-
                                  currentStarAllele %>%
                                  gsub(
                                    pattern = "CYP3A4",
                                    replacement = "",
                                    x = .
                                    )
                              }
                                ) %>%
                            unlist(x = .)
                        
                          # if the haplotypeToReturn array is longer than length 1,
                          # meaning there is a novel combination of star alleles present, combine the star alleles
                          # into a single string and separate them with a plus sign ("+")
                          if(length(x = haplotypeToReturn)>1)
                          {
                            haplotypeToReturn <-
                            haplotypeToReturn %>%
                            paste0(.,collapse = "+")
                            
                          }
                          
                        }
                         
                         return(haplotypeToReturn)
                      
                       },
                       `CYP3A4*1G`,
                       `CYP3A4*3`,
                       `CYP3A4*10`,
                       `CYP3A4*22`,
                       SIMPLIFY = FALSE
                         ) %>%
                      unlist(x = .)
    )

HaplotypeFrequencyPlot <-
#### create a frequency plot from 0 to 1 for the haplotypes
CYP3A4hapSampleMatrix %>%
  # select the haplotype column
  select(
    .data = .,
    Haplotype
    ) %>%
  # group by the haplotype
  group_by(
    .data = .,
    Haplotype
    ) %>%
  # count each haplotype
  summarise(
    .data = .,
    haplotypeCount = n()
    ) %>%
  # create a columnn with a label for the gene
  mutate(
    .data = .,
    Gene = "CYP3A4"
    ) %>%
  ggplot(
    data = .,
    mapping = aes(
                  x = Gene,
                  y = haplotypeCount,
                  fill = Haplotype
                  )
      ) +
    geom_col(position = "fill") +
    theme_classic() +
    ylab(label = "Allele Frequency") +
    labs(fill = "Star Allele") +
    ggtitle(
      label = paste0(
                    "CYP3A4",
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
    scale_y_continuous(expand = c(0,0))
 
# create a table with haplotypes and haplotype counts 
HaplotypeFrequencyTable <-
  CYP3A4hapSampleMatrix %>%
  # select the haplotype column
  select(
    .data = .,
    Haplotype
  ) %>%
  # group by the haplotype
  group_by(
    .data = .,
    Haplotype
  ) %>%
  # count each haplotype
  summarise(
    .data = .,
    "count" = n()
  ) %>%
  # ungroup the table
  ungroup(x = .) %>%
  # compute the frequency of each haplotype and
  # compute the total sample size
  mutate(
    .data = .,
    HaplotypeFrequency = round(x = count/sum(count),digits = 4),
    SampleSize = sum(count)
    )
  

# create a table of sample IDs and diplotypes
DiplotypeData <-
CYP3A4hapSampleMatrix %>%
  # select the sampleIDs and haplotypes
  select(
    .data = .,
    SampleID,
    Haplotype
    ) %>%
  # remove "_hap1" and "_hap2" from the sample ID column
  mutate(
    .data = .,
    SampleID = SampleID %>%
               str_replace(
                 string = .,
                 pattern = "(_hap1|_hap2)",
                 replacement = ""
                   )
    ) %>%
  # group by the sampleID
  group_by(
    .data = .,
    SampleID
    ) %>%
  # split into a list of dataframes based on the sample ID
  group_split(.tbl = .) %>%
  lapply(
    X = .,
    FUN = function(currentDataFrame)
    {
      # create a string of each haplotype split with a "|" to create a diplotype
      Diplotype <-
        c(
          currentDataFrame$Haplotype[1],
          currentDataFrame$Haplotype[2]
        ) %>%
        # sort before pasting together so the diplotype with the lowest number is always in front
        sort(x = .) %>%
        paste0(
          .,
          collapse = "|"
          )
      # return a dataframe of the diplotype and the sample ID
      DataToReturn <-
        data.frame(
          "Diplotype" = Diplotype,
          "SampleID" = currentDataFrame$SampleID %>% unique(x = .)
        )
      
      return(DataToReturn)
    }
      ) %>%
  # bind each dataframe by row
  do.call(
    what = "rbind",
    args = .
    ) 

##### create a dipltoype frequency plot
DiplotypeFrequencyPlot <-
DiplotypeData %>%
  # create a gene column
  mutate(
    .data = .,
    Gene = "CYP3A4"
    ) %>%
  group_by(
    .data = .,
    Diplotype
    ) %>%
  # count each haplotype
  summarise(
    .data = .,
    DiplotypeCount = n()
    ) %>%
  # create a columnn with a label for the gene
  mutate(
    .data = .,
    Gene = "CYP3A4"
    ) %>%
  ggplot(
    data = .,
    mapping = aes(
                  x = Gene,
                  y = DiplotypeCount,
                  fill = Diplotype
                  )
      ) +
    geom_col(position = "fill") +
    theme_classic() +
    ylab(label = "Allele Frequency") +
    labs(fill = "Star Allele") +
    ggtitle(
      label = paste0(
                    "CYP3A4",
                    " ",
                    "Diplotypes"
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
    scale_y_continuous(expand = c(0,0))
  
### create a diplotype count table
DiplotypeFrequencyTable <-
DiplotypeData %>%
  group_by(
    .data = .,
    Diplotype
    ) %>%
  # count each diplotype
  summarise(
    .data = .,
    "count" = n()
    ) %>%
  # ungroup the table
  ungroup(x = .) %>%
  # compute the frequency of each diplotype
  # and the total sample size
  mutate(
    .data = .,
    DiplotypeFrequency = round(x = count/sum(count),digits = 4),
    SampleSize = sum(count)
    )

# create a directory for saving cyp3a4 star allele calling data
if(!dir.exists(paths = "./custom_cyp3A4starAllele_calls/"))
{
  dir.create(path = "./custom_cyp3A4starAllele_calls/")
}

# save the haplotype and diplotype frequency plots
PlotsToSave <-
list(
  "HaplotypeFrequencyPlot" = HaplotypeFrequencyPlot,
  "DiplotypeFrequencyPlot" = DiplotypeFrequencyPlot
)

mapply(
  FUN = function(currentPlot,currentPlotName)
  {
      png(
        filename = paste0("./custom_cyp3A4starAllele_calls/",currentPlotName,".png"),
        width = 1000,
        height = 1000
        )

    print(x = currentPlot)

    dev.off()
  },
  PlotsToSave,
  names(x = PlotsToSave)
    )

HaplotypeFrequencyTable <-
# reformat the haplotype frequency table
HaplotypeFrequencyTable %>%
  # add a gene column
  mutate(
    .data = .,
    Gene = "CYP3A4"
    ) %>%
  # reorder and rename columns
  select(
    .data = .,
    Gene,
    Haplotype,
    "Haplotype Frequency (N = 938)" = HaplotypeFrequency
    ) %>%
  # arrange the rows by the haplotype
  arrange(
    .data = .,
    Haplotype
    )

# reformat the diplotype frequency table
DiplotypeFrequencyTable <-
  DiplotypeFrequencyTable %>%
  # add a gene column
  mutate(
    .data = .,
    Gene = "CYP3A4"
    ) %>%
  # create a phenotype column
  mutate(
    .data = .,
    Phenotype = "not assigned"
    ) %>%
  # reorder and rename columns
  select(
    .data = .,
    Gene,
    Diplotype,
    Phenotype,
    "Frequency (N = 469)" = DiplotypeFrequency
    )

# save the hapltoype and diplotype frequency tables
TablesToSave <-
  list(
    "HaplotypeFrequencyTable" = HaplotypeFrequencyTable,
    "DiplotypeFrequencyTable" = DiplotypeFrequencyTable
  )

mapply(
  FUN = function(currentTable,currentTableName)
  {
    currentTable %>%
      write_tsv(
        x = .,
        file = paste0("./custom_cyp3A4starAllele_calls/",currentTableName,".tsv")
          )
      
  },
  TablesToSave,
  names(x = TablesToSave)
    )

print(x = "~~~~~~~~~~~~~ CYP3A4 star alleles ~~~~~~~~~~~~~~~~~")
print(x = HaplotypeFrequencyTable)
print(x = DiplotypeFrequencyTable)
print(x = "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")

