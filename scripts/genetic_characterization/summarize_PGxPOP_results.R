# this is a script for parsing the output of a star allele haplotype, diplotype, and phenotype assignment file generated with PGx-POP software

# load necesssary packages
source("../scripts/load_R_packages.R")

# load the UGT1A1 star allele assignment results without the header
StarAlleleData <-
  read.delim(
    file = "./PGx_POP_results_allVariants_UGT1A1.txt",
    header = FALSE,
    sep = ",",
    skip = 1
    ) %>%
  # deselect the empty column that is created
  select(
    .data = .,
    -`V16`
    )

# load the header for the StarAlleleData
StarAlleleDataHeader <-
  # read in the PGxPOP results file
  read.delim(
    file = "./PGx_POP_results_allVariants_UGT1A1.txt",
    header = TRUE
    ) %>%
  # extract the header
  names(x = .) %>%
  # convert the header to a string
  toString(x = .) %>%
  # split the string on dots (.)
  str_split(
    string = .,
    pattern = "\\."
    ) %>%
  # unlist the result
  unlist(x = .)

# add the header to the StarAlleleData
names(x = StarAlleleData) <- StarAlleleDataHeader

# change the phenotype column to all lower case letters
StarAlleleData <-
  StarAlleleData %>%
  mutate(
    .data = .,
    phenotype = phenotype %>% tolower(x = .)
  )

#### compute haplotype counts and haplotype frequencies
HaplotypeSummary <-
  # put both haplotype columns into an array
  c(
    StarAlleleData$hap_1,
    StarAlleleData$hap_2
    ) %>%
  # put the haplotypes into a single column of a tibble
  tibble(
    "Haplotype" = .
  ) %>%
  # group by the haplotype
  group_by(
    .data = .,
    Haplotype
    ) %>%
  # count the haplotypes
  summarise(
    .data = .,
    "HaplotypeCount" = n()
    ) %>%
  # ungroup the results
  ungroup(x = .) %>%
  # calculate haplotype frequencies and total sample size
  mutate(
    .data = .,
    HaplotypeFrequency = round(x = HaplotypeCount/sum(HaplotypeCount),digits = 4),
    SampleSize = sum(HaplotypeCount)
    ) %>%
  # create a gene column
  mutate(
    .data = .,
    Gene = "UGT1A1"
    ) %>%
  # arrange and relabel columns
  select(
    .data = .,
    Gene,
    Haplotype,
    "Haplotype Frequency (N = 938)" = HaplotypeFrequency
    )

### compute diplotype counts, diplotype frequencies, and phenotype count/frequency
DiplotypeSummary <-
StarAlleleData %>%
  select(
    .data = .,
    diplotype,
    phenotype
    ) %>%
  # combine equivalent diplotypes by sorting the diplotype permutations
  mutate(
    .data = .,
    diplotype = diplotype %>%
                lapply(
                  X = .,
                  FUN = function(currentDiplotype)
                  {
                    
                    DiplotypeToReturn <-
                    currentDiplotype %>%
                      # split the current diplotype based on the verticle pipe ("|")
                      str_split(
                        string = .,
                        pattern = "\\|"
                        ) %>%
                      unlist(x = .) %>%
                      # sort the haplotypes
                      sort(x = .) %>%
                      # recombine back into a diplotype now that it is sorted
                      paste0(.,collapse = "|")
                      
                      return(DiplotypeToReturn)
                  }
                    ) %>%
                  unlist(x = .)
    ) %>%
  # group by diplotype and phenotype
  group_by(
    .data = .,
    diplotype,
    phenotype
    ) %>%
  # count the diplotypes and phenotypes
  summarise(
    .data = .,
    "Count" = n()
    ) %>%
  # ungroup the result
  ungroup(x = .) %>%
  # compute the frequencies of the diplotypes and phenotypes
  mutate(
    .data = .,
    Frequency = round(x = Count/sum(Count),digits = 4)
    ) %>%
  # compute the sample size
  mutate(
    .data = .,
    SampleSize = sum(Count)
    ) %>%
  # make a gene column
  mutate(
    .data = .,
    Gene = "UGT1A1"
  ) %>%
  # rearrange and relabel columns
  select(
    .data = .,
    Gene,
    "Diplotype" = diplotype,
    "Phenotype" = phenotype,
    "Frequency (N = 469)" = Frequency
    )

# save the haplotype and diplotype tables to the file system
HaplotypeSummary %>%
  write_tsv(
    x = .,
    file = "./ugt1a1_HaplotypeFrequencyTable.tsv"
    )

DiplotypeSummary %>%
  write_tsv(
    x = .,
    file = "./ugt1a1_DiplotypeFrequencyTable.tsv"
    )

print(x = "~~~~~~~~~~~~~ PGx-POP summary ~~~~~~~~~~~~~~~~~")
print(x = HaplotypeSummary)
print(x = DiplotypeSummary)


print(x = "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
