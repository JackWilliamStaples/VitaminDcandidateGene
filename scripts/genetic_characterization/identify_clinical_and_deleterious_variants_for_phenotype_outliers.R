# this is a script for filtering annotated variant calling data of vitamin D phenotype outliers
# to variants that have clinical significance or a deleterious CADD score of at least 20

# load necessary packages
source("../scripts/load_R_packages.R")

# list the files in the ./phenotype_outliers/ directory
FilesToLoad <-
  list.files(path = "./phenotype_outliers/")

# combine the annotated variants for each individual into one dataframe
AnnotatedVariants <-
  # load each of the files into memory
  FilesToLoad %>%
  lapply(
    X = .,
    FUN = function(currentFile)
    {
      
      FileToReturn <-
        currentFile %>%
          paste0("./phenotype_outliers/",.) %>%
          # load the current file
          read.delim(
            file = .,
            header = TRUE
            ) %>%
          # rename the allele frequency column to make it more informative
          dplyr::rename(
            .data = .,
            "allele_frequency" = Otherinfo1
            ) %>%
          # filter to variants with an allele frequency above 0 for the current individual
          filter(
            .data = .,
            allele_frequency > 0
            ) %>%
          # create a column with a label for the current individual sequencing code
          mutate(
            .data = .,
            Sample_code = currentFile %>%
                          gsub(
                            pattern = "annotated_variants_",
                            replacement = "",
                            x = .
                            ) %>%
                          gsub(
                            pattern = ".hg19_multianno.txt",
                            replacement = "",
                            x = .
                            ) %>%
                          str_split(
                            string = .,
                            pattern = "_"
                            ) %>%
                          unlist(x = .) %>%
                          unique(x = .)
            ) %>%
        # convert missing phred scaled CADD values to NA
        # then convert the phred scaled CADD score column to a numeric
        mutate(
          .data = .,
          CADD_phred = if_else(
                               condition = CADD_phred == ".",
                               true = as.character(x = NA),
                               false = CADD_phred
                               )
          ) %>%
        mutate(
          .data = .,
          CADD_phred = CADD_phred %>% as.numeric(x = .)
          ) %>%
        # filter to variants with a pathogenic clinVar phenotype 
        # or a Phred scaled CADD score of at least 20 with no clinVar phenotype
        filter(
          .data = .,
          (
          (CLNSIG == "Pathogenic") |
          (CLNSIG == "Likely_pathogenic") |
          (CLNSIG == "Pathogenic/Likely_pathogenic") |
          (CLNSIG == "Pathogenic/Likely_pathogenic,_other")
          ) |
          (
          (CADD_phred >= 20) & (CLNSIG == ".")
          )
          ) %>%
        # reorder and rename columns
        select(
          .data = .,
          Sample_code,
          CADD_phred,
          allele_frequency,
          everything()
          ) %>%
        # group by the sample code and
        group_by(
          .data = .,
          Sample_code
          ) %>%
        # arrange rows based on CADD score in descending order
        arrange(
          .data = .,
          CADD_phred
          )
      
      return(FileToReturn)
      
    }
      ) %>%
  do.call(
    what = "rbind",
    args = .
    ) %>%
  # save the results to the file system
  write_tsv(
    x = .,
    file = "../annotation_output/variant_annotation_for_phenotype_outliers.txt"
      )
