# this is a script for parsing the results of calling star alleles with pypgx software

# load necesssary packages
source("../scripts/load_R_packages.R")

  # list the files that contain pypgx results
  # by searching for the string "pipeline"
  list.files(path = "./") %>%
    grepl(
      pattern = "pipeline",
      x = .
      ) %>%
  list.files(path = "./")[.] %>%
  # load each of the pypgx results files and summarize the results
  lapply(
    X = .,
    FUN = function(currentFile)
    {
      fileToReturn <-
        # load the current file
        read.delim(
          file = glue("{currentFile}/data.tsv"),
          header = TRUE
            ) %>%
        # sort the genotype column so all permutations are sorted the same way
        mutate(
          .data = .,
          Genotype = Genotype %>%
                     lapply(
                       X = .,
                       FUN = function(currentGenotype)
                       {
                         ResultToReturn <-
                           currentGenotype %>%
                           str_split(
                             string = .,
                             pattern = "/"
                           ) %>%
                             unlist(x = .) %>%
                             sort(x = .) %>%
                             paste(.,collapse = "|")
                         
                         return(ResultToReturn)
                       }
                         ) %>%
                      unlist(x = .)
                     
          ) %>%
        # select the genotype column
        select(
          .data = .,
          Genotype,
          Phenotype
          ) %>%
        # group by the genotype and phenotype
        group_by(
          .data = .,
          Genotype,
          Phenotype
          ) %>%
        # count the number of occurances of the genotype/phenotype pair
        summarise(
          .data = .,
          Count = n()
          ) %>%
        # ungroup the result
        ungroup(x = .) %>%
        # calculate the frequency of the genotype/phenotype pair by dividing the count by the total count
        mutate(
          .data = .,
          Frequency = Count/sum(Count)
          ) %>%
        # add a column for the current gene based on the current file name
        mutate(
          .data = .,
          Gene = gsub(
                      pattern = "-pipeline",
                      replacement = "",
                      x = currentFile
                      ) %>%
                  gsub(
                    pattern = "grch37-",
                    replacement = "",
                    x = .
                  ) 
          ) %>%
        # move the gene column to the front of the table
        select(
          .data = .,
          Gene,
          everything()
          )
      
      return(fileToReturn)
    }
      ) %>%
    do.call(
      what = "rbind",
      args = .
        ) %>%
    # save the results to the file system
    write_delim(
      x = .,
      file = "./pypgx_results_allGenes.txt"
        )
  
  
  