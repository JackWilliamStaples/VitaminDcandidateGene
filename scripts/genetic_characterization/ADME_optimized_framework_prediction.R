# load necessary packages
source("../scripts/load_R_packages.R")

# This is a script for performing the ADME optimized framework predictions for deleterious ADME variants 
# using an annotation file that is output from ANNOVAR annotation of a VCF file.

# The ADME optimized prediction framework is described in Zhou, et al., 2018, The Pharmacogenomics Journal, 19:115-126.

# define the relative paths to the files to load for the the annotated variant calling data
FileNamesToLoad <- 
  c(
    # all variant calling data
    "../annotation_output/merge_annotation.hg19_multianno.txt",
    # singleton variants
    "../singleton_annotation_output/singleton_annotation.hg19_multianno.txt",
    # INDEL variants
    "../indel_annotation_output/merge_indel_annotation.hg19_multianno.txt"
  )

# load each of the files by looping through the filenames to load and store as a list
AnnotationFilesWithADMEoptimizedPrediction <- 
  FileNamesToLoad %>%
    lapply(
      X = .,
      FUN = function(fileName)
      {
        
        file <- 
          # read in the file
        read.delim(
          file = fileName,
          header = FALSE
            )
        
        # select the columns of the file by parsing the first row of the table
        file <- 
          # select the first row
          file[1,] %>%
          # transpose to create a matrix with a single column
          t(x = .) %>%
          # convert the matrix to a dataframe of a single column
          data.frame(stringsAsFactors = FALSE) %>%
          # pull out the single column as an array
          pull(.data = .) %>%
          # remove the names of the array elements
          unname(obj = .) %>%
          # identify the array elements that are not NA or empty strings
          lapply(
            X = .,
            FUN = function(currentElement)
            {
              !(is.na(x = currentElement) == TRUE | currentElement == "")
              # the "!" is intentional here
            }
          ) %>%
          # unlist the elements
          unlist(x = .) %>%
          # identify the elements that are true
          which(x = .) %>%
          # use the true elements to select the desired columns from the file
          file[,.]
        
          # name the columns of the dataframe with the first row of the dataframe
          names(x = file) <-
            # take the first row
            file[1,] %>%
            # transpose to a single column matrix
            t(x = .) %>%
            # convert the single column matrix to a dataframe
            data.frame(stringsAsFactors = FALSE) %>%
            # pull out the first column into an array
            pull(.data = .) %>%
            # remove the names of the array
            unname(obj = .)
          
          # now that the columns have been named with the first row of the table, remove the first row from the table
          # also rename some of the "Otherinfo" columns
          file <-
            file[-1,] %>%
          dplyr::rename(
            .data = .,
            CSKT_AF = Otherinfo1,
            CHROM_fromInputFile = Otherinfo4,
            POS_fromInputFile = Otherinfo5,
            VariantID_fromInputFile = Otherinfo6,
            variant_summary_stats = Otherinfo11
            )
        
          # return a table of the annovar annotations with the ADME optimized framework prediction added
          TableToReturn <-
            file %>%
            # arrange the table by chromosome number then by starting nucleotide position
            arrange(
              .data = .,
              Chr,
              Start
              ) %>%
            # convert missing (indicated with a ".") LRT, mutation assessor, provean, vest, and phred-scaled CADD scores to NA
            mutate(
              .data = .,
              across(
                .cols = c(
                          LRT_score,
                          MutationAssessor_score,
                          PROVEAN_score,
                          VEST3_score,
                          CADD_phred
                          ),
                .fns = ~ case_when(
                                 .x == "." ~ as.character(x = NA),
                                 TRUE ~ as.character(x = .x)
                                )
                  )
              ) %>%
            # now convert LRT, mutation assessor, provean, vest and CADD scores to numeric column vectors
            mutate(
              .data = .,
              across(
                .cols = c(LRT_score,MutationAssessor_score,PROVEAN_score,VEST3_score,CADD_phred),
                .fns = ~ as.numeric(x = .x)
                )
              ) %>%
            # create columns for performing the ADME optimized deletrious predicition according to the methods of Zhou, et al., 2018, The Pharmacogenomics Journal, 19:115-126
            # The ADME optimized criteria are:
                # (a) LRT_score < 0.0025
                # (b) MutationAssessor_score > 2.0566
                # (c) PROVEAN_score < -3.286
                # (d) VEST3_score > 0.4534
                # (e) CADD_phred > 19.19
            # if the value for the score listed above is available, compute if the value satisfies the ADME-optimized framework condition,
              # otherwise return an NA
            mutate(
              .data = .,
              LRT_ADME =
                if_else(
                condition = is.na(x = LRT_score),
                true = as.numeric(x = NA),
                false = as.numeric(x = as.numeric(x = LRT_score) < 0.0025)
                ),
              MutationAssessor_ADME =
                if_else(
                condition = is.na(x = MutationAssessor_score),
                true = as.numeric(x = NA),
                false = as.numeric(x = as.numeric(x = MutationAssessor_score) > 2.0566)
                ),
              PROVEAN_ADME =
                if_else(
                  condition = is.na(x = PROVEAN_score),
                  true = as.numeric(x = NA),
                  false = as.numeric(x = as.numeric(x = PROVEAN_score) < -3.286)
                  ),
              VEST3_ADME =
                if_else(
                  condition = is.na(x = VEST3_score),
                  true = as.numeric(x = NA),
                  false = as.numeric(x = as.numeric(x = VEST3_score) > 0.4534)
                  ),
              CADD_ADME =
                if_else(
                condition = is.na(x = CADD_phred),
                true = as.numeric(x = NA),
                false = as.numeric(x = as.numeric(x = CADD_phred) > 19.19)
                )
              ) %>%
              # compute the sum of the binary ADME optimized prediction scores if they are not all NA, otherwise leave as NA
              mutate(
                .data = .,
                ADME_optimized_prediction = if_else(
                                                    condition = is.na(x = LRT_ADME) & is.na(x = MutationAssessor_ADME) & is.na(x = PROVEAN_ADME) & is.na(x = VEST3_ADME) & is.na(x = CADD_ADME),
                                                    true = as.numeric(x = NA),
                                                    false = rowSums(x = across(.cols = c(LRT_ADME,MutationAssessor_ADME,PROVEAN_ADME,VEST3_ADME,CADD_ADME)),na.rm = TRUE)/5
                                                      )
                                            
                )
          
          return(TableToReturn)
          
      }
      )

# save the annotated file of all variants as a tsv file
AnnotationFilesWithADMEoptimizedPrediction[[1]] %>%
write_tsv(
  x = .,
  file = #"../annotation_output/merge_annovar_and_pharmGKB.tsv"
    "../annotation_output/merge_annovar_and_ADMEoptimizedPredictions.tsv"
  )


# save the annotated singleton data as a tsv file
AnnotationFilesWithADMEoptimizedPrediction[[2]] %>%
write_tsv(
  x = .,
  file = #"../singleton_annotation_output/merge_singleton_annovar_and_pharmGKB.tsv"
    "../singleton_annotation_output/merge_singleton_annovar_and_ADMEoptimizedPredictions.tsv"
    )

# save the annotated INDEL data as a tsv file
AnnotationFilesWithADMEoptimizedPrediction[[3]] %>%
  write_tsv(
    x = .,
    file = #"../indel_annotation_output/merge_indel_annovar_and_pharmGKB.tsv"
      "../indel_annotation_output/merge_indel_annovar_and_ADMEoptimizedPredictions.tsv"
      )

print(x = "***************************** ADME optimized framework predictions are complete :) ***************************************")


