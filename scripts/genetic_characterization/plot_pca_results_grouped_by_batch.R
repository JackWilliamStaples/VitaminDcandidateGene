# this is a script for plotting the first two principal components of principal component analysis performed on the vitamin D genetic data
# that is output from the PLINK quality control, colored based on whether the sample was sequenced in round 1 or round 2 of sequencing
# if there are two distinct clusters that are formed in this plot, then there MUST be a batch effect that is present and 
# if there ARE NOT distinct clusters (i.e., the data is randomly scattered) then there is likely not a batch effect

# load necessary packages
source("../scripts/load_R_packages.R")

# read in the principal components analysis that is output from plink based on a linkage disequilibrium pruned set of SNPs
pca_table <- 
  read.table(
    file = "pca_results.eigenvec",
    header = FALSE, 
    comment.char = ""
    )

pca_df <- 
  pca_table %>% 
  as.data.frame(x = .) 

# read in the control_case_batch1_batch2.txt file that was created earlier in the vitamin D quality control pipeline 
# to use for labeling samples based on sequencing round in the PCA plot below
samples_with_batch_numbers <-
  # read in the table with sample IDs and batch numbers
  read.delim(
    file = "control_case_batch1_batch2.txt",
    header = FALSE,
    stringsAsFactors = FALSE
    ) %>%
   # convert to a dataframe
   data.frame() %>%
   # rename columns
   select(
     .data = .,
     "sampleID" = V1,
     "batch" = V2
     ) %>%
   # remove the duplicate sample ID
   mutate(
     .data = .,
     sampleID = sampleID %>%
                lapply(
                  X = .,
                  FUN = function(currentSampleID)
                  {
                    ValueToReturn <-
                      currentSampleID %>%
                      str_split(string = .,pattern = " ") %>%
                      unlist(x = .) %>%
                      unique(x = .)
                    
                    return(ValueToReturn)
                  }
                    )
     ) %>%
  # make sure the sample ID is a character vector
  mutate(
    .data = .,
    sampleID = sampleID %>% as.character(x = .),
    batch = batch %>% as.character(x = .)
    )

# create a plot of the first two principal components with points colored by batch number to see if a batch sequencing effect exists
pca_plot <-
  # change the name of the V1 column in the pca results to sampleID
  pca_df %>%
    dplyr::rename(
      .data = .,
      "sampleID" = V1
      ) %>%
    # deselect the duplicate sample ID column
    select(
      .data = .,
      -V2
      ) %>%
    # change the sample ID column from a numeric to a character vector
    mutate(
      .data = .,
      sampleID = sampleID %>%
                 as.character(x = .)
      ) %>% 
    # change the unnamed columns of the first and second principal component to "PC1" and "PC2"
    select(
      .data = .,
      sampleID,
      "PC1" = V3,
      "PC2" = V4
      ) %>%
    left_join(
      x = .,
      y = samples_with_batch_numbers,
      by = "sampleID"
      ) %>%
  # plot the results as a scatter plot of the second principal component versus the first, with points colored by batch number
   ggplot(
     data = .,
     mapping = 
       aes(
         x = PC1,
         y = PC2, 
         color = batch
         )
   ) +
   geom_point(
     size = 0.8, 
     alpha = 1,
     colour = "black",
     pch = 21,
     aes(fill = batch)
     ) +
    # use the classic theme for the plot
   theme_classic() +
   scale_fill_discrete(name = "Sequencing round") +
   theme(
     legend.title = element_text(family = "Arial",hjust = 0.5,face = "bold",size = 11),
     axis.title = element_text(family = "Arial",face = "bold",size = 11),
     axis.text = element_text(family = "Arial",size = 11,face = "bold"),
     legend.text = element_text(family = "Arial",face = "bold",size = 11),
     legend.position = c(0.8,0.6),
     legend.box.background = element_rect(colour = "black",size = 1.5)
     )
  


# save the result as a tiff file
tiff(
  filename = "./pca_plot.tiff",
  compression = "none",
  res = 300,
  units = "in",
  width = 7,
  height = 5
  )

pca_plot
dev.off()



