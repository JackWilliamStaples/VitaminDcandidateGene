# this is a script for creating a plink range file for extracting nucleotide positions 
# in a plink file for a specified chromosome and range of nucleotide positions

# load necessary packages
source("../scripts/load_R_packages.R")

GeneListAndPositions <-
# load the gene list bed file that was used for the sequencing
read.delim(
  file = "../invoice_gene_list_BED_and_coverage_summaries_batch_1/gene_list.bed.txt",
  stringsAsFactors = FALSE,
  header = FALSE
  )

# change the names of the columns of the file
names(x = GeneListAndPositions) <- c("chrom","start","end","gene")

# split the GeneListAndPositions dataframe into a list of dataframes for each individual gene
GeneListAndPositions %>%
  group_by(
    .data = .,
    gene
    ) %>%
  group_split(
    .tbl = .
  ) %>%
  # loop through each individual gene-specific data frame,
  # select the relevant columns to create a PLINK range file,
  # and save the table to the file system as a text file named according to the gene
  lapply(
    X = .,
    FUN = function(currentGenePositionRange)
    {
      currentGenePositionRange %>%
        select(
          .data = .,
          "CHR" = chrom,
          "BP1" = start,
          "BP2" = end,
          "LABEL" = gene
          ) %>%
        write_tsv(
          x = .,
          file = paste0(
                    unique(x = currentGenePositionRange$gene),
                    "_",
                    "regionFile.txt"
                    ),
          col_names = FALSE
          )
    }
    )

print(x = "$$$$$$$$$$$$$$$$$$$$$ gene region files created for filtering VCF data $$$$$$$$$$$$$$$$$$")