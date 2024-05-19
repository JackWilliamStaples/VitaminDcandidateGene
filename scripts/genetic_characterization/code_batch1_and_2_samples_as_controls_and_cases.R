# load necessary packages
source("../scripts/load_R_packages.R")

# read in the list of samples from batch 1 
sample_IDs_batch1 <- read.delim(file = "sample_IDs_batch1.txt",header = FALSE)
# add a column of 1's to the sample IDs from batch 1, which specifies samples as "controls"
sample_IDs_batch1 <-
sample_IDs_batch1 %>% 
  mutate(
    .data = .,
    V2 = rep_len(
      x = 1,
      length.out = nrow(x = sample_IDs_batch1)
      )
    )
# read in the list of samples from batch 2
sample_IDs_batch2 <- read.delim(file = "sample_IDs_batch2.txt",header = FALSE)
# add a column of 2's to the sample IDs from batch 2, which specifies samples as "cases"
sample_IDs_batch2 <- 
  sample_IDs_batch2 %>%
  mutate(
    .data = .,
    V2 = rep_len(
      x = 2,
      length.out = nrow(x = sample_IDs_batch2)
      )
    )

# create a data frame of all samples from batches 1 and 2 with control = 1 and case = 2 codes
FileToReturn <-
rbind(
  sample_IDs_batch1,
  sample_IDs_batch2
  )

# save to the file system as a .txt file without column names
write_tsv(x = FileToReturn,file = "control_case_batch1_batch2.txt",col_names = FALSE)
