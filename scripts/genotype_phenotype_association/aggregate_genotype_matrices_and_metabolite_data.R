# this is a script for joining genotype matrices with the metabolite data
# the resulting data from this is used for plotting the metabolite data grouped by genotype for each variant for each gene for each metabolite
# the resulting data from this is also used for performing multiple regression with the genotypes and demographic factors (Age, Gender, and BMI)

# load necessary packages
source("../scripts/load_R_packages.R")

# load functions from genetic characterization code to reduce existing VEP variant IDs to rs-numbers only
source("../scripts/genetic_characterization/genetic_characterization_functions.R")

# load other necessary functions
source("../scripts/genotype_phenotype_association/genetic_association_analysis_functions.R")

### create an rsID, snpID lookup table for variants that were significant (P <= 0.05) in linear and logistic regressions

PLINKvariantID_dbSNP_rsID_lookupTable <-
  c(
    "./linear_regression_results_file_and_annotation_common_variants.txt",
    "./logistic_regression_results_file_and_annotation_common_variants.txt"
    ) %>%
  lapply(
    X = .,
    FUN = function(currentFileName)
    {
      FileToReturn <-
        # load the metabolite regression results that are joined with the variant annotation data
        read.delim(
          file = currentFileName,
          header = TRUE,
          stringsAsFactors = FALSE
          ) %>%
          # filter to variants with a P-value <= 0.05
          filter(
            .data = .,
            Pvalue <= 0.05
            ) %>%
          # select only the variant ID, and dbSNP annotation columns
          select(
            SNP,
            avsnp138,
            avsnp142,
            snp138,
            existing_variant_VEP
            ) %>%
          # remove extra variant IDs that are not rs-numbers from the existing_variant_VEP column if an rs-number is present
          RemoveAlternateVariantIDsAndReplaceMissingRSnumbers(
            InputDataFrame = .,
            rsIDcolumnName = "existing_variant_VEP"
            ) %>%
          # for any dbsnp annotation column labeled with a ".", change it to an NA
          # for any VEP variant ID column with a "*rsNA", change it to an NA
          mutate_at(
            .tbl = .,
            .vars = c("avsnp138","avsnp142","snp138","existing_variant_VEP"),
            .funs = function(currentColumnVector)
                    {
                      VectorToReturn <-
                      currentColumnVector %>%
                      lapply(
                        X = .,
                        FUN = function(currentElement)
                        {
                          if(isTRUE(x = currentElement==as.character(x = ".")) | isTRUE(x = currentElement==as.character(x = "*rsNA")))
                          {
                            ValueToReturn <- as.character(x = NA)
                          } else {
                            ValueToReturn <- currentElement
                          }

                          return(ValueToReturn)
                        }
                        ) %>%
                      unlist(x = .)

                      return(VectorToReturn)
                    }
            )
      
      
        return(FileToReturn)
    }
      ) %>%
    do.call(
      what = "rbind",
      args = .
      ) %>%
    # obtain distinct rows based on the SNP column
    distinct(
      .data = .,
      SNP,
      .keep_all = TRUE
      ) %>%
    # create a column with the PLINK variant ID converted to the format that the genotype matrices are in
    # the PLINK variant ID column has the format chrom:pos:ref:alt
    # the genotype matrices have the format Xchrom.pos.ref.alt_alt
    ConvertPLINKvariantIDtoGenotypeMatrixFormat(DataFrameToUse = .)
  

#save the lookup table of PLINK variant IDs and dbSNP rsIDs to the file system
PLINKvariantID_dbSNP_rsID_lookupTable %>%
write_tsv(
  x = .,
  file = "PLINKvariantID_dbSNP_rsID_lookupTable_commonVariants.txt"
    )

# load the plink continuous phenotype file with metabolite data
PhenotypeData <-
  read.delim(
    file = "VitaminDplinkPhenotypeFile.txt",
    header = TRUE
  )

PhenotypeData <-
  PhenotypeData %>%
  # deselect unnecessary columns
  select(
    .data = .,
    -IID
  ) %>%
  # rename the sample ID column
  select(
    .data = .,
    "sampleID" = FID,
    everything()
  ) %>%
  # change any phenotype values that are -9 to NA for the metabolite data columns
  mutate_at(
    .tbl = .,
    # identify metabolite columns based on the ngperml or the ratio suffix
    .vars = names(x = PhenotypeData) %>%
      grepl(pattern = "(ngperml|ratio)",x = .) %>%
      which(x = .) %>%
      names(x = PhenotypeData)[.],
    # if the current column value is -9, change it to an NA
    .funs = function(currentColumn)
    {
      ArrayToReturn <-
        # loop through all the values of the current column
        # if the value is a -9, change it to an NA,
        # if the value is not a -9, leave it as is
        currentColumn %>%
        as.numeric(x = .) %>%
        lapply(
          X = .,
          FUN = function(currentValue)
          {
            if(currentValue==as.numeric(x = -9))
            {
              ValueToReturn <- as.numeric(x = NA)
            } else
            {
              ValueToReturn <- currentValue
            }

            return(ValueToReturn)
          }
        ) %>%
        unlist(x = .)

      return(ArrayToReturn)
    }
  )


FilesInDirectory <-
# list the files in the genotype matrices temporary directory
list.files(path = "./genotype_matrices_temp/")

FilesToLoad <-
  FilesInDirectory %>%
  # identify indexes of files that end with the prefix raw
  grepl(
    pattern = "raw",
    x = .
      ) %>%
  which(x = .) %>%
  # extract the files with the suffix raw
  FilesInDirectory[.]

GenotypeMatricesWithMetaboliteDataAdded <-
# load all of the genotype matrices and join each one with the metabolite data
FilesToLoad %>%
  lapply(
    X = .,
    FUN = function(currentFileToLoad)
    {
      LoadedFile <-
      # load the current file
        read.table(
          file = paste0("./genotype_matrices_temp/",currentFileToLoad),
          header = TRUE
          ) %>%
        # create a column with a gene label
        mutate(
          .data = .,
          Gene = currentFileToLoad %>%
                 str_extract(
                   string = .,
                   pattern = "(CYP3A4|DHCR7|CYP2R1|CYP27B1|SULT2A1|UGT1A4|GC|VDR|CASR|CUBN|CYP24A1|RXRG|RXRB|RXRA)"
                     )
          ) %>%
        # deselect unnecessary columns
        select(
          .data = .,
          -IID,
          -PAT,
          -MAT,
          -SEX,
          -PHENOTYPE
          ) %>%
        # rename the sample ID column
        select(
          .data = .,
          "sampleID" = FID,
          everything()
          ) %>%
        # join the genotype matrix file with the PLINK phenotype file by sample ID
        left_join(
          x = .,
          y = PhenotypeData,
          by = "sampleID"
          )

      return(LoadedFile)
    }
      )


# name the list of genotype matrices with metabolite data added according to the file name
names(x = GenotypeMatricesWithMetaboliteDataAdded) <- FilesToLoad

# make a directory for saving the genotype matrices with metabolite data added if it does not exist
if(!dir.exists(paths = "./genotype_matrices_with_metabolite_data/"))
{
  dir.create(path = "./genotype_matrices_with_metabolite_data/")
}

# save each of the genotype matrices with metabolite data added to the file system
mapply(
  FUN = function(GenotypeMatrix,NameOfGenotypeMatrix)
  {
    GenotypeMatrix %>%
    write_tsv(
      x = .,
      file = glue("./genotype_matrices_with_metabolite_data/{NameOfGenotypeMatrix}")
        )
  },
  GenotypeMatricesWithMetaboliteDataAdded,
  names(x = GenotypeMatricesWithMetaboliteDataAdded)
    )


# load the plink binary phenotype file with vitamin D status phenotypes
BinaryPhenotypeData <-
  read.delim(
    file = "VitaminDplinkPhenotypeFileBinary.txt",
    header = TRUE
  ) %>%
  # deselect unnecessary columns
  select(
    .data = .,
    -IID
  ) %>%
  # rename the sample ID column
  select(
    .data = .,
    "sampleID" = FID,
    everything()
  ) %>%
  # change any phenotype values that are -9 to NA for the binary phenotype columns
  mutate_at(
    .tbl = .,
    c("VitDsufficiency","VitDdeficiency"),
    # if the current column value is -9, change it to an NA
    .funs = function(currentColumn)
    {
      ArrayToReturn <-
        # loop through all the values of the current column
        # if the value is a -9, change it to an NA,
        # if the value is not a -9, leave it as is
        currentColumn %>%
        as.numeric(x = .) %>%
        lapply(
          X = .,
          FUN = function(currentValue)
          {
            if(currentValue==as.numeric(x = -9))
            {
              ValueToReturn <- as.numeric(x = NA)
            } else
            {
              ValueToReturn <- currentValue
            }
  
            return(ValueToReturn)
          }
        ) %>%
        unlist(x = .)
  
      return(ArrayToReturn)
    }
  )


GenotypeMatricesWithBinaryPhenotypes <-
# load all of the genotype matrices and join each one with the binary vitamin D status phenotypes
FilesToLoad %>%
  lapply(
    X = .,
    FUN = function(currentFileToLoad)
    {
      LoadedFile <-
      # load the current file
        read.table(
          file = paste0("./genotype_matrices_temp/",currentFileToLoad),
          header = TRUE
          ) %>%
        # create a column with a gene label
        mutate(
          .data = .,
          Gene = currentFileToLoad %>%
                 str_extract(
                   string = .,
                   pattern = "(CYP3A4|DHCR7|CYP2R1|CYP27B1|SULT2A1|UGT1A4|GC|VDR|CASR|CUBN|CYP24A1|RXRG|RXRB|RXRA)"
                     )
          ) %>%
        # deselect unnecessary columns
        select(
          .data = .,
          -IID,
          -PAT,
          -MAT,
          -SEX,
          -PHENOTYPE
          ) %>%
        # rename the sample ID column
        select(
          .data = .,
          "sampleID" = FID,
          everything()
          ) %>%
        # join the genotype matrix file with the PLINK phenotype file with binary phenotypes by sample ID
        left_join(
          x = .,
          y = BinaryPhenotypeData,
          by = "sampleID"
          )

      return(LoadedFile)
    }
      )

# name the list of genotype matrices with binary phenotype data added according to the file name
names(x = GenotypeMatricesWithBinaryPhenotypes) <- FilesToLoad

# make a directory for saving the genotype matrices with binary phenotype data added if it does not exist
if(!dir.exists(paths = "./genotype_matrices_with_binary_phenotypes/"))
{
  dir.create(path = "./genotype_matrices_with_binary_phenotypes/")
}

# save each of the genotype matrices with binary phenotype data added to the file system
mapply(
  FUN = function(GenotypeMatrix,NameOfGenotypeMatrix)
  {
    GenotypeMatrix %>%
    write_tsv(
      x = .,
      file = glue("./genotype_matrices_with_binary_phenotypes/{NameOfGenotypeMatrix}")
        )
  },
  GenotypeMatricesWithBinaryPhenotypes,
  names(x = GenotypeMatricesWithBinaryPhenotypes)
    )
 