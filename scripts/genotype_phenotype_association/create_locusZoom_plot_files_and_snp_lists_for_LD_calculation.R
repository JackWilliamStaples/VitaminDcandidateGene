# load necessary packages
source("../scripts/load_R_packages.R")
# load necessary functions
source("../scripts/genotype_phenotype_association/genetic_association_analysis_functions.R")

# read in the linear additive allele association results for all variants
PLINKalleleDoseRegressionResults <-
  read.delim(
    file = "./linear_regression_results_file_and_annotation_common_variants.txt",
    header = TRUE
  ) %>%
  # ensure that the base pair column is a numeric
  mutate(
    .data = .,
    BP = BP %>% as.numeric(x = .)
  )

# read in the logistic regression association results for all variants
PLINKlogisticRegressionResults <-
  read.delim(
    file = "./logistic_regression_results_file_and_annotation_common_variants.txt",
    header = TRUE
    ) %>%
  # ensure that the base pair column is a numeric
  mutate(
    .data = .,
    BP = BP %>% as.numeric(x = .)
    )

# save a text file with the variant name (MarkerName) and p values (P-value) for LocusZoom
# for each metabolite and metabolite ratio and for each gene
PLINKalleleDoseRegressionResults %>%
  CreateMarkerNameAndPvalueFilesForLocusZoom(
    DataFrameToUse = .,
    PhenotypeColumnName = "Metabolite",
    GeneColumnName = "Gene"
    )

# save a text file with the variant name (MarkerName) and p values (P-value) for LocusZoom
# for each binary phenotype and for each gene
PLINKlogisticRegressionResults %>%
  CreateMarkerNameAndPvalueFilesForLocusZoom(
    DataFrameToUse = .,
    PhenotypeColumnName = "Phenotype",
    GeneColumnName = "Gene"
      )

# save a list of variant IDs to filter the VCF file used for calculating CSKT linkage disequilibrium
# that are the same exact variants that are plottted for the LocusZoom association results
# for multi-allelic sites, the variant is selected based on the lowest P-value for the multi-allelic site
# do this for each metabolite and metabolite ratio, for each gene
PLINKalleleDoseRegressionResults %>%
  CreateLocusZoomSNPlists(
    DataFrameToUse = .,
    PhenotypeColumnName = "Metabolite",
    GeneColumnName = "Gene"
    )

# save a list of variant IDs to filter the VCF file used for calculating CSKT linkage disequilibrium
# that are the same exact variants that are plottted for the LocusZoom association results
# for multi-allelic sites, the variant is selected based on the lowest P-value for the multi-allelic site
# do this for each binary phenotype and for each gene
PLINKlogisticRegressionResults %>%
  CreateLocusZoomSNPlists(
    DataFrameToUse = .,
    PhenotypeColumnName = "Phenotype",
    GeneColumnName = "Gene"
    )
