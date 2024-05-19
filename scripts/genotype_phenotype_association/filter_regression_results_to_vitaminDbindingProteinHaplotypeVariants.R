# load necessary packages
source("../scripts/load_R_packages.R")
# load functions
source(file = "../scripts/genotype_phenotype_association/genetic_association_analysis_functions.R")

RegressionResultsAndAnnotations <-
# load the regression results with variant annotation data added
read.delim(
  file = "./linear_regression_results_file_and_annotation_common_variants.txt",
  header = TRUE
)

#### filter the RegressionResultsAndAnnotations to relevant metabolites and common vitamin D binding protein polymorphism variants
RegressionResultsAndAnnotations %>%
  # filter to the primary 25(OH)D3 metabolite
  filter(
    .data = .,
    Metabolite == "primary_25OHD3_ngperml"
    ) %>%
  # filter to the common vitamin D binding protein polymorphisms (rs7041 and rs4588)
  filter(
    .data = .,
    (avsnp142 == "rs7041") | 
    (avsnp142 == "rs4588")
    ) %>%
  # add astericks to p-values with levels of significance
  mutate(
    .data = .,
    Pvalue = Pvalue %>%
             AddAstericksToPvalues(columnVector = .)
    ) %>%
  # select relevant columns
  select(
    .data = .,
    avsnp142,
    Gene,
    SNP,
    BETA,
    SE,
    "Pvalue_unadjust_OLS" = Pvalue,
    AAChange.refGeneWithVer,
    CSKT_AF,
    X1000g2015aug_all,
    X1000g2015aug_afr,
    X1000g2015aug_amr,
    X1000g2015aug_eas,
    X1000g2015aug_eur,
    X1000g2015aug_sas
    ) %>%
  # save to the file system
  write_tsv(
    x = .,
    file = "./VDBP_variants_for_presentation.txt"
      ) 


