# load necessary packages
source("./scripts/load_R_packages.R")
source("./scripts/genetic_characterization/genetic_characterization_functions.R")

# load the ANNOVAR, pharmGKB, and VEP annotation data
FileToLoad <- "./annotation_output/allVariants_annovar_pharmGKB_and_VEP_data.tsv"

ANNOVARpharmGKBandVEPdata <-
  FileToLoad %>%
  # load the annotation file
  read.delim(
    file = .,
    header = TRUE
      ) %>%
  # select relevant columns
  select(
    .data = .,
    Chr:variant_summary_stats,
    LRT_ADME:starAllele
    ) %>%
  # Replace any occurrence of "UGT1A1" in the gene column with "UGT1A" since UGT1A1 was not truly sequenced.
  # UGT1A4 was actually sequenced but the UGT1A genes (UGT1A1,3-10) have overlapping transcripts so the whole region is annotated as UGT1A1
  UpdateUGT1AfamilyGeneNames(
    InputDataFrame = .,
    VEPgeneColumnName = "Gene_VEP"
  ) %>%
  # create an HGVSg column for the genomic location and ref/alt alleles 
  mutate(
    .data = .,
    NucleotideID = paste0(
                          "GRCh37.p13:",
                          Location_VEP,
                          USED_REF,
                          ">",
                          Allele_VEP
                          )
    ) %>%
  # replace any amino acids described as "%3D" with an equal sign for synonymous variants
  mutate(
    .data = .,
    HGVSp = if_else(
                    condition = grepl(x = HGVSp,pattern = "%3D")==TRUE,
                    true = sub(pattern = "%3D",replacement = "=",x = HGVSp),
                    false = HGVSp
                      )
    ) %>%
  # add a space before the exon and intron number so the number does not
  # convert to a date when the table is opened with Microsoft Excel
  mutate(
    .data = .,
    exon_number_VEP = exon_number_VEP %>%
                      paste0(" ",.),
    intron_number_VEP = intron_number_VEP %>%
                        paste0(" ",.)
    ) %>%
  # sort by chromosome then by nucleotide position
  arrange(
    .data = .,
    Chr,
    Start
  ) 


##### identify how many variants are known and novel based on the absence of an rs-identification number or existing variant ID
# identify number of novel variants based on not having a dbSNP annotation
NumberOfNovelVariants <-
ANNOVARpharmGKBandVEPdata %>%
  select(
    .data = .,
    avsnp142,
    avsnp138,
    snp138,
    existing_variant_VEP
    ) %>%
  filter(
    .data = .,
      (avsnp142 == ".") & 
      (avsnp138 == ".") & 
      (snp138 == ".") & 
      (existing_variant_VEP == "-")
    ) %>%
  nrow(x = .)

  NumberOfNovelVariants %>%
  paste0(
    "########################### There are ",
    .,
    " novel variants out of ",
    nrow(x = ANNOVARpharmGKBandVEPdata),
    " total variants ",
    "in the vitamin D variant call data ##############################"
    ) %>%
  print(x = .)
      
# identify how many variants are known to dbSNP based on the presence of an rs identification number or existing variant ID
NumberOfKnownVariants <-
ANNOVARpharmGKBandVEPdata %>%
  select(
    .data = .,
    avsnp142,
    avsnp138,
    snp138,
    existing_variant_VEP
  ) %>%
  filter(
    .data = .,
    (avsnp142 != ".") | 
    (avsnp138 != ".") | 
    (snp138 != ".") | 
    (existing_variant_VEP != "-")
  ) %>%
  nrow(x = .)

  NumberOfKnownVariants %>%
  paste0(
    "########################### There are ",
    .,
    " known variants out of ",
    nrow(x = ANNOVARpharmGKBandVEPdata),
    " total variants",
    " ##############################"
    ) %>%
  print(x = .)

###### identify how many variants have minor allele frequency (MAF) < 0.05 and how many variants are 
     # both novel (not in dbSNP database or has existing variant ID) AND have MAF < 0.05
     # I am classifying a novel variant as a variant with no dbSNP rs-identification number or existing variant ID

VariantsWithMAFlessThan0.05 <-
  ANNOVARpharmGKBandVEPdata %>%
  # Filter the annotation data to variants with an alternate allele frequency < 0.05 OR a (1 - alternate allele frequency) < 0.05,
    # This will filter to variants with a MINOR allele frequency (MAF < 0.05)
  filter(
    .data = .,
      (CSKT_AF < 0.05) | 
      ((1 - CSKT_AF) < 0.05)
    ) %>%
  nrow(x = .)

print(
  x = paste0(
             "####################### There are ",
             VariantsWithMAFlessThan0.05,
             " variants with minor allele frequency (MAF) < 0.05"
             )
  )

print(
  x = paste0(
             "####################### There are ",
             nrow(x = ANNOVARpharmGKBandVEPdata)-VariantsWithMAFlessThan0.05,
             " variants with minor allele frequency (MAF) ≥ 0.05"
             )
  )

VariantsThatAreNovelAndHaveMAFlessThan0.05 <-
  ANNOVARpharmGKBandVEPdata %>%
  # filter the annotation data to variants with a minor allele frequency (MAF) of less than 0.05
  # and no rs-identification number from any of the dbSNP columns or an existing variant ID
  filter(
    .data = .,
    (
      (CSKT_AF < 0.05) | 
      ((1 - CSKT_AF) < 0.05)
    ) & 
    (avsnp142 == ".") & 
    (avsnp138 == ".") & 
    (snp138 == ".") & 
    (existing_variant_VEP == "-")
  ) %>%
  nrow(x = .)

print(
  x = paste0(
             "####################### There are ",
             VariantsThatAreNovelAndHaveMAFlessThan0.05,
             " novel (no dbsnp rsID or existing variant ID) variants with minor allele frequency (MAF) < 0.05"
             )
  )


###### identify how many variants have minor allele count (MAC) < 5 and a minor allele count (MAC) >= 1
###### as well as variants with minor allele count of at least 5
###### and additionally, variants with a minor allele count of zero

# load the variants with a MAC < 5 and a MAC >= 1
variantsWithMACofAtMost4 <-
  read.table(
    file = "./annovar_annotation/variantsWithMACofAtMost4.bim",
    header = FALSE
    ) %>%
  nrow(x = .)

# load the variants with a MAC >= 5
variantsWithMACofAtLeast5 <-
  read.table(
    file = "./annovar_annotation/variantsWithMACofAtLeast5.bim",
    header = FALSE
    ) %>%
  nrow(x = .)

# load the variants with a MAC = 0
variantsWithMACofZero <-
  read.table(
    file = "./annovar_annotation/variantsWithMACofZero.bim",
    header = FALSE
    ) %>%
  nrow(x = .)

print(
  x = paste0(
            "############################ There are ",
            variantsWithMACofAtLeast5,
            " variants with minor allele count (MAC) ≥ 5 out of ",
            variantsWithMACofAtLeast5 + variantsWithMACofAtMost4,
            "total variants with a minor allele count of at least 1"
            )
  )

print(
  x = paste0(
            "############################ There are ",
            variantsWithMACofAtMost4,
            " variants with minor allele count (MAC) < 5 and MAC >= 1 out of ",
            variantsWithMACofAtLeast5 + variantsWithMACofAtMost4,
            " total variants with a minor allele count (MAC) of at least 1"
            )
  )

print(
  x = paste0(
            "########################## There are ",
            variantsWithMACofZero,
            "variants with minor allele count (MAC) = 0"
            )
    )

##################################################################################

##### identify the number of EXONIC variants that are novel
NumberOfNovelExonicVariants <-
    ANNOVARpharmGKBandVEPdata %>%
      # filter to variants that are both novel (no rs-identification number or existing variant ID) AND exonic 
      # exonic variants include: 
        # (1) frameshift_variant, 
        # (2) [frameshift_variant,splice_region_variant], 
        # (3) inframe_deletion, 
        # (4) missense_variant, 
        # (5) [missense_variant,splice_region_variant], 
        # (6) [splice_region_variant,synonymous_variant], 
        # (7) start_lost, 
        # (8) stop_gained, 
        # (9) synonymous_variant
      filter(
        .data = .,
        ((avsnp142 == ".") & (avsnp138 == ".") & (snp138 == ".") & (existing_variant_VEP == "-")) &
          ( 
            (Consequence_VEP == "frameshift_variant") |
            (Consequence_VEP == "frameshift_variant,splice_region_variant") |
            (Consequence_VEP == "inframe_deletion") |
            (Consequence_VEP == "missense_variant") |
            (Consequence_VEP == "missense_variant,splice_region_variant") |
            (Consequence_VEP == "splice_region_variant,synonymous_variant") |
            (Consequence_VEP == "start_lost") |
            (Consequence_VEP == "stop_gained") |
            (Consequence_VEP == "synonymous_variant")
          )
      ) %>%
    nrow(x = .)

print(
  x = paste0(
             "####################### There are ",
             NumberOfNovelExonicVariants,
             " novel exonic variants"
             )
  )

# identify variants in the merged dataset that have  associated, pathogenic, likely pathogenic, or drug response phenotypes in clinVar
# or a pharmGKB clinical annotation with a non-missing evidence level
NoteworthyClinVarAndPharmGKBannnotations <-
ANNOVARpharmGKBandVEPdata %>%
  filter(
    .data = .,
          (CLNSIG == "Pathogenic") | 
          (CLNSIG == "Likely_pathogenic") | 
          (CLNSIG == "Pathogenic/Likely_pathogenic") |
          (CLNSIG == "Pathogenic/Likely_pathogenic,_other") |
          (CLNSIG == "drug_response") |
          (CLNSIG == "association") |
          !is.na(x = pharmGKB_EvidenceLevel)
  ) 

# print the number of variants with ClinVar associated, pathogenic, likely pathogenic, 
# or drug response phenotypes or pharmGKB clinical annotations with a non-missing level of evidence
print(
  x = paste0(
            "################# There are ",
            nrow(x = NoteworthyClinVarAndPharmGKBannnotations),
            " variants with ClinVar (associated, pathogenic, likely pathogenic, or drug response phenotypes) or pharmGKB clinical annotations with a non-missing level of evidence ##########################"
            )
  )

# create a directory for saving the data if it doesn't exist
if(!dir.exists(paths = "./ClinicalVariantsForPresentations/"))
{
  dir.create(path = "./ClinicalVariantsForPresentations/")
}

CADDphredGreaterThan20Table <-
# identify variants with a Phred-scaled CADD score of at least 20 from the VEP annotation
ANNOVARpharmGKBandVEPdata %>%
  mutate(
    .data = .,
    # change empty CADD scores to an NA
    CADD_PHRED_VEP = if_else(
                            condition = CADD_PHRED_VEP == "-",
                            true = as.character(x = NA),
                            false = as.character(x = CADD_PHRED_VEP)
                            )
  ) %>%
  # make sure the phred-scaled CADD score is a numeric
  mutate(
    .data = .,
    CADD_PHRED_VEP = as.numeric(x = CADD_PHRED_VEP)
    ) %>%
  filter(
    .data = .,
    CADD_PHRED_VEP >= 20
    ) 

# print the number of variants with a Phred-scaled CADD score of at least 20
print(
  x = paste0(
            "################# There are ",
            nrow(x = CADDphredGreaterThan20Table),
            " variants with a phred-scaled CADD score of at least 20 ##########################"
            )
  )

# load the file of ANNOVAR and pharmGKB annotated singleton data to count the number of singletons that are present
SingletonData <-
"./annotation_output/singletons_annovar_pharmGKB_and_VEP_data.tsv" %>%
  read.delim(
    file = .,
    header = TRUE
    )  
  
NumberOfSingletons <-
  SingletonData %>%
  nrow(x = .)

  NumberOfSingletons %>%
  paste0("%%%%%%%%%%%%%%%%%%%%% There are ",.," total singletons %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%") %>%
  print(x = .)
  
# load the file of ANNOVAR and pharmGKB annotated INDEL data to count the number of INDELS that are present
IndelData <-
  "./annotation_output/indels_annovar_pharmGKB_and_VEP_data.tsv" %>%
  read.delim(
    file = .,
    header = TRUE
    )
  
NumberOfIndels <-
  IndelData %>%
  nrow(x = .)

  NumberOfIndels %>%
    paste0("%%%%%%%%%%%%%%%%%%%%% There are ",.," total INDELS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%") %>%
    print(x = .)
  
# create a table with number of total variants (known and novel), novel variants, known variants, 
# variants minor allele frequency (MAF) < 0.05, novel AND minor allele frequency (MAF) < 0.05 variants, variants with minor allele count (MAC) ≥ 5, 
# variants with minor allele count (MAC) < 5 and MAC >= 1, variants iwth minor allele count (MAC) = 0, number of INDELS, and number of singletons
VariantCountTable <-
  tibble(
    "Category" = c(
                  "TotalVariants",
                  "NovelVariants",
                  "KnownVariants",
                  "Variants with minor allele frequency (MAF) < 0.05",
                  "Novel variants with minor allele frequency (MAF) < 0.05",
                  "Variants with minor allele count (MAC) ≥ 5",
                  "Variants with minor allele count (MAC) < 5 and minor allele count (MAC) >= 1",
                  "Variants with minor allele count (MAC) = 0",
                  "NovelANDexonic",
                  "Singletons",
                  "INDELs"
                  ),
    "VariantCount" = c(
                       nrow(x = ANNOVARpharmGKBandVEPdata),
                       NumberOfNovelVariants,
                       NumberOfKnownVariants,
                       VariantsWithMAFlessThan0.05,
                       VariantsThatAreNovelAndHaveMAFlessThan0.05,
                       variantsWithMACofAtLeast5,
                       variantsWithMACofAtMost4,
                       variantsWithMACofZero,
                       NumberOfNovelExonicVariants,
                       NumberOfSingletons,
                       NumberOfIndels
                       )
  )

# save the table of variant counts to the file system
VariantCountTable %>%
  write_tsv(
    x = .,
    file = "./ClinicalVariantsForPresentations/VariantCounts.txt"
    )

print(x = "########################## Variant count table computed ################################")
print(x = VariantCountTable)

### create a table of all clinically relevant variants:
    # identify variants in the merged dataset that have associated, pathogenic, likely pathogenic, or drug response phenotypes in clinVar
    # OR a pharmGKB annotation with a non-missing level of evidence
    # OR a CADD score of at least 20
    # OR an ADME optimized prediction score > 0.5 for missense SNVs only ("missense_variant" or "missense_variant,splice_region_variant") 
AllUniqueClinicalVariants <-
    ANNOVARpharmGKBandVEPdata %>%
      mutate(
        .data = .,
        # change empty CADD scores to an NA
        CADD_PHRED_VEP = if_else(
                                condition = CADD_PHRED_VEP == "-",
                                true = as.character(x = NA),
                                false = as.character(x = CADD_PHRED_VEP)
                                  ),
        # make sure the ADME optimized prediction score is a numeric
        ADME_optimized_prediction = ADME_optimized_prediction %>% as.numeric(x = .)
        ) %>%
      # change the CADD scores to a numeric
      mutate(
        .data = .,
        CADD_PHRED_VEP = CADD_PHRED_VEP %>% as.numeric(x = .)
        ) %>%
      # identify variants in the merged dataset that have  associated, pathogenic, likely pathogenic, or drug response phenotypes in clinVar
      # OR a pharmGKB clinical annnotation with a non-missing level of evidence
      # OR a CADD score of at least 20
      # OR an ADME optimized prediction score > 0.5 for missense SNVs ("missense_variant" or "missense_variant,splice_region_variant")
      filter(
        .data = .,
          (CLNSIG == "Pathogenic") | 
          (CLNSIG == "Likely_pathogenic") | 
          (CLNSIG == "Pathogenic/Likely_pathogenic") |
          (CLNSIG == "Pathogenic/Likely_pathogenic,_other") |
          (CLNSIG == "drug_response") |
          (CLNSIG == "association") |
          !is.na(x = pharmGKB_EvidenceLevel) |
          (CADD_PHRED_VEP >= 20) |
          (
            (ADME_optimized_prediction > 0.5) & 
            ( Consequence_VEP == "missense_variant" | Consequence_VEP == "missense_variant,splice_region_variant" )
            )
        ) %>%
      # convert all table entries with a "." to a "-" or an NA to a "-"
      mutate(
        .data = .,
        across(
               .cols = c(
                        starAllele,
                        Chr,
                        Start,
                        CLNDN,
                        CLNSIG,
                        X1000g2015aug_all,
                        X1000g2015aug_afr,
                        X1000g2015aug_amr,
                        X1000g2015aug_eas,
                        X1000g2015aug_eur,
                        X1000g2015aug_sas,
                        ADME_optimized_prediction,
                        pharmGKB_EvidenceLevel
                         ),
               .fns = ~ case_when(
                                  (.x == ".") | (is.na(x = .x)) ~ "-",
                                  TRUE ~ as.character(x = .x)
                                 )
                 )
        ) %>%
      # replace the commas separating the pubmed IDs with semi-colons so they can be viewed with Microsoft Excel
      mutate(
        .data = .,
        Pubmed_ID = Pubmed_ID %>%
                    gsub(
                      pattern = ",",
                      replacement = ";",
                      x = .
                      )
        ) %>%
      # remove any cosmic variant IDs (COSV) or CM variant IDs if there is an rs-number present
      # if there is no rs-id in the existing_variant_VEP, change it to a "*rsNA" for no rs-number
      RemoveAlternateVariantIDsAndReplaceMissingRSnumbers(
        InputDataFrame = .,
        rsIDcolumnName = "existing_variant_VEP"
      ) %>%
      # select and rename useful columns
      select(
        .data = .,
        "Variant ID" = existing_variant_VEP,
        "Star Allele" = starAllele,
        Chr,
        Start,
        "Nucleotide" = NucleotideID,
        "cDNA" = HGVSc,
        "Exon" = exon_number_VEP,
        "Intron" = intron_number_VEP,
        "Gene" = Gene_VEP,
        "Consequence" = Consequence_VEP,
        "Amino Acid" = HGVSp,
        "CADD (Phred scale)" = CADD_PHRED_VEP,
        "ADME Framework Score" = ADME_optimized_prediction,
        "PharmGKB Evidence" = pharmGKB_EvidenceLevel,
        "AAF**" = CSKT_AF,
        "ALL" = X1000g2015aug_all,
        "AFR" = X1000g2015aug_afr,
        "AMR" = X1000g2015aug_amr,
        "EAS" = X1000g2015aug_eas,
        "EUR" = X1000g2015aug_eur,
        "SAS" = X1000g2015aug_sas,
        "ClinVar Descript." = CLNDN,
        "ClinVar Signif." = CLNSIG,
        "Pubmed ID" = Pubmed_ID
        )

print(
  x = paste0(
             "There are ",
             nrow(x = AllUniqueClinicalVariants),
             " total variants with a ClinVar associated, pathogenic, likely pathogenic, or drug response phenotypes a pharmGKB clinical annotation with non-missing level of evidence or CADD score >= 20 or an ADME optimized prediction score > 0.5"
             )
  )

# save the table of clinical variants to the file system
AllUniqueClinicalVariants %>%
  write_tsv(
    x = .,
    file = "./ClinicalVariantsForPresentations/AllUniqueClinicalVariants.txt"
      )


##### create a table with clinical annotation tallies for the categories computed above
ClinicalAnnotationTallyTable <-
# create a table with variant tallies for each clinical variant annotation category that were used to create the table directly above
AllUniqueClinicalVariants %>%
  # create a column for identifying whether or not a variant has a associated, pathogenic, likely pathogenic, or drug response phenotypes in clinVar
  mutate(
    .data = .,
    ClinVarTally = (                     
                    (`ClinVar Signif.` == "Pathogenic") |
                    (`ClinVar Signif.` == "Likely_pathogenic") |
                    (`ClinVar Signif.` == "Pathogenic/Likely_pathogenic") |
                    (`ClinVar Signif.` == "Pathogenic/Likely_pathogenic,_other") |
                    (`ClinVar Signif.` == "drug_response") |
                    (`ClinVar Signif.` == "association")
                    )
    ) %>%
  # create a column for identifying whether or not a variant has a clinical pharmGKB annotation with non-missing level of evidence
  mutate(
    .data = .,
    PharmGKBClinicalTally = (`PharmGKB Evidence` != "-")
    ) %>%
  # create a column for identifying whether or not a missense variant has an ADME optimized prediction score greater than 0.5
  mutate(
    .data = .,
    ADMEpredictionTally =  (
                             (`ADME Framework Score` > 0.5) & 
                             ( Consequence == "missense_variant" | Consequence == "missense_variant,splice_region_variant" )
                            )
    ) %>%
  # create a column for identifying whether or not a variant has a phred-scaled CADD score of at least 20
  mutate(
    .data = .,
    CADDtally = `CADD (Phred scale)` >= 20
    ) %>%
  select(
    .data = .,
    ClinVarTally,
    PharmGKBClinicalTally,
    ADMEpredictionTally,
    CADDtally
    ) %>%
  colSums(
    x = .,
    na.rm = TRUE
    )
  print(x = "############### Clinical and in silico prediction annotation tallies computed ###################")
  print(x = ClinicalAnnotationTallyTable)
  
  # save the ClinicalAnnotationTallyTable to the file system
  ClinicalAnnotationTallyTable %>%
    enframe(
      x = .,
      name = "Category",
      value = "Variant count"
        ) %>%
    write_tsv(
      x = .,
      file = "./ClinicalVariantsForPresentations/ClinicalAnnotationTallyTable.txt"
        )

  
########################################################
# create a table with general variant tallies grouped by gene:
# for each gene, create a table with tallies of:

      # (1) total variants 
      # (2) number that are novel (no rsID or existing variant ID from dbSNP from ANNOVAR or VEP annotations) 
      # (3) number in  non-exonic functional classes (intronic; 3_prime_UTR_variant; 5_prime_UTR_variant; downstream_gene_variant; upstream_gene_variant; intronic splicing; splice acceptor; splice donor) 
      # (4) number that are exonic (includes frameshift_variant; frameshift_variant,splice_region_variant; inframe_deletion; missense_variant; missense_variant,splice_region_variant; splice_region_variant,synonymous_variant; start_lost; stop_gained; synonymous_variant)
      # (5) number in each exonic functional class  
      # (6) number that have minor allele frequency (MAF) < 0.05 
      # (7) number that are novel & have minor allele frequency (MAF) < 0.05 
      # (8) number that have minor allele frequency (MAF) < 0.01
      # (9) number that are novel & have minor allele frequency (MAF) < 0.01
    AllvariantTalliesPerGene <-
        ANNOVARpharmGKBandVEPdata %>%
        ###### create a column with a 1=yes if the variant is novel, otherwise 0=no
          mutate(
            .data = .,
            NovelVariant =  ( 
                              (avsnp142 == ".") & 
                              (avsnp138 == ".") & 
                              (snp138 == ".") & 
                              (existing_variant_VEP == "-") 
                             ),
            ###### create a column with a 1=yes for the main variant functional classes:
             # exonic (includes frameshift_variant; frameshift_variant,splice_region_variant; inframe_deletion; missense_variant; missense_variant,splice_region_variant; splice_region_variant,synonymous_variant; start_lost; stop_gained; synonymous_variant)
            Exonic = (
                       (Consequence_VEP == "frameshift_variant") |
                       (Consequence_VEP == "frameshift_variant,splice_region_variant") |
                       (Consequence_VEP == "inframe_deletion") |
                       (Consequence_VEP == "missense_variant") |
                       (Consequence_VEP == "missense_variant,splice_region_variant") |
                       (Consequence_VEP == "splice_region_variant,synonymous_variant") |
                       (Consequence_VEP == "start_lost") |
                       (Consequence_VEP == "stop_gained") |
                       (Consequence_VEP == "synonymous_variant")
                     ),
            # exonic splicing
            ExonicSplicing = (
                              (Consequence_VEP == "frameshift_variant,splice_region_variant") |
                              (Consequence_VEP == "missense_variant,splice_region_variant") |
                              (Consequence_VEP == "splice_region_variant,synonymous_variant")
                              ),
            # intronic (includes variants in intronic regions that influence splicing)
            Intronic = (
                        (Consequence_VEP == "intron_variant") | 
                        (Consequence_VEP == "splice_region_variant,intron_variant") |
                        (Consequence_VEP == "splice_acceptor_variant") |
                        (Consequence_VEP == "splice_donor_variant") |
                        (Consequence_VEP == "splice_donor_5th_base_variant,intron_variant") |
                        (Consequence_VEP == "splice_donor_region_variant,intron_variant") |
                        (Consequence_VEP == "splice_polypyrimidine_tract_variant,intron_variant") |
                        (Consequence_VEP == "splice_region_variant,splice_polypyrimidine_tract_variant,intron_variant") 
                        ),
            # intronic splicing
            IntronicSplicing = (
                               (Consequence_VEP == "splice_acceptor_variant") | 
                               (Consequence_VEP == "splice_donor_variant") |
                               (Consequence_VEP == "splice_region_variant,intron_variant") |
                               (Consequence_VEP == "splice_donor_5th_base_variant,intron_variant") |
                               (Consequence_VEP == "splice_donor_region_variant,intron_variant") |
                               (Consequence_VEP == "splice_polypyrimidine_tract_variant,intron_variant") |
                               (Consequence_VEP == "splice_region_variant,splice_polypyrimidine_tract_variant,intron_variant") 
                               ),
            # UTR3
            UTR3 = Consequence_VEP == "3_prime_UTR_variant",
            # UTR5
            UTR5 = Consequence_VEP == "5_prime_UTR_variant",
            # downstream gene variant
            DownstreamGene = Consequence_VEP == "downstream_gene_variant",
            # upstream gene variant
            UpstreamGene = Consequence_VEP == "upstream_gene_variant",
            ###### create a column with a 1=yes for the EXONIC variant functional classes:
              # frameshift_variant (includes frameshift variants in splicing regions)
            frameshift_variant = ( 
                                  (Consequence_VEP == "frameshift_variant") | 
                                  (Consequence_VEP == "frameshift_variant,splice_region_variant")
                                  ),
              # frameshift_splicing
            frameshift_splicing = Consequence_VEP == "frameshift_variant,splice_region_variant",
            # inframe_deletion 
            inframe_deletion = Consequence_VEP == "inframe_deletion",
            # missense_variant (includes missense variants in splicing regions)
            missense_variant = (
                                 (Consequence_VEP == "missense_variant") | 
                                 (Consequence_VEP == "missense_variant,splice_region_variant")
                                ),
            # missense_splicing 
            missense_splicing = Consequence_VEP == "missense_variant,splice_region_variant",
            # start_lost
            start_lost = Consequence_VEP == "start_lost",
            # stop_gained
            stop_gained = Consequence_VEP == "stop_gained",
            # synonymous_variant (includes synonymous variants in splice regions)
            synonymous_variant = (
                                  (Consequence_VEP == "synonymous_variant") | 
                                  (Consequence_VEP == "splice_region_variant,synonymous_variant")
                                  ),
            # synonymous_splicing
            synonymous_splicing = Consequence_VEP == "splice_region_variant,synonymous_variant",
            ###### create a column with a 1=yes if the minor allele frequency (MAF) is less than 5 %
            MAFbelow0.05 = (CSKT_AF < 0.05) | ((1 - CSKT_AF) < 0.05),
            ##### create a column with a 1=yes if the minor allele frequency (MAF) is less than 1 %
            MAFbelow0.01 = (CSKT_AF < 0.01) | ((1 - CSKT_AF) < 0.01),
            ###### create a column with a 1=yes if the variant is novel AND the minor allele frequency (MAF) is less than 5 %
            Novel_and_MAFbelow0.05 =  (
                                        ( 
                                          (avsnp142 == ".") & 
                                          (avsnp138 == ".") & 
                                          (snp138 == ".") & 
                                          (existing_variant_VEP == "-") 
                                        ) & 
                                        (
                                          (CSKT_AF < 0.05) | 
                                          ((1 - CSKT_AF) < 0.05)
                                        )
                                      ),
            ###### create a column with a 1=yes if the variant is novel AND the minor allele frequency (MAF) is less than 1 %
            Novel_and_MAFbelow0.01 = (
                                        (
                                          (avsnp142 == ".") & 
                                          (avsnp138 == ".") & 
                                          (snp138 == ".") & 
                                          (existing_variant_VEP == "-")
                                        ) & 
                                        (
                                          (CSKT_AF < 0.01) | 
                                          ((1 - CSKT_AF) < 0.01)
                                        )
                                      ),
          ) %>%
        # if the gene column contains the string "UGT1A" change the gene label to UGT1A
        mutate(
          .data = .,
          Gene_VEP = if_else(
                             condition = grepl(
                                               pattern = "UGT1A",
                                               x = Gene_VEP
                                               ),
                             true = "UGT1A",
                             false = Gene_VEP
                               )
          ) %>%
        # group by the gene
        group_by(
          .data = .,
          Gene_VEP
          ) %>%
        summarise(
          .data = .,
          # count total variants per gene
          "Total_Count" = n(),
          # count variant functional classes per gene
          "Intronic_Count" = sum(Intronic,na.rm = TRUE),
          "IntronicSplicing_Count" = sum(IntronicSplicing,na.rm = TRUE),
          "UTR3_Count" = sum(UTR3,na.rm = TRUE),
          "UTR5_Count" = sum(UTR5,na.rm = TRUE),
          "DownstreamGene_Count" = sum(DownstreamGene,na.rm = TRUE),
          "UpstreamGene_Count" = sum(UpstreamGene,na.rm = TRUE),
          "Exonic_Count" = sum(Exonic,na.rm = TRUE),
          "ExonicSplicing_Count" = sum(ExonicSplicing,na.rm = TRUE),
          "frameshift_variant_Count" = sum(frameshift_variant,na.rm = TRUE),
          "frameshift_splicing_Count" = sum(frameshift_splicing,na.rm = TRUE),
          "inframe_deletion_Count" = sum(inframe_deletion,na.rm = TRUE),
          "missense_variant_Count" = sum(missense_variant,na.rm = TRUE),
          "missense_splicing_Count" = sum(missense_splicing,na.rm = TRUE),
          "start_lost_Count" = sum(start_lost,na.rm = TRUE),
          "stop_gained_Count" = sum(stop_gained,na.rm = TRUE),
          "synonymous_variant_Count" = sum(synonymous_variant,na.rm = TRUE),
          "synonymous_splicing_Count" = sum(synonymous_splicing,na.rm = TRUE),
          # count the number of novel variants per gene by taking the column sum of the novel variants
          "Novel_Count" = sum(NovelVariant,na.rm = TRUE),
          # count the number of variants per gene with minor allele frequency (MAF) < 0.05
          "MAFbelow0.05_Count" = sum(MAFbelow0.05,na.rm = TRUE),
          # count the number of variants per gene with minor allele frequency (MAF) < 0.01
          "MAFbelow0.01_Count" = sum(MAFbelow0.01,na.rm = TRUE),
          # count the number of variants that are novel AND have minor allele frequency (MAF) < 0.05
          "MAFbelow0.05andNovel_Count" = sum(Novel_and_MAFbelow0.05,na.rm = TRUE),
          # count the number of variants that are novel AND have minor allele frequency (MAF) < 0.01
          "Novel_and_MAFbelow0.01_Count" = sum(Novel_and_MAFbelow0.01,na.rm = TRUE)
          )
      
      
      # (8) number of singletons
      # (9) number that are novel & singleton
      # (10) number of exonic singletons
      
      SingletonVariantTalliesPerGene <-
        SingletonData %>%
        mutate(
          .data = .,
          # create a column with a 1=yes if the singleton is novel (has no rsID or existing variant ID from ANNOVAR or VEP annotations)
          NovelSingleton = (
                              (avsnp142 == ".") & 
                              (avsnp138 == ".") & 
                              (snp138 == ".") & 
                              (existing_variant_VEP == "-")
                            ),
          # create a column with 1=yes if the singleton is exonic 
          ExonicSingleton = (
                              (Consequence_VEP == "frameshift_variant") |
                              (Consequence_VEP == "frameshift_variant,splice_region_variant") |
                              (Consequence_VEP == "inframe_deletion") |
                              (Consequence_VEP == "missense_variant") |
                              (Consequence_VEP == "missense_variant,splice_region_variant") |
                              (Consequence_VEP == "splice_region_variant,synonymous_variant") |
                              (Consequence_VEP == "start_lost") |
                              (Consequence_VEP == "stop_gained") |
                              (Consequence_VEP == "synonymous_variant")
                            )
          ) %>%
        # if the gene column contains the string "UGT1A" change the gene label to UGT1A
        mutate(
          .data = .,
          Gene_VEP = if_else(
                             condition = grepl(
                                               pattern = "UGT1A",
                                               x = Gene_VEP
                                               ),
                             true = "UGT1A",
                             false = Gene_VEP
                               )
          ) %>%
        # group by the gene
        group_by(
          .data = .,
          Gene_VEP
          ) %>%
        summarise(
          .data = .,
          # count the number of singletons per gene
          "TotalSingleton_Count" = n(),
          # count the number of novel singletons per gene
          "NovelSingleton_Count" = sum(NovelSingleton,na.rm = TRUE),
          # count the number of exonic singletons
          "ExonicSingleton_Count" = sum(ExonicSingleton,na.rm = TRUE)
          )
      
      # (11) number of INDELs
      # (12) number of exonic INDELs
      # (13) number of singleton INDELs
      
      IndelTalliesPerGene <-
        IndelData %>%
        mutate(
          .data = .,
          # create a column with a 1=yes if the INDEL is exonic 
          ExonicINDEL = (
                          (Consequence_VEP == "frameshift_variant") |
                          (Consequence_VEP == "frameshift_variant,splice_region_variant") |
                          (Consequence_VEP == "inframe_deletion") |
                          (Consequence_VEP == "missense_variant") |
                          (Consequence_VEP == "missense_variant,splice_region_variant") |
                          (Consequence_VEP == "splice_region_variant,synonymous_variant") |
                          (Consequence_VEP == "start_lost") |
                          (Consequence_VEP == "stop_gained") |
                          (Consequence_VEP == "synonymous_variant")
                        ),
          # create a column with a 1=yes if the INDEL is a singleton [ minor allele frequency == 0.001066 == (1/(469*2)) ]
          INDEL_and_singleton = (CSKT_AF > 0) & 
                                (
                                 (CSKT_AF <= (1/(469*2))) |
                                 ((1 - CSKT_AF) <= (1/(469*2)))
                                )
          ) %>%
        # if the gene column contains the string "UGT1A" change the gene label to UGT1A
        mutate(
          .data = .,
          Gene_VEP = if_else(
                             condition = grepl(
                                               pattern = "UGT1A",
                                               x = Gene_VEP
                                               ),
                             true = "UGT1A",
                             false = Gene_VEP
                               )
          ) %>%
        # group by the gene
        group_by(
          .data = .,
          Gene_VEP
        ) %>%
        summarise(
          .data = .,
          # count the number of INDELs per gene
          "TotalINDEL_Count" = n(),
          # count the number of exonic INDELS 
          "ExonicINDEL_Count" = sum(ExonicINDEL,na.rm = TRUE),
          # count the number of singleton INDELS
          "INDELandSingleton_Count" = sum(INDEL_and_singleton,na.rm = TRUE)
          )
      
    # join the tables of variant tallies for all genes, singletons, and indels together by gene 
    # and save the final table to the file system
    VariantTallyTablePerGene <-
      AllvariantTalliesPerGene %>%
      left_join(
        x = .,
        y = SingletonVariantTalliesPerGene,
        by = "Gene_VEP"
          ) %>%
      left_join(
        x = .,
        y = IndelTalliesPerGene,
        by = "Gene_VEP"
        ) %>%
      # order the columns as they will be specified in the manuscript for this data
      select(
        .data = .,
        "Gene" = Gene_VEP,
        "Variants" = Total_Count,
        "Novel" = Novel_Count,
        "*MAF < 0.05" = MAFbelow0.05_Count,
        "Novel and *MAF < 0.05" = MAFbelow0.05andNovel_Count,
        "*MAF < 0.01" = MAFbelow0.01_Count,
        "Novel and *MAF < 0.01" = Novel_and_MAFbelow0.01_Count,
        "Singleton" = TotalSingleton_Count,
        "Novel singleton" = NovelSingleton_Count,
        "Exonic singleton" = ExonicSingleton_Count,
        "INDEL" = TotalINDEL_Count,
        "Exonic INDEL" = ExonicINDEL_Count,
        "Singleton INDEL" = INDELandSingleton_Count,
        "Intronic" = Intronic_Count,
        "Intronic splicing" = IntronicSplicing_Count,
        "UTR3" = UTR3_Count,
        "UTR5" = UTR5_Count,
        "Downstream gene" = DownstreamGene_Count,
        "Upstream gene" = UpstreamGene_Count,
        "Exonic" = Exonic_Count,
        "Exonic splicing" = ExonicSplicing_Count,
        "Frameshift" = frameshift_variant_Count,
        "Frameshift splicing" = frameshift_splicing_Count,
        "Inframe del." = inframe_deletion_Count,
        "Missense" = missense_variant_Count,
        "Missense splicing" = missense_splicing_Count,
        "Start lost" = start_lost_Count,
        "Stop gained" = stop_gained_Count,
        "Synonymous" = synonymous_variant_Count,
        "Synonymous splicing" = synonymous_splicing_Count
        ) %>%
      # add column totals to the last row of the table
      bind_rows(
                summarise(
                          .data = .,
                          across(where(is.numeric), ~ sum(.x,na.rm = TRUE)),
                          across(where(is.character), ~"Total")
                          )
                ) 
    
    VariantTallyTablePerGene %>%
      write_tsv(
        x = .,
        file = "./ClinicalVariantsForPresentations/VariantTallyTablePerGene.txt"
          )
    
    print(x = "################################# Variant tallies per gene computed ##################################")
    print(x = VariantTallyTablePerGene)
    
# create a table of all variants that are either novel or exonic to be included as a supplementary table in the vitamin D publication
# the table is formatted analogous to the table in Fohner, A., et al. (2013). Pharmacogenetics in American Indian populations: analysis of CYP2D6, CYP3A4, CYP3A5, and CYP2C9 in the Confederated Salish and Kootenai Tribes. Pharmacogenetics and genomics, 23(8), 403–414. https://doi.org/10.1097/FPC.0b013e3283629ce9
ExonicORnovelVariants <-
  ANNOVARpharmGKBandVEPdata %>%
  # filter to variants that are either novel (no rs-identification number or existing variant ID) or exonic 
  filter(
    .data = .,
    ((avsnp142 == ".") & (avsnp138 == ".") & (snp138 == ".") & (existing_variant_VEP == "-")) | 
    (
      (Consequence_VEP == "frameshift_variant") |
      (Consequence_VEP == "frameshift_variant,splice_region_variant") |
      (Consequence_VEP == "inframe_deletion") |
      (Consequence_VEP == "missense_variant") |
      (Consequence_VEP == "missense_variant,splice_region_variant") |
      (Consequence_VEP == "splice_region_variant,synonymous_variant") |
      (Consequence_VEP == "start_lost") |
      (Consequence_VEP == "stop_gained") |
      (Consequence_VEP == "synonymous_variant")
    )
    ) %>%
  # remove any cosmic variant IDs (COSV) or CM variant IDs if there is an rs-number present
  # if there is no rs-id in the existing_variant_VEP, change it to a "*rsNA" for no rs-number
  RemoveAlternateVariantIDsAndReplaceMissingRSnumbers(
    InputDataFrame = .,
    rsIDcolumnName = "existing_variant_VEP"
  ) %>%
  # rename and select relevant columns
  select(
    .data = .,
    "Variant ID" = existing_variant_VEP,
    "Star Allele" = starAllele,
    Chr,
    Start,
    "Nucleotide" = NucleotideID,
    "cDNA" = HGVSc,
    "Exon" = exon_number_VEP,
    "Intron" = intron_number_VEP,
    "Gene" = Gene_VEP,
    "Consequence" = Consequence_VEP,
    "Amino Acid" = HGVSp,
    "CADD (Phred scale)" = CADD_PHRED_VEP,
    "ADME Framework Score" = ADME_optimized_prediction,
    "PharmGKB Evidence" = pharmGKB_EvidenceLevel,
    "AAF**" = CSKT_AF,
    "ALL" = X1000g2015aug_all,
    "AFR" = X1000g2015aug_afr,
    "AMR" = X1000g2015aug_amr,
    "EAS" = X1000g2015aug_eas,
    "EUR" = X1000g2015aug_eur,
    "SAS" = X1000g2015aug_sas,
    "ClinVar Descript." = CLNDN,
    "ClinVar Signif." = CLNSIG,
    "Pubmed ID" = Pubmed_ID
    ) %>%
  # convert all table entries with a "." to a "-" or an NA to a "-"
  mutate(
    .data = .,
    across(
            .cols = c(
                      `Variant ID`,
                      `Star Allele`,
                      Chr,
                      Start,
                      Nucleotide,
                      cDNA,
                      Gene,
                      Consequence,
                      `Amino Acid`,
                      `ClinVar Descript.`,
                      `ClinVar Signif.`,
                      `ALL`,
                      `AFR`,
                      `AMR`,
                      `EAS`,
                      `EUR`,
                      `SAS`,
                      `ADME Framework Score`,
                      `PharmGKB Evidence`
                      ),
            .fns = ~ case_when(
                               (.x == ".") | (is.na(x = .x)) ~ "-",
                                TRUE ~ as.character(x = .x)
                              )
    )
  )

# save the exonic/novel variant table to the file system
ExonicORnovelVariants %>%
  write_tsv(
    x = .,
    file = "./ClinicalVariantsForPresentations/ExonicOrNovelVariantsForPaper.txt"
    )
 

# combine the clinical variants table and the exonic/novel variant table and save to the file system
  rbind(
    AllUniqueClinicalVariants,
    ExonicORnovelVariants
    ) %>%
    # obtain unique rows only based on the variant ID
    distinct(
      .data = .,
      Nucleotide,
      .keep_all = TRUE
      ) %>%
    # # obtain unique rows only
    # unique(x = .) %>%
    # sort by chromosome number and nucleotide position in ascending order
    arrange(
      .data = .,
      as.numeric(x = Chr),
      as.numeric(x = Start)
      ) %>%
    # deselect the Chr and Start columns 
    select(
      .data = .,
      -Chr,
      -Start
      ) %>%
    # save the results to the file system
    write_tsv(
      x = .,
      file = "./ClinicalVariantsForPresentations/AllUniqueClinicalVariantsAndExonicOrNovelVariantsForPaper.txt"
      )
    
