# this is a script for combining ANNOVAR and pharmGKB annotated data with useful data from the VEP annotation output
# this script adds pharmGKB levels of evidence to the ANNOVAR annotated data based on rs-number
# this script also adds any star alleles to the annotation data for cyp2r1, cyp3a4, ugt1a1, and ugt1a4 based on a table of rs-numbers and star alleles that was manually created

# load necessary packages
source("../scripts/load_R_packages.R")
source("../scripts/genetic_characterization/genetic_characterization_functions.R")
# load the ANNOVAR annotated data with ADME optimized prediction framework scores computed
ANNOVARfileNames <-
  list(
    "allVariants" = "../annotation_output/merge_annovar_and_ADMEoptimizedPredictions.tsv",
    "singletons" = "../singleton_annotation_output/merge_singleton_annovar_and_ADMEoptimizedPredictions.tsv",
    "indels" = "../indel_annotation_output/merge_indel_annovar_and_ADMEoptimizedPredictions.tsv"
  )

ANNOVARdata <-
ANNOVARfileNames %>%
  lapply(
    X = .,
    FUN = function(currentFile)
    {
      LoadedFile <-
        read.delim(
          file = currentFile,
          header = TRUE
          )
      
      return(LoadedFile)
    }
      )

# load the VEP annotation data
VEPfileNames <-
  list(
    "allVariants" = "../VEP_annotation_output/allVariants_VEP_results.txt",
    "singletons" = "../VEP_annotation_output/singletons_VEP_results.txt",
    "indels" = "../VEP_annotation_output/indels_VEP_results.txt"
  )

VEPannotationData <-
  VEPfileNames %>%
  lapply(
    X = .,
    FUN = function(currentFile)
    {
      
      # load the VEP file
      LoadedFile <-
      read.delim(
        file = currentFile
        ) %>%
        # filter to refSeq transcripts ONLY
        # also filter to the relevant gene symbols that were used in this study
        filter(
          .data = .,
          (SOURCE == "RefSeq") &
            (
              (SYMBOL == "CASR") |
              (SYMBOL == "CUBN") |
              (SYMBOL == "CYP24A1") |
              (SYMBOL == "CYP27B1") |
              (SYMBOL == "CYP2R1") |
              (SYMBOL == "CYP3A4") |
              (SYMBOL == "DHCR7") |
              (SYMBOL == "GC") |
              (SYMBOL == "RXRA") |
              (SYMBOL == "RXRB") |
              (SYMBOL == "RXRG") |
              (SYMBOL == "SULT2A1") |
                # all ugt1a family genes except for the UGT1A2 pseudo-gene
              (grepl(pattern = "UGT1A",x = SYMBOL)==TRUE & grepl(pattern = "UGT1A2P",x = SYMBOL)==FALSE) |
              (SYMBOL == "VDR")
            )
          ) %>%
        # select useful columns
        select(
          .data = .,
          "VariantID_VEP" = X.Uploaded_variation,
          "Location_VEP" = Location,
          GIVEN_REF,
          USED_REF,
          "Allele_VEP" = Allele,
          "Consequence_VEP" = Consequence,
          "IMPACT_rating_VEP" = IMPACT,
          "Gene_VEP" = SYMBOL,
          "Gene_code_VEP" = Gene,
          Feature_type:BIOTYPE,
          "exon_number_VEP" = EXON,
          "intron_number_VEP" = INTRON,
          HGVSc:CDS_position,
          "Protein_position_VEP" = Protein_position,
          "Amino_acids_VEP" = Amino_acids,
          Codons,
          "existing_variant_VEP" = Existing_variation,
          STRAND,
          SOURCE,
          "VAR_SYNONYMS_VEP" = VAR_SYNONYMS,
          "Pubmed_ID" = PUBMED,
          "PHENOTYPES_VEP" = PHENOTYPES,
          MOTIF_NAME,
          DOMAINS,
          "TRANSCRIPTION_FACTORS_VEP" = TRANSCRIPTION_FACTORS,
          "CADD_RAW_VEP" = CADD_RAW,
          "CADD_PHRED_VEP" = CADD_PHRED#,
          #"LoFtool_VEP" = LoFtool
          ) %>%
        # obtain one row for each variant ID
        distinct(
          .data = .,
          VariantID_VEP,
          .keep_all = TRUE
          ) %>%
        # create a duplicate copy of the existing variant ID column
        mutate(
          .data = .,
          rsNumber_VEP = existing_variant_VEP
          ) %>%
        # but only with the rs-numbers in the existing variant ID column, if an rs-number exists
        RemoveAlternateVariantIDsAndReplaceMissingRSnumbers(
          InputDataFrame = .,
          rsIDcolumnName = "rsNumber_VEP"
          )
    }
  )

# create a file of all novel variants (variants with no rsIDs or existing variant ID) and reformat the novel variants to a format that is compatible with SNP-nexus (https://www.snp-nexus.org/v4/)
NovelVariants <-
ANNOVARdata$allVariants %>%
  mutate(
    .data = .,
    VariantID_fromInputFile_copy = VariantID_fromInputFile
  ) %>%
  left_join(
    x = .,
    y = VEPannotationData$allVariants,
    by = c("VariantID_fromInputFile" = "VariantID_VEP")
      ) %>%
  select(
    .data = .,
    avsnp142,
    avsnp138,
    snp138,
    existing_variant_VEP,
    VariantID_fromInputFile
  ) %>%
  # filter to novel variants only
  filter(
    .data = .,
    (avsnp142 == ".") & (avsnp138 == ".") & (snp138 == ".") & (existing_variant_VEP == "-")
  ) %>%
  # select the snp ID only
  select(
    .data = .,
    "SNP" = VariantID_fromInputFile
  ) %>%
  # convert to a tibble
  tibble() %>%
  # split the snp ID into a list on the ":"
  mutate(
    .data = .,
    SNP = SNP %>% str_split(string = .,pattern = ":")
  ) %>%
  # pull the column out to create a list of snp IDs
  pull(.data = .,SNP) %>%
  # create a dataframe with each split snpID as a row in the dataframe
  do.call(
    what = "rbind",
    args = .
  ) %>%
  as.data.frame(
    x = .,
    row.names = NULL
  ) %>%
  # create a chromosome column,
  # and replace "*" with "-" for indels,
  # and add a strand column with 1 for the plus strand
  mutate(
    .data = .,
    chrom = "Chromosome",
    chrom_number = V1,
    pos = V2,
    ref = if_else(
      condition = V3 == "*",
      true = "-",
      false = V3
    ),
    alt = if_else(
      condition = V4 == "*",
      true = "-",
      false = V4
    ),
    strand = 1
  ) %>%
  # select relevant columns
  select(
    .data = .,
    chrom,
    chrom_number,
    pos,
    ref,
    alt,
    strand
  )

# create a directory for saving the novel variants if it does not exist
if(!dir.exists(paths = "../novel_variants_for_snp_nexus_annotation/"))
{
  dir.create(path = "../novel_variants_for_snp_nexus_annotation/")
}

# save the resulting file to the file system
write_tsv(
  x = NovelVariants,
  file = "../novel_variants_for_snp_nexus_annotation/novel_variants.txt",
  col_names = FALSE
)

# (1) annotate the novel variants file with SNP-nexus using the GRCh37 genome build manually on the snp-nexus website: https://www.snp-nexus.org/v4/ with ALL annotation categories selected
# (2) download and save the resulting annotation text file from the SNP-nexus interface to the ../novel_variants_for_snp_nexus_annotation/ directory and name as SNP_nexus_annotation.txt

# identify variants that have rs numbers from the snp-nexus annotation
SNPnexusRSids <-
  # load the SNP-nexus annotation file
  read.delim(
    file = "../novel_variants_for_snp_nexus_annotation/SNP_nexus_annotation.txt",
    header = TRUE
  ) %>%
  # filter to variants that have rs identification numbers (rs is not none)
  filter(
    .data = .,
    dbSNP != "None"
  ) %>%
  # convert "-" to a "*" for indel ref/alt alleles to match the format that was output from ANNOVAR
  mutate(
    .data = .,
    REF.Allele = if_else(
      condition = REF.Allele == "-",
      true = "*",
      false = REF.Allele
    )
  ) %>%
  mutate(
    .data = .,
    ALT.Allele..IUPAC. = if_else(
      condition = ALT.Allele..IUPAC. == "-",
      true = "*",
      false = ALT.Allele..IUPAC.
    )
  ) %>%
  # create a SNP ID column that matches that variant ID format from the ANNOVAR annotation data
  mutate(
    .data = .,
    VariantID_fromInputFile = paste(Chromosome,Position,REF.Allele,ALT.Allele..IUPAC.,sep = ":")
  ) %>%
  # select the snp ID and rs ID columns
  select(
    .data = .,
    VariantID_fromInputFile,
    dbSNP
  ) %>%
  # make sure rows are unique
  distinct(.data = .)

# load the lookup table of rs-numbers and star alleles
StarAlleleTable <-
  read.delim(
    file = "../starAllele_rsNumber_lookupTable/starAllele_rsNumber_map.txt",
    header = TRUE
  )

# load the table of clinical pharmGKB annotations with clinical levels of evidence for a variant influencing drug metabolism or drug response
pharmGKBevidenceData <-
  read.delim(
    file = "./pharmGKB_annotation_file/clinical_annotations.tsv",
    header = TRUE
      ) %>%
  # separate multiple entries in the Variant.Haplotypes column into multiple rows based on a comma
  separate_rows(
    data = .,
    Variant.Haplotypes,
    sep = ", "
    ) %>%
  # select the Variant.Haplotypes and Level.of.Evidence columns only
  select(
    .data = .,
    Variant.Haplotypes,
    Level.of.Evidence
    ) %>%
  # replace UGT1A1*36 and UGT1A1*37 star alleles with their rs-number pasted together with the star allele since they have the same rs-number
  mutate(
    .data = .,
    Variant.Haplotypes = if_else(
                                 condition = (Variant.Haplotypes == "UGT1A1*36") | (Variant.Haplotypes == "UGT1A1*37"),
                                 true = paste0(Variant.Haplotypes,":","rs3064744"),
                                 false = Variant.Haplotypes
                                 )
  ) %>%
  # replace any star alleles in the Variant.Haplotypes with an rs-number if the star allele is present in the StarAlleleTable
  mutate(
    .data = .,
    Variant.Haplotypes = Variant.Haplotypes %>%
                         lapply(
                           X = .,
                           FUN = function(currentVariant)
                           {
                               # if the current variant is not an rs-number, meaning it is a star allele,
                               # replace the star allele with an rs-number, if the star allele is present in the StarAlleleTable,
                               # if the current variant is an rs-number, leave it as is
                               if(grepl(pattern = "rs",x = currentVariant)==FALSE)
                               {
                                   # if the current variant is in the star allele table, replace the star allele with the rs-number,
                                   # otherwise, leave the star allele as is
                                   if(any(StarAlleleTable$starAllele == currentVariant)==TRUE)
                                   {
                                     VariantToReturn <-
                                     StarAlleleTable %>%
                                       filter(
                                         .data = .,
                                         starAllele == currentVariant
                                         ) %>%
                                       pull(
                                         .data = .,
                                         rsNumber
                                         )
                                   } else {
                                     VariantToReturn <- currentVariant
                                   }

                               } else {
                                 VariantToReturn <- currentVariant
                               }

                             return(VariantToReturn)
                           }
                           ) %>%
                          unlist(x = .)
    ) %>%
    # group by the Variant.Haplotypes and split into a list of dataframes based on Variant.Haplotypes
    group_by(
      .data = .,
      Variant.Haplotypes
    ) %>%
    group_split(.tbl = .) %>%
    # collapse the Variant.Haplotypes and levels of evidence into a single row
    # after sorting the level(s) of evidence in numeric order
    lapply(
      X = .,
      FUN = function(currentDataFrame)
      {
        TableToReturn <-
          tibble(
            "Variant.Haplotypes" = unique(x = currentDataFrame$Variant.Haplotypes),
            "pharmGKB_EvidenceLevel" = currentDataFrame$Level.of.Evidence %>%
              unique(x = .) %>%
              sort(x = .) %>%
              paste(.,collapse = " | ")
          )
        
        return(TableToReturn)
      }
    ) %>%
    # bind the list of dataframes back together by row
    do.call(
      what = "rbind",
      args = .
    ) %>% 
  # remove the rs-number from the UGT1A1*36 and *37 alleles
  mutate(
    .data = .,
    Variant.Haplotypes = if_else(
                                condition = (Variant.Haplotypes == "UGT1A1*36:rs3064744") | (Variant.Haplotypes == "UGT1A1*37:rs3064744"),
                                true = str_replace(
                                                    string = Variant.Haplotypes,
                                                    pattern = ":rs3064744",
                                                    replacement = ""
                                                    ),
                                false = Variant.Haplotypes
                                )
  )

# for allVariants, singletons, and indels join the ANNOVAR and pharmGKB files with the VEP annotations
# then join the data with the star allele rs-number table to obtain star alleles
# and also update any novel variants that were found to have rs-numbers with SNP-nexus
# save the resulting file to the file system
mapply(
  FUN = function(currentANNOVARfile,currentVEPfile,currentFileName)
    {
        # create a directory for storing data if it doesn't exist
        if(!dir.exists("../annotation_output/"))
        {
          dir.create(path = "../annotation_output/")
        }
    
      FileToSave <-
        currentANNOVARfile %>%
          # join the annovar and vep files, keeping all of the variant IDs in the VEP file
          right_join(
            x = .,
            y = currentVEPfile,
            by = c("VariantID_fromInputFile" = "VariantID_VEP")
            ) %>%
          # join the resulting file with the pharmGKB evidence table based on the rs-numbers present in the VEP data
          left_join(
            x = .,
            y = pharmGKBevidenceData,
            by = c("rsNumber_VEP" = "Variant.Haplotypes")
            ) %>%
          # join the resulting file with the starAllele table based on the variant ID in the VEP data
          left_join(
            x = .,
            y = StarAlleleTable %>%
                select(
                  .data = .,
                  VariantID_fromInputFile,
                  starAllele
                  ),
            by = "VariantID_fromInputFile"
            ) %>%
          # if the star allele is UGT1A1*36 or UGT1A1*37, update the pharmGKB level of evidence
          # to match the level of evidence in the pharmGKBevidenceData
          mutate(
            .data = .,
            pharmGKB_EvidenceLevel = mapply(
                                            FUN = function(currentStarAllele,currentpharmGKB_EvidenceLevel)
                                            {
                                                    # if the current star allele is NA, return the currentpharmGKB_EvidenceLevel
                                                    if(is.na(x = currentStarAllele)==TRUE)
                                                    {
                                                      EvidenceLevelToReturn <- currentpharmGKB_EvidenceLevel
                                                      
                                                      # if the current star allele is UGT1A1*36 or UGT1A1*37, return the evidence level from the pharmGKBevidenceData
                                                    } else if ((currentStarAllele == "UGT1A1*36") | (currentStarAllele == "UGT1A1*37")) {
                                                      
                                                      # filter the pharmGKBevidenceData to the currentStarAllele
                                                      # and select the pharmGKB evidence level
                                                      EvidenceLevelToReturn <-
                                                        pharmGKBevidenceData %>%
                                                        select(
                                                          .data = .,
                                                          Variant.Haplotypes,
                                                          "pharmGKB_EvidenceLevel*" = pharmGKB_EvidenceLevel
                                                          ) %>%
                                                        filter(
                                                          .data = .,
                                                          Variant.Haplotypes == currentStarAllele
                                                          ) %>%
                                                          pull(
                                                            .data = .,
                                                            `pharmGKB_EvidenceLevel*`
                                                            )
                                                      
                                                    } else {
                                                      # otherwise, leave the pharmGKB evidence level as is
                                                      EvidenceLevelToReturn <- currentpharmGKB_EvidenceLevel
                                                    }
                                              
                                                    return(EvidenceLevelToReturn)
                                            },
                                            starAllele,
                                            pharmGKB_EvidenceLevel,
                                            SIMPLIFY = FALSE
                                            ) %>%
                                            unlist(x = .)
            )
      
      # If the current file name is allVariants, identify the variant IDs that were removed when joining
      # the ANNOVAR and VEP annotation data above. These are variants that do not have RefSeq transcripts.
      # Save a file of these variants to the file system for filtering these variants out of the genetic association analysis
      if(currentFileName=="allVariants")
      {
        currentANNOVARfile %>%
          # anti join the annovar and vep files, to identify variants that do not have refSeq transcripts
          anti_join(
            x = .,
            y = currentVEPfile,
            by = c("VariantID_fromInputFile" = "VariantID_VEP")
            ) %>%
          select(
            .data = .,
            VariantID_fromInputFile
            ) %>%
          # save the variant IDs to the file system with out a header
          write_tsv(
            x = .,
            file = "../annotation_output/variants_without_RefSeq_annotations.txt",
            col_names = FALSE
            )
      }
      
      # add any additional rsIDs from the SNPnexus annotation to the file before it is saved to the file system
        # update the ANNOVARandPharmGKBdata with the rsID(s) identified using SNP-nexus if there are any [i.e., nrow(SNPnexusRSids > 0)]
         if(nrow(x = SNPnexusRSids) > 0)
         {
            FileToSave <-
              FileToSave %>%
              mutate(
                .data = .,
                snp138 = if_else(
                  condition = VariantID_fromInputFile == SNPnexusRSids$VariantID_fromInputFile,
                  true = SNPnexusRSids$dbSNP,
                  false = snp138
                ),
                avsnp138 = if_else(
                  condition = VariantID_fromInputFile == SNPnexusRSids$VariantID_fromInputFile,
                  true = SNPnexusRSids$dbSNP,
                  false = avsnp138
                ),
                avsnp142 = if_else(
                  condition = VariantID_fromInputFile == SNPnexusRSids$VariantID_fromInputFile,
                  true = SNPnexusRSids$dbSNP,
                  false = avsnp142
                )
              )
         }
      
      FileToSave %>%
          write_tsv(
            x = .,
            file = paste0("../annotation_output/",currentFileName,"_","annovar_pharmGKB_and_VEP_data.tsv")
              )
          
    },
  ANNOVARdata,
  VEPannotationData,
  names(x = ANNOVARdata),
  SIMPLIFY = FALSE
  )