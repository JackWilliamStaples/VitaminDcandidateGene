####### this is a wrapper script for running the entire vitamin D genetic analysis pipeline
setwd(dir = "~/Vitamin_D_project_MCCF/Vitamin_D_project_MCCF/")

# print current time at start to the screen
system(command = 'date +"%r"')

### grant permission to execute all the scripts in the pipeline
system(command = "chmod +x ./scripts/genetic_characterization/Vitamin_D_quality_control.sh")
system(command = "chmod +x ./scripts/genetic_characterization/Vitamin_D_vcf_summary_stats.sh")
system(command = "chmod +x ./scripts/genetic_characterization/ANNOVAR_annotation.sh")
system(command = "chmod +x ./scripts/genetic_characterization/pgx_pop_star_allele_caller.sh")
system(command = "chmod +x ./scripts/genetic_characterization/stargazer_vitaminD.sh")
system(command = "chmod +x ./scripts/genetic_characterization/pypgx_star_allele_caller.sh")
system(command = "chmod +x ./scripts/genotype_phenotype_association/candidate_gene_association_analysis.sh")

### execute each script

# perform quality control and phasing of the vitamin D sequencing data of batch 1 and batch 2 of sequencing
# and combine all variant calling data into a single, phased VCF file
system(command = "./scripts/genetic_characterization/Vitamin_D_quality_control.sh")

# compute VCFtools and BCFtools summary stats on the vitamin D sequencing data after quality control and phasing
# additionally, identify singleton variants and INDELS
system(command = "./scripts/genetic_characterization/Vitamin_D_vcf_summary_stats.sh")

# annotate the genetic data (all variants, indels, and singletons) that have been phased and undergone quality control with 
# ANNOVAR, pharmGKB, and the Ensemble Variant Effect Predictor (VEP)
# Additionally, use SNP-nexus to identify any variants that were not given existing variants IDs by ANNOVAR or VEP
system(command = "./scripts/genetic_characterization/ANNOVAR_annotation.sh")

# perform PGx-POP star allele calling on the VCF file of all variants for available genes (UGT1A1 only)
system(command = "./scripts/genetic_characterization/pgx_pop_star_allele_caller.sh")

# Call star alleles on the variants that have undergone quality control with starGazer for genes CYP2R1 and UGT1A4 (needs python version 3 to run)
# Additionally call star alleles for CYP3A4 variants with a custom method
system(command = "./scripts/genetic_characterization/stargazer_vitaminD.sh")

# call star alleles with pypgx for genes cyp2r1, cyp3a4, ugt1a1, and ugt1a4
system(command = "./scripts/genetic_characterization/pypgx_star_allele_caller.sh")

# consolidate the star allele haplotype, diplotype, and phenotype calling results from
# starGazer, PGxPOP, and custom star allele calling into tables for manuscripts
source(file = "./scripts/genetic_characterization/consolidate_starAllele_data.R")

# source the R script for creating variant summary statistics tables
source(file = "./scripts/genetic_characterization/compute_variant_annotation_summaries.R")

# perform candidate gene variant-phenotype association analysis
system(command = "./scripts/genotype_phenotype_association/candidate_gene_association_analysis.sh")

# print current time at the end to the screen
system(command = 'date +"%r"')



