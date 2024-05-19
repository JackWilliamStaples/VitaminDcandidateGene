# this is a script for computing summary stats of vitamin D variant calling data
PATH=$PATH:~/bioinformatics_software/

# create a working directory for computing vcf summary statistics if it does not exist
if [[ -d ~/Vitamin_D_project_MCCF/Vitamin_D_project_MCCF/vcf_summary_stats ]]
then
  echo "Vcf summary statistics directory exists"
else
  echo "Created working directory for computing vcf summary statistics"
  mkdir ~/Vitamin_D_project_MCCF/Vitamin_D_project_MCCF/vcf_summary_stats/
fi

# switch to the working directory
cd ~/Vitamin_D_project_MCCF/Vitamin_D_project_MCCF/vcf_summary_stats/

####### generated bcftools stats file on the clean, phased VCF file:
bcftools stats "../quality_control_final_output_directory/vitamin_D_merged_clean.gt.vcf.gz" > "vitaminD_bcftools_stat.vchk"
		
##### filter the VCF down to singletons only by filtering sites with a minor allele count (MAC) of exactly 1
bcftools +fill-tags "../quality_control_final_output_directory/vitamin_D_merged_clean.gt.vcf.gz" -o temp_allele_frequency_update.vcf -- -t all
vcftools --vcf temp_allele_frequency_update.vcf --mac 1 --max-mac 1 --recode --recode-INFO-all --out vitaminD-merge_clean_updateINFO_singletons.vcf
rm -r temp_allele_frequency_update.vcf

##### filter VCF files down to INDELS only using BCFtools
bcftools +fill-tags "../quality_control_final_output_directory/vitamin_D_merged_clean.gt.vcf.gz" -o temp_allele_frequency_update.vcf -- -t all
bcftools view "temp_allele_frequency_update.vcf" --types indels -o "vitaminD-merge_clean_updateINFO.INDELS.vcf"
rm -r temp_allele_frequency_update.vcf

##### filter VCF files to individuals with vitamin D phenotype outliers with bcftools
# the phenotype outliers are:
  # (1) visit code A-150 with VCF file ID 302172 with 25(OH)D3 concentration of 170 ng/mL
  # (2) visit code A-020 with VCF file ID 302061 with vitamin D2 concentration of 41982.5 pg/mL
  # (3) visit code A-025 with VCF file ID 302066 with vitamin D2 concentration of 40790 pg/mL
  # (4) visit code A-333 with VCF file ID 336170 with 25(OH)D2 concentration of 298453 pg/mL
  # (5) visit code A-297 with VCF file ID 302297 with 1,25(OH)2D3 concentration of 129 pg/mL
  # (6) visit code A-363 with VCF file ID 336192 with 1,25(OH)2D3 concentration of 117 pg/mL
  # (7) visit code A-267 with VCF file ID 302269 with 1,25(OH)2D3 concentration of 109 pg/mL
  # (8) visit code A-022 with VCF file ID 302063 with 4-beta,25(OH)2D3 concentration of 421.4 pg/mL
  # (9) visit code A-323 with VCF file ID 302322 with 4-beta,25(OH)2D3 concentration of 383 pg/mL
bcftools +fill-tags "../quality_control_final_output_directory/vitamin_D_merged_clean.gt.vcf.gz" -o temp_allele_frequency_update.vcf -- -t all
for sampleID in {"302172_302172","302061_302061","302066_302066","336170_336170","302297_302297","336192_336192","302269_302269","302063_302063","302322_302322"}
do
  bcftools view -s "${sampleID}" "temp_allele_frequency_update.vcf" -o "phenotype_outlier_${sampleID}.vcf"
done
rm -r temp_allele_frequency_update.vcf
  
  