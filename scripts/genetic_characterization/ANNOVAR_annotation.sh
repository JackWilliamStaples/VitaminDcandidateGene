# add the bioinformatics software directory to the path
PATH=$PATH:~/bioinformatics_software/

# if a working directory for annotating genetic data does not exist, create it
if [[ -d ~/Vitamin_D_project_MCCF/Vitamin_D_project_MCCF/annovar_annotation ]]
then
  echo "Variant annotation working directory exists"
else
  echo "Created variant annnotation working directory"
  mkdir ~/Vitamin_D_project_MCCF/Vitamin_D_project_MCCF/annovar_annotation/
fi

# working directory for this is:
cd ~/Vitamin_D_project_MCCF/Vitamin_D_project_MCCF/annovar_annotation/

##### code for downloading ANNOVAR data for annotation, only run this code if necessary
#	#Downloaded refGeneWithVer hg19 annotation information
#		perl ~/bioinformatics_software/annovar/annotate_variation.pl --buildver hg19 --hgvs --downdb --webfrom annovar refGeneWithVer ~/bioinformatics_software/annovar/humandb/
#	#Downloaded hg19 cytogenetic band identifiers
#		perl ~/bioinformatics_software/annovar/annotate_variation.pl --downdb --buildver hg19 cytoBand ~/bioinformatics_software/annovar/humandb/
#	#Downloaded current 1000 genomes annotation for whole genome data
#		perl ~/bioinformatics_software/annovar/annotate_variation.pl --buildver hg19 --downdb --webfrom annovar 1000g2015aug ~/bioinformatics_software/annovar/humandb/
#	#Downloaded current exac annotation for whole exome data
#		perl ~/bioinformatics_software/annovar/annotate_variation.pl --buildver hg19 --downdb --webfrom annovar exac03 ~/bioinformatics_software/annovar/humandb/
#	#Downloaded dbnsfp30a: includes SIFT, PolyPhen2 HDIV, PolyPhen2 HVAR, LRT, MutationTaster, MutationAssessor, FATHMM, MetaSVM, MetaLR, VEST, CADD, GERP++, DANN, fitCons, PhyloP and SiPhy scores, but ONLY on coding variants
#		perl ~/bioinformatics_software/annovar/annotate_variation.pl --buildver hg19 --downdb --webfrom annovar dbnsfp35a ~/bioinformatics_software/annovar/humandb/
#	#Downloaded dbscsnv11: dbscSNV version 1.1 for splice site prediction by AdaBoost and Random Forest, which score how likely that the variant may affect splicing
#		perl ~/bioinformatics_software/annovar/annotate_variation.pl --buildver hg19 --downdb --webfrom annovar dbscsnv11 ~/bioinformatics_software/annovar/humandb/
#	#Downloaded dbsnp 138 database and avsnp150 to see if there were any differences in the version/to see if one was better than the other for recognizing indels (the avsnp data bases are dbsnp but with left normalization)
#		perl ~/bioinformatics_software/annovar/annotate_variation.pl --buildver hg19 --downdb --webfrom annovar snp138 ~/bioinformatics_software/annovar/humandb/
#		perl ~/bioinformatics_software/annovar/annotate_variation.pl --buildver hg19 --downdb --webfrom annovar avsnp138 ~/bioinformatics_software/annovar/humandb/
#		perl ~/bioinformatics_software/annovar/annotate_variation.pl --buildver hg19 --downdb --webfrom annovar avsnp142 ~/bioinformatics_software/annovar/humandb/
#	# Downloaded clinVar registry
#		perl ~/bioinformatics_software/annovar/annotate_variation.pl --buildver hg19 --downdb --webfrom annovar clinvar_20190305 ~/bioinformatics_software/annovar/humandb/
#	# Downloaded catalog of variants previously reported in GWAS studies
#		perl ~/bioinformatics_software/annovar/annotate_variation.pl -build hg19 -downdb gwasCatalog ~/bioinformatics_software/annovar/humandb/
#	# Download transcription factor binding site motif annotation
#		perl ~/bioinformatics_software/annovar/annotate_variation.pl -buildver hg19 -downdb tfbsConsSites ~/bioinformatics_software/annovar/humandb/
#	#Downloaded GWAVA scores (tentative) [too many GB]
#		#perl ~/bioinformatics_software/annovar/annotate_variation.pl --buildver hg19 --downdb --webfrom annovar gwava ~/bioinformatics_software/annovar/humandb/
#	# Downloaded CADD scores (tentative) [350GB]
#		#perl ~/bioinformatics_software/annovar/annotate_variation.pl --buildver hg19 --downdb --webfrom annovar cadd13 ~/bioinformatics_software/annovar/humandb/
#	# Downloaded DANN scores (tentative) [200GB]
#		#perl ~/bioinformatics_software/annovar/annotate_variation.pl --buildver hg19 --downdb --webfrom annovar dann ~/bioinformatics_software/annovar/humandb/
#	# Downloaded COSMIC database
#		perl ~/bioinformatics_software/annovar/annotate_variation.pl --buildver hg19 --downdb --webfrom annovar cosmic70 ~/bioinformatics_software/annovar/humandb/
	
# Annotate the clean VCF file for all samples with INFO column recomputed
bcftools +fill-tags "../quality_control_final_output_directory/vitamin_D_merged_clean.gt.vcf.gz" -o vitaminD-merge_clean_updateINFO_allVariants.vcf -- -t all
perl ~/bioinformatics_software/annovar/table_annovar.pl "vitaminD-merge_clean_updateINFO_allVariants.vcf" ~/bioinformatics_software/annovar/humandb/ --outfile "merge_annotation" --buildver hg19 --onetranscript --nastring '.' --protocol refGeneWithVer,cytoBand,gwasCatalog,tfbsConsSites,1000g2015aug_all,1000g2015aug_afr,1000g2015aug_amr,1000g2015aug_eas,1000g2015aug_eur,1000g2015aug_sas,exac03,dbnsfp35a,dbscsnv11,snp138,avsnp138,avsnp142,clinvar_20190305,cosmic70 --operation g,r,r,r,f,f,f,f,f,f,f,f,f,f,f,f,f,f --remove --vcfinput 

# create an annotation output directory if it does not exist and
# move all of the final annotated files to the annotation_output directory
if [[ -d ../annotation_output ]]
then
  echo "Annotation output directory exists"
else
  echo "Created annotation output directory"
  mkdir ../annotation_output/
fi
mv *multianno* ../annotation_output/
  
# make a directory for storing a pharmGKB annotation file from the pharmGKB website (https://www.pharmgkb.org/downloads) if it does not exist
if [[ -d ./pharmGKB_annotation_file ]]
then
  echo "Directory for pharmGKB file exists"
else
  echo "Created directory for pharmGKB annotation file"
  mkdir ./pharmGKB_annotation_file/
fi

# change directories to the pharmGKB_annotation_file directory
cd ./pharmGKB_annotation_file/
# download the pharmGKB clinicalAnnotations.zip  file to this directory
wget https://api.pharmgkb.org/v1/download/file/data/clinicalAnnotations.zip 
# unzip the clinicalAnnotations.zip file
unzip clinicalAnnotations.zip
# change directories back to the annovar_annotation directory
cd ../

# Annotate the VCF file that was filtered down to singletons only
perl ~/bioinformatics_software/annovar/table_annovar.pl "../vcf_summary_stats/vitaminD-merge_clean_updateINFO_singletons.vcf.recode.vcf" ~/bioinformatics_software/annovar/humandb/ --outfile "singleton_annotation" --buildver hg19 --onetranscript --nastring '.' --protocol refGeneWithVer,cytoBand,gwasCatalog,tfbsConsSites,1000g2015aug_all,1000g2015aug_afr,1000g2015aug_amr,1000g2015aug_eas,1000g2015aug_eur,1000g2015aug_sas,exac03,dbnsfp35a,dbscsnv11,snp138,avsnp138,avsnp142,clinvar_20190305,cosmic70 --operation g,r,r,r,f,f,f,f,f,f,f,f,f,f,f,f,f,f --remove --vcfinput 

# move the singleton annotation data to the newly made singleton_annotation_output directory
if [[ -d ./singleton_annotation_output ]]
then
  echo "Directory for singleton annotation output exists"
else
  echo "Created directory for singleton annotation output"
  mkdir ../singleton_annotation_output/
fi  
mv *singleton* ../singleton_annotation_output/

# Annotate the VCF file that was filtered down to INDELS only
perl ~/bioinformatics_software/annovar/table_annovar.pl "../vcf_summary_stats/vitaminD-merge_clean_updateINFO.INDELS.vcf" ~/bioinformatics_software/annovar/humandb/ --outfile "merge_indel_annotation" --buildver hg19 --onetranscript --nastring '.' --protocol refGeneWithVer,cytoBand,gwasCatalog,tfbsConsSites,1000g2015aug_all,1000g2015aug_afr,1000g2015aug_amr,1000g2015aug_eas,1000g2015aug_eur,1000g2015aug_sas,exac03,dbnsfp35a,dbscsnv11,snp138,avsnp138,avsnp142,clinvar_20190305,cosmic70 --operation g,r,r,r,f,f,f,f,f,f,f,f,f,f,f,f,f,f --remove --vcfinput 

# move the INDEL annotation data to the newly made indel_annotation_directory
if [[ -d ./indel_annotation_output ]]
then
  echo "Directory for indel annotation output exists"
else
  echo "Created directory for indel annotation output"
  mkdir ../indel_annotation_output/
fi
mv *indel* ../indel_annotation_output/

# Annotate the annovar text files additionally with ADME optimized prediction framework scrores as 
# described in Zhou, et al., 2018, The Pharmacogenomics Journal, 19:115-126.
# the tidyverse package must be installed for this script to run
Rscript --no-save ../scripts/genetic_characterization/ADME_optimized_framework_prediction.R

# create a directory for storing variant effect predictor (VEP) annotation output if it doesn't exist
if [[ -d ../VEP_annotation_output ]]
then
  echo "VEP annotation directory exists"
else
  echo "Created VEP annotation directory"
  mkdir ../VEP_annotation_output/
fi
    
# manually perform an ANNOTATION of VCF files with all variants, singletons, and INDELS using the VEP web interface
# see VEP_settings.pdf screenshot in ../VEP_annotation_output for the settings used for the annotation
# the input files used were:
  #./vitaminD-merge_clean_updateINFO_allVariants.vcf
  #../vcf_summary_stats/vitaminD-merge_clean_updateINFO_singletons.vcf.recode.vcf
  #../vcf_summary_stats/vitaminD-merge_clean_updateINFO.INDELS.vcf
# the results were manually downloaded from the VEP website interface and saved to ../VEP_annotation_output with the file names:
  #allVariants_VEP_results.txt
  #singletons_VEP_results.txt
  #indels_VEP_results.txt
      
# Parse the VEP results with an R script for any interesting annotations that are not included as part of the ANNOVAR annotation 
# and join the useful columns from the VEP annotation with ANNOVAR and pharmGKB annotated data.
# Additionally, add pharmGKB clinical annotations and star alleles (based on rs-number matching) to the annotation file
Rscript --no-save ../scripts/genetic_characterization/combine_ANNOVAR_and_pharmGKBdata_with_VEP_data.R 


# create various plots of the annotation results 
Rscript --no-save ../scripts/genetic_characterization/plot_annotation_results.R

# create a VCF file of variants with a minor allele count (MAC) of at least 5 
# and remove variants that DO NOT have refSeq transcripts
vcftools --vcf ./vitaminD-merge_clean_updateINFO_allVariants.vcf --exclude ../annotation_output/variants_without_RefSeq_annotations.txt --mac 5 --recode --recode-INFO-all --out temp_vcf
# create a PLINK file with the VCF file from the previous step
plink --vcf temp_vcf.recode.vcf --keep-allele-order --allow-no-sex --make-bed --out variantsWithMACofAtLeast5
rm temp_vcf*

# create a VCF file of variants with a minor allele count (MAC) of at most 4 and a MAC of at least 1
# and remove variants that DO NOT have refSeq transcripts
vcftools --vcf ./vitaminD-merge_clean_updateINFO_allVariants.vcf --exclude ../annotation_output/variants_without_RefSeq_annotations.txt --mac 1 --max-mac 4 --recode --recode-INFO-all --out temp_vcf
# create a PLINK file with the VCF file from the previous step
plink --vcf temp_vcf.recode.vcf --keep-allele-order --allow-no-sex --make-bed --out variantsWithMACofAtMost4
rm temp_vcf*

# create a VCF file of variants with a minor allele count (MAC) equal to 0 that have refSeq transcripts
vcftools --vcf ./vitaminD-merge_clean_updateINFO_allVariants.vcf --exclude ../annotation_output/variants_without_RefSeq_annotations.txt --mac 0 --max-mac 0 --recode --recode-INFO-all --out temp_vcf
# create a PLINK file with the VCF file from the previous step
plink --vcf temp_vcf.recode.vcf --keep-allele-order --allow-no-sex --make-bed --out variantsWithMACofZero
rm temp_vcf*

# annotate the VCF files of vitamin D phenotype outliers
for currentSampleID in {"302172_302172","302061_302061","302066_302066","336170_336170","302297_302297","336192_336192","302269_302269","302063_302063","302322_302322"}
do
  perl ~/bioinformatics_software/annovar/table_annovar.pl "../vcf_summary_stats/phenotype_outlier_${currentSampleID}.vcf" ~/bioinformatics_software/annovar/humandb/ --outfile "annotated_variants_${currentSampleID}" --buildver hg19 --onetranscript --nastring '.' --protocol refGeneWithVer,cytoBand,gwasCatalog,tfbsConsSites,1000g2015aug_all,1000g2015aug_afr,1000g2015aug_amr,1000g2015aug_eas,1000g2015aug_eur,1000g2015aug_sas,exac03,dbnsfp35a,dbscsnv11,snp138,avsnp138,avsnp142,clinvar_20190305,cosmic70 --operation g,r,r,r,f,f,f,f,f,f,f,f,f,f,f,f,f,f --remove --vcfinput
done

# move the annotated variants of vitamin D phenotype outliers to a new folder
mkdir ./phenotype_outliers/
mv *hg19_multianno.txt ./phenotype_outliers/

# filter the variants for the phenotype outliers to clinical and deleterious variants
Rscript --no-save ../scripts/genetic_characterization/identify_clinical_and_deleterious_variants_for_phenotype_outliers.R
          