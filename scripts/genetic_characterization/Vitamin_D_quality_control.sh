# make sure the path points to the bioinformatics software directory
PATH=$PATH:~/bioinformatics_software/

# create a working directory for the quality control if it does not exist
if [[ -d ~/Vitamin_D_project_MCCF/Vitamin_D_project_MCCF/quality_control_working_directory ]]
then
  echo "Quality control working directory exists"
else
  echo "Created working directory for quality control"
  mkdir ~/Vitamin_D_project_MCCF/Vitamin_D_project_MCCF/quality_control_working_directory/
fi

# working directory for this is ~/Vitamin_D_project_MCCF/Vitamin_D_project_MCCF/quality_control_working_directory/
cd ~/Vitamin_D_project_MCCF/Vitamin_D_project_MCCF/quality_control_working_directory/

#################################################### quality control of VCF file from the first batch of Univ. of Washington sequencing #########################################################################

# filter the variants in the VCF from the first batch of sequencing down to variants that PASS the University of Washington quality control criteria using bcftools
bcftools view -i 'FILTER="PASS"' ../VCF_batch1_sequencing/woodahl_grc_1.HF.final.vcf > VCF_batch1_PASS_only.vcf

# zip the VCF file of passing variants from batch 1 with bgzip and index with tabix
bgzip -c VCF_batch1_PASS_only.vcf > VCF_batch1_PASS_only.vcf.gz
tabix -p vcf VCF_batch1_PASS_only.vcf.gz

#split the multiallelic sites into seperate rows in the vcf file from batch 1 using bcftools
bcftools norm -m-any VCF_batch1_PASS_only.vcf.gz -o VCF_batch1_PASS_only_multiallelic_split.vcf
# Lines   total/split/realigned/skipped:	6035/213/0/0

### make a directory and download the  hs37d5.fa fasta reference genome file if it does not exist already

# test to see if the reference genome file exists
if [[ -f ../hs37d5_fasta_reference_genome/hs37d5.fa.gz ]]
then
  echo "The hs37d5_fasta_reference_genome already exists, and does not need to be downloaded."
else
  echo "The hs37d5_fasta_reference_genome file does not exist."
  echo "Starting download..................."
  # make a directory for storing the hs37d5.fa fasta reference genome file
  mkdir ../hs37d5_fasta_reference_genome/
  # change directories to the newly created directory
  cd ../hs37d5_fasta_reference_genome/
  # download the hs37d5.fa fasta reference genome file from 1000 genomes
  wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz
  echo " Reference genome has been downloaded. "
  # change directories back to the quality control working directory
  cd ../quality_control_working_directory/
fi
  
# normalize the indels using the hs37d5.fa reference genome so each indel in the indel split is represented uniquely in the file (see http://apol1.blogspot.com/2014/11/best-practice-for-converting-vcf-files.html for explanation)
bcftools norm VCF_batch1_PASS_only_multiallelic_split.vcf -f ../hs37d5_fasta_reference_genome/hs37d5.fa.gz -o VCF_batch1_PASS_only_multiallelic_split_indel_fix.vcf
# total/split/realigned/skipped:	6538/0/368/0

# convert all the existing variant IDs to unique variant IDs with the format chromosome number, position, reference allele, alternate allele. the variant IDs MUST be unique before loading into PLINK
bcftools annotate VCF_batch1_PASS_only_multiallelic_split_indel_fix.vcf -x ID -I +'%CHROM:%POS:%REF:%ALT' -o VCF_batch1_ready_for_plink.vcf
	
# load the VCF file of passing variants split at multiallelic sites into PLINK and name the plink files as "vitaminD_batch1_plink_1" 
plink --vcf VCF_batch1_ready_for_plink.vcf --keep-allele-order --make-bed --out vitaminD_batch1_plink_1
# 6538 variants and 286 people loaded into PLINK
 
# check to see that all site IDs are unique (this step should print nothing to the console)
cut -f2 vitaminD_batch1_plink_1.bim | sort | uniq -c | awk '$1>1'
	
# calculate missingness per individual sample and per SNP in batch 1 using plink version 1.9
plink --bfile vitaminD_batch1_plink_1 --keep-allele-order --missing

# delete individuals with missingness > 0.1 based on Anderson et. al. Nature paper on candidate gene study QC
plink --bfile vitaminD_batch1_plink_1 --keep-allele-order --mind 0.1 --make-bed --out vitaminD_batch1_plink_2
# 0 people removed due to missing genotype data (--mind)

# Delete SNPs with missingness >0.05 based on Anderson et. al. Nature paper on candidate gene study QC
plink --bfile vitaminD_batch1_plink_2 --keep-allele-order --geno 0.05 --make-bed --out vitaminD_batch1_plink_3
# 218 variants removed due to missing genotype data (--geno)
# 6320 variants and 286 people pass filters and QC
 
# obtain a new individual missingness report (plink.imiss) after removing the 218 variants in the previous step
plink --bfile vitaminD_batch1_plink_3 --keep-allele-order --missing
	
	  # lookup the individual missingness in the plink.imiss file for the sample that was sequenced twice on accident (the duplicate IDs are 302313 with visit code A-313 and 336348 with visit code A-552)
    #FID      IID         MISS_PHENO  N_MISS  N_GENO   F_MISS
    #302313   302313          Y        6     6320 0.0009494
    
# remove sample 302313 with visit code of A-313 since A-313 has the greater proportion of missing genotypes
# sample 302313 must be removed
# create a text file with the sample ID to remove if it does not exist
if [[ -f ./sample_to_remove.txt ]]
then
  echo "sample_to_remove.txt file already exists"
else 
  echo "sample_to_remove.txt file does not exist"
  echo "302313 302313" > sample_to_remove.txt
  echo "sample_to_remove.txt file created"
fi
 
# remove sample A-313 because sample A-313 has the greater proportion of missing genotypes
plink --bfile vitaminD_batch1_plink_3 --keep-allele-order --remove sample_to_remove.txt --make-bed --out vitaminD_batch1_plink_4
# --remove: 285 people remaining.
	
# compute HWE p-values for all SNPS
plink --bfile vitaminD_batch1_plink_4 --keep-allele-order --hardy

#create a file for zooming in on strongly deviating HWE SNPs
awk '{ if ($9 <0.00001) print $0 }' plink.hwe > plinkzoomhwe.hwe

# remove snps with HWE p-value < 0.00001 based on Anderson et. al. Nature paper on candidate gene study QC
plink --bfile vitaminD_batch1_plink_4 --keep-allele-order --hwe include-nonctrl 0.00001 --make-bed --out vitaminD_batch1_plink_5
#--hwe: 88 variants removed due to Hardy-Weinberg exact test.
# 6232 variants and 285 people pass filters and QC.

# keeping alleles at all alternate allele frequencies for now

# convert the QC'd plink file back to a VCF file 
plink --bfile vitaminD_batch1_plink_5 --keep-allele-order --recode vcf --out vitaminD-batch1_clean_vcf
# 6232 variants and 285 people pass filters and QC.
	
##################################### quality control of VCF file from the second batch of Univ. of Washington sequencing ############################################################################################

# filter the variants in the VCF from the second batch of sequencing down to variants that PASS the University of Washington quality control criteria using bcftools
bcftools view -i 'FILTER="PASS"' ../VCF_batch2_sequencing/woodahl_grc_custom_2.HF.final.vcf > VCF_batch2_PASS_only.vcf 

# zip the VCF file of passing variants from batch 2 with bgzip and index with tabix
bgzip -c VCF_batch2_PASS_only.vcf > VCF_batch2_PASS_only.vcf.gz
tabix -p vcf VCF_batch2_PASS_only.vcf.gz

#split the multiallelic sites into seperate rows in the vcf file from batch 2 using bcftools
bcftools norm -m-any VCF_batch2_PASS_only.vcf.gz -o VCF_batch2_PASS_only_multiallelic_split.vcf
# Lines   total/split/realigned/skipped:	5746/201/0/0

# normalize the indels using the hs37d5.fa reference genome so each indel in the indel split is represented uniquely in the file (see http://apol1.blogspot.com/2014/11/best-practice-for-converting-vcf-files.html for explanation)
bcftools norm VCF_batch2_PASS_only_multiallelic_split.vcf -f ../hs37d5_fasta_reference_genome/hs37d5.fa.gz -o VCF_batch2_PASS_only_multiallelic_split_indel_fix.vcf
# total/split/realigned/skipped:	6227/0/343/0

# convert all the existing variant IDs to unique variant IDs with the format chromosome number, position, reference allele, alternate allele. the variant IDs MUST be unique before loading into PLINK
bcftools annotate VCF_batch2_PASS_only_multiallelic_split_indel_fix.vcf -x ID -I +'%CHROM:%POS:%REF:%ALT' -o VCF_batch2_ready_for_plink.vcf

# load the VCF file of passing variants split at multiallelic sites into PLINK and name the plink files as "vitaminD_batch2_plink_1" 
plink --vcf VCF_batch2_ready_for_plink.vcf --keep-allele-order --make-bed --out vitaminD_batch2_plink_1
# 6227 variants and 184 people pass filters and QC.
	 
# check to see that all site IDs are unique (this step should print nothing to the console)
cut -f2 vitaminD_batch2_plink_1.bim | sort | uniq -c | awk '$1>1'
	
# calculate missingness per individual sample and per SNP in batch 2 using plink version 1.9
plink --bfile vitaminD_batch2_plink_1 --keep-allele-order --missing

# delete individuals with missingness > 0.1 based on Anderson et. al. Nature paper on candidate gene study QC
plink --bfile vitaminD_batch2_plink_1 --keep-allele-order --mind 0.1 --make-bed --out vitaminD_batch2_plink_2
# 0 people removed due to missing genotype data (--mind)
 
# Delete SNPs with missingness >0.05 based on Anderson et. al. Nature paper on candidate gene study QC
plink --bfile vitaminD_batch2_plink_2 --keep-allele-order --geno 0.05 --make-bed --out vitaminD_batch2_plink_3
# 170 variants removed due to missing genotype data (--geno)
# 6057 variants and 184 people pass filters and QC.

# compute individual missingness again after removing 170 variants
# this will create a file called plink.imiss
plink --bfile vitaminD_batch2_plink_3 --keep-allele-order --missing

# lookup the individual missingness in the plink.imiss file for the sample that was sequenced twice on accident (the duplicate IDs are 302313 with visit code A-313  and 336348 with visit code A-552)
#FID      IID         MISS_PHENO  N_MISS  N_GENO   F_MISS
#336348   336348          Y        2     6057 0.0003302
# sample A-313 has more missing data and was removed from the batch 1 sequencing data in an earlier step

# compute HWE p-values for all SNPS
plink --bfile vitaminD_batch2_plink_3 --keep-allele-order --hardy

#create a file for zooming in on strongly deviating HWE SNPs
awk '{ if ($9 <0.00001) print $0 }' plink.hwe > plinkzoomhwe.hwe

# remove snps with HWE p-value < 0.00001 based on Anderson et. al. Nature paper on candidate gene study QC
plink --bfile vitaminD_batch2_plink_3 --keep-allele-order --hwe include-nonctrl 0.00001 --make-bed --out vitaminD_batch2_plink_4
#--hwe: 104 variants removed due to Hardy-Weinberg exact test.
#5953 variants and 184 people pass filters and QC.

# keeping alleles at all alternate allele frequencies for now
	
# convert the QC'd plink file back to a VCF file
plink --bfile vitaminD_batch2_plink_4 --keep-allele-order --recode vcf --out vitaminD-batch2_clean_vcf
# 5953 variants and 184 people pass filters and QC.
    
################################ merging the files from the two separate batches of sequencing and performing quality control on the merged set ##################################################

# zip (with bgzip) and index (with tabix) the clean QC'd files from batches 1 and 2
for batch in {batch1,batch2}
do 
  bgzip -c "vitaminD-${batch}_clean_vcf.vcf" > "vitaminD-${batch}_clean_ready_for_merging_vcf.vcf.gz" && tabix -p vcf "vitaminD-${batch}_clean_ready_for_merging_vcf.vcf.gz"
done

#merge the zipped VCF files from sequencing batches 1 and 2 of vitamin D sequencing using bcftools
#use the "none" parameter to ensure no new multiallelic sites are created, they are just left split into individual snps
# use the -0 parameter to set the variants that are absent from one file (batch1 or batch2) and present in the other (batch1 or batch2) to the reference allele for the missing genotypes that are created upon merging
bcftools merge -m none -0 vitaminD-batch1_clean_ready_for_merging_vcf.vcf.gz vitaminD-batch2_clean_ready_for_merging_vcf.vcf.gz -o vitaminD-batch1_batch2_initial_merge_vcf.vcf
	
# zip the merged VCF file with bgzip and index with tabix
bgzip -c vitaminD-batch1_batch2_initial_merge_vcf.vcf > vitaminD-batch1_batch2_initial_merge_vcf.vcf.gz
tabix -p vcf vitaminD-batch1_batch2_initial_merge_vcf.vcf.gz

# split the rows containing multi-allelic sites if there are any present with bcftools
bcftools norm -m-any vitaminD-batch1_batch2_initial_merge_vcf.vcf.gz -o vitaminD-batch1_batch2_initial_merge_vcf.multiallelic_split.vcf
# total/split/realigned/skipped:	7532/0/0/0

# normalize the indels using the hs37d5.fa reference genome so each indel in the indel split is represented uniquely in the file (see http://apol1.blogspot.com/2014/11/best-practice-for-converting-vcf-files.html for explanation)
bcftools norm vitaminD-batch1_batch2_initial_merge_vcf.multiallelic_split.vcf -f ../hs37d5_fasta_reference_genome/hs37d5.fa.gz -o vitaminD-batch1_batch2_initial_merge_vcf_indel_realign.multiallelic_split.vcf
# Line total/split/realigned/skipped:	7532/0/0/0

# convert all the existing variant IDs to unique variant IDs with the format chromosome number, position, reference allele, alternate allele. the variant IDs MUST be unique before loading into PLINK
bcftools annotate vitaminD-batch1_batch2_initial_merge_vcf_indel_realign.multiallelic_split.vcf -x ID -I +'%CHROM:%POS:%REF:%ALT' -o vitaminD-batch1_batch2_initial_merge_updated_variant_IDs.vcf

# load the merged VCF file into PLINK and name the plink files as "vitamin_D_merged_plink_1"
plink --vcf vitaminD-batch1_batch2_initial_merge_updated_variant_IDs.vcf --keep-allele-order --make-bed --out vitamin_D_merged_plink_1
# 7532 variants and 469 people pass filters and QC.
 
# check to see that all site IDs are unique (this step should print nothing to the console)
cut -f2 vitamin_D_merged_plink_1.bim | sort | uniq -c | awk '$1>1'
	
# create a text file of samples IDs from the QC'd data from batch 1
awk '{print $1,$2}' vitaminD_batch1_plink_5.fam > sample_IDs_batch1.txt

# create a text file of sample IDs from the QC'd data from batch 2 
awk '{print $1,$2}' vitaminD_batch2_plink_4.fam > sample_IDs_batch2.txt

# use the file of sample IDs from batch 2 to code the samples from batch 2 as "cases=2" and the samples from batch 1 as "controls=1"
Rscript --no-save ../scripts/genetic_characterization/code_batch1_and_2_samples_as_controls_and_cases.R

# use the file of sample IDs from batch 2 to code the samples from batch 2 as "cases=2" and the samples from batch 1 as "controls=1"
plink --bfile vitamin_D_merged_plink_1 --keep-allele-order --allow-no-sex --pheno control_case_batch1_batch2.txt --make-bed --out vitamin_D_merged_plink_case_control_to_check_for_missingness

#test for significant differences in genotype calls between batches 1 and 2
plink --bfile vitamin_D_merged_plink_case_control_to_check_for_missingness --allow-no-sex --keep-allele-order --test-missing midp --out vitamin_D_merged_plink_case_control_to_check_for_missingness

#identify variants that have a P-value of < 0.05 based on fisher's exact test for no differences in genotype missingness for individual snps sequenced in batch 1 and batch 2
Rscript --no-save ../scripts/genetic_characterization/alleles_with_batch_genotyping_effects.R

#remove the variants with P-value < 0.05 from the analysis from the file that was created before coding individual as cases and controls, individuals do not need to be coded as cases and controls for further analysis
plink --bfile vitamin_D_merged_plink_1 --allow-no-sex --keep-allele-order --exclude SNPsWithSignificantDifferencesInGenotypeMissingness.txt --make-bed --out vitamin_D_merged_plink_2
# 7449 variants and 469 people pass filters and QC.
# --exclude: 7449 variants remaining.

# do a final check for deviations from hardy weinberg equilibrium once again
# compute HWE p-values for all SNPS
plink --bfile vitamin_D_merged_plink_2 --keep-allele-order --hardy

# remove snps with HWE p-value < 0.00001 based on Anderson et. al. Nature paper on candidate gene study QC
plink --bfile vitamin_D_merged_plink_2 --keep-allele-order --hwe include-nonctrl 0.00001 --make-bed --out vitamin_D_merged_plink_3
# --hwe: 56 variants removed due to Hardy-Weinberg exact test.
# 7393 variants and 469 people pass filters and QC.

# recode the final file to a VCF
plink --bfile vitamin_D_merged_plink_3 --keep-allele-order --recode vcf --out vitamin_D_merged_clean
# 7393 variants and 469 people pass filters and QC.

# zip with bgzip and index with tabix
bgzip -c vitamin_D_merged_clean.vcf > vitamin_D_merged_clean.vcf.gz
tabix -p vcf vitamin_D_merged_clean.vcf.gz
 
# make a temporary directory for storing all files that end in .bed, .bim, .fam, .vcf, .vcf.gz, and .vcf.gz.tbi, .log, .nosex
mkdir ../temp_dir_for_intermediate_files/
# move file types listed above to the temp_dir_for_intermediate_files
mv *.bed ../temp_dir_for_intermediate_files/
mv *.bim ../temp_dir_for_intermediate_files/
mv *.fam ../temp_dir_for_intermediate_files/
mv *.vcf ../temp_dir_for_intermediate_files/
mv *.vcf* ../temp_dir_for_intermediate_files/
mv *.log* ../temp_dir_for_intermediate_files/
mv *.nosex* ../temp_dir_for_intermediate_files/
	
######################################################################## phasing and imputation with beagle version 4.1 created by Brian Browning (see https://faculty.washington.edu/browning/beagle/b4_1.html#download) #######################################

# working directory for this is:
mkdir ../phasing_and_imputation_working_directory/
cd ../phasing_and_imputation_working_directory/

# perform phasing with beagle using the DEFAULT parameters on the clean, QC'd VCF file with all 469 participants included
# with a random number seed for reproducibility
java -jar ~/bioinformatics_software/beagle.27Jan18.7e1.jar gt=../temp_dir_for_intermediate_files/vitamin_D_merged_clean.vcf.gz seed=777 out=vitamin_D_merged_clean.gt

# move all of the phased and imputed VCF files to the working directory for quality control
mv *gt* ../quality_control_working_directory/
 
# move back to the quality_control_working_directory
cd ../quality_control_working_directory/

# delete the phasing_and_imputation_working_directory now that it is no longer needed
rm -r ../phasing_and_imputation_working_directory/

#index the phased VCF file from the merge of batches 1 and 2 with tabix, the file is already zipped upon output from BEAGLE
tabix -p vcf vitamin_D_merged_clean.gt.vcf.gz

# perform principal components analysis using an LD-pruned set of snps with the final clean data set to create the plot pca_plot.pdf to see if there is a batch effect, 
# indicated by two distinct clusters on the PCA plot
plink --vcf vitamin_D_merged_clean.gt.vcf.gz --keep-allele-order --indep-pairwise 200 1 0.5 --out vitamin_D_LD_pruned_SNPs
plink --vcf vitamin_D_merged_clean.gt.vcf.gz --keep-allele-order --extract vitamin_D_LD_pruned_SNPs.prune.in --make-bed --out vitamin_D_no_linkage
plink --bfile vitamin_D_no_linkage --keep-allele-order --pca 10 --out pca_results
Rscript --no-save ../scripts/genetic_characterization/plot_pca_results_grouped_by_batch.R

# make a directory for storing all clean, phased data if it does not exist
if [[ -d ../quality_control_final_output_directory ]]
then
  echo "final output directory exists"
else
  echo "final output directory created"
  mkdir ../quality_control_final_output_directory/
fi

### move the clean, phased data to the final output directory for the quality control
mv vitamin_D_merged_clean.gt.vcf.gz ../quality_control_final_output_directory/
mv vitamin_D_merged_clean.gt.vcf.gz.tbi ../quality_control_final_output_directory/
  
# move all files that end in .bed, .bim, .fam, .vcf, .vcf.gz, and .vcf.gz.tbi, .log, and .nosex
# move file types listed above to the temp_dir_for_intermediate_files
mv *.bed ../temp_dir_for_intermediate_files/
mv *.bim ../temp_dir_for_intermediate_files/
mv *.fam ../temp_dir_for_intermediate_files/
mv *.log* ../temp_dir_for_intermediate_files/
mv *.nosex* ../temp_dir_for_intermediate_files/
mv *.prune.in ../temp_dir_for_intermediate_files/
mv *.prune.out ../temp_dir_for_intermediate_files/
  
# delete the temp_dir_for_intermediate_files and all of its contents now that it is no longer needed
rm -r ../temp_dir_for_intermediate_files/
    


