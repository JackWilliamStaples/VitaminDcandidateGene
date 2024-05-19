#### this is a wrapper script for candidate gene association analysis of vitamin D endocrinology gene variants with vitamin D and vitamin D metabolite levels

# add the bioinformatics software to the path
PATH=$PATH:~/bioinformatics_software/
# add the locusZoom executable to the path
PATH=$PATH:~/bioinformatics_software/locuszoom/bin/

# create a working directory for the genotype-phenotype association analysis if it does not exist
if [[ -d ~/Vitamin_D_project_MCCF/Vitamin_D_project_MCCF/genotype_phenotype_assoc_working_directory/ ]]
then
  echo "Genetic association analysis working directory already exists"
else
  echo "Created working directory for genetic association analysis"
  mkdir ~/Vitamin_D_project_MCCF/Vitamin_D_project_MCCF/genotype_phenotype_assoc_working_directory/
fi

# working directory for this is 
cd ~/Vitamin_D_project_MCCF/Vitamin_D_project_MCCF/genotype_phenotype_assoc_working_directory/

### create continuous and binary phenotype files for common and rare variant analysis as well as demographic/seasonal covariate files:
  # (1) this will create a report of missing values in the vitamin D and metabolite quantitation data
  # (2) this will identify study participants that have metabolite quantitation duplicates
  # (3) this will create a file called VitaminDplinkPhenotypeFile.txt with continuous phenotypes for linear regressions
  # (4) this will create a file with demographic and seasonal covariates for adjusting regressions
  # (5) this will create a file called VitaminDplinkPhenotypeFileBinary.txt with binary phenotypes for logistic regressions
  # (6) this will also create a PLINK phenotype file for each individual vitamin D phenotype (continuous and binary) to be used for the SKAT-O rare variant analysis
  # (7) this will create histograms and box and whisker plots for each metabolite
  # (8) this will also perform Shapiro Wilke's tests for deviations from normality on untransformed and log10() transformed metabolite values
  # (9) this will create correlation plots of all pairwise D3 and D2 metabolite correlations
  # (10) this will also test for heteroscedasticity for appropriate D3 and D2 metabolite correlations
# the tidyverse package must be installed for this script to work
Rscript --no-save ../scripts/genotype_phenotype_association/create_plink_phenotype_file_for_all_metabolites.R

# create a file containing the chromosome and nucleotide position range for each gene being analyzed in the study
# the tidyverse package must be installed for this to run 
# this will create multiple files with file names {gene name}_regionFile.txt
Rscript --no-save ../scripts/genotype_phenotype_association/create_gene_region_files.R

# make a directory for storing all of the intermediate regression association results
mkdir ./regression_results_temp/

# create individual plink files for each gene for variants with a minor allele count of at least 5 (MAC >= 5)
# and test the variants in each gene for association with each individual CONTINUOUS phenotype in the VitaminDplinkPhenotypeFile.txt plink phenotype file 
# with an additive allele dose association model while adjusting for multiple testing 
# and creating files with the mean and standard deviation of each phenotype grouped by genotype
for geneName in {CYP3A4,DHCR7,CYP2R1,CYP27B1,SULT2A1,UGT1A4,GC,VDR,CASR,CUBN,CYP24A1,RXRG,RXRB,RXRA}
do 
  vcftools --vcf ../annovar_annotation/vitaminD-merge_clean_updateINFO_allVariants.vcf --exclude ../annotation_output/variants_without_RefSeq_annotations.txt --mac 5 --recode --recode-INFO-all --out temp_vcf
  plink --vcf temp_vcf.recode.vcf --keep-allele-order --extract range "${geneName}_regionFile.txt" --pheno VitaminDplinkPhenotypeFile.txt --allow-no-sex --all-pheno --assoc qt-means --adjust --out "./regression_results_temp/${geneName}_additive_assoc"
  rm temp_vcf*
done

# create individual plink files for each gene for variants with a minor allele count of at least 5 (MAC >= 5)
# and test the variants in each gene for association with each individual BINARY phenotype in the VitaminDplinkPhenotypeFileBinary.txt plink phenotype file 
# with a logistic regression while adjusting for multiple testing
# the --ci flag is used to add 95% confidence intervals for odds ratios to the results file
# A random number seed must be set for this for reproducibility
# The --1 flag specifies case = 1, case = 0 phenotype
for geneName in {CYP3A4,DHCR7,CYP2R1,CYP27B1,SULT2A1,UGT1A4,GC,VDR,CASR,CUBN,CYP24A1,RXRG,RXRB,RXRA}
do 
  vcftools --vcf ../annovar_annotation/vitaminD-merge_clean_updateINFO_allVariants.vcf --exclude ../annotation_output/variants_without_RefSeq_annotations.txt --mac 5 --recode --recode-INFO-all --out temp_vcf
  plink --vcf temp_vcf.recode.vcf --keep-allele-order --extract range "${geneName}_regionFile.txt" --pheno VitaminDplinkPhenotypeFileBinary.txt --1 --allow-no-sex --all-pheno --logistic --adjust --ci 0.95 --seed 777 --out "./regression_results_temp/${geneName}_logistic_assoc"
  rm temp_vcf*
done

# Perform SNP clumping with PLINK for CONTINUOUS phenotypes (for more information on clumping, see: https://www.cog-genomics.org/plink/1.9/postproc#clump)
# to reduce the SNPs down to the most associated SNPs with phenotype (i.e., removing any SNPs that are strongly correlated with the causative SNP based on linkage disequilibrium).
# Do this with each gene and metabolite/metabolite ratio
for currentMetabolite in {active_1alpha25OH2D2_ngperml,active_1alpha25OH2D3_ngperml,beta_4beta25OH2D3_ngperml,Overall_primary_25OHD_ngperml,primary_25OHD2_ngperml,primary_25OHD3_ngperml,R_24R25OH2D3_ngperml,VitaminD2_ngperml,VitaminD3_ngperml,primary_25OHD3_VitaminD3_ratio,active_1alpha25OH2D3_primary_25OHD3_ratio,R_24R25OH2D3_primary_25OHD3_ratio,beta_4beta25OH2D3_primary_25OHD3_ratio,primary_25OHD2_VitaminD2_ratio,active_1alpha25OH2D2_primary_25OHD2_ratio}
do 
  for geneName in {CYP3A4,DHCR7,CYP2R1,CYP27B1,SULT2A1,UGT1A4,GC,VDR,CASR,CUBN,CYP24A1,RXRG,RXRB,RXRA} 
  do 
    vcftools --vcf ../annovar_annotation/vitaminD-merge_clean_updateINFO_allVariants.vcf --exclude ../annotation_output/variants_without_RefSeq_annotations.txt --mac 5 --recode --recode-INFO-all --out temp_vcf
    plink --vcf temp_vcf.recode.vcf --keep-allele-order --allow-no-sex --extract range "${geneName}_regionFile.txt" --clump-p1 1 --clump-r2 0.1 --clump-kb 250 --clump "./regression_results_temp/${geneName}_additive_assoc.${currentMetabolite}.qassoc" --clump-snp-field SNP --clump-field P --out "./regression_results_temp/clump_file${geneName}_additive_assoc.${currentMetabolite}"
    rm temp_vcf*
  done
done
    
# Perform SNP clumping with PLINK for BINARY phenotypes
# Do this with each gene and VitDsufficiency and VitDdeficiency phenotypes
for phenotype in {VitDsufficiency,VitDdeficiency}
do
  for geneName in {CYP3A4,DHCR7,CYP2R1,CYP27B1,SULT2A1,UGT1A4,GC,VDR,CASR,CUBN,CYP24A1,RXRG,RXRB,RXRA}
  do 
    vcftools --vcf ../annovar_annotation/vitaminD-merge_clean_updateINFO_allVariants.vcf --exclude ../annotation_output/variants_without_RefSeq_annotations.txt --mac 5 --recode --recode-INFO-all --out temp_vcf
    plink --vcf temp_vcf.recode.vcf --keep-allele-order --allow-no-sex --extract range "${geneName}_regionFile.txt" --clump-p1 1 --clump-r2 0.1 --clump-kb 250 --clump "./regression_results_temp/${geneName}_logistic_assoc.${phenotype}.assoc.logistic" --clump-snp-field SNP --clump-field P --out "./regression_results_temp/clump_file${geneName}_logistic_assoc.${phenotype}"
    rm temp_vcf*
  done
done

# make a directory for storing genotype matrices
mkdir ./genotype_matrices_temp/

# convert the VCF file to gentoype matrices with genotypes represented as (0,1,2) for each individual gene for variants with a minor allele count of at least 5 (MAC >= 5) only
for geneName in {CYP3A4,DHCR7,CYP2R1,CYP27B1,SULT2A1,UGT1A4,GC,VDR,CASR,CUBN,CYP24A1,RXRG,RXRB,RXRA}
do 
  vcftools --vcf ../annovar_annotation/vitaminD-merge_clean_updateINFO_allVariants.vcf --exclude ../annotation_output/variants_without_RefSeq_annotations.txt --mac 5 --recode --recode-INFO-all --out temp_vcf
  plink --vcf temp_vcf.recode.vcf --keep-allele-order --allow-no-sex --extract range "${geneName}_regionFile.txt" --recode A --out "./genotype_matrices_temp/${geneName}_geno_matrix"
  rm temp_vcf*
done

# use BCFtools to filter the VCF file of all variants to the vitamin D binding protein region (chr 4: 72607410-72671237)
# and create a sample x haplotype matrix for vitamin D binding protein (GC) only for a haplotype-trait association analysis with the variants rs7041 and rs4588
bcftools view "../quality_control_final_output_directory/vitamin_D_merged_clean.gt.vcf.gz" --regions 4:72607410-72671237 -o temp.vcf
bcftools convert --hapsample --vcf-ids temp.vcf -o "GC_haplotypes"
rm -r temp.vcf
gunzip GC_haplotypes.hap.gz

# make a file with variant IDs for the vitamin D haplotype defining SNPs: rs7041 (4:72618334:A:C) and rs4588 (4:72618323:G:T)
Rscript --no-save ../scripts/genotype_phenotype_association/create_VDBP_haplotype_snps_file.R

# calculate the linkage disequilibrium between rs7041 and rs4588 in the CSKT population
vcftools --gzvcf "../quality_control_final_output_directory/vitamin_D_merged_clean.gt.vcf.gz" --snps VDBP_haplotype_SNPs_forLDcalc.txt --recode --recode-INFO-all --out temp
plink --vcf temp.recode.vcf --r2 --ld-snp-list VDBP_haplotype_reference_SNP_forLDcalc.txt --ld-window 4000 --ld-window-kb 1000000 --ld-window-r2 0
rm -r temp.recode.vcf
# based on the resulting plink.ld file, the linkage disequilibrium is r2 = 0.291804

# make a temporary directory for deleting unnecessary files
mkdir ./temp_directory/

# move files that aren't needed to the temp_directory
mv *.log ./temp_directory/
mv *.nosex ./temp_directory/

# delete the temp_directory and all of its contents
rm -r ./temp_directory/

# combine all of the regression results for each gene and for each metabolite all into one file for variants that have a minor allele count of at least 5 (MAC >= 5)
# and combine this data with the variant annotation data for easy look-up of significant variants
Rscript --no-save ../scripts/genotype_phenotype_association/combine_regression_results_and_variant_annotation_data.R

# combine all genotype matrices for common variants tested with the additive allele association model with the metabolite data
# and also create a lookup table of rs identification numbers and variants IDs in the PLINK format (chrom:pos:ref:alt)
Rscript --no-save ../scripts/genotype_phenotype_association/aggregate_genotype_matrices_and_metabolite_data.R

# remove the temporary directory for storing raw genotype matrices now that it is no longer needed
rm -r ./genotype_matrices_temp/

# create a file of the linear regression results for the common vitamin D binding protein haplotype defining variants (rs7041 and rs4588)
Rscript --no-save ../scripts/genotype_phenotype_association/filter_regression_results_to_vitaminDbindingProteinHaplotypeVariants.R

# create a plot of UVB radiation data for each month in Polson to see if it matches the 6 week lag time that
# is expected between peak UVB exposure and peak 25(OH)D3 concentration
Rscript --no-save ../scripts/genotype_phenotype_association/plot_UVB_radiation_data.R

# prepare phenotype and covariate files for linear and logistic regressions in the next script
Rscript --no-save ../scripts/genotype_phenotype_association/prepare_files_for_linear_and_logistic_regressions.R

# This script performs analysis of continuous metabolite concentrations and ratios, the analysis includes:
  # (1) D3 and D2 metabolite phenotype principal components analysis to visualize metabolite quantitation batch effects
  # (2) Cosinor seasonal regression analysis to identify seasonal 25(OH)D3 trends with and without demographic covariates
  # (3) Grand summary statistics and summary statistics stratified by season and gender for demographic factors
  # (4) Summary statistics for clinical vitamin D status
  # (5) Univariate analysis of metabolite concentration and ratio versus individual seasonal and demographic factors
  # (6) Association analysis for significant relationships among seasonal and demographic factors (age, BMI, gender, and season)
  # (7) Tests for confounding relationships between demographic and seasonal factors
  # (8) Tests for confounding of genetic associations by seasonal and demographic factors
  # (9) multiple regressions of significant variants identified in primary analysis with adjustment for all demographic/seasonal covariates
  # (10) multiple regression with all significant causative variants identified in primary analysis to estimate polygenic contribution to phenotypic variability for each phenotype
  # (11) multiple regression with all significant causative variants identified in primary analysis and seasonal/demographic covariates together to estimate genetic and environmental contribution to variability of each phenotype
Rscript --no-save ../scripts/genotype_phenotype_association/linear_regressions_continuous_phenotypes.R

# use the genotype matrices with metabolite data added to create boxplots of the metabolite data grouped by genotype for significant variants only
Rscript --no-save ../scripts/genotype_phenotype_association/create_genotype_metabolite_boxplots.R

# This script performs analysis of binary vitamin D status phenotypes [sufficiency and deficiency], the analysis includes:
  # (1) Univariate and multiple logistic regressions of demographic/seasonal covariates with binary vitamin D status phenotypes
  # (2) multiple logistic regressions with individual genetic variants included with adjustment for all seasonal and demographic covariates
  # (3) Logistic regression of binary phenotypes versus all significant causative variants with and without seasonal/demographic covariates included
  # (4) Logistic regression of binary phenotypes versus common vitamin D binding protein diplotypes
Rscript --no-save ../scripts/genotype_phenotype_association/logistic_regressions_binary_phenotypes.R

# delete the temporary regression results directory
rm -r ./regression_results_temp/

# Create files for LocusZoom plots and lists of variants for calculating linkage disequilibrium.
# Create individual files with the reference SNP with the lowest P-value so linkage disequilibrium can be computed with reference to the reference SNP for each metabolite and gene.
# The calculated linkage disequilibrium is used to color the points in the locusZoom plots that are computed in the code below. 
Rscript --no-save ../scripts/genotype_phenotype_association/create_locusZoom_plot_files_and_snp_lists_for_LD_calculation.R

# create a temporary directory for storing locusZoom intermediate files
mkdir ./locusZoom_temp/

for currentPhenotype in {active_1alpha25OH2D2_ngperml,active_1alpha25OH2D3_ngperml,beta_4beta25OH2D3_ngperml,Overall_primary_25OHD_ngperml,primary_25OHD2_ngperml,primary_25OHD3_ngperml,R_24R25OH2D3_ngperml,VitaminD2_ngperml,VitaminD3_ngperml,primary_25OHD3_VitaminD3_ratio,active_1alpha25OH2D3_primary_25OHD3_ratio,R_24R25OH2D3_primary_25OHD3_ratio,beta_4beta25OH2D3_primary_25OHD3_ratio,primary_25OHD2_VitaminD2_ratio,active_1alpha25OH2D2_primary_25OHD2_ratio,VitDsufficiency,VitDdeficiency}
do
  for geneName in {CYP3A4,DHCR7,CYP2R1,CYP27B1,SULT2A1,UGT1A4,GC,VDR,CASR,CUBN,CYP24A1,RXRG,RXRB,RXRA}
  do
    # Filter the variant IDs in the VCF file that were used for the candidate gene association tests to variant IDs
    # that will be plotted with LocusZoom with the lowest P-value used for multi-allelic sites
    # using the SNPsForLDcalculation_${currentPhenotype}_${geneName}.txt files created in the previous script with VCFtools,
    vcftools --vcf "../annovar_annotation/vitaminD-merge_clean_updateINFO_allVariants.vcf" --snps "./Results_for_Calculating_LD_for_LocusZoom/SNPsForLDcalculation_${currentPhenotype}_${geneName}.txt" --recode --recode-INFO-all --out "./locusZoom_temp/${currentPhenotype}_${geneName}_VariantForCalculatingLocusZoom_linkageDisequilibrium"
    # Calculate linkage disequilibrium for the SNPs that will be included in LocusZoom plots using PLINK with reference to the SNP with the lowest P-value for that phenotype/gene,
    # a file with the reference SNP for the LD calculation was created in the R script run directly above 
    plink --vcf "./locusZoom_temp/${currentPhenotype}_${geneName}_VariantForCalculatingLocusZoom_linkageDisequilibrium.recode.vcf" --exclude ../annotation_output/variants_without_RefSeq_annotations.txt --r2 --ld-snp-list "ReferenceSNPs_for_LinkDis/referenceSNP_forLDcalc_${currentPhenotype}_${geneName}.txt" --ld-window 4000 --ld-window-kb 1000000 --ld-window-r2 0
    mv plink.ld "./locusZoom_temp/${currentPhenotype}_${geneName}_linkage_disequilibrium.txt"
    # Reformat the variant ID column in the VCF file using BCFtools to match the format required for locusZoom,
    bcftools annotate "./locusZoom_temp/${currentPhenotype}_${geneName}_VariantForCalculatingLocusZoom_linkageDisequilibrium.recode.vcf" -x ID -I +'%CHROM:%POS' -o "./locusZoom_temp/${currentPhenotype}_${geneName}_VariantForCalculatingLocusZoom_linkageDisequilibrium_updatedVariantID.vcf"
    # Zip and index the vcf files created with bgzip and tabix so they take up less memory
    bgzip -c "./locusZoom_temp/${currentPhenotype}_${geneName}_VariantForCalculatingLocusZoom_linkageDisequilibrium_updatedVariantID.vcf" > "./locusZoom_temp/${currentPhenotype}_${geneName}_VariantForCalculatingLocusZoom_linkageDisequilibrium_updatedVariantID.vcf.gz"
    tabix -p vcf "./locusZoom_temp/${currentPhenotype}_${geneName}_VariantForCalculatingLocusZoom_linkageDisequilibrium_updatedVariantID.vcf.gz"
  done
done

# Clean the files with linkage disequilibrium calculations that were output from PLINK in the step above into the format required by LocusZoom and combine the linkage disequilibrium data with the P-value data
Rscript --no-save ../scripts/genotype_phenotype_association/filter_LDresults_to_the_Variant_with_the_lowest_Pvalue.R

# create a file for storing locusZoom plots
mkdir ./locusZoom_plots/

# create the locusZoom plots for each phenotype and gene
for currentPhenotype in {active_1alpha25OH2D2_ngperml,active_1alpha25OH2D3_ngperml,beta_4beta25OH2D3_ngperml,Overall_primary_25OHD_ngperml,primary_25OHD2_ngperml,primary_25OHD3_ngperml,R_24R25OH2D3_ngperml,VitaminD2_ngperml,VitaminD3_ngperml,primary_25OHD3_VitaminD3_ratio,active_1alpha25OH2D3_primary_25OHD3_ratio,R_24R25OH2D3_primary_25OHD3_ratio,beta_4beta25OH2D3_primary_25OHD3_ratio,primary_25OHD2_VitaminD2_ratio,active_1alpha25OH2D2_primary_25OHD2_ratio,VitDsufficiency,VitDdeficiency}
do
  for geneName in {CYP3A4,DHCR7,CYP2R1,CYP27B1,SULT2A1,UGT1A4,GC,VDR,CASR,CUBN,CYP24A1,RXRG,RXRB,RXRA}
  do
    locuszoom --ld "./locusZoom_temp/linkageDisequilibriumAndPvalues_forLocusZoom_${currentPhenotype}_${geneName}.txt" --metal "./Results_for_LocusZoom/readyToPlot_${currentPhenotype}_${geneName}_Results_for_LocusZoom.txt" --build hg19 --refgene "${geneName}" showRecomb=F theme="publication" title="${geneName} Locus" --prefix "./locusZoom_plots/LocusZoom_plot_${currentPhenotype}_${geneName}" --plotonly
  done
done

# delete the locusZoom temporary directory now that it is no longer needed
rm -r ./locusZoom_temp/

# make a file with variant IDs for calculating linkage disequilibrium between:
    # (1) cyp3a4*1G (rs2242480, 7:99361466:C:T) and rs4646440 (7:99360870:G:A)
          # rs4646440 (7:99360870:G:A) was the cyp3a4 variant that was most significantly associated with 25(OH)D3 concentration and cyp3a4*1G has some evidence of influencing 25(OH)D3 metabolism in other studies.
    # (2) CUBN rs58278501 (10:17008452:A:AT) and [10:16935394:CGTGTGTGT:C (rs34859798), 10:17047864:G:T (rs77520345), 10:16890650:A:G (rs7089377), 10:16888811:CA:C (rs59819645)]
          # rs58278501 (10:17008452:A:AT) was significantly associated with 25(OH)D sufficiency in logistic regression and the other 4 variants were significantly associated with 25(OH)D3 concentration in linear regression.
    # (3) CUBN rs201370313 (10:17144251:C:CTTAT) and [10:16935394:CGTGTGTGT:C (rs34859798), 10:17047864:G:T (rs77520345), 10:16890650:A:G (rs7089377), 10:16888811:CA:C (rs59819645)]
          # rs201370313 (10:17144251:C:CTTAT) was significantly associated with 25(OH)D sufficiency in logistic regression and the other 4 variants were significantly associated with 25(OH)D3 concentration in linear regression.
Rscript --no-save ../scripts/genotype_phenotype_association/create_cubn_and_cyp3a4_linkage_snps_file.R
# calculate the linkage disequilibrium between:
  # (1) cyp3a4*1G (rs2242480, 7:99361466:C:T) and rs4646440 (7:99360870:G:A)
  # (2) CUBN rs58278501 (10:17008452:A:AT) and [10:16935394:CGTGTGTGT:C (rs34859798), 10:17047864:G:T (rs77520345), 10:16890650:A:G (rs7089377), 10:16888811:CA:C (rs59819645)]
  # (3) CUBN rs201370313 (10:17144251:C:CTTAT) and [10:16935394:CGTGTGTGT:C (rs34859798), 10:17047864:G:T (rs77520345), 10:16890650:A:G (rs7089377), 10:16888811:CA:C (rs59819645)]
vcftools --gzvcf "../quality_control_final_output_directory/vitamin_D_merged_clean.gt.vcf.gz" --snps CUBN_and_CYP3A4_SNPs_forLDcalc.txt --recode --recode-INFO-all --out temp
plink --vcf temp.recode.vcf --r2 --ld-snp-list CUBN_and_CYP3A4_reference_SNP_forLDcalc.txt --ld-window 4000 --ld-window-kb 1000000 --ld-window-r2 0
rm -r temp.recode.vcf
# rename the resulting linkage disequilibrium file
mv plink.ld significant_cubn_and_cyp3a4_variant_LD_calculations.txt
# Add rs-numbers to the variants in the significant_cubn_and_cyp3a4_variant_LD_calculations.txt file
Rscript --no-save ../scripts/genotype_phenotype_association/add_rsNumbers_to_significant_cubn_and_cyp3a4_variant_LD_calculations.R

# make a directory for storing rare variant analysis intermediates
mkdir ./rareVariantAnalysisIntermediates/

# create a file for updating the gender of study participants in the PLINK file
# this will create a file called ParticipantGenderPLINKfile.txt
Rscript --no-save ../scripts/genotype_phenotype_association/create_file_for_updating_gender.R

# calculate allele frequencies for all rare variants that are included in the rare variant analysis with PLINK (variants with a minimum minor allele count of at least 1 and a maximum minor allele count of 4)
vcftools --vcf ../annovar_annotation/vitaminD-merge_clean_updateINFO_allVariants.vcf --exclude ../annotation_output/variants_without_RefSeq_annotations.txt --mac 1 --max-mac 4 --recode --recode-INFO-all --out temp_vcf
plink --vcf temp_vcf.recode.vcf --keep-allele-order --freq
rm temp_vcf*

# filter the VCF file to a file with a max minor allele count equal to 4 and a minimum minor allele count of at least 1
# update the gender of the participants in the file also
# and update the phenotypes in the file with continuous phenotype data for each individual phenotype
for currentPhenotype in {active_1alpha25OH2D2_ngperml,active_1alpha25OH2D3_ngperml,beta_4beta25OH2D3_ngperml,Overall_primary_25OHD_ngperml,primary_25OHD2_ngperml,primary_25OHD3_ngperml,R_24R25OH2D3_ngperml,VitaminD2_ngperml,VitaminD3_ngperml,primary_25OHD3_VitaminD3_ratio,active_1alpha25OH2D3_primary_25OHD3_ratio,R_24R25OH2D3_primary_25OHD3_ratio,beta_4beta25OH2D3_primary_25OHD3_ratio,primary_25OHD2_VitaminD2_ratio,active_1alpha25OH2D2_primary_25OHD2_ratio}
do
  vcftools --vcf ../annovar_annotation/vitaminD-merge_clean_updateINFO_allVariants.vcf --exclude ../annotation_output/variants_without_RefSeq_annotations.txt --mac 1 --max-mac 4 --recode --recode-INFO-all --out temp_vcf
  plink --vcf temp_vcf.recode.vcf --keep-allele-order --update-sex ParticipantGenderPLINKfile.txt --pheno "./pheno_files_for_rare_variant_analysis/${currentPhenotype}_plink_pheno_for_rareVariantAnalysis.txt" --make-bed --out "./rareVariantAnalysisIntermediates/rare_variants_for_SKAT_analysis_${currentPhenotype}"
  rm temp_vcf*
done

# filter the VCF file to a file with a max minor allele count equal to 4 and a minimum minor allele count of at least 1
# update the gender of the participants in the file also
# and update the phenotypes in the file with binary phenotype data 
# The --1 flag specifies case=1 control=0 phenotype
for currentPhenotype in {VitDsufficiency,VitDdeficiency}
do
  vcftools --vcf ../annovar_annotation/vitaminD-merge_clean_updateINFO_allVariants.vcf --exclude ../annotation_output/variants_without_RefSeq_annotations.txt --mac 1 --max-mac 4 --recode --recode-INFO-all --out temp_vcf
  plink --vcf temp_vcf.recode.vcf --keep-allele-order --update-sex ParticipantGenderPLINKfile.txt --pheno "./pheno_files_for_rare_variant_analysis/${currentPhenotype}_plink_pheno_for_rareVariantAnalysis.txt" --1 --make-bed --out "./rareVariantAnalysisIntermediates/rare_variants_for_SKAT_analysis_${currentPhenotype}"
  rm temp_vcf*
done

# using the plink files for rare variant analysis computed in the previous step, run the SKAT-O rare variant association model for all phenotypes
Rscript --no-save ../scripts/genotype_phenotype_association/SKATO_analysis.R

# delete the rare variant analysis intermediate directory now that it is no longer needed
rm -r ./rareVariantAnalysisIntermediates/



