# this is a script for calling star alleles with Stargazer (ref: Lee, et al., 2018) on the QC'd VCF with all variants for cyp2r1 and ugt1a4
# this script also calls star alleles for cyp3a4 with a custom method

# make a working directory for using stargazer if it does not exist
if [[ -d ~/Vitamin_D_project_MCCF/Vitamin_D_project_MCCF/starGazer ]]
then
  echo "stargazer haplotype assignment directory exists"
else
  echo "Created working directory for star allele haplotype assignment with stargazer"
  mkdir ~/Vitamin_D_project_MCCF/Vitamin_D_project_MCCF/starGazer/
fi

# working directory for this is:
cd ~/Vitamin_D_project_MCCF/Vitamin_D_project_MCCF/starGazer/

# make sure that python version 3 is running (this is required for starGazer)
# add the pyenv executable to the PATH
export PATH="$HOME/.pyenv/bin:$PATH"
eval "$(pyenv init -)"
eval "$(pyenv virtualenv-init -)"

# set the current version of python as python 3
pyenv global 3.8.2

#### running stargazer with the VCF file from the quality control output
# run stagazer genotyping tool for genes cyp2r1 and ugt1a4 with the vcf with the targeted sequencing flag
gunzip ../quality_control_final_output_directory/vitamin_D_merged_clean.gt.vcf.gz
for currentGene in {cyp2r1,cyp3a4,ugt1a1,ugt1a4}
do 
  python ~/bioinformatics_software/Stargazer_v1.0.8/stargazer.py genotype -o "${currentGene}_results" -d ts -t "${currentGene}" --vcf "../quality_control_final_output_directory/vitamin_D_merged_clean.gt.vcf"
done
		
# set the current version of python as the masOS system version (python 2)
pyenv global system 

# plot (make haplotype and genotype-predicted phenotype bar charts and tables for each gene) with the stargazer star allele calling results 
Rscript --no-save ../scripts/genetic_characterization/summarize_stargazer_results.R

# use BCFtools to filter the VCF file of all variants to the CYP3A4 protein region (chr 7: 99354582-99381811)
# and create a sample x haplotype matrix for CYP3A4 only 
bgzip -c ../quality_control_final_output_directory/vitamin_D_merged_clean.gt.vcf > ../quality_control_final_output_directory/vitamin_D_merged_clean.gt.vcf.gz
tabix -p vcf ../quality_control_final_output_directory/vitamin_D_merged_clean.gt.vcf.gz
bcftools view "../quality_control_final_output_directory/vitamin_D_merged_clean.gt.vcf.gz" --regions 7:99354582-99381811 -o temp.vcf; bcftools convert --hapsample --vcf-ids temp.vcf -o "CYP3A4_haplotypes"
rm -r temp.vcf
gunzip CYP3A4_haplotypes.hap.gz

# assign CYP3A4 star alleles to the vitamin D variant calling data using a custom script 
# since stargazer does not include important star alleles such as CYP3A4*1G
Rscript --no-save ../scripts/genetic_characterization/assign_CYP3A4_star_allele_haplotypes.R



