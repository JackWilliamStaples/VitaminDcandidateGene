# make a working directory for using pypgx if it does not exist
if [[ -d ~/Vitamin_D_project_MCCF/Vitamin_D_project_MCCF/pypgx ]]
then
  echo "pypgx haplotype assignment directory exists"
else
  echo "Created working directory for star allele haplotype assignment with pypgx"
  mkdir ~/Vitamin_D_project_MCCF/Vitamin_D_project_MCCF/pypgx/
fi

# working directory for this is:
cd ~/Vitamin_D_project_MCCF/Vitamin_D_project_MCCF/pypgx/

# make sure that python version 3 is running (this is required for pypgx)
# add the pyenv executable to the PATH
export PATH="$HOME/.pyenv/bin:$PATH"
eval "$(pyenv init -)"
eval "$(pyenv virtualenv-init -)"

# set the current version of python as python 3
pyenv global 3.8.2

# move vcf files to the current working directory
mv ../pgx_pop_annotation/vitaminD-merge_clean_updateINFO_allVariants.vcf.gz ./vitaminD-merge_clean_updateINFO_allVariants.vcf.gz
mv ../pgx_pop_annotation/vitaminD-merge_clean_updateINFO_allVariants.vcf.gz.tbi ./vitaminD-merge_clean_updateINFO_allVariants.vcf.gz.tbi

# call star alleles with pypgx for all possible genes
for currentGene in {CYP2R1,CYP3A4,UGT1A1,UGT1A4}
do
  pypgx run-ngs-pipeline "${currentGene}" "grch37-${currentGene}-pipeline" --variants vitaminD-merge_clean_updateINFO_allVariants.vcf.gz
done

# set the current version of python as the masOS system version (python 2)
pyenv global system 

# uncompress the pypgx results file for each gene
for currentGene in {CYP2R1,CYP3A4,UGT1A1,UGT1A4}
do
  cd "grch37-${currentGene}-pipeline/"
  unzip results.zip 
  rm *zip
  mv */data.tsv ./data.tsv
  cd ../
done

# summarize the pypgx results
Rscript --no-save ../scripts/genetic_characterization/summarize_pypgx_results.R

# move vcf files back to their original directory
mv vitaminD-merge_clean_updateINFO_allVariants.vcf.gz ../pgx_pop_annotation/vitaminD-merge_clean_updateINFO_allVariants.vcf.gz
mv vitaminD-merge_clean_updateINFO_allVariants.vcf.gz.tbi ../pgx_pop_annotation/vitaminD-merge_clean_updateINFO_allVariants.vcf.gz.tbi


