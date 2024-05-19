# this is a script for calling star alleles with pgx-pop

# make a working directory for using pgx-pop if it does not exist
if [[ -d ~/Vitamin_D_project_MCCF/Vitamin_D_project_MCCF/pgx_pop_annotation ]]
then
  echo "PGx-POP directory exists"
else
  echo "Created working directory for star allele haplotype assignment with PGx-POP"
  mkdir ~/Vitamin_D_project_MCCF/Vitamin_D_project_MCCF/pgx_pop_annotation/
fi

# working directory for this is:
cd ~/Vitamin_D_project_MCCF/Vitamin_D_project_MCCF/pgx_pop_annotation/
# add the directory with pgx-pop python scripts to the path
PATH=$PATH:~/bioinformatics_software/PGxPOP/bin/

# add the pyenv executable to the PATH for ensuring that the python version 3 is running for PGx-POP
export PATH="$HOME/.pyenv/bin:$PATH"
eval "$(pyenv init -)"
eval "$(pyenv virtualenv-init -)"

# set the current version of python as python 3
pyenv global 3.8.2

# compress and zip the VCF file with all variants to be used for PGx-POP
bgzip -c ../annovar_annotation/vitaminD-merge_clean_updateINFO_allVariants.vcf > vitaminD-merge_clean_updateINFO_allVariants.vcf.gz
tabix -p vcf vitaminD-merge_clean_updateINFO_allVariants.vcf.gz

# run PGxPOP with python 3 (python 3.6 or higher required for PGx-POP) for UGT1A1
python ~/bioinformatics_software/PGxPOP/bin/PGxPOP.py --vcf vitaminD-merge_clean_updateINFO_allVariants.vcf.gz -g UGT1A1 --phased --build hg19 -o ./PGx_POP_results_allVariants_UGT1A1.txt

# set the current version of python as the masOS system version (python 2)
pyenv global system 

# summarize the PGx-POP results
Rscript --no-save ../scripts/genetic_characterization/summarize_PGxPOP_results.R
