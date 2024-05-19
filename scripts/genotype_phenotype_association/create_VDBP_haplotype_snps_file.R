# this is a script for creating a text file with variant IDs for the vitamin D haplotype defining SNPs:
# rs7041 (4:72618334:A:C) and rs4588 (4:72618323:G:T)

# load necessary packages
source("../scripts/load_R_packages.R")


VDBP_haplotype_SNPs_forLDcalc <-
data.frame("snps" = c("4:72618334:A:C","4:72618323:G:T"))

VDBP_haplotype_reference_SNP_forLDcalc <-
  data.frame("referenceSNP" = c("4:72618334:A:C"))

# save both files to the file system
write.table(
  x = VDBP_haplotype_SNPs_forLDcalc,
  file = "./VDBP_haplotype_SNPs_forLDcalc.txt",
  col.names = FALSE,
  row.names = FALSE,
  quote = FALSE
  )

write.table(
  x = VDBP_haplotype_reference_SNP_forLDcalc,
  file = "./VDBP_haplotype_reference_SNP_forLDcalc.txt",
  col.names = FALSE,
  row.names = FALSE,
  quote = FALSE
  )

print(x = "created VDBP_haplotype_SNPs_forLDcalc.txt")
print(x = "created VDBP_haplotype_reference_SNP_forLDcalc.txt")