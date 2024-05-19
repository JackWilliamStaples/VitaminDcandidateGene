# This is a script for creating a text file with variant IDs for calculating linkage disequilibrium between:
  # (1) cyp3a4*1G (rs2242480, 7:99361466:C:T) and rs4646440 (7:99360870:G:A)
  # (2) CUBN rs58278501 (10:17008452:A:AT) and [10:16935394:CGTGTGTGT:C (rs34859798), 10:17047864:G:T (rs77520345), 10:16890650:A:G (rs7089377), 10:16888811:CA:C (rs59819645)]
  # (3) CUBN rs201370313 (10:17144251:C:CTTAT) and [10:16935394:CGTGTGTGT:C (rs34859798), 10:17047864:G:T (rs77520345), 10:16890650:A:G (rs7089377), 10:16888811:CA:C (rs59819645)]

# load necessary packages
source("../scripts/load_R_packages.R")

# define the variants that will be used for calculating pairwise linkage disequilibrium
CUBN_and_CYP3A4_SNPs_forLDcalc <-
  data.frame(
    "snps" = c(
               # CYP3A4 variants
               "7:99361466:C:T",
               "7:99360870:G:A",
               # CUBN variants
               "10:17008452:A:AT",
               "10:17144251:C:CTTAT",
               "10:16935394:CGTGTGTGT:C",
               "10:17047864:G:T",
               "10:16890650:A:G",
               "10:16888811:CA:C",
               # UGT1A variants
               "2:234652740:T:C",
               "2:234627608:T:G"
               )
    )

# define the reference variants for calculating pairwise linkage disequilibrium
CUBN_and_CYP3A4_reference_SNP_forLDcalc <-
  data.frame(
    "referenceSNP" = c(
                      # CYP3A4 variant rs4646440 (7:99360870:G:A)
                      "7:99360870:G:A",
                      # CUBN variants rs58278501 (10:17008452:A:AT) and rs201370313 (10:17144251:C:CTTAT)
                      "10:17008452:A:AT",
                      "10:17144251:C:CTTAT",
                      # UGT1A variant rs139116240 (2:234652740:T:C)
                      "2:234652740:T:C"
                      )
    )

# save both files to the file system
write.table(
  x = CUBN_and_CYP3A4_SNPs_forLDcalc,
  file = "./CUBN_and_CYP3A4_SNPs_forLDcalc.txt",
  col.names = FALSE,
  row.names = FALSE,
  quote = FALSE
  )

write.table(
  x = CUBN_and_CYP3A4_reference_SNP_forLDcalc,
  file = "./CUBN_and_CYP3A4_reference_SNP_forLDcalc.txt",
  col.names = FALSE,
  row.names = FALSE,
  quote = FALSE
  )

print(x = "created CUBN_and_CYP3A4_SNPs_forLDcalc.txt")
print(x = "created CUBN_and_CYP3A4_reference_SNP_forLDcalc.txt")
