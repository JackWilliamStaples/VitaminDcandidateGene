# Vitamin D candidate-gene study pipeline

##### **Author:** Jack Staples

**Overview:**

This is a repository for the code for a candidate-gene study of how genetic variation in vitamin D regulatory genes (*CASR*, *CUBN*, *CYP2R1*, *CYP3A4*, *CYP24A1*, *CYP27B1*, *DHCR7*, *GC*, *RXRA*, *RXRB*, *RXRG*, *SULT2A1*, *UGT1A*, and *VDR*), season, and demographic factors (age, body mass index [BMI], and gender) influence vitamin D phenotypes.

The vitamin D phenotypes included in this analysis were concentration of vitamin D [vitamin D3 and D2], vitamin D metabolites [25(OH)D3 and D2; 1-alpha,25(OH)2D3 and D2; 4-beta,25(OH)2D3; and 24R,25(OH)2D3], metabolite ratios [25(OH)D3/Vitamin D3; 1-alpha,25(OH)2D3/25(OH)D3; 4-beta,25(OH)2D3/25(OH)D3; and 24R,25(OH)2D3/25(OH)D3], and vitamin D status [25(OH)D sufficiency (20 ng/mL) and 25(OH)D deficiency (12 ng/mL)].

This pipeline consists of a combination of unix bioinformatics command line tools and R scripts wrapped in bash scripts so the pipeline can be run on a bash shell terminal.

**Required software:**

-   ANNOVAR version date 06/08/2020 (<https://annovar.openbioinformatics.org/en/latest/>)

-   BEAGLE version 4.1 (<https://faculty.washington.edu/browning/beagle/b4_1.html>)

-   BCFtools version 1.11 from htslib version 1.11 (<http://www.htslib.org/download/>)

-   bgzip version 1.11 (included with htslib)

-   Ensembl variant effect predictor (VEP) web interface (<http://grch37.ensembl.org/Homo_sapiens/Tools/VEP>)

-   LocusZoom version 1.3 (<https://genome.sph.umich.edu/wiki/LocusZoom_Standalone>)

-   PGx-POP version 1.0 (<https://github.com/PharmGKB/PGxPOP>)

-   PLINK version 1.90 (<https://www.cog-genomics.org/plink/>)

-   Pyenv version 1.2.23 for switching between python version 2 and 3 (<https://github.com/pyenv/pyenv>)

-   pypgx version 0.16.0 (<https://github.com/sbslee/pypgx>)

-   Python version 2.7.18 and 3.8.2

-   R version 4.0.4

-   SNPnexus web interface (<https://www.snp-nexus.org/v4/>)

-   stargazer version 1.0.8 (<https://stargazer.gs.washington.edu/stargazerweb/index.html>)

-   tabix (included with htslib)

-   VCFtools version 0.1.16 (<https://vcftools.github.io/index.html>)

**Required R packages:**

-   aod version 1.3.1

-   cosinor version 1.1

-   data.table version 1.14.0

-   estimatr version 0.30.2

-   factoextra version 1.0.7

-   GGally version 2.1.2

-   ggfortify version 0.4.11

-   ggpubr version 0.4.0

-   glue version 1.4.2

-   mltools version 0.3.5

-   olsrr version 0.5.3

-   parallel version 4.0.4

-   readxl version 1.3.1

-   reshape2 version 1.4.4

-   SKAT version 2.0.1

-   tidyverse version 1.3.0

-   viridis version 0.5.1

**The input files for this pipeline are:**

-   Variant call format (VCF) file from sequencing of the first subset of research participants (woodahl_grc_1.HF.final.vcf)

-   VCF file from sequencing of the second subset of research participants (woodahl_grc_custom_2.HF.final.vcf)

-   A BED file with genomic position ranges of the genes that were sequenced in this study (gene_list.bed.txt)

-   A 1000 genomes reference genome FASTA file (hs37d5.fa) downloaded from the 1000 genomes ftp (<ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz>)

-   Phenotype files with concentrations of vitamin D and vitamin D metabolites for all study participants (Vit_D\_Master_clean.csv and UM_Vit_D\_data_summary_A417_toA560_clean.csv)

-   De-identified covariates files with demographic information (Encounters-Grid view.csv and Participant info-Grid view_deidentified.csv)

-   Files to match research participant visit codes with gene sequencing IDs (lookup_woodahl_grc_1\_updated.xlsx and lookup_woodahl_grc_custom_2.xlsx)

-   Star allele haplotype data for *CYP3A4* downloaded from the pharmacogene variation consortium (PharmVar) for the GRCh37 reference genome sequence (CYP3A4.NC_000007.13.haplotypes.tsv) (link: <https://www.pharmvar.org/gene/CYP3A4>)

-   Pharmacogenomics Knowledgebase (PharmGKB) clinical variant annotation data (downloaded from: <https://api.pharmgkb.org/v1/download/file/data/clinicalAnnotations.zip>)

-   Estimated ultraviolet light B (UVB) radiation data in Polson, MT for a 12 month year (UVB_radiation_per_month.csv) with the online Quick TUV calculator (<https://www.acom.ucar.edu/Models/TUV/Interactive_TUV/>)

**Pipeline scripts and descriptions:**

-   **Vitamin_D\_genetic_analysis_pipeline_master.R** lists the order that pipeline scripts need to be run in order to complete the analysis

    This runs the following scripts:

    -   **Vitamin_D\_quality_control.sh** performs quality control specific for candidate-gene studies and statistical phasing of variant calling data.

    -   **Vitamin_D\_vcf_summary_stats.sh** calculates basic summary statistics for variant calling data, creates VCF files of singleton and INDEL variants, and creates VCF files for individuals that were phenotypic outliers.

    -   **ANNOVAR_annotation.sh** annotates variant calling data with ANNOVAR, the Ensembl variant effect predictor, and SNP-nexus (includes variant annotation from multiple databases); clinical variant annotation data from the pharmGKB; and star allele haplotypes from pharmVAR. This also performs *in-silico* predictions of deleterious variants using the ADME optimized prediction framework and identifies deleterious variants that individuals who are phenotypic outliers have. This also runs a script for plotting some of the results of the annotation.

    -   **pgx_pop_star_allele_caller.sh** assigns star allele haplotypes to variant calling data with PGx-POP software.

    -   **stargazer_vitaminD.sh** assigns star allele haplotypes to variant calling data with stargazer software and assigns star allele haplotypes for the gene *CYP3A4* using a custom method.

    -   **pypgx_star_allele_caller.sh** assigns star allele haplotypes to variant calling data with pypgx software.

    -   **consolidate_starAllele_data.R** this aggregates star allele haplotype calling data from the previous 3 scripts into a results table.

    -   **compute_variant_annotation_summaries.R** computes variant summary statistics tables.

    -   **candidate_gene_association_analysis.sh** performs all the genotype-phenotype statistical association analysis, which includes:

        -   Linear regressions to test individual genetic variants, season, and demographic factors for associations with continuous vitamin D phenotypes.

        -   Logistic regressions to test individual genetic variants, season, and demographic factors for associations with dichotomous vitamin D phenotypes.

        -   Linkage disequilibrium pruning of association results.

        -   Vitamin D binding protein haplotype/diplotype assignment and testing for association with vitamin D phenotypes.

        -   Plots of estimated UVB radiation data in Polson, MT throughout the year.

        -   Boxplots of continous phenotypes stratified by genotype for significant variants.

        -   LocusZoom plots of association results for each variant tested.

        -   Pooled rare variant analysis using the optimal combination of the rare variant burden and sequence kernel association tests (SKAT-O method).

    -   **load_R\_packages.R** is run at the beginning of every R script for loading necessary R packages.
