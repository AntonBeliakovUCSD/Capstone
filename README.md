# How does Heritability Correlate With Causal Variance?

This project aims to answer the question how heritability correlates with causal variance. The initial hypothesis is that there is a positive linear relationship and as the number of causal snps increases for a gene id, the higher it's heritability is likely to be. This hypothesis would align with the assumption in existing literature that gene ids or traits with few causal variants are also not highly conserved. The linear regression, grouping, and visualizations align with the hypothesis. The GERP score analysis offers a key interesting result about gene ids with high heritability and low causal snps. The TWAS analysis is still underway, but the results have been calculated for one trait: breast carcinoma. The next steps are to scale the TWAS analysis and interpret its results in relation with our overall research question. 

## File Descriptions

### Data 
- **LDREF/**: Genotype Data
- **BreastCarcinoma.22.dat**: TWAS Output
- **GEUVADIS_EUR_covariates.txt**: Covariates File
- **causal_variance_snps.csv**: Causal Variance Scores and Number of Causal SNPs per gene id
- **gene_annotation.txt**: Supplementary gene id information
- **heritability_scores.csv**: Scores for each gene id
- **Summary Statistics**: [http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST007001-GCST008000/GCST007236/harmonised/](Link to Summary Statistics - Too large to be attatched)
- **GERP Scores**: [https://genome.ucsc.edu/cgi-bin/hgTables?db=hg19&hgta_group=compGeno&hgta_track=allHg19RS_BW&hgta_table=allHg19RS_BW&hgta_doSchema=describe+table+schema](Link to GERP Scores - Too large to be attatched)

### Scripts 
- **Data Viz, GERP Analysis, TWAS Interpretation.ipynb**: A working notebook that has exploratory data visualizations, GERP score analysis, and the beginning of a TWAS interpretation. 
- **FUSION.assoc_test.R**: TWAS Requirement
- **FUSION.compute_weights.R**: TWAS Requirement
- **FUSION.post_process.R**: TWAS Requirement
- **Heritability Scores and Regression.py**: Script that calculated heritability scores, cleaned them, and did the statistical analysis. 
- **Scaling-Susie.Rmd**: A R script to fine map all gene ids. 
- **susie_function.Rmd**: A R script that calculated the causal variance for one gene id. 

## To Recreate Analysis

This project used multiple languages and tools including Python, R, the GCTA command line tool, and plink2. Ensure that all languages are set up before running the notebooks or attempting to recreate the files from scratch. Once the analysis is completely finished, the team will work on integrating the different components to make the re-creation more seamless for a third party. 
