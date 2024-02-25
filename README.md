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
- **Summary Statistics**: [Link to Summary Statistics - Too large to be attatched - rename it to bc](http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST007001-GCST008000/GCST007236/harmonised/)
- **GERP Scores**: [Link to GERP Scores - Too large to be attatched](https://genome.ucsc.edu/cgi-bin/hgTables?db=hg19&hgta_group=compGeno&hgta_track=allHg19RS_BW&hgta_table=allHg19RS_BW&hgta_doSchema=describe+table+schema)
- **BRCA**: These are the weights for the TWAS analysis, please download the appropriate ones from this link. These could not be pushed onto github. [Link - Download the Breast Invasive Carcinoma folder and rename it BRCA] (http://gusevlab.org/projects/fusion/) 

### Scripts 
- **Data Viz, GERP Analysis, TWAS Interpretation.ipynb**: A working notebook that has exploratory data visualizations, GERP score analysis, and the beginning of a TWAS interpretation. 
- **FUSION.assoc_test.R**: TWAS Requirement
- **FUSION.compute_weights.R**: TWAS Requirement
- **FUSION.post_process.R**: TWAS Requirement
- **Heritability Scores and Regression.py**: Script that calculated heritability scores, cleaned them, and did the statistical analysis. 
- **Scaling-Susie.Rmd**: A R script to fine map all gene ids. 
- **susie_function.Rmd**: A R script that calculated the causal variance for one gene id. 

## To Recreate TWAS Analysis

The bulk of the reproducable analysis in our project will be the TWAS, where the user can specify a trait and see a plot with significant SNPs. 

This project used multiple languages and tools including Python, R, the GCTA command line tool, and plink2. Ensure that all languages are set up before running the notebooks or attempting to recreate the files from scratch. Use the following command: 
```
git clone https://github.com/AntonBeliakovUCSD/Capstone.git

cd Capstone 

pip3 install -r requirements.txt
```

To Run The TWAS for all chromosomes for all the cancers in our analysis, follow the following steps. 

1. Unzip the Weight Files and Prepare Summary Statistics
The weight files for this analysis include TCGA-BRCA.TUMOR.tar.bz2, TCGA-OV.TUMOR.tar.bz2, TCGA-PRAD.TUMOR.tar.bz2, and TCGA-SKCM.TUMOR.tar.bz2. 

Ensure that your summary statistics are in the data folder. FIX THIS LATER - how to include summary statistics zipped. 
2. Run the TWAS Command in CMD or Terminal

This is the general format of the command in case you want to run it for another cancer not specified in this analysis. 

```
./src/twas_all_chromosomes.sh path/to/sumstats path/to/weights path/to/weights_dir desired_name path/to/output
```

These are the commands to rerun for this particular analysis.

Skin Cutaneous Melanoma
```
./src/twas_all_chromosomes.sh data/UKB_460K.cancer_MELANOMA.sumstats data/TCGA-SKCM.TUMOR/TCGA-SKCM.TUMOR.pos data/TCGA-SKCM.TUMOR Melanoma output 
```
Breast Invasive Carcinoma
```
./src/twas_all_chromosomes.sh data/PASS_BreastCancer.sumstats data/TCGA-BRCA.TUMOR/TCGA-BRCA.TUMOR.pos data/TCGA-BRCA.TUMOR breast_carcinoma output 
```
Prostate Adenocarcinoma
```
./src/twas_all_chromosomes.sh data/PASS_ProstateCancer.sumstats data/TCGA-PRAD.TUMOR/TCGA-PRAD.TUMOR.pos data/TCGA-PRAD.TUMOR prostate_cancer output 
```

Note that many genes were skipped, especially on chromosome 6,9,
12 for prostate cancer. 

Ovarian Serous Cystadenocarcinoma
```
./src/twas_all_chromosomes.sh data/PASS_OvarianCancer.sumstats data/TCGA-OV.TUMOR/TCGA-OV.TUMOR.pos data/TCGA-OV.TUMOR ovarian_cancer output 
```
3. Analyze Results in R

Run the script "Analyzing TWAS Results.Rmd" and Manhattan plots and Miami plots will be saved for each cancer in the output folder. 

How should we improve this?

Make the script executable on it's own? o
