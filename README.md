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

Please note that you have to download the summary statistics such as bc.h.tsv from the GWAS Catelog. They cannot be uploaded to GitHub due to size. 

This step will take in the raw summary statistics and output a prepared file for TWAS analysis. This step can take minutes. 
```
python3 src/clean_summary_stats.py data/bc.h.tsv "['variant_id', 'effect_allele', 'other_allele', 'Z']" data/BreastCarcinoma_GWAS.txt
```

Run the TWAS command in Terminal
```
Rscript src/FUSION.assoc_test.R \
--sumstats data/BreastCarcinoma_GWAS.txt_summary_stats.txt \
--weights ./data/BRCA/TCGA-BRCA.TUMOR.pos \
--weights_dir ./data/BRCA/ \
--ref_ld_chr ./data/LDREF/1000G.EUR. \
--chr 22 \
--out data/BreastCarcinoma.22.dat
```

Clean the TWAS Output. This output can be used for a variety of interpretations and significance testing. A Manhatten plot with labeled SNPs will also be shown! 

```
python3 src/twas_scripts.py data/BreastCarcinoma.22.dat
```

Once the TWAS analysis is completely finished, the team will work on integrating the different components to make the re-creation more seamless for a third party. Currently, one should be able to run the code and follow the reading structure. However, this will be more generalized so that one can run a TWAS by inserting their own summary statistics file and specifying a few parameters. 
