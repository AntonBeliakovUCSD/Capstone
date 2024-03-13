# Dissecting Heritability and Causal Variants in Cancer Genomics

This project dissects the genetic architecture, specifically heritability and causal variants, of common cancers to further understanding of significant genes that cause these cancers and the relationship between causal variants and heritability more broadly. Three primary methods were used to 1) significantly test the impact of heritability and the number of causal single nucleotide polymorphisms (SNPs) across all genes available 2) integrate genomic evolutionary rate profiling scores into fine mapping results to further evolutionary detail for all genes 3) conduct a transcriptome wide association study on particular cancers. The results of the first step were used as inputs into the final step to corroborate the results and provide additional insights into the significance of the causal genes. Our results built upon existing empirical studies by specifying which particular genes in a gene family are causal and by investigating how fine mapping results could corroborate existing findings. 

## File Descriptions

### Data 

Many data files are zipped due to their size. Please unzip all before running any analysis to avoid errors. Download all data and make any appropriate name changes according to your own preferences. 

- **LDREF/**: Genotype Data
- **GEUVADIS_EUR_covariates.txt**: Covariates File
- **causal_variance_snps.csv**: Causal Variance Scores and Number of Causal SNPs per gene id
- **gene_annotation.txt**: Supplementary gene id information
- **PASS_{disease}.sumstats.zip**: Summary statistics for breast cancer, ovarian cancer, and prostate cancer in 3 files. 
- **TCGA_{disease}.tar.bz2**: 4 Weight Tissue files for cancers. 
- **UKB_460K.cancer_MELANOMA.sumstats.zip**: Summary statistics for Melanoma. 
- **Updated_{disease}.TUMOR.zip**: 4 Weight Tissue files for cancers with SuSiE model. 
- **heritability_scores.csv**: Scores for each gene id
- **Summary Statistics**: [Link to Summary Statistics - Too large to be attatched](http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST007001-GCST008000/GCST007236/harmonised/)
- **GERP Scores**: [Link to GERP Scores - Too large to be attatched](https://genome.ucsc.edu/cgi-bin/hgTables?db=hg19&hgta_group=compGeno&hgta_track=allHg19RS_BW&hgta_table=allHg19RS_BW&hgta_doSchema=describe+table+schema)
- **Gene Expression Data**: Phenotype Data Hosted on Google Drive Since it is too Large to be Attatched. [Drive Link](https://drive.google.com/file/d/1lJwqtsWGMn9kLe_YIRjj5V8_7Nj3UlPt/view?usp=share_link)


### Scripts 
- **Analyzing TWAS Results.Rmd**: A R script that analyzes TWAS results with different models. This saves all the plots for all the cancers in this analysis. 
- **EDA and Visualization + GERP.ipynb**: A working notebook that has exploratory data visualizations, GERP score analysis, and the beginning of a TWAS interpretation. 
- **FUSION.assoc_test.R**: TWAS Requirement
- **FUSION.compute_weights.R**: TWAS Requirement
- **FUSION.post_process.R**: TWAS Requirement
- **Heritability Scores and Regression.py**: Script that calculated heritability scores, cleaned them, and did the statistical analysis. 
- **Scaling-Susie.Rmd**: A R script to fine map all gene ids. 
- **Weight Incorporation.Rmd**: A R script that creates updated weight files for all cancers. 
- **clean_summary_stats.py**: A Python script that cleans summary statistics for a particular trait if you wish to run this for a different disease not analyzed in this project. 
- **susie_function**: A R script for fine mapping any one particular gene. 
- **twas_all_chromosome.sh**: A script to run the chromosome level TWAS command line prompt for all chromosomes in this analysis. 
- **twas_scripts.py**: A Python script to generate run TWAS and generate its corresponding plots for any particular disease. Can be run on command line with appropriate paramaters for flexability for traits not in this analysis. 

## Output

In the output folder, there are multiple TWAS results file specific to each disease, chromosome, and what models were used in it. These results were compilied into 2 aggregrate files "full_twas_susie" and "full_twas_default." There are also Maimi and Manhatten plots for each cancer. Here is an example of the Manhatten plot that shows the significant genes for ovarian cancer: 

![Example](/output/ovarian_cancer.png 

## To Recreate TWAS Analysis

The bulk of the reproducable analysis in our project will be the TWAS, where the user can specify a trait and see a plot with significant SNPs. 

You can also recreate the regression and GERP analysis that drove the TWAS by running "EDA and Visualization + GERP.ipynb."

This project used multiple languages and tools including Python, R, the [GCTA](https://yanglab.westlake.edu.cn/software/gcta/#Overview) command line tool, and [Plink](https://www.cog-genomics.org/plink/2.0/). Ensure that all languages are set up before running the notebooks or attempting to recreate the files from scratch. Note that some packages and dependencies vary by operating system. For example, the R scripts used a package that is not compatible with Windows operating systems. 

Use the following command to clone the repository: 
```
git clone https://github.com/AntonBeliakovUCSD/Capstone.git

cd Capstone 

pip3 install -r requirements.txt
```

To Run The TWAS for all chromosomes for all the cancers in our analysis, follow the following steps. 

1. Unzip the Weight Files and Prepare Summary Statistics
The weight files for this analysis include TCGA-BRCA.TUMOR.tar.bz2, TCGA-OV.TUMOR.tar.bz2, TCGA-PRAD.TUMOR.tar.bz2, and TCGA-SKCM.TUMOR.tar.bz2. The summary statistics start with the prefix "PASS" or "UKB."


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

4. Optional: Rerun the TWAS with a SuSiE Model to integrate Fine Mapping Results

In this step, all the weight files were updated to include only the weights from the genes that had appeared in a credible set in the fine mapping portion of the project. You do not need to rerun the weights portion of this analysis (the weights were updated in the "Weight_Incorporation.Rmd"). The updated weights are in the data folder as: 

- data/Updated_BRCA.TUMOR
- data/Updated_SKCM.TUMOR
- data/Updated_PRAD.TUMOR
- data/Updated_OV.TUMOR

You can rerun each TWAS command above and replace the path/to/weights_dir paramter with an updated version. It's also advisable to change the output name to not override previous results.  

Here is an updated example for reruning the TWAS on Melanoma with the susie model. 

```
./src/twas_all_chromosomes.sh data/UKB_460K.cancer_MELANOMA.sumstats data/TCGA-SKCM.TUMOR/TCGA-SKCM.TUMOR.pos data/Updated_SKCM.TUMOR Melanoma_susie output 
```
5. Running TWAS for Your Own Trait! 

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

## R Shiny App

Check out interactive Manhatten plots and gene searcher on our app (which is also integrated on our website) [here](https://dissecting-cancer-genomics.shinyapps.io/capstone/)! The code for this is in the "capstone" folder. 
