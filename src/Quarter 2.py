#!/usr/bin/env python
# coding: utf-8

# ## Generate Heritability Scores and Perform Linear Regression Between the Scores and the Number of Causal SNPs

# In[1]:


import pandas_plink as pp
import statsmodels.api as sm
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import glob
import os
import statsmodels.api as sm


# In[14]:


main_dir = os.getcwd()[:-3]
dire = main_dir + "data/"


# ### Get the genes we want to analyze

# In[2]:


genes = pd.read_csv(dire+"causal_variance_cv_10.csv")
genes = genes.drop(genes.columns[0], axis=1)
genes = genes[genes['CausalVariance'].notna()]
gene_id = genes["GeneID"]
genes


# In[3]:


gene_id


# In[4]:


causal_snp = pd.read_csv(dire+"causal_variance_snps.csv").drop(columns = "Unnamed: 0")
causal_snp = causal_snp[causal_snp['CausalVariance'].notna()]
causal_snp


# ### Read in the Genotype files (.bim, .bed, .fam) and the Phenotype file (gene expression data for practice)

# In[5]:


gene_expression_data = pd.read_csv(dire+'GD462.GeneQuantRPKM.50FN.samplename.resk10.txt', delimiter='\t')
gene_expression_data['TargetID'] = gene_expression_data['TargetID'].str.split('.').str[0]
gene_expression_data


# #### Filter out the genes we are interested in

# In[6]:


gene_expression_data = gene_expression_data[gene_expression_data["TargetID"].isin(list(gene_id))]
gene_expression_data


# ### Read in the genotype data for further filtering

# In[7]:


bfile = dire+"LDREF/1000G.EUR.1"
(bim, fam, bed) = pp.read_plink(bfile)


# In[8]:


fam


# #### Find the common individuals and create feasible phenotype data files for GCTA to process

# In[8]:


fam_ids = fam['fid'].tolist()

gene_expression_filtered = gene_expression_data[['TargetID'] + [col for col in gene_expression_data.columns if col in fam_ids]]
gene_expression_filtered


# ### Change the order of columns in the gene_expression_data to adhere to the fam file

# In[9]:


# Ensure 'TargetID' is the first column in gene_expression_data
first_col = gene_expression_filtered.pop('TargetID')

# Get the order of 'fid' from the fam DataFrame
order = [fid for fid in fam_ids if fid in gene_expression_data.columns]

# Reorder the columns in gene_expression_data to match the order in 'fid'
gene_expression_filtered = gene_expression_filtered[order]

# Reinsert 'TargetID' as the first column
gene_expression_filtered.insert(0, 'TargetID', first_col)
gene_expression_filtered


# In[ ]:





# ## Use the h2 values from the FUSION website

# In[12]:


# Directory containing your .hsq files
hsq_directory = dire+"GTExv8.ALL.Cells_EBV-transformed_lymphocytes"  # Change this to the path where your .hsq files are stored

# Initialize a list to store your data
data = []

# Loop through each file in the directory
for filename in os.listdir(hsq_directory):
    if filename.endswith(".hsq"):
        gene_id = filename.split(".")[0]  # Assumes the gene ID is the part of the filename before the first dot
        filepath = os.path.join(hsq_directory, filename)
        
        # Open and read the .hsq file
        with open(filepath, 'r') as file:
            lines = file.readlines()
            # Assuming the heritability score is on the first line, following the first space
            # Modify this line if your file format is different
            hsq_score = lines[0].split()[1]  # This gets the second element from the first line
            
            # Append the gene ID and heritability score to the data list
            data.append([gene_id, hsq_score])

# Convert the list to a pandas DataFrame
df_hsq = pd.DataFrame(data, columns=['Gene_ID', 'Heritability_Score'])
df_hsq


# In[13]:


lst2 = df_hsq['Gene_ID'].tolist()
common_values = list(set(lst2) & set(gene_expression_filtered["TargetID"].tolist()))
len(common_values)


# In[14]:


df_hsq_renamed = df_hsq.rename(columns={"Gene_ID": "TargetID"})

# Perform an inner join to filter df_hsq to only include genes present in gene_expression_filtered
filtered_hsq = pd.merge(df_hsq_renamed, gene_expression_filtered[['TargetID']], on="TargetID", how="inner")
filtered_hsq


# In[15]:


# map the chromosome to each gene
chromosome_mapping = gene_expression_data.set_index('TargetID')['Chr'].to_dict()

filtered_hsq['Chromosome'] = filtered_hsq['TargetID'].map(chromosome_mapping)
filtered_hsq


# ### Change all negative heritability scores to 0

# In[16]:


filtered_hsq['Heritability_Score'] = pd.to_numeric(filtered_hsq['Heritability_Score'])
filtered_hsq['Heritability_Score'] = filtered_hsq['Heritability_Score'].apply(lambda x: max(x, 0))
filtered_hsq


# In[17]:


# save the dataframe to local
filtered_hsq.to_csv("heritability_scores.csv", index=False)


# In[ ]:





# In[18]:


filtered_hsq['Heritability_Score'].value_counts()


# ## Linear Regression

# #### Merge the dataframe first

# In[22]:


merged_df = pd.merge(causal_snp, filtered_hsq, left_on='GeneID', right_on='TargetID')
merged_df


# In[23]:


# add the log of NumCausalSNPs and heritability scores, adding 1e-6 to avoid log of 0
merged_df['Log_Heritability_Score'] = np.log(merged_df['Heritability_Score'] + 1e-6)
merged_df['Log_NumCausalSNPs'] = np.log(merged_df['NumCausalSNPs'] + 1e-6)

merged_df


# In[61]:


X = merged_df['NumCausalSNPs']  # Independent variable
y = merged_df['Heritability_Score']    # Dependent variable

X = sm.add_constant(X)

model = sm.OLS(y, X).fit()


# In[62]:


print(model.summary())


# In[ ]:





# In[43]:


# Plot it on a graph
plt.figure(figsize=(10, 6))
plt.scatter(merged_df['NumCausalSNPs'], merged_df['Heritability_Score'], alpha=0.5)
plt.title('Scatter Plot of Heritability Scores vs NumCausalSNPs')
plt.xlabel('NumCausalSNPs')
plt.ylabel('Heritability')
plt.grid(True)
plt.savefig("h2_NumCausalSNPs")
plt.show()


# In[ ]:





# #### Let's see if we change to "given that there is at least 1 causal variants by Susie, what is the h2"

# In[40]:


filtered_df = merged_df[merged_df['NumCausalSNPs'] != 0]
filtered_df


# In[41]:


X2 = filtered_df['NumCausalSNPs']  # Independent variable
y2 = filtered_df['Heritability_Score']    # Dependent variable

X2 = sm.add_constant(X2)

model2 = sm.OLS(y2, X2).fit()
print(model2.summary())


# In[42]:


# Plot it on a graph
plt.figure(figsize=(10, 6))
plt.scatter(filtered_df['NumCausalSNPs'], filtered_df['Heritability_Score'], alpha=0.5)
plt.title('Scatter Plot of Heritability Scores vs NumCausalSNPs')
plt.xlabel('NumCausalSNPs')
plt.ylabel('Heritability')
plt.grid(True)
plt.show()


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# # Appendix

# # Codes below belong to the failed approach due to potential incompatibility of R package and Windows

# #### Find the common individuals between Genoype data and Gene_expression data, for further filter of 22 genitype plink files

# In[11]:


idlist = gene_expression_filtered.columns.tolist()[1:]

with open('common_individuals.txt', 'w') as file:
    file.write('FID IID\n') 
    for iid in idlist:
        file.write(f'{iid} {iid}\n')


# In[ ]:





# ### Filter the 22 files
#  .\plink.exe --bfile 1000G.EUR.[chromosome_number] --keep common_individuals.txt --make-bed --out 1000G.EUR.[chromosome_number]_filtered
# 

# ### Filter Minor Allele Frequency < 0.05
# 1..22 | ForEach-Object {
#   .\plink.exe --bfile 1000G.EUR.${_}_filtered --maf 0.05 --make-bed --out 1000G.EUR.${_}_maf
# }

# ### Combine the genotype data first
# 
# .\plink.exe --bfile 1000G.EUR.1_maf --merge-list all_files_to_merge.txt --make-bed --out combined_genotype

# ### Then Generate the Combined GRM (without having to combine the 22 GRMs, since genotype data are already merged)
# 
# .\gcta-win-1.94.1.exe --bfile combined_genotype --autosome --make-grm --out combined_GRM

# ### Generate GRMs and Filter Individuals with Threshold of 0.025
# .\gcta-win-1.94.1.exe --grm combined_GRM --grm-cutoff 0.025 --make-grm --out filtered_combined_GRM
# 
# 

# #### After filtering with 0.025, 44 individuals are removed, so filter the common individuals again adn create the phenotype files for the remaining 3000 individuals

# In[11]:


with open('filtered_combined_GRM.grm.id', 'r') as file:
    grm_individuals = file.read().splitlines()

# Extract FID and IID from the GRM file
grm_ids = [line.split('\t')[0] for line in grm_individuals] 
grm_ids


# In[12]:


filtered_gene_expression = gene_expression_filtered[['TargetID'] + grm_ids]
filtered_gene_expression


# In[13]:


# create phenotype data
individual_ids = filtered_gene_expression.columns.tolist()[1:]


# In[14]:


for index, row in filtered_gene_expression.iterrows():
    gene_id = row['TargetID']

    phenotype_data = pd.DataFrame({
        'FID': individual_ids,
        'IID': individual_ids,
        'Phenotype': row[individual_ids].values
    })

    # Save the phenotype data without a header and index
    phenotype_data.to_csv(f'phenotypes/{gene_id}.txt', sep='\t', header=False, index=False)


# ## Filter the covariates file and then change it to the format that GCTA can manipulate

# In[15]:


co = pd.read_csv('tiffany_covariates.txt', sep='\t') 
co_filtered = co[['ID'] + [col for col in co.columns if col in individual_ids]]
co_filtered


# #### Change the Format of the dataframe

# In[16]:


melted_df = co_filtered.melt(id_vars=['ID'], var_name='IID', value_name='Covariate_Value')

pivoted_df = melted_df.pivot(index='IID', columns='ID', values='Covariate_Value').reset_index()

co_res = pivoted_df.rename_axis(None, axis=1)
co_res.insert(0, 'FID', co_res['IID'])
co_res


# #### Drop columns that are not needed

# In[17]:


co_res = co_res[["FID","IID","genoPC1","genoPC2","genoPC3","genoPC4","genoPC5","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10","sex"]]
co_res


# In[18]:


# Separate the discrete and quantitative covariates
discrete_co = co_res[['FID', 'IID', 'sex']]
quantitative_co = co_res.drop('sex', axis=1)
discrete_co


# In[19]:


quantitative_co


# In[67]:


# Save to files
discrete_co.to_csv('discrete_covariates.txt', sep='\t', index=False, header=False)
quantitative_co.to_csv('quantitative_covariates.txt', sep='\t', index=False, header=False)
co_res.to_csv('total_covariates.txt', sep='\t', index=False, header=False)


# In[ ]:





# #### See the proportion of variance explained by the covariates

# In[21]:


cc = co_res.copy()
cc = cc.drop(columns = "FID")
cc = cc.set_index("IID")
cc


# In[22]:


pp = pd.read_csv("phenotypes/ENSG00000000419.txt",sep = "\t",header = None)
pp = pp.drop(columns = 0)
pp = pp.set_index(1)
pp


# In[23]:


X = sm.add_constant(cc)  # independent variables (covariates)
y = pp[2]  # dependent variable (phenotype)

model = sm.OLS(y, X).fit()

r_squared = model.rsquared
print(model.summary())
print(f"Proportion of variance explained by covariates: {r_squared}")


# ### Try generating phenotypes data but regressing out covariates first

# In[61]:


# Define paths
covariates_file = 'total_covariates.txt'
phenotypes_folder = 'phenotypes'
adjusted_phenotypes_folder = 'phenouse'

covariates = pd.read_csv(covariates_file, sep='\t')
covariates.set_index('IID', inplace=True)
covariates.drop(columns = "FID", inplace=True)


# In[62]:


# Function to adjust phenotype
def adjust_phenotype(phenotype_file):
    phenotype_path = os.path.join(phenotypes_folder, phenotype_file)
    pp = pd.read_csv(phenotype_path, sep='\t', header=None, names=['FID', 'IID', 'Phenotype'])
    pp.set_index('IID', inplace=True)
    pp.drop(columns = "FID", inplace=True)

    # Merge phenotype and covariate data
    merged_data = pp.join(covariates, how='inner')

    # Prepare model data
    X = sm.add_constant(merged_data.drop(columns=['Phenotype']))  # Covariates with constant
    y = merged_data['Phenotype']  # Phenotype

    # Fit model and calculate residuals
    model = sm.OLS(y, X).fit()
    residuals = model.resid

    # Prepare adjusted phenotype DataFrame
    adjusted_phenotypes = pd.DataFrame({'FID': merged_data.index, 'IID': merged_data.index, 'AdjustedPhenotype': residuals})
    adjusted_phenotypes.reset_index(drop=True, inplace=True)

    # Save adjusted phenotype
    adjusted_file_path = os.path.join(adjusted_phenotypes_folder, phenotype_file)
    adjusted_phenotypes.to_csv(adjusted_file_path, sep='\t', index=False, header=False)

# Iterate over phenotype files and adjust
for phenotype_file in os.listdir(phenotypes_folder):
    adjust_phenotype(phenotype_file)


# In[ ]:





# In[ ]:





# ### Generate the Heritability scores
# Get-ChildItem -Path "phenouse\*.txt" | ForEach-Object {
#     $gene_id = $_.BaseName
#     .\gcta-win-1.94.1.exe --grm filtered_combined_GRM --pheno $_.FullName --reml --reml-maxit 1000 --out "results\$gene_id`_heritability_covar"
# }
# 
# 
# 
# #### Past Approach:
# Get-ChildItem -Path "phenotypes\*.txt" | ForEach-Object {
#     $gene_id = $_.BaseName
#     .\gcta-win-1.94.1.exe --grm filtered_combined_GRM --pheno $_.FullName --covar discrete_covariates.txt --qcovar quantitative_covariates.txt --reml --reml-maxit 1000 --out "results\$gene_id`_heritability_covar"
# }
# 

# In[ ]:





# ## Interpret the results

# In[63]:


results_lst = []

for file in os.listdir('results'):
    if file.endswith('_heritability_covar.hsq'):
        gene_id = file.split('_')[0]

        hsq_data = pd.read_csv(f'results/{file}', delim_whitespace=True)
        hsq_score = hsq_data.loc[hsq_data['Source'] == 'V(G)/Vp', 'Variance'].values[0]
        results_lst.append({'gene_id': gene_id, 'heritability': hsq_score})

results = pd.DataFrame(results_lst)


# In[64]:


# map the chromosome to each gene
chromosome_mapping = gene_expression_data.set_index('TargetID')['Chr'].to_dict()

results['Chromosome'] = results['gene_id'].map(chromosome_mapping)


# In[65]:


results


# #### Save the results to csv file

# In[58]:


results.to_csv('old_heritability_scores.csv', index=False)


# In[66]:


#5th
results["heritability"].value_counts()


# In[35]:


#4th
results["heritability"].value_counts()


# In[43]:


#3rd
results["heritability"].value_counts()


# In[57]:


#2nd
results["heritability"].value_counts()


# In[60]:


#1st
df = pd.read_csv("scores.csv")
df["heritability"].value_counts()


# After applying the covariates data, there is an increase in the number of 0.000001 but decrease in 0.999999.

# ## Linear Regression

# #### Merge the dataframes first

# In[61]:


merged_df = pd.merge(genes, results, left_on='GeneID', right_on='gene_id')
merged_df


# In[63]:


X = merged_df['CausalVariance']  # Independent variable
y = merged_df['heritability']    # Dependent variable

X = sm.add_constant(X)

model = sm.OLS(y, X).fit()


# In[64]:


print(model.summary())


# The R-squared value is 0.000, indicating that the model explains none of the variability of the response data around its mean. In other words, the causal variance does not seem to explain any variation in heritability in your dataset.
# 
# The Adjusted R-squared is also -0.000. This is a modified version of R² that has been adjusted for the number of predictors in the model. The negative value suggests that the model does not fit the data well and does not improve over a simpler model without predictors.
# 
# The F-statistic value is very low (0.04584), and the probability of the F-statistic (Prob (F-statistic)) is high (0.830). This suggests that the model is not statistically significant. The high p-value indicates that the observed R² is likely due to chance.
# 
# The constant (const) coefficient is 0.3737, which is the intercept of the regression line. It represents the expected mean value of the heritability when the causal variance is zero.
# The coefficient for CausalVariance is -0.0054. This means that for each unit increase in causal variance, the heritability decreases by 0.0054 units. However, this relationship is not statistically significant, as indicated by the p-value.
# 
# The p-value for the CausalVariance coefficient is 0.830, which is much higher than the typical alpha level of 0.05. This high p-value indicates that there is no significant relationship between causal variance and heritability.
# 
# The 95% confidence interval for the CausalVariance coefficient ranges from -0.055 to 0.044, which includes zero. This further suggests that the coefficient is not significantly different from zero.
# 
# The low R-squared value, coupled with the high p-value for the F-statistic and the coefficients, suggests that the linear model does not fit the data well. The causal variance does not appear to be a good predictor of heritability in this case.
# 
# In summary, ther linear regression analysis indicates that there is no significant relationship between causal variance and heritability in your dataset. The causal variance does not seem to be a good predictor of heritability for the genes you have analyzed.

# In[65]:


# Plot it on a graph
plt.figure(figsize=(10, 6))
plt.scatter(merged_df['CausalVariance'], merged_df['heritability'], alpha=0.5)
plt.title('Scatter Plot of Heritability vs Causal Variance')
plt.xlabel('Causal Variance')
plt.ylabel('Heritability')
plt.grid(True)
plt.show()


# In[ ]:





# In[ ]:





# In[43]:


### Pruning Linkage disequilibrium

#1..22 | ForEach-Object {
#  .\plink.exe --bfile 1000G.EUR.${_}_maf --indep-pairwise 50 5 0.2 --out 1000G.EUR.${_}_pruned
#}
### Extract SNPs
#
#1..22 | ForEach-Object {
#  .\plink.exe --bfile 1000G.EUR.${_}_maf --extract 1000G.EUR.${_}_pruned.prune.in --make-bed --out 1000G.EUR.${_}_final
#}


# In[ ]:


### Then we will merge the 22 GRMs into a combined GRMs 

# .\gcta-win-1.94.1.exe --mgrm grm_list.txt --make-grm --out GRM_combined


# In[ ]:


### 
#foreach ($i in 1..22) {
#    .\gcta-win-1.94.1.exe --bfile "1000G.EUR.${i}_maf" --autosome --make-grm --out "GRM_chr$i"
#    .\gcta-win-1.94.1.exe --grm "GRM_chr$i" --grm-cutoff 0.025 --make-grm --out "GRM_chr${i}_final"
#}


# In[ ]:


#Command lines to use the Rscript tp regress out covariates and generate new phenotype data fiels (for 1 phenotype)
Rscript compute_weights.R --bfile combined_genotype --pheno phenotypes/ENSG00000000419.txt --covar total_covariates.txt --out phenouse/ENSG00000000419 --tmp tmp
Rscript compute_weights.R --bfile combined_genotype --pheno phenotypes/ENSG00000000419.txt --covar total_covariates.txt --out phenouse/ENSG00000000419 --tmp tmp --save_hsq TRUE

# Command lines to use the Rscript tp regress out covariates and generate new phenotype data fiels (for all)
$RScriptPath = "Rscript"
$ComputeWeightsScriptPath = "compute_weights.R"
$BfilePath = "combined_genotype"
$CovarPath = "total_covariates.txt"
$OutputDir = "phenouse"
$TmpDir = "tmp" 

if (-not (Test-Path $TmpDir)) {
    New-Item -ItemType Directory -Force -Path $TmpDir
}

Get-ChildItem "phenotypes\*.txt" | ForEach-Object {
    $PhenoFilePath = $_.FullName
    $OutputPrefix = Join-Path $OutputDir $_.BaseName

    $Cmd = "$RScriptPath $ComputeWeightsScriptPath --bfile $BfilePath --pheno $PhenoFilePath --covar $CovarPath --out $OutputPrefix --tmp $TmpDir"

    Invoke-Expression $Cmd

    # Optional: Output the command for logging
    #Write-Host "Processed: $PhenoFilePath"
}

