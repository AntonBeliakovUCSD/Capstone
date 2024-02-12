import pandas as pd 
import numpy as np
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler
import seaborn as sns
import pyBigWig
from statsmodels.stats.multitest import multipletests
import os 
import sys 
import argparse


#clean twas results
def clean_twas_results(filepath):
    """
    filepath: this is the filepath for raw TWAS results
    
    output: The output is the cleaned up TWAS with appropriate column types and 
    it is corrected for significance testing. This can now be used for visualization
    and further significance testing. 
    """
    twas_results = pd.read_csv(filepath, sep='\t') 
    def try_convert_to_float(x):
        try:
            return float(x)
        except ValueError:  
            return np.nan
        
    float_cols = ['HSQ', 'BEST.GWAS.Z', 'EQTL.R2', 'EQTL.Z', 'EQTL.GWAS.Z', 'TWAS.Z', 'TWAS.P']
    for col in float_cols:
        twas_results[col] = twas_results[col].apply(lambda x: try_convert_to_float(x))
    twas_results['CHR'] = twas_results['CHR'].astype(str)
    # correct for multiple testin
    pvals_corrected = multipletests(twas_results['TWAS.P'].dropna(), method='fdr_bh')[1]
    twas_results.loc[twas_results['TWAS.P'].notnull(), 'TWAS.P.adj'] = pvals_corrected
    return twas_results

#create manhatten plot 
def plot_manhattan(df):
    """
    df: This is a cleaned TWAS results dataframe

    output: A Manhatten Plot is saved in the current directory. The signicant
    gene ids are labeled. 
    """
    df['-log10(TWAS.P)'] = -np.log10(df['TWAS.P'].replace(0, np.nan))
    df.sort_values(['CHR', 'P0'], inplace=True)
    df['CHR_num'] = df['CHR'].map(dict(zip(df['CHR'].unique(), range(1, len(df['CHR'].unique()) + 1))))
    num_chromosomes = len(df['CHR'].unique())
    colors = ['red', 'blue'] * (num_chromosomes // 2 + 1)  

    plt.figure(figsize=(12, 6))
    significance_threshold = -np.log10(0.05 / len(df))
    
    for chrom in df['CHR'].unique():
        chrom_df = df[df['CHR'] == chrom]
        plt.scatter(chrom_df['P0'], chrom_df['-log10(TWAS.P)'], label=chrom, s=10)
        
        significant_snps = chrom_df[chrom_df['-log10(TWAS.P)'] > significance_threshold]
        for _, row in significant_snps.iterrows():
            plt.text(row['P0'], row['-log10(TWAS.P)'], row['ID'], fontsize=8, rotation=45)

    plt.xlabel('Genomic Position')
    plt.ylabel('-log10(TWAS.P)')
    plt.title('Manhattan Plot of TWAS Results')
    plt.axhline(y=significance_threshold, color='grey', linestyle='--')  # Draw significance threshold line
    plt.legend(title='Chromosome', bbox_to_anchor=(1.05, 1), loc='upper left', ncol=2)
    plt.tight_layout()
    plt.savefig("manhatten_plot.png")
    plt.show()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="TWAS data processing")
    parser.add_argument('twas_file', help='Path to the TWAS results file')
    args = parser.parse_args()

    cleaned_data = clean_twas_results(args.twas_file)

    plot_manhattan(cleaned_data)