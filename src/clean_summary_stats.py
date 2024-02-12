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
import ast 

# Read in GWAS Summary Statistics and Prepare them for TWAS
def clean_summary_stats(filepath, key_col_names, trait):
    """
    filepath: filepath to summary statistics file 
    key_col_names: a list of column names that represent 'SNP', 'A1', 'A2', 'Z'. 
    This also assumes that the columns 'beta' and 'standard_error' are present. 
    
    output: summary statistics file prepared for TWAS
    """
    key_col_names = ast.literal_eval( key_col_names)
    sum_stats = pd.read_csv(filepath, sep='\t')
    sum_stats['Z'] = sum_stats['beta'] / sum_stats['standard_error']
    sum_stats = sum_stats[key_col_names]
    sum_stats.columns = ['SNP', 'A1', 'A2', 'Z']
    sum_stats.to_csv(f'{trait}_summary_stats.txt', sep='\t', index=False)
    return sum_stats

if __name__ == "__main__":
    input_file_path = sys.argv[1]
    cols = sys.argv[2]
    trait = sys.argv[3]
    clean_summary_stats(input_file_path, cols, trait)