#!/usr/bin/env python
# coding: utf-8

# Importing libraries
import pandas as pd
import sys
import numpy as np

# Getting files from sys.argv
stats_file = sys.argv[1]
nextclade_file = sys.argv[2]

# Reading files 
stats = pd.read_csv(stats_file, sep='\t', dtype={'Genome': str})
nextclade = pd.read_csv(nextclade_file, sep='\t', dtype={'seqName': str})

# Processing
nextclade = nextclade.rename(columns = {'seqName':'Genome'}) # Renaming column in nextclade output

# Ensure both Genome columns are strings
stats['Genome'] = stats['Genome'].astype(str)
nextclade['Genome'] = nextclade['Genome'].astype(str)

nextclade_useful = nextclade[['Genome','clade','qc.overallStatus', 'deletions','insertions','aaSubstitutions', 'aaDeletions', 'totalFrameShifts', 'frameShifts', 'substitutions']] # Selecting only necessary information

df = pd.merge(stats, nextclade_useful, on=['Genome'], how="outer")
df['Passed_QC'] = np.where(df['Npos_Depth>=10'] < 20000, 'R', np.where(df['Coverage'] < 85, 'R', np.where(df['totalFrameShifts'] > 0, 'M', np.where(df['SNPs_count'] >= 10, 'R', 'A'))))

output_fileName = stats_file.replace('_sorted.tsv', '.csv')
df.to_csv(output_fileName, index=False) # Writing results to file