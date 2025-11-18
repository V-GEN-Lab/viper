#!/usr/bin/env python
# coding: utf-8

# Importing libraries
import pandas as pd
import sys
import numpy as np

# Getting files from sys.argv
stats_file = sys.argv[1]
nextclade_file = sys.argv[2]
lcs_file = sys.argv[3]

# Reading files 
stats = pd.read_csv(stats_file, sep='\t')
nextclade = pd.read_csv(nextclade_file, sep='\t')
lcs = pd.read_csv(lcs_file, sep='\t')

# Processing

lcs_grouped = lcs.sort_values(['sample','proportion'],ascending=False).groupby('sample').head(2)
first = lcs_grouped.drop_duplicates(subset=['sample'], keep='first')
last = lcs_grouped.drop_duplicates(subset=['sample'], keep='last')
lcs_joined = pd.merge(first, last, on=['sample'], how='left')
lcs_joined = lcs_joined[['sample', 'variant_group_x', 'proportion_x', 'variant_group_y', 'proportion_y']]

nextclade = nextclade.rename(columns = {'seqName':'Genome'}) # Renaming column in nextclade output
#nextclade['Genome'] = ["Genoma_" + str(x) + ".fasta" for x in nextclade["Genome"]] # Adding "Genoma_" prefix and ".fasta" sufix 
nextclade_useful = nextclade[['Genome','clade','qc.overallStatus', 'deletions','insertions','aaSubstitutions', 'aaDeletions', 'substitutions']] # Selecting only necessary information
#df = pd.merge(stats, nextclade_useful, on=['Genome']) # Merging data
df = pd.merge(stats, nextclade_useful,on=['Genome'], how="outer")
#df = df[df['clade'].notna()] # Dropping data without clade information
#df['Passed QC?'] = np.where(df['Npos_Depth>=10'] < 20000, 'N', 'S')
df['Passed QC?'] = np.where(df['Npos_Depth>=10'] < 20000, 'N', np.where(df['Coverage'] < 85, 'N', np.where(df['SNPs_count'] >= 5, 'R', 'S')))

df['sample']= df['Genome'].str.split('__', expand=True)[1]
df = pd.merge(df, lcs_joined,on=['sample'], how="outer")
df = df.drop(columns=['sample'])
df['Recombination suspect?'] = np.where(df['proportion_y'] > 0.1, 'S', 'N')

output_fileName = stats_file[0:-11] + ".xlsx"
df.to_excel(output_fileName, index=False) # Writing results to file 
