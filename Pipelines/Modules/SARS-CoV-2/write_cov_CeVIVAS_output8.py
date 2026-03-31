#!/usr/bin/env python
# coding: utf-8

# Importando bibliotecas
import pandas as pd
import sys
import os
import re

# Entrada de arquivos
statistics_file = sys.argv[1]  # All_Statistics.csv
assembly_path = sys.argv[2]  # Caminho da montagem
raw_data_path = sys.argv[3]  # Caminho dos dados brutos
folder_name = sys.argv[4]  # Nome da pasta

# Tente ler o arquivo como CSV ou TSV
try:
    statistics_df = pd.read_csv(statistics_file, sep=',')  # ou sep=',' para CSV
    print(f"Arquivo {statistics_file} lido com sucesso como CSV/TSV.")
except Exception as e:
    print(f"Erro ao ler o arquivo {statistics_file} como CSV/TSV: {e}")
    sys.exit(1)

# Lendo os dados brutos
try:
    path = pd.read_csv(raw_data_path, sep='\t', header=None)
    if path.empty:
        print(f"O arquivo {raw_data_path} está vazio.")
        path = pd.DataFrame(columns=[0])  # Cria um DataFrame vazio com cabeçalho
except Exception as e:
    print(f"Erro ao ler o arquivo CSV: {e}")
    path = pd.DataFrame(columns=[0])  # Cria um DataFrame vazio caso o arquivo não exista

# Ajustando ID_metadata para fusão
if 'Genome' in statistics_df.columns:
    statistics_df['ID_metadata'] = statistics_df['Genome'].str.replace('Genoma_Indels__', '', regex=True)
    statistics_df['ID_metadata'] = statistics_df['ID_metadata'].str.replace(r'_S\d+.+', '', regex=True)
    
    # Verificando se ID_metadata foi criado corretamente
    print("Pré-merge - ID_metadata:")
    print(statistics_df[['Genome', 'ID_metadata']].head())
else:
    print("Erro: A coluna 'Genome' não existe no DataFrame!")
    print("Colunas disponíveis:", statistics_df.columns)
    sys.exit(1)

# Confirmar se 'ID_metadata' foi criada corretamente
if 'ID_metadata' not in statistics_df.columns:
    print("Erro crítico: 'ID_metadata' não foi criada corretamente!")
    sys.exit(1)

# Depurar também os dados do arquivo read_path.tsv
print("Pré-merge - Dados de read_path.tsv:")
print(path.head())

# Ajustar a coluna 0 do path:
# 1. Remover espaços em branco no final
# 2. Remover sufixos como _S2, _S5, etc.
path[0] = path[0].str.strip()  # Remove espaços em branco
path[0] = path[0].str.replace(r'_S\d+$', '', regex=True)  # Remove sufixos _S<número>

# Verificar valores únicos nas colunas de junção
print("Valores únicos em statistics_df['ID_metadata']:")
print(statistics_df['ID_metadata'].unique())

print("Valores únicos em path[0]:")
print(path[0].unique())

# Verificar se há valores em comum
common_values = set(statistics_df['ID_metadata']).intersection(set(path[0]))
print("Valores em comum entre statistics_df['ID_metadata'] e path[0]:")
print(common_values)

if not common_values:
    print("Erro: Não há valores em comum entre statistics_df['ID_metadata'] e path[0]. O merge não pode ser realizado.")
    sys.exit(1)

# Fazer o merge
try:
    cevivas_output = pd.merge(statistics_df, path, left_on='ID_metadata', right_on=0)
    print("Merge realizado com sucesso!")
    print("Colunas disponíveis após o merge:")
    print(cevivas_output.columns)
except Exception as e:
    print(f"Erro ao mesclar os DataFrames: {e}")
    sys.exit(1)

# Removendo a coluna duplicada do merge
cevivas_output = cevivas_output.drop(columns=[0])

# Renomeando colunas
cevivas_output = cevivas_output.rename(columns={
    'Genome': 'Genoma',
    'clade': 'Clade',
    'SNPs_count': 'SNPs',  # Renomeia SNPs_count para SNPs
    1: 'Caminho_dados_brutos_vital'
})

# Adicionando colunas vazias para preenchimento manual
cevivas_output['CEVIVAS_ID'] = ''
cevivas_output['Repeat'] = ''
cevivas_output['Estudo'] = ''
cevivas_output['Analysis'] = ''

# Adicionando caminho da montagem e data de montagem
PATHWAY = assembly_path
assembly_basename = os.path.basename(os.path.normpath(PATHWAY))
date_match = re.search(r'(\d{2}-\d{2}-\d{4})', assembly_basename)
if date_match:
    Data_Montagem = date_match.group(1)
else:
    pathway_tokens = [token for token in assembly_basename.split('_') if token]
    Data_Montagem = pathway_tokens[1] if len(pathway_tokens) > 1 else assembly_basename

cevivas_output['Caminho_montagem_vital'] = assembly_path
cevivas_output['Data_Montagem'] = Data_Montagem

# Reordenando colunas conforme necessário
final_columns = [
    'ID_metadata', 'CEVIVAS_ID', 'Genoma', 'Pangolin_lineage', 'Clade', 'Coverage', 'Passed_QC', 'Repeat',
    'Estudo', 'Analysis', 'SNPs', 'Number_of_Ns', 'N_Reads', 'Reads_mapped', 'Percent_mapped', 'Mean_depth',
    'Median_depth', 'Npos_Depth>=10', 'Npos_Depth>=25', 'scorpio_call', 'S_Mean_depth', 'ORF1ab_mean_depth',
    'N_mean_depth', 'PCR1_mean_depth', 'PCR2_mean_depth', 'PCR1_median_depth',
    'PCR2_median_depth', 'Intrahost', 'qc.overallStatus', 'deletions', 'insertions', 'aaSubstitutions',
    'aaDeletions', 'totalFrameShifts', 'frameShifts', 'substitutions', 'Caminho_montagem_vital',
    'Caminho_dados_brutos_vital', 'Data_Montagem'
]

# Verificar se todas as colunas em final_columns existem no DataFrame
missing_columns = [col for col in final_columns if col not in cevivas_output.columns]
if missing_columns:
    print(f"Aviso: As seguintes colunas não existem no DataFrame e serão ignoradas: {missing_columns}")
    final_columns = [col for col in final_columns if col in cevivas_output.columns]

# Selecionar apenas as colunas presentes no DataFrame
cevivas_output = cevivas_output[final_columns]

# Verifique os dados antes de salvar o arquivo CSV
print("Dados finais antes de salvar:")
print(cevivas_output.head())  # Isso vai imprimir as primeiras linhas do dataframe

# Verifique se chegou até o ponto de salvar o arquivo
output_file_name = f'CeVIVAS_CoV_{folder_name}.csv'
print(f"Salvando o arquivo {output_file_name}")

# Salvando o arquivo consolidado
try:
    cevivas_output.to_csv(output_file_name, index=False)
    print(f"Arquivo {output_file_name} salvo com sucesso!")
except Exception as e:
    print(f"Erro ao salvar o arquivo {output_file_name}: {e}")
    sys.exit(1)
