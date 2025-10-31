#!/usr/bin/env python
# coding: utf-8

# Importing libraries
import pandas as pd
import sys
import re

# ------------------------------------------------------------------------
# 1) Ler argumentos de linha de comando
# ------------------------------------------------------------------------
statistics_file = sys.argv[1]    # Arquivo *_complete.Statistics
assembly_path   = sys.argv[2]    # Caminho absoluto para a pasta de montagem
raw_data_path   = sys.argv[3]    # TSV com colunas: <sample_id> <caminho_R1> ...
folder_name     = sys.argv[4]    # Nome/prefixo (ex: data ou nome do projeto)

# ------------------------------------------------------------------------
# 2) Ler arquivo de estatísticas completas (tab-separated)
#    Este arquivo contém todas as informações de coverage, profundidade etc.
# ------------------------------------------------------------------------
statistics_df = pd.read_csv(statistics_file, sep='\t')

# ------------------------------------------------------------------------
# 3) Selecionar as colunas necessárias para o relatório final
#    Ajuste se o seu arquivo tiver nomes de colunas diferentes.
# ------------------------------------------------------------------------
stats_columns = [
    # Básico
    'Genome',
    'N_Reads',
    'Reads_assembled',
    'Percent_assembled',
    'Segments_Assembled',
    'Type',
    'H_genotype',
    'N_genotype',
    'Subtype',
    'Clade',

    # Profundidades e cobertura por segmento
    'segment_1_Mean_depth',
    'segment_1_Median_depth',
    'segment_1_Npos_Depth>=10',
    'segment_1_Npos_Depth>=25',
    'segment_1_Coverage',
    'segment_1_Number_of_Ns',
    'segment_1_SNPs',

    'segment_2_Mean_depth',
    'segment_2_Median_depth',
    'segment_2_Npos_Depth>=10',
    'segment_2_Npos_Depth>=25',
    'segment_2_Coverage',
    'segment_2_Number_of_Ns',
    'segment_2_SNPs',

    'segment_3_Mean_depth',
    'segment_3_Median_depth',
    'segment_3_Npos_Depth>=10',
    'segment_3_Npos_Depth>=25',
    'segment_3_Coverage',
    'segment_3_Number_of_Ns',
    'segment_3_SNPs',

    'segment_4_Mean_depth',
    'segment_4_Median_depth',
    'segment_4_Npos_Depth>=10',
    'segment_4_Npos_Depth>=25',
    'segment_4_Coverage',
    'segment_4_Number_of_Ns',
    'segment_4_SNPs',

    'segment_5_Mean_depth',
    'segment_5_Median_depth',
    'segment_5_Npos_Depth>=10',
    'segment_5_Npos_Depth>=25',
    'segment_5_Coverage',
    'segment_5_Number_of_Ns',
    'segment_5_SNPs',

    'segment_6_Mean_depth',
    'segment_6_Median_depth',
    'segment_6_Npos_Depth>=10',
    'segment_6_Npos_Depth>=25',
    'segment_6_Coverage',
    'segment_6_Number_of_Ns',
    'segment_6_SNPs',

    'segment_7_Mean_depth',
    'segment_7_Median_depth',
    'segment_7_Npos_Depth>=10',
    'segment_7_Npos_Depth>=25',
    'segment_7_Coverage',
    'segment_7_Number_of_Ns',
    'segment_7_SNPs',

    'segment_8_Mean_depth',
    'segment_8_Median_depth',
    'segment_8_Npos_Depth>=10',
    'segment_8_Npos_Depth>=25',
    'segment_8_Coverage',
    'segment_8_Number_of_Ns',
    'segment_8_SNPs',
]

# Garante que as colunas existam no DataFrame (evita erro se faltar algo).
for col in stats_columns:
    if col not in statistics_df.columns:
        statistics_df[col] = None  # Cria coluna vazia se não existir

statistics_df_useful = statistics_df[stats_columns].copy()

# ------------------------------------------------------------------------
# 4) Criar ID_metadata a partir da coluna "Genome"
#    Exemplo: "Genoma_FLU__<ID>_S1..." => "<ID>"
# ------------------------------------------------------------------------
statistics_df_useful['ID_metadata'] = (
    statistics_df_useful['Genome']
    .str.replace('Genoma_FLU__', '', regex=False)
    .str.replace(r'_S\d+.+', '', regex=True)
)

# ------------------------------------------------------------------------
# 5) Determinar quais segmentos foram realmente montados (coluna "Segments")
#    Baseado na Coverage != "not_detected"/"could_not_be_assembled"
# ------------------------------------------------------------------------
segments_series = {}
for index, row in statistics_df_useful.iterrows():
    seg_list = []
    for seg_col in [
        'segment_1_Coverage','segment_2_Coverage','segment_3_Coverage',
        'segment_4_Coverage','segment_5_Coverage','segment_6_Coverage',
        'segment_7_Coverage','segment_8_Coverage'
    ]:
        val = str(row[seg_col])
        if ('not_detected' not in val) and ('could_not_be_assembled' not in val):
            # Ex: 'segment_1_Coverage' => '1'
            seg_num = seg_col.split('_')[1]
            seg_list.append(seg_num)
    # Junta com pipe, ex: "|1|2|3"
    segments_series[index] = '|' + '|'.join(seg_list) if seg_list else ''

statistics_df_useful['Segments'] = pd.Series(segments_series)

# ----------------------------------------------------------------
# 6) Definir Passed_QC com base em segment_4_Coverage
#    - "A" se coverage >= 85
#    - "E" se 50 <= coverage < 85
#    - "R" se not_detected, could_not_be_assembled ou < 50
# ----------------------------------------------------------------
Passed_QC = {}
for index, row in statistics_df_useful.iterrows():
    seg4_value = str(row['segment_4_Coverage'])
    try:
        if seg4_value in ['not_detected', 'could_not_be_assembled']:
            Passed_QC[index] = 'R'
        else:
            val = float(seg4_value)
            if val >= 85:
                Passed_QC[index] = 'A'
            elif 50 <= val < 85:
                Passed_QC[index] = 'E'
            else:
                Passed_QC[index] = 'R'
    except ValueError:
        # Caso a conversão float falhe (por ex. valor vazio), marca como 'R'
        Passed_QC[index] = 'R'

statistics_df_useful['Passed_QC'] = pd.Series(Passed_QC)

# ----------------------------------------------------------------
# 7) Ler arquivo de caminhos de dados brutos (raw_data_path)
#    e fazer merge para obter 'Caminho_dados_brutos_vital'
# ----------------------------------------------------------------
path_df = pd.read_csv(raw_data_path, sep='\t', header=None)
# Ajusta a primeira coluna para remover sufixos, ex: "_S1_L001" etc.
path_df[0] = path_df[0].str.replace(r'_S\d+.+', '', regex=True)

# Faz merge com base em 'ID_metadata'
joined = pd.merge(statistics_df_useful, path_df, left_on='ID_metadata', right_on=0, how='left')

# ----------------------------------------------------------------
# 8) Ajustes finais de colunas
#    - Renomear "Genome" -> "Genoma"
#    - Renomear a coluna 1 (path_df) para "Caminho_dados_brutos_vital"
# ----------------------------------------------------------------
joined.rename(columns={
    'Genome': 'Genoma',
    1: 'Caminho_dados_brutos_vital'
}, inplace=True)

# Cria coluna 'folder' para identificar subpasta da montagem
joined['folder'] = joined['Genoma'].str.replace('Genoma_', '', regex=False)
joined['folder'] = joined['folder'].str.replace('.fasta', '', regex=False)

# Extrai a data do sequenciamento a partir do assembly_path
# Exemplo: se assembly_path = "28-02-2025"
# Ajuste se seu padrão de pasta for diferente
parts = assembly_path.rstrip().split('/')
last_part = parts[-1]

# Tenta capturar algo no formato dd-mm-aaaa
match = re.search(r'(\d{2}-\d{2}-\d{4})', last_part)
if match:
    # Se encontrar a data no padrão "28-02-2024", pega ela
    Data_Montagem = match.group(1)
else:
    # Caso não encontre, use o nome inteiro da pasta (ou outro critério)
    Data_Montagem = last_part

joined['Data_Montagem'] = Data_Montagem

# Define caminhos
joined['Caminho_montagem_vital'] = joined.apply(
    lambda row: f"{assembly_path}/{row['folder']}/segments/", axis=1
)

# Cria colunas vazias
joined['Notes'] = ''
joined['CEVIVAS_ID'] = ''
joined['Repeat'] = ''
joined['Estudo'] = ''
joined['Analysis'] = ''

# ----------------------------------------------------------------
# 9) Reordenar colunas para a forma final solicitada
# ----------------------------------------------------------------
final_cols = [
    'ID_metadata',
    'CEVIVAS_ID',
    'Genoma',
    'Type',
    'Subtype',
    'Clade',
    'Segments_Assembled',
    'Segments',
    'Passed_QC',
    'Repeat',
    'Estudo',
    'Analysis',

    'segment_1_Coverage',
    'segment_2_Coverage',
    'segment_3_Coverage',
    'segment_4_Coverage',
    'segment_5_Coverage',
    'segment_6_Coverage',
    'segment_7_Coverage',
    'segment_8_Coverage',

    'N_Reads',
    'Reads_assembled',
    'Percent_assembled',

    'segment_1_Mean_depth',
    'segment_1_Median_depth',
    'segment_1_Npos_Depth>=10',
    'segment_1_Npos_Depth>=25',
    'segment_1_Number_of_Ns',
    'segment_1_SNPs',

    'segment_2_Mean_depth',
    'segment_2_Median_depth',
    'segment_2_Npos_Depth>=10',
    'segment_2_Npos_Depth>=25',
    'segment_2_Number_of_Ns',
    'segment_2_SNPs',

    'segment_3_Mean_depth',
    'segment_3_Median_depth',
    'segment_3_Npos_Depth>=10',
    'segment_3_Npos_Depth>=25',
    'segment_3_Number_of_Ns',
    'segment_3_SNPs',

    'segment_4_Mean_depth',
    'segment_4_Median_depth',
    'segment_4_Npos_Depth>=10',
    'segment_4_Npos_Depth>=25',
    'segment_4_Number_of_Ns',
    'segment_4_SNPs',

    'segment_5_Mean_depth',
    'segment_5_Median_depth',
    'segment_5_Npos_Depth>=10',
    'segment_5_Npos_Depth>=25',
    'segment_5_Number_of_Ns',
    'segment_5_SNPs',

    'segment_6_Mean_depth',
    'segment_6_Median_depth',
    'segment_6_Npos_Depth>=10',
    'segment_6_Npos_Depth>=25',
    'segment_6_Number_of_Ns',
    'segment_6_SNPs',

    'segment_7_Mean_depth',
    'segment_7_Median_depth',
    'segment_7_Npos_Depth>=10',
    'segment_7_Npos_Depth>=25',
    'segment_7_Number_of_Ns',
    'segment_7_SNPs',

    'segment_8_Mean_depth',
    'segment_8_Median_depth',
    'segment_8_Npos_Depth>=10',
    'segment_8_Npos_Depth>=25',
    'segment_8_Number_of_Ns',
    'segment_8_SNPs',

    'Caminho_dados_brutos_vital',
    'Caminho_montagem_vital',
    'Data_Montagem',
    'Notes'
]

# Garante que todas as colunas existam (se faltarem, preenche com None)
for col in final_cols:
    if col not in joined.columns:
        joined[col] = None

# Ordena de fato
final_df = joined[final_cols]

# ----------------------------------------------------------------
# 10) Salvar a planilha final em CSV
# ----------------------------------------------------------------
output_file_name = f'CeVIVAS_FLU_{folder_name}.csv'
final_df.to_csv(output_file_name, index=False)

print(f"[OK] Planilha final gerada: {output_file_name}")
