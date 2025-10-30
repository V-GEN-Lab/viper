#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import sys
from openpyxl import Workbook

def main():
    # Parâmetros de entrada
    statistics_file = sys.argv[1]  # Arquivo .Statistics (ex: All_Statistics__teste_06-02-2025_SE0000.tsv)
    assembly_path = sys.argv[2]    # Caminho da montagem
    raw_data_path = sys.argv[3]    # Caminho dos dados brutos
    folder_name = sys.argv[4]      # Nome da pasta para output

    # Ler o arquivo de estatísticas com ponto como decimal (ajuste se necessário)
    stats_full = pd.read_csv(
        statistics_file,
        sep='\t',
        decimal='.',     # Se o arquivo usa ponto como separador decimal
        na_values=['NA','']
    )

    # Converter colunas numéricas, se necessário
    for col in ['Coverage', 'E_Coverage', 'SNPs']:
        if col in stats_full.columns:
            stats_full[col] = pd.to_numeric(stats_full[col], errors='coerce')

    # Usar uma cópia para outras operações (se necessário)
    statistics_df = stats_full.copy()

    # Ler o arquivo com os caminhos dos dados brutos
    path = pd.read_csv(raw_data_path, sep='\t', header=None)
    path[0] = path[0].str.replace(r'_S\d+.+', '', regex=True)

    # Extrair data de montagem do assembly_path (ajuste conforme sua convenção)
    PATHWAY = assembly_path
    Data_Montagem = PATHWAY.rstrip().split('/')[-1].split('_')[1]

    # Criar DataFrame inicial com as colunas de interesse
    statistics_df_useful = statistics_df[['Genome', 'Serotype', 'clade', 'SNPs', 'Coverage', 'E_Coverage']].copy()

    # Ajustar ID_metadata
    statistics_df_useful.loc[:, 'ID_metadata'] = (
        statistics_df_useful['Genome']
        .str.replace('Genoma_DENV__', '', regex=False)
        .str.replace(r'_S\d+.+', '', regex=True)
    )

    # Dicionário para mapear clade para genotype
    clade_to_genotype = {
        '1I': 'DENV1/I', '1II': 'DENV1/II', '1III': 'DENV1/III', '1IV': 'DENV1/IV', '1V': 'DENV1/V', '1VI': 'DENV1/VI',
        '2I': 'DENV2/I', '2II': 'DENV2/II', '2III': 'DENV2/III', '2IV': 'DENV2/IV', '2V': 'DENV2/V', '2VI': 'DENV2/VI',
        '3I': 'DENV3/I', '3II': 'DENV3/II', '3III': 'DENV3/III', '3IV': 'DENV3/IV', '3V': 'DENV3/V',
        '4I': 'DENV4/I', '4II': 'DENV4/II', '4III': 'DENV4/III', '4IV': 'DENV4/IV'
    }
    statistics_df_useful.loc[:, 'Genotype'] = (
        statistics_df_useful['clade']
        .str.split('_').str[0]
        .map(clade_to_genotype)
    )

    # Função para classificar QC
    def classify_qc(cov, ecov):
        if pd.isna(cov):
            return 'R'
        if cov >= 85.0:
            return 'A'
        elif not pd.isna(ecov) and ecov >= 55.0:
            return 'E'
        else:
            return 'R'

    statistics_df_useful['Passed_QC'] = [
        classify_qc(row['Coverage'], row['E_Coverage'])
        for _, row in statistics_df_useful.iterrows()
    ]

    # Merge com o arquivo de caminhos
    statistics_df_useful2 = pd.merge(
        statistics_df_useful,
        path,
        left_on='ID_metadata',
        right_on=0,
        how='left'
    )
    DF_copy = statistics_df_useful2.copy()

    # Adicionar colunas extras
    DF_copy['Caminho_montagem_vital'] = PATHWAY
    DF_copy['Data_Montagem'] = Data_Montagem

    # Renomear colunas e criar CEVIVAS_ID
    cevivas_output = DF_copy.rename(columns={
        'Genome': 'Genoma',
        'clade': 'Lineage',
        1: 'Caminho_dados_brutos_vital'
    })
    cevivas_output['CEVIVAS_ID'] = ''

    # Merge final com stats_full usando suffixes para evitar conflito
    merged = pd.merge(
        cevivas_output,
        stats_full,
        left_on='Genoma',
        right_on='Genome',
        how='left',
        suffixes=('_left','_right')
    )

    # Escolher explicitamente as colunas desejadas
    merged['Coverage'] = merged['Coverage_right']
    merged['E_Coverage'] = merged['E_Coverage_right']
    merged['Serotype'] = merged['Serotype_right']
    merged['SNPs'] = merged['SNPs_right']
    
    # Para a coluna Genotype: se "Genotype_left" existir, use-a; senão, mantenha "Genotype" (do lado esquerdo)
    if 'Genotype_left' in merged.columns:
        merged['Genotype'] = merged['Genotype_left']
    else:
        merged['Genotype'] = merged['Genotype']

    # Definir colunas finais desejadas
    final_columns = [
        'ID_metadata', 'CEVIVAS_ID', 'Genoma', 'Serotype', 'Genotype', 'Lineage',
        'Coverage', 'E_Coverage', 'Passed_QC',
        'Repeat', 'Estudo', 'Analysis',  # ficarão vazias
        'SNPs', 'Number_of_Ns', 'N_Reads', 'Reads_mapped', 'Percent_mapped',
        'Mean_depth', 'Median_depth', 'Npos_Depth>=10', 'Npos_Depth>=25',
        'E_Depth', 'NS1_Coverage', 'NS1_Depth', 'NS3_Coverage', 'NS3_Depth',
        'NS5_Coverage', 'NS5_Depth',
        'Caminho_montagem_vital', 'Caminho_dados_brutos_vital',
        'Data_Montagem', 'Notes'
    ]

    # Para as colunas que devem ficar vazias, atribuir string vazia se não existirem
    for col in final_columns:
        if col not in merged.columns:
            if col in ['Repeat','Estudo','Analysis','Notes']:
                merged[col] = ''
            else:
                merged[col] = 'NA'

    final_df = merged[final_columns]

    # Salvar o resultado em Excel
    output_file = f'CeVIVAS_DENV_{folder_name}.xlsx'
    with pd.ExcelWriter(output_file, engine='openpyxl') as writer:
        final_df.to_excel(writer, index=False, sheet_name='Consolidado')

if __name__ == "__main__":
    main()
