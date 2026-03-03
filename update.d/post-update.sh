#!/usr/bin/env bash
set -euo pipefail

# Determina o caminho raiz do repo a partir deste script
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

# Recria os links simbólicos dos pipelines principais
sudo ln -sf "$ROOT/Pipelines/Modules/SARS-CoV-2/Exec_assembly_pipeline_Illumina_v8_bowtie2_ref_iVar_CeVIVAS.sh" /usr/local/bin/VIPER_CoV.sh
sudo ln -sf "$ROOT/Pipelines/Modules/DENV/Exec_assembly_pipeline_Illumina_v5_bwa_mem_ref_iVar.sh" /usr/local/bin/VIPER_DENV.sh
sudo ln -sf "$ROOT/Pipelines/Modules/Influenza/Exec_assembly_pipeline_Illumina_v3_Vapor_SPAdes.sh" /usr/local/bin/VIPER_FLU.sh

# Garante permissao de execucao nos alvos
sudo chmod +x /usr/local/bin/VIPER_CoV.sh /usr/local/bin/VIPER_DENV.sh /usr/local/bin/VIPER_FLU.sh

echo "Links simbolicos atualizados em /usr/local/bin (VIPER_CoV.sh, VIPER_DENV.sh, VIPER_FLU.sh)"
