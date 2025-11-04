#!/bin/bash

# Obtém o caminho absoluto do diretório atual
PIPELINE_PATH=$(realpath .)

# Adiciona o diretório ao PATH
echo "export PATH=\$PATH:$PIPELINE_PATH" >> ~/.bashrc
echo "export PIPELINE=$PIPELINE_PATH" >> ~/.bashrc

# Carrega as alterações no shell atual
source ~/.bashrc

echo "Current directory added to PATH as PIPELINE"

sudo ln -s $PIPELINE_PATH/SARS-CoV-2/Exec_assembly_pipeline_Illumina_v8_bowtie2_ref_iVar_CeVIVAS.sh /usr/local/bin/VIPER_CoV.sh
sudo ln -s $PIPELINE_PATH/DENV/Exec_assembly_pipeline_Illumina_v5_bwa_mem_ref_iVar.sh /usr/local/bin/VIPER_DENV.sh
sudo ln -s $PIPELINE_PATH/Influenza/Exec_assembly_pipeline_Illumina_v3_Vapor_SPAdes.sh /usr/local/bin/VIPER_Influenza.sh

chmod +x $PIPELINE_PATH/SARS-CoV-2/Exec_assembly_pipeline_Illumina_v8_bowtie2_ref_iVar_CeVIVAS.sh
chmod +x $PIPELINE_PATH/DENV/Exec_assembly_pipeline_Illumina_v5_bwa_mem_ref_iVar.sh
chmod +x $PIPELINE_PATH/Influenza/Exec_assembly_pipeline_Illumina_v3_Vapor_SPAdes.sh

echo "All requirements are done!"