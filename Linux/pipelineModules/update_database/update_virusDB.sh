#!/bin/bash

# Updating pangolin
echo "Updating pangolin data"
pangolin --update-data

# Update SARS-CoV-2 Nextclade data
echo "Updating SARS-CoV-2 Nextclade data"
nextclade dataset get --name='sars-cov-2' --output-dir="$PIPELINE/SARS-CoV-2/nextstrain_files/"

# Updating Influenza Nextclade data

echo "Updating Influenza Nextclade data"
nextclade dataset get --name='flu_h1n1pdm_ha' --output-dir="$PIPELINE/Influenza/nextclade_files/H1"
nextclade dataset get --name='flu_h3n2_ha' --output-dir="$PIPELINE/Influenza/nextclade_files/H3"
nextclade dataset get --name='flu_vic_ha' --output-dir="$PIPELINE/Influenza/nextclade_files/Vic"
nextclade dataset get --name='flu_yam_ha' --output-dir="$PIPELINE/Influenza/nextclade_files/Yam"

# Updating Dengue Nextclade data

echo "Updating Dengue Nextclade data"
nextclade dataset get --name='community/v-gen-lab/dengue/denv1' --output-dir="$PIPELINE/DENV/nextclade_files/denv1"
nextclade dataset get --name='community/v-gen-lab/dengue/denv2' --output-dir="$PIPELINE/DENV/nextclade_files/denv2"
nextclade dataset get --name='community/v-gen-lab/dengue/denv3' --output-dir="$PIPELINE/DENV/nextclade_files/denv3"
nextclade dataset get --name='community/v-gen-lab/dengue/denv4' --output-dir="$PIPELINE/DENV/nextclade_files/denv4"

echo "Done"
