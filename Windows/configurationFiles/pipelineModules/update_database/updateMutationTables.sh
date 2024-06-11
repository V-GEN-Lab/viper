#!/bin/bash
# Update mutation tables
echo "Updating mutation tables..."
wget https://raw.githubusercontent.com/V-GEN-Lab/DyDEVILS/main/mutation_tables/DENV1_augur_clades.tsv -O $PIPELINE/phylogeny/mutation_tables/DENV1_augur_clades.tsv
wget https://raw.githubusercontent.com/V-GEN-Lab/DyDEVILS/main/mutation_tables/DENV2_augur_clades.tsv -O $PIPELINE/phylogeny/mutation_tables/DENV2_augur_clades.tsv
wget https://raw.githubusercontent.com/V-GEN-Lab/DyDEVILS/main/mutation_tables/DENV3_augur_clades.tsv -O $PIPELINE/phylogeny/mutation_tables/DENV3_augur_clades.tsv
wget https://raw.githubusercontent.com/V-GEN-Lab/DyDEVILS/main/mutation_tables/DENV4_augur_clades.tsv -O $PIPELINE/phylogeny/mutation_tables/DENV4_augur_clades.tsv
wget https://raw.githubusercontent.com/nextstrain/seasonal-flu/master/config/h1n1pdm/ha/clades-long.tsv -O $PIPELINE/phylogeny/mutation_tables/h1n1pdm09_ha_clades-long.tsv
wget https://raw.githubusercontent.com/nextstrain/seasonal-flu/master/config/h3n2/ha/clades-long.tsv -O $PIPELINE/phylogeny/mutation_tables/h3n2_ha_clades-long.tsv
wget https://raw.githubusercontent.com/nextstrain/seasonal-flu/master/config/nextstrain_clades_vic_ha.tsv -O $PIPELINE/phylogeny/mutation_tables/vic_ha_clades.tsv
wget https://raw.githubusercontent.com/nextstrain/seasonal-flu/master/config/nextstrain_clades_yam_ha.tsv -O $PIPELINE/phylogeny/mutation_tables/yam_ha_clades.tsv
wget https://raw.githubusercontent.com/nextstrain/ncov/master/defaults/clades.tsv -O $PIPELINE/phylogeny/mutation_tables/sars-cov-2_clades.tsv
echo "Done!"