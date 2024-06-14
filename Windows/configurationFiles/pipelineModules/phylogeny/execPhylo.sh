#!/bin/bash

# Modified version for VIPER
# Script for conducting phylogenetic sequence analysis
# Author: James Siqueira Pereira
# Supervisors: Alex Ranieri Jerônimo Lima; Gabriela Ribeiro; Vinicius Cairus;
# Funding Institutions: FAPESP; Instituto Butantan;
# Date: 15-03-2024

# Description: This script performs phylogenetic analyses on viral genomic data, covering sequential steps for sequence alignment, tree construction, temporal data-based refinement, ancestral reconstruction, translation, and export for subsequent interactive visualization using the auspice.us program.

# Note: This script can be adapted according to the user's needs. Its use and distribution are subject to crediting the authors and funding institutions.

# Set default variables
SEQUENCES=""
OUTPUT=""
NTHREADS="2"
METADATA="null"
CONFIG="$PIPELINE/phylogeny/auspice_config.json"
MODEL="TEST"
UFBOOT="1000"
LAT_LONGS="$PIPELINE/phylogeny/lat_longs.tsv"
# Process command line arguments
# Function to handle options
process_args() {
    local key="$1"
    local value="$2"

    case $key in
        -i|--sequences)
            SEQUENCES="$value"
            ;;
        -v|--virus)
            VIRUS="$value"
            ;;
        -j|--jobname)
            OUTPUT="$value"
            ;;
        -n|--nthreads)
            NTHREADS="$value"
            ;;
        -metadata|--metadata)
            METADATA="$value"
            ;;
        -m)
            MODEL="$value"
            ;;
        -b|--ufboot)
            UFBOOT="$value"
            ;;
        *)
            echo "Unknown option: $key"
            exit 1
            ;;
    esac
}

# Process command line arguments
while [[ $# -gt 0 ]]; do
    process_args "$1" "$2"
    shift 2
done

if [[ $VIRUS == "SARS-CoV-2" ]]; then
    CLADES="$PIPELINE/phylogeny/mutation_tables/sars-cov-2_clades.tsv"
    REFERENCE_SEQUENCE="$PIPELINE/phylogeny/references/Ref_Wuhan.fasta"
    REFERENCE_ANNOTATION="$PIPELINE/phylogeny/references/Ref_Wuhan.gb"
    ROOT="MN908947"
    #TITLE='"Phylogenetic tree built using VIPER for SARS-CoV-2"'
    if [ "$METADATA" != "null" ]; then
        # Check if $METADATA has the date column
        if awk -F'\t' 'NR==1 {for (i=1; i<=NF; i++) if ($i == "date") { found=1; break } } END { exit !found }' "${METADATA}"; then
            # Read the header to find the column positions
            IFS=$'\t' read -r -a headers < "${METADATA}"
            for i in "${!headers[@]}"; do
            if [ "${headers[$i]}" == "strain" ]; then
                strain_pos=$i
            elif [ "${headers[$i]}" == "date" ]; then
                date_pos=$i
            fi
            done
            # Prepare the new line with the values
            num_columns=${#headers[@]}
            new_line=()

            for (( i=0; i<num_columns; i++ )); do
            if [ $i -eq $strain_pos ]; then
                new_line+=("MN908947")
            elif [ $i -eq $date_pos ]; then
                new_line+=("2019-12-31")
            else
                new_line+=("")
            fi
            done

            # Join the new line with tabs
            new_line_joined=$(IFS=$'\t'; echo "${new_line[*]}")

            # Copy the input file to the output file
            cp "${METADATA}" "${METADATA:0:-4}_formated.tsv"

            # Append the new line to the output file
            echo -e "$new_line_joined" >> "${METADATA:0:-4}_formated.tsv"
            METADATA2="${METADATA:0:-4}_formated.tsv"
            METADATA=$METADATA2
            ROOT="oldest"
        fi
    fi
elif [[ $VIRUS == "DENV-1" ]]; then
    CLADES="$PIPELINE/phylogeny/mutation_tables/DENV1_augur_clades.tsv"
    REFERENCE_SEQUENCE="$PIPELINE/phylogeny/references/denv1.fasta"
    REFERENCE_ANNOTATION="$PIPELINE/phylogeny/references/denv1.gb"
    ROOT="mid_point"
    #TITLE='"Phylogenetic tree built using VIPER for DENV-1"'
elif [[ $VIRUS == "DENV-2" ]]; then
    CLADES="$PIPELINE/phylogeny/mutation_tables/DENV2_augur_clades.tsv"
    REFERENCE_SEQUENCE="$PIPELINE/phylogeny/references/denv2.fasta"
    REFERENCE_ANNOTATION="$PIPELINE/phylogeny/references/denv2.gb"
    ROOT="mid_point"
    #TITLE='"Phylogenetic tree built using VIPER for DENV-2"'
elif [[ $VIRUS == "DENV-3" ]]; then
    CLADES="$PIPELINE/phylogeny/mutation_tables/DENV3_augur_clades.tsv"
    REFERENCE_SEQUENCE="$PIPELINE/phylogeny/references/denv3.fasta"
    REFERENCE_ANNOTATION="$PIPELINE/phylogeny/references/denv3.gb"
    ROOT="mid_point"
    #TITLE='"Phylogenetic tree built using VIPER for DENV-3"'
elif [[ $VIRUS == "DENV-4" ]]; then
    CLADES="$PIPELINE/phylogeny/mutation_tables/DENV4_augur_clades.tsv"
    REFERENCE_SEQUENCE="$PIPELINE/phylogeny/references/denv4.fasta"
    REFERENCE_ANNOTATION="$PIPELINE/phylogeny/references/denv4.gb"
    ROOT="mid_point"
    #TITLE='"Phylogenetic tree built using VIPER for DENV-4"'
elif [[ $VIRUS == "FLUA_H1" ]]; then
    CLADES="$PIPELINE/phylogeny/mutation_tables/h1n1pdm09_ha_clades-long.tsv"
    REFERENCE_SEQUENCE="$PIPELINE/phylogeny/references/h1n1pdm09_ha_reference.fasta"
    REFERENCE_ANNOTATION="$PIPELINE/phylogeny/references/h1n1pdm09_ha.gb"
    ROOT="best"
    #TITLE='"Phylogenetic tree built using VIPER for Influenza A/H1N1 - Hemagglutinin (segment 4)"'
    if [ "$METADATA" != "null" ]; then
        # Check if $METADATA has the date column
        if awk -F'\t' 'NR==1 {for (i=1; i<=NF; i++) if ($i == "date") { found=1; break } } END { exit !found }' "${METADATA}"; then
            ROOT="oldest"
        fi
    fi
elif [[ $VIRUS == "FLUA_N1" ]]; then
    CLADES=""
    REFERENCE_SEQUENCE="$PIPELINE/phylogeny/references/h1n1pdm09_na_reference.fasta"
    REFERENCE_ANNOTATION="$PIPELINE/phylogeny/references/h1n1pdm09_na.gb"
    ROOT="best"
    #TITLE='"Phylogenetic tree built using VIPER for Influenza A/H1N1 - Neuraminidase (segment 6)"'
    if [  "$METADATA" != "null" ]; then
        # Check if $METADATA has the date column
        if awk -F'\t' 'NR==1 {for (i=1; i<=NF; i++) if ($i == "date") { found=1; break } } END { exit !found }' "${METADATA}"; then
            ROOT="oldest"
        fi
    fi
elif [[ $VIRUS == "FLUA_N2" ]]; then
    CLADES=""
    REFERENCE_SEQUENCE="$PIPELINE/phylogeny/references/h3n2_na_reference.fasta"
    REFERENCE_ANNOTATION="$PIPELINE/phylogeny/references/h3n2_na.gb"
    ROOT="best"
    #TITLE='"Phylogenetic tree built using VIPER for Influenza A/H3N2 - Neuraminidase (segment 6)"'
    if [ "$METADATA" != "null" ]; then
        # Check if $METADATA has the date column
        if awk -F'\t' 'NR==1 {for (i=1; i<=NF; i++) if ($i == "date") { found=1; break } } END { exit !found }' "${METADATA}"; then
            ROOT="oldest"
        fi
    fi
elif [[ $VIRUS == "FLUA_H3" ]]; then
    CLADES="$PIPELINE/phylogeny/mutation_tables/h3n2_ha_clades-long.tsv"
    REFERENCE_SEQUENCE="$PIPELINE/phylogeny/references/h3n2_ha_reference.fasta"
    REFERENCE_ANNOTATION="$PIPELINE/phylogeny/references/h3n2_ha.gb"
    ROOT="best"
    #TITLE='"Phylogenetic tree built using VIPER for Influenza A/H3N2 - Hemagglutinin (segment 4)"'
    if [ "$METADATA" != "null" ]; then
        # Check if $METADATA has the date column
        if awk -F'\t' 'NR==1 {for (i=1; i<=NF; i++) if ($i == "date") { found=1; break } } END { exit !found }' "${METADATA}"; then
            ROOT="oldest"
        fi
    fi
elif [[ $VIRUS == "FLUB_VIC_HA" ]]; then
    CLADES="$PIPELINE/phylogeny/mutation_tables/vic_ha_clades.tsv"
    REFERENCE_SEQUENCE="$PIPELINE/phylogeny/references/vic_ha_reference.fasta"
    REFERENCE_ANNOTATION="$PIPELINE/phylogeny/references/vic_ha.gb"
    ROOT="best"
    #TITLE='"Phylogenetic tree built using VIPER for Influenza B/Victoria - Hemagglutinin (segment 4)"'
    if [ "$METADATA" != "null" ]; then
        # Check if $METADATA has the date column
        if awk -F'\t' 'NR==1 {for (i=1; i<=NF; i++) if ($i == "date") { found=1; break } } END { exit !found }' "${METADATA}"; then
            ROOT="oldest"
        fi
    fi
elif [[ $VIRUS == "FLUB_VIC_NA" ]]; then
    CLADES=""
    REFERENCE_SEQUENCE="$PIPELINE/phylogeny/references/vic_na_reference.fasta"
    REFERENCE_ANNOTATION="$PIPELINE/phylogeny/references/vic_na.gb"
    ROOT="best"
    #TITLE='"Phylogenetic tree built using VIPER for Influenza B/Victoria - Neuraminidase (segment 6)"'
    if [ "$METADATA" != "null" ]; then
        # Check if $METADATA has the date column
        if awk -F'\t' 'NR==1 {for (i=1; i<=NF; i++) if ($i == "date") { found=1; break } } END { exit !found }' "${METADATA}"; then
            ROOT="oldest"
        fi
    fi
elif [[ $VIRUS == "FLUB_YAM_HA" ]]; then
    CLADES="$PIPELINE/phylogeny/mutation_tables/yam_ha_clades.tsv"
    REFERENCE_SEQUENCE="$PIPELINE/phylogeny/references/yam_ha_reference.fasta"
    REFERENCE_ANNOTATION="$PIPELINE/phylogeny/references/yam_ha.gb"
    ROOT="best"
    #TITLE='"Phylogenetic tree built using VIPER for Influenza B/Yamagata - Hemagglutinin (segment 4)"'
    if [ "$METADATA" != "null" ]; then
        # Check if $METADATA has the date column
        if awk -F'\t' 'NR==1 {for (i=1; i<=NF; i++) if ($i == "date") { found=1; break } } END { exit !found }' "${METADATA}"; then
            ROOT="oldest"
        fi
    fi
elif [[ $VIRUS == "FLUB_YAM_NA" ]]; then
    CLADES=""
    REFERENCE_SEQUENCE="$PIPELINE/phylogeny/references/yam_na_reference.fasta"
    REFERENCE_ANNOTATION="$PIPELINE/phylogeny/references/yam_na.gb"
    ROOT="best"
    #TITLE='"Phylogenetic tree built using VIPER for Influenza B/Yamagata - Neuraminidase (segment 6)"'
    if [ "$METADATA" != "null" ]; then
        # Check if $METADATA has the date column
        if awk -F'\t' 'NR==1 {for (i=1; i<=NF; i++) if ($i == "date") { found=1; break } } END { exit !found }' "${METADATA}"; then
            ROOT="oldest"
        fi
    fi
else
    echo "Unknown virus: $VIRUS"
    # Handle the case when $VIRUS is not recognized
    exit 1
fi

# Create alignment folder
mkdir ${OUTPUT}_alignment

echo "Sequence alignment in progress..."
# Command 1: nextalign
echo "Executing command:"
echo "augur align --sequences $SEQUENCES --reference-sequence $REFERENCE_SEQUENCE --output ${OUTPUT}_alignment/${OUTPUT}_aligned.fasta --fill-gaps --nthreads $NTHREADS"
augur align --sequences $SEQUENCES --reference-sequence $REFERENCE_SEQUENCE --output ${OUTPUT}_alignment/${OUTPUT}_aligned.fasta --fill-gaps --nthreads $NTHREADS

echo "Tree construction using IQTree2..."

# Create output folder
mkdir ${OUTPUT}_tree

# Copy alignment to folder where the tree will be generated
cp ${OUTPUT}_alignment/${OUTPUT}_aligned.fasta ${OUTPUT}_tree/

# Command 2: IQtree2
echo "Executing command:"
echo "iqtree2 -s ${OUTPUT}_tree/${OUTPUT}_aligned.fasta -m $MODEL -nt $NTHREADS -bb $UFBOOT"
iqtree2 -s ${OUTPUT}_tree/${OUTPUT}_aligned.fasta -m $MODEL -nt $NTHREADS -bb $UFBOOT

# Remove duplicate alignment
rm -rf ${OUTPUT}_tree/${OUTPUT}_aligned.fasta

# Rename output file
mv ${OUTPUT}_tree/${OUTPUT}_aligned.fasta.treefile ${OUTPUT}_tree/${OUTPUT}.treefile

# Command 3: augur refine
if [  "$METADATA" != "null" ]; then
  echo "Refining tree with TreeTime..."
  echo "Executing command:"
  echo "augur refine --alignment ${OUTPUT}_alignment/${OUTPUT}_aligned.fasta --tree ${OUTPUT}_tree/${OUTPUT}.treefile --metadata $METADATA --output-tree ${OUTPUT}_tree/${OUTPUT}_refinedTree.nwk --root $ROOT --output-node-data ${OUTPUT}_tree/${OUTPUT}_branch-lengths.json --timetree --coalescent const --date-confidence --stochastic-resolve --date-inference marginal --clock-filter-iqd 4"
  augur refine --alignment ${OUTPUT}_alignment/${OUTPUT}_aligned.fasta --tree ${OUTPUT}_tree/${OUTPUT}.treefile --metadata $METADATA --output-tree ${OUTPUT}_tree/${OUTPUT}_refinedTree.nwk --root $ROOT --output-node-data ${OUTPUT}_tree/${OUTPUT}_branch-lengths.json --timetree --coalescent const --date-confidence --stochastic-resolve --date-inference marginal --clock-filter-iqd 4
else
  echo "Refining undated tree..."
  echo "Executing command:"
  echo "augur refine --tree ${OUTPUT}_tree/${OUTPUT}.treefile --root $ROOT --output-tree ${OUTPUT}_tree/${OUTPUT}_refinedTree.nwk --output-node-data ${OUTPUT}_tree/${OUTPUT}_branch-lengths.json"
  augur refine --tree ${OUTPUT}_tree/${OUTPUT}.treefile --root $ROOT --output-tree ${OUTPUT}_tree/${OUTPUT}_refinedTree.nwk --output-node-data ${OUTPUT}_tree/${OUTPUT}_branch-lengths.json
fi

# Command 4: augur ancestral
echo "Reconstructing ancestral..."
echo "Executing command:"
echo "augur ancestral --tree ${OUTPUT}_tree/${OUTPUT}_refinedTree.nwk --alignment ${OUTPUT}_alignment/${OUTPUT}_aligned.fasta --output-node-data ${OUTPUT}_tree/${OUTPUT}_nt_muts.json --inference joint --root-sequence $REFERENCE_SEQUENCE"
augur ancestral --tree ${OUTPUT}_tree/${OUTPUT}_refinedTree.nwk --alignment ${OUTPUT}_alignment/${OUTPUT}_aligned.fasta --output-node-data ${OUTPUT}_tree/${OUTPUT}_nt_muts.json --inference joint --root-sequence $REFERENCE_SEQUENCE

# Command 5: augur translate
echo "Defining mutations in amino acid sequences..."
echo "Executing command:"
echo "augur translate --tree ${OUTPUT}_tree/${OUTPUT}_refinedTree.nwk --ancestral-sequences ${OUTPUT}_tree/${OUTPUT}_nt_muts.json --reference-sequence $REFERENCE_ANNOTATION --output-node-data ${OUTPUT}_tree/${OUTPUT}_aa_muts.json"
augur translate --tree ${OUTPUT}_tree/${OUTPUT}_refinedTree.nwk --ancestral-sequences ${OUTPUT}_tree/${OUTPUT}_nt_muts.json --reference-sequence $REFERENCE_ANNOTATION --output-node-data ${OUTPUT}_tree/${OUTPUT}_aa_muts.json

# Command 6: augur clades for lineages
if [ -n "$CLADES" ]; then
  echo "Annotating clades/lineages..."
  echo "Executing command:"
  echo "augur clades --tree ${OUTPUT}_tree/${OUTPUT}_refinedTree.nwk --mutations ${OUTPUT}_tree/${OUTPUT}_aa_muts.json ${OUTPUT}_tree/${OUTPUT}_nt_muts.json ${OUTPUT}_tree/${OUTPUT}_branch-lengths.json --clade ${CLADES} --output-node-data ${OUTPUT}_tree/lineages.json --membership-name lineage --label-name lineages"
  augur clades --tree ${OUTPUT}_tree/${OUTPUT}_refinedTree.nwk --mutations ${OUTPUT}_tree/${OUTPUT}_aa_muts.json ${OUTPUT}_tree/${OUTPUT}_nt_muts.json ${OUTPUT}_tree/${OUTPUT}_branch-lengths.json --clade ${CLADES} --output-node-data ${OUTPUT}_tree/lineages.json --membership-name clade_or_lineage --label-name clade_or_lineage
fi

# Extract metadata columns to color by
if [  "$METADATA" != "null" ]; then
    COLOR_METADATA=$(awk -F'\t' 'NR==1 {for (i=1; i<=NF; i++) if ($i != "strain" && $i != "date") printf "%s ", $i; printf "\n"; exit}' $METADATA)
fi

# Command 7: augur export v2
if [  "$METADATA" != "null" ] && [ -n "$CLADES" ]; then
  echo "Exporting tree with metadata..."
  echo "Executing command:"
  echo "augur export v2 --tree ${OUTPUT}_tree/${OUTPUT}_refinedTree.nwk --metadata $METADATA --node-data ${OUTPUT}_tree/${OUTPUT}_nt_muts.json ${OUTPUT}_tree/${OUTPUT}_aa_muts.json ${OUTPUT}_tree/${OUTPUT}_branch-lengths.json ${OUTPUT}_tree/lineages.json --output ${OUTPUT}_tree/${OUTPUT}.json --color-by-metadata ${COLOR_METADATA} --auspice-config $CONFIG --lat-longs $LAT_LONGS"
  augur export v2 --tree ${OUTPUT}_tree/${OUTPUT}_refinedTree.nwk --metadata $METADATA --node-data ${OUTPUT}_tree/${OUTPUT}_nt_muts.json ${OUTPUT}_tree/${OUTPUT}_aa_muts.json ${OUTPUT}_tree/${OUTPUT}_branch-lengths.json ${OUTPUT}_tree/lineages.json --output ${OUTPUT}_tree/${OUTPUT}.json --color-by-metadata ${COLOR_METADATA} --auspice-config $CONFIG --lat-longs $LAT_LONGS --skip-validation
elif [  "$METADATA" != "null" ] && [ -z "$CLADES" ]; then
  echo "Exporting tree with metadata..."
  echo "Executing command:"
  echo "augur export v2 --tree ${OUTPUT}_tree/${OUTPUT}_refinedTree.nwk --metadata $METADATA --node-data ${OUTPUT}_tree/${OUTPUT}_nt_muts.json ${OUTPUT}_tree/${OUTPUT}_aa_muts.json ${OUTPUT}_tree/${OUTPUT}_branch-lengths.json --output ${OUTPUT}_tree/${OUTPUT}.json --color-by-metadata ${COLOR_METADATA} --auspice-config $CONFIG --lat-longs $LAT_LONGS"
  augur export v2 --tree ${OUTPUT}_tree/${OUTPUT}_refinedTree.nwk --metadata $METADATA --node-data ${OUTPUT}_tree/${OUTPUT}_nt_muts.json ${OUTPUT}_tree/${OUTPUT}_aa_muts.json ${OUTPUT}_tree/${OUTPUT}_branch-lengths.json --output ${OUTPUT}_tree/${OUTPUT}.json --color-by-metadata ${COLOR_METADATA} --auspice-config $CONFIG --lat-longs $LAT_LONGS --skip-validation
elif [ "$METADATA" == "null" ] && [ -n "$CLADES" ]; then
  echo "Exporting tree without metadata..."
  echo "Executing command:"
  echo "augur export v2 --tree ${OUTPUT}_tree/${OUTPUT}_refinedTree.nwk --node-data ${OUTPUT}_tree/${OUTPUT}_nt_muts.json ${OUTPUT}_tree/${OUTPUT}_aa_muts.json ${OUTPUT}_tree/${OUTPUT}_branch-lengths.json ${OUTPUT}_tree/lineages.json --output ${OUTPUT}_tree/${OUTPUT}.json --auspice-config $CONFIG"
  augur export v2 --tree ${OUTPUT}_tree/${OUTPUT}_refinedTree.nwk --node-data ${OUTPUT}_tree/${OUTPUT}_nt_muts.json ${OUTPUT}_tree/${OUTPUT}_aa_muts.json ${OUTPUT}_tree/${OUTPUT}_branch-lengths.json ${OUTPUT}_tree/lineages.json --output ${OUTPUT}_tree/${OUTPUT}.json --auspice-config $CONFIG --skip-validation
else
  echo "Exporting tree without metadata..."
  echo "Executing command:"
  echo "augur export v2 --tree ${OUTPUT}_tree/${OUTPUT}_refinedTree.nwk --node-data ${OUTPUT}_tree/${OUTPUT}_nt_muts.json ${OUTPUT}_tree/${OUTPUT}_aa_muts.json ${OUTPUT}_tree/${OUTPUT}_branch-lengths.json ${OUTPUT}_tree/lineages.json --output ${OUTPUT}_tree/${OUTPUT}.json --auspice-config $CONFIG"
  augur export v2 --tree ${OUTPUT}_tree/${OUTPUT}_refinedTree.nwk --node-data ${OUTPUT}_tree/${OUTPUT}_nt_muts.json ${OUTPUT}_tree/${OUTPUT}_aa_muts.json ${OUTPUT}_tree/${OUTPUT}_branch-lengths.json ${OUTPUT}_tree/lineages.json --output ${OUTPUT}_tree/${OUTPUT}.json --auspice-config $CONFIG --skip-validation
  
fi

echo "Finish!!"
exit 0
