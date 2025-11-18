S=$1; # Arquivo All Stats
F=$2; # Arquivo All.fasta
H=$3;  # data

# Correção: método mais robusto para criar nome do arquivo sorted
BASENAME=$(basename "${S}" .tsv)
SORTED_FILE="${BASENAME}_sorted.tsv"

# Sort All Stats file for latter
sort -ru "${S}" > "${SORTED_FILE}"

nextclade dataset get --name='sars-cov-2' --output-dir='$PIPELINE/SARS-CoV-2/nextstrain_files'

# Run NextClade 

mkdir -p ${H}__nextClade # Create NextClade output folder
nextclade run  --input-dataset $PIPELINE/SARS-CoV-2/nextstrain_files --output-all ${H}__nextClade --output-basename nextclade --jobs 4 ${F}

# Join NextClade results with Sorted All Stats using python
python $PIPELINE/SARS-CoV-2/nextstrain_files/process_AllStats_NextClade.py ${S:0:-4}_sorted.tsv ${H}__nextClade/nextclade.tsv

# Cleaning intermediate files
rm "${SORTED_FILE}"
