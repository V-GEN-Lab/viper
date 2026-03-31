FOLDER=${PWD##*/}
THREADS=${1:-1}

# Create list of samples based on fastq files
L=$(ls *.fastq.gz); for f in $L; do g=${f%.*}; echo ${g%_L00*};  done | uniq > Lista.txt

RAW_DATA_PATH="${PWD}/read_path.tsv"
: > "${RAW_DATA_PATH}"
shopt -s nullglob
while IFS= read -r SAMPLE; do
  [ -z "$SAMPLE" ] && continue
  sample_files=( "${PWD}"/"${SAMPLE}"*.fastq.gz )
  if [ ${#sample_files[@]} -gt 0 ]; then
    sample_path=$(printf '%s;' "${sample_files[@]}")
    sample_path=${sample_path%;}
  else
    sample_path="${PWD}"
  fi
  printf "%s\t%s\n" "$SAMPLE" "$sample_path" >> "${RAW_DATA_PATH}"
done < Lista.txt
shopt -u nullglob

# Update nextclade datasets
nextclade dataset get --name='flu_h1n1pdm_ha' --output-dir="$PIPELINE/Influenza/nextclade_files/H1"
nextclade dataset get --name='flu_h3n2_ha' --output-dir="$PIPELINE/Influenza/nextclade_files/H3"
nextclade dataset get --name='flu_vic_ha'     --output-dir="$PIPELINE/Influenza/nextclade_files/Vic"
nextclade dataset get --name='flu_yam_ha'     --output-dir="$PIPELINE/Influenza/nextclade_files/Yam"

# Init assemblies in folder
while IFS= read -r SAMPLE; do
  [ -z "$SAMPLE" ] && continue
  bash "$PIPELINE/Influenza/Influenza_assembly_v4.sh" "$SAMPLE" "${FOLDER}" "${THREADS}"
done < Lista.txt

# Join all fasta and statistics files
cat *_complete.Statistics | sort -ru > All_Statistics__${FOLDER}.tsv
cat Genoma_FLU* > All_FLU_${FOLDER}.fasta

# ------------------------------------------------------------
# Add QC status column (Influenza): based on segment_4_Coverage
# A: >=85 | E: 50-84.999 | R: <50 or not_detected/could_not_be_assembled/invalid
# ------------------------------------------------------------
STATS="All_Statistics__${FOLDER}.tsv"
TMP="${STATS}.tmp"

gawk -v FS='\t' -v OFS='\t' '
function is_missing(x) { return (x=="" || x=="NA" || x=="NaN") }

NR==1 {
  seg4_i = 0
  for (i=1; i<=NF; i++) {
    if ($i == "segment_4_Coverage") seg4_i = i
  }
  if (seg4_i == 0) { print "ERRO: coluna segment_4_Coverage não encontrada no header" > "/dev/stderr"; exit 2 }

  print $0, "Passed_QC"
  next
}

{
  seg4 = $(seg4_i)
  gsub(",", ".", seg4)  # se vier decimal com vírgula

  status = "R"

  if (seg4 == "not_detected" || seg4 == "could_not_be_assembled" || is_missing(seg4)) {
    status = "R"
  } else {
    # se não for número, vira R
    if (seg4 ~ /^-?[0-9]+([.][0-9]+)?$/) {
      val = seg4 + 0
      if (val >= 85)       status = "A"
      else if (val >= 50)  status = "E"
      else                 status = "R"
    } else {
      status = "R"
    }
  }

  print $0, status
}
' "$STATS" > "$TMP"

mv -f "$TMP" "$STATS"
python "$PIPELINE/Influenza/write_flu_CeVIVAS_output_v5.py" "${STATS}" "${PWD}" "${RAW_DATA_PATH}" "${FOLDER}"
