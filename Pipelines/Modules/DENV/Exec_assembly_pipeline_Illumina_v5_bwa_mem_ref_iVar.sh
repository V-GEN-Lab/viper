FOLDER=${PWD##*/}
THREADS=${1:-1}

# Create list of samples based on fastq files
L=$(ls *.fastq.gz)
for f in $L; do
  g=${f%.*}
  echo "${g%_L00*}"
done | uniq > Lista.txt

# Init assemblies in folder
eval "$(conda shell.bash hook)"
while IFS= read -r SAMPLE; do
  [ -z "$SAMPLE" ] && continue
  bash "$PIPELINE/DENV/Dengue_assembly_v5.1.sh" "$SAMPLE" "${FOLDER}" "${THREADS}"
done < Lista.txt
conda deactivate

# Join all fasta and statistics files
cat *.fasta > "All_Fastas__${FOLDER}.fas"
cat *.Statistics | sort -ru > "All_Statistics__${FOLDER}.tsv"

# Add QC status column (heurística: Coverage>=85 -> A; senão E_Coverage>=55 -> E; senão R)
STATS="All_Statistics__${FOLDER}.tsv"
TMP="${STATS}.tmp"

gawk -v FS='\t' -v OFS='\t' '
function is_missing(x) { return (x=="" || x=="NA" || x=="NaN") }

NR==1 {
  cov_i = ecov_i = 0
  for (i=1; i<=NF; i++) {
    if ($i == "Coverage")   cov_i = i
    if ($i == "E_Coverage") ecov_i = i
  }
  if (cov_i == 0)  { print "ERRO: coluna Coverage não encontrada no header" > "/dev/stderr"; exit 2 }
  if (ecov_i == 0) { print "ERRO: coluna E_Coverage não encontrada no header" > "/dev/stderr"; exit 2 }

  print $0, "Passed_QC"
  next
}

{
  cov  = $(cov_i)
  ecov = $(ecov_i)

  # se vier decimal com vírgula, normaliza
  gsub(",", ".", cov)
  gsub(",", ".", ecov)

  status = "R"
  if (is_missing(cov)) {
    status = "R"
  } else if ((cov + 0) >= 85.0) {
    status = "A"
  } else if (!is_missing(ecov) && (ecov + 0) >= 55.0) {
    status = "E"
  } else {
    status = "R"
  }

  print $0, status
}
' "$STATS" > "$TMP"

mv -f "$TMP" "$STATS"