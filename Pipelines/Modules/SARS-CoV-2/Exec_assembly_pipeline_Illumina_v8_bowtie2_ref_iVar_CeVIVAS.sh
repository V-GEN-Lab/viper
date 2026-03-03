FOLDER=${PWD##*/}
THREADS=${1:-1}

# Create list of samples based on fastq files
L=$(ls *.fastq.gz); for f in $L; do g=${f%.*}; echo ${g%_L00*};  done | uniq > Lista.txt
pangolin --update-data

while IFS= read -r SAMPLE; do
  [ -z "$SAMPLE" ] && continue
  bash "$PIPELINE/SARS-CoV-2/Pipeline_Illumina_v8_bowtie2_ref_iVar.sh" "$SAMPLE" "${FOLDER}" "${THREADS}"
done < Lista.txt

# Join all fasta and statistics files
cat *.fasta > All_Fastas__${FOLDER}.fas
cat *.Statistics > All_Statistics__${FOLDER}.tsv 


# Run nextClade 
bash $PIPELINE/SARS-CoV-2/nextClade_v1.sh All_Statistics__${FOLDER}.tsv All_Fastas__${FOLDER}.fas ${FOLDER} ${THREADS}
rm -rf $PIPELINE/SARS-CoV-2/nextstrain_files/nextClade_v1.sh All_Statistics__${FOLDER}.tsv
