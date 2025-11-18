FOLDER=${PWD##*/}
THREADS=$1

# Create list of samples based on fastq files
L=$(ls *.fastq.gz); for f in $L; do g=${f%.*}; echo ${g%_L00*};  done | uniq > Lista.txt
pangolin --update 

cat Lista.txt | xargs -P ${THREADS} -I {} sh -c "bash $PIPELINE/SARS-CoV-2/Pipeline_Illumina_v8_bowtie2_ref_iVar.sh {} ${FOLDER}";

# Join all fasta and statistics files
cat *.fasta > All_Fastas__${FOLDER}.fas
cat *.Statistics > All_Statistics__${FOLDER}.tsv 


# Run nextClade 
bash $PIPELINE/SARS-CoV-2/nextstrain_files/nextClade_v1.sh All_Statistics__${FOLDER}.tsv All_Fastas__${FOLDER}.fas ${FOLDER} ${THREADS}


