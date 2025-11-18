FOLDER=${PWD##*/}
THREADS=$1

# Create list of samples based on fastq files
L=$(ls *.fastq.gz); for f in $L; do g=${f%.*}; echo ${g%_L00*};  done | uniq > Lista.txt

# Init assemblies in folder
eval "$(conda shell.bash hook)"
conda activate $PIPELINE/DENV/denvAssembly; cat Lista.txt | xargs -P ${THREADS} -I {} sh -c "bash $PIPELINE/DENV/Dengue_assembly_v5.1.sh {} ${FOLDER}"; conda deactivate
 
# Join all fasta and statistics files
cat *.fasta > All_Fastas__${FOLDER}.fas
cat *.Statistics | sort -ru > All_Statistics__${FOLDER}.tsv 

