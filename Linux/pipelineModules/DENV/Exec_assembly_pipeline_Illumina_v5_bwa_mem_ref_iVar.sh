FOLDER=${PWD##*/}
THREADS=$1

# Create list of samples based on fastq files
L=$(ls *.fastq.gz); for f in $L; do g=${f%.*}; echo ${g%_L00*};  done | uniq > Lista.txt

# Init assemblies in folder
eval "$(conda shell.bash hook)"
conda activate /project/carol/dengue/pipeline/denvAssembly; cat Lista.txt | xargs -P ${THREADS} -I {} sh -c "bash /project/carol/dengue/pipeline/Dengue_assembly_v5.1.sh {} ${FOLDER}"; conda deactivate
 
# Join all fasta and statistics files
cat *.fasta > All_Fastas__${FOLDER}.fas
cat *.Statistics | sort -ru > All_Statistics__${FOLDER}.tsv 

# Process CeVIVAS output
assemlby_path=$(realpath .)

while read value; do read_path=$(realpath $value*R1*); echo -e $value '\t' $read_path; done < Lista.txt > read_path.tsv

/usr/bin/python /project/carol/dengue/pipeline/write_dengue_CeVIVAS_output5.py All_Statistics__${FOLDER}.tsv $assemlby_path read_path.tsv $FOLDER
