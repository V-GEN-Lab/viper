FOLDER=${PWD##*/}
THREADS=$1

# Create list of samples based on fastq files
L=$(ls *.fastq.gz); for f in $L; do g=${f%.*}; echo ${g%_L00*};  done | uniq > Lista.txt

# Update nextclade datasets
nextclade dataset get --name='flu_h1n1pdm_ha' --output-dir='/storage/vital/volume1/carol/influenza/pipeline/nextclade_files/H1'
nextclade dataset get --name='flu_h3n2_ha' --output-dir='/storage/vital/volume1/carol/influenza/pipeline/nextclade_files/H3'
nextclade dataset get --name='flu_vic_ha' --output-dir='/storage/vital/volume1/carol/influenza/pipeline/nextclade_files/Vic'
nextclade dataset get --name='flu_yam_ha' --output-dir='/storage/vital/volume1/carol/influenza/pipeline/nextclade_files/Yam'

# Init assemblies in folder
eval "$(conda shell.bash hook)"
conda activate /storage/vital/volume1/carol/influenza/pipeline/fluAssembly; cat Lista.txt | xargs -P ${THREADS} -I {} sh -c "bash /storage/vital/volume1/carol/influenza/pipeline/Influenza_assembly_v4.sh {} ${FOLDER}"; conda deactivate

# Join all fasta and statistics files
#cat *.fasta > All_Fastas__${FOLDER}.fas
cat *_complete.Statistics | sort -ru > All_Statistics__${FOLDER}.tsv

# Process CeVIVAS output
assemlby_path=$(realpath .)
while read value; do read_path=$(realpath $value*R1*); echo -e $value '\t' $read_path; done < Lista.txt > read_path.tsv

/usr/bin/python /storage/vital/volume1/carol/influenza/pipeline/write_flu_CeVIVAS_output_v4.py All_Statistics__${FOLDER}.tsv $assemlby_path read_path.tsv $FOLDER

cat Genoma_FLU* > All_FLU_${FOLDER}.fasta
