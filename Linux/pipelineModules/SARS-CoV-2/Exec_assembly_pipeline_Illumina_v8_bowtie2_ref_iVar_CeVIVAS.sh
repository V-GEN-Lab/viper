FOLDER=${PWD##*/}
THREADS=$1

# Create list of samples based on fastq files
L=$(ls *.fastq.gz); for f in $L; do g=${f%.*}; echo ${g%_L00*};  done | uniq > Lista.txt

# Init assemblies in folder
eval "$(conda shell.bash hook)"
conda activate pangolin;
pangolin --update

cat Lista.txt | xargs -P ${THREADS} -I {} sh -c "bash /storage/zuleika/volume1/project/carol/sarsCov2/CeVIVAS/pipeline/Pipeline_Illumina_v8_bowtie2_ref_iVar.sh {} ${FOLDER}"; conda deactivate

# Join all fasta and statistics files
cat *.fasta > All_Fastas__${FOLDER}.fas
cat *.Statistics > All_Statistics__${FOLDER}.tsv 

#Update nextclade (tree.json)
#curl -OL https://raw.githubusercontent.com/nextstrain/nextclade/master/data/sars-cov-2/tree.json

# Esse repositório é antigo. Atualização de árvore via nextclade dataset get em nextClade_v1.sh
#curl -OL https://github.com/nextstrain/nextclade/blob/feat/otuput-json-dataset-info/data2/new/cov-sc2/tree.json
#mv -f tree.json /storage/zuleika/volume1/project/carol/sarsCov2/CeVIVAS/pipeline/nextstrain_files/tree.json # move file to suport folder

# Run nextClade 
bash /storage/zuleika/volume1/project/carol/sarsCov2/CeVIVAS/pipeline/nextstrain_files/nextClade_v1.sh All_Statistics__${FOLDER}.tsv All_Fastas__${FOLDER}.fas ${FOLDER} ${THREADS}

# Process CeVIVAS output
assemlby_path=$(realpath .)
while read value; do read_path=$(realpath $value*R1*); echo -e $value '\t' $read_path; done < Lista.txt > read_path.tsv

conda activate pangolin
python /storage/zuleika/volume1/project/carol/sarsCov2/CeVIVAS/pipeline/write_cov_CeVIVAS_output8.py All_Statistics__${FOLDER}.csv $assemlby_path read_path.tsv $FOLDER
conda deactivate

