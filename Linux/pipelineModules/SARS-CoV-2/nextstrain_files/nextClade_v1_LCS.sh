S=$1; # Arquivo All Stats
F=$2; # Arquivo All.fasta
H=$3;  # data
P=$4; # Number of processors to use
L=$5 # LCS output file

# Sort All Stats file for latter
sort -ru ${S} > ${S:0:-4}_sorted.tsv

# Update tree.json
#curl -OL https://raw.githubusercontent.com/nextstrain/nextclade/master/data/sars-cov-2/tree.json
#mv tree.json /project/carol/sarsCov2/nexstrain_files/tree.json # move file to suport folder 
nextclade dataset get --name='sars-cov-2' --output-dir='/project/carol/sarsCov2/nextstrain_files/'

# Run NextClade 

mkdir -p ${H}__nextClade # Create NextClade output folder
nextclade run  --input-dataset /project/carol/sarsCov2/nextstrain_files/ --genes E,M,N,ORF1a,ORF1b,ORF3a,ORF6,ORF7a,ORF7b,ORF8,ORF9b,S --output-all ${H}__nextClade --output-basename nextclade --jobs ${P} ${F}

# Join NextClade results with Sorted All Stats using python
python3 /project/carol/sarsCov2/nextstrain_files/process_AllStats_NextClade_LCS.py ${S:0:-4}_sorted.tsv ${H}__nextClade/nextclade.tsv ${L}

# Cleaning intermediate files
rm ${S:0:-4}_sorted.tsv
