#!/usr/bin/env bash

#SBATCH --job-name=Rodar_Montagens_Covid
#SBATCH --output=Rodar_Montagens_Covid_%A_%a.out
#SBATCH --error=Rodar_Montagens_Covid_%A_%a.err


rm -rf output/ || true
rm *.err *.out  || true



# Rodar a partir do diretorio "/project/jspatane/sarsCov2/project/Mendelics_[*]".
# Os genomas já foram des-zipados pra esse teste. Se nao forem, usar linha abaixo:

# Command to submmit job.
# sbatch -n4 --mem=100G -w "vital" Rodar_Montagens_Covid.sh 2> log/sbatch.err > log/sbatch.out

# If you want to run again, remove contents from log and output directories.
#rm -rf log/* output/*

# Unzip all files.
#for F in *.zip; do unzip $F; done

START_TIME=`date +%s`


Lista_fileNames="Lista_Genomas.txt"
n=`wc -l ${Lista_fileNames} | cut -d' ' -f1`


#mkdir log 
mkdir output



for LINE in `seq ${n}`; do
        F=$(sed -n ${LINE}p ${Lista_fileNames})
        mkdir output/${F}
        bwa mem Ref_Wuhan.fasta ${F}.fastq | samtools sort -o output/${F}/${F}.sorted.bam &
done
wait


for LINE in `seq ${n}`; do
	F=$(sed -n ${LINE}p ${Lista_fileNames})
	samtools index output/${F}/${F}.sorted.bam &
done
wait

for LINE in `seq ${n}`; do
	F=$(sed -n ${LINE}p ${Lista_fileNames})
	samtools mpileup -uf Ref_Wuhan.fasta output/${F}/${F}.sorted.bam | bcftools call -c --ploidy 1 | vcfutils.pl vcf2fq > output/${F}/${F}.mpileup &
done
wait


for LINE in `seq ${n}`; do
	F=$(sed -n ${LINE}p ${Lista_fileNames})
	seqtk seq -aQ64 -q20 -n N output/${F}/${F}.mpileup > output/${F}/Genoma_${F}.fasta &
done
wait

END_TIME=`date +%s`
runtime=$( echo "${END_TIME} - ${START_TIME}" | bc -l )
hours=$((runtime / 3600))
minutes=$(( (runtime % 3600) / 60 ))
seconds=$(( (runtime % 3600) % 60 ))
echo "O script demorou: $hours h:$minutes min:$seconds s"

exit 0
