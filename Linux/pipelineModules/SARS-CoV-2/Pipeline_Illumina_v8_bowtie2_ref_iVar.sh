# 1. Definição de Variáveis e Preparação Inicial

G=$1;  # Nome base da amostra
H=$2;  # Data
F=Indels__${G}__${H};  # Nome do diretório de saída


# Criar diretório para cada genoma a ser processado
mkdir -p ${F};

# Descompactar e concatenar sequências
cat ${G}*R1*.fastq.gz  > ${F}/${F}_R1.fq.gz
cat ${G}*R2*.fastq.gz  > ${F}/${F}_R2.fq.gz

# Entrar no diretório de análise
cd ${F}

# 1) FastQC dos dados brutos
mkdir -p fastqc_raw
fastqc *_R1.fq.gz *_R2.fq.gz -o fastqc_raw -t 10

# 2. Limpeza das Reads
# Limpar reads usando Trimmomatic ou fastp
trimmomatic PE -phred33 ${F}_R1.fq.gz ${F}_R2.fq.gz  ${F}_R1_paired.fq.gz  ${F}_R1_unpaired.fq.gz  ${F}_R2_paired.fq.gz  ${F}_R2_unpaired.fq.gz  ILLUMINACLIP:$PIPELINE/SARS-CoV-2/primer_and_adapter_colection_v5.fasta:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36  -threads 10 > ${F}_trimmomatic.log

# 3) FastQC dos dados pós-trim (use apenas os *_paired*)
mkdir -p fastqc_trimmed
fastqc ${F}_R1_paired.fq.gz ${F}_R2_paired.fq.gz -o fastqc_trimmed -t 10

# 4) (Opcional, mas altamente recomendado) MultiQC para agregar tudo
# instala se não tiver: pip install --user multiqc
python3 -m multiqc . -o ${F}_multiqc_report

# Juntar reads unpaired
cat ${F}_R1_unpaired.fq.gz ${F}_R2_unpaired.fq.gz > ${F}_unpaired.fq.gz
rm -rf ${F}_R1_unpaired.fq.gz ${F}_R2_unpaired.fq.gz

# 3. Mapeamento Inicial e Correção de Genoma
# Mapeamento inicial com Bowtie2, posso usar também BWA-MEM ou o Minimap2
bowtie2 -p 5 -x $PIPELINE/SARS-CoV-2/Reference_Genomes/BA.2_ref -1 ${F}_R1_paired.fq.gz -2 ${F}_R2_paired.fq.gz -U ${F}_unpaired.fq.gz | samtools view -S -b -q 30 > ${F}_Wuhan_mapped.bam

# Ordenar e indexar o BAM
samtools sort ${F}_Wuhan_mapped.bam -o ${F}_Wuhan_mapped.sorted.bam
samtools index ${F}_Wuhan_mapped.sorted.bam

# Talvez adicionar Picard MarkDuplicates
#ferramenta que identifica e marca (ou remove, se configurado) leituras duplicadas em arquivos SAM/BAM. Essas duplicatas normalmente surgem durante a amplificação por PCR na preparação da biblioteca de sequenciamento (ou devido a artefatos ópticos) e podem inflar artificialmente a cobertura, levando a vieses em análises subsequentes, como a chamada de variantes. Ao marcar essas leituras duplicadas (geralmente definindo uma flag específica no arquivo BAM), os programas de downstream podem ignorá-las, garantindo que apenas leituras independentes sejam consideradas nas análises.

# java -jar picard.jar MarkDuplicates I=input.bam O=marked_duplicates.bam M=marked_dup_metrics.txt

# Correção do genoma com Pilon
java -jar /usr/local/src/pilon/1.24/pilon-1.24.jar --genome $PIPELINE/SARS-CoV-2/Reference_Genomes/BA.2_ref.fasta --frags ${F}_Wuhan_mapped.sorted.bam --minmq 30 --output ${F}.Pilon --fix "gaps,indels"  --threads 10 --mindepth 10;


# Correção do genoma com Pilon usando a versão do Picard
#java -jar /usr/local/src/pilon/1.24/pilon-1.24.jar --genome $PIPELINE/SARS-CoV-2/Reference_Genomes/BA.2_ref.fasta --frags marked_duplicates.bam --minmq 30 --output ${F}.Pilon --fix "gaps,indels" --threads 10 --mindepth 10;

# 4. Remapeamento e Chamada de Variantes
# Remapeamento com o genoma corrigido usando bowtie2 como mapeados e ara correção o ivar (mas pode ser usado LoFreq)
bowtie2-build ${F}.Pilon.fasta ${F}_pilon_index
bowtie2 -p 1 -x ${F}_pilon_index -1 ${F}_R1_paired.fq.gz -2 ${F}_R2_paired.fq.gz -U ${F}_unpaired.fq.gz | samtools view -S -b -q 30 > ${F}.Pilon.bam

# Ordenar e indexar o BAM
samtools sort ${F}.Pilon.bam -o ${F}.Pilon.sorted.bam
samtools index ${F}.Pilon.sorted.bam

# Chamada de consenso com iVar (acho que essa parte esta errada)
# -q 30 é a qualidade da base ...talvez abaixar para 20 ou 15, assim teremos uma cobertura mais completa em regiões de baixa profundidade
samtools mpileup -aa -A -d 0 -Q 0 ${F}.Pilon.sorted.bam | ivar consensus -q 30 -p ${F}_ivar -i ${F}.Pilon.fasta

# Ajuste de remapeamento
bowtie2-build ${F}_ivar.fa ${F}_ivar.fa

if [[ ! -s ${F}_ivar.fa ]]; then
    echo "Erro: Arquivo ${F}_ivar.fa não foi gerado."
    exit 1
fi

bowtie2 -p 1 -x ${F}_ivar.fa  -1 ${F}_R1_paired.fq.gz -2 ${F}_R2_paired.fq.gz -U ${F}_unpaired.fq.gz | samtools view -S -b -q 30 > ${F}.Pilon2.bam

# Ordenar e indexar o BAM
samtools sort ${F}.Pilon2.bam -o ${F}.Pilon2.sorted.bam
samtools index ${F}.Pilon2.sorted.bam

# Chamada de variantes com iVar
samtools mpileup -aa -A -d 0 -B -Q 0 ${F}.Pilon2.sorted.bam | ivar variants -p ${F}_iVar_variants -q 20 -t 0.25 -m 10 -r ${F}_ivar.fa
grep 'TRUE' ${F}_iVar_variants.tsv | cut -f 3,4 | grep -v '+' | grep -v '-' | grep -vP 'N\t' | grep -vP '\tN' | wc -l > ${F}.SNPsCount

# Remover bases degeneradas e gerar o Genoma fasta
python $PIPELINE/SARS-CoV-2/substitute_degenarate_bases.py ${F}_ivar.fa Genoma_${F}.fasta

# 5. Análise de Cobertura e Estatísticas

# Diretório onde estão os arquivos BED
BED_DIR="/storage/zuleika/volume1/project/carol/sarsCov2/CeVIVAS/pipeline"

# Copiar os arquivos BED para o diretório de trabalho
cp ${BED_DIR}/v3_check_*.bed ./
cp ${BED_DIR}/pcr_mix*.bed ./

# Obter o nome da referência do arquivo BAM
REF_NAME=$(samtools view -H ${F}.Pilon2.sorted.bam | awk -F'\t' '$1=="@SQ" {print $2}' | cut -d':' -f2)

# Se REF_NAME estiver vazio, exibir erro e parar o pipeline
if [[ -z "$REF_NAME" ]]; then
    echo "Erro: Não foi possível encontrar o nome da referência no BAM."
    exit 1
fi


# Atualizar os cabeçalhos dos arquivos BED
for BED in v3_check_*.bed pcr_mix*.bed; do
    sed -i "s/^NC_045512.2_pilon/${REF_NAME}/" "$BED"
done

echo "Arquivos BED corrigidos e prontos para análise!"

# Verificar cobertura em regiões chave
samtools depth -b v3_check_S.bed ${F}.Pilon2.sorted.bam | awk '{sum+=$3} END {print sum/NR}' > ${F}.S_MeanDepth;
samtools depth -b v3_check_ORF1ab.bed ${F}.Pilon2.sorted.bam | awk '{sum+=$3} END {print sum/NR}' > ${F}.ORF1ab_MeanDepth;
samtools depth -b v3_check_N.bed ${F}.Pilon2.sorted.bam | awk '{sum+=$3} END {print sum/NR}' > ${F}.N_MeanDepth;

# Verificar cobertura em regiões de PCR (MÉDIA)
samtools depth -b pcr_mix1.bed ${F}.Pilon2.sorted.bam | awk '{sum+=$3} END {print sum/NR}' > ${F}.PCR1_MeanDepth;
samtools depth -b pcr_mix2.bed ${F}.Pilon2.sorted.bam | awk '{sum+=$3} END {print sum/NR}' > ${F}.PCR2_MeanDepth;

# Verificar cobertura em regiões de PCR (MEDIANA)
samtools depth -b pcr_mix1.bed ${F}.Pilon2.sorted.bam | awk '{print $3}' | sort -n | awk 'NF{a[NR]=$1;c++}END {print (c%2==0)?(a[int(c/2)]+a[int(c/2)+1])/2:a[int(c/2)]}' > ${F}.PCR1_MedianDepth;
samtools depth -b pcr_mix2.bed ${F}.Pilon2.sorted.bam | awk '{print $3}' | sort -n | awk 'NF{a[NR]=$1;c++}END {print (c%2==0)?(a[int(c/2)]+a[int(c/2)+1])/2:a[int(c/2)]}' > ${F}.PCR2_MedianDepth;


# Plotar profundidade de cobertura
samtools depth -a ${F}.Pilon2.sorted.bam > ${F}.GenomeDepth
python $PIPELINE/SARS-CoV-2/plot_depth_per_sample.py ${F}.GenomeDepth ${F}_Depth_plot.png ${F}


# 6. Classificação de Linhagem com Pangolin
pangolin Genoma_${F}.fasta --outfile Pangolin__${F}.tsv
sed -i "s/,/\t/g" Pangolin__${F}.tsv

# 7. Geração de Estatísticas Finais
#zcat ${F}_R1.fq.gz ${F}_R2.fq.gz | echo $((`wc -l`/4)) > ${F}.ReadCount;
echo $(( $(wc -l < <(zcat ${F}_R1.fq.gz ${F}_R2.fq.gz)) / 4 )) > ${F}.ReadCount;

# Criar estatísticas básicas, assumindo que não há mapeamento
if [[ ! -s ${F}.Pilon2.sorted.bam ]]; then
    echo "0" > ${F}.ReadsMapped
    echo "0.00" > ${F}.PercentMapped
    echo "0.00" > ${F}.MeanDepth
    echo "0.00" > ${F}.MedianDepth
    echo "0" > ${F}.Depth10
    echo "0" > ${F}.Depth25
    echo "0.00" > ${F}.CoverageBlastn
    echo "0" > ${F}.CountNs
    echo "NA" > ${F}.Pangolin
    echo "NA" > ${F}.scorpio

    # Criar um genoma fictício de SARS-CoV-2 com 29903 Ns
    echo ">Genoma_${F}" > Genoma_${F}.fasta
    printf 'N%.0s' {1..29903} >> Genoma_${F}.fasta
    echo "" >> Genoma_${F}.fasta
else
    # Calcular estatísticas normalmente
    samtools view -c -F 260 ${F}.Pilon2.sorted.bam > ${F}.ReadsMapped
    x=$(cat ${F}.ReadsMapped); y=$(cat ${F}.ReadCount); python -c "print(round(float(${x}/${y}*100), 2))" > ${F}.PercentMapped
    samtools depth -a ${F}.Pilon2.sorted.bam | awk '{sum+=$3} END {print (NR>0) ? sum/NR : 0}' > ${F}.MeanDepth
    samtools depth -a ${F}.Pilon2.sorted.bam | awk '{print $3}' | sort -n | awk '{a[NR]=$1} END {if (NR>0) {print (NR%2==0)?(a[int(NR/2)]+a[int(NR/2)+1])/2:a[int(NR/2)+1]} else {print 0}}' > ${F}.MedianDepth
    samtools depth ${F}.Pilon2.sorted.bam  |  awk '{print $3 >= 10}' | grep '1' | wc -l > ${F}.Depth10
    samtools depth ${F}.Pilon2.sorted.bam  |  awk '{print $3 >= 25}' | grep '1' | wc -l > ${F}.Depth25
    makeblastdb -in $PIPELINE/SARS-CoV-2/Reference_Genomes/BA.2_ref.fasta -dbtype nucl -out $PIPELINE/SARS-CoV-2/Reference_Genomes/BA.2_ref
    blastn -query Genoma_${F}.fasta -db $PIPELINE/SARS-CoV-2/Reference_Genomes/BA.2_ref -outfmt 6 > ${F}.blastn
    cat ${F}.blastn | awk '{x = $8-$7; print x < 0 ? -x+1 : x+1}' | awk '{sum+=$1} END {coverage = sum/29903 * 100"%"; printf "%0.2f\n", coverage}' > ${F}.CoverageBlastn
    seqtk comp Genoma_${F}.fasta | awk '{x+=$9}END{print x}' > ${F}.CountNs
    cat Pangolin__${F}.tsv | awk NR==2 | awk '{print $2}' > ${F}.Pangolin
    cut -f 5  Pangolin__${F}.tsv | tail -n+2 > ${F}.scorpio
fi

# Gerar o arquivo GenomeName
ls Genoma_${F}.fasta > ${F}.GenomeName;

# Colocar o mesmo nome do aquivo no cabeçalho 
K=Genoma_${F}.fasta # New header
sed -i "1s/^.*$/>${K}/" Genoma_${F}.fasta;

# Processar variantes intrahost
eval "$(conda shell.bash hook)"
conda activate /storage/zuleika/volume1/project/carol/sarsCov2/conda_envs/alexranieri/intrahost

cat Genoma_${F}.fasta | sed 's/__/_/g' | sed 's/.fasta//g' > Genoma_${F}_minor_check.fasta
mHeader=$(grep '>' Genoma_${F}_minor_check.fasta | sed 's/>//g')
mafft --keeplength --add Genoma_${F}_minor_check.fasta $PIPELINE/SARS-CoV-2/Reference_Genomes/BA.2_ref.fasta > ${F}_ref.algn
bam-readcount -w 0 -f $PIPELINE/SARS-CoV-2/Reference_Genomes/BA.2_ref.fasta -q 30 ${F}_Wuhan_mapped.sorted.bam  > ${mHeader}.fa.bc
python $PIPELINE/SARS-CoV-2/intrahost_script.py -in ${mHeader}.fa.bc -dp 100 -al ${F}_ref.algn

conda deactivate

if test -f ${F}_ref.algn.minor.fa; then
    mv ${F}_ref.algn.minor.fa Genoma_${F}_minor.fasta
    pangolin -t 1 Genoma_${F}_minor.fasta --outfile Pangolin__${F}_minor.tsv
    sed -i "s/,/\t/g" Pangolin__${F}_minor.tsv
    cat Pangolin__${F}_minor.tsv | awk NR==2 | awk '{print $2}' > ${F}_pango_minor;
    echo 'Yes' > minor_check.txt
else
    echo 'No' > minor_check.txt
fi

# Imprimir estatísticas em um arquivo .Statistics
printf "Genome\tN_Reads\tReads_mapped\tPercent_mapped\tMean_depth\tMedian_depth\tNpos_Depth>=10\tNpos_Depth>=25\tCoverage\tNumber_of_Ns\tPangolin_lineage\tscorpio_call\tS_Mean_depth\tORF1ab_mean_depth\tN_mean_depth\tSNPs_count\tPCR1_mean_depth\tPCR2_mean_depth\tPCR1_median_depth\tPCR2_median_depth\tIntrahost\n" > ${F}.Statistics;

paste -d "\t" \
*GenomeName* *ReadC* *ReadsM* *Percent* \
*.MeanDepth *.MedianDepth *Depth10 *Depth25 *CoverageBlastn *CountN* \
*.Pangolin *.scorpio *.S_MeanDepth *.ORF1ab_MeanDepth *.N_MeanDepth *.SNPsCount \
*.PCR1_MeanDepth *.PCR2_MeanDepth *.PCR1_MedianDepth *.PCR2_MedianDepth minor_check.txt >> ${F}.Statistics;

# Copiar genoma e estatísticas para o diretório superior
cp Genoma_${F}.fasta ../;
cp ${F}.Statistics ../;

# 8. Limpeza e Finalização
# Remover arquivos temporários
#rm *.gz
#rm  ${F}.Pilon.bam ${F}.Pilon2.bam ${F}.Pilon.sorted.bam ${F}_Wuhan_mapped.sorted.bam ${F}_Wuhan_mapped.bam ${F}_minor_mapped.bam
#rm ${F}_Wuhan_mapped.paired.bam ${F}_Wuhan_mapped.paired.sorted.bam

# Voltar ao diretório superior
cd ..;
