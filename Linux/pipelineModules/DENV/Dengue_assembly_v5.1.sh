#!/bin/bash

# Parâmetros de entrada:
# $1 = Nome base da amostra (ex: BR-01)
# $2 = Data no formato AAAAMMDD (ex: 20231015)
G=$1;  # Armazena o primeiro argumento como nome base da amostra
H=$2;  # Armazena o segundo argumento como data
F=DENV__${G}__${H}; # Cria nome padronizado para arquivos (ex: DENV__BR-01__20231015)


# Configuração inicial do diretório
s=$(pwd); cd $s # Armazena e mantém o diretório atual (evita mudanças acidentais)

# Cria pasta para análise se não existir (-p evita erros se já existir)
mkdir -p ${F};

# Processamento inicial das reads
# Concatena todos arquivos R1 (forward) e R2 (reverse) em arquivos únicos
cat ${G}*R1*.fastq.gz  > ${F}/${F}_R1.fq.gz
cat ${G}*R2*.fastq.gz  > ${F}/${F}_R2.fq.gz

# Entra na pasta da análise
cd ${F}

########################################
# Etapa 1: Remoção de primers e adaptadores
########################################
# Arquivo com sequências de adaptadores/primer (forward)
# Arquivo com sequências de adaptadores/primer (reverse)
# Usa 1 núcleo de processamento
cutadapt -b file:/project/carol/dengue/pipeline/primer_and_adapter_colection_dengue.fasta -B file:/project/carol/dengue/pipeline/primer_and_adapter_colection_dengue.fasta -j 1 -o ${F}_R1_cutadapt.fq.gz -p ${F}_R2_cutadapt.fq.gz ${F}_R1.fq.gz ${F}_R2.fq.gz

########################################
# Etapa 2: Controle de qualidade das reads
########################################
# trimmomatic PE Modo Paired-End
# Encoding da qualidade
# LEADING:3 Remove bases com qualidade <3 no início
# TRAILING:3 Remove bases com qualidade <3 no final
# SLIDINGWINDOW:5:20 Janela deslizante de 5 bases, qualidade média mínima 20
# MINLEN:35 Descarta reads com menos de 35 bases
# TOPHRED33 Mantém encoding Phred33
#-threads 1 Usa 1 thread
trimmomatic PE -phred33 ${F}_R1_cutadapt.fq.gz ${F}_R2_cutadapt.fq.gz ${F}_R1_paired.fq.gz ${F}_R1_unpaired.fq.gz ${F}_R2_paired.fq.gz ${F}_R2_unpaired.fq.gz LEADING:3 TRAILING:3 SLIDINGWINDOW:5:20 MINLEN:35 TOPHRED33 -threads 10

# Remove arquivos intermediários
rm -rf ${F}_R1_cutadapt.fq.gz ${F}_R2_cutadapt.fq.gz

########################################
# Etapa 3: Mapeamento contra genomas de referência dos 4 sorotipos
########################################
sorotipo=(1 2 3 4) # Array com os 4 sorotipos de dengue

for i in "${sorotipo[@]}" 
do
    # Mapeamento com BWA MEM
    bwa mem -t 1 /project/carol/dengue/pipeline/reference_genomes/dengue_${i}.fasta ${F}_R1_paired.fq.gz ${F}_R2_paired.fq.gz > ${F}_dengue_${i}_mapped.sam
    
    # Filtra reads mapeadas com qualidade mínima (MAPQ >=10)
    samtools view -S -b -q 10 ${F}_dengue_${i}_mapped.sam > ${F}_dengue_${i}_mapped.bam
    
    # Ordena e indexa o arquivo BAM para análises posteriores
    samtools sort ${F}_dengue_${i}_mapped.bam -o ${F}_dengue_${i}_mapped.sorted.bam
    samtools index ${F}_dengue_${i}_mapped.sorted.bam
    
    # Remove arquivo SAM intermediário
    rm ${F}_dengue_${i}_mapped.sam
done

########################################
# Etapa 4: Seleção do sorotipo com maior mapeamento
########################################
# Conta reads mapeadas para cada sorotipo
for i in "${sorotipo[@]}" 
do
    samtools view -c -F 260 ${F}_dengue_${i}_mapped.sorted.bam > ${F}_${i}.ReadsMapped
done 

# Combina resultados e seleciona sorotipo com maior contagem
awk '{print $0"\t"FILENAME}' *.ReadsMapped > stats_ReadsMapped.tsv
selected_serotype=$(sort -k1 -rn stats_ReadsMapped.tsv | head -n1 | cut -f 2 | tail -c 14 | head -c 1)

########################################
# Etapa 5: Refinamento do genoma com Pilon
########################################
java -jar /usr/local/src/pilon/1.24/pilon-1.24.jar --genome /project/carol/dengue/pipeline/reference_genomes/dengue_${selected_serotype}.fasta --frags ${F}_dengue_${selected_serotype}_mapped.sorted.bam --output ${F}_${selected_serotype}.Pilon --fix "gaps,indels" --threads 1 --mindepth 5 --minqual 20 --minmq 10

########################################
# Etapa 6: Remapeamento contra o genoma corrigido
########################################
# Indexa novo genoma e remapeia
bwa index ${F}_${selected_serotype}.Pilon.fasta
bwa mem -t 1 ${F}_${selected_serotype}.Pilon.fasta ${F}_R1_paired.fq.gz ${F}_R2_paired.fq.gz > ${F}_dengue_${selected_serotype}_Pilon_mapped.sam

# Processamento do SAM para BAM ordenado
samtools view -S -b -q 10 ${F}_dengue_${selected_serotype}_Pilon_mapped.sam > ${F}_dengue_${selected_serotype}_Pilon_mapped.bam
samtools sort ${F}_dengue_${selected_serotype}_Pilon_mapped.bam -o ${F}_dengue_${selected_serotype}_Pilon_mapped.sorted.bam
samtools index ${F}_dengue_${selected_serotype}_Pilon_mapped.sorted.bam
rm ${F}_dengue_${selected_serotype}_Pilon_mapped.sam

########################################
# Etapa 7: Geração de consenso com iVar
########################################
# -m 5 Frequência mínima 5% para chamar variante
# -aa Gera pilha para todas as posições do genoma, mesmo aquelas sem cobertura (importante para manter o tamanho do genoma completo).
# -A Inclui reads marcados como "não primários" (por exemplo, secundários, duplicatas). Pode ser útil em genomas virais onde não se quer perder cobertura.
# -d 0 Remove o limite de profundidade de cobertura (por padrão o samtools usa limite de 8000). Aqui você deixa o programa usar toda a cobertura disponível.
# -Q 0 Inclui todas as bases, independentemente da qualidade de base. O filtro de qualidade será aplicado depois, pelo iVar.

samtools mpileup -aa -A -d 0 -Q 0 ${F}_dengue_${selected_serotype}_Pilon_mapped.sorted.bam | ivar consensus -m 5 -p ${F}_ivar -i ${F}_pilon

########################################
# Etapa 8: Ajuste final de mapeamento
########################################
bwa index ${F}_ivar.fa # Indexa genoma consenso
bwa mem ${F}_ivar.fa ${F}_R1_paired.fq.gz ${F}_R2_paired.fq.gz -t 1 | samtools view -S -b -q 10 > ${F}.Pilon2.bam # Filtra por qualidade
samtools sort -o ${F}.Pilon2.sorted.bam ${F}.Pilon2.bam # Ordena
samtools index ${F}.Pilon2.sorted.bam # Indexa

########################################
# Etapa 9: Chamada de variantes
########################################
samtools mpileup -aa -A -d 0 -B -Q 0 ${F}.Pilon2.sorted.bam | ivar variants -p ${F}_iVar_variants -q 20 -t 0.25 -m 5 -r ${F}_ivar.fa # Parâmetros de filtro

# Conta SNPs significativos
grep 'TRUE' ${F}_iVar_variants.tsv | grep -vP 'N\t' | grep -vP '\tN' | wc -l > ${F}.SNPsCount

# Substitui bases degeneradas no genoma consenso
python /project/carol/dengue/pipeline/substitute_degenarate_bases.py ${F}_ivar.fa Genoma_${F}.fasta
rm -rf ${F}_ivar.fa

########################################
# Etapa 10: Cálculo de estatísticas
########################################
# Contagem total de reads
echo $(zcat ${F}_R1.fq.gz ${F}_R2.fq.gz | wc -l)/4|bc > ${F}.ReadCount

# Porcentagem de reads mapeadas
samtools view -c -F 260 ${F}.Pilon2.sorted.bam > ${F}.ReadsMappedFinal
x=$(cat ${F}.ReadsMappedFinal); y=$(cat ${F}.ReadCount)
python -c "print(round(float(${x}/${y}*100), 2))" > ${F}.PercentMapped

# Profundidade média e mediana
samtools depth -a ${F}.Pilon2.sorted.bam | awk '{sum+=$3} END {print sum/NR}' > ${F}.MeanDepth
samtools depth -a ${F}.Pilon2.sorted.bam | awk '{print $3}' | sort -n | awk 'NF {a[++i]=$1} END {if (i%2==1) print a[(i+1)/2]; else print (a[i/2] + a[i/2+1])/2}' > ${F}.MedianDepth

# Cobertura em diferentes profundidades
samtools depth ${F}.Pilon2.sorted.bam | awk '$3 >= 10 {count++} END {print count}' > ${F}.Depth10
samtools depth ${F}.Pilon2.sorted.bam | awk '$3 >= 25 {count++} END {print count}' > ${F}.Depth25

# Cobertura geral do genoma
bases=$(seqtk comp Genoma_${F}.fasta | awk '{print $1 "\t" ($3+$4+$5+$6)}' | cut -f 2)
total=$(seqtk comp Genoma_${F}.fasta | cut -f2)
echo "scale=4; ($bases / $total) * 100" | bc > ${F}.CoverageBlastn

# Pegar o Serotype
echo $selected_serotype > ${F}.Serotype 

########################################
# Etapa 11: Análise específica por gene
########################################
# Extrai o nome do contig do cabeçalho do BAM (campo SN)
CONTIG=$(samtools view -H ${F}.Pilon2.sorted.bam | grep "^@SQ" | sed -E 's/.*SN:([^[:space:]]+).*/\1/' | head -n 1)
echo "Contig extraído: ${CONTIG}"

# Cria um FASTA temporário com o nome do contig atualizado
sed "s/>.*/>${CONTIG}/" Genoma_${F}.fasta > Genoma_${F}.temp.fasta

# Define genes e diretório dos BEDs
GENES=("E" "NS1" "NS3" "NS5")
BED_DIR="/project/carol/dengue/pipeline/gene_check"

# Decide prefixo do BED com base no sorotipo
case ${selected_serotype} in
    1) PREFIX="denv1" ;;
    2) PREFIX="denv2" ;;
    3) PREFIX="denv3" ;;
    4) PREFIX="denv4" ;;
    *) echo "Sorotipo inválido: ${selected_serotype}"; exit 1 ;;
esac

# Loop por genes
for gene in "${GENES[@]}"
do
    # Copia BED e ajusta contig
    sed "s/gene/${F}_pilon/g" ${BED_DIR}/${PREFIX}_${gene}.bed > ${PREFIX}_${gene}_tmp.bed

    # Profundidade média
    samtools depth -a ${F}.Pilon2.sorted.bam -b ${PREFIX}_${gene}_tmp.bed | \
        awk '{sum+=$3} END {print sum/NR}' > ${F}.Coverage_${gene}_MeanDepth

    # Porcentagem de bases com cobertura (coverage)
    bases=$(seqtk comp -r ${PREFIX}_${gene}_tmp.bed Genoma_${F}.fasta | awk '{print $1 "\t" ($7+$4+$5+$6)}' | cut -f2)
    total=$(seqtk comp -r ${PREFIX}_${gene}_tmp.bed Genoma_${F}.fasta | awk '{print $1 "\t" ($3-$2)}' | cut -f2)
    echo "scale=4; ($bases / $total) * 100" | bc > ${F}.Coverage_gene_${gene}
done

# Cobertura global com blastn (opcional, se quiser manter)
# cat ${F}.blastn | awk '{x = $8-$7; print x < 0 ? -x+1 : x+1}' | \
#     awk -v L=${GENOME_LEN} '{sum+=$1} END {printf "%.2f\n", (sum/L)*100}' > ${F}.CoverageBlastn

# Corrige o cabeçalho do arquivo FASTA
sed -i "s/^>.*/>Genoma_${F}.fasta/" Genoma_${F}.fasta

# Pegar a contagem de Ns
seqtk comp Genoma_${F}.fasta | awk '{x+=$9}END{print x}' > ${F}.Number_of_Ns

########################################
# Etapa 12: Classificação filogenética
########################################
nextclade run -D /project/carol/dengue/pipeline/nextclade_files2/denv${selected_serotype}/2024-08-31--20-44-06Z/ -j 1 -t nextclade.tsv Genoma_${F}.fasta # Arquivo de saída

# Extrai clado atribuído
csvcut -t -c clade nextclade.tsv | tail -n+2 > ${F}.clade

########################################
# Etapa Final: Consolidação de resultados
########################################
# Gerar nomes consistentes
mv ${F}.CoverageBlastn ${F}.Coverage

# Cabeçalho
printf "Genome\tN_Reads\tReads_mapped\tPercent_mapped\tMean_depth\tMedian_depth\tNpos_Depth>=10\tNpos_Depth>=25\tCoverage\tE_Coverage\tE_Depth\tNS1_Coverage\tNS1_Depth\tNS3_Coverage\tNS3_Depth\tNS5_Coverage\tNS5_Depth\tNumber_of_Ns\tSNPs\tSerotype\tclade\n" > ${F}.Statistics

# Corpo (usando nomes explícitos)
printf "Genoma_${F}.fasta\t$(cat ${F}.ReadCount)\t$(cat ${F}.ReadsMappedFinal)\t$(cat ${F}.PercentMapped)\t$(cat ${F}.MeanDepth)\t$(cat ${F}.MedianDepth)\t$(cat ${F}.Depth10)\t$(cat ${F}.Depth25)\t$(cat ${F}.Coverage)\t$(cat ${F}.Coverage_gene_E)\t$(cat ${F}.Coverage_E_MeanDepth)\t$(cat ${F}.Coverage_gene_NS1)\t$(cat ${F}.Coverage_NS1_MeanDepth)\t$(cat ${F}.Coverage_gene_NS3)\t$(cat ${F}.Coverage_NS3_MeanDepth)\t$(cat ${F}.Coverage_gene_NS5)\t$(cat ${F}.Coverage_NS5_MeanDepth)\t$(cat ${F}.Number_of_Ns)\t$(cat ${F}.SNPsCount)\t$(cat ${F}.Serotype)\t$(cat ${F}.clade)\n" >> ${F}.Statistics


# Move resultados para diretório principal
cp Genoma_${F}.fasta ../;
cp ${F}.Statistics ../;

# Limpeza opcional (descomentar se necessário)
# rm *.gz
# rm  *_Pilon_mapped.sorted.bam ${F}.Pilon2.bam

cd ..  # Retorna ao diretório principal
