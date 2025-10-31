#!/usr/bin/env bash

# Recebe dois argumentos:
#   1) Nome base da amostra (ex: ID da amostra, sem extensões)
#   2) Data (ou identificador de data/pasta)
G=$1;  # nome base da amostra
H=$2;  # data ou identificador de pasta
F=FLU__${G}__${H};  # Nome padrão para o diretório de saída e prefixos

# Armazena o diretório atual e navega até ele (garantindo que o script rode de onde foi chamado)
s=$(pwd); cd $s

# Cria um diretório para cada amostra, caso não exista.
mkdir -p ${F};

# Agrupa e concatena (cat) todos os arquivos R1*.fastq.gz e R2*.fastq.gz
# correspondentes à amostra para formar um único R1.fq.gz e R2.fq.gz.
cat ${G}*R1*.fastq.gz  > ${F}/${F}_R1.fq.gz
cat ${G}*R2*.fastq.gz  > ${F}/${F}_R2.fq.gz

# Entra no diretório da amostra.
cd ${F}
echo $F > sample_name.txt

# ========================================================================
# 1) Remoção de primers e adaptadores usando cutadapt
#    Utiliza um arquivo fasta com sequências de primers e adaptadores
# ========================================================================
cutadapt -b file:$PIPELINE/Influenza/primer_and_adapter_colection.fasta -B file:$PIPELINE/Influenza/primer_and_adapter_colection.fasta -j 10 -o ${F}_R1_cutadapt.fq.gz -p ${F}_R2_cutadapt.fq.gz ${F}_R1.fq.gz ${F}_R2.fq.gz

# ========================================================================
# 2) Limpeza de reads com trimmomatic
#    Remove bases de baixa qualidade e reads muito curtas
# ========================================================================
trimmomatic PE -phred33 ${F}_R1_cutadapt.fq.gz ${F}_R2_cutadapt.fq.gz ${F}_R1_paired.fq.gz ${F}_R1_unpaired.fq.gz ${F}_R2_paired.fq.gz ${F}_R2_unpaired.fq.gz LEADING:3 TRAILING:3 SLIDINGWINDOW:5:20 MINLEN:35 TOPHRED33 -threads 10

# Remove os arquivos intermediários de cutadapt para economizar espaço.
rm -rf ${F}_R1_cutadapt.fq.gz ${F}_R2_cutadapt.fq.gz

# Concatena os reads não pareados em um único arquivo .fq.gz
cat ${F}_R1_unpaired.fq.gz ${F}_R2_unpaired.fq.gz > ${F}_unpaired.fq.gz
rm -rf ${F}_R1_unpaired.fq.gz ${F}_R2_unpaired.fq.gz

# ========================================================================
# 3) Identificação de melhor referência para cada segmento usando vapor
#    "vapor.py" compara reads vs. um conjunto de referências
# ========================================================================
cat $PIPELINE/Influenza/segment_list.txt \
| xargs -P 8 -I {} sh -c 'N=$(cat sample_name.txt); /project/carol/influenza/pipeline/software/vapor/vapor.py -fa $PIPELINE/Influenza/segments_database/{}.fasta -fq "$N"_R1_paired.fq.gz "$N"_R2_paired.fq.gz | cut -f 6 | sed "s/>//g" > {}_result.txt'

# ========================================================================
# 4) Determina o tipo (A ou B), o H (H1, H3, Victoria, Yamagata, etc) e o N
#    usando os arquivos resultantes de vapor
# ========================================================================
if [ -s segment_7_result.txt ]; then
    cut -f 1 -d '_' segment_7_result.txt > A_or_B_genotype.txt
else
    echo "not_detected" > A_or_B_genotype.txt
fi

if [ -s segment_4_result.txt ]; then
    cut -f 1 -d '_' segment_4_result.txt > H_genotype.txt
else
    echo "not_detected" > H_genotype.txt
    echo "not_detected" > clade.txt
fi

if [ -s segment_6_result.txt ]; then
    cut -f 1 -d '_' segment_6_result.txt > N_genotype.txt
else
    echo "not_detected" > N_genotype.txt
fi

# ========================================================================
# 5) Extrai a sequência de referência escolhida (para cada segmento)
#    e gera arquivos segment_i_ref.fasta
# ========================================================================
cat $PIPELINE/Influenza/segment_list.txt | xargs -P 8 -I {} sh -c 'seqtk subseq $PIPELINE/Influenza/segments_database/{}.fasta {}_result.txt> {}_ref.fasta'

# ========================================================================
# 6) Mapeamento, montagem e refinamento de cada segmento
#    - Usa bowtie2 para mapear contra a referência
#    - Usa spades.py para montar
#    - Usa minimap2 + samtools + ivar para gerar consenso
#    - Faz remapeamentos sucessivos para refinamento
#    - Avalia a clade (Nextclade) se for H1, H3, Victoria ou Yamagata
#    - Gera estatísticas por segmento
# ========================================================================
mkdir segments
for i in 1 2 3 4 5 6 7 8; do

    # Verifica se o arquivo segment_i_result.txt existe e tem conteúdo
    # (ou seja, se o segmento foi detectado pelo vapor)
    if [ -s segment_${i}_result.txt ]; then

        # Constrói índice bowtie2 para a referência escolhida do segmento
        bowtie2-build segment_${i}_ref.fasta segment_${i}_ref.fasta

        # Mapeia reads paired contra a ref
        bowtie2 --very-sensitive -p 10 -x segment_${i}_ref.fasta -1 ${F}_R1_paired.fq.gz -2 ${F}_R2_paired.fq.gz | samtools view -S -b -q 20 > segment_${i}_mapped.bam

        # Ordena o BAM
        samtools sort segment_${i}_mapped.bam -o segment_${i}_mapped.sorted.bam
        rm -rf segment_${i}_mapped.bam
        samtools index segment_${i}_mapped.sorted.bam

        # Recupera apenas os pares mapeados
        samtools view  -u -f 1 -F 12 segment_${i}_mapped.sorted.bam > segment_${i}_mapped.paired.bam
        rm -rf segment_${i}_mapped.sorted.bam

        samtools sort -n -o segment_${i}_mapped.sorted.paired.bam segment_${i}_mapped.paired.bam
        rm -rf segment_${i}_mapped.paired.bam

        bedtools bamtofastq -i segment_${i}_mapped.sorted.paired.bam -fq segment_${i}_R1.fq -fq2 segment_${i}_R2.fq

        gzip segment_${i}_R1.fq segment_${i}_R2.fq

        # Montagem por SPAdes
        spades.py -1 segment_${i}_R1.fq.gz -2 segment_${i}_R2.fq.gz -o SPAdes_segment${i} -t 10 --only-assembler --careful --cov-cutoff 10.0

        # Se SPAdes gerar scaffolds.fasta, faz o refinamento
        if test -f SPAdes_segment${i}/scaffolds.fasta; then
            minimap2 -a segment_${i}_ref.fasta SPAdes_segment${i}/scaffolds.fasta | samtools view -Sb > minimap_segment_${i}.bam
            samtools sort -o minimap_segment_${i}.sorted.bam minimap_segment_${i}.bam
            samtools index minimap_segment_${i}.sorted.bam
            rm -rf minimap_segment_${i}.bam

            # Gera consenso com ivar
            samtools mpileup -aa -A -d 1000 -C 50 -Q 0 minimap_segment_${i}.sorted.bam | ivar consensus -p minimap_segment_${i}_preconsensus -m 1 -q 0 -t 0

            # Substitui bases degeneradas (script Python externo)
            python $PIPELINE/Influenza/substitute_degenarate_bases.py minimap_segment_${i}_preconsensus.fa minimap_segment_${i}_preconsensus_fixed.fa
            rm -rf minimap_segment_${i}_preconsensus.fa

            # Move arquivo final para segments/
            mv minimap_segment_${i}_preconsensus_fixed.fa segments/segment_preconsensus_${i}_${F}.fasta

            # Reconstrói índice e remapeia novamente
            bowtie2-build segments/segment_preconsensus_${i}_${F}.fasta segments/segment_preconsensus_${i}_${F}.fasta
            bowtie2 --very-sensitive -p 10 -x segments/segment_preconsensus_${i}_${F}.fasta -1 segment_${i}_R1.fq.gz -2 segment_${i}_R2.fq.gz | samtools view -S -b -q 20 > segments/segment_preconsensus_${i}_${F}.bam

            samtools sort segments/segment_preconsensus_${i}_${F}.bam -o segments/segment_preconsensus_${i}_${F}.sorted.bam
            rm -rf segments/segment_preconsensus_${i}_${F}.bam segments/segment_preconsensus_${i}_${F}.fasta

            # Novo consenso (ivar)
            samtools mpileup -aa -A -d 0 -Q 0 segments/segment_preconsensus_${i}_${F}.sorted.bam | ivar consensus -p segment_${i}_${F}_ivar -i segment${i}_${F}.fasta

            # Substitui bases degeneradas novamente
            python $PIPELINE/Influenza/substitute_degenarate_bases.py segment_${i}_${F}_ivar.fa segments/segment_${i}_${F}.fasta
            rm -rf segment_${i}_${F}_ivar.fa
            rm -rf segments/segment_preconsensus_${i}_${F}.sorted.bam

            # Se o segmento for 4 (HA), faz análise de clade via nextclade
            if [ $(cat H_genotype.txt) == "H1" ] && [ $i -eq 4 ]; then
                nextclade run -D $PIPELINE/Influenza/nextclade_files/H1 -j 10 -t nextclade.tsv segments/segment_${i}_${F}.fasta
                csvcut -t -c clade nextclade.tsv | tail -n+2 > clade.txt
            elif [ $(cat H_genotype.txt) == "H3" ] && [ $i -eq 4 ]; then
                nextclade run -D $PIPELINE/Influenza/nextclade_files/H3 -j 10 -t nextclade.tsv segments/segment_${i}_${F}.fasta
                csvcut -t -c clade nextclade.tsv | tail -n+2 > clade.txt
            elif [ $(cat H_genotype.txt) == "Victoria" ] && [ $i -eq 4 ]; then
                nextclade run -D $PIPELINE/Influenza/nextclade_files/Vic -j 10 -t nextclade.tsv segments/segment_${i}_${F}.fasta
                csvcut -t -c clade nextclade.tsv | tail -n+2 > clade.txt
            elif [ $(cat H_genotype.txt) == "Yamagata" ] && [ $i -eq 4 ]; then
                nextclade run -D $PIPELINE/Influenza/nextclade_files/Yam -j 10 -t nextclade.tsv segments/segment_${i}_${F}.fasta
                csvcut -t -c clade nextclade.tsv | tail -n+2 > clade.txt
            elif [ $i -eq 4 ]; then
                echo "nextclade_not_available" > clade.txt
            fi

            # Faz um novo remapeamento final (ajuste do remapeamento)
            bowtie2-build segments/segment_${i}_${F}.fasta segments/segment_${i}_${F}.fasta
            bowtie2 --very-sensitive -p 10 -x segments/segment_${i}_${F}.fasta -1 segment_${i}_R1.fq.gz -2 segment_${i}_R2.fq.gz | samtools view -S -b -q 20 > segments/segment_${i}_${F}.bam

            samtools sort segments/segment_${i}_${F}.bam -o segments/segment_${i}_${F}.sorted.bam
            rm -rf segments/segment_${i}_${F}.bam

            # Usa ivar variants para identificar SNPs
            samtools mpileup -aa -A -d 0 -B -Q 0 segments/segment_${i}_${F}.sorted.bam | ivar variants -p segments/segment_${i}_${F}_iVar_variants -t 0.25 -m 10 -r segments/segment_${i}_${F}.fasta

            # Gera estatísticas de profundidade, etc.
            grep 'TRUE' segments/segment_${i}_${F}_iVar_variants.tsv | grep -vP 'N\t' | grep -vP '\tN' | wc -l > segments/segment_${i}.SNPsCount

            samtools depth -a segments/segment_${i}_${F}.sorted.bam | awk '{sum+=$3} END {print sum/NR}' > segments/segment_${i}.MeanDepth

            samtools depth -a segments/segment_${i}_${F}.sorted.bam | awk '{print $3}' | sort -n | awk 'NF{a[NR]=$1;c++}END {print (c%2==0)?(a[int(c/2)+1]+a[int(c/2)])/2:a[int(c/2)+1]}' > segments/segment_${i}.MedianDepth

            samtools depth -a segments/segment_${i}_${F}.sorted.bam | awk '{print $3 >= 10}' | grep '1' | wc -l > segments/segment_${i}.Depth10

            samtools depth -a segments/segment_${i}_${F}.sorted.bam | awk '{print $3 >= 25}' | grep '1' | wc -l > segments/segment_${i}.Depth25

            seqtk comp segments/segment_${i}_${F}.fasta | awk '{x+=$9}END{print x}' > segments/segment_${i}.CountNs

            total_segment=$(seqtk comp segments/segment_${i}_${F}.fasta | awk '{print $1 "\t" ($3+$4+$5+$6)}' |  cut -f2)
            total_ref=$(seqtk comp segment_${i}_ref.fasta | cut -f2)

            echo "scale=2; ($total_segment / $total_ref) * 100" | bc > segments/segment_${i}.coverage

            # Cria cabeçalho e junta estatísticas num arquivo .Statistics
            printf "segment_${i}_Mean_depth\tsegment_${i}_Median_depth\tsegment_${i}_Npos_Depth>=10\tsegment_${i}_Npos_Depth>=25\tsegment_${i}_Coverage\tsegment_${i}_Number_of_Ns\tsegment_${i}_SNPs\n" > segments/segment_${i}.Statistics

            paste -d "\t" segments/segment_${i}.MeanDepth segments/segment_${i}.MedianDepth segments/segment_${i}.Depth10 segments/segment_${i}.Depth25 segments/segment_${i}.coverage segments/segment_${i}.CountNs segments/segment_${i}.SNPsCount >> segments/segment_${i}.Statistics

            rm -rf segment_${i}_R1.fq.gz segment_${i}_R2.fq.gz

        else
            # Caso não seja possível montar (scaffolds.fasta ausente)
            printf "segment_${i}_Mean_depth\tsegment_${i}_Median_depth\tsegment_${i}_Npos_Depth>=10\tsegment_${i}_Npos_Depth>=25\tsegment_${i}_Coverage\tsegment_${i}_Number_of_Ns\tsegment_${i}_SNPs\n" > segments/segment_${i}.Statistics

            printf "could_not_be_assembled\tcould_not_be_assembled\tcould_not_be_assembled\tcould_not_be_assembled\tcould_not_be_assembled\tcould_not_be_assembled\tcould_not_be_assembled" >> segments/segment_${i}.Statistics

            # Se for o segmento 4 e não montou, define clade como não detectado
            if [ $i -eq 4 ]; then
                echo "not_detected" > clade.txt
            fi
        fi

    else
        # Caso o segmento não tenha sido detectado pelo vapor (arquivo vazio ou inexistente)
        printf "segment_${i}_Mean_depth\tsegment_${i}_Median_depth\tsegment_${i}_Npos_Depth>=10\tsegment_${i}_Npos_Depth>=25\tsegment_${i}_Coverage\tsegment_${i}_Number_of_Ns\tsegment_${i}_SNPs\n" > segments/segment_${i}.Statistics

        printf "not_detected\tnot_detected\tnot_detected\tnot_detected\tnot_detected\tnot_detected\tnot_detected" >> segments/segment_${i}.Statistics
    fi
done

# Junta as estatísticas de todos os segmentos em um único arquivo
paste -d "\t" segments/*.Statistics > segments.Statistics

# ========================================================================
# 7) Estatísticas gerais para o genoma inteiro
#    - Cria Genoma_${F}.fasta com todos os segmentos montados
#    - Mapeia novamente para estatísticas globais
# ========================================================================
cat segments/*.fasta > Genoma_${F}.fasta

if [ -s Genoma_${F}.fasta ]; then
    bowtie2-build Genoma_${F}.fasta Genoma_${F}.fasta
    bowtie2 --very-sensitive -p 10 -x Genoma_${F}.fasta -1 ${F}_R1_paired.fq.gz -2 ${F}_R2_paired.fq.gz | samtools view -S -b -q 20 > ${F}.bam

    samtools sort ${F}.bam -o ${F}.sorted.bam
    samtools index ${F}.sorted.bam
    rm -rf ${F}.bam

    # Conta quantos reads mapearam
    samtools view -c -F 260 ${F}.sorted.bam > ${F}.ReadsMappedFinal

    # Conta total de reads de entrada (paired)
    echo $(zcat ${F}_R1_paired.fq.gz ${F}_R2_paired.fq.gz | wc -l)/4 | bc > ${F}.ReadCount

    # Calcula % de reads mapeados
    x=$(cat ${F}.ReadsMappedFinal)
    y=$(cat ${F}.ReadCount)
    python -c "print(round(float(${x}/${y}*100), 2))" > ${F}.PercentMapped

    # Conta quantos segmentos foram montados
    grep -c '>' Genoma_${F}.fasta > ${F}.SegmentsAssembled

    # Gera um arquivo 'subtype.txt' unindo H e N
    paste H_genotype.txt N_genotype.txt | sed 's/\t//g' > subtype.txt
else
    # Caso não tenha conseguido gerar nenhum segmento
    echo $(zcat ${F}_R1_paired.fq.gz ${F}_R2_paired.fq.gz | wc -l)/4 | bc > ${F}.ReadCount
    echo "0" > ${F}.ReadsMappedFinal
    echo "0" > ${F}.PercentMapped
    echo "0" > ${F}.SegmentsAssembled
    echo "not_detected" > subtype.txt
fi

# Cria um arquivo Genome_basic.Statistics com colunas básicas
ls Genoma_${F}.fasta > ${F}.GenomeName;
printf "Genome\tN_Reads\tReads_assembled\tPercent_assembled\tSegments_Assembled\tType\tH_genotype\tN_genotype\tSubtype\tClade\n" > Genome_basic.Statistics

paste -d "\t" ${F}.GenomeName ${F}.ReadCount ${F}.ReadsMappedFinal ${F}.PercentMapped ${F}.SegmentsAssembled A_or_B_genotype.txt H_genotype.txt N_genotype.txt  subtype.txt clade.txt >> Genome_basic.Statistics

# Junta Genome_basic.Statistics com segments.Statistics
paste -d "\t" Genome_basic.Statistics segments.Statistics > ${F}_complete.Statistics

# ========================================================================
# 9) Voltando para o diretório pai
# ========================================================================

# Copia o genoma montado e as estatísticas finais para o diretório pai
cp Genoma_${F}.fasta ../;
cp ${F}_complete.Statistics ../;

# Remove arquivos .fq.gz para liberar espaço
rm -rf *.fq.gz

# Retorna ao diretório pai
find . -type d | xargs -i chmod g+w {}
cd ..
