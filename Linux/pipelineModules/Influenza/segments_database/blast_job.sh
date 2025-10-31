#!/bin/bash
#SBATCH --job-name=blastn_parallel
#SBATCH --output=blastn_parallel_%j.out
#SBATCH --error=blastn_parallel_%j.err
#SBATCH --cpus-per-task=8
#SBATCH --mem=40G

echo "Iniciando BLASTN paralelo para todos os segmentos"
date

# Função para processar um segmento
process_segment() {
    segment=$1
    if [ -f "$segment" ]; then
        echo "Processando: $segment (PID: $$)"
        output_file="blast_results_${segment%.fasta}.txt"
        
        blastn -query $segment \
               -db /database/ncbi-blast+/nt_viruses \
               -out $output_file \
               -outfmt "6 qseqid sseqid pident qcovs evalue bitscore" \
               -max_target_seqs 4 \
               -evalue 1e-10
        
        echo "Concluído: $segment"
    else
        echo "AVISO: Arquivo $segment não encontrado!"
    fi
}

# Exportar função para subprocessos
export -f process_segment

# Processar todos os segmentos em paralelo
for segment in segment_{1..8}.fasta; do
    process_segment $segment &
done

# Aguardar todos terminarem
wait

echo "==============================================="
echo "TODOS OS SEGMENTOS CONCLUÍDOS!"
echo "Arquivos de resultado criados:"
ls -la blast_results_segment_*.txt
date
