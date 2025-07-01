#!/bin/bash

GENE_NAME="rpoB"
CONSENSO="rpoB_consensus.fasta"
GENOMAS="bacillus_genomes.fasta"
SALIDA="${GENE_NAME}_reference.fasta"

# Crear base BLAST temporal
makeblastdb -in "$GENOMAS" -dbtype nucl -out tempdb

# Ejecutar blastn local y guardar hits con ≥70% identidad
blastn -query "$CONSENSO" -db tempdb -outfmt "6 sseqid sstart send pident length" -perc_identity 70 -max_target_seqs 1 > hits.tsv

# Extraer secuencias con seqkit (requiere seqkit)
cut -f1 hits.tsv | sort | uniq | seqkit grep -f - "$GENOMAS" > "$SALIDA"

# Limpiar
rm tempdb.* hits.tsv

echo "Archivo generado: $SALIDA"
