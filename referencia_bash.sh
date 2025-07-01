#!/bin/bash

echo "🔍 BUSCADOR DE AMPLICONES POR IDENTIDAD (70%)"

# 1. Pedir inputs al usuario
read -p "📂 Nombre del archivo FASTA de entrada: " INPUT_FASTA
read -p "🧬 Nombre del amplicón (gen): " GENE_NAME
read -p "🧬 Ingresa la secuencia consenso (una sola línea): " CONSENSO

# 2. Crear script temporal en Python
TEMP_SCRIPT=$(mktemp)

cat << EOF > "$TEMP_SCRIPT"
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

fasta_input = "$INPUT_FASTA"
gene_name = "$GENE_NAME"
consenso = "$CONSENSO".upper()
umbral_identidad = 0.7

def porcentaje_identidad(seq1, seq2):
    min_len = min(len(seq1), len(seq2))
    matches = sum(a == b for a, b in zip(seq1[:min_len], seq2[:min_len]))
    return matches / min_len

resultados = []
for record in SeqIO.parse(fasta_input, "fasta"):
    secuencia = str(record.seq).upper()
    for i in range(len(secuencia) - len(consenso) + 1):
        fragmento = secuencia[i:i+len(consenso)]
        identidad = porcentaje_identidad(consenso, fragmento)
        if identidad >= umbral_identidad:
            nuevo_record = record[:0]
            nuevo_record.seq = Seq(fragmento)
            resultados.append(nuevo_record)
            break

salida = f"{gene_name}_reference.fasta"
if resultados:
    SeqIO.write(resultados, salida, "fasta")
    print(f"✅ Se guardaron {len(resultados)} secuencias en '{salida}'")
else:
    print("❌ No se encontraron secuencias con ≥70% identidad.")
EOF

# 3. Ejecutar el script Python
python3 "$TEMP_SCRIPT"

# 4. Limpiar
rm "$TEMP_SCRIPT"
