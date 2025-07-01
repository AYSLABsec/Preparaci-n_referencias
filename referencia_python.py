#!/usr/bin/env python3

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# --------------------------
# FUNCIONES
# --------------------------

def porcentaje_identidad(seq1, seq2):
    min_len = min(len(seq1), len(seq2))
    matches = sum(a == b for a, b in zip(seq1[:min_len], seq2[:min_len]))
    return matches / min_len

def buscar_amplicones(fasta_input, consenso, gene_name, umbral_identidad=0.7):
    resultados = []
    consenso = consenso.upper()

    for record in SeqIO.parse(fasta_input, "fasta"):
        secuencia = str(record.seq).upper()
        for i in range(len(secuencia) - len(consenso) + 1):
            fragmento = secuencia[i:i+len(consenso)]
            identidad = porcentaje_identidad(consenso, fragmento)
            if identidad >= umbral_identidad:
                nuevo_record = record[:0]
                nuevo_record.seq = Seq(fragmento)
                resultados.append(nuevo_record)
                break  # solo una coincidencia por secuencia
    return resultados

# --------------------------
# INTERFAZ DE CONSOLA
# --------------------------

def main():
    print("🔍 BUSCADOR DE AMPLICONES POR IDENTIDAD\n")

    # Entrada de usuario
    fasta_input = input("📂 Nombre del archivo FASTA de entrada: ").strip()
    gene_name = input("🧬 Nombre del amplicón (gen): ").strip()
    consenso = input("🧬 Ingresa la secuencia consenso (una sola línea): ").strip()

    # Buscar amplicones
    print("\n⏳ Buscando coincidencias con ≥70% identidad...")
    resultados = buscar_amplicones(fasta_input, consenso, gene_name)

    # Guardar resultados
    salida = f"{gene_name}_reference.fasta"
    if resultados:
        SeqIO.write(resultados, salida, "fasta")
        print(f"✅ Se guardaron {len(resultados)} fragmentos en: {salida}")
    else:
        print("❌ No se encontraron secuencias con suficiente similitud.")

if __name__ == "__main__":
    main()
