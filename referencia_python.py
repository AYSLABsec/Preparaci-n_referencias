#!/usr/bin/env python3

from Bio import SeqIO
from Bio.Seq import Seq
import argparse
import os

def porcentaje_identidad(seq1, seq2):
    """Calcula el porcentaje de identidad simple (por alineamiento local básico)"""
    l1, l2 = len(seq1), len(seq2)
    min_len = min(l1, l2)
    matches = sum(a == b for a, b in zip(seq1[:min_len], seq2[:min_len]))
    return matches / min_len

def buscar_similares(consenso, fasta_input, umbral=0.70):
    resultados = []
    for record in SeqIO.parse(fasta_input, "fasta"):
        secuencia = str(record.seq).upper()
        encontrado = False
        for i in range(len(secuencia) - len(consenso) + 1):
            fragmento = secuencia[i:i+len(consenso)]
            identidad = porcentaje_identidad(consenso, fragmento)
            if identidad >= umbral:
                nuevo_record = record[:0]  # copia del encabezado
                nuevo_record.seq = Seq(fragmento)
                resultados.append(nuevo_record)
                encontrado = True
                break  # solo toma la primera coincidencia
        if not encontrado:
            continue
    return resultados

def main():
    parser = argparse.ArgumentParser(description="Extrae fragmentos similares a una secuencia consenso")
    parser.add_argument("-i", "--input", required=True, help="Archivo FASTA con las secuencias (genomas/CDS)")
    parser.add_argument("-s", "--seq", required=True, help="Secuencia consenso (texto o archivo)")
    parser.add_argument("-g", "--gene", required=True, help="Nombre del gen para nombrar el archivo de salida")
    parser.add_argument("-t", "--threshold", type=float, default=0.7, help="Umbral de identidad (default: 0.7)")
    args = parser.parse_args()

    if os.path.isfile(args.seq):
        with open(args.seq) as f:
            consenso = f.read().strip().upper()
    else:
        consenso = args.seq.strip().upper()

    resultados = buscar_similares(consenso, args.input, umbral=args.threshold)

    salida = f"{args.gene}_reference.fasta"
    if resultados:
        SeqIO.write(resultados, salida, "fasta")
        print(f"Se guardaron {len(resultados)} secuencias en {salida}")
    else:
        print("No se encontraron secuencias con suficiente similitud.")

if __name__ == "__main__":
    main()
