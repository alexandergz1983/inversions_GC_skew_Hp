#!/usr/bin/env python
# -*- coding: utf-8 -*-
# CREATED BY: Alexander Garcia-Zea
# jagarcia@eafit.edu.co

import os
import argparse
import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import subprocess
import logging

# Configurar logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def process_clusters(gbk_dir, cog_file, cdho_file, output_dir, csv_output):
    # Leer los archivos
    cog_df = pd.read_csv(cog_file)
    cdho_df = pd.read_csv(cdho_file, sep='\t')

    # Lista para el resumen
    gene_summary = []

    # Procesar cada fila de la columna Gene en el archivo --cog
    for _, row in cog_df.iterrows():
        group_name = row['Gene']  # Nombre del grupo (columna Gene)
        genome_annotations = row.iloc[3:]  # Columnas después de 'Annotation'
        logging.info(f"Procesando grupo: {group_name}")

        # Diccionarios para almacenar secuencias por orientación
        sequences = {'CD_nt': [], 'HO_nt': [], 'CD_aa': [], 'HO_aa': []}
        cd_count = 0
        ho_count = 0

        # Procesar cada columna genoma
        for genome, annotations in genome_annotations.items():
            if pd.isna(annotations):
                continue  # Saltar si no hay anotaciones en esta celda

            # Separar múltiples anotaciones por ';'
            for annotation in map(str.strip, annotations.split(';')):
                # Buscar la anotación en --cdho
                gene_info = cdho_df[cdho_df['FEATURE_NAME'].str.lower() == annotation.lower()]

                if gene_info.empty:
                    logging.warning(f"Anotación {annotation} no encontrada en CDHO.")
                    continue

                # Procesar cada fila correspondiente en CDHO
                for _, gene_row in gene_info.iterrows():
                    seq_name = gene_row['SEQ_NAME']
                    start = int(gene_row['START'])
                    end = int(gene_row['END'])
                    strand = gene_row['STRAND']
                    orientation = 'CD' if strand == '+' else 'HO'

                    logging.info(f"Anotación {annotation} - SEQ_NAME: {seq_name} - Coordenadas: {start}-{end} - Strand: {strand}")

                    # Localizar el archivo GenBank correspondiente
                    gbk_file = f"{seq_name}.gbk"
                    gbk_path = os.path.join(gbk_dir, gbk_file)

                    if not os.path.exists(gbk_path):
                        logging.warning(f"Archivo {gbk_path} no encontrado.")
                        continue

                    # Leer y procesar el archivo GenBank
                    for record in SeqIO.parse(gbk_path, 'genbank'):
                        # Extraer la secuencia nucleotídica
                        nucleotide_seq = record.seq[start - 1:end] if start < end else record.seq[end - 1:start].reverse_complement()
                        seq_id = f"{seq_name}_{annotation}"
                        seq_record_nt = SeqRecord(nucleotide_seq, id=seq_id, description="")

                        # Agregar la secuencia nucleotídica al diccionario correspondiente
                        sequences[f"{orientation}_nt"].append(seq_record_nt)

                        # Guardar la secuencia nucleotídica en un archivo temporal
                        temp_nt_fasta = f"{seq_id}_temp_nt.fasta"
                        with open(temp_nt_fasta, 'w') as temp_file:
                            SeqIO.write(seq_record_nt, temp_file, 'fasta')

                        # Definir archivo temporal para la traducción a proteínas
                        temp_aa_fasta = f"{seq_id}_temp_aa.fasta"

                        # Ejecutar transeq (EMBOSS) para traducir la secuencia nucleotídica a proteínas
                        try:
                            subprocess.run(
                                ['transeq', '-sequence', temp_nt_fasta, '-outseq', temp_aa_fasta, '-table', '1', '-frame', '1'],
                                check=True
                            )
                        except subprocess.CalledProcessError as e:
                            logging.error(f"Error al ejecutar transeq: {e}")
                            continue

                        # Leer la secuencia de proteínas traducida desde el archivo de salida de transeq
                        try:
                            protein_seq_record = next(SeqIO.parse(temp_aa_fasta, 'fasta'))
                            sequences[f"{orientation}_aa"].append(protein_seq_record)
                        except Exception as e:
                            logging.error(f"Error al leer la secuencia traducida desde {temp_aa_fasta}: {e}")
                            continue
                        finally:
                            # Eliminar archivos temporales
                            if os.path.exists(temp_nt_fasta):
                                os.remove(temp_nt_fasta)
                            if os.path.exists(temp_aa_fasta):
                                os.remove(temp_aa_fasta)

                        # Incrementar el contador según la orientación
                        if orientation == 'CD':
                            cd_count += 1
                        elif orientation == 'HO':
                            ho_count += 1

        # Generar archivos FASTA
        for seq_type, seq_list in sequences.items():
            if seq_list:
                orientation, seq_format = seq_type.split('_')
                output_file = os.path.join(output_dir, f"{group_name}_{orientation}_{seq_format}.fasta")
                os.makedirs(os.path.dirname(output_file), exist_ok=True)
                write_fasta(seq_list, output_file)
                logging.info(f"Archivo {output_file} generado con {len(seq_list)} secuencias.")
            else:
                logging.warning(f"No se encontraron secuencias para {seq_type} en el grupo {group_name}.")

        # Agregar al resumen
        gene_summary.append({
            'Gene': group_name,
            'CD': cd_count,
            'HO': ho_count
        })

    # Guardar el resumen
    save_summary_csv(gene_summary, csv_output)

def write_fasta(records, output_path):
    with open(output_path, 'w') as fasta_file:
        SeqIO.write(records, fasta_file, 'fasta')
    logging.info(f"Archivo FASTA escrito en {output_path}.")

def save_summary_csv(summary_data, output_path):
    summary_df = pd.DataFrame(summary_data)
    if os.path.exists(output_path):
        summary_df.to_csv(output_path, mode='a', header=False, index=False)
    else:
        summary_df.to_csv(output_path, index=False)
    logging.info(f"Resumen guardado en {output_path}.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Procesar clusters y generar archivos FASTA según COG y CDHO.")
    parser.add_argument("--dir", required=True, help="Directorio con archivos .gbk.")
    parser.add_argument("--cog", required=True, help="Archivo CSV con información de COG.")
    parser.add_argument("--cdho", required=True, help="Archivo .tsv con información de CDHO.")
    parser.add_argument("--output", required=True, help="Directorio de salida para los archivos FASTA.")
    parser.add_argument("--csv_output", required=True, help="Archivo de salida CSV con el resumen de conteo de genes.")

    args = parser.parse_args()
    os.makedirs(args.output, exist_ok=True)
    process_clusters(args.dir, args.cog, args.cdho, args.output, args.csv_output)
