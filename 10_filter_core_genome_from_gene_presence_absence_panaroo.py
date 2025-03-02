#!/usr/bin/env python
# -*- coding: utf-8 -*-
# CREATED BY: Alexander Garcia-Zea
# jagarcia@eafit.edu.co

import pandas as pd
import argparse

def filter_core_genome(input_file, list_file, output_file):
    # Cargar el archivo gene_presence_absence.csv
    df = pd.read_csv(input_file)
    
    # Cargar la lista de genes centrales
    with open(list_file, 'r') as f:
        core_genes = {line.strip() for line in f}
    
    # Verificar que la columna 'Gene' exista
    if 'Gene' not in df.columns:
        raise ValueError("El archivo de entrada no contiene la columna 'Gene'.")
    
    # Filtrar filas donde el 'Gene' esté en la lista de genes centrales
    filtered_df = df[df['Gene'].isin(core_genes)]
    
    # Guardar el resultado en un archivo CSV
    filtered_df.to_csv(output_file, index=False)
    print(f"Archivo de genoma central guardado en {output_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Filtrar el genoma central basado en una lista externa y un archivo gene_presence_absence.csv generado por Panaroo.")
    parser.add_argument("--input", required=True, help="Ruta al archivo gene_presence_absence.csv de entrada.")
    parser.add_argument("--list", required=True, help="Ruta al archivo de texto que contiene la lista del genoma central (un gen por línea).")
    parser.add_argument("--output", required=True, help="Ruta al archivo CSV de salida que contendrá solo los genes del genoma central.")
    
    args = parser.parse_args()
    filter_core_genome(args.input, args.list, args.output)

