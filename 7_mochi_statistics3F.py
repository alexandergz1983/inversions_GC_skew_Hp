#!/usr/bin/env python
# -*- coding: utf-8 -*-
# CREATED BY: Alexander Garcia-Zea
# jagarcia@eafit.edu.co

import os
import pandas as pd
import matplotlib.pyplot as plt
import argparse

def parse_arguments():
    parser = argparse.ArgumentParser(description='Analizar GC skew de archivos tabulares en un directorio.')
    parser.add_argument('--input', required=True, help='Directorio que contiene los archivos tabulares.')
    return parser.parse_args()

def load_data(directory):
    # Combinar todos los archivos en un solo DataFrame
    data_frames = []
    for file in os.listdir(directory):
        if file.endswith('.mochi'):  # Asegúrate de que los archivos tengan la extensión correcta
            file_path = os.path.join(directory, file)
            df = pd.read_csv(file_path, sep='\t')  # Asume que los archivos están delimitados por tabulaciones
            data_frames.append(df)
    combined_data = pd.concat(data_frames, ignore_index=True)
    return combined_data

def plot_gc_skew(data):
    # Filtrar los genes según el GC skew positivo (>0) y negativo (<=0)
    genes_positive_gc_skew = data[data['GCskew(rel. to LS)'] > 0]
    genes_negative_gc_skew = data[data['GCskew(rel. to LS)'] <= 0]

    # Separar los genes por orientación
    cd_genes_positive = genes_positive_gc_skew[genes_positive_gc_skew['OrientationOrientation'] == 'CD']
    ho_genes_positive = genes_positive_gc_skew[genes_positive_gc_skew['OrientationOrientation'] == 'HO']

    cd_genes_negative = genes_negative_gc_skew[genes_negative_gc_skew['OrientationOrientation'] == 'CD']
    ho_genes_negative = genes_negative_gc_skew[genes_negative_gc_skew['OrientationOrientation'] == 'HO']

    # Crear la figura A
    plt.figure(figsize=(10, 5))
    plt.subplot(1, 2, 1)
    plt.bar(range(len(cd_genes_positive)), cd_genes_positive['GCskew(rel. to LS)'], color='gray', label='CD genes (retained)')
    plt.bar(range(len(ho_genes_positive)), ho_genes_positive['GCskew(rel. to LS)'], color='black', label='HO genes (retained)')
    plt.axhline(0, color='black', linestyle='--')
    plt.xlabel('Genes individuales')
    plt.ylabel('GC skew')
    plt.title('Genes Retenidos')
    plt.legend()

    # Crear la figura para los genes afectados por inversión
    plt.subplot(1, 2, 2)
    plt.bar(range(len(cd_genes_negative)), cd_genes_negative['GCskew(rel. to LS)'], color='lightcoral', label='CD genes (inverted)')
    plt.bar(range(len(ho_genes_negative)), ho_genes_negative['GCskew(rel. to LS)'], color='darkred', label='HO genes (inverted)')
    plt.axhline(0, color='black', linestyle='--')
    plt.xlabel('Genes individuales')
    plt.ylabel('GC skew')
    plt.title('Genes Afectados por Inversión')
    plt.legend()

    plt.tight_layout()
    plt.show()

    # Calcular las medias y desviaciones estándar
    mean_retained_cd = cd_genes_positive['GCskew(rel. to LS)'].mean()
    mean_retained_ho = ho_genes_positive['GCskew(rel. to LS)'].mean()
    mean_inverted_cd = cd_genes_negative['GCskew(rel. to LS)'].mean()
    mean_inverted_ho = ho_genes_negative['GCskew(rel. to LS)'].mean()

    std_retained_cd = cd_genes_positive['GCskew(rel. to LS)'].std()
    std_retained_ho = ho_genes_positive['GCskew(rel. to LS)'].std()
    std_inverted_cd = cd_genes_negative['GCskew(rel. to LS)'].std()
    std_inverted_ho = ho_genes_negative['GCskew(rel. to LS)'].std()

    # Crear la figura B
    plt.figure(figsize=(8, 4))
    plt.bar(['CD+ Retained', 'HO+ Retained', 'CD- Inverted', 'HO- Inverted'],
            [mean_retained_cd, mean_retained_ho, mean_inverted_cd, mean_inverted_ho],
            yerr=[std_retained_cd, std_retained_ho, std_inverted_cd, std_inverted_ho],
            color=['gray', 'black', 'lightcoral', 'darkred'])
    plt.axhline(0, color='black')
    plt.ylabel('Average GC skew')
    plt.title('Media de GC skew por Grupo de Genes')
    
    plt.tight_layout()
    plt.show()

def main():
    args = parse_arguments()
    input_directory = args.input

    # Cargar y combinar los datos
    data = load_data(input_directory)

    # Generar los gráficos
    plot_gc_skew(data)

if __name__ == "__main__":
    main()

