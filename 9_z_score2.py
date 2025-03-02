#!/usr/bin/env python
# -*- coding: utf-8 -*-
# CREATED BY: Alexander Garcia-Zea
# jagarcia@eafit.edu.co

import os
import pandas as pd
import matplotlib.pyplot as plt
import argparse
from scipy.stats import norm
import numpy as np

# use: z.py --input _GC.mochi

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
    # Separar los genes CD y HO
    cd_genes = data[data['OrientationOrientation'] == 'CD']
    ho_genes = data[data['OrientationOrientation'] == 'HO']

    # Ordenar los datos por GC skew
    cd_genes_sorted = cd_genes.sort_values(by='GCskew(rel. to LS)', ascending=False)
    ho_genes_sorted = ho_genes.sort_values(by='GCskew(rel. to LS)', ascending=False)

    # Create figure A
    plt.figure(figsize=(8, 4))
    plt.subplot(1, 2, 1)
    plt.bar(range(len(cd_genes_sorted)), cd_genes_sorted['GCskew(rel. to LS)'], color='gray', label='CD genes')
    plt.bar(range(len(ho_genes_sorted)), ho_genes_sorted['GCskew(rel. to LS)'], color='red', label='HO genes')
    plt.axhline(0, color='black', linestyle='--')
    plt.xlabel('Individual genes')
    plt.ylabel('GC skew')
    plt.legend()
    
    # Figure B
    # Filtrar genes CD con GC skew > 0 (retenidos) y HO con GC skew <= 0 (invertidos) para la Figura B
    cd_retained = cd_genes[cd_genes['GCskew(rel. to LS)'] > 0]
    ho_inverted = ho_genes[ho_genes['GCskew(rel. to LS)'] <= 0]

    # Calcular medias y errores estándar
    mean_cd = cd_retained['GCskew(rel. to LS)'].mean()
    mean_ho = ho_inverted['GCskew(rel. to LS)'].mean()
    std_cd = cd_retained['GCskew(rel. to LS)'].std()
    std_ho = ho_inverted['GCskew(rel. to LS)'].std()
    sem_cd = std_cd / np.sqrt(len(cd_retained))
    sem_ho = std_ho / np.sqrt(len(ho_inverted))

    # Prueba Z para comparar los dos grupos
    z_score = (mean_cd - mean_ho) / np.sqrt(sem_cd**2 + sem_ho**2)
    p_value = 2 * (1 - norm.cdf(abs(z_score)))

    # Guardar los estadísticos en un archivo CSV
    stats_df = pd.DataFrame({
        'Category': ['CD+ (Retained)', 'HO- (Inverted)'],
        'Mean GC Skew': [mean_cd, mean_ho],
        'Standard Error': [sem_cd, sem_ho],
        'Z-score': [z_score, np.nan],  # Solo hay un Z-score para la comparación
        'P-value': [p_value, np.nan]   # Solo hay un p-valor para la comparación
    })
    stats_df.to_csv('z_FigB.csv', index=False)

    # Crear la figura B: Medias y errores estándar
    plt.subplot(1, 2, 2)
    plt.bar(['CD+ (Retained)', 'HO- (Inverted)'], [mean_cd, mean_ho], yerr=[sem_cd, sem_ho], color=['gray', 'red'])
    plt.axhline(0, color='black')
    plt.ylabel('Average GC skew')
    plt.title(f'**p < 0.0001 (p = {p_value:.5f})')
    
    plt.tight_layout()
    plt.show()

def main():
    args = parse_arguments()
    input_directory = args.input

    # Cargar y combinar los datos
    data = load_data(input_directory)

    # Generar los gráficos y guardar estadísticas
    plot_gc_skew(data)

if __name__ == "__main__":
    main()


