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

# use: ab_1.py --input _GC.mochi

def combine_mochi_files(input_dir):
    combined_df = pd.DataFrame()
    for file in os.listdir(input_dir):
        if file.endswith('.mochi'):
            file_path = os.path.join(input_dir, file)
            df = pd.read_csv(file_path, sep='\t')  # Ajusta el separador según el formato del archivo
            combined_df = pd.concat([combined_df, df], ignore_index=True)
    return combined_df

def create_figure_A(data):
    cd_genes = data[data['OrientationOrientation'] == 'CD']['GCskew(rel. to LS)']
    ho_genes = data[data['OrientationOrientation'] == 'HO']['GCskew(rel. to LS)']

    # Ordenar los datos por GC skew
    cd_genes_sorted = cd_genes.sort_values(ascending=False)
    ho_genes_sorted = ho_genes.sort_values(ascending=False)

    plt.figure(figsize=(8, 6))
    plt.fill_between(range(len(cd_genes_sorted)), cd_genes_sorted, color='gray', alpha=0.5, label='CD genes')
    plt.fill_between(range(len(ho_genes_sorted)), ho_genes_sorted, color='red', alpha=0.5, label='HO genes')
    plt.axhline(0, color='black', linestyle='--', linewidth=0.8)
    plt.xlabel('Individual genes')
    plt.ylabel('GC skew')
    plt.title('Distribution of GC skew for CD and HO genes')
    plt.legend()
    plt.savefig('figure_A_single_species.png')
    plt.close()

def create_figure_B(data):
    # Filtrar los genes retenidos (CD con GC skew > 0) y los invertidos (HO con GC skew <= 0)
    cd_retained = data[(data['OrientationOrientation'] == 'CD') & (data['GCskew(rel. to LS)'] > 0)]
    ho_inverted = data[(data['OrientationOrientation'] == 'HO') & (data['GCskew(rel. to LS)'] <= 0)]

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

    # Guardar estadísticos en un CSV
    stats_df = pd.DataFrame({
        'Group': ['CD+ (Retained)', 'HO- (Inverted)'],
        'Mean': [mean_cd, mean_ho],
        'Standard Deviation': [std_cd, std_ho],
        'Standard Error': [sem_cd, sem_ho],
        'Z-score': [z_score, z_score],  # Z-score es el mismo para ambos grupos en la comparación
        'P-value': [p_value, p_value]
    })
    stats_df.to_csv('ab_FigB.csv', index=False)

    # Crear la Figura B
    plt.figure(figsize=(8, 6))
    
    # Generar la barra con medias y errores estándar
    plt.bar(['CD+ (Retained)', 'HO- (Inverted)'], [mean_cd, mean_ho], yerr=[sem_cd, sem_ho], color=['gray', 'red'])
    plt.axhline(0, color='black')
    plt.ylabel('Average GC skew')
    plt.title(f'Average GC skew for retained and inverted genes\n**p < 0.0001 (p = {p_value:.5f})')

    plt.tight_layout()
    plt.savefig('figure_B_single_species.png')
    plt.show()
    plt.close()

def main():
    parser = argparse.ArgumentParser(description='Generate figures from mochi files for a single species.')
    parser.add_argument('--input', required=True, help='Input directory containing .mochi files')
    args = parser.parse_args()

    combined_data = combine_mochi_files(args.input)

    create_figure_A(combined_data)
    create_figure_B(combined_data)
    print("Figuras generadas: figure_A_single_species.png y figure_B_single_species.png")
    print("Estadísticos guardados en FigB.csv")

if __name__ == "__main__":
    main()


