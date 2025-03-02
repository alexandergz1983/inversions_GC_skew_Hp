#!/usr/bin/env python
# -*- coding: utf-8 -*-
# CREATED BY: Jerson Alexander Garcia-Zea
# jagarcia@eafit.edu.co

# use: python mochi_statisticsV2F.py --input gc_mochi/ --output summary_mochi_statisticsV2.csv

import os
import pandas as pd
import matplotlib.pyplot as plt
import argparse

def parse_arguments():
    parser = argparse.ArgumentParser(description='Analizar GC skew de archivos tabulares en un directorio.')
    parser.add_argument('--input', required=True, help='Directorio que contiene los archivos tabulares.')
    parser.add_argument('--output', required=True, help='Archivo de salida CSV para estadísticas.')
    return parser.parse_args()

def load_data(directory):
    # Combinar todos los archivos en un solo DataFrame
    data_frames = []
    for file in os.listdir(directory):
        if file.endswith('.mochi'):
            file_path = os.path.join(directory, file)
            df = pd.read_csv(file_path, sep='\t')
            df['Genoma'] = file.replace('.mochi', '')  # Agregar una columna con el nombre del genoma
            data_frames.append(df)
    combined_data = pd.concat(data_frames, ignore_index=True)
    return combined_data

def calculate_frequencies_and_proportions(data):
    # Calcula frecuencias y proporciones por cada condición
    conditions = {
        'CD <= 0': (data['OrientationOrientation'] == 'CD') & (data['GCskew(rel. to LS)'] <= 0),
        'CD > 0': (data['OrientationOrientation'] == 'CD') & (data['GCskew(rel. to LS)'] > 0),
        'HO <= 0': (data['OrientationOrientation'] == 'HO') & (data['GCskew(rel. to LS)'] <= 0),
        'HO > 0': (data['OrientationOrientation'] == 'HO') & (data['GCskew(rel. to LS)'] > 0)
    }

    frequencies = {label: data[cond].shape[0] for label, cond in conditions.items()}
    total_genes = sum(frequencies.values())
    proportions = {label: freq / total_genes for label, freq in frequencies.items()}
    
    return frequencies, proportions

def calculate_mean_std(data):
    # Calcula la media y desviación estándar por condición
    conditions = {
        'CD <= 0': (data['OrientationOrientation'] == 'CD') & (data['GCskew(rel. to LS)'] <= 0),
        'CD > 0': (data['OrientationOrientation'] == 'CD') & (data['GCskew(rel. to LS)'] > 0),
        'HO <= 0': (data['OrientationOrientation'] == 'HO') & (data['GCskew(rel. to LS)'] <= 0),
        'HO > 0': (data['OrientationOrientation'] == 'HO') & (data['GCskew(rel. to LS)'] > 0)
    }

    stats = {}
    for label, cond in conditions.items():
        subset = data[cond]
        mean = subset['GCskew(rel. to LS)'].mean()
        std = subset['GCskew(rel. to LS)'].std()
        stats[label] = (mean, std)
    
    return stats

def save_statistics_to_csv(data, output_file):
    # Generar estadísticas para cada genoma y una fila TOTAL
    summary = []
    
    total_data = pd.DataFrame()  # Para acumular datos de todos los genomas
    
    for genoma, group in data.groupby('Genoma'):
        frequencies, proportions = calculate_frequencies_and_proportions(group)
        stats = calculate_mean_std(group)
        
        summary_row = {
            'Genoma': genoma,
            'CD frequencies <= 0': frequencies['CD <= 0'],
            'CD frequencies > 0': frequencies['CD > 0'],
            'HO frequencies <= 0': frequencies['HO <= 0'],
            'HO frequencies > 0': frequencies['HO > 0'],
            'CD Proportion <= 0': proportions['CD <= 0'],
            'CD Proportion > 0': proportions['CD > 0'],
            'HO Proportion <= 0': proportions['HO <= 0'],
            'HO Proportion > 0': proportions['HO > 0'],
            'Mean CD <= 0': stats['CD <= 0'][0],
            'STD CD <= 0': stats['CD <= 0'][1],
            'Mean CD > 0': stats['CD > 0'][0],
            'STD CD > 0': stats['CD > 0'][1],
            'Mean HO <= 0': stats['HO <= 0'][0],
            'STD HO <= 0': stats['HO <= 0'][1],
            'Mean HO > 0': stats['HO > 0'][0],
            'STD HO > 0': stats['HO > 0'][1]
        }
        
        summary.append(summary_row)
        total_data = pd.concat([total_data, group], ignore_index=True)  # Acumular todos los datos
    
    # Calcular estadísticas totales
    frequencies_total, proportions_total = calculate_frequencies_and_proportions(total_data)
    stats_total = calculate_mean_std(total_data)
    
    # Crear fila TOTAL
    total_row = {
        'Genoma': 'TOTAL',
        'CD frequencies <= 0': frequencies_total['CD <= 0'],
        'CD frequencies > 0': frequencies_total['CD > 0'],
        'HO frequencies <= 0': frequencies_total['HO <= 0'],
        'HO frequencies > 0': frequencies_total['HO > 0'],
        'CD Proportion <= 0': proportions_total['CD <= 0'],
        'CD Proportion > 0': proportions_total['CD > 0'],
        'HO Proportion <= 0': proportions_total['HO <= 0'],
        'HO Proportion > 0': proportions_total['HO > 0'],
        'Mean CD <= 0': stats_total['CD <= 0'][0],
        'STD CD <= 0': stats_total['CD <= 0'][1],
        'Mean CD > 0': stats_total['CD > 0'][0],
        'STD CD > 0': stats_total['CD > 0'][1],
        'Mean HO <= 0': stats_total['HO <= 0'][0],
        'STD HO <= 0': stats_total['HO <= 0'][1],
        'Mean HO > 0': stats_total['HO > 0'][0],
        'STD HO > 0': stats_total['HO > 0'][1]
    }
    
    # Añadir fila TOTAL a las estadísticas
    summary.append(total_row)
    
    # Convertir a DataFrame y guardar en CSV
    summary_df = pd.DataFrame(summary)
    summary_df.to_csv(output_file, sep=';', index=False)
    return summary_df

def plot_pie_chart(data, title):
    labels = data.keys()
    sizes = data.values()
    plt.figure(figsize=(8, 8))
    plt.pie(sizes, labels=labels, autopct='%1.1f%%', startangle=140)
    plt.title(title)
    plt.axis('equal')
    plt.show()

import matplotlib.pyplot as plt

# Modificación de la función plot_pie_chart para aceptar colores
def plot_pie_chart(data, title, colors=None):
    labels = data.keys()
    sizes = data.values()
    
    # Si no se proporcionan colores, se usan los predeterminados
    if colors is None:
        colors = plt.get_cmap('tab10').colors[:len(data)]  # Usa una paleta predeterminada
    
    plt.figure(figsize=(6, 6))
    plt.pie(sizes, labels=labels, autopct='%1.1f%%', colors=colors)
    plt.title(title)
    plt.show()

# Función plot_comparisons modificada para incluir colores
def plot_comparisons(frequencies, proportions, colors=None):
    # Graficar las tortas solicitadas
    plot_pie_chart(frequencies, 'Frecuencia total de genes')
    plot_pie_chart(proportions, 'Proporción total de genes')

    # Comparación entre genes CD y HO
    plot_pie_chart(
        {'CD <= 0': proportions['CD <= 0'], 'HO <= 0': proportions['HO <= 0']},
        'Proporción de Genes CD vs HO con GC skew <= 0',
        colors=colors
    )
    plot_pie_chart(
        {'CD > 0': proportions['CD > 0'], 'HO > 0': proportions['HO > 0']},
        'Proporción de Genes CD vs HO con GC skew > 0',
        colors=colors
    )
    
    # Comparación dentro de los mismos grupos
    plot_pie_chart(
        {'CD <= 0': proportions['CD <= 0'], 'CD > 0': proportions['CD > 0']},
        'Proportion of Genes CD with GC skew <= 0 vs > 0',
        colors=['red', 'gray']  # Rojo para <= 0 y gris para > 0
    )

    plot_pie_chart(
        {'HO <= 0': proportions['HO <= 0'], 'HO > 0': proportions['HO > 0']},
        'Proportion of HO Genes with GC skew <= 0 vs > 0',
        colors=['red', 'gray']  # Rojo para <= 0 y gris para > 0
    )

def main():
    args = parse_arguments()
    data = load_data(args.input)
    frequencies, proportions = calculate_frequencies_and_proportions(data)
    plot_comparisons(frequencies, proportions)
    save_statistics_to_csv(data, args.output)

if __name__ == '__main__':
    main()
