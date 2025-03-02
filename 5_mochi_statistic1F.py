#!/usr/bin/env python
# -*- coding: utf-8 -*-
# CREATED BY: Jerson Alexander Garcia-Zea
# jagarcia@eafit.edu.co

# use: python mochi_statisticsV1F.py --input gc_mochi/

import os
import pandas as pd
import numpy as np
import argparse
import matplotlib.pyplot as plt

def calculate_statistics(file_path):
    df = pd.read_csv(file_path, sep='\t')

    # Calculate CD and HO frequency
    cd_frequencies = (df['OrientationOrientation'] == 'CD').sum()
    ho_frequencies = (df['OrientationOrientation'] == 'HO').sum()
    total_genes = cd_frequencies + ho_frequencies

    # Calculate A, T, C, G frequency
    A_frequencies = df['#A'].sum()
    T_frequencies = df['#T'].sum()
    C_frequencies = df['#C'].sum()
    G_frequencies = df['#G'].sum()
    total_bases = A_frequencies + T_frequencies + C_frequencies + G_frequencies

    # Calculate AT and GC frequencies
    AT_frequencies = A_frequencies + T_frequencies
    GC_frequencies = G_frequencies + C_frequencies

    # Calculate mean and standard deviation values
    mean_A = df['#A'].mean()
    std_A = df['#A'].std()
    mean_T = df['#T'].mean()
    std_T = df['#T'].std()
    mean_C = df['#C'].mean()
    std_C = df['#C'].std()
    mean_G = df['#G'].mean()
    std_G = df['#G'].std()
    mean_AT = (df['#A'] + df['#T']).mean()
    std_AT = (df['#A'] + df['#T']).std()
    mean_GC = (df['#G'] + df['#C']).mean()
    std_GC = (df['#G'] + df['#C']).std()

    mean_gccontent = df['GCcontent'].mean()
    mean_gcskew = df['GCskew(rel. to LS)'].mean()

    # Calculate proportions
    cd_proportion = cd_frequencies / total_genes if total_genes > 0 else 0
    ho_proportion = ho_frequencies / total_genes if total_genes > 0 else 0

    A_proportion = A_frequencies / total_bases if total_bases > 0 else 0
    T_proportion = T_frequencies / total_bases if total_bases > 0 else 0
    C_proportion = C_frequencies / total_bases if total_bases > 0 else 0
    G_proportion = G_frequencies / total_bases if total_bases > 0 else 0

    AT_proportion = AT_frequencies / total_bases if total_bases > 0 else 0
    GC_proportion = GC_frequencies / total_bases if total_bases > 0 else 0

    return (cd_frequencies, ho_frequencies, mean_gccontent, mean_gcskew, cd_proportion, ho_proportion,
            A_frequencies, T_frequencies, C_frequencies, G_frequencies, 
            AT_frequencies, GC_frequencies, 
            mean_A, std_A, mean_T, std_T, mean_C, std_C, mean_G, std_G, mean_AT, std_AT, mean_GC, std_GC,
            A_proportion, T_proportion, C_proportion, G_proportion, AT_proportion, GC_proportion)


def process_directory(input_dir):
    results = []

    for filename in os.listdir(input_dir):
        if filename.endswith('.mochi'):
            file_path = os.path.join(input_dir, filename)
            (cd_frequencies, ho_frequencies, mean_gccontent, mean_gcskew, cd_proportion, ho_proportion,
             A_frequencies, T_frequencies, C_frequencies, G_frequencies,
             AT_frequencies, GC_frequencies, 
             mean_A, std_A, mean_T, std_T, mean_C, std_C, mean_G, std_G, mean_AT, std_AT, mean_GC, std_GC,
             A_proportion, T_proportion, C_proportion, G_proportion, AT_proportion, GC_proportion) = calculate_statistics(file_path)

            results.append({
                'Genoma': filename,
                'CD frequencies': cd_frequencies,
                'HO frequencies': ho_frequencies,
                'Mean GC Content': mean_gccontent,
                'Mean GC Skew': mean_gcskew,
                'CD Proportion': cd_proportion,
                'HO Proportion': ho_proportion,
                'A Frequencies': A_frequencies,
                'T Frequencies': T_frequencies,
                'C Frequencies': C_frequencies,
                'G Frequencies': G_frequencies,
                'AT Frequencies': AT_frequencies,
                'GC Frequencies': GC_frequencies,
                'Mean A ± STD': mean_A,
                'Std A': std_A,
                'Mean T ± STD': mean_T,
                'Std T': std_T,
                'Mean C ± STD': mean_C,
                'Std C': std_C,
                'Mean G ± STD': mean_G,
                'Std G': std_G,
                'Mean AT ± STD': mean_AT,
                'Std AT': std_AT,
                'Mean GC ± STD': mean_GC,
                'Std GC': std_GC,
                'A Proportion': A_proportion,
                'T Proportion': T_proportion,
                'C Proportion': C_proportion,
                'G Proportion': G_proportion,
                'AT Proportion': AT_proportion,
                'GC Proportion': GC_proportion
            })

    results_df = pd.DataFrame(results)

    # Cálculo de medias y desviaciones estándar totales
    summary_row = {
        'Genoma': 'TOTAL',
        'CD frequencies': results_df['CD frequencies'].sum(),
        'HO frequencies': results_df['HO frequencies'].sum(),
        'Mean GC Content': results_df['Mean GC Content'].mean(),
        'Mean GC Skew': results_df['Mean GC Skew'].mean(),
        'CD Proportion': results_df['CD frequencies'].sum() / (results_df['CD frequencies'].sum() + results_df['HO frequencies'].sum()),
        'HO Proportion': results_df['HO frequencies'].sum() / (results_df['CD frequencies'].sum() + results_df['HO frequencies'].sum()),
        'A Frequencies': results_df['A Frequencies'].sum(),
        'T Frequencies': results_df['T Frequencies'].sum(),
        'C Frequencies': results_df['C Frequencies'].sum(),
        'G Frequencies': results_df['G Frequencies'].sum(),
        'AT Frequencies': results_df['AT Frequencies'].sum(),
        'GC Frequencies': results_df['GC Frequencies'].sum(),
        'Mean A ± STD': f"{results_df['Mean A ± STD'].mean():.2f} ± {results_df['Std A'].std():.2f}",
        'Mean T ± STD': f"{results_df['Mean T ± STD'].mean():.2f} ± {results_df['Std T'].std():.2f}",
        'Mean C ± STD': f"{results_df['Mean C ± STD'].mean():.2f} ± {results_df['Std C'].std():.2f}",
        'Mean G ± STD': f"{results_df['Mean G ± STD'].mean():.2f} ± {results_df['Std G'].std():.2f}",
        'Mean AT ± STD': f"{results_df['Mean AT ± STD'].mean():.2f} ± {results_df['Std AT'].std():.2f}",
        'Mean GC ± STD': f"{results_df['Mean GC ± STD'].mean():.2f} ± {results_df['Std GC'].std():.2f}",
        'A Proportion': results_df['A Frequencies'].sum() / (results_df['A Frequencies'].sum() + results_df['T Frequencies'].sum() + results_df['C Frequencies'].sum() + results_df['G Frequencies'].sum()),
        'T Proportion': results_df['T Frequencies'].sum() / (results_df['A Frequencies'].sum() + results_df['T Frequencies'].sum() + results_df['C Frequencies'].sum() + results_df['G Frequencies'].sum()),
        'C Proportion': results_df['C Frequencies'].sum() / (results_df['A Frequencies'].sum() + results_df['T Frequencies'].sum() + results_df['C Frequencies'].sum() + results_df['G Frequencies'].sum()),
        'G Proportion': results_df['G Frequencies'].sum() / (results_df['A Frequencies'].sum() + results_df['T Frequencies'].sum() + results_df['C Frequencies'].sum() + results_df['G Frequencies'].sum()),
        'AT Proportion': results_df['AT Frequencies'].sum() / (results_df['A Frequencies'].sum() + results_df['T Frequencies'].sum() + results_df['C Frequencies'].sum() + results_df['G Frequencies'].sum()),
        'GC Proportion': results_df['GC Frequencies'].sum() / (results_df['A Frequencies'].sum() + results_df['T Frequencies'].sum() + results_df['C Frequencies'].sum() + results_df['G Frequencies'].sum())
    }

    results_df = pd.concat([results_df, pd.DataFrame([summary_row])], ignore_index=True)
    return results_df

def create_pie_charts(results_df):
    total_result = results_df[results_df['Genoma'] == 'TOTAL'].iloc[0]

    if not total_result.empty:

        # Gráfico de Proporción Total entre Genes CD y HO
        fig1, ax1 = plt.subplots()
        ax1.pie([total_result['CD frequencies'], total_result['HO frequencies']], 
                labels=['CD', 'HO'], autopct='%1.1f%%', colors=['red', 'gray'])
        ax1.axis('equal')
        plt.title('Proporción Total entre Genes CD y HO')
        plt.savefig('CD_HO_pie_chart.png')

        # Gráfico de Proporción Total entre Bases
        fig2, ax2 = plt.subplots()
        ax2.pie([total_result['A Frequencies'], total_result['T Frequencies'], 
                 total_result['C Frequencies'], total_result['G Frequencies']],
                labels=['A', 'T', 'C', 'G'], autopct='%1.1f%%', 
                colors=['yellowgreen', 'gold', 'lightskyblue', 'lightcoral'])
        ax2.axis('equal')
        plt.title('Proporción Total entre Bases A, T, C, G')
        plt.savefig('Bases_pie_chart.png')

        # Gráfico de Proporción Total de AT vs GC
        fig3, ax3 = plt.subplots()
        ax3.pie([total_result['AT Frequencies'], total_result['GC Frequencies']], 
                labels=['AT', 'GC'], autopct='%1.1f%%', colors=['plum', 'lightblue'])
        ax3.axis('equal')
        plt.title('Proporción Total de AT vs GC')
        plt.savefig('AT_GC_pie_chart.png')

def main():
    parser = argparse.ArgumentParser(description='Process .mochi files for CD and HO frequency analysis')
    parser.add_argument('--input', required=True, help='Input directory containing .mochi files')
    args = parser.parse_args()

    results_df = process_directory(args.input)

    # Save results to an output file
    results_df.to_csv('Summary_Results_MochiStatisticV1F.csv', sep='\t', index=False)
    print("Analysis completed. Results saved to gc_proportions_mochi.txt.")

    # Create pie charts
    create_pie_charts(results_df)

if __name__ == '__main__':
    main()

