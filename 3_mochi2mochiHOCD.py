#!/usr/bin/env python
# -*- coding: utf-8 -*-
# CREATED BY: Alexander Garcia-Zea
# jagarcia@eafit.edu.co

import pandas as pd
import argparse
import os
import glob

# Argumentos de entrada
parser = argparse.ArgumentParser(description="Procesa todos los archivos .mochi en un directorio")
parser.add_argument('--indir', required=True, help="Directorio de entrada con archivos .mochi")
args = parser.parse_args()

# Directorio de entrada y salida
input_dir = args.indir
output_dir = os.path.join(input_dir, "mochi_f")

# Crear el directorio de salida si no existe
os.makedirs(output_dir, exist_ok=True)

# Procesar todos los archivos .mochi en el directorio de entrada
for mochi_file in glob.glob(os.path.join(input_dir, "*.mochi")):
    try:
        # Leer el archivo .mochi
        df = pd.read_csv(mochi_file, sep="\t")
        
        # Crear el DataFrame sin la columna 'AAseq'
        df_final = df.drop(['AAseq'], axis=1)

        # Determinar la orientación basada en la hebra
        for i in range(len(df_final)):
            if df_final.loc[i, "STRAND"] == '-':
                df_final.loc[i, "Orientation"] = "HO"
            else:
                df_final.loc[i, "Orientation"] = "CD"
        
        # Generar el nombre del archivo de salida con el sufijo "_HOCD"
        base_name = os.path.basename(mochi_file).replace(".mochi", "")
        output_file = os.path.join(output_dir, f"{base_name}_HOCD.mochi")
        
        # Guardar el archivo procesado en el directorio de salida
        df_final.to_csv(output_file, sep="\t", index=False)
        print(f"Archivo procesado: {output_file}")

    except Exception as e:
        print(f"Error al procesar {mochi_file}: {e}")


