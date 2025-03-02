#!/usr/bin/env python
# -*- coding: utf-8 -*-
# CREATED BY: Alexander Garcia-Zea
# jagarcia@eafit.edu.co

import os
import subprocess
import argparse

def ejecutar_prokka(input_dir, output_dir):
    # Verificar que el directorio de entrada existe
    if not os.path.isdir(input_dir):
        print(f"El directorio {input_dir} no existe.")
        return

    # Verificar que el directorio de salida existe o crearlo
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Listar archivos en el directorio que terminen en .fa o .fasta
    for filename in os.listdir(input_dir):
        if filename.endswith(".fa") or filename.endswith(".fasta"):
            # Crear el prefijo basado en el nombre original del archivo
            base_name = os.path.splitext(filename)[0]  # Nombre sin extensión
            prefix = base_name  # Mantener el nombre original como prefijo

            # Construir el comando para ejecutar Prokka con las opciones necesarias
            input_file = os.path.join(input_dir, filename)
            cmd = [
                "prokka", input_file, 
                "--outdir", output_dir,  # Un solo directorio de salida
                "--prefix", prefix,      # Mantener el nombre original como prefijo
                "--force",                # Sobrescribe archivos si ya existen
                "--kingdom", "Bacteria",   # Especificar el reino
                "--genus", "Helicobacter",
                "--compliant",            # Cumplir con las reglas de anotación
            ]

            try:
                # Ejecutar el comando Prokka
                print(f"Ejecutando Prokka en {filename}...")
                subprocess.run(cmd, check=True)
                
                # Verificar si los archivos .gff y .faa fueron creados correctamente
                gff_file = os.path.join(output_dir, f"{prefix}.gff")
                faa_file = os.path.join(output_dir, f"{prefix}.faa")
                
                if os.path.exists(gff_file) and os.path.exists(faa_file):
                    print(f"Anotación completada para {filename}. Resultados en {output_dir}/")
                else:
                    print(f"Error: Los archivos {gff_file} o {faa_file} no fueron generados correctamente.")
                    
            except subprocess.CalledProcessError as e:
                print(f"Error al ejecutar Prokka en {filename}: {e}")
            except Exception as e:
                print(f"Error inesperado: {e}")

if __name__ == "__main__":
    # Configuración de argumentos para la línea de comandos
    parser = argparse.ArgumentParser(description="Ejecutar Prokka en múltiples archivos .fasta en un directorio.")
    parser.add_argument("--input_dir", required=True, help="Directorio que contiene archivos .fa o .fasta")
    parser.add_argument("--output_dir", required=True, help="Directorio donde se guardarán los archivos anotados")
    
    args = parser.parse_args()
    
    # Llamar a la función para ejecutar Prokka en los archivos del directorio
    ejecutar_prokka(args.input_dir, args.output_dir)
