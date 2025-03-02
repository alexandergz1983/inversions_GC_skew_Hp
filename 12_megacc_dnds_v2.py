#!/usr/bin/env python
# -*- coding: utf-8 -*-
# CREATED BY: Alexander Garcia-Zea
# jagarcia@eafit.edu.co

import os
import csv
import re
import concurrent.futures
import argparse
from subprocess import check_output, CalledProcessError

# Configuración
model_ds = "modelos_mega/distance_estimation_overall_mean_Syn_syn-nonsynonymous.mao"
model_dn = "modelos_mega/distance_estimation_overall_mean_Non_syn-nonsynonymous.mao"

def run_megacc(model_file, input_fasta, output_prefix):
    megacc_path = "megacc"
    output_csv = f"{output_prefix}.csv"
    
    if os.path.exists(output_csv):
        print(f"Archivo {output_csv} ya existe. Saltando procesamiento.")
        return output_csv
    
    command = [megacc_path, "-a", model_file, "-d", input_fasta, "-o", output_csv]
    try:
        check_output(command)
        return output_csv
    except CalledProcessError as e:
        print(f"Error en MEGACC para {input_fasta}: {e}")
        return None

def process_csv(csv_file, target):
    if not csv_file or not os.path.exists(csv_file):
        return None, None
    
    try:
        with open(csv_file, 'r') as f:
            # Leer todas las líneas no vacías
            lines = [line.strip() for line in f.readlines() if line.strip()]
            
            if len(lines) < 2:
                print(f"Archivo {csv_file} incompleto. Se reprocesará.")
                os.remove(csv_file)
                return None, None
            
            # Dividir encabezado y datos usando espacios múltiples como delimitador
            header = re.split(r'\s{2,}', lines[0])
            data = re.split(r'\s{2,}', lines[1])
            
            # Validar estructura
            if len(header) < 2 or len(data) < 2:
                print(f"Formato inválido en {csv_file}. Se reprocesará.")
                os.remove(csv_file)
                return None, None
            
            # Obtener índices según el target (dS o dN)
            value_idx = 0 if target in header[0] else None
            se_idx = 1 if "S.E." in header[1] else None
            
            if value_idx is None or se_idx is None:
                print(f"Encabezado incorrecto en {csv_file}. Se reprocesará.")
                os.remove(csv_file)
                return None, None
            
            try:
                value = float(data[value_idx])
                se = float(data[se_idx])
                return (value, se)
            except (ValueError, IndexError):
                print(f"Valores inválidos en {csv_file}. Se reprocesará.")
                os.remove(csv_file)
                return None, None
            
    except Exception as e:
        print(f"Error crítico leyendo {csv_file}: {str(e)}")
        return None, None

def create_summary_table(input_dir, model_ds, model_dn, output_table, max_workers=8):
    results = {}
    
    def process_file(file):
        if file.endswith(".fasta"):
            match = re.match(r"(.*?)_(CD|HO)_nt\.fasta", file)
            if not match:
                return
            cluster_gene = match.group(1)
            condition = match.group(2)
            input_file = os.path.join(input_dir, file)
            
            # Generar nombres de archivos CSV
            base_name = f"{cluster_gene}_{condition}"
            ds_csv = f"{base_name}_ds.csv"
            dn_csv = f"{base_name}_dn.csv"
            
            # Procesar dS
            if not os.path.exists(ds_csv):
                ds_csv = run_megacc(model_ds, input_file, base_name + "_ds")
            ds, ds_se = process_csv(ds_csv, "dS") if ds_csv else (None, None)
            
            # Procesar dN
            if not os.path.exists(dn_csv):
                dn_csv = run_megacc(model_dn, input_file, base_name + "_dn")
            dn, dn_se = process_csv(dn_csv, "dN") if dn_csv else (None, None)
            
            return (cluster_gene, condition, ds, ds_se, dn, dn_se)
    
    with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = [executor.submit(process_file, file) for file in os.listdir(input_dir) if file.endswith(".fasta")]
        for future in concurrent.futures.as_completed(futures):
            result = future.result()
            if result:
                cluster_gene, condition, ds, ds_se, dn, dn_se = result
                if cluster_gene not in results:
                    results[cluster_gene] = {}
                results[cluster_gene][condition] = (ds, ds_se, dn, dn_se)
    
    # Escribir resultados
    with open(output_table, 'w', newline='') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(["Cluster", "CD_dS", "CD_SE", "CD_dN", "CD_SE", "HO_dS", "HO_SE", "HO_dN", "HO_SE"])
        for cluster_gene in sorted(results.keys()):
            conditions = results[cluster_gene]
            cd_data = conditions.get("CD", (None, None, None, None))
            ho_data = conditions.get("HO", (None, None, None, None))
            
            writer.writerow([
                cluster_gene,
                cd_data[0] or "",
                cd_data[1] or "",
                cd_data[2] or "",
                cd_data[3] or "",
                ho_data[0] or "",
                ho_data[1] or "",
                ho_data[2] or "",
                ho_data[3] or ""
            ])

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Automatizar ejecución de MEGACC para dS y dN en alineamientos FASTA.")
    parser.add_argument("--input", required=True, help="Directorio con alineamientos FASTA.")
    parser.add_argument("--output", default="resultados_dn_ds.csv", help="Archivo CSV de salida.")
    parser.add_argument("--force", action="store_true", help="Forzar reprocesamiento de todos los archivos")
    parser.add_argument("--max-workers", type=int, default=8, help="Número de núcleos a usar.")
    args = parser.parse_args()
    
    if args.force:
        print("Modo forzado: Se reprocesarán todos los archivos")
        for f in os.listdir(args.input):
            if f.endswith(("_ds.csv", "_dn.csv")):
                os.remove(os.path.join(args.input, f))
    
    create_summary_table(args.input, model_ds, model_dn, args.output, max_workers=args.max_workers)
    print(f"Tabla de resultados generada: {args.output}")
