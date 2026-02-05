#!/usr/bin/env python3
"""
Análisis de Codones ATG en E. coli K-12 MG1655
==============================================
Este script identifica y cuantifica todos los codones ATG presentes
en el genoma completo, calculando su densidad y comparándolos con
los genes oficialmente anotados.

Autor: Proyecto Bioinformática
Fecha: 2024
"""

import os
import json
from Bio import SeqIO
from collections import Counter
import pandas as pd

# ============================================================================
# CONFIGURACIÓN
# ============================================================================

INPUT_FILE = os.path.expanduser('~/data/raw/ecoli_k12_mg1655.gbk')
OUTPUT_DIR = os.path.expanduser('~/projects/bioinfo/results/tables')
os.makedirs(OUTPUT_DIR, exist_ok=True)

# ============================================================================
# FUNCIONES
# ============================================================================

def contar_codon(secuencia, codon):
    """
    Cuenta todas las apariciones de un codón en una secuencia.
    
    Args:
        secuencia (str): Secuencia de ADN
        codon (str): Codón a buscar (ej: 'ATG')
    
    Returns:
        list: Lista de posiciones donde aparece el codón (0-indexed)
    """
    secuencia = str(secuencia).upper()
    codon = codon.upper()
    posiciones = []
    
    # Buscar en las 3 fases de lectura
    for i in range(len(secuencia) - len(codon) + 1):
        if secuencia[i:i+len(codon)] == codon:
            posiciones.append(i)
    
    return posiciones

def extraer_genes_anotados(record):
    """
    Extrae información de genes anotados del archivo GenBank.
    
    Args:
        record: Objeto SeqRecord de BioPython
    
    Returns:
        dict: Diccionario con estadísticas de genes
    """
    genes = []
    cds_features = []
    
    for feature in record.features:
        if feature.type == 'gene':
            genes.append({
                'tipo': 'gene',
                'inicio': int(feature.location.start),
                'fin': int(feature.location.end),
                'strand': feature.location.strand,
                'locus_tag': feature.qualifiers.get('locus_tag', ['NA'])[0]
            })
        elif feature.type == 'CDS':
            cds_features.append({
                'tipo': 'CDS',
                'inicio': int(feature.location.start),
                'fin': int(feature.location.end),
                'strand': feature.location.strand,
                'product': feature.qualifiers.get('product', ['Unknown'])[0],
                'translation': feature.qualifiers.get('translation', [''])[0]
            })
    
    return {
        'total_genes': len(genes),
        'total_cds': len(cds_features),
        'genes': genes,
        'cds': cds_features
    }

def calcular_densidad(count, tamano_genoma):
    """
    Calcula densidad por kilobase.
    
    Args:
        count (int): Número de elementos
        tamano_genoma (int): Tamaño del genoma en bp
    
    Returns:
        float: Densidad por kb
    """
    return (count / tamano_genoma) * 1000

# ============================================================================
# ANÁLISIS PRINCIPAL
# ============================================================================

def main():
    print("="*70)
    print("ANÁLISIS DE CODONES ATG - E. coli K-12 MG1655")
    print("="*70)
    
    # Cargar genoma
    print(f"\n[1/5] Cargando genoma desde: {INPUT_FILE}")
    record = SeqIO.read(INPUT_FILE, 'genbank')
    
    print(f"✓ Genoma cargado: {record.id}")
    print(f"✓ Descripción: {record.description[:60]}...")
    print(f"✓ Tamaño: {len(record.seq):,} bp")
    
    # Contar codones ATG
    print("\n[2/5] Buscando codones ATG en el genoma completo...")
    posiciones_atg = contar_codon(record.seq, 'ATG')
    total_atg = len(posiciones_atg)
    
    print(f"✓ Total de codones ATG encontrados: {total_atg:,}")
    
    # Calcular densidad
    print("\n[3/5] Calculando densidad de ATG...")
    densidad_atg = calcular_densidad(total_atg, len(record.seq))
    print(f"✓ Densidad de ATG: {densidad_atg:.2f} por kb")
    
    # Extraer genes anotados
    print("\n[4/5] Extrayendo genes oficialmente anotados...")
    info_genes = extraer_genes_anotados(record)
    
    print(f"✓ Total de genes anotados: {info_genes['total_genes']:,}")
    print(f"✓ Total de CDS anotados: {info_genes['total_cds']:,}")
    
    # Análisis comparativo
    print("\n[5/5] Análisis comparativo...")
    ratio_atg_genes = total_atg / info_genes['total_cds'] if info_genes['total_cds'] > 0 else 0
    
    print(f"✓ Ratio ATG/CDS: {ratio_atg_genes:.2f}")
    print(f"  → Hay ~{ratio_atg_genes:.1f} codones ATG por cada gen anotado")
    
    # ========================================================================
    # GUARDAR RESULTADOS
    # ========================================================================
    
    print("\n" + "="*70)
    print("GUARDANDO RESULTADOS")
    print("="*70)
    
    # Resultado en JSON
    resultado_json = {
        'genoma': {
            'id': record.id,
            'descripcion': record.description,
            'tamano_bp': len(record.seq)
        },
        'analisis_atg': {
            'total_atg': total_atg,
            'densidad_por_kb': round(densidad_atg, 2),
            'primeras_10_posiciones': posiciones_atg[:10],
            'ultimas_10_posiciones': posiciones_atg[-10:]
        },
        'genes_anotados': {
            'total_genes': info_genes['total_genes'],
            'total_cds': info_genes['total_cds']
        },
        'comparacion': {
            'ratio_atg_vs_cds': round(ratio_atg_genes, 2),
            'atg_no_codificantes_estimados': total_atg - info_genes['total_cds']
        }
    }
    
    json_file = os.path.join(OUTPUT_DIR, 'atg_analysis.json')
    with open(json_file, 'w') as f:
        json.dump(resultado_json, f, indent=2)
    print(f"✓ JSON guardado: {json_file}")
    
    # Tabla resumen en CSV
    df_resumen = pd.DataFrame([{
        'Genoma_ID': record.id,
        'Tamaño_bp': len(record.seq),
        'Total_ATG': total_atg,
        'Densidad_ATG_por_kb': round(densidad_atg, 2),
        'Total_Genes': info_genes['total_genes'],
        'Total_CDS': info_genes['total_cds'],
        'Ratio_ATG_CDS': round(ratio_atg_genes, 2)
    }])
    
    csv_file = os.path.join(OUTPUT_DIR, 'atg_summary.csv')
    df_resumen.to_csv(csv_file, index=False)
    print(f"✓ CSV guardado: {csv_file}")
    
    # Tabla de posiciones ATG (primeras 100)
    df_posiciones = pd.DataFrame({
        'Posicion': posiciones_atg[:100],
        'Codon': ['ATG'] * min(100, len(posiciones_atg))
    })
    
    pos_file = os.path.join(OUTPUT_DIR, 'atg_positions_sample.csv')
    df_posiciones.to_csv(pos_file, index=False)
    print(f"✓ Muestra de posiciones guardada: {pos_file}")
    
    # ========================================================================
    # RESUMEN FINAL
    # ========================================================================
    
    print("\n" + "="*70)
    print("RESUMEN FINAL")
    print("="*70)
    print(f"Tamaño del genoma:        {len(record.seq):,} bp")
    print(f"Total codones ATG:        {total_atg:,}")
    print(f"Densidad ATG:             {densidad_atg:.2f} /kb")
    print(f"Total genes anotados:     {info_genes['total_genes']:,}")
    print(f"Total CDS anotados:       {info_genes['total_cds']:,}")
    print(f"Ratio ATG/CDS:            {ratio_atg_genes:.2f}x")
    print("="*70)
    print("\n✅ ANÁLISIS COMPLETADO EXITOSAMENTE\n")

if __name__ == '__main__':
    main()
