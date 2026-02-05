#!/usr/bin/env python3
"""
Análisis de Codones de Terminación en E. coli K-12 MG1655
=========================================================
Este script identifica y cuantifica los tres codones de terminación
(TAA, TAG, TGA) presentes en el genoma completo, calculando sus
proporciones relativas y densidad.

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

STOP_CODONS = ['TAA', 'TAG', 'TGA']

# ============================================================================
# FUNCIONES
# ============================================================================

def contar_codon(secuencia, codon):
    """
    Cuenta todas las apariciones de un codón en una secuencia.
    
    Args:
        secuencia (str): Secuencia de ADN
        codon (str): Codón a buscar (ej: 'TAA')
    
    Returns:
        list: Lista de posiciones donde aparece el codón (0-indexed)
    """
    secuencia = str(secuencia).upper()
    codon = codon.upper()
    posiciones = []
    
    # Buscar en toda la secuencia
    for i in range(len(secuencia) - len(codon) + 1):
        if secuencia[i:i+len(codon)] == codon:
            posiciones.append(i)
    
    return posiciones

def extraer_stop_codons_cds(record):
    """
    Extrae los codones de terminación de los CDS anotados.
    
    Args:
        record: Objeto SeqRecord de BioPython
    
    Returns:
        dict: Diccionario con stop codons por CDS
    """
    stop_codons_cds = []
    
    for feature in record.features:
        if feature.type == 'CDS':
            # Obtener secuencia del CDS
            cds_seq = feature.extract(record.seq)
            
            # El stop codon debería estar al final (últimos 3 nucleótidos)
            if len(cds_seq) >= 3:
                stop_codon = str(cds_seq[-3:]).upper()
                
                stop_codons_cds.append({
                    'producto': feature.qualifiers.get('product', ['Unknown'])[0],
                    'inicio': int(feature.location.start),
                    'fin': int(feature.location.end),
                    'stop_codon': stop_codon,
                    'es_stop_valido': stop_codon in STOP_CODONS
                })
    
    return stop_codons_cds

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

def calcular_proporciones(conteos):
    """
    Calcula las proporciones relativas de cada stop codon.
    
    Args:
        conteos (dict): Diccionario con conteos de cada stop codon
    
    Returns:
        dict: Proporciones en porcentaje
    """
    total = sum(conteos.values())
    if total == 0:
        return {k: 0.0 for k in conteos.keys()}
    
    return {k: (v / total) * 100 for k, v in conteos.items()}

# ============================================================================
# ANÁLISIS PRINCIPAL
# ============================================================================

def main():
    print("="*70)
    print("ANÁLISIS DE CODONES DE TERMINACIÓN - E. coli K-12 MG1655")
    print("="*70)
    
    # Cargar genoma
    print(f"\n[1/6] Cargando genoma desde: {INPUT_FILE}")
    record = SeqIO.read(INPUT_FILE, 'genbank')
    
    print(f"✓ Genoma cargado: {record.id}")
    print(f"✓ Descripción: {record.description[:60]}...")
    print(f"✓ Tamaño: {len(record.seq):,} bp")
    
    # Contar cada stop codon en todo el genoma
    print(f"\n[2/6] Buscando codones de terminación en el genoma completo...")
    print(f"   Codones a buscar: {', '.join(STOP_CODONS)}")
    
    resultados_stop = {}
    posiciones_stop = {}
    
    for stop_codon in STOP_CODONS:
        print(f"\n   → Buscando {stop_codon}...", end=" ")
        posiciones = contar_codon(record.seq, stop_codon)
        total = len(posiciones)
        densidad = calcular_densidad(total, len(record.seq))
        
        resultados_stop[stop_codon] = {
            'total': total,
            'densidad_por_kb': densidad,
            'posiciones': posiciones
        }
        posiciones_stop[stop_codon] = posiciones
        
        print(f"✓ {total:,} encontrados ({densidad:.2f} /kb)")
    
    # Calcular totales
    print(f"\n[3/6] Calculando estadísticas globales...")
    total_stops = sum(r['total'] for r in resultados_stop.values())
    densidad_total = calcular_densidad(total_stops, len(record.seq))
    
    print(f"✓ Total de stop codons: {total_stops:,}")
    print(f"✓ Densidad total: {densidad_total:.2f} /kb")
    
    # Calcular proporciones
    print(f"\n[4/6] Calculando proporciones relativas...")
    conteos = {k: v['total'] for k, v in resultados_stop.items()}
    proporciones = calcular_proporciones(conteos)
    
    for stop_codon in STOP_CODONS:
        print(f"✓ {stop_codon}: {proporciones[stop_codon]:.2f}% ({conteos[stop_codon]:,} codones)")
    
    # Analizar stop codons en CDS anotados
    print(f"\n[5/6] Analizando stop codons en CDS anotados...")
    stop_cds = extraer_stop_codons_cds(record)
    
    # Contar stop codons en CDS
    stop_cds_counter = Counter([s['stop_codon'] for s in stop_cds if s['es_stop_valido']])
    stop_cds_proporciones = calcular_proporciones(stop_cds_counter)
    
    print(f"✓ Total CDS analizados: {len(stop_cds)}")
    print(f"✓ Stop codons en CDS:")
    for stop_codon in STOP_CODONS:
        count = stop_cds_counter.get(stop_codon, 0)
        prop = stop_cds_proporciones.get(stop_codon, 0.0)
        print(f"   → {stop_codon}: {count:,} ({prop:.2f}%)")
    
    # Comparación genoma vs CDS
    print(f"\n[6/6] Comparación: Genoma completo vs CDS anotados...")
    print("\n   Codón  | Genoma (%) | CDS (%)  | Diferencia")
    print("   " + "-"*50)
    for stop_codon in STOP_CODONS:
        gen_prop = proporciones[stop_codon]
        cds_prop = stop_cds_proporciones.get(stop_codon, 0.0)
        diff = gen_prop - cds_prop
        print(f"   {stop_codon}    |   {gen_prop:5.2f}   |  {cds_prop:5.2f}  |  {diff:+6.2f}")
    
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
        'analisis_stop_codons': {
            'total_stop_codons': total_stops,
            'densidad_total_por_kb': round(densidad_total, 2),
            'por_tipo': {
                stop: {
                    'total': resultados_stop[stop]['total'],
                    'densidad_por_kb': round(resultados_stop[stop]['densidad_por_kb'], 2),
                    'proporcion_pct': round(proporciones[stop], 2),
                    'primeras_10_posiciones': resultados_stop[stop]['posiciones'][:10],
                    'ultimas_10_posiciones': resultados_stop[stop]['posiciones'][-10:]
                }
                for stop in STOP_CODONS
            }
        },
        'analisis_cds': {
            'total_cds': len(stop_cds),
            'stop_codons_por_tipo': {
                stop: {
                    'total': stop_cds_counter.get(stop, 0),
                    'proporcion_pct': round(stop_cds_proporciones.get(stop, 0.0), 2)
                }
                for stop in STOP_CODONS
            }
        },
        'comparacion': {
            'genoma_vs_cds': {
                stop: {
                    'genoma_pct': round(proporciones[stop], 2),
                    'cds_pct': round(stop_cds_proporciones.get(stop, 0.0), 2),
                    'diferencia': round(proporciones[stop] - stop_cds_proporciones.get(stop, 0.0), 2)
                }
                for stop in STOP_CODONS
            }
        }
    }
    
    json_file = os.path.join(OUTPUT_DIR, 'stop_codons_analysis.json')
    with open(json_file, 'w') as f:
        json.dump(resultado_json, f, indent=2)
    print(f"✓ JSON guardado: {json_file}")
    
    # Tabla resumen en CSV
    df_resumen = pd.DataFrame([
        {
            'Stop_Codon': stop,
            'Total_Genoma': resultados_stop[stop]['total'],
            'Densidad_por_kb': round(resultados_stop[stop]['densidad_por_kb'], 2),
            'Proporcion_Genoma_%': round(proporciones[stop], 2),
            'Total_CDS': stop_cds_counter.get(stop, 0),
            'Proporcion_CDS_%': round(stop_cds_proporciones.get(stop, 0.0), 2)
        }
        for stop in STOP_CODONS
    ])
    
    csv_file = os.path.join(OUTPUT_DIR, 'stop_codons_summary.csv')
    df_resumen.to_csv(csv_file, index=False)
    print(f"✓ CSV guardado: {csv_file}")
    
    # Tabla detallada de stop codons en CDS
    df_cds = pd.DataFrame(stop_cds[:100])  # Primeros 100
    cds_file = os.path.join(OUTPUT_DIR, 'stop_codons_cds_sample.csv')
    df_cds.to_csv(cds_file, index=False)
    print(f"✓ Muestra CDS guardada: {cds_file}")
    
    # ========================================================================
    # RESUMEN FINAL
    # ========================================================================
    
    print("\n" + "="*70)
    print("RESUMEN FINAL - STOP CODONS")
    print("="*70)
    print(f"Tamaño del genoma:           {len(record.seq):,} bp")
    print(f"Total stop codons (genoma):  {total_stops:,}")
    print(f"Densidad total:              {densidad_total:.2f} /kb")
    print(f"\nProporciones en genoma completo:")
    for stop in STOP_CODONS:
        print(f"  {stop}: {proporciones[stop]:5.2f}% ({conteos[stop]:,} codones)")
    print(f"\nProporciones en CDS anotados:")
    for stop in STOP_CODONS:
        count = stop_cds_counter.get(stop, 0)
        prop = stop_cds_proporciones.get(stop, 0.0)
        print(f"  {stop}: {prop:5.2f}% ({count:,} CDS)")
    print("="*70)
    print("\n✅ ANÁLISIS DE STOP CODONS COMPLETADO EXITOSAMENTE\n")

if __name__ == '__main__':
    main()
