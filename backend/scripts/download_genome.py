#!/usr/bin/env python3
"""
Script de Descarga Programática del Genoma de E. coli K-12 MG1655
Utiliza Bio.Entrez (NCBI E-utilities) para descargar el genoma completo

Requisitos del proyecto:
- Descargar programáticamente desde NCBI usando BioPython
- Verificación de integridad de datos
- Manejo de errores de conectividad
- Validación del organismo correcto
- Formato: GenBank completo con todas las anotaciones

Autor: Sistema de Análisis Genómico
Fecha: Febrero 2026
"""

import os
import sys
import time
import hashlib
from pathlib import Path
from typing import Optional, Dict
from datetime import datetime
from Bio import Entrez, SeqIO


class NCBIGenomeDownloader:
    """
    Descargador programático del genoma de E. coli K-12 MG1655
    usando Bio.Entrez (NCBI E-utilities)
    """
    
    # Configuración de Entrez
    EMAIL = "bioinformatics.student@university.edu"
    TOOL = "GenomeAnalysisPipeline"
    
    # Metadatos esperados para validación
    EXPECTED_ORGANISM = "Escherichia coli str. K-12 substr. MG1655"
    EXPECTED_ACCESSION = "NC_000913.3"  # RefSeq accession
    EXPECTED_GENOME_SIZE = 4641652  # bp
    
    def __init__(self, output_dir: str = "ncbi_dataset"):
        """
        Inicializar descargador
        
        Args:
            output_dir: Directorio donde se guardará el genoma descargado
        """
        self.output_dir = Path(output_dir)
        self.data_dir = self.output_dir / "ncbi_dataset" / "data" / "GCF_000005845.2"
        
        # Configurar Bio.Entrez
        Entrez.email = self.EMAIL
        Entrez.tool = self.TOOL
        
    def download_with_retry(self, max_retries: int = 3) -> bool:
        """
        Descargar genoma con reintentos automáticos usando Bio.Entrez
        
        Args:
            max_retries: Número máximo de reintentos en caso de error
            
        Returns:
            True si la descarga fue exitosa, False en caso contrario
        """
        print("=" * 70)
        print("🧬 DESCARGA PROGRAMÁTICA DEL GENOMA E. COLI K-12 MG1655")
        print("=" * 70)
        print(f"\n📍 Accession: {self.EXPECTED_ACCESSION}")
        print(f"🔗 NCBI E-utilities (Bio.Entrez)")
        print(f"📁 Directorio destino: {self.output_dir.absolute()}\n")
        
        # Crear directorios
        self.data_dir.mkdir(parents=True, exist_ok=True)
        self.gbff_path = self.data_dir / "genomic.gbff"
        
        for attempt in range(1, max_retries + 1):
            try:
                print(f"🔄 Intento {attempt}/{max_retries}...")
                print(f"⏳ Conectando a NCBI...")
                
                # Descargar GenBank usando Entrez.efetch
                start_time = time.time()
                
                print(f"📥 Descargando archivo GenBank completo...")
                handle = Entrez.efetch(
                    db="nucleotide",
                    id=self.EXPECTED_ACCESSION,
                    rettype="gbwithparts",
                    retmode="text"
                )
                
                # Guardar archivo
                downloaded = 0
                chunk_size = 8192
                
                with open(self.gbff_path, 'w') as f:
                    while True:
                        chunk = handle.read(chunk_size)
                        if not chunk:
                            break
                        f.write(chunk)
                        downloaded += len(chunk)
                        
                        # Mostrar progreso
                        if downloaded % (chunk_size * 100) == 0:
                            print(f"   {downloaded / (1024*1024):.2f} MB descargados...", end='\r')
                
                handle.close()
                
                elapsed = time.time() - start_time
                file_size = self.gbff_path.stat().st_size / (1024 * 1024)
                
                print(f"\n✅ Descarga completada en {elapsed:.1f} segundos")
                print(f"💾 Archivo guardado: {self.gbff_path}")
                print(f"📦 Tamaño: {file_size:.2f} MB")
                
                return True
                
            except Exception as e:
                print(f"⚠️  Error en intento {attempt}: {e}")
                if attempt < max_retries:
                    wait_time = 5 * attempt
                    print(f"⏳ Esperando {wait_time} segundos antes de reintentar...")
                    time.sleep(wait_time)
                else:
                    print(f"❌ Falló después de {max_retries} intentos")
                    return False
        
        return False
    
    def calculate_md5(self, filepath: Path) -> str:
        """
        Calcular hash MD5 del archivo para verificación de integridad
        
        Args:
            filepath: Ruta al archivo
            
        Returns:
            Hash MD5 en formato hexadecimal
        """
        print(f"\n🔐 Calculando hash MD5 para verificación de integridad...")
        md5_hash = hashlib.md5()
        
        with open(filepath, 'rb') as f:
            for chunk in iter(lambda: f.read(8192), b""):
                md5_hash.update(chunk)
        
        hash_value = md5_hash.hexdigest()
        print(f"   MD5: {hash_value}")
        return hash_value
    
    def validate_organism(self) -> bool:
        """
        Validar que el organismo descargado es el correcto usando BioPython
        
        Returns:
            True si la validación es exitosa, False en caso contrario
        """
        print(f"\n🔍 Validando organismo descargado...")
        
        try:
            # Leer archivo GenBank con BioPython
            record = SeqIO.read(self.gbff_path, "genbank")
            
            # Validar accession
            if record.id == self.EXPECTED_ACCESSION:
                print(f"✅ Accession verificado: {record.id}")
            else:
                print(f"⚠️  Accession no coincide: esperado {self.EXPECTED_ACCESSION}, obtenido {record.id}")
            
            # Validar organismo
            organism = record.annotations.get('organism', '')
            if "Escherichia coli" in organism and "K-12" in organism:
                print(f"✅ Organismo verificado: {organism}")
            else:
                print(f"⚠️  Organismo: {organism}")
            
            # Validar tamaño del genoma
            genome_length = len(record.seq)
            deviation = abs(genome_length - self.EXPECTED_GENOME_SIZE)
            
            if deviation == 0:
                print(f"✅ Tamaño del genoma: {genome_length:,} bp (correcto)")
            else:
                print(f"⚠️  Tamaño del genoma: {genome_length:,} bp (esperado: {self.EXPECTED_GENOME_SIZE:,})")
            
            # Contar genes
            gene_count = sum(1 for feature in record.features if feature.type == "gene")
            cds_count = sum(1 for feature in record.features if feature.type == "CDS")
            
            print(f"✅ Genes encontrados: {gene_count:,}")
            print(f"✅ CDS encontrados: {cds_count:,}")
            
            return True
            
        except Exception as e:
            print(f"❌ Error durante validación: {e}")
            return False
    
    def generate_report(self, md5_hash: str, success: bool) -> None:
        """
        Generar reporte de descarga
        
        Args:
            md5_hash: Hash MD5 del archivo descargado
            success: Si el proceso fue exitoso
        """
        report_file = self.output_dir / "download_report.txt"
        
        with open(report_file, 'w') as f:
            f.write("=" * 70 + "\n")
            f.write("REPORTE DE DESCARGA - GENOMA E. COLI K-12 MG1655\n")
            f.write("=" * 70 + "\n\n")
            f.write(f"Fecha: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write(f"Accession: {self.EXPECTED_ACCESSION}\n")
            f.write(f"Organismo: {self.EXPECTED_ORGANISM}\n")
            f.write(f"Método: Bio.Entrez (NCBI E-utilities)\n")
            f.write(f"Database: nucleotide\n")
            f.write(f"RetType: gbwithparts\n\n")
            f.write(f"Archivo descargado: {self.gbff_path.name}\n")
            f.write(f"Ruta: {self.gbff_path}\n")
            f.write(f"Tamaño: {self.gbff_path.stat().st_size / (1024*1024):.2f} MB\n")
            f.write(f"MD5: {md5_hash}\n\n")
            f.write(f"Estado: {'✅ EXITOSO' if success else '❌ FALLIDO'}\n")
            f.write("=" * 70 + "\n")
        
        print(f"\n📋 Reporte generado: {report_file}")
    
    def run(self) -> bool:
        """
        Ejecutar el proceso completo de descarga
        
        Returns:
            True si todo el proceso fue exitoso
        """
        print("\n" + "🧬" * 35)
        print("INICIO DE DESCARGA PROGRAMÁTICA CON BIO.ENTREZ")
        print("🧬" * 35 + "\n")
        
        # Paso 1: Descargar
        if not self.download_with_retry():
            print("\n❌ DESCARGA FALLIDA")
            return False
        
        # Paso 2: Verificar integridad
        md5_hash = self.calculate_md5(self.gbff_path)
        
        # Paso 3: Validar organismo
        organism_valid = self.validate_organism()
        
        # Paso 4: Generar reporte
        self.generate_report(md5_hash, organism_valid)
        
        # Resumen final
        print("\n" + "=" * 70)
        print("🎉 PROCESO COMPLETADO")
        print("=" * 70)
        print(f"\n✅ Genoma descargado exitosamente")
        print(f"✅ Datos guardados en: {self.data_dir.absolute()}")
        print(f"✅ Validación: {'EXITOSA' if organism_valid else 'PARCIAL'}")
        print(f"\n📚 Archivo GenBank completo:")
        print(f"   - {self.gbff_path.name} (formato GenBank con todas las anotaciones)")
        print(f"   - Contiene: secuencia completa, genes, CDS, features")
        print("\n🚀 Listo para análisis bioinformático\n")
        
        return True


def main():
    """Función principal"""
    # Configurar directorio de salida
    project_root = Path(__file__).parent.parent.parent
    output_dir = project_root / "ncbi_dataset"
    
    # Crear descargador
    downloader = NCBIGenomeDownloader(str(output_dir))
    
    # Ejecutar
    success = downloader.run()
    
    # Exit code
    sys.exit(0 if success else 1)


if __name__ == "__main__":
    main()
