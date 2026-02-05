"""
Genome Comparator Module
Compara múltiples genomas secuencialmente para análisis comparativo

Funcionalidades:
- Comparación de genomas descargados
- Identificación de genes únicos y compartidos
- Análisis de genes de mayor y menor tamaño por genoma
- Comparación de estadísticas entre cepas
- Búsqueda individual y por grupos funcionales
- Soporte para cepas bacterianas emparentadas

Autor: Sistema de Análisis Genómico
Fecha: Febrero 2026
"""

import os
import json
from typing import List, Dict, Optional, Any, Tuple
from dataclasses import dataclass, asdict, field
from datetime import datetime
import numpy as np
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction

from ..genome_parser import GenomeParser, GeneData


@dataclass
class GenomeStats:
    """Estadísticas de un genoma individual"""
    accession: str
    organism_name: str
    genome_length: int
    gene_count: int
    cds_count: int
    gc_content: float
    gene_density: float  # genes/Mb
    avg_gene_length: float
    min_gene_length: int
    max_gene_length: int


@dataclass
class GeneGroup:
    """Grupo de genes por categoría funcional"""
    name: str
    description: str
    genes: List[Dict]
    count: int
    avg_length: float


@dataclass
class GenomeComparisonResult:
    """Resultado de comparación entre genomas"""
    genomes: List[GenomeStats]
    comparison_date: str
    total_genomes_compared: int
    
    # Estadísticas comparativas
    largest_genome: str  # accession
    smallest_genome: str
    highest_gc: str
    lowest_gc: str
    highest_gene_density: str
    lowest_gene_density: str
    
    # Genes extremos globales
    longest_genes_global: List[Dict]  # Top 10 genes más largos de todos los genomas
    shortest_genes_global: List[Dict]  # Top 10 genes más cortos de todos los genomas
    
    # Comparación de métricas
    metrics_comparison: Dict[str, Dict]  # métrica -> {genoma: valor}


@dataclass
class GeneFilterResult:
    """Resultado de filtrado de genes"""
    genes: List[Dict]
    total_count: int
    filter_applied: str
    genome_source: str


# Cepas de E. coli emparentadas para comparación rápida
RELATED_ECOLI_STRAINS = [
    {
        "accession": "GCF_000005845.2",
        "organism": "Escherichia coli str. K-12 substr. MG1655",
        "strain": "K-12 MG1655",
        "description": "Cepa de laboratorio de referencia, genoma completamente secuenciado",
        "category": "laboratory",
        "genome_size_mb": 4.6,
        "is_reference": True
    },
    {
        "accession": "GCF_000008865.2",
        "organism": "Escherichia coli O157:H7 str. Sakai",
        "strain": "O157:H7",
        "description": "Cepa patógena causante de colitis hemorrágica",
        "category": "pathogenic",
        "genome_size_mb": 5.5,
        "is_reference": True
    },
    {
        "accession": "GCF_000009565.2",
        "organism": "Escherichia coli BL21(DE3)",
        "strain": "BL21(DE3)",
        "description": "Cepa de expresión de proteínas, usada en biotecnología",
        "category": "industrial",
        "genome_size_mb": 4.6,
        "is_reference": False
    },
    {
        "accession": "GCF_000019425.1",
        "organism": "Escherichia coli CFT073",
        "strain": "CFT073",
        "description": "Cepa uropatógena, causa infecciones del tracto urinario",
        "category": "pathogenic",
        "genome_size_mb": 5.2,
        "is_reference": True
    },
    {
        "accession": "GCF_000026305.1",
        "organism": "Escherichia coli O104:H4 str. 2011C-3493",
        "strain": "O104:H4",
        "description": "Cepa del brote europeo 2011, altamente virulenta",
        "category": "pathogenic",
        "genome_size_mb": 5.3,
        "is_reference": False
    },
    {
        "accession": "GCF_000750555.1",
        "organism": "Escherichia coli Nissle 1917",
        "strain": "Nissle 1917",
        "description": "Cepa probiótica, usada en tratamientos médicos",
        "category": "probiotic",
        "genome_size_mb": 5.0,
        "is_reference": False
    },
    {
        "accession": "GCF_000007445.1",
        "organism": "Escherichia coli W3110",
        "strain": "W3110",
        "description": "Cepa derivada de K-12, usada en investigación",
        "category": "laboratory",
        "genome_size_mb": 4.6,
        "is_reference": False
    }
]

# Grupos funcionales de genes para filtrado
GENE_FUNCTIONAL_GROUPS = {
    "metabolism": {
        "name": "Metabolismo",
        "keywords": ["dehydrogenase", "synthase", "kinase", "transferase", "oxidase", "reductase"],
        "description": "Genes involucrados en procesos metabólicos"
    },
    "transport": {
        "name": "Transporte",
        "keywords": ["transporter", "permease", "channel", "pump", "ABC", "import", "export"],
        "description": "Genes de transporte de membrana"
    },
    "regulation": {
        "name": "Regulación",
        "keywords": ["regulator", "repressor", "activator", "transcription", "response"],
        "description": "Genes reguladores de expresión"
    },
    "dna_rna": {
        "name": "DNA/RNA",
        "keywords": ["polymerase", "helicase", "ligase", "primase", "topoisomerase", "nuclease"],
        "description": "Genes de replicación, reparación y transcripción"
    },
    "protein": {
        "name": "Proteínas",
        "keywords": ["ribosomal", "chaperone", "protease", "peptidase", "aminoacyl"],
        "description": "Genes de síntesis y procesamiento de proteínas"
    },
    "cell_structure": {
        "name": "Estructura Celular",
        "keywords": ["membrane", "wall", "flagell", "pilus", "fimbr", "capsule"],
        "description": "Genes de estructuras celulares"
    },
    "stress_response": {
        "name": "Respuesta a Estrés",
        "keywords": ["stress", "shock", "heat", "cold", "oxidative", "SOS"],
        "description": "Genes de respuesta a condiciones adversas"
    },
    "hypothetical": {
        "name": "Hipotéticos",
        "keywords": ["hypothetical", "unknown", "uncharacterized", "predicted"],
        "description": "Genes de función desconocida"
    }
}


class GenomeComparator:
    """
    Comparador de genomas para análisis secuencial
    
    Permite comparar múltiples genomas descargados,
    identificar genes extremos y agrupar por función.
    """
    
    def __init__(self, cache_dir: str = ".cache"):
        self.cache_dir = cache_dir
        self._genome_cache: Dict[str, Dict] = {}
        self._parser = GenomeParser(cache_dir)
    
    def get_related_strains(self, category: Optional[str] = None) -> List[Dict]:
        """
        Obtener lista de cepas de E. coli emparentadas
        
        Args:
            category: Filtrar por categoría (laboratory, pathogenic, industrial, probiotic)
            
        Returns:
            Lista de cepas relacionadas
        """
        if category:
            return [s for s in RELATED_ECOLI_STRAINS if s["category"] == category]
        return RELATED_ECOLI_STRAINS
    
    def load_genome_from_directory(self, genome_dir: str) -> Optional[Dict]:
        """
        Cargar datos de un genoma desde su directorio
        
        Args:
            genome_dir: Directorio del genoma (ej: ncbi_dataset/ncbi_dataset/data/GCF_000005845.2)
            
        Returns:
            Datos del genoma o None si no se encuentra
        """
        accession = os.path.basename(genome_dir)
        
        # Buscar archivo GenBank
        gbff_file = None
        for f in os.listdir(genome_dir):
            if f.endswith('.gbff'):
                gbff_file = os.path.join(genome_dir, f)
                break
        
        if not gbff_file:
            return None
        
        # Parsear genoma
        try:
            genome_data = self._parser.parse_genbank(gbff_file)
            
            return {
                "accession": accession,
                "organism_name": genome_data.organism,
                "genome_length": genome_data.length,
                "gc_content": genome_data.gc_content,
                "gene_count": len(genome_data.genes),
                "cds_count": genome_data.cds_count,
                "genes": genome_data.genes,
                "sequence": genome_data.sequence
            }
        except Exception as e:
            print(f"Error al cargar genoma {accession}: {e}")
            return None
    
    def scan_downloaded_genomes(self, base_path: str) -> List[str]:
        """
        Escanear genomas descargados en el directorio
        
        Args:
            base_path: Ruta base del proyecto
            
        Returns:
            Lista de rutas a directorios de genomas
        """
        genomes = []
        ncbi_data_path = os.path.join(base_path, "ncbi_dataset", "ncbi_dataset", "data")
        
        if not os.path.exists(ncbi_data_path):
            return genomes
        
        for item in os.listdir(ncbi_data_path):
            item_path = os.path.join(ncbi_data_path, item)
            if os.path.isdir(item_path) and (item.startswith("GCF_") or item.startswith("GCA_")):
                genomes.append(item_path)
        
        return genomes
    
    def compare_genomes(self, genome_dirs: List[str]) -> GenomeComparisonResult:
        """
        Comparar múltiples genomas secuencialmente
        
        Args:
            genome_dirs: Lista de directorios de genomas a comparar
            
        Returns:
            Resultado de la comparación
        """
        genome_stats = []
        all_genes = []
        
        for genome_dir in genome_dirs:
            genome_data = self.load_genome_from_directory(genome_dir)
            if not genome_data:
                continue
            
            # Calcular estadísticas
            genes = genome_data["genes"]
            lengths = [g.length for g in genes]
            
            stats = GenomeStats(
                accession=genome_data["accession"],
                organism_name=genome_data["organism_name"],
                genome_length=genome_data["genome_length"],
                gene_count=genome_data["gene_count"],
                cds_count=genome_data["cds_count"],
                gc_content=genome_data["gc_content"],
                gene_density=round((genome_data["gene_count"] / genome_data["genome_length"]) * 1_000_000, 2),
                avg_gene_length=round(np.mean(lengths), 2) if lengths else 0,
                min_gene_length=min(lengths) if lengths else 0,
                max_gene_length=max(lengths) if lengths else 0
            )
            genome_stats.append(stats)
            
            # Agregar genes con referencia al genoma
            for gene in genes:
                all_genes.append({
                    "locus_tag": gene.locus_tag,
                    "product": gene.product,
                    "length": gene.length,
                    "gc_content": gene.gc_content,
                    "start": gene.start,
                    "end": gene.end,
                    "strand": gene.strand,
                    "genome": genome_data["accession"]
                })
        
        if not genome_stats:
            return None
        
        # Ordenar genes por longitud
        sorted_by_length = sorted(all_genes, key=lambda g: g["length"], reverse=True)
        
        # Encontrar extremos
        largest = max(genome_stats, key=lambda s: s.genome_length)
        smallest = min(genome_stats, key=lambda s: s.genome_length)
        highest_gc = max(genome_stats, key=lambda s: s.gc_content)
        lowest_gc = min(genome_stats, key=lambda s: s.gc_content)
        highest_density = max(genome_stats, key=lambda s: s.gene_density)
        lowest_density = min(genome_stats, key=lambda s: s.gene_density)
        
        # Construir comparación de métricas
        metrics_comparison = {
            "genome_length": {s.accession: s.genome_length for s in genome_stats},
            "gene_count": {s.accession: s.gene_count for s in genome_stats},
            "gc_content": {s.accession: s.gc_content for s in genome_stats},
            "gene_density": {s.accession: s.gene_density for s in genome_stats},
            "avg_gene_length": {s.accession: s.avg_gene_length for s in genome_stats}
        }
        
        return GenomeComparisonResult(
            genomes=genome_stats,
            comparison_date=datetime.now().isoformat(),
            total_genomes_compared=len(genome_stats),
            largest_genome=largest.accession,
            smallest_genome=smallest.accession,
            highest_gc=highest_gc.accession,
            lowest_gc=lowest_gc.accession,
            highest_gene_density=highest_density.accession,
            lowest_gene_density=lowest_density.accession,
            longest_genes_global=sorted_by_length[:10],
            shortest_genes_global=sorted_by_length[-10:] if len(sorted_by_length) >= 10 else sorted_by_length,
            metrics_comparison=metrics_comparison
        )
    
    def get_genes_by_size(
        self,
        genes: List[GeneData],
        order: str = "largest",
        count: int = 10,
        genome_accession: str = ""
    ) -> GeneFilterResult:
        """
        Obtener genes ordenados por tamaño
        
        Args:
            genes: Lista de genes
            order: "largest" o "smallest"
            count: Número de genes a retornar
            genome_accession: Accession del genoma fuente
            
        Returns:
            Resultado con genes filtrados
        """
        if order == "largest":
            sorted_genes = sorted(genes, key=lambda g: g.length, reverse=True)[:count]
            filter_desc = f"Top {count} genes más largos"
        else:
            sorted_genes = sorted(genes, key=lambda g: g.length)[:count]
            filter_desc = f"Top {count} genes más cortos"
        
        return GeneFilterResult(
            genes=[{
                "locus_tag": g.locus_tag,
                "product": g.product or "Sin descripción",
                "length": g.length,
                "gc_content": g.gc_content,
                "start": g.start,
                "end": g.end,
                "strand": "+" if g.strand == 1 else "-"
            } for g in sorted_genes],
            total_count=len(sorted_genes),
            filter_applied=filter_desc,
            genome_source=genome_accession
        )
    
    def filter_genes_by_group(
        self,
        genes: List[GeneData],
        group_id: str,
        genome_accession: str = ""
    ) -> GeneFilterResult:
        """
        Filtrar genes por grupo funcional
        
        Args:
            genes: Lista de genes
            group_id: ID del grupo funcional (metabolism, transport, etc.)
            genome_accession: Accession del genoma fuente
            
        Returns:
            Resultado con genes del grupo
        """
        if group_id not in GENE_FUNCTIONAL_GROUPS:
            return GeneFilterResult(
                genes=[],
                total_count=0,
                filter_applied=f"Grupo desconocido: {group_id}",
                genome_source=genome_accession
            )
        
        group = GENE_FUNCTIONAL_GROUPS[group_id]
        keywords = group["keywords"]
        
        # Filtrar genes por palabras clave en el producto
        filtered = []
        for gene in genes:
            if gene.product:
                product_lower = gene.product.lower()
                if any(kw.lower() in product_lower for kw in keywords):
                    filtered.append(gene)
        
        return GeneFilterResult(
            genes=[{
                "locus_tag": g.locus_tag,
                "product": g.product or "Sin descripción",
                "length": g.length,
                "gc_content": g.gc_content,
                "start": g.start,
                "end": g.end,
                "strand": "+" if g.strand == 1 else "-"
            } for g in filtered],
            total_count=len(filtered),
            filter_applied=f"Grupo: {group['name']} - {group['description']}",
            genome_source=genome_accession
        )
    
    def search_genes_advanced(
        self,
        genes: List[GeneData],
        query: str = "",
        min_length: Optional[int] = None,
        max_length: Optional[int] = None,
        min_gc: Optional[float] = None,
        max_gc: Optional[float] = None,
        strand: Optional[str] = None,
        genome_accession: str = ""
    ) -> GeneFilterResult:
        """
        Búsqueda avanzada de genes con múltiples filtros
        
        Args:
            genes: Lista de genes
            query: Texto de búsqueda (locus_tag o producto)
            min_length: Longitud mínima
            max_length: Longitud máxima
            min_gc: GC mínimo
            max_gc: GC máximo
            strand: Hebra ('+' o '-')
            genome_accession: Accession del genoma fuente
            
        Returns:
            Resultado con genes filtrados
        """
        filtered = genes.copy()
        filters_applied = []
        
        # Filtrar por texto
        if query:
            query_lower = query.lower()
            filtered = [
                g for g in filtered
                if query_lower in g.locus_tag.lower() or 
                   (g.product and query_lower in g.product.lower())
            ]
            filters_applied.append(f"Texto: '{query}'")
        
        # Filtrar por longitud
        if min_length is not None:
            filtered = [g for g in filtered if g.length >= min_length]
            filters_applied.append(f"Longitud ≥ {min_length}")
        
        if max_length is not None:
            filtered = [g for g in filtered if g.length <= max_length]
            filters_applied.append(f"Longitud ≤ {max_length}")
        
        # Filtrar por GC
        if min_gc is not None:
            filtered = [g for g in filtered if g.gc_content >= min_gc]
            filters_applied.append(f"GC ≥ {min_gc}%")
        
        if max_gc is not None:
            filtered = [g for g in filtered if g.gc_content <= max_gc]
            filters_applied.append(f"GC ≤ {max_gc}%")
        
        # Filtrar por hebra
        if strand:
            strand_value = 1 if strand == "+" else -1
            filtered = [g for g in filtered if g.strand == strand_value]
            filters_applied.append(f"Hebra: {strand}")
        
        return GeneFilterResult(
            genes=[{
                "locus_tag": g.locus_tag,
                "product": g.product or "Sin descripción",
                "length": g.length,
                "gc_content": g.gc_content,
                "start": g.start,
                "end": g.end,
                "strand": "+" if g.strand == 1 else "-"
            } for g in filtered],
            total_count=len(filtered),
            filter_applied=" | ".join(filters_applied) if filters_applied else "Sin filtros",
            genome_source=genome_accession
        )
    
    def get_gene_groups_summary(
        self,
        genes: List[GeneData],
        genome_accession: str = ""
    ) -> List[Dict]:
        """
        Obtener resumen de genes por grupo funcional
        
        Args:
            genes: Lista de genes
            genome_accession: Accession del genoma fuente
            
        Returns:
            Lista de resúmenes por grupo
        """
        summaries = []
        
        for group_id, group_info in GENE_FUNCTIONAL_GROUPS.items():
            keywords = group_info["keywords"]
            matching_genes = []
            
            for gene in genes:
                if gene.product:
                    product_lower = gene.product.lower()
                    if any(kw.lower() in product_lower for kw in keywords):
                        matching_genes.append(gene)
            
            if matching_genes:
                lengths = [g.length for g in matching_genes]
                avg_length = round(np.mean(lengths), 2)
            else:
                avg_length = 0
            
            summaries.append({
                "id": group_id,
                "name": group_info["name"],
                "description": group_info["description"],
                "count": len(matching_genes),
                "avg_length": avg_length,
                "percentage": round((len(matching_genes) / len(genes)) * 100, 2) if genes else 0
            })
        
        # Ordenar por cantidad de genes
        summaries.sort(key=lambda x: x["count"], reverse=True)
        return summaries
    
    def get_functional_groups(self) -> Dict:
        """Obtener definiciones de grupos funcionales"""
        return GENE_FUNCTIONAL_GROUPS


def get_genome_comparator(cache_dir: str = ".cache") -> GenomeComparator:
    """Factory function para obtener instancia del comparador"""
    return GenomeComparator(cache_dir)
