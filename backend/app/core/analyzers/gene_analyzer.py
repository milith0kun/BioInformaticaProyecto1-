"""
Gene Analyzer Module
Analyzes gene statistics from genomic data using BioPython

Módulos BioPython utilizados:
- Bio.SeqUtils.gc_fraction: Cálculo de contenido GC por gen
- Bio.SeqFeature: Extracción de información de genes (via genome_parser)
- Bio.SeqRecord: Gestión de registros genómicos

Este módulo implementa el algoritmo de extracción de genes
especificado en la sección 7.2 del informe técnico.

Funcionalidades:
- Número total de genes anotados
- Número de CDS (coding sequences)  
- Longitud de cada gen
- Contenido GC global y por gen
- Posiciones genómicas (start, end, strand)
- Densidad génica (genes/Mb)
- Distribución estadística de tamaños
"""
import numpy as np
from typing import List, Dict, Optional
from dataclasses import dataclass, asdict
from ..genome_parser import GeneData


@dataclass
class GeneStatistics:
    """Estadísticas de tamaño de genes"""
    mean: float
    median: float
    min: int
    max: int
    std: float


@dataclass
class GeneAnalysisResult:
    """Resultado completo del análisis de genes"""
    total_genes: int
    total_cds: int
    genome_length: int
    gc_content: float
    gene_density: float
    size_statistics: GeneStatistics
    genes: List[Dict]


class GeneAnalyzer:
    """
    Analyzer for gene statistics and properties
    
    Implementa el algoritmo 7.2 del informe:
    - Parsear features tipo "gene" y "CDS" del GenBank
    - Extraer atributos: locus_tag, product, location
    - Calcular longitud: end - start + 1
    - Calcular GC: (G + C) / total_bases * 100 usando BioPython
    - Almacenar en estructura tabular
    """
    
    def __init__(self):
        self._last_result: Optional[GeneAnalysisResult] = None
    
    def analyze(self, genes: List[GeneData], genome_length: int, 
                cds_count: int, genome_gc: float) -> GeneAnalysisResult:
        """
        Perform comprehensive gene analysis
        
        Los genes ya vienen procesados por genome_parser.py que usa
        BioPython para extraer la información del archivo GenBank.
        
        Args:
            genes: List of GeneData objects (parsed with BioPython)
            genome_length: Total genome length in bp
            cds_count: Number of CDS features
            genome_gc: Genome-wide GC content (calculated with Bio.SeqUtils)
            
        Returns:
            GeneAnalysisResult with all metrics
        """
        total_genes = len(genes)
        
        # Calculate gene density (genes per megabase)
        gene_density = round((total_genes / genome_length) * 1_000_000, 2) if genome_length > 0 else 0
        
        # Calculate size statistics
        if genes:
            lengths = [g.length for g in genes]
            size_stats = GeneStatistics(
                mean=round(float(np.mean(lengths)), 2),
                median=round(float(np.median(lengths)), 2),
                min=int(np.min(lengths)),
                max=int(np.max(lengths)),
                std=round(float(np.std(lengths)), 2)
            )
        else:
            size_stats = GeneStatistics(mean=0, median=0, min=0, max=0, std=0)
        
        # Convert genes to dictionaries
        genes_dict = [
            {
                "locus_tag": g.locus_tag,
                "start": g.start,
                "end": g.end,
                "length": g.length,
                "strand": g.strand,
                "product": g.product,
                "gc_content": g.gc_content,
                "gene_name": getattr(g, 'gene_name', None),
                "protein_id": getattr(g, 'protein_id', None),
                "start_codon": getattr(g, 'start_codon', None),
                "stop_codon": getattr(g, 'stop_codon', None),
                "has_introns": getattr(g, 'has_introns', False)
            }
            for g in genes
        ]
        
        result = GeneAnalysisResult(
            total_genes=total_genes,
            total_cds=cds_count,
            genome_length=genome_length,
            gc_content=genome_gc,
            gene_density=gene_density,
            size_statistics=size_stats,
            genes=genes_dict
        )
        
        self._last_result = result
        return result
    
    def get_size_distribution(self, genes: List[GeneData], bins: int = 20) -> Dict:
        """
        Get gene size distribution for histogram
        
        Args:
            genes: List of GeneData objects
            bins: Number of histogram bins
            
        Returns:
            Dictionary with bin edges and counts
        """
        if not genes:
            return {"bin_edges": [], "counts": []}
        
        lengths = [g.length for g in genes]
        counts, bin_edges = np.histogram(lengths, bins=bins)
        
        return {
            "bin_edges": [int(b) for b in bin_edges],
            "counts": [int(c) for c in counts]
        }
    
    def get_gc_distribution(self, genes: List[GeneData]) -> Dict:
        """
        Get GC content distribution across genes
        
        Args:
            genes: List of GeneData objects
            
        Returns:
            Dictionary with GC statistics
        """
        if not genes:
            return {"mean": 0, "std": 0, "min": 0, "max": 0, "values": []}
        
        gc_values = [g.gc_content for g in genes]
        
        return {
            "mean": round(float(np.mean(gc_values)), 2),
            "std": round(float(np.std(gc_values)), 2),
            "min": round(float(np.min(gc_values)), 2),
            "max": round(float(np.max(gc_values)), 2),
            "values": gc_values
        }
    
    def get_strand_distribution(self, genes: List[GeneData]) -> Dict:
        """
        Get distribution of genes across strands
        
        Args:
            genes: List of GeneData objects
            
        Returns:
            Dictionary with strand counts
        """
        forward = sum(1 for g in genes if g.strand == 1)
        reverse = sum(1 for g in genes if g.strand == -1)
        
        return {
            "forward": forward,
            "reverse": reverse,
            "forward_percentage": round((forward / len(genes)) * 100, 1) if genes else 0,
            "reverse_percentage": round((reverse / len(genes)) * 100, 1) if genes else 0
        }
    
    def get_genes_by_length_range(self, genes: List[GeneData], 
                                   min_length: int = 0, 
                                   max_length: int = None) -> List[GeneData]:
        """
        Filter genes by length range
        
        Args:
            genes: List of GeneData objects
            min_length: Minimum gene length
            max_length: Maximum gene length (None for no limit)
            
        Returns:
            Filtered list of genes
        """
        if max_length is None:
            return [g for g in genes if g.length >= min_length]
        
        return [g for g in genes if min_length <= g.length <= max_length]
    
    def get_longest_genes(self, genes: List[GeneData], n: int = 10) -> List[GeneData]:
        """Get the n longest genes"""
        return sorted(genes, key=lambda g: g.length, reverse=True)[:n]
    
    def get_shortest_genes(self, genes: List[GeneData], n: int = 10) -> List[GeneData]:
        """Get the n shortest genes"""
        return sorted(genes, key=lambda g: g.length)[:n]
    
    def search_genes(self, genes: List[GeneData], query: str) -> List[GeneData]:
        """
        Search genes by locus tag or product
        
        Args:
            genes: List of GeneData objects
            query: Search query string
            
        Returns:
            Matching genes
        """
        query = query.lower()
        return [
            g for g in genes 
            if query in g.locus_tag.lower() or 
               (g.product and query in g.product.lower())
        ]
    
    def to_dict(self) -> Dict:
        """Convert last result to dictionary"""
        if not self._last_result:
            return {}
        
        return {
            "total_genes": self._last_result.total_genes,
            "total_cds": self._last_result.total_cds,
            "genome_length": self._last_result.genome_length,
            "gc_content": self._last_result.gc_content,
            "gene_density": self._last_result.gene_density,
            "size_statistics": asdict(self._last_result.size_statistics),
            "genes": self._last_result.genes
        }
