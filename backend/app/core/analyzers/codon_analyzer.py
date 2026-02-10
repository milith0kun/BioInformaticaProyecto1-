"""
Codon Analyzer Module
Analyzes start and stop codons in genomic sequences using BioPython

Módulos BioPython utilizados:
- Bio.Seq: Representación de secuencias como objetos Seq
- Bio.Seq.count(): Conteo de patrones (codones ATG, TAA, TAG, TGA)
- Bio.Seq.find(): Búsqueda de posiciones de codones

Este módulo implementa el algoritmo de conteo de codones especificado
en la sección 7.1 del informe técnico.
"""
from Bio.Seq import Seq
from typing import Dict, Tuple, List, Union
from dataclasses import dataclass
import numpy as np


@dataclass
class CodonAnalysisResult:
    genome_length: int
    atg_count: int
    atg_density: float
    stop_codons: Dict[str, Dict]
    gene_comparison: Dict[str, int]
    spatial_distribution: Dict = None  # Nueva: distribución espacial
    statistical_quality: Dict = None   # Nueva: métricas de calidad


class CodonAnalyzer:
    """
    Analyzer for start (ATG) and stop (TAA, TAG, TGA) codons
    
    Utiliza BioPython Bio.Seq para:
    - Representación eficiente de secuencias de DNA
    - Conteo optimizado de codones usando seq.count()
    - Búsqueda de patrones en secuencias
    """
    
    # Codon patterns
    START_CODON = "ATG"
    STOP_CODONS = ["TAA", "TAG", "TGA"]
    
    def __init__(self):
        self._last_result: CodonAnalysisResult = None
        self._seq_object: Seq = None
        self._atg_positions: List[int] = []  # Nueva: almacenar posiciones
    
    def _analyze_spatial_distribution(self, positions: List[int], genome_length: int, num_windows: int = 100) -> Dict:
        """
        Analiza la distribución espacial de codones a lo largo del genoma.
        Calcula desviación estándar, varianza y coeficiente de variación.
        
        Args:
            positions: Lista de posiciones de codones
            genome_length: Tamaño total del genoma
            num_windows: Número de ventanas para dividir el genoma
        
        Returns:
            Dict con estadísticas de distribución espacial
        """
        if not positions or genome_length == 0:
            return {
                'mean_per_window': 0,
                'median_per_window': 0,
                'std_deviation': 0,
                'variance': 0,
                'coefficient_of_variation': 0,
                'uniformity_score': 0
            }
        
        # Dividir genoma en ventanas
        window_size = max(1, genome_length // num_windows)
        count_per_window = [0] * num_windows
        
        for pos in positions:
            window_idx = min(pos // window_size, num_windows - 1)
            count_per_window[window_idx] += 1
        
        # Calcular estadísticas
        counts_array = np.array(count_per_window)
        mean = float(np.mean(counts_array))
        std_dev = float(np.std(counts_array))
        
        # Coeficiente de variación (menor = más uniforme)
        cv = (std_dev / mean * 100) if mean > 0 else 0
        
        # Score de uniformidad (100 = perfectamente uniforme)
        uniformity = max(0, 100 - cv)
        
        # Calcular distancias entre codones consecutivos
        if len(positions) > 1:
            sorted_pos = sorted(positions)
            distances = np.diff(sorted_pos)
            mean_distance = float(np.mean(distances))
            std_distance = float(np.std(distances))
        else:
            mean_distance = 0
            std_distance = 0
        
        return {
            'mean_per_window': round(mean, 2),
            'median_per_window': round(float(np.median(counts_array)), 2),
            'std_deviation': round(std_dev, 3),
            'variance': round(float(np.var(counts_array)), 3),
            'coefficient_of_variation': round(cv, 2),
            'uniformity_score': round(uniformity, 1),
            'min_per_window': int(np.min(counts_array)),
            'max_per_window': int(np.max(counts_array)),
            'mean_distance_between_codons': round(mean_distance, 1),
            'std_distance_between_codons': round(std_distance, 1),
            'num_windows': num_windows,
            'window_size_bp': window_size
        }
    
    def _calculate_statistical_quality(self, atg_count: int, total_stops: int, genome_length: int, annotated_genes: int) -> Dict:
        """
        Calcula métricas de calidad estadística del análisis.
        Compara con valores esperados para minimizar desviación.
        
        Args:
            atg_count: Cantidad de ATG encontrados
            total_stops: Total de stop codons
            genome_length: Tamaño del genoma
            annotated_genes: Genes anotados oficialmente
        
        Returns:
            Dict con métricas de calidad
        """
        # Ratio ATG/Stop (ideal ~1.0 para genomas bacterianos)
        atg_stop_ratio = atg_count / total_stops if total_stops > 0 else 0
        
        # Desviación de genes esperados vs encontrados
        if annotated_genes > 0:
            gene_deviation_pct = abs(atg_count - annotated_genes) / annotated_genes * 100
            accuracy = max(0, 100 - gene_deviation_pct)
        else:
            gene_deviation_pct = 0
            accuracy = 0
        
        # Densidad esperada para E. coli: ~0.8-1.2 ATG per kb
        expected_density_min = 0.8
        expected_density_max = 1.2
        actual_density = (atg_count / genome_length) * 1000
        
        density_within_expected = expected_density_min <= actual_density <= expected_density_max
        
        return {
            'atg_stop_ratio': round(atg_stop_ratio, 3),
            'gene_count_accuracy': round(accuracy, 1),
            'gene_deviation_percent': round(gene_deviation_pct, 2),
            'density_within_expected_range': density_within_expected,
            'expected_density_range': f"{expected_density_min}-{expected_density_max}",
            'actual_density': round(actual_density, 3),
            'analysis_quality_score': round((accuracy + (100 if density_within_expected else 50)) / 2, 1)
        }
    
    def analyze(self, sequence: Union[str, Seq], annotated_genes: int = 0) -> CodonAnalysisResult:
        """
        Perform comprehensive codon analysis on a sequence using BioPython
        
        Implementación del algoritmo 7.1 del informe con mejoras estadísticas:
        - Cargar genoma usando Bio.Seq
        - Buscar ATG usando seq.count("ATG")
        - Buscar codones stop usando seq.count()
        - Calcular densidad: (count / genome_length) * 1000
        - Calcular distribución espacial con desviación estándar
        - Calcular métricas de calidad estadística
        
        Args:
            sequence: The DNA sequence (str or Bio.Seq object)
            annotated_genes: Number of annotated genes for comparison
            
        Returns:
            CodonAnalysisResult with all metrics including spatial distribution and quality
        """
        # Convert to Bio.Seq object if string
        if isinstance(sequence, str):
            self._seq_object = Seq(sequence.upper())
        else:
            self._seq_object = sequence
        
        genome_length = len(self._seq_object)
        
        # Count ATG (start codons) using BioPython's Seq.count()
        atg_count = self._count_codon_biopython(self.START_CODON)
        
        # Encontrar posiciones de ATG para análisis espacial (como en los scripts)
        self._atg_positions = self.find_codon_positions(self.START_CODON)
        
        # Calculate ATG density (per kilobase)
        atg_density = round((atg_count / genome_length) * 1000, 3) if genome_length > 0 else 0
        
        # Count stop codons using BioPython
        stop_counts = {}
        stop_positions_all = {}  # Almacenar posiciones de cada stop codon
        total_stops = 0
        
        for codon in self.STOP_CODONS:
            count = self._count_codon_biopython(codon)
            positions = self.find_codon_positions(codon)
            stop_counts[codon] = count
            stop_positions_all[codon] = positions
            total_stops += count
        
        # Calculate percentages for stop codons (como en analyze_stop_codons.py)
        stop_codons_data = {}
        for codon, count in stop_counts.items():
            percentage = round((count / total_stops) * 100, 2) if total_stops > 0 else 0
            density = round((count / genome_length) * 1000, 3) if genome_length > 0 else 0
            
            # Calcular distribución espacial para este stop codon
            spatial_dist = self._analyze_spatial_distribution(
                stop_positions_all[codon], 
                genome_length,
                num_windows=100
            )
            
            stop_codons_data[codon] = {
                "count": count,
                "percentage": percentage,
                "density_per_kb": density,
                "spatial_distribution": spatial_dist,
                "first_10_positions": stop_positions_all[codon][:10],
                "last_10_positions": stop_positions_all[codon][-10:] if len(stop_positions_all[codon]) >= 10 else []
            }
        
        # Gene comparison (más detallado)
        ratio_atg_genes = round(atg_count / annotated_genes, 2) if annotated_genes > 0 else 0
        non_coding_atg = atg_count - annotated_genes if annotated_genes > 0 else 0
        
        gene_comparison = {
            "annotated_genes": annotated_genes,
            "atg_found": atg_count,
            "difference": atg_count - annotated_genes,
            "ratio_atg_to_genes": ratio_atg_genes,
            "estimated_non_coding_atg": non_coding_atg,
            "percent_non_coding": round((non_coding_atg / atg_count * 100), 2) if atg_count > 0 else 0
        }
        
        # Análisis de distribución espacial de ATG (como en analyze_atg.py)
        spatial_distribution = self._analyze_spatial_distribution(
            self._atg_positions,
            genome_length,
            num_windows=100
        )
        
        # Métricas de calidad estadística
        statistical_quality = self._calculate_statistical_quality(
            atg_count,
            total_stops,
            genome_length,
            annotated_genes
        )
        
        result = CodonAnalysisResult(
            genome_length=genome_length,
            atg_count=atg_count,
            atg_density=atg_density,
            stop_codons=stop_codons_data,
            gene_comparison=gene_comparison,
            spatial_distribution=spatial_distribution,
            statistical_quality=statistical_quality
        )
        
        self._last_result = result
        return result
    
    def _count_codon_biopython(self, codon: str) -> int:
        """
        Count occurrences of a codon using BioPython's Seq.count()
        
        Bio.Seq.count() es más eficiente que regex para secuencias grandes
        y maneja correctamente secuencias de DNA.
        
        Args:
            codon: The codon pattern to count (e.g., "ATG")
            
        Returns:
            Number of non-overlapping occurrences
        """
        if self._seq_object is None:
            return 0
        return self._seq_object.count(codon)
    
    def _count_pattern(self, sequence: str, pattern: str) -> int:
        """
        Fallback: Count occurrences using string count
        Kept for compatibility
        """
        return sequence.upper().count(pattern)
    
    def find_codon_positions(self, codon: str) -> List[int]:
        """
        Find all positions of a codon in the sequence using BioPython
        
        Utiliza Bio.Seq.find() iterativamente para encontrar
        todas las posiciones de un codón.
        
        Args:
            codon: The codon to find (e.g., "ATG")
            
        Returns:
            List of positions (0-indexed)
        """
        if self._seq_object is None:
            return []
        
        positions = []
        start = 0
        while True:
            pos = self._seq_object.find(codon, start)
            if pos == -1:
                break
            positions.append(pos)
            start = pos + 1
        return positions
    
    def analyze_reading_frames(self, sequence: Union[str, Seq]) -> Dict[int, Dict]:
        """
        Analyze codons in each of the three reading frames
        
        Args:
            sequence: The DNA sequence to analyze
            
        Returns:
            Dictionary with analysis for each reading frame (0, 1, 2)
        """
        sequence = sequence.upper()
        results = {}
        
        for frame in range(3):
            frame_seq = sequence[frame:]
            # Trim to multiple of 3
            trim_length = len(frame_seq) - (len(frame_seq) % 3)
            frame_seq = frame_seq[:trim_length]
            
            atg_count = 0
            stop_counts = {codon: 0 for codon in self.STOP_CODONS}
            
            # Iterate through codons
            for i in range(0, len(frame_seq), 3):
                codon = frame_seq[i:i+3]
                if codon == self.START_CODON:
                    atg_count += 1
                elif codon in self.STOP_CODONS:
                    stop_counts[codon] += 1
            
            results[frame] = {
                "atg_count": atg_count,
                "stop_codons": stop_counts,
                "total_codons": len(frame_seq) // 3
            }
        
        return results
    
    def find_orfs(self, sequence: str, min_length: int = 100) -> list:
        """
        Find Open Reading Frames (ORFs) in the sequence
        
        Args:
            sequence: DNA sequence to analyze
            min_length: Minimum ORF length in nucleotides
            
        Returns:
            List of ORFs with their positions
        """
        sequence = sequence.upper()
        orfs = []
        
        # Search in all three reading frames
        for frame in range(3):
            frame_seq = sequence[frame:]
            i = 0
            
            while i < len(frame_seq) - 2:
                codon = frame_seq[i:i+3]
                
                if codon == self.START_CODON:
                    # Look for stop codon
                    for j in range(i + 3, len(frame_seq) - 2, 3):
                        stop_codon = frame_seq[j:j+3]
                        if stop_codon in self.STOP_CODONS:
                            orf_length = j - i + 3
                            if orf_length >= min_length:
                                orfs.append({
                                    "start": frame + i,
                                    "end": frame + j + 3,
                                    "length": orf_length,
                                    "frame": frame,
                                    "stop_codon": stop_codon
                                })
                            break
                i += 3
        
        return sorted(orfs, key=lambda x: x["length"], reverse=True)
    
    def get_codon_usage(self, sequence: str) -> Dict[str, int]:
        """
        Calculate usage of all 64 possible codons
        
        Args:
            sequence: DNA sequence to analyze
            
        Returns:
            Dictionary mapping each codon to its count
        """
        sequence = sequence.upper()
        codon_counts = {}
        
        # Initialize all possible codons
        bases = ['A', 'T', 'G', 'C']
        for b1 in bases:
            for b2 in bases:
                for b3 in bases:
                    codon_counts[b1 + b2 + b3] = 0
        
        # Count codons (reading frame 0)
        for i in range(0, len(sequence) - 2, 3):
            codon = sequence[i:i+3]
            if codon in codon_counts:
                codon_counts[codon] += 1
        
        return codon_counts
    
    def to_dict(self) -> Dict:
        """
        Convert last result to dictionary with full statistical data.
        Includes all metrics from analyze_atg.py and analyze_stop_codons.py
        """
        if not self._last_result:
            return {}
        
        result = {
            "genome_length": self._last_result.genome_length,
            "atg_count": self._last_result.atg_count,
            "atg_density": self._last_result.atg_density,
            "stop_codons": self._last_result.stop_codons,
            "gene_comparison": self._last_result.gene_comparison
        }
        
        # Agregar distribución espacial si está disponible
        if self._last_result.spatial_distribution:
            result["spatial_distribution"] = self._last_result.spatial_distribution
        
        # Agregar calidad estadística si está disponible
        if self._last_result.statistical_quality:
            result["statistical_quality"] = self._last_result.statistical_quality
        
        return result
