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


@dataclass
class CodonAnalysisResult:
    genome_length: int
    atg_count: int
    atg_density: float
    stop_codons: Dict[str, Dict]
    gene_comparison: Dict[str, int]


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
    
    def analyze(self, sequence: Union[str, Seq], annotated_genes: int = 0) -> CodonAnalysisResult:
        """
        Perform comprehensive codon analysis on a sequence using BioPython
        
        Implementación del algoritmo 7.1 del informe:
        - Cargar genoma usando Bio.Seq
        - Buscar ATG usando seq.count("ATG")
        - Buscar codones stop usando seq.count()
        - Calcular densidad: (count / genome_length) * 1000
        
        Args:
            sequence: The DNA sequence (str or Bio.Seq object)
            annotated_genes: Number of annotated genes for comparison
            
        Returns:
            CodonAnalysisResult with all metrics
        """
        # Convert to Bio.Seq object if string
        if isinstance(sequence, str):
            self._seq_object = Seq(sequence.upper())
        else:
            self._seq_object = sequence
        
        genome_length = len(self._seq_object)
        
        # Count ATG (start codons) using BioPython's Seq.count()
        atg_count = self._count_codon_biopython(self.START_CODON)
        
        # Calculate ATG density (per kilobase)
        atg_density = round((atg_count / genome_length) * 1000, 3) if genome_length > 0 else 0
        
        # Count stop codons using BioPython
        stop_counts = {}
        total_stops = 0
        
        for codon in self.STOP_CODONS:
            count = self._count_codon_biopython(codon)
            stop_counts[codon] = count
            total_stops += count
        
        # Calculate percentages for stop codons
        stop_codons_data = {}
        for codon, count in stop_counts.items():
            percentage = round((count / total_stops) * 100, 1) if total_stops > 0 else 0
            stop_codons_data[codon] = {
                "count": count,
                "percentage": percentage
            }
        
        # Gene comparison
        gene_comparison = {
            "annotated_genes": annotated_genes,
            "atg_found": atg_count,
            "difference": atg_count - annotated_genes
        }
        
        result = CodonAnalysisResult(
            genome_length=genome_length,
            atg_count=atg_count,
            atg_density=atg_density,
            stop_codons=stop_codons_data,
            gene_comparison=gene_comparison
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
        """Convert last result to dictionary"""
        if not self._last_result:
            return {}
        
        return {
            "genome_length": self._last_result.genome_length,
            "atg_count": self._last_result.atg_count,
            "atg_density": self._last_result.atg_density,
            "stop_codons": self._last_result.stop_codons,
            "gene_comparison": self._last_result.gene_comparison
        }
