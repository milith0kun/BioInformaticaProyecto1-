"""
Genome Parser Module
Parses GenBank and other genomic file formats using BioPython

Módulos BioPython utilizados:
- Bio.SeqIO: Parsing de archivos GenBank y FASTA
- Bio.Seq: Representación y manipulación de secuencias
- Bio.SeqRecord: Gestión completa de información genómica
- Bio.SeqFeature: Extracción de features (genes, CDS)
- Bio.SeqUtils: Cálculo de contenido GC y otras utilidades
"""
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature
from Bio.SeqUtils import gc_fraction  # GC content calculation
from typing import List, Dict, Optional, Tuple, Generator
import os
import json
from dataclasses import dataclass, asdict


@dataclass
class GeneData:
    locus_tag: str
    start: int
    end: int
    length: int
    strand: int
    product: Optional[str]
    gc_content: float
    gene_name: Optional[str] = None
    protein_id: Optional[str] = None


@dataclass 
class GenomeData:
    sequence: str
    length: int
    gc_content: float
    description: str
    organism: str
    accession: str
    genes: List[GeneData]
    cds_count: int
    

class GenomeParser:
    """Parser for GenBank and FASTA genomic files"""
    
    def __init__(self, cache_dir: Optional[str] = None):
        self.cache_dir = cache_dir
        self._cached_genome: Optional[GenomeData] = None
        self._cached_filepath: Optional[str] = None
        
    def parse_genbank(self, filepath: str) -> GenomeData:
        """Parse a GenBank (.gbff) file"""
        # Check cache
        if self._cached_filepath == filepath and self._cached_genome:
            return self._cached_genome
            
        print(f"\n🧬 [PARSER] Parseando archivo: {filepath}")
        genes = []
        cds_count = 0
        
        # Parse GenBank file - may contain multiple records
        records = list(SeqIO.parse(filepath, "genbank"))
        
        if not records:
            raise ValueError(f"No records found in {filepath}")
        
        print(f"📊 [PARSER] Total de records en archivo: {len(records)}")
        
        # Concatenate all records (chromosomes + plasmids or contigs)
        full_sequence = ""
        all_genes = []
        all_cds_count = 0
        
        organism = records[0].annotations.get('organism', 'Unknown')
        description = records[0].description
        accession = records[0].id
        
        # Process each record (chromosome, plasmid, or contig)
        for idx, record in enumerate(records):
            sequence_part = str(record.seq)
            full_sequence += sequence_part
            
            if idx == 0:
                print(f"📏 [PARSER] Record 1: {len(sequence_part):,} bp - {record.description[:80]}")
            elif idx < 5:
                print(f"📏 [PARSER] Record {idx+1}: {len(sequence_part):,} bp - {record.description[:80]}")
            elif idx == 5:
                print(f"📏 [PARSER] ... y {len(records)-5} records más")
            
            # Extract genes and CDS from this record
            for feature in record.features:
                if feature.type == "gene":
                    gene_data = self._extract_gene_data(feature, record.seq)
                    if gene_data:
                        all_genes.append(gene_data)
                        
                elif feature.type == "CDS":
                    all_cds_count += 1
        
        genes = all_genes
        cds_count = all_cds_count
        
        print(f"📏 [PARSER] Longitud total del genoma: {len(full_sequence):,} bp ({len(records)} records)")
        print(f"🔬 [PARSER] Organismo: {organism}")
        print(f"🆔 [PARSER] Accession: {accession}")
        
        # Calculate genome GC content using BioPython on the full concatenated sequence
        from Bio.Seq import Seq
        full_seq_obj = Seq(full_sequence)
        gc_content = self._calculate_gc_biopython(full_seq_obj)
        
        print(f"🧮 [PARSER] Genes encontrados: {len(genes)}")
        print(f"🧮 [PARSER] CDS encontrados: {cds_count}")
        print(f"🧮 [PARSER] GC content: {gc_content:.2f}%\n")
        
        genome_data = GenomeData(
            sequence=full_sequence,
            length=len(full_sequence),
            gc_content=gc_content,
            description=description,
            organism=organism,
            accession=accession,
            genes=genes,
            cds_count=cds_count
        )
        
        # Cache result
        self._cached_genome = genome_data
        self._cached_filepath = filepath
        
        return genome_data
    
    def _extract_gene_data(self, feature: SeqFeature, sequence) -> Optional[GeneData]:
        """
        Extract gene data from a SeqFeature using BioPython
        
        Utiliza Bio.SeqFeature para:
        - Obtener ubicación genómica (start, end, strand)
        - Extraer qualifiers (locus_tag, product, protein_id)
        - Extraer la secuencia del gen usando feature.extract()
        """
        try:
            # Get location using BioPython's SeqFeature.location
            start = int(feature.location.start)
            end = int(feature.location.end)
            strand = feature.location.strand if feature.location.strand else 1
            
            # Get qualifiers from GenBank annotations
            qualifiers = feature.qualifiers
            locus_tag = qualifiers.get('locus_tag', ['unknown'])[0]
            gene_name = qualifiers.get('gene', [None])[0]
            product = qualifiers.get('product', [None])[0]
            protein_id = qualifiers.get('protein_id', [None])[0]
            
            # Extract gene sequence using BioPython's feature.extract()
            # This properly handles strand and compound locations
            gene_seq = feature.extract(sequence)
            
            # Calculate GC content using Bio.SeqUtils.gc_fraction
            gc_content = self._calculate_gc_biopython(gene_seq)
            
            return GeneData(
                locus_tag=locus_tag,
                start=start,
                end=end,
                length=len(gene_seq),
                strand=strand,
                product=product,
                gc_content=gc_content,
                gene_name=gene_name,
                protein_id=protein_id
            )
        except Exception as e:
            print(f"Error extracting gene data: {e}")
            return None
    
    def _calculate_gc_biopython(self, sequence) -> float:
        """
        Calculate GC content using BioPython's Bio.SeqUtils.gc_fraction
        
        Esta función usa la implementación optimizada de BioPython
        para calcular el porcentaje de G+C en una secuencia.
        """
        if not sequence or len(sequence) == 0:
            return 0.0
        
        try:
            # gc_fraction returns a value between 0 and 1
            gc = gc_fraction(sequence) * 100
            return round(gc, 2)
        except Exception:
            # Fallback to manual calculation
            return self._calculate_gc(str(sequence))
    
    def _calculate_gc(self, sequence: str) -> float:
        """Calculate GC content percentage"""
        if not sequence:
            return 0.0
        
        sequence = sequence.upper()
        g_count = sequence.count('G')
        c_count = sequence.count('C')
        total = len(sequence)
        
        if total == 0:
            return 0.0
            
        return round((g_count + c_count) / total * 100, 2)
    
    def parse_fasta(self, filepath: str) -> Tuple[str, int]:
        """Parse a FASTA file and return sequence and length"""
        records = list(SeqIO.parse(filepath, "fasta"))
        
        if not records:
            raise ValueError(f"No records found in {filepath}")
        
        # Concatenate all sequences (usually just one for chromosome)
        full_sequence = "".join(str(record.seq) for record in records)
        
        return full_sequence, len(full_sequence)
    
    def get_sequence_iterator(self, filepath: str, format: str = "genbank") -> Generator:
        """Get an iterator over sequences for memory-efficient processing"""
        return SeqIO.parse(filepath, format)
    
    def save_to_cache(self, data: Dict, cache_key: str) -> None:
        """Save parsed data to cache"""
        if not self.cache_dir:
            return
            
        os.makedirs(self.cache_dir, exist_ok=True)
        cache_path = os.path.join(self.cache_dir, f"{cache_key}.json")
        
        with open(cache_path, 'w') as f:
            json.dump(data, f)
    
    def load_from_cache(self, cache_key: str) -> Optional[Dict]:
        """Load data from cache if exists"""
        if not self.cache_dir:
            return None
            
        cache_path = os.path.join(self.cache_dir, f"{cache_key}.json")
        
        if os.path.exists(cache_path):
            with open(cache_path, 'r') as f:
                return json.load(f)
        
        return None
