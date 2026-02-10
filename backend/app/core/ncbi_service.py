"""
NCBI Enhanced Service Module
Provides advanced NCBI API integration for:
- Protein retrieval from genomes
- Gene location/position lookup
- Full sequence retrieval
- Literature search via PubMed/E-Utilities
- Enhanced codon usage from CDS sequences

Uses NCBI Datasets API v2 + E-Utilities (Entrez)
"""
import os
import re
import json
import requests
import time
from typing import Dict, Optional, List, Any
from dataclasses import dataclass, field
from Bio import SeqIO, Entrez
from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction


# Configure Entrez
Entrez.email = "bioinformatics@proyecto1.edu"
Entrez.tool = "BioInformatica-Proyecto1"


@dataclass
class ProteinInfo:
    protein_id: str
    locus_tag: str
    gene_name: str
    product: str
    length: int
    sequence: str
    molecular_weight_approx: float
    start: int
    end: int
    strand: int


@dataclass
class GeneDetail:
    locus_tag: str
    gene_name: str
    start: int
    end: int
    strand: int
    length: int
    product: str
    protein_id: str
    gc_content: float
    sequence: str
    translation: str
    codon_start: int
    note: str
    db_xref: List[str]


@dataclass
class CodonUsageComplete:
    """Complete codon usage analysis with RSCU and CAI"""
    codon_counts: Dict[str, int]
    codon_frequencies: Dict[str, float]
    rscu: Dict[str, float]  # Relative Synonymous Codon Usage
    amino_acid_usage: Dict[str, Dict[str, Any]]
    total_codons: int
    cai_reference: Dict[str, float]  # CAI reference values
    gc3_content: float  # GC at 3rd position
    effective_number_of_codons: float  # Nc


# Genetic code table (standard)
GENETIC_CODE = {
    'TTT': 'Phe', 'TTC': 'Phe',
    'TTA': 'Leu', 'TTG': 'Leu', 'CTT': 'Leu', 'CTC': 'Leu', 'CTA': 'Leu', 'CTG': 'Leu',
    'ATT': 'Ile', 'ATC': 'Ile', 'ATA': 'Ile',
    'ATG': 'Met',
    'GTT': 'Val', 'GTC': 'Val', 'GTA': 'Val', 'GTG': 'Val',
    'TCT': 'Ser', 'TCC': 'Ser', 'TCA': 'Ser', 'TCG': 'Ser', 'AGT': 'Ser', 'AGC': 'Ser',
    'CCT': 'Pro', 'CCC': 'Pro', 'CCA': 'Pro', 'CCG': 'Pro',
    'ACT': 'Thr', 'ACC': 'Thr', 'ACA': 'Thr', 'ACG': 'Thr',
    'GCT': 'Ala', 'GCC': 'Ala', 'GCA': 'Ala', 'GCG': 'Ala',
    'TAT': 'Tyr', 'TAC': 'Tyr',
    'TAA': 'Stop', 'TAG': 'Stop', 'TGA': 'Stop',
    'CAT': 'His', 'CAC': 'His',
    'CAA': 'Gln', 'CAG': 'Gln',
    'AAT': 'Asn', 'AAC': 'Asn',
    'AAA': 'Lys', 'AAG': 'Lys',
    'GAT': 'Asp', 'GAC': 'Asp',
    'GAA': 'Glu', 'GAG': 'Glu',
    'TGT': 'Cys', 'TGC': 'Cys',
    'TGG': 'Trp',
    'CGT': 'Arg', 'CGC': 'Arg', 'CGA': 'Arg', 'CGG': 'Arg', 'AGA': 'Arg', 'AGG': 'Arg',
    'GGT': 'Gly', 'GGC': 'Gly', 'GGA': 'Gly', 'GGG': 'Gly',
}

AMINO_ACID_NAMES = {
    'Ala': 'Alanina', 'Arg': 'Arginina', 'Asn': 'Asparagina', 'Asp': 'Ácido aspártico',
    'Cys': 'Cisteína', 'Gln': 'Glutamina', 'Glu': 'Ácido glutámico', 'Gly': 'Glicina',
    'His': 'Histidina', 'Ile': 'Isoleucina', 'Leu': 'Leucina', 'Lys': 'Lisina',
    'Met': 'Metionina', 'Phe': 'Fenilalanina', 'Pro': 'Prolina', 'Ser': 'Serina',
    'Thr': 'Treonina', 'Trp': 'Triptófano', 'Tyr': 'Tirosina', 'Val': 'Valina',
    'Stop': 'Stop'
}

AMINO_ACID_1LETTER = {
    'Ala': 'A', 'Arg': 'R', 'Asn': 'N', 'Asp': 'D', 'Cys': 'C',
    'Gln': 'Q', 'Glu': 'E', 'Gly': 'G', 'His': 'H', 'Ile': 'I',
    'Leu': 'L', 'Lys': 'K', 'Met': 'M', 'Phe': 'F', 'Pro': 'P',
    'Ser': 'S', 'Thr': 'T', 'Trp': 'W', 'Tyr': 'Y', 'Val': 'V',
    'Stop': '*'
}


class NCBIService:
    """Enhanced NCBI service for bioinformatics analysis"""

    DATASETS_BASE = "https://api.ncbi.nlm.nih.gov/datasets/v2"
    EUTILS_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

    def __init__(self, api_key: Optional[str] = None):
        self.api_key = api_key or os.getenv("NCBI_API_KEY")
        self.session = requests.Session()
        self.session.headers.update({
            "Accept": "application/json",
            "User-Agent": "BioInformatica-Proyecto1/2.0"
        })
        if self.api_key:
            self.session.headers["api-key"] = self.api_key

    # ==================== PROTEIN METHODS ====================

    def get_proteins_from_genbank(self, genbank_path: str, limit: int = 500) -> List[Dict]:
        """
        Extract proteins from a parsed GenBank file
        Returns list of protein info dicts
        """
        proteins = []
        try:
            records = SeqIO.parse(genbank_path, "genbank")
            for record in records:
                for feature in record.features:
                    if feature.type == "CDS":
                        qualifiers = feature.qualifiers
                        translation = qualifiers.get('translation', [''])[0]
                        if not translation:
                            continue

                        protein_id = qualifiers.get('protein_id', [''])[0]
                        locus_tag = qualifiers.get('locus_tag', [''])[0]
                        gene_name = qualifiers.get('gene', [''])[0]
                        product = qualifiers.get('product', [''])[0]

                        start = int(feature.location.start)
                        end = int(feature.location.end)
                        strand = feature.location.strand or 1

                        # Approximate molecular weight (average amino acid MW ~110 Da)
                        mw_approx = len(translation) * 110.0

                        proteins.append({
                            "protein_id": protein_id,
                            "locus_tag": locus_tag,
                            "gene_name": gene_name,
                            "product": product,
                            "length": len(translation),
                            "sequence": translation[:100] + ("..." if len(translation) > 100 else ""),
                            "full_sequence": translation,
                            "molecular_weight_approx": round(mw_approx, 1),
                            "start": start,
                            "end": end,
                            "strand": strand
                        })

                        if len(proteins) >= limit:
                            break
                if len(proteins) >= limit:
                    break

        except Exception as e:
            print(f"Error extracting proteins: {e}")

        return proteins

    # ==================== GENE DETAIL METHODS ====================

    def get_gene_at_position(self, genbank_path: str, position: int) -> Optional[Dict]:
        """
        Find which gene is at a specific genomic position
        """
        try:
            records = list(SeqIO.parse(genbank_path, "genbank"))
            cumulative_offset = 0

            for record in records:
                record_length = len(record.seq)

                if position < cumulative_offset + record_length:
                    # Position is in this record
                    local_pos = position - cumulative_offset

                    for feature in record.features:
                        if feature.type in ("gene", "CDS"):
                            start = int(feature.location.start)
                            end = int(feature.location.end)

                            if start <= local_pos <= end:
                                qualifiers = feature.qualifiers
                                gene_seq = feature.extract(record.seq)

                                result = {
                                    "found": True,
                                    "position": position,
                                    "feature_type": feature.type,
                                    "locus_tag": qualifiers.get('locus_tag', [''])[0],
                                    "gene_name": qualifiers.get('gene', [''])[0],
                                    "product": qualifiers.get('product', [''])[0],
                                    "start": start + cumulative_offset,
                                    "end": end + cumulative_offset,
                                    "strand": feature.location.strand or 1,
                                    "length": len(gene_seq),
                                    "gc_content": round(gc_fraction(gene_seq) * 100, 2),
                                    "relative_position": round((local_pos - start) / (end - start) * 100, 1) if end > start else 0,
                                    "record_id": record.id,
                                    "record_description": record.description[:100]
                                }

                                # Add protein info if CDS
                                if feature.type == "CDS":
                                    result["protein_id"] = qualifiers.get('protein_id', [''])[0]
                                    result["translation"] = qualifiers.get('translation', [''])[0][:50] + "..."

                                return result

                    # Position found but not in any gene
                    return {
                        "found": False,
                        "position": position,
                        "message": "Posición en región intergénica",
                        "record_id": record.id,
                        "nearest_genes": self._find_nearest_genes(record, local_pos, cumulative_offset)
                    }

                cumulative_offset += record_length

            return {"found": False, "position": position, "message": "Posición fuera del genoma"}

        except Exception as e:
            return {"found": False, "position": position, "error": str(e)}

    def _find_nearest_genes(self, record, position: int, offset: int = 0) -> List[Dict]:
        """Find the nearest genes to a given position"""
        genes = []
        for feature in record.features:
            if feature.type == "gene":
                start = int(feature.location.start)
                end = int(feature.location.end)
                distance = min(abs(position - start), abs(position - end))
                qualifiers = feature.qualifiers
                genes.append({
                    "locus_tag": qualifiers.get('locus_tag', [''])[0],
                    "gene_name": qualifiers.get('gene', [''])[0],
                    "start": start + offset,
                    "end": end + offset,
                    "distance": distance,
                    "direction": "upstream" if start > position else "downstream"
                })

        genes.sort(key=lambda x: x["distance"])
        return genes[:4]  # Return 4 nearest genes

    def get_gene_detail(self, genbank_path: str, locus_tag: str) -> Optional[Dict]:
        """
        Get detailed information about a specific gene by locus_tag
        """
        try:
            records = list(SeqIO.parse(genbank_path, "genbank"))
            cumulative_offset = 0

            for record in records:
                for feature in record.features:
                    if feature.type in ("gene", "CDS"):
                        qualifiers = feature.qualifiers
                        lt = qualifiers.get('locus_tag', [''])[0]

                        if lt == locus_tag:
                            gene_seq = feature.extract(record.seq)
                            start = int(feature.location.start)
                            end = int(feature.location.end)

                            result = {
                                "locus_tag": lt,
                                "gene_name": qualifiers.get('gene', [''])[0],
                                "start": start + cumulative_offset,
                                "end": end + cumulative_offset,
                                "strand": feature.location.strand or 1,
                                "length": len(gene_seq),
                                "product": qualifiers.get('product', ['unknown'])[0],
                                "gc_content": round(gc_fraction(gene_seq) * 100, 2),
                                "sequence": str(gene_seq),
                                "codon_start": int(qualifiers.get('codon_start', [1])[0]),
                                "note": qualifiers.get('note', [''])[0] if 'note' in qualifiers else '',
                                "db_xref": qualifiers.get('db_xref', []),
                                "record_id": record.id
                            }

                            # Add translation if CDS
                            if feature.type == "CDS":
                                result["protein_id"] = qualifiers.get('protein_id', [''])[0]
                                result["translation"] = qualifiers.get('translation', [''])[0]
                                result["feature_type"] = "CDS"
                            else:
                                result["feature_type"] = "gene"

                            return result

                cumulative_offset += len(record.seq)

            return None

        except Exception as e:
            print(f"Error getting gene detail: {e}")
            return None

    # ==================== SEQUENCE METHODS ====================

    def get_genome_sequence_segment(self, genbank_path: str, start: int = 0, end: int = 1000) -> Dict:
        """
        Get a segment of the genome sequence with gene annotations
        """
        try:
            records = list(SeqIO.parse(genbank_path, "genbank"))
            full_sequence = "".join(str(r.seq) for r in records)
            genome_length = len(full_sequence)

            # Clamp values
            start = max(0, min(start, genome_length - 1))
            end = min(end, genome_length)

            segment = full_sequence[start:end]

            # Find genes in this segment
            genes_in_segment = []
            cumulative = 0
            for record in records:
                for feature in record.features:
                    if feature.type == "gene":
                        g_start = int(feature.location.start) + cumulative
                        g_end = int(feature.location.end) + cumulative

                        # Check overlap
                        if g_start < end and g_end > start:
                            qualifiers = feature.qualifiers
                            genes_in_segment.append({
                                "locus_tag": qualifiers.get('locus_tag', [''])[0],
                                "gene_name": qualifiers.get('gene', [''])[0],
                                "start": g_start,
                                "end": g_end,
                                "strand": feature.location.strand or 1,
                                "product": qualifiers.get('product', [''])[0],
                                "segment_start": max(0, g_start - start),
                                "segment_end": min(end - start, g_end - start)
                            })
                cumulative += len(record.seq)

            # Calculate local GC content
            gc = round(gc_fraction(Seq(segment)) * 100, 2) if segment else 0

            return {
                "sequence": segment,
                "start": start,
                "end": end,
                "length": len(segment),
                "genome_length": genome_length,
                "gc_content": gc,
                "genes": genes_in_segment,
                "formatted": self._format_sequence(segment, start)
            }

        except Exception as e:
            return {"error": str(e)}

    def _format_sequence(self, sequence: str, offset: int = 0, line_width: int = 60) -> List[Dict]:
        """Format sequence with position numbers like FASTA display"""
        lines = []
        for i in range(0, len(sequence), line_width):
            chunk = sequence[i:i + line_width]
            lines.append({
                "position": offset + i + 1,
                "sequence": chunk,
                "gc_count": chunk.count('G') + chunk.count('C'),
                "at_count": chunk.count('A') + chunk.count('T'),
            })
        return lines

    # ==================== CODON USAGE METHODS ====================

    def calculate_complete_codon_usage(self, genbank_path: str) -> Dict:
        """
        Calculate complete codon usage from all CDS in the genome.
        Includes RSCU, amino acid frequencies, GC3 content.
        """
        codon_counts = {codon: 0 for codon in GENETIC_CODE.keys()}
        gc3_count = 0
        total_third_positions = 0
        total_codons = 0

        try:
            records = SeqIO.parse(genbank_path, "genbank")
            for record in records:
                for feature in record.features:
                    if feature.type == "CDS":
                        try:
                            cds_seq = str(feature.extract(record.seq)).upper()
                            # Count codons in reading frame
                            for i in range(0, len(cds_seq) - 2, 3):
                                codon = cds_seq[i:i + 3]
                                if codon in codon_counts:
                                    codon_counts[codon] += 1
                                    total_codons += 1

                                    # GC at 3rd position
                                    if codon[2] in ('G', 'C'):
                                        gc3_count += 1
                                    total_third_positions += 1
                        except Exception:
                            continue

        except Exception as e:
            print(f"Error calculating codon usage: {e}")

        if total_codons == 0:
            return {"error": "No codons found"}

        # Calculate frequencies
        codon_frequencies = {
            codon: round(count / total_codons * 1000, 3)
            for codon, count in codon_counts.items()
        }

        # Calculate RSCU
        rscu = self._calculate_rscu(codon_counts)

        # Group by amino acid
        amino_acid_usage = self._group_by_amino_acid(codon_counts, rscu, total_codons)

        # GC3 content
        gc3 = round(gc3_count / total_third_positions * 100, 2) if total_third_positions > 0 else 0

        # Effective number of codons (Nc) - Wright's formula approximation
        nc = self._calculate_nc(codon_counts)

        return {
            "codon_counts": codon_counts,
            "codon_frequencies": codon_frequencies,
            "rscu": rscu,
            "amino_acid_usage": amino_acid_usage,
            "total_codons": total_codons,
            "gc3_content": gc3,
            "effective_number_of_codons": round(nc, 2),
            "codon_table": self._build_codon_table(codon_counts, codon_frequencies, rscu)
        }

    def _calculate_rscu(self, codon_counts: Dict[str, int]) -> Dict[str, float]:
        """
        Calculate Relative Synonymous Codon Usage (RSCU)
        RSCU = observed / expected (if all synonymous codons used equally)
        """
        # Group codons by amino acid
        aa_groups = {}
        for codon, aa in GENETIC_CODE.items():
            if aa not in aa_groups:
                aa_groups[aa] = []
            aa_groups[aa].append(codon)

        rscu = {}
        for aa, codons in aa_groups.items():
            total = sum(codon_counts.get(c, 0) for c in codons)
            n_synonymous = len(codons)

            for codon in codons:
                observed = codon_counts.get(codon, 0)
                expected = total / n_synonymous if n_synonymous > 0 else 0

                if expected > 0:
                    rscu[codon] = round(observed / expected, 3)
                else:
                    rscu[codon] = 0.0

        return rscu

    def _group_by_amino_acid(self, codon_counts: Dict, rscu: Dict, total_codons: int) -> Dict:
        """Group codon usage by amino acid"""
        aa_usage = {}

        for codon, aa in GENETIC_CODE.items():
            if aa == 'Stop':
                continue

            if aa not in aa_usage:
                aa_usage[aa] = {
                    "name": AMINO_ACID_NAMES.get(aa, aa),
                    "letter": AMINO_ACID_1LETTER.get(aa, '?'),
                    "codons": [],
                    "total_count": 0,
                    "fraction": 0
                }

            count = codon_counts.get(codon, 0)
            aa_usage[aa]["codons"].append({
                "codon": codon,
                "count": count,
                "rscu": rscu.get(codon, 0),
                "frequency_per_thousand": round(count / total_codons * 1000, 2) if total_codons > 0 else 0
            })
            aa_usage[aa]["total_count"] += count

        # Calculate fraction for each amino acid
        all_aa_total = sum(v["total_count"] for v in aa_usage.values())
        for aa in aa_usage:
            aa_usage[aa]["fraction"] = round(
                aa_usage[aa]["total_count"] / all_aa_total * 100, 2
            ) if all_aa_total > 0 else 0

        return aa_usage

    def _calculate_nc(self, codon_counts: Dict[str, int]) -> float:
        """
        Calculate effective number of codons (Nc)
        Simplified Wright (1990) method
        """
        aa_groups = {}
        for codon, aa in GENETIC_CODE.items():
            if aa == 'Stop':
                continue
            if aa not in aa_groups:
                aa_groups[aa] = []
            aa_groups[aa].append(codon)

        f_values = {1: [], 2: [], 3: [], 4: [], 6: []}

        for aa, codons in aa_groups.items():
            k = len(codons)
            total = sum(codon_counts.get(c, 0) for c in codons)

            if total <= 1 or k < 2:
                continue

            # Calculate F (homozygosity)
            f = sum((codon_counts.get(c, 0) / total) ** 2 for c in codons)

            if k in f_values:
                f_values[k].append(f)

        # Nc formula
        nc = 2.0  # Met + Trp (single codon amino acids)

        for k, f_list in f_values.items():
            if f_list:
                avg_f = sum(f_list) / len(f_list)
                if avg_f > 0:
                    nc += k / avg_f
                else:
                    nc += k

        return min(nc, 61)  # Cap at 61

    def _build_codon_table(self, counts: Dict, freqs: Dict, rscu: Dict) -> List[Dict]:
        """Build a structured codon table for UI display"""
        table = []
        for codon in sorted(GENETIC_CODE.keys()):
            aa = GENETIC_CODE[codon]
            table.append({
                "codon": codon,
                "amino_acid": aa,
                "amino_acid_name": AMINO_ACID_NAMES.get(aa, aa),
                "letter": AMINO_ACID_1LETTER.get(aa, '?'),
                "count": counts.get(codon, 0),
                "frequency": freqs.get(codon, 0),
                "rscu": rscu.get(codon, 0)
            })
        return table

    # ==================== GENOME MAP DATA ====================

    def get_genome_map_data(self, genbank_path: str) -> Dict:
        """
        Get all data needed for a circular genome map visualization
        """
        try:
            records = list(SeqIO.parse(genbank_path, "genbank"))
            full_sequence = "".join(str(r.seq) for r in records)
            genome_length = len(full_sequence)

            # Collect genes for the map
            genes = []
            cumulative = 0
            for record in records:
                for feature in record.features:
                    if feature.type == "gene":
                        qualifiers = feature.qualifiers
                        start = int(feature.location.start) + cumulative
                        end = int(feature.location.end) + cumulative
                        strand = feature.location.strand or 1

                        gene_seq = feature.extract(record.seq)

                        genes.append({
                            "locus_tag": qualifiers.get('locus_tag', [''])[0],
                            "gene_name": qualifiers.get('gene', [''])[0],
                            "product": qualifiers.get('product', [''])[0],
                            "start": start,
                            "end": end,
                            "strand": strand,
                            "length": end - start,
                            "gc_content": round(gc_fraction(gene_seq) * 100, 2),
                            "angle_start": round(start / genome_length * 360, 2),
                            "angle_end": round(end / genome_length * 360, 2)
                        })
                cumulative += len(record.seq)

            # GC content by window (for the outer ring)
            window_size = genome_length // 500 if genome_length > 500 else 1
            gc_windows = []
            for i in range(0, genome_length, window_size):
                window_seq = full_sequence[i:i + window_size]
                gc = gc_fraction(Seq(window_seq)) * 100 if window_seq else 0
                gc_windows.append({
                    "position": i,
                    "gc": round(gc, 2),
                    "angle": round(i / genome_length * 360, 2)
                })

            # GC skew
            gc_skew = []
            for i in range(0, genome_length, window_size):
                window_seq = full_sequence[i:i + window_size].upper()
                g = window_seq.count('G')
                c = window_seq.count('C')
                skew = (g - c) / (g + c) if (g + c) > 0 else 0
                gc_skew.append({
                    "position": i,
                    "skew": round(skew, 4),
                    "angle": round(i / genome_length * 360, 2)
                })

            return {
                "genome_length": genome_length,
                "organism": records[0].annotations.get('organism', 'Unknown') if records else 'Unknown',
                "accession": records[0].id if records else '',
                "description": records[0].description[:100] if records else '',
                "total_genes": len(genes),
                "genes": genes,
                "gc_windows": gc_windows,
                "gc_skew": gc_skew,
                "window_size": window_size,
                "average_gc": round(gc_fraction(Seq(full_sequence)) * 100, 2)
            }

        except Exception as e:
            return {"error": str(e)}

    # ==================== PUBMED / LITERATURE ====================

    def search_literature(self, query: str, max_results: int = 10) -> List[Dict]:
        """
        Search PubMed for relevant literature using NCBI E-Utilities
        """
        try:
            # Search PubMed
            handle = Entrez.esearch(db="pubmed", term=query, retmax=max_results, sort="relevance")
            results = Entrez.read(handle)
            handle.close()

            pmids = results.get("IdList", [])
            if not pmids:
                return []

            # Fetch details
            handle = Entrez.efetch(db="pubmed", id=",".join(pmids), rettype="abstract", retmode="xml")
            articles = Entrez.read(handle)
            handle.close()

            references = []
            for article in articles.get("PubmedArticle", []):
                medline = article.get("MedlineCitation", {})
                article_data = medline.get("Article", {})

                title = article_data.get("ArticleTitle", "")
                abstract = article_data.get("Abstract", {}).get("AbstractText", [""])
                if isinstance(abstract, list):
                    abstract = " ".join(str(a) for a in abstract)

                # Authors
                authors_list = article_data.get("AuthorList", [])
                authors = []
                for auth in authors_list[:3]:
                    last = auth.get("LastName", "")
                    fore = auth.get("ForeName", "")
                    if last:
                        authors.append(f"{last} {fore}".strip())
                if len(authors_list) > 3:
                    authors.append("et al.")

                # Journal
                journal = article_data.get("Journal", {})
                journal_title = journal.get("Title", "")
                pub_date = journal.get("JournalIssue", {}).get("PubDate", {})
                year = pub_date.get("Year", "")

                pmid = str(medline.get("PMID", ""))

                references.append({
                    "pmid": pmid,
                    "title": title,
                    "authors": ", ".join(authors),
                    "journal": journal_title,
                    "year": year,
                    "abstract": abstract[:500] + ("..." if len(abstract) > 500 else ""),
                    "url": f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/",
                    "doi": ""
                })

            return references

        except Exception as e:
            print(f"Error searching PubMed: {e}")
            return []

    def get_genbank_reference(self, accession: str) -> Dict:
        """Get GenBank record reference info"""
        try:
            url = f"{self.DATASETS_BASE}/genome/accession/{accession}/dataset_report"
            response = self.session.get(url, timeout=30)
            response.raise_for_status()
            data = response.json()

            if "reports" in data and data["reports"]:
                report = data["reports"][0]
                return {
                    "accession": report.get("accession", accession),
                    "organism": report.get("organism", {}).get("organism_name", ""),
                    "assembly_info": report.get("assembly_info", {}),
                    "annotation_info": report.get("annotation_info", {}),
                    "bioproject": report.get("assembly_info", {}).get("bioproject_accession", ""),
                    "url": f"https://www.ncbi.nlm.nih.gov/datasets/genome/{accession}/"
                }

            return {"accession": accession, "error": "Not found"}

        except Exception as e:
            return {"accession": accession, "error": str(e)}

    # ==================== NCBI GENE / NUCLEOTIDE SEARCH ====================

    def search_ncbi_gene(self, query: str, organism: str = "", max_results: int = 5) -> List[Dict]:
        """
        Search NCBI Gene database for gene information.
        Returns gene records with links to NCBI.
        """
        try:
            search_term = query
            if organism:
                search_term = f"{query} AND {organism}[Organism]"

            handle = Entrez.esearch(db="gene", term=search_term, retmax=max_results, sort="relevance")
            results = Entrez.read(handle)
            handle.close()

            gene_ids = results.get("IdList", [])
            if not gene_ids:
                return []

            # Fetch gene details
            handle = Entrez.efetch(db="gene", id=",".join(gene_ids), rettype="docsum", retmode="xml")
            records = Entrez.read(handle)
            handle.close()

            genes = []
            doc_sums = records.get("DocumentSummarySet", {}).get("DocumentSummary", [])
            if not doc_sums:
                doc_sums = records if isinstance(records, list) else []

            for doc in doc_sums:
                gene_id = str(doc.get("Id", doc.get("uid", "")))
                name = str(doc.get("Name", doc.get("name", "")))
                description = str(doc.get("Description", doc.get("description", "")))
                organism_name = str(doc.get("Organism", {}).get("ScientificName", "")) if isinstance(doc.get("Organism"), dict) else str(doc.get("Organism", ""))
                chromosome = str(doc.get("Chromosome", ""))
                map_location = str(doc.get("MapLocation", ""))

                genes.append({
                    "gene_id": gene_id,
                    "name": name,
                    "description": description,
                    "organism": organism_name,
                    "chromosome": chromosome,
                    "map_location": map_location,
                    "url": f"https://www.ncbi.nlm.nih.gov/gene/{gene_id}",
                    "genbank_url": f"https://www.ncbi.nlm.nih.gov/gene/?term={name}"
                })

            return genes

        except Exception as e:
            print(f"Error searching NCBI Gene: {e}")
            return []

    def search_ncbi_nucleotide(self, query: str, max_results: int = 3) -> List[Dict]:
        """
        Search NCBI Nucleotide database.
        Returns nucleotide records with GenBank links.
        """
        try:
            handle = Entrez.esearch(db="nucleotide", term=query, retmax=max_results, sort="relevance")
            results = Entrez.read(handle)
            handle.close()

            nuc_ids = results.get("IdList", [])
            if not nuc_ids:
                return []

            handle = Entrez.efetch(db="nucleotide", id=",".join(nuc_ids), rettype="docsum", retmode="xml")
            records = Entrez.read(handle)
            handle.close()

            nucleotides = []
            for record in records:
                acc = str(record.get("AccessionVersion", record.get("Caption", "")))
                title = str(record.get("Title", ""))
                length = int(record.get("Length", 0))

                nucleotides.append({
                    "accession": acc,
                    "title": title,
                    "length": length,
                    "url": f"https://www.ncbi.nlm.nih.gov/nuccore/{acc}",
                    "genbank_url": f"https://www.ncbi.nlm.nih.gov/nuccore/{acc}?report=genbank"
                })

            return nucleotides

        except Exception as e:
            print(f"Error searching NCBI Nucleotide: {e}")
            return []

    # ==================== CENTRAL DOGMA DATA ====================

    def get_central_dogma_data(self, genbank_path: str, locus_tag: str) -> Optional[Dict]:
        """
        Get complete central dogma data for a gene:
        DNA (5'→3' template + coding strand) → mRNA → Protein
        Shows transcription and translation process
        """
        try:
            records = list(SeqIO.parse(genbank_path, "genbank"))

            for record in records:
                for feature in record.features:
                    if feature.type == "CDS":
                        qualifiers = feature.qualifiers
                        lt = qualifiers.get('locus_tag', [''])[0]

                        if lt == locus_tag:
                            # Extract CDS sequence (coding strand, 5'→3')
                            coding_seq = str(feature.extract(record.seq)).upper()
                            strand = feature.location.strand or 1

                            # DNA: coding strand (sense) 5'→3'
                            dna_sense_5to3 = coding_seq

                            # DNA: template strand (antisense) 3'→5'
                            complement_map = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
                            dna_template_3to5 = ''.join(complement_map.get(b, 'N') for b in dna_sense_5to3)

                            # mRNA: same as coding strand but U instead of T (5'→3')
                            mrna_5to3 = dna_sense_5to3.replace('T', 'U')

                            # Protein: translate mRNA
                            translation = qualifiers.get('translation', [''])[0]
                            if not translation:
                                # Translate manually
                                protein = []
                                for i in range(0, len(coding_seq) - 2, 3):
                                    codon = coding_seq[i:i+3]
                                    aa = GENETIC_CODE.get(codon, '?')
                                    if aa == 'Stop':
                                        protein.append('*')
                                        break
                                    protein.append(AMINO_ACID_1LETTER.get(aa, '?'))
                                translation = ''.join(protein)

                            # Build codon-by-codon alignment
                            codons = []
                            for i in range(0, min(len(coding_seq), 150), 3):
                                dna_codon = coding_seq[i:i+3]
                                rna_codon = dna_codon.replace('T', 'U')
                                aa_idx = i // 3
                                aa = translation[aa_idx] if aa_idx < len(translation) else '?'
                                aa_name = GENETIC_CODE.get(dna_codon, 'Unknown')

                                codons.append({
                                    "position": i + 1,
                                    "dna_codon": dna_codon,
                                    "template_codon": ''.join(complement_map.get(b, 'N') for b in dna_codon),
                                    "rna_codon": rna_codon,
                                    "amino_acid": aa,
                                    "amino_acid_name": aa_name,
                                    "is_start": dna_codon == 'ATG',
                                    "is_stop": dna_codon in ('TAA', 'TAG', 'TGA'),
                                })

                            return {
                                "gene_name": qualifiers.get('gene', [''])[0],
                                "locus_tag": lt,
                                "product": qualifiers.get('product', [''])[0],
                                "strand": strand,
                                "strand_label": "Forward (+)" if strand == 1 else "Reverse (-)",
                                "start": int(feature.location.start),
                                "end": int(feature.location.end),
                                "length_bp": len(coding_seq),
                                "length_aa": len(translation),

                                # Full sequences (first 300 bases shown)
                                "dna_sense_5to3": dna_sense_5to3[:300],
                                "dna_template_3to5": dna_template_3to5[:300],
                                "mrna_5to3": mrna_5to3[:300],
                                "protein": translation[:100],
                                "full_protein": translation,

                                # Codon-by-codon alignment (first 50 codons)
                                "codons": codons,
                                "total_codons": len(coding_seq) // 3,

                                # Process info
                                "gc_content": round(gc_fraction(Seq(coding_seq)) * 100, 2),
                                "molecular_weight_kda": round(len(translation) * 110 / 1000, 1),

                                # NCBI references
                                "protein_id": qualifiers.get('protein_id', [''])[0],
                                "db_xref": qualifiers.get('db_xref', []),
                                "ncbi_protein_url": f"https://www.ncbi.nlm.nih.gov/protein/{qualifiers.get('protein_id', [''])[0]}" if qualifiers.get('protein_id', [''])[0] else None,
                                "ncbi_gene_url": f"https://www.ncbi.nlm.nih.gov/gene/?term={qualifiers.get('gene', [''])[0]}" if qualifiers.get('gene', [''])[0] else None,
                            }

            return None

        except Exception as e:
            print(f"Error getting central dogma data: {e}")
            return None

    def get_gc_sliding_window(self, genbank_path: str, window_size: int = 5000, step: int = 1000):
        """
        Calculate GC content across the genome using a sliding window.
        Returns data points for plotting GC content along the genome.
        """
        try:
            records = list(SeqIO.parse(genbank_path, "genbank"))
            if not records:
                return None
            record = max(records, key=lambda r: len(r.seq))
            seq = str(record.seq).upper()
            genome_length = len(seq)
            avg_gc = round(gc_fraction(record.seq) * 100, 2)

            windows = []
            for start in range(0, genome_length - window_size + 1, step):
                window = seq[start:start + window_size]
                gc = (window.count('G') + window.count('C')) / len(window) * 100
                windows.append({
                    "position": start + window_size // 2,
                    "gc": round(gc, 2),
                    "start": start,
                    "end": start + window_size
                })

            # Also calculate GC skew (G-C)/(G+C) for each window
            gc_skew = []
            for start in range(0, genome_length - window_size + 1, step):
                window = seq[start:start + window_size]
                g = window.count('G')
                c = window.count('C')
                skew = (g - c) / (g + c) if (g + c) > 0 else 0
                gc_skew.append({
                    "position": start + window_size // 2,
                    "skew": round(skew, 4)
                })

            return {
                "genome_length": genome_length,
                "average_gc": avg_gc,
                "window_size": window_size,
                "step": step,
                "total_windows": len(windows),
                "gc_windows": windows,
                "gc_skew": gc_skew,
                "organism": record.annotations.get("organism", ""),
                "accession": record.id,
                "gc_stats": {
                    "min_gc": round(min(w["gc"] for w in windows), 2) if windows else 0,
                    "max_gc": round(max(w["gc"] for w in windows), 2) if windows else 0,
                    "std_gc": round(float(__import__('numpy').std([w["gc"] for w in windows])), 2) if windows else 0,
                }
            }
        except Exception as e:
            print(f"Error calculating GC sliding window: {e}")
            return None

    def submit_blast_search(self, sequence: str, program: str = "blastn",
                           database: str = "nt", max_hits: int = 10):
        """
        Submit a BLAST search to NCBI's BLAST API.
        Returns request ID (RID) for polling results.
        """
        try:
            url = "https://blast.ncbi.nlm.nih.gov/blast/Blast.cgi"
            params = {
                "CMD": "Put",
                "PROGRAM": program,
                "DATABASE": database,
                "QUERY": sequence[:10000],  # Limit to 10kb
                "FORMAT_TYPE": "JSON2",
                "HITLIST_SIZE": max_hits,
                "EXPECT": "0.01"
            }
            response = requests.post(url, data=params, timeout=30)
            response.raise_for_status()

            # Extract RID from NCBI response
            text = response.text
            rid = ""
            
            # Method 1: Look for "RID = XXXXX" pattern in the QBlastInfo block
            # NCBI format is typically: <!-- QBlastInfoBegin RID = RID_VALUE ... -->
            rid_match = re.search(r'RID\s*=\s*([A-Z0-9_-]+)', text, re.IGNORECASE)
            if rid_match:
                candidate = rid_match.group(1).strip()
                # RID is typically uppercase alphanumeric, sometimes with dashes
                if 8 <= len(candidate) <= 25 and candidate.replace("-", "").replace("_", "").isalnum():
                    rid = candidate
            
            # Method 2: Fallback to line-by-line search if regex was too loose
            if not rid:
                for line in text.split('\n'):
                    if 'RID' in line.upper() and '=' in line:
                        # Strip HTML and extract value
                        clean_line = re.sub(r'<[^>]+>', '', line)
                        if '=' in clean_line:
                            parts = clean_line.split('=')
                            candidate = parts[-1].strip()
                            if 8 <= len(candidate) <= 25 and candidate.replace("-", "").replace("_", "").isalnum():
                                rid = candidate
                                break
            
            # Final validation: RID must NOT contain HTML and must match NCBI pattern
            if not rid or '<' in rid or '>' in rid or '"' in rid:
                print(f"⚠️ [BLAST] Error parsing RID. Response sample: {text[:300]}")
                return {
                    "error": "No se pudo obtener un ID de rastreo (RID) válido de NCBI. Es posible que el servicio esté saturado o la secuencia sea inválida.",
                    "status": "ERROR"
                }

            return {
                "rid": rid,
                "program": program,
                "database": database,
                "status": "SUBMITTED",
                "url": f"https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Get&RID={rid}&FORMAT_TYPE=HTML"
            }
        except Exception as e:
            print(f"Error submitting BLAST search: {e}")
            return {"error": str(e), "status": "ERROR"}

    def get_blast_results(self, rid: str):
        """
        Check BLAST results for a given RID.
        """
        try:
            url = f"https://blast.ncbi.nlm.nih.gov/blast/Blast.cgi"
            params = {
                "CMD": "Get",
                "RID": rid,
                "FORMAT_TYPE": "JSON2"
            }
            response = requests.get(url, params=params, timeout=30)
            response.raise_for_status()

            text = response.text

            # Check if still running
            if "WAITING" in text:
                return {"status": "WAITING", "rid": rid}
            if "READY" not in text and "BlastOutput2" not in text:
                return {"status": "WAITING", "rid": rid}

            # Try to parse JSON results
            try:
                import json
                data = json.loads(text)
                results = []
                if "BlastOutput2" in data:
                    for report in data["BlastOutput2"]:
                        search = report.get("report", {}).get("results", {}).get("search", {})
                        hits = search.get("hits", [])
                        for hit in hits[:10]:
                            desc = hit.get("description", [{}])[0]
                            hsps = hit.get("hsps", [{}])
                            best_hsp = hsps[0] if hsps else {}
                            results.append({
                                "accession": desc.get("accession", ""),
                                "title": desc.get("title", "")[:200],
                                "taxid": desc.get("taxid", ""),
                                "score": best_hsp.get("bit_score", 0),
                                "evalue": best_hsp.get("evalue", 0),
                                "identity_pct": round(best_hsp.get("identity", 0) / max(best_hsp.get("align_len", 1), 1) * 100, 1),
                                "align_len": best_hsp.get("align_len", 0),
                                "query_from": best_hsp.get("query_from", 0),
                                "query_to": best_hsp.get("query_to", 0),
                                "url": f"https://www.ncbi.nlm.nih.gov/nucleotide/{desc.get('accession', '')}"
                            })
                return {
                    "status": "COMPLETE",
                    "rid": rid,
                    "total_hits": len(results),
                    "hits": results,
                    "blast_url": f"https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Get&RID={rid}&FORMAT_TYPE=HTML"
                }
            except:
                return {
                    "status": "COMPLETE",
                    "rid": rid,
                    "total_hits": 0,
                    "hits": [],
                    "blast_url": f"https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Get&RID={rid}&FORMAT_TYPE=HTML",
                    "note": "Results are available at the BLAST URL"
                }
        except Exception as e:
            print(f"Error getting BLAST results: {e}")
            return {"status": "ERROR", "error": str(e), "rid": rid}

    def get_trna_rrna_analysis(self, genbank_path: str):
        """
        Extract and analyze tRNA and rRNA features from GenBank file.
        """
        try:
            records = list(SeqIO.parse(genbank_path, "genbank"))
            if not records:
                return None
            record = max(records, key=lambda r: len(r.seq))
            trna_list = []
            rrna_list = []

            for feature in record.features:
                qualifiers = feature.qualifiers
                start = int(feature.location.start)
                end = int(feature.location.end)
                strand = feature.location.strand

                if feature.type == "tRNA":
                    trna_list.append({
                        "product": qualifiers.get("product", ["unknown tRNA"])[0],
                        "locus_tag": qualifiers.get("locus_tag", [""])[0],
                        "anticodon": qualifiers.get("note", [""])[0] if "anticodon" not in qualifiers else qualifiers["anticodon"][0],
                        "start": start,
                        "end": end,
                        "strand": strand,
                        "length": end - start,
                        "amino_acid": qualifiers.get("product", [""])[0].replace("tRNA-", "")
                    })
                elif feature.type == "rRNA":
                    rrna_list.append({
                        "product": qualifiers.get("product", ["unknown rRNA"])[0],
                        "locus_tag": qualifiers.get("locus_tag", [""])[0],
                        "start": start,
                        "end": end,
                        "strand": strand,
                        "length": end - start,
                    })

            # Count amino acid coverage
            aa_coverage = {}
            for trna in trna_list:
                aa = trna["amino_acid"]
                aa_coverage[aa] = aa_coverage.get(aa, 0) + 1

            # rRNA types
            rrna_types = {}
            for rrna in rrna_list:
                prod = rrna["product"]
                rrna_types[prod] = rrna_types.get(prod, 0) + 1

            return {
                "organism": record.annotations.get("organism", ""),
                "total_trna": len(trna_list),
                "total_rrna": len(rrna_list),
                "trna_genes": trna_list,
                "rrna_genes": rrna_list,
                "amino_acid_coverage": aa_coverage,
                "rrna_types": rrna_types,
                "trna_by_strand": {
                    "forward": len([t for t in trna_list if t["strand"] == 1]),
                    "reverse": len([t for t in trna_list if t["strand"] == -1])
                },
                "rrna_by_strand": {
                    "forward": len([r for r in rrna_list if r["strand"] == 1]),
                    "reverse": len([r for r in rrna_list if r["strand"] == -1])
                }
            }
        except Exception as e:
            print(f"Error analyzing tRNA/rRNA: {e}")
            return None

    def get_functional_categories(self, genbank_path: str):
        """
        Classify genes into COG-like functional categories based on product descriptions.
        Uses keyword matching against known functional category terms.
        """
        try:
            records = list(SeqIO.parse(genbank_path, "genbank"))
            if not records:
                return None
            record = max(records, key=lambda r: len(r.seq))

            # COG functional categories with keywords
            COG_CATEGORIES = {
                "J": {"name": "Traducción, ribosomas", "keywords": ["ribosom", "trna", "rrna", "translation", "aminoacyl", "elongation factor", "initiation factor"], "color": "#ef4444"},
                "K": {"name": "Transcripción", "keywords": ["transcription", "rna polymerase", "sigma factor", "anti-sigma", "transcriptional"], "color": "#f97316"},
                "L": {"name": "Replicación y reparación DNA", "keywords": ["dna polymerase", "helicase", "ligase", "topoisomerase", "recombina", "repair", "replicat", "gyrase", "primase"], "color": "#f59e0b"},
                "D": {"name": "División celular", "keywords": ["cell division", "ftsz", "ftsa", "septum", "chromosome partition", "min"], "color": "#eab308"},
                "V": {"name": "Defensa", "keywords": ["restriction", "modification", "crispr", "toxin", "antitoxin", "defense", "immunity"], "color": "#84cc16"},
                "T": {"name": "Transducción de señales", "keywords": ["sensor", "signal", "kinase", "phosphatase", "response regulator", "two-component", "chemotaxis"], "color": "#22c55e"},
                "M": {"name": "Pared celular / membrana", "keywords": ["membrane", "lipopolysaccharide", "peptidoglycan", "outer membrane", "murein", "lps", "porin", "lipoprotein"], "color": "#10b981"},
                "N": {"name": "Motilidad celular", "keywords": ["flagell", "motility", "pilus", "fimbri", "chemotaxis receptor"], "color": "#14b8a6"},
                "U": {"name": "Secreción / transporte vesicular", "keywords": ["secretion", "type ii", "type iii", "type iv", "sec", "tat"], "color": "#06b6d4"},
                "O": {"name": "Chaperonas / modificación post-traduccional", "keywords": ["chaperone", "protease", "peptidase", "heat shock", "groel", "dnak", "clp", "lon protease", "grpe"], "color": "#0ea5e9"},
                "C": {"name": "Metabolismo energético", "keywords": ["dehydrogenase", "oxidoreductase", "cytochrome", "nadh", "atp synthase", "electron transport", "oxidase", "reductase", "ferredoxin"], "color": "#3b82f6"},
                "G": {"name": "Metabolismo de carbohidratos", "keywords": ["sugar", "glucose", "galactose", "phosphotransferase", "glycolysis", "gluconate", "maltose", "lactose", "xylose", "fructose"], "color": "#6366f1"},
                "E": {"name": "Metabolismo de aminoácidos", "keywords": ["amino acid", "aminotransferase", "synthase", "deaminase", "biosynthesis of", "tryptophan", "leucine", "histidine", "arginine"], "color": "#8b5cf6"},
                "F": {"name": "Metabolismo de nucleótidos", "keywords": ["purine", "pyrimidine", "nucleotide", "nucleoside", "thymidylate"], "color": "#a855f7"},
                "H": {"name": "Metabolismo de coenzimas", "keywords": ["coenzyme", "cofactor", "biotin", "thiamine", "folate", "riboflavin", "nad", "cobalamin", "pantothenate"], "color": "#d946ef"},
                "I": {"name": "Metabolismo de lípidos", "keywords": ["fatty acid", "lipid", "acyl", "phospholipid", "lipase"], "color": "#ec4899"},
                "P": {"name": "Transporte de iones inorgánicos", "keywords": ["iron", "zinc", "copper", "magnesium", "phosphate transport", "sulfate transport", "potassium", "sodium"], "color": "#f43f5e"},
                "Q": {"name": "Biosíntesis metabolitos secundarios", "keywords": ["siderophore", "enterobactin", "polyketide", "secondary metabol"], "color": "#fb7185"},
                "R": {"name": "Función general predicha", "keywords": ["predicted", "putative", "probable", "uncharacterized protein"], "color": "#94a3b8"},
                "S": {"name": "Función desconocida", "keywords": ["hypothetical", "unknown function", "duf", "domain of unknown"], "color": "#cbd5e1"},
            }

            categories = {code: {"code": code, **info, "genes": [], "count": 0} for code, info in COG_CATEGORIES.items()}
            uncategorized = []
            total_cds = 0

            for feature in record.features:
                if feature.type != "CDS":
                    continue
                total_cds += 1
                qualifiers = feature.qualifiers
                product = qualifiers.get("product", ["hypothetical protein"])[0].lower()
                locus_tag = qualifiers.get("locus_tag", [""])[0]
                gene_name = qualifiers.get("gene", [""])[0]

                assigned = False
                for code, info in COG_CATEGORIES.items():
                    for kw in info["keywords"]:
                        if kw in product:
                            categories[code]["genes"].append({
                                "locus_tag": locus_tag,
                                "gene_name": gene_name,
                                "product": qualifiers.get("product", [""])[0],
                            })
                            categories[code]["count"] += 1
                            assigned = True
                            break
                    if assigned:
                        break

                if not assigned:
                    uncategorized.append({
                        "locus_tag": locus_tag,
                        "gene_name": gene_name,
                        "product": qualifiers.get("product", [""])[0],
                    })

            # Sort by count
            sorted_cats = sorted(
                [v for v in categories.values() if v["count"] > 0],
                key=lambda x: x["count"],
                reverse=True
            )

            # Limit genes list to 20 per category for response size
            for cat in sorted_cats:
                cat["sample_genes"] = cat["genes"][:20]
                del cat["genes"]

            return {
                "organism": record.annotations.get("organism", ""),
                "total_cds": total_cds,
                "categorized": sum(c["count"] for c in sorted_cats),
                "uncategorized": len(uncategorized),
                "categories": sorted_cats,
                "uncategorized_sample": uncategorized[:20],
                "top_categories": [
                    {"code": c["code"], "name": c["name"], "count": c["count"], "color": c["color"]}
                    for c in sorted_cats[:10]
                ]
            }
        except Exception as e:
            print(f"Error classifying genes: {e}")
            return None

    def calculate_cai_per_gene(self, genbank_path: str, top_n: int = 50):
        """
        Calculate Codon Adaptation Index (CAI) for each gene.
        Uses the codon usage of highly expressed genes (ribosomal proteins)
        as the reference set.
        """
        try:
            import math
            records = list(SeqIO.parse(genbank_path, "genbank"))
            if not records:
                return None
            record = max(records, key=lambda r: len(r.seq))

            # Step 1: Build reference codon usage from highly expressed genes
            # (ribosomal proteins are highly expressed in all organisms)
            ref_codons = {}
            total_ref_codons = 0
            all_genes_data = []

            for feature in record.features:
                if feature.type != "CDS":
                    continue
                qualifiers = feature.qualifiers
                product = qualifiers.get("product", [""])[0].lower()
                locus_tag = qualifiers.get("locus_tag", [""])[0]
                gene_name = qualifiers.get("gene", [""])[0]

                try:
                    seq = str(feature.location.extract(record.seq)).upper()
                except:
                    continue

                if len(seq) < 30 or len(seq) % 3 != 0:
                    continue

                # Count codons for this gene
                gene_codons = {}
                for i in range(0, len(seq) - 2, 3):
                    codon = seq[i:i+3]
                    if 'N' not in codon:
                        gene_codons[codon] = gene_codons.get(codon, 0) + 1

                all_genes_data.append({
                    "locus_tag": locus_tag,
                    "gene_name": gene_name,
                    "product": qualifiers.get("product", [""])[0],
                    "length": len(seq),
                    "strand": feature.location.strand,
                    "codons": gene_codons
                })

                # If ribosomal protein, add to reference set
                is_ref = any(kw in product for kw in [
                    "ribosomal protein", "elongation factor", "translation"
                ])
                if is_ref:
                    for codon, count in gene_codons.items():
                        ref_codons[codon] = ref_codons.get(codon, 0) + count
                        total_ref_codons += count

            # If not enough reference genes, use all genes as reference
            if total_ref_codons < 1000:
                ref_codons = {}
                total_ref_codons = 0
                for g in all_genes_data:
                    for codon, count in g["codons"].items():
                        ref_codons[codon] = ref_codons.get(codon, 0) + count
                        total_ref_codons += count

            # Step 2: Calculate relative adaptiveness (w) for each codon
            # Group codons by amino acid
            CODON_TABLE = {
                'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
                'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
                'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
                'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
                'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
                'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
                'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
                'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
                'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
                'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
                'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
                'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
                'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
                'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
                'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
                'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
            }

            # Group by amino acid
            aa_codons = {}
            for codon, aa in CODON_TABLE.items():
                if aa == '*':
                    continue
                if aa not in aa_codons:
                    aa_codons[aa] = []
                aa_codons[aa].append(codon)

            # Relative adaptiveness w(codon) = freq(codon) / max_freq(synonyms)
            w_values = {}
            for aa, codons in aa_codons.items():
                freqs = {c: ref_codons.get(c, 0) for c in codons}
                max_freq = max(freqs.values()) if freqs else 1
                if max_freq == 0:
                    max_freq = 1
                for c in codons:
                    w_values[c] = freqs[c] / max_freq if freqs[c] > 0 else 0.01

            # Step 3: Calculate CAI for each gene
            results = []
            for gene in all_genes_data:
                log_sum = 0
                n = 0
                for codon, count in gene["codons"].items():
                    if codon in w_values and CODON_TABLE.get(codon, '*') != '*':
                        w = w_values[codon]
                        if w > 0:
                            log_sum += count * math.log(w)
                            n += count

                cai = math.exp(log_sum / n) if n > 0 else 0

                results.append({
                    "locus_tag": gene["locus_tag"],
                    "gene_name": gene["gene_name"],
                    "product": gene["product"],
                    "length": gene["length"],
                    "strand": gene["strand"],
                    "cai": round(cai, 4),
                })

            # Sort by CAI descending
            results.sort(key=lambda x: x["cai"], reverse=True)

            # Statistics
            cai_values = [r["cai"] for r in results if r["cai"] > 0]
            import numpy as np

            return {
                "organism": record.annotations.get("organism", ""),
                "total_genes": len(results),
                "reference_genes": total_ref_codons,
                "cai_stats": {
                    "mean": round(float(np.mean(cai_values)), 4) if cai_values else 0,
                    "median": round(float(np.median(cai_values)), 4) if cai_values else 0,
                    "min": round(min(cai_values), 4) if cai_values else 0,
                    "max": round(max(cai_values), 4) if cai_values else 0,
                    "std": round(float(np.std(cai_values)), 4) if cai_values else 0,
                },
                "top_expressed": results[:top_n],
                "low_expressed": results[-top_n:] if len(results) > top_n else [],
                "all_genes": results,
                "cai_distribution": {
                    "0.0-0.2": len([c for c in cai_values if c < 0.2]),
                    "0.2-0.4": len([c for c in cai_values if 0.2 <= c < 0.4]),
                    "0.4-0.6": len([c for c in cai_values if 0.4 <= c < 0.6]),
                    "0.6-0.8": len([c for c in cai_values if 0.6 <= c < 0.8]),
                    "0.8-1.0": len([c for c in cai_values if c >= 0.8]),
                }
            }
        except Exception as e:
            print(f"Error calculating CAI: {e}")
            import traceback
            traceback.print_exc()
            return None

    def calculate_phylogenetic_distances(self, genbank_paths: list):
        """
        Calculate pairwise distances between multiple genomes for phylogenetic tree.
        Uses Average Nucleotide Identity approximation based on:
        - GC content difference
        - Genome size ratio
        - Shared gene content (Jaccard)
        Returns distance matrix for dendrogram visualization.
        """
        try:
            import numpy as np
            genomes = []

            for path in genbank_paths:
                try:
                    records = list(SeqIO.parse(path, "genbank"))
                    if not records:
                        continue
                    # Pick largest record (chromosome)
                    record = max(records, key=lambda r: len(r.seq))
                    gc = round(gc_fraction(record.seq) * 100, 2)
                    genes = set()
                    gene_products = set()
                    for f in record.features:
                        if f.type == "CDS":
                            gene = f.qualifiers.get("gene", [""])[0]
                            product = f.qualifiers.get("product", [""])[0].lower()
                            if gene:
                                genes.add(gene.lower())
                            if product and "hypothetical" not in product:
                                # Use simplified product name for matching
                                gene_products.add(product[:50])

                    genomes.append({
                        "name": record.annotations.get("organism", record.id),
                        "accession": record.id,
                        "gc": gc,
                        "length": len(record.seq),
                        "gene_count": len(genes),
                        "genes": genes,
                        "products": gene_products,
                    })
                except Exception as e:
                    print(f"Error parsing {path}: {e}")

            n = len(genomes)
            if n < 2:
                return {"error": "Se necesitan al menos 2 genomas para construir el árbol"}

            # Calculate distance matrix
            dist_matrix = np.zeros((n, n))

            for i in range(n):
                for j in range(i + 1, n):
                    # GC distance (0-50 range normalized to 0-1)
                    gc_dist = abs(genomes[i]["gc"] - genomes[j]["gc"]) / 50.0

                    # Size distance (ratio, 0-1)
                    max_len = max(genomes[i]["length"], genomes[j]["length"])
                    min_len = min(genomes[i]["length"], genomes[j]["length"])
                    size_dist = 1 - (min_len / max_len) if max_len > 0 else 1

                    # Gene content distance (1 - Jaccard similarity)
                    union = genomes[i]["genes"] | genomes[j]["genes"]
                    intersection = genomes[i]["genes"] & genomes[j]["genes"]
                    jaccard_genes = len(intersection) / len(union) if union else 0

                    # Product similarity
                    prod_union = genomes[i]["products"] | genomes[j]["products"]
                    prod_intersect = genomes[i]["products"] & genomes[j]["products"]
                    jaccard_products = len(prod_intersect) / len(prod_union) if prod_union else 0

                    # Combined distance (weighted)
                    distance = (
                        0.3 * gc_dist +
                        0.1 * size_dist +
                        0.3 * (1 - jaccard_genes) +
                        0.3 * (1 - jaccard_products)
                    )

                    dist_matrix[i][j] = round(distance, 4)
                    dist_matrix[j][i] = round(distance, 4)

            # UPGMA clustering for dendrogram
            labels = [g["name"] for g in genomes]
            tree = self._upgma_cluster(dist_matrix.tolist(), labels)

            return {
                "genomes": [
                    {
                        "name": g["name"],
                        "accession": g["accession"],
                        "gc": g["gc"],
                        "length": g["length"],
                        "gene_count": g["gene_count"],
                    }
                    for g in genomes
                ],
                "distance_matrix": dist_matrix.tolist(),
                "labels": labels,
                "tree": tree,
                "method": "UPGMA (GC + Genome Size + Gene Content + Product Similarity)"
            }
        except Exception as e:
            print(f"Error calculating phylogenetic distances: {e}")
            import traceback
            traceback.print_exc()
            return None

    def _upgma_cluster(self, dist_matrix, labels):
        """Simple UPGMA clustering returning Newick-like tree structure."""
        import copy
        n = len(labels)
        if n <= 1:
            return {"name": labels[0] if labels else "", "children": []}

        # Copy matrix and labels
        matrix = copy.deepcopy(dist_matrix)
        nodes = [{"name": label, "height": 0} for label in labels]
        sizes = [1] * n

        while len(nodes) > 1:
            # Find minimum distance
            min_dist = float('inf')
            mi, mj = 0, 1
            for i in range(len(matrix)):
                for j in range(i + 1, len(matrix)):
                    if matrix[i][j] < min_dist:
                        min_dist = matrix[i][j]
                        mi, mj = i, j

            # Create new node
            height = min_dist / 2
            new_node = {
                "name": f"({nodes[mi]['name']}, {nodes[mj]['name']})",
                "children": [
                    {**nodes[mi], "branch_length": round(height - nodes[mi]["height"], 4)},
                    {**nodes[mj], "branch_length": round(height - nodes[mj]["height"], 4)},
                ],
                "height": height,
                "distance": round(min_dist, 4)
            }

            # Update matrix (UPGMA formula)
            new_row = []
            size_mi = sizes[mi]
            size_mj = sizes[mj]
            
            for k in range(len(matrix)):
                if k == mi or k == mj:
                    continue
                d = (matrix[mi][k] * size_mi + matrix[mj][k] * size_mj) / (size_mi + size_mj)
                new_row.append(round(d, 4))

            # Remove old entries (higher index first)
            for idx in sorted([mi, mj], reverse=True):
                matrix.pop(idx)
                nodes.pop(idx)
                sizes.pop(idx)

            # Remove columns
            for row in matrix:
                for idx in sorted([mi, mj], reverse=True):
                    row.pop(idx)

            # Add new entries
            for i, row in enumerate(matrix):
                row.append(new_row[i])
            matrix.append(new_row + [0])
            nodes.append(new_node)
            sizes.append(size_mi + size_mj)

        return nodes[0] if nodes else {}


# Singleton
_ncbi_service: Optional[NCBIService] = None


def get_ncbi_service(api_key: Optional[str] = None) -> NCBIService:
    global _ncbi_service
    if _ncbi_service is None:
        _ncbi_service = NCBIService(api_key)
    return _ncbi_service
