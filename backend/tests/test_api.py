"""
Tests Completos — Sistema de Análisis Genómico
Verifica backend + nuevas funcionalidades contra genomas reales.

Cubre:
- Tests básicos (health, files, genome)
- Activar genomas y verificar
- GC Sliding Window
- tRNA/rRNA Analysis
- Categorías Funcionales (COG)
- CAI per Gene
- BLAST (submit only, no polling)
- Phylogenetic Tree (multi-genome)
- Genome Summary
- Proteins, Codon Usage, Genome Map

Usa los 3 genomas disponibles:
  1. GCA_000005845.2 — E. coli K-12 MG1655
  2. GCA_000008865.2 — E. coli O157:H7 Sakai
  3. GCF_020097475.1 — E. fergusonii
"""
import pytest
import sys
import os
import time

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from fastapi.testclient import TestClient
from app.main import app

# === FIXTURE FOR APP STATE ===
@pytest.fixture(autouse=True)
def setup_app_state():
    """Initialize app state (FileDetector) for tests"""
    from app.core.file_detector import FileDetector
    project_root = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    app.state.project_root = project_root
    app.state.file_detector = FileDetector(project_root)
    
    # Pre-scan the genomes directory
    genomes_dir = os.path.join(project_root, "genomes")
    if os.path.exists(genomes_dir):
        app.state.file_detector.scan_directory(genomes_dir)
        
    yield

client = TestClient(app)

# === GenBank file paths ===
PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
GENOMES = {
    "GCA_000005845.2": os.path.join(PROJECT_ROOT, "genomes/GCA_000005845.2/extracted/ncbi_dataset/data/GCA_000005845.2/genomic.gbff"),
    "GCA_000008865.2": os.path.join(PROJECT_ROOT, "genomes/GCA_000008865.2/extracted/ncbi_dataset/data/GCA_000008865.2/genomic.gbff"),
    "GCF_020097475.1": os.path.join(PROJECT_ROOT, "genomes/GCF_020097475.1/extracted/ncbi_dataset/data/GCF_020097475.1/genomic.gbff"),
}


# ==============================================================
#  BASIC TESTS
# ==============================================================

class TestBasicHealth:
    """Health and root endpoints"""

    def test_root(self):
        r = client.get("/")
        assert r.status_code == 200
        assert "E. coli" in r.json()["message"]

    def test_health(self):
        r = client.get("/health")
        assert r.status_code == 200
        assert r.json()["status"] == "healthy"


class TestFileEndpoints:

    def test_detect_files(self):
        r = client.get("/api/files/detect")
        assert r.status_code == 200
        assert "files" in r.json()


class TestGenomeManagement:

    def test_popular_genomes(self):
        r = client.get("/api/genome/popular")
        assert r.status_code == 200
        assert len(r.json()["genomes"]) > 0

    def test_downloaded_genomes(self):
        r = client.get("/api/genome/downloaded")
        assert r.status_code == 200
        data = r.json()
        assert "genomes" in data
        # Verify we have some genomes
        print(f"  ✓ {len(data['genomes'])} genomas descargados")

    def test_analysis_status(self):
        r = client.get("/api/analysis/status")
        assert r.status_code == 200


# ==============================================================
#  ACTIVATE EACH GENOME AND TEST CORE ENDPOINTS
# ==============================================================

class TestGenomeActivation:
    """Activate each genome and verify basic functionality"""

    @pytest.mark.parametrize("accession", list(GENOMES.keys()))
    def test_activate_genome(self, accession):
        """Activate a genome and verify it's active"""
        r = client.post(f"/api/genome/activate/{accession}")
        assert r.status_code in [200, 404], f"Activate {accession}: {r.status_code} - {r.text}"
        if r.status_code == 200:
            data = r.json()
            assert data.get("accession") == accession or "activated" in str(data).lower()
            print(f"  ✓ Activated {accession}")


# ==============================================================
#  NCBI ENHANCED ENDPOINTS — PER GENOME
# ==============================================================

class TestPerGenomeAnalysis:
    """Test analysis endpoints on E. coli K-12 (primary genome)"""

    @classmethod
    def setup_class(cls):
        """Activate E. coli K-12 for testing"""
        r = client.post("/api/genome/activate/GCA_000005845.2")
        if r.status_code != 200:
            pytest.skip("Cannot activate E. coli K-12")

    # ---------- Proteins ----------
    def test_proteins(self):
        r = client.get("/api/ncbi/proteins?page=1&page_size=10")
        assert r.status_code == 200
        data = r.json()
        assert "proteins" in data
        assert len(data["proteins"]) > 0
        assert "locus_tag" in data["proteins"][0]
        print(f"  ✓ Proteins: {data.get('total', len(data['proteins']))} total")

    # ---------- Genome Map ----------
    def test_genome_map(self):
        r = client.get("/api/ncbi/genome-map")
        assert r.status_code == 200
        data = r.json()
        assert "genes" in data
        assert "genome_length" in data
        assert data["genome_length"] > 4_000_000  # E. coli ~4.6 Mb
        print(f"  ✓ Genome Map: {data['genome_length']:,} bp, {len(data['genes'])} genes")

    # ---------- Codon Usage ----------
    def test_codon_usage(self):
        r = client.get("/api/ncbi/codon-usage")
        assert r.status_code == 200
        data = r.json()
        assert "codon_table" in data or "codons" in data
        print(f"  ✓ Codon Usage OK")

    # ---------- GC Sliding Window ----------
    def test_gc_window(self):
        r = client.get("/api/ncbi/gc-window?window_size=5000&step=2500")
        assert r.status_code == 200
        data = r.json()
        assert "gc_windows" in data
        assert "gc_skew" in data
        assert "genome_length" in data
        assert "gc_stats" in data
        assert len(data["gc_windows"]) > 100
        assert 40 < data["average_gc"] < 60  # E. coli GC ~50.8%
        print(f"  ✓ GC Window: {len(data['gc_windows'])} windows, avg GC={data['average_gc']}%")
        print(f"    Stats: min={data['gc_stats']['min_gc']}%, max={data['gc_stats']['max_gc']}%, std={data['gc_stats']['std_gc']}%")

    # ---------- tRNA/rRNA Analysis ----------
    def test_rna_analysis(self):
        r = client.get("/api/ncbi/rna-analysis")
        assert r.status_code == 200
        data = r.json()
        assert "total_trna" in data
        assert "total_rrna" in data
        assert "trna_genes" in data
        assert "rrna_genes" in data
        # E. coli has ~86 tRNAs and ~22 rRNAs
        assert data["total_trna"] > 50
        assert data["total_rrna"] > 10
        assert "amino_acid_coverage" in data
        print(f"  ✓ RNA Analysis: {data['total_trna']} tRNAs, {data['total_rrna']} rRNAs")
        print(f"    AA coverage: {len(data.get('amino_acid_coverage', {}))} amino acids")

    # ---------- Functional Categories (COG) ----------
    def test_functional_categories(self):
        r = client.get("/api/ncbi/functional-categories")
        assert r.status_code == 200
        data = r.json()
        assert "total_cds" in data
        assert "categories" in data
        assert "categorized" in data
        assert len(data["categories"]) > 5  # Should have multiple categories
        assert data["total_cds"] > 3000  # E. coli has >4000 CDS
        assert data["categorized"] > 0
        total_cat = data["categorized"]
        total_cds = data["total_cds"]
        pct = round(total_cat / total_cds * 100, 1)
        print(f"  ✓ COG Categories: {len(data['categories'])} categories, {total_cat}/{total_cds} classified ({pct}%)")
        for cat in data["categories"][:3]:
            print(f"    {cat['code']} ({cat['name']}): {cat['count']} genes")

    # ---------- CAI Analysis ----------
    def test_cai_analysis(self):
        r = client.get("/api/ncbi/cai-analysis?top_n=30")
        assert r.status_code == 200
        data = r.json()
        assert "total_genes" in data
        assert "cai_stats" in data
        assert "top_expressed" in data
        assert "low_expressed" in data
        assert "cai_distribution" in data
        assert data["total_genes"] > 3000
        assert 0 < data["cai_stats"]["mean"] < 1
        assert 0 < data["cai_stats"]["max"] <= 1
        # Top expressed should have high CAI
        if data["top_expressed"]:
            top = data["top_expressed"][0]
            assert top["cai"] > 0.5, f"Top gene CAI too low: {top['cai']}"
        print(f"  ✓ CAI Analysis: {data['total_genes']} genes")
        print(f"    Stats: mean={data['cai_stats']['mean']}, max={data['cai_stats']['max']}, min={data['cai_stats']['min']}")
        print(f"    Distribution: {data['cai_distribution']}")
        if data["top_expressed"]:
            t = data["top_expressed"][0]
            print(f"    Top expressed: {t['locus_tag']} ({t.get('gene_name','')}) CAI={t['cai']} — {t['product']}")

    # ---------- Genome Summary ----------
    def test_genome_summary(self):
        r = client.get("/api/ncbi/genome-summary")
        assert r.status_code == 200
        data = r.json()
        assert "organism" in data
        assert "genome_length" in data
        print(f"  ✓ Genome Summary: {data['organism']} — {data['genome_length']:,} bp")

    # ---------- Gene at Position ----------
    def test_gene_at_position(self):
        r = client.get("/api/ncbi/gene-at-position/1000000")
        assert r.status_code == 200
        data = r.json()
        assert "nearest_gene" in data or "gene" in data or "locus_tag" in data
        print(f"  ✓ Gene at position 1M: {data}")

    # ---------- Central Dogma (pick first protein) ----------
    def test_central_dogma(self):
        # First get a protein locus_tag
        pr = client.get("/api/ncbi/proteins?page=1&page_size=1")
        if pr.status_code != 200:
            pytest.skip("Cannot get proteins")
        proteins = pr.json().get("proteins", [])
        if not proteins:
            pytest.skip("No proteins available")
        locus = proteins[0]["locus_tag"]

        r = client.get(f"/api/ncbi/central-dogma/{locus}")
        assert r.status_code == 200
        data = r.json()
        assert "dna_sense_5to3" in data or "codon_table" in data
        print(f"  ✓ Central Dogma for {locus}: DNA len={len(data.get('dna_sense_5to3', ''))}")


# ==============================================================
#  TEST ON SECOND GENOME (E. coli O157:H7)
# ==============================================================

class TestSecondGenome:
    """Switch to E. coli O157:H7 and verify key endpoints work"""

    @classmethod
    def setup_class(cls):
        r = client.post("/api/genome/activate/GCA_000008865.2")
        if r.status_code != 200:
            pytest.skip("Cannot activate O157:H7")

    def test_proteins_o157(self):
        r = client.get("/api/ncbi/proteins?page=1&page_size=5")
        assert r.status_code == 200
        data = r.json()
        assert len(data.get("proteins", [])) > 0
        print(f"  ✓ O157:H7 Proteins: {data.get('total', '?')} total")

    def test_gc_window_o157(self):
        r = client.get("/api/ncbi/gc-window?window_size=10000&step=5000")
        assert r.status_code == 200
        data = r.json()
        assert len(data.get("gc_windows", [])) > 50
        print(f"  ✓ O157:H7 GC Window: {len(data['gc_windows'])} windows, avg={data['average_gc']}%")

    def test_rna_o157(self):
        r = client.get("/api/ncbi/rna-analysis")
        assert r.status_code == 200
        data = r.json()
        assert data["total_trna"] > 0
        print(f"  ✓ O157:H7 RNA: {data['total_trna']} tRNAs, {data['total_rrna']} rRNAs")

    def test_cog_o157(self):
        r = client.get("/api/ncbi/functional-categories")
        assert r.status_code == 200
        data = r.json()
        assert data["total_cds"] > 4000  # O157:H7 has ~5113 CDS
        print(f"  ✓ O157:H7 COG: {data['categorized']}/{data['total_cds']} classified")

    def test_cai_o157(self):
        r = client.get("/api/ncbi/cai-analysis?top_n=10")
        assert r.status_code == 200
        data = r.json()
        assert data["total_genes"] > 3000
        print(f"  ✓ O157:H7 CAI: mean={data['cai_stats']['mean']}, max={data['cai_stats']['max']}")


# ==============================================================
#  TEST ON THIRD GENOME (E. fergusonii)
# ==============================================================

class TestThirdGenome:
    """Switch to E. fergusonii and verify"""

    @classmethod
    def setup_class(cls):
        r = client.post("/api/genome/activate/GCF_020097475.1")
        if r.status_code != 200:
            pytest.skip("Cannot activate E. fergusonii")

    def test_proteins_fergusonii(self):
        r = client.get("/api/ncbi/proteins?page=1&page_size=5")
        assert r.status_code == 200
        print(f"  ✓ E. fergusonii Proteins: {r.json().get('total', '?')} total")

    def test_gc_window_fergusonii(self):
        r = client.get("/api/ncbi/gc-window?window_size=5000&step=2500")
        assert r.status_code == 200
        data = r.json()
        assert len(data.get("gc_windows", [])) > 50
        print(f"  ✓ E. fergusonii GC: avg={data['average_gc']}%")

    def test_cai_fergusonii(self):
        r = client.get("/api/ncbi/cai-analysis?top_n=10")
        assert r.status_code == 200
        data = r.json()
        assert data["total_genes"] > 2000
        print(f"  ✓ E. fergusonii CAI: mean={data['cai_stats']['mean']}")

    def test_cog_fergusonii(self):
        r = client.get("/api/ncbi/functional-categories")
        assert r.status_code == 200
        data = r.json()
        assert data["total_cds"] > 2000
        print(f"  ✓ E. fergusonii COG: {data['categorized']}/{data['total_cds']} classified")


# ==============================================================
#  PHYLOGENETIC TREE (MULTI-GENOME)
# ==============================================================

class TestPhylogeneticTree:
    """Test phylogenetic tree with all available genomes"""

    def test_phylo_with_paths(self):
        """Build tree with explicit GenBank paths"""
        available_paths = [p for p in GENOMES.values() if os.path.exists(p)]
        if len(available_paths) < 2:
            pytest.skip("Need >=2 genomes for phylogeny")

        r = client.post("/api/ncbi/phylogenetic-tree", json={
            "genbank_paths": available_paths
        })
        assert r.status_code == 200, f"Phylo tree failed: {r.status_code} - {r.text}"
        data = r.json()
        assert "tree" in data
        assert "distance_matrix" in data
        assert "labels" in data
        assert "genomes" in data
        assert len(data["genomes"]) >= 2
        assert len(data["distance_matrix"]) == len(data["genomes"])
        print(f"  ✓ Phylogenetic Tree: {len(data['genomes'])} genomes")
        for g in data["genomes"]:
            print(f"    — {g['name']} ({g['accession']}) GC={g['gc']}%, {g['length']:,} bp")
        print(f"    Method: {data['method']}")
        print(f"    Distance matrix ({len(data['labels'])}x{len(data['labels'])}): {data['distance_matrix']}")

    def test_phylo_tree_structure(self):
        """Verify tree structure has proper children"""
        available_paths = [p for p in GENOMES.values() if os.path.exists(p)]
        if len(available_paths) < 2:
            pytest.skip("Need >=2 genomes")

        r = client.post("/api/ncbi/phylogenetic-tree", json={
            "genbank_paths": available_paths
        })
        data = r.json()
        tree = data["tree"]
        # Tree should have children
        assert "children" in tree or "name" in tree
        # Verify distances are within [0, 1]
        for row in data["distance_matrix"]:
            for val in row:
                assert 0 <= val <= 1.5, f"Distance out of range: {val}"
        print(f"  ✓ Tree structure valid, distances in range")


# ==============================================================
#  BLAST (Submit only — no waiting)
# ==============================================================

class TestBLASTSubmit:
    """Test BLAST submission endpoint (not results, which require waiting)"""

    @classmethod
    def setup_class(cls):
        # Re-activate K-12
        client.post("/api/genome/activate/GCA_000005845.2")

    def test_blast_submit_format(self):
        """Verify BLAST submit endpoint accepts proper format"""
        # Use a short test sequence from E. coli
        test_seq = "ATGAAACGCATTAGCACCACCATTACCACCACCATCACCATTACCACAGGT"

        r = client.post("/api/ncbi/blast/submit", json={
            "sequence": test_seq,
            "program": "blastn",
            "database": "nt",
            "max_hits": 5
        })
        # May fail due to network but should be properly formatted
        if r.status_code == 200:
            data = r.json()
            assert "rid" in data
            print(f"  ✓ BLAST submitted: RID={data['rid']}")
        else:
            print(f"  ⚠ BLAST submit returned {r.status_code} (network dependent)")
            # Accept 422, 500 for network issues but not 404 (missing endpoint)
            assert r.status_code != 404, "BLAST endpoint not found!"


# ==============================================================
#  UNIT TESTS — Core Analyzers
# ==============================================================

class TestCoreAnalyzers:

    def test_codon_analyzer(self):
        from app.core.analyzers.codon_analyzer import CodonAnalyzer
        analyzer = CodonAnalyzer()
        result = analyzer.analyze("ATGATGATG", 0)
        assert result.atg_count == 3

    def test_gene_analyzer(self):
        from app.core.analyzers.gene_analyzer import GeneAnalyzer
        analyzer = GeneAnalyzer()
        assert analyzer is not None

    def test_validator(self):
        from app.core.analyzers.validator import Validator
        validator = Validator()
        calculated = {
            "total_genes": 4651,
            "genome_length": 4641652,
            "gc_content": 50.79,
            "total_cds": 4318,
            "gene_density": 1002.0
        }
        result = validator.validate(calculated)
        assert result.overall_status.value == "PASS"


# ==============================================================
#  UNIT TESTS — NCBI Service Methods
# ==============================================================

class TestNCBIServiceDirect:
    """Test NCBI service methods directly"""

    def _get_k12_path(self):
        path = GENOMES["GCA_000005845.2"]
        if not os.path.exists(path):
            pytest.skip("K-12 GenBank not found")
        return path

    def test_functional_categories_direct(self):
        from app.core.ncbi_service import NCBIService
        svc = NCBIService()
        result = svc.get_functional_categories(self._get_k12_path())
        assert result is not None
        assert result["total_cds"] > 3000
        assert len(result["categories"]) > 5
        print(f"  ✓ Direct COG: {result['categorized']}/{result['total_cds']}")

    def test_cai_per_gene_direct(self):
        from app.core.ncbi_service import NCBIService
        svc = NCBIService()
        result = svc.calculate_cai_per_gene(self._get_k12_path(), top_n=20)
        assert result is not None
        assert result["total_genes"] > 3000
        assert 0 < result["cai_stats"]["mean"] < 1
        # Top expressed gene should be ribosomal or housekeeping
        top = result["top_expressed"][0]
        assert top["cai"] > 0.5
        print(f"  ✓ Direct CAI: top gene {top['locus_tag']} ({top.get('gene_name','')}) CAI={top['cai']}")

    def test_gc_sliding_window_direct(self):
        from app.core.ncbi_service import NCBIService
        svc = NCBIService()
        result = svc.get_gc_sliding_window(self._get_k12_path(), window_size=10000, step=5000)
        assert result is not None
        assert len(result["gc_windows"]) > 100
        assert 40 < result["average_gc"] < 60
        print(f"  ✓ Direct GC: {len(result['gc_windows'])} windows, avg={result['average_gc']}%")

    def test_trna_rrna_direct(self):
        from app.core.ncbi_service import NCBIService
        svc = NCBIService()
        result = svc.get_trna_rrna_analysis(self._get_k12_path())
        assert result is not None
        assert result["total_trna"] > 50
        assert result["total_rrna"] > 10
        print(f"  ✓ Direct RNA: {result['total_trna']} tRNAs, {result['total_rrna']} rRNAs")

    def test_phylo_distances_direct(self):
        from app.core.ncbi_service import NCBIService
        svc = NCBIService()
        paths = [p for p in GENOMES.values() if os.path.exists(p)]
        if len(paths) < 2:
            pytest.skip("Need >=2 genomes")
        result = svc.calculate_phylogenetic_distances(paths)
        assert result is not None
        assert "tree" in result
        assert "distance_matrix" in result
        print(f"  ✓ Direct Phylo: {len(result['genomes'])} genomes clustered")

    def test_upgma_simple(self):
        """Test UPGMA with a simple distance matrix"""
        from app.core.ncbi_service import NCBIService
        svc = NCBIService()
        # Simple 3x3 matrix
        matrix = [
            [0, 0.1, 0.5],
            [0.1, 0, 0.4],
            [0.5, 0.4, 0],
        ]
        labels = ["A", "B", "C"]
        tree = svc._upgma_cluster(matrix, labels)
        assert "children" in tree
        assert tree["children"] is not None
        print(f"  ✓ UPGMA: {tree['name']}")


if __name__ == "__main__":
    pytest.main(["-v", "-s", "--tb=short", __file__])
