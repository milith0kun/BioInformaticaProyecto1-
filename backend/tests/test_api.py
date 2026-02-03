"""
Tests para el Sistema de Análisis Genómico E. coli K-12
"""
import pytest
from fastapi.testclient import TestClient
import sys
import os

# Add parent to path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from app.main import app

client = TestClient(app)


class TestHealthEndpoints:
    """Tests para endpoints de salud"""
    
    def test_root(self):
        """Test endpoint raíz"""
        response = client.get("/")
        assert response.status_code == 200
        data = response.json()
        assert "message" in data
        assert "E. coli" in data["message"]
    
    def test_health(self):
        """Test health check"""
        response = client.get("/health")
        assert response.status_code == 200
        assert response.json()["status"] == "healthy"


class TestFileEndpoints:
    """Tests para endpoints de archivos"""
    
    def test_detect_files(self):
        """Test detección de archivos"""
        response = client.get("/api/files/detect")
        assert response.status_code == 200
        data = response.json()
        assert "files" in data
        assert "count" in data
    
    def test_extract_zip(self):
        """Test extracción de ZIP"""
        response = client.post("/api/files/extract")
        # Puede fallar si no hay ZIP, pero debe responder
        assert response.status_code in [200, 404]


class TestGenomeEndpoints:
    """Tests para endpoints de genoma"""
    
    def test_popular_genomes(self):
        """Test lista de genomas populares"""
        response = client.get("/api/genome/popular")
        assert response.status_code == 200
        data = response.json()
        assert "genomes" in data
        assert len(data["genomes"]) > 0
        # Verificar que E. coli K-12 está en la lista
        accessions = [g["accession"] for g in data["genomes"]]
        assert "GCF_000005845.2" in accessions
    
    def test_downloaded_genomes(self):
        """Test lista de genomas descargados"""
        response = client.get("/api/genome/downloaded")
        assert response.status_code == 200
        data = response.json()
        assert "genomes" in data
    
    def test_search_genomes_valid(self):
        """Test búsqueda con término válido"""
        response = client.get("/api/genome/search/Escherichia%20coli?limit=5")
        # Puede ser 200 o 404 si no hay conexión a NCBI
        assert response.status_code in [200, 404, 500]
        if response.status_code == 200:
            data = response.json()
            assert "results" in data
    
    def test_genome_info_valid(self):
        """Test info de genoma conocido"""
        response = client.get("/api/genome/info/GCF_000005845.2")
        # Puede fallar si no hay conexión a NCBI
        assert response.status_code in [200, 404, 500]


class TestAnalysisEndpoints:
    """Tests para endpoints de análisis"""
    
    def test_analysis_status(self):
        """Test estado del análisis"""
        response = client.get("/api/analysis/status")
        assert response.status_code == 200
        data = response.json()
        assert "status" in data
    
    def test_clear_cache(self):
        """Test limpiar caché"""
        response = client.post("/api/analysis/clear-cache")
        assert response.status_code == 200
        data = response.json()
        assert "message" in data


class TestExportEndpoints:
    """Tests para endpoints de exportación"""
    
    def test_export_json_no_data(self):
        """Test exportar JSON sin datos previos"""
        response = client.get("/api/export/json")
        # Debe fallar si no hay datos
        assert response.status_code in [200, 404]
    
    def test_export_csv_no_data(self):
        """Test exportar CSV sin datos"""
        response = client.get("/api/export/csv/genes")
        assert response.status_code in [200, 404]


class TestCodonAnalyzer:
    """Tests unitarios para el analizador de codones"""
    
    def test_codon_analyzer_import(self):
        """Test importación del módulo"""
        from app.core.analyzers.codon_analyzer import CodonAnalyzer
        analyzer = CodonAnalyzer()
        assert analyzer is not None
    
    def test_find_atg_simple(self):
        """Test búsqueda de ATG en secuencia simple"""
        from app.core.analyzers.codon_analyzer import CodonAnalyzer
        analyzer = CodonAnalyzer()
        
        sequence = "ATGATGATG"  # 3 ATG
        result = analyzer.analyze(sequence, 0)
        assert result.atg_count == 3
    
    def test_stop_codon_detection(self):
        """Test detección de codones de parada"""
        from app.core.analyzers.codon_analyzer import CodonAnalyzer
        analyzer = CodonAnalyzer()
        
        # TAA, TAG, TGA
        sequence = "TAATAGTGA"
        result = analyzer.analyze(sequence, 0)
        
        assert result.stop_codons["TAA"]["count"] >= 1
        assert result.stop_codons["TAG"]["count"] >= 1
        assert result.stop_codons["TGA"]["count"] >= 1


class TestGeneAnalyzer:
    """Tests unitarios para el analizador de genes"""
    
    def test_gene_analyzer_import(self):
        """Test importación del módulo"""
        from app.core.analyzers.gene_analyzer import GeneAnalyzer
        analyzer = GeneAnalyzer()
        assert analyzer is not None


class TestValidator:
    """Tests para el validador de resultados"""
    
    def test_validator_import(self):
        """Test importación del módulo"""
        from app.core.analyzers.validator import Validator
        validator = Validator()
        assert validator is not None
    
    def test_validation_ecoli_reference(self):
        """Test validación contra valores de referencia de E. coli"""
        from app.core.analyzers.validator import Validator
        validator = Validator()
        
        # Valores cercanos a E. coli K-12
        calculated = {
            "total_genes": 4651,
            "genome_length": 4641652,
            "gc_content": 50.79,
            "total_cds": 4318,
            "gene_density": 1002.0
        }
        
        result = validator.validate(calculated)
        assert result.overall_status.value == "PASS"


class TestNCBIDownloader:
    """Tests para el módulo de descarga de NCBI"""
    
    def test_downloader_import(self):
        """Test importación del módulo"""
        from app.core.ncbi_downloader import NCBIDownloader
        import tempfile
        
        with tempfile.TemporaryDirectory() as tmpdir:
            downloader = NCBIDownloader(tmpdir)
            assert downloader is not None


if __name__ == "__main__":
    pytest.main(["-v", __file__])
