"""
NCBI Datasets API Downloader Module
Handles automatic genome download from NCBI Datasets REST API v2
"""
import os
import zipfile
import requests
import json
import shutil
import time
from typing import Dict, Optional, List
from dataclasses import dataclass, field
from datetime import datetime
import hashlib


@dataclass
class GenomeInfo:
    """Information about a genome from NCBI"""
    accession: str
    organism_name: str = ""
    organism_common_name: str = ""
    tax_id: int = 0
    assembly_name: str = ""
    assembly_level: str = ""
    genome_size: int = 0
    gene_count: int = 0
    gc_percent: float = 0.0
    submission_date: str = ""
    refseq_category: str = ""
    bioproject: str = ""
    biosample: str = ""
    strain: str = ""
    is_reference: bool = False
    annotation_info: Dict = field(default_factory=dict)


class NCBIDownloader:
    """
    Downloads genome data packages from NCBI Datasets API v2
    
    API Documentation: https://www.ncbi.nlm.nih.gov/datasets/docs/v2/api/rest-api/
    """
    
    BASE_URL = "https://api.ncbi.nlm.nih.gov/datasets/v2"
    
    # Available annotation types for download
    ANNOTATION_TYPES = [
        "GENOME_FASTA",      # Genomic sequence (FASTA)
        "GENOME_GFF",        # Annotation (GFF3)
        "GENOME_GBFF",       # GenBank format (full annotations)
        "RNA_FASTA",         # RNA sequences
        "CDS_FASTA",         # Coding sequences
        "PROT_FASTA",        # Protein sequences
        "SEQUENCE_REPORT",   # Sequence report (JSONL)
    ]
    
    def __init__(self, download_dir: str, api_key: Optional[str] = None):
        """
        Initialize the downloader
        
        Args:
            download_dir: Directory to save downloaded files
            api_key: Optional NCBI API key for higher rate limits
        """
        self.download_dir = download_dir
        self.api_key = api_key
        self.current_genome: Optional[GenomeInfo] = None
        self.session = requests.Session()
        
        # Set headers
        self.session.headers.update({
            "Accept": "application/json",
            "User-Agent": "BioInformatica-Proyecto1/1.0"
        })
        
        if api_key:
            self.session.headers["api-key"] = api_key
        
        # Create download directory
        os.makedirs(download_dir, exist_ok=True)
    
    def get_genome_info(self, accession: str) -> Optional[GenomeInfo]:
        """
        Get metadata about a genome from NCBI
        
        Args:
            accession: Assembly accession (e.g., GCF_000005845.2)
            
        Returns:
            GenomeInfo object with genome metadata
        """
        url = f"{self.BASE_URL}/genome/accession/{accession}/dataset_report"
        
        try:
            response = self.session.get(url, timeout=30)
            response.raise_for_status()
            data = response.json()
            
            if "reports" not in data or len(data["reports"]) == 0:
                return None
            
            report = data["reports"][0]
            assembly_info = report.get("assembly_info", {})
            assembly_stats = report.get("assembly_stats", {})
            annotation_info = report.get("annotation_info", {})
            organism = report.get("organism", {})
            
            # Safe type conversion helper
            def safe_int(value, default=0):
                try:
                    return int(value) if value is not None else default
                except (ValueError, TypeError):
                    return default
            
            def safe_float(value, default=0.0):
                try:
                    return float(value) if value is not None else default
                except (ValueError, TypeError):
                    return default
            
            genome_info = GenomeInfo(
                accession=report.get("accession", accession),
                organism_name=organism.get("organism_name", "Unknown"),
                organism_common_name=organism.get("organism_common_name", ""),
                tax_id=safe_int(organism.get("tax_id")),
                assembly_name=assembly_info.get("assembly_name", ""),
                assembly_level=assembly_info.get("assembly_level", ""),
                genome_size=safe_int(assembly_stats.get("total_sequence_length")),
                gene_count=safe_int(annotation_info.get("stats", {}).get("gene_counts", {}).get("total")),
                gc_percent=safe_float(assembly_stats.get("gc_percent")),
                submission_date=assembly_info.get("submission_date", ""),
                refseq_category=assembly_info.get("refseq_category", ""),
                bioproject=assembly_info.get("bioproject_accession", ""),
                biosample=assembly_info.get("biosample", {}).get("accession", "") if isinstance(assembly_info.get("biosample"), dict) else "",
                strain=organism.get("infraspecific_names", {}).get("strain", "") if isinstance(organism.get("infraspecific_names"), dict) else "",
                is_reference=assembly_info.get("refseq_category", "") == "reference genome",
                annotation_info=annotation_info
            )
            
            self.current_genome = genome_info
            return genome_info
            
        except requests.exceptions.RequestException as e:
            print(f"Error fetching genome info: {e}")
            return None
        except (KeyError, json.JSONDecodeError) as e:
            print(f"Error parsing genome info: {e}")
            return None
    
    def search_genomes(self, query: str, limit: int = 10) -> List[Dict]:
        """
        Search for genomes by organism name or taxon
        
        Args:
            query: Search term (organism name, taxon ID, etc.)
            limit: Maximum number of results
            
        Returns:
            List of matching genome summaries
        """
        # Safe type conversion helpers
        def safe_int(value, default=0):
            try:
                return int(value) if value is not None else default
            except (ValueError, TypeError):
                return default
        
        def safe_float(value, default=0.0):
            try:
                return float(value) if value is not None else default
            except (ValueError, TypeError):
                return default
        
        # If query is very short, try common organisms
        search_queries = [query]
        
        # Add common expansions for short queries
        if len(query) <= 10:
            common_expansions = {
                "coli": ["Escherichia coli", "562"],  # 562 is E. coli tax ID
                "salmonella": ["Salmonella enterica", "28901"],
                "bacillus": ["Bacillus subtilis", "1423"],
                "staph": ["Staphylococcus aureus", "1280"],
                "pseudo": ["Pseudomonas aeruginosa", "287"],
                "tuberculosis": ["Mycobacterium tuberculosis", "1773"],
                "strepto": ["Streptococcus pneumoniae", "1313"],
            }
            
            query_lower = query.lower()
            for key, expansions in common_expansions.items():
                if key in query_lower:
                    search_queries.extend(expansions)
                    break
        
        # Try each query until we get results
        all_results = []
        for search_query in search_queries:
            # Try to search by taxon name or ID
            url = f"{self.BASE_URL}/genome/taxon/{search_query}/dataset_report"
            params = {
                "page_size": limit,
                "filters.reference_only": "false"
            }
            
            try:
                response = self.session.get(url, params=params, timeout=30)
                response.raise_for_status()
                data = response.json()
                
                for report in data.get("reports", [])[:limit]:
                    assembly_info = report.get("assembly_info", {})
                    organism = report.get("organism", {})
                    assembly_stats = report.get("assembly_stats", {})
                    
                    # Safe extraction of strain
                    infraspecific = organism.get("infraspecific_names")
                    strain = ""
                    if isinstance(infraspecific, dict):
                        strain = infraspecific.get("strain", "")
                    
                    genome_size = safe_int(assembly_stats.get("total_sequence_length"))
                    genome_size_mb = round(genome_size / 1_000_000, 2) if genome_size > 0 else 0
                    
                    all_results.append({
                        "accession": report.get("accession", ""),
                        "organism_name": organism.get("organism_name", "Unknown"),
                        "strain": strain,
                        "assembly_name": assembly_info.get("assembly_name", ""),
                        "assembly_level": assembly_info.get("assembly_level", ""),
                        "genome_size_mb": genome_size_mb,
                        "is_reference": assembly_info.get("refseq_category", "") == "reference genome",
                        "submission_date": assembly_info.get("submission_date", "")
                    })
                
                # If we got results, stop trying other queries
                if all_results:
                    break
                    
            except requests.exceptions.RequestException as e:
                print(f"Error searching genomes with query '{search_query}': {e}")
                continue
        
        return all_results[:limit]
    
    def download_genome(
        self, 
        accession: str, 
        include_gbff: bool = True,
        include_gff: bool = True,
        include_fasta: bool = True,
        include_protein: bool = False,
        include_cds: bool = False,
        include_rna: bool = False,
        callback=None
    ) -> Dict:
        """
        Download genome data package from NCBI
        
        Args:
            accession: Assembly accession (e.g., GCF_000005845.2)
            include_gbff: Include GenBank format file
            include_gff: Include GFF3 annotation
            include_fasta: Include genomic FASTA
            include_protein: Include protein sequences
            include_cds: Include CDS sequences
            include_rna: Include RNA sequences
            callback: Optional callback function for progress updates
            
        Returns:
            Dictionary with download result and file paths
        """
        # Build annotation types list
        annotation_types = []
        if include_fasta:
            annotation_types.append("GENOME_FASTA")
        if include_gbff:
            annotation_types.append("GENOME_GBFF")
        if include_gff:
            annotation_types.append("GENOME_GFF")
        if include_protein:
            annotation_types.append("PROT_FASTA")
        if include_cds:
            annotation_types.append("CDS_FASTA")
        if include_rna:
            annotation_types.append("RNA_FASTA")
        
        # Always include sequence report
        annotation_types.append("SEQUENCE_REPORT")
        
        # Build download URL
        url = f"{self.BASE_URL}/genome/accession/{accession}/download"
        params = {
            "include_annotation_type": annotation_types,
            "hydrated": "FULLY_HYDRATED"
        }
        
        # Create accession-specific directory
        genome_dir = os.path.join(self.download_dir, accession)
        os.makedirs(genome_dir, exist_ok=True)
        
        zip_path = os.path.join(genome_dir, f"{accession}_dataset.zip")
        
        # Intentar hasta 3 veces si hay errores
        max_retries = 3
        last_error = None
        
        for attempt in range(max_retries):
            try:
                if attempt > 0:
                    if callback:
                        callback("downloading", f"Reintentando descarga (intento {attempt + 1}/{max_retries})...")
                    # Limpiar archivo corrupto si existe
                    if os.path.exists(zip_path):
                        os.remove(zip_path)
                
                if callback:
                    callback("downloading", f"Descargando {accession} desde NCBI...")
                
                # Download the ZIP file
                response = self.session.get(url, params=params, stream=True, timeout=300)
                response.raise_for_status()
                
                # Verificar que la respuesta es un ZIP válido
                content_type = response.headers.get('content-type', '')
                if 'html' in content_type.lower() or 'text' in content_type.lower():
                    raise Exception(f"NCBI devolvió {content_type} en lugar de un archivo ZIP. Posible error de API o accession inválido.")
                
                # Get total size if available
                total_size = int(response.headers.get('content-length', 0))
                
                # Verificar tamaño mínimo (un ZIP válido tiene al menos algunos KB)
                if total_size > 0 and total_size < 1000:
                    raise Exception(f"Archivo muy pequeño ({total_size} bytes). Probablemente no sea un genoma válido.")
                
                downloaded = 0
                
                with open(zip_path, 'wb') as f:
                    for chunk in response.iter_content(chunk_size=8192):
                        if chunk:
                            f.write(chunk)
                            downloaded += len(chunk)
                            if callback and total_size > 0:
                                progress = (downloaded / total_size) * 100
                                callback("downloading", f"Descargando... {progress:.1f}%")
                
                # Verificar que el archivo descargado es un ZIP válido
                if not zipfile.is_zipfile(zip_path):
                    file_size = os.path.getsize(zip_path)
                    raise Exception(f"El archivo descargado ({file_size} bytes) no es un ZIP válido. Reintentando...")
                
                # Si llegamos aquí, la descarga fue exitosa
                break
                
            except Exception as e:
                last_error = e
                if attempt < max_retries - 1:
                    print(f"⚠️ Error en descarga de {accession} (intento {attempt + 1}): {str(e)}")
                    time.sleep(2)  # Esperar 2 segundos antes de reintentar
                else:
                    # Último intento falló
                    raise Exception(f"No se pudo descargar {accession} después de {max_retries} intentos. Error: {str(last_error)}")
        
        try:
            if callback:
                callback("extracting", "Extrayendo archivos...")
            
            # Extract the ZIP
            extracted_dir = os.path.join(genome_dir, "extracted")
            if os.path.exists(extracted_dir):
                shutil.rmtree(extracted_dir)
            
            with zipfile.ZipFile(zip_path, 'r') as zip_ref:
                zip_ref.extractall(extracted_dir)
            
            # Find extracted files
            files = self._scan_extracted_files(extracted_dir)
            
            # Get genome info if not already loaded
            if not self.current_genome or self.current_genome.accession != accession:
                self.get_genome_info(accession)
            
            # Create download report
            report = self._create_download_report(accession, genome_dir, files)
            
            if callback:
                callback("completed", f"Descarga completada: {len(files)} archivos")
            
            return {
                "success": True,
                "accession": accession,
                "download_dir": genome_dir,
                "extracted_dir": extracted_dir,
                "zip_path": zip_path,
                "files": files,
                "genome_info": self.current_genome.__dict__ if self.current_genome else None,
                "report": report
            }
            
        except requests.exceptions.RequestException as e:
            error_msg = f"Error downloading genome: {e}"
            print(error_msg)
            if callback:
                callback("error", error_msg)
            return {
                "success": False,
                "error": str(e),
                "accession": accession
            }
        except zipfile.BadZipFile as e:
            error_msg = f"Invalid ZIP file: {e}"
            print(error_msg)
            if callback:
                callback("error", error_msg)
            return {
                "success": False,
                "error": error_msg,
                "accession": accession
            }
    
    def _scan_extracted_files(self, directory: str) -> List[Dict]:
        """Scan extracted directory for genomic files"""
        files = []
        
        for root, dirs, filenames in os.walk(directory):
            for filename in filenames:
                filepath = os.path.join(root, filename)
                ext = os.path.splitext(filename.lower())[1]
                
                # Handle compound extensions
                if filename.lower().endswith('.gbff'):
                    ext = '.gbff'
                elif filename.lower().endswith('.gff'):
                    ext = '.gff'
                elif filename.lower().endswith('.fna'):
                    ext = '.fna'
                elif filename.lower().endswith('.faa'):
                    ext = '.faa'
                
                file_type = self._get_file_type(ext)
                
                files.append({
                    "filename": filename,
                    "filepath": filepath,
                    "extension": ext,
                    "size_bytes": os.path.getsize(filepath),
                    "size_mb": round(os.path.getsize(filepath) / (1024 * 1024), 2),
                    "file_type": file_type
                })
        
        return files
    
    def _get_file_type(self, extension: str) -> str:
        """Get human-readable file type from extension"""
        type_map = {
            '.gbff': 'GenBank Full Flat File',
            '.gff': 'GFF3 Annotation',
            '.fna': 'FASTA Nucleotide',
            '.faa': 'FASTA Amino Acid',
            '.fasta': 'FASTA',
            '.jsonl': 'JSON Lines Report',
            '.json': 'JSON',
            '.txt': 'Text File',
            '.md': 'Markdown'
        }
        return type_map.get(extension, 'Unknown')
    
    def _create_download_report(self, accession: str, genome_dir: str, files: List[Dict]) -> str:
        """Create a download report file"""
        report_path = os.path.join(genome_dir, "download_report.txt")
        
        report_content = f"""======================================================================
REPORTE DE DESCARGA - NCBI Datasets API v2
======================================================================

Fecha: {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}
Accession: {accession}
"""
        
        if self.current_genome:
            report_content += f"""Organismo: {self.current_genome.organism_name}
Cepa: {self.current_genome.strain}
Tamaño del genoma: {self.current_genome.genome_size:,} bp
Contenido GC: {self.current_genome.gc_percent}%
Genes anotados: {self.current_genome.gene_count:,}
Nivel de ensamblaje: {self.current_genome.assembly_level}
Categoría RefSeq: {self.current_genome.refseq_category}
BioProject: {self.current_genome.bioproject}
"""
        
        report_content += f"""
Archivos descargados: {len(files)}
"""
        
        for f in files:
            report_content += f"  - {f['filename']} ({f['size_mb']} MB) - {f['file_type']}\n"
        
        report_content += """
Estado: ✅ EXITOSO
======================================================================
"""
        
        with open(report_path, 'w', encoding='utf-8') as f:
            f.write(report_content)
        
        return report_path
    
    def get_available_genomes_dir(self) -> List[Dict]:
        """List all downloaded genomes"""
        genomes = []
        
        if not os.path.exists(self.download_dir):
            return genomes
        
        for entry in os.listdir(self.download_dir):
            genome_path = os.path.join(self.download_dir, entry)
            if os.path.isdir(genome_path):
                extracted_dir = os.path.join(genome_path, "extracted")
                if os.path.exists(extracted_dir):
                    files = self._scan_extracted_files(extracted_dir)
                    
                    # Check for genome info
                    info_file = os.path.join(genome_path, "genome_info.json")
                    genome_info = None
                    if os.path.exists(info_file):
                        with open(info_file, 'r') as f:
                            genome_info = json.load(f)
                    
                    genomes.append({
                        "accession": entry,
                        "path": genome_path,
                        "extracted_dir": extracted_dir,
                        "file_count": len(files),
                        "files": files,
                        "genome_info": genome_info
                    })
        
        return genomes
    
    def set_active_genome(self, accession: str) -> bool:
        """
        Set a downloaded genome as the active genome for analysis
        
        Returns True if successful, False otherwise
        """
        genome_path = os.path.join(self.download_dir, accession)
        extracted_dir = os.path.join(genome_path, "extracted")
        
        if not os.path.exists(extracted_dir):
            return False
        
        # Load genome info
        self.get_genome_info(accession)
        
        return True
    
    def delete_genome(self, accession: str) -> bool:
        """Delete a downloaded genome"""
        genome_path = os.path.join(self.download_dir, accession)
        
        if os.path.exists(genome_path):
            shutil.rmtree(genome_path)
            return True
        
        return False
