"""
File Detector Module
Handles detection, extraction and cataloging of genomic files
"""
import os
import zipfile
from typing import List, Dict, Optional
from dataclasses import dataclass, field
import hashlib


@dataclass
class DetectedFile:
    filename: str
    filepath: str
    extension: str
    size_bytes: int
    file_type: str
    is_primary: bool = False
    file_hash: Optional[str] = None


FILE_TYPE_MAP = {
    '.gbff': 'GenBank Full Flat File',
    '.gb': 'GenBank',
    '.genbank': 'GenBank',
    '.fna': 'FASTA Nucleotide',
    '.fasta': 'FASTA',
    '.fa': 'FASTA',
    '.faa': 'FASTA Amino Acid',
    '.gff': 'General Feature Format',
    '.gff3': 'General Feature Format v3',
    '.gtf': 'Gene Transfer Format',
    '.jsonl': 'JSON Lines',
    '.json': 'JSON',
    '.tsv': 'Tab-Separated Values',
    '.csv': 'Comma-Separated Values',
    '.txt': 'Text File',
}

# Priority order for primary file selection
PRIMARY_FILE_PRIORITY = ['.gbff', '.gb', '.genbank', '.gff', '.gtf', '.fna', '.fasta']


class FileDetector:
    """Detects and catalogs genomic files in a project directory"""
    
    def __init__(self, project_root: str):
        self.project_root = project_root
        self.extracted_dir = os.path.join(project_root, "backend", "extracted")
        self.detected_files: List[DetectedFile] = []
        self.primary_file: Optional[DetectedFile] = None
        
    def extract_zip(self, zip_path: str) -> bool:
        """Extract ZIP file to extracted directory"""
        try:
            os.makedirs(self.extracted_dir, exist_ok=True)
            
            with zipfile.ZipFile(zip_path, 'r') as zip_ref:
                zip_ref.extractall(self.extracted_dir)
            
            # Scan the extracted directory
            self.scan_directory(self.extracted_dir)
            return True
        except Exception as e:
            print(f"Error extracting ZIP: {e}")
            return False
    
    def scan_directory(self, directory: str) -> List[DetectedFile]:
        """Recursively scan directory for genomic files with performance tracking"""
        import time
        start_time = time.time()
        print(f"🔍 [SCAN] Iniciando escaneo en: {directory}")
        
        self.detected_files = []
        
        # EXCLUDE heavy directories to prevent timeouts
        EXCLUDE_DIRS = {'.git', '.venv', 'node_modules', '__pycache__', 'dist', 'venv'}
        
        try:
            for root, dirs, files in os.walk(directory):
                # Prune excluded directories in-place to prevent os.walk from entering them
                dirs[:] = [d for d in dirs if d not in EXCLUDE_DIRS]
                
                for filename in files:
                    filepath = os.path.join(root, filename)
                    ext = self._get_extension(filename)
                    
                    if ext in FILE_TYPE_MAP:
                        try:
                            size = os.path.getsize(filepath)
                            file_type = FILE_TYPE_MAP.get(ext, 'Unknown')
                            
                            detected = DetectedFile(
                                filename=filename,
                                filepath=filepath,
                                extension=ext,
                                size_bytes=size,
                                file_type=file_type
                            )
                            self.detected_files.append(detected)
                        except Exception as e:
                            print(f"⚠️ [SCAN] Error en archivo {filename}: {e}")
        except Exception as e:
            print(f"🔥 [SCAN] Error crítico durante os.walk: {e}")
        
        # Determine primary file
        self._select_primary_file()
        
        duration = time.time() - start_time
        print(f"✅ [SCAN] Escaneo completado en {duration:.2f}s. Encontrados: {len(self.detected_files)} archivos.")
        
        return self.detected_files
    
    def _get_extension(self, filename: str) -> str:
        """Get file extension, handling compound extensions"""
        lower = filename.lower()
        
        # Check for compound extensions
        if lower.endswith('.gbff'):
            return '.gbff'
        if lower.endswith('.gff3'):
            return '.gff3'
            
        # Standard extension
        _, ext = os.path.splitext(lower)
        return ext
    
    def _select_primary_file(self) -> None:
        """Select the primary file based on priority"""
        for ext in PRIMARY_FILE_PRIORITY:
            for file in self.detected_files:
                if file.extension == ext:
                    # Prefer GCF (RefSeq) over GCA (GenBank) if available
                    if 'GCF_' in file.filepath:
                        file.is_primary = True
                        self.primary_file = file
                        return
            
            # If no GCF found, use first match
            for file in self.detected_files:
                if file.extension == ext:
                    file.is_primary = True
                    self.primary_file = file
                    return
    
    def get_file_by_extension(self, extension: str) -> Optional[DetectedFile]:
        """Get first file with given extension"""
        for file in self.detected_files:
            if file.extension == extension:
                return file
        return None
    
    def get_genbank_file(self) -> Optional[DetectedFile]:
        """Get the GenBank file (.gbff)"""
        return self.get_file_by_extension('.gbff')
    
    def get_file_hash(self, filepath: str) -> str:
        """Calculate MD5 hash of file for caching"""
        hash_md5 = hashlib.md5()
        with open(filepath, "rb") as f:
            for chunk in iter(lambda: f.read(4096), b""):
                hash_md5.update(chunk)
        return hash_md5.hexdigest()
    
    def get_files_summary(self) -> Dict:
        """Get summary of detected files"""
        return {
            "files": [
                {
                    "filename": f.filename,
                    "filepath": f.filepath,
                    "extension": f.extension,
                    "size_bytes": f.size_bytes,
                    "file_type": f.file_type,
                    "is_primary": f.is_primary
                }
                for f in self.detected_files
            ],
            "total_count": len(self.detected_files),
            "primary_file": self.primary_file.filename if self.primary_file else None
        }
