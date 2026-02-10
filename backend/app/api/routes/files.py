"""
Files API Routes
Handles file detection, extraction, and information endpoints
"""
from fastapi import APIRouter, HTTPException, Request
from typing import Optional
import os

from app.models.schemas import FileListResponse, FileInfo

router = APIRouter()


@router.get("/detect", response_model=FileListResponse)
async def detect_files(request: Request):
    """
    Detect and list all available genomic files from ALL downloaded genomes
    """
    file_detector = request.app.state.file_detector
    project_root = file_detector.project_root
    
    # Buscar en la carpeta genomes/ todos los genomas descargados
    genomes_path = os.path.join(project_root, "genomes")
    all_files = []
    
    if os.path.exists(genomes_path):
        for accession_dir in os.listdir(genomes_path):
            if not (accession_dir.startswith("GCF_") or accession_dir.startswith("GCA_")):
                continue
            
            # Buscar archivos en extracted/ncbi_dataset/data/
            extracted_path = os.path.join(
                genomes_path, 
                accession_dir, 
                "extracted", 
                "ncbi_dataset", 
                "data",
                accession_dir
            )
            
            if os.path.exists(extracted_path):
                for filename in os.listdir(extracted_path):
                    filepath = os.path.join(extracted_path, filename)
                    if os.path.isfile(filepath):
                        file_size = os.path.getsize(filepath)
                        extension = os.path.splitext(filename)[1]
                        
                        # Determinar tipo de archivo
                        file_type = "Unknown"
                        if extension == ".gbff":
                            file_type = "GenBank Full Flat File"
                        elif extension == ".fna":
                            file_type = "FASTA Nucleotide"
                        elif extension == ".gff":
                            file_type = "GFF3 Annotation"
                        elif extension == ".gtf":
                            file_type = "GTF Annotation"
                        elif extension == ".faa":
                            file_type = "FASTA Protein"
                        elif extension == ".jsonl":
                            file_type = "JSON Lines Report"
                        elif extension == ".txt":
                            file_type = "Text File"
                        
                        all_files.append({
                            "filename": filename,
                            "filepath": filepath,
                            "extension": extension,
                            "size_bytes": file_size,
                            "size_mb": round(file_size / (1024 * 1024), 2),
                            "file_type": file_type,
                            "is_primary": extension == ".gbff",
                            "accession": accession_dir
                        })
    
    # También buscar en ncbi_dataset/ (carpeta antigua)
    ncbi_folder = os.path.join(project_root, "ncbi_dataset")
    if os.path.exists(ncbi_folder) and not file_detector.detected_files:
        file_detector.scan_directory(ncbi_folder)
        summary = file_detector.get_files_summary()
        for f in summary["files"]:
            all_files.append({
                **f,
                "accession": "ncbi_dataset"
            })
    
    files = [
        FileInfo(
            filename=f["filename"],
            filepath=f["filepath"],
            extension=f["extension"],
            size_bytes=f["size_bytes"],
            size_mb=round(f["size_bytes"] / (1024 * 1024), 2),
            file_type=f["file_type"],
            is_primary=f.get("is_primary", False),
            accession=f.get("accession", "")
        )
        for f in all_files
    ]
    
    primary_file = next((f["filepath"] for f in all_files if f.get("is_primary")), None)
    
    return FileListResponse(
        files=files,
        total_count=len(all_files),
        primary_file=primary_file
    )


@router.post("/extract")
async def extract_zip(request: Request):
    """
    Extract ncbi_dataset.zip if not already extracted
    """
    file_detector = request.app.state.file_detector
    project_root = file_detector.project_root
    
    zip_path = os.path.join(project_root, "ncbi_dataset.zip")
    
    if not os.path.exists(zip_path):
        # Check if already extracted
        ncbi_folder = os.path.join(project_root, "ncbi_dataset")
        if os.path.exists(ncbi_folder):
            file_detector.scan_directory(ncbi_folder)
            return {
                "message": "Los datos ya están extraídos",
                "extracted": True,
                "file_count": len(file_detector.detected_files)
            }
        raise HTTPException(status_code=404, detail="ncbi_dataset.zip no encontrado")
    
    success = file_detector.extract_zip(zip_path)
    
    if not success:
        raise HTTPException(status_code=500, detail="Error al extraer el archivo ZIP")
    
    return {
        "message": "Archivo ZIP extraído correctamente",
        "extracted": True,
        "file_count": len(file_detector.detected_files)
    }


@router.get("/{filename}/info")
async def get_file_info(filename: str, request: Request):
    """
    Get detailed information about a specific file
    """
    file_detector = request.app.state.file_detector
    
    for file in file_detector.detected_files:
        if file.filename == filename:
            return {
                "filename": file.filename,
                "filepath": file.filepath,
                "extension": file.extension,
                "size_bytes": file.size_bytes,
                "size_mb": round(file.size_bytes / (1024 * 1024), 2),
                "file_type": file.file_type,
                "is_primary": file.is_primary,
                "hash": file_detector.get_file_hash(file.filepath)
            }
    
    raise HTTPException(status_code=404, detail=f"Archivo {filename} no encontrado")


@router.get("/primary")
async def get_primary_file(request: Request):
    """
    Get the primary GenBank file for analysis
    """
    file_detector = request.app.state.file_detector
    
    if not file_detector.primary_file:
        raise HTTPException(
            status_code=404, 
            detail="No se encontró archivo principal. Ejecute /detect primero."
        )
    
    pf = file_detector.primary_file
    return {
        "filename": pf.filename,
        "filepath": pf.filepath,
        "extension": pf.extension,
        "size_bytes": pf.size_bytes,
        "file_type": pf.file_type
    }
