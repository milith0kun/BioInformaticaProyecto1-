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
    Detect and list all available genomic files
    """
    file_detector = request.app.state.file_detector
    
    if not file_detector.detected_files:
        # Try to scan ncbi_dataset folder
        project_root = file_detector.project_root
        ncbi_folder = os.path.join(project_root, "ncbi_dataset")
        
        if os.path.exists(ncbi_folder):
            file_detector.scan_directory(ncbi_folder)
    
    summary = file_detector.get_files_summary()
    
    files = [
        FileInfo(
            filename=f["filename"],
            filepath=f["filepath"],
            extension=f["extension"],
            size_bytes=f["size_bytes"],
            file_type=f["file_type"],
            is_primary=f["is_primary"]
        )
        for f in summary["files"]
    ]
    
    return FileListResponse(
        files=files,
        total_count=summary["total_count"],
        primary_file=summary["primary_file"]
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
