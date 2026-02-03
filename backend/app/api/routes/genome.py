"""
Genome API Routes
Handles genome search, download, and management from NCBI Datasets API
"""
from fastapi import APIRouter, HTTPException, Request, BackgroundTasks
from fastapi.responses import StreamingResponse
from typing import Optional, List
from pydantic import BaseModel
import os
import json
import asyncio

from app.core.ncbi_downloader import NCBIDownloader, GenomeInfo


router = APIRouter()


# Pydantic models for request/response
class GenomeSearchResult(BaseModel):
    accession: str
    organism_name: str
    strain: str
    assembly_name: str
    assembly_level: str
    genome_size_mb: float
    is_reference: bool
    submission_date: str


class GenomeInfoResponse(BaseModel):
    accession: str
    organism_name: str
    organism_common_name: str
    tax_id: int
    assembly_name: str
    assembly_level: str
    genome_size: int
    genome_size_mb: float
    gene_count: int
    gc_percent: float
    submission_date: str
    refseq_category: str
    bioproject: str
    biosample: str
    strain: str
    is_reference: bool


class DownloadRequest(BaseModel):
    accession: str
    include_gbff: bool = True
    include_gff: bool = True
    include_fasta: bool = True
    include_protein: bool = False
    include_cds: bool = False
    include_rna: bool = False


class DownloadedGenome(BaseModel):
    accession: str
    path: str
    file_count: int
    files: List[dict]
    genome_info: Optional[dict] = None


# Store download status
download_status = {}


def get_downloader(request: Request) -> NCBIDownloader:
    """Get or create NCBIDownloader instance"""
    if not hasattr(request.app.state, 'ncbi_downloader'):
        project_root = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
        download_dir = os.path.join(project_root, "genomes")
        
        # Get API key from environment if available
        api_key = os.environ.get("NCBI_API_KEY")
        
        request.app.state.ncbi_downloader = NCBIDownloader(download_dir, api_key)
    
    return request.app.state.ncbi_downloader


@router.get("/search/{query}")
async def search_genomes(query: str, request: Request, limit: int = 10):
    """
    Search for genomes by organism name or taxon
    
    Examples:
    - /api/genome/search/Escherichia%20coli
    - /api/genome/search/562 (tax ID)
    - /api/genome/search/Bacillus%20subtilis
    """
    downloader = get_downloader(request)
    
    results = downloader.search_genomes(query, limit)
    
    if not results:
        raise HTTPException(status_code=404, detail=f"No se encontraron genomas para: {query}")
    
    return {
        "query": query,
        "count": len(results),
        "results": results
    }


@router.get("/info/{accession}", response_model=GenomeInfoResponse)
async def get_genome_info(accession: str, request: Request):
    """
    Get detailed information about a specific genome
    
    Example: /api/genome/info/GCF_000005845.2
    """
    downloader = get_downloader(request)
    
    info = downloader.get_genome_info(accession)
    
    if not info:
        raise HTTPException(status_code=404, detail=f"Genoma no encontrado: {accession}")
    
    return GenomeInfoResponse(
        accession=info.accession,
        organism_name=info.organism_name,
        organism_common_name=info.organism_common_name,
        tax_id=info.tax_id,
        assembly_name=info.assembly_name,
        assembly_level=info.assembly_level,
        genome_size=info.genome_size,
        genome_size_mb=round(info.genome_size / 1_000_000, 2),
        gene_count=info.gene_count,
        gc_percent=info.gc_percent,
        submission_date=info.submission_date,
        refseq_category=info.refseq_category,
        bioproject=info.bioproject,
        biosample=info.biosample,
        strain=info.strain,
        is_reference=info.is_reference
    )


@router.post("/download")
async def download_genome(download_req: DownloadRequest, request: Request, background_tasks: BackgroundTasks):
    """
    Download a genome from NCBI Datasets API
    
    This initiates a background download task.
    Use /api/genome/download-status/{accession} to check progress.
    """
    downloader = get_downloader(request)
    accession = download_req.accession
    
    # Check if already downloading
    if accession in download_status and download_status[accession].get("status") == "downloading":
        return {
            "message": "La descarga ya está en progreso",
            "accession": accession,
            "status": "downloading"
        }
    
    # Initialize status
    download_status[accession] = {
        "status": "starting",
        "message": "Iniciando descarga...",
        "accession": accession
    }
    
    def update_status(status: str, message: str):
        download_status[accession] = {
            "status": status,
            "message": message,
            "accession": accession
        }
    
    def download_task():
        result = downloader.download_genome(
            accession=accession,
            include_gbff=download_req.include_gbff,
            include_gff=download_req.include_gff,
            include_fasta=download_req.include_fasta,
            include_protein=download_req.include_protein,
            include_cds=download_req.include_cds,
            include_rna=download_req.include_rna,
            callback=update_status
        )
        
        if result["success"]:
            download_status[accession] = {
                "status": "completed",
                "message": f"Descarga completada: {len(result['files'])} archivos",
                "accession": accession,
                "result": result
            }
            
            # Save genome info
            if result.get("genome_info"):
                info_path = os.path.join(result["download_dir"], "genome_info.json")
                with open(info_path, 'w') as f:
                    json.dump(result["genome_info"], f, indent=2)
        else:
            download_status[accession] = {
                "status": "error",
                "message": result.get("error", "Error desconocido"),
                "accession": accession
            }
    
    # Add to background tasks
    background_tasks.add_task(download_task)
    
    return {
        "message": "Descarga iniciada",
        "accession": accession,
        "status": "starting"
    }


@router.get("/download-status/{accession}")
async def get_download_status(accession: str):
    """
    Get the status of a genome download
    """
    if accession not in download_status:
        return {
            "status": "not_found",
            "message": f"No hay descarga pendiente para {accession}",
            "accession": accession
        }
    
    return download_status[accession]


@router.get("/downloaded")
async def list_downloaded_genomes(request: Request):
    """
    List all downloaded genomes
    """
    downloader = get_downloader(request)
    
    genomes = downloader.get_available_genomes_dir()
    
    return {
        "count": len(genomes),
        "genomes": genomes
    }


@router.post("/activate/{accession}")
async def activate_genome(accession: str, request: Request):
    """
    Set a downloaded genome as the active genome for analysis
    
    This updates the file detector to use this genome's files
    """
    downloader = get_downloader(request)
    
    # Check if genome is downloaded
    genome_path = os.path.join(downloader.download_dir, accession)
    extracted_dir = os.path.join(genome_path, "extracted")
    
    if not os.path.exists(extracted_dir):
        raise HTTPException(status_code=404, detail=f"Genoma no descargado: {accession}")
    
    # Update the file detector to point to this genome
    file_detector = request.app.state.file_detector
    file_detector.scan_directory(extracted_dir)
    
    # Get genome info
    info = downloader.get_genome_info(accession)
    
    return {
        "message": f"Genoma {accession} activado",
        "accession": accession,
        "files_detected": len(file_detector.detected_files),
        "primary_file": file_detector.primary_file.filename if file_detector.primary_file else None,
        "genome_info": info.__dict__ if info else None
    }


@router.get("/current")
async def get_current_genome(request: Request):
    """
    Get information about the currently active genome
    """
    downloader = get_downloader(request)
    file_detector = request.app.state.file_detector
    
    if not file_detector.primary_file:
        return {
            "active": False,
            "message": "No hay genoma activo"
        }
    
    # Try to determine which genome is active
    current_path = file_detector.primary_file.filepath
    
    return {
        "active": True,
        "primary_file": file_detector.primary_file.filename,
        "primary_file_path": current_path,
        "files_count": len(file_detector.detected_files),
        "genome_info": downloader.current_genome.__dict__ if downloader.current_genome else None
    }


@router.delete("/{accession}")
async def delete_genome(accession: str, request: Request):
    """
    Delete a downloaded genome
    """
    downloader = get_downloader(request)
    
    success = downloader.delete_genome(accession)
    
    if not success:
        raise HTTPException(status_code=404, detail=f"Genoma no encontrado: {accession}")
    
    return {
        "message": f"Genoma {accession} eliminado",
        "accession": accession
    }


@router.get("/popular")
async def get_popular_genomes():
    """
    Get a list of popular/example genomes for quick access
    """
    popular = [
        {
            "accession": "GCF_000005845.2",
            "organism": "Escherichia coli str. K-12 substr. MG1655",
            "description": "Modelo clásico para estudios de genética bacteriana",
            "genome_size_mb": 4.6,
            "is_reference": True
        },
        {
            "accession": "GCF_000009045.1",
            "organism": "Bacillus subtilis subsp. subtilis str. 168",
            "description": "Bacteria modelo Gram-positiva",
            "genome_size_mb": 4.2,
            "is_reference": True
        },
        {
            "accession": "GCF_000146045.2",
            "organism": "Saccharomyces cerevisiae S288C",
            "description": "Levadura modelo eucariota",
            "genome_size_mb": 12.1,
            "is_reference": True
        },
        {
            "accession": "GCF_000195955.2",
            "organism": "Mycobacterium tuberculosis H37Rv",
            "description": "Agente causante de la tuberculosis",
            "genome_size_mb": 4.4,
            "is_reference": True
        },
        {
            "accession": "GCF_000006765.1",
            "organism": "Pseudomonas aeruginosa PAO1",
            "description": "Patógeno oportunista importante",
            "genome_size_mb": 6.3,
            "is_reference": True
        },
        {
            "accession": "GCF_000196035.1",
            "organism": "Staphylococcus aureus subsp. aureus NCTC 8325",
            "description": "Bacteria patógena común",
            "genome_size_mb": 2.8,
            "is_reference": True
        }
    ]
    
    return {
        "count": len(popular),
        "genomes": popular
    }
