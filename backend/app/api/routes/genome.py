"""
Genome API Routes
Handles genome search, download, and management from NCBI Datasets API
"""
from fastapi import APIRouter, HTTPException, Request, BackgroundTasks, Query
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
        # Usar el project_root del app.state que se configura en main.py
        if hasattr(request.app.state, 'project_root'):
            project_root = request.app.state.project_root
        else:
            # Fallback: calcular desde la ubicación actual
            project_root = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))))
            print(f"⚠️ [DOWNLOADER] project_root no encontrado en app.state, usando fallback")
        
        download_dir = os.path.join(project_root, "genomes")
        
        print(f"🔧 [DOWNLOADER] Inicializando con project_root: {project_root}")
        print(f"🔧 [DOWNLOADER] download_dir: {download_dir}")
        
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
    
    # Always return 200, even if no results
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
        print(f"🔽 [DOWNLOAD] Iniciando descarga de {accession}")
        try:
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
            
            print(f"🔽 [DOWNLOAD] Resultado para {accession}: success={result.get('success')}")
            print(f"🔽 [DOWNLOAD] download_dir: {result.get('download_dir')}")
            print(f"🔽 [DOWNLOAD] files: {result.get('files', [])}")
            
            if result["success"]:
                download_status[accession] = {
                    "status": "completed",
                    "message": f"Descarga completada: {len(result['files'])} archivos",
                    "accession": accession,
                    "result": result
                }
                print(f"✅ [DOWNLOAD] {accession} completado exitosamente")
                
                # Save genome info
                if result.get("genome_info"):
                    info_path = os.path.join(result["download_dir"], "genome_info.json")
                    with open(info_path, 'w') as f:
                        json.dump(result["genome_info"], f, indent=2)
            else:
                error_msg = result.get("error", "Error desconocido")
                print(f"❌ [DOWNLOAD] {accession} falló: {error_msg}")
                download_status[accession] = {
                    "status": "error",
                    "message": error_msg,
                    "accession": accession
                }
        except Exception as e:
            print(f"❌ [DOWNLOAD] Error en descarga de {accession}: {str(e)}")
            import traceback
            print(f"❌ [DOWNLOAD] Traceback: {traceback.format_exc()}")
            download_status[accession] = {
                "status": "error",
                "message": f"Error en descarga: {str(e)}",
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


# ================== COMPARACIÓN GENÓMICA ==================

from app.core.analyzers.genome_comparator import (
    get_genome_comparator, 
    RELATED_ECOLI_STRAINS,
    GENE_FUNCTIONAL_GROUPS
)


@router.get("/related-strains")
async def get_related_strains(category: Optional[str] = None):
    """
    Obtener cepas de E. coli emparentadas para comparación
    
    Args:
        category: Filtrar por categoría (laboratory, pathogenic, industrial, probiotic)
    """
    comparator = get_genome_comparator()
    strains = comparator.get_related_strains(category)
    
    return {
        "count": len(strains),
        "strains": strains,
        "categories": ["laboratory", "pathogenic", "industrial", "probiotic"]
    }


@router.get("/functional-groups")
async def get_functional_groups():
    """
    Obtener grupos funcionales disponibles para filtrado de genes
    """
    groups = []
    for group_id, info in GENE_FUNCTIONAL_GROUPS.items():
        groups.append({
            "id": group_id,
            "name": info["name"],
            "description": info["description"],
            "keywords": info["keywords"]
        })
    
    return {
        "count": len(groups),
        "groups": groups
    }


@router.post("/compare")
async def compare_genomes(
    request: Request, 
    accessions: Optional[str] = Query(None, description="Comma-separated list of accession numbers")
):
    """
    Comparar múltiples genomas descargados
    
    Si no se especifican accessions, compara todos los genomas descargados
    """
    print(f"🔍 [COMPARE] Iniciando comparación con accessions: {accessions}")
    project_root = request.app.state.project_root
    print(f"🔍 [COMPARE] project_root: {project_root}")
    
    # Convertir string separado por comas a lista
    accession_list = None
    if accessions:
        accession_list = [acc.strip() for acc in accessions.split(',')]
        print(f"🔍 [COMPARE] accession_list: {accession_list}")
    
    comparator = get_genome_comparator()
    
    # Obtener directorios de genomas
    all_genome_dirs = comparator.scan_downloaded_genomes(project_root)
    print(f"🔍 [COMPARE] all_genome_dirs encontrados: {len(all_genome_dirs)}")
    print(f"🔍 [COMPARE] Directorios: {all_genome_dirs}")
    
    if not all_genome_dirs:
        print("❌ [COMPARE] No hay genomas descargados")
        raise HTTPException(
            status_code=404,
            detail="No hay genomas descargados. Primero debes descargar los genomas seleccionados."
        )
    
    # Filtrar por accessions si se especifican
    if accession_list:
        print(f"🔍 [COMPARE] Filtrando genomas por accession_list: {accession_list}")
        genome_dirs = [
            d for d in all_genome_dirs 
            if os.path.basename(d) in accession_list
        ]
        print(f"🔍 [COMPARE] genome_dirs después del filtro: {genome_dirs}")
        print(f"🔍 [COMPARE] Basenames: {[os.path.basename(d) for d in genome_dirs]}")
        
        # Verificar cuáles genomas faltan
        found_accessions = {os.path.basename(d) for d in genome_dirs}
        missing_accessions = set(accession_list) - found_accessions
        
        if missing_accessions:
            raise HTTPException(
                status_code=404,
                detail=f"Los siguientes genomas no están descargados: {', '.join(missing_accessions)}. Descárgalos primero antes de comparar."
            )
    else:
        genome_dirs = all_genome_dirs
    
    if len(genome_dirs) < 1:
        raise HTTPException(
            status_code=400,
            detail="Se necesita al menos 1 genoma para analizar"
        )
    
    # Realizar comparación
    result = comparator.compare_genomes(genome_dirs)
    
    if not result:
        raise HTTPException(
            status_code=500,
            detail="Error al comparar genomas"
        )
    
    # Store genome data for dynamic validation
    from app.api.routes.analysis import _analysis_cache
    _analysis_cache["compared_genomes"] = [
        {
            "accession": g.accession,
            "organism_name": g.organism_name,
            "genome_length": g.genome_length,
            "total_genes": g.gene_count,
            "total_cds": g.cds_count,
            "gc_content": g.gc_content,
            "gene_density": g.gene_density
        }
        for g in result.genomes
    ]
    
    # Convertir a diccionario serializable
    return {
        "comparison_date": result.comparison_date,
        "total_genomes_compared": result.total_genomes_compared,
        "genomes": [
            {
                "accession": g.accession,
                "organism_name": g.organism_name,
                "genome_length": g.genome_length,
                "gene_count": g.gene_count,
                "cds_count": g.cds_count,
                "gc_content": g.gc_content,
                "gene_density": g.gene_density,
                "avg_gene_length": g.avg_gene_length,
                "min_gene_length": g.min_gene_length,
                "max_gene_length": g.max_gene_length
            }
            for g in result.genomes
        ],
        "extremes": {
            "largest_genome": result.largest_genome,
            "smallest_genome": result.smallest_genome,
            "highest_gc": result.highest_gc,
            "lowest_gc": result.lowest_gc,
            "highest_gene_density": result.highest_gene_density,
            "lowest_gene_density": result.lowest_gene_density
        },
        "longest_genes_global": result.longest_genes_global,
        "shortest_genes_global": result.shortest_genes_global,
        "metrics_comparison": result.metrics_comparison
    }


@router.get("/genes/by-size/{order}")
async def get_genes_by_size(
    request: Request,
    order: str = "largest",
    count: int = 10
):
    """
    Obtener genes ordenados por tamaño (mayor o menor)
    
    Args:
        order: 'largest' o 'smallest'
        count: Número de genes a retornar (default 10)
    """
    from app.api.routes.analysis import _analysis_cache
    
    if _analysis_cache["genome_data"] is None:
        raise HTTPException(
            status_code=404,
            detail="No hay análisis ejecutado. Ejecute el análisis primero."
        )
    
    genome_data = _analysis_cache["genome_data"]
    comparator = get_genome_comparator()
    
    result = comparator.get_genes_by_size(
        genes=genome_data.genes,
        order=order,
        count=count,
        genome_accession=genome_data.accession
    )
    
    return {
        "genes": result.genes,
        "total_count": result.total_count,
        "filter_applied": result.filter_applied,
        "genome_source": result.genome_source
    }


@router.get("/genes/by-group/{group_id}")
async def get_genes_by_group(request: Request, group_id: str):
    """
    Obtener genes filtrados por grupo funcional
    
    Args:
        group_id: ID del grupo (metabolism, transport, regulation, etc.)
    """
    from app.api.routes.analysis import _analysis_cache
    
    if _analysis_cache["genome_data"] is None:
        raise HTTPException(
            status_code=404,
            detail="No hay análisis ejecutado. Ejecute el análisis primero."
        )
    
    genome_data = _analysis_cache["genome_data"]
    comparator = get_genome_comparator()
    
    result = comparator.filter_genes_by_group(
        genes=genome_data.genes,
        group_id=group_id,
        genome_accession=genome_data.accession
    )
    
    return {
        "genes": result.genes,
        "total_count": result.total_count,
        "filter_applied": result.filter_applied,
        "genome_source": result.genome_source,
        "group_info": GENE_FUNCTIONAL_GROUPS.get(group_id, {})
    }


@router.get("/genes/search-advanced")
async def search_genes_advanced(
    request: Request,
    query: Optional[str] = None,
    min_length: Optional[int] = None,
    max_length: Optional[int] = None,
    min_gc: Optional[float] = None,
    max_gc: Optional[float] = None,
    strand: Optional[str] = None,
    page: int = 1,
    page_size: int = 50
):
    """
    Búsqueda avanzada de genes con múltiples filtros
    """
    from app.api.routes.analysis import _analysis_cache
    
    if _analysis_cache["genome_data"] is None:
        raise HTTPException(
            status_code=404,
            detail="No hay análisis ejecutado. Ejecute el análisis primero."
        )
    
    genome_data = _analysis_cache["genome_data"]
    comparator = get_genome_comparator()
    
    result = comparator.search_genes_advanced(
        genes=genome_data.genes,
        query=query or "",
        min_length=min_length,
        max_length=max_length,
        min_gc=min_gc,
        max_gc=max_gc,
        strand=strand,
        genome_accession=genome_data.accession
    )
    
    # Aplicar paginación
    total = result.total_count
    total_pages = (total + page_size - 1) // page_size
    start = (page - 1) * page_size
    end = start + page_size
    
    return {
        "genes": result.genes[start:end],
        "total_count": total,
        "page": page,
        "page_size": page_size,
        "total_pages": total_pages,
        "filter_applied": result.filter_applied,
        "genome_source": result.genome_source
    }


@router.get("/genes/groups-summary")
async def get_gene_groups_summary(request: Request):
    """
    Obtener resumen de genes por grupo funcional
    """
    from app.api.routes.analysis import _analysis_cache
    
    if _analysis_cache["genome_data"] is None:
        raise HTTPException(
            status_code=404,
            detail="No hay análisis ejecutado. Ejecute el análisis primero."
        )
    
    genome_data = _analysis_cache["genome_data"]
    comparator = get_genome_comparator()
    
    summaries = comparator.get_gene_groups_summary(
        genes=genome_data.genes,
        genome_accession=genome_data.accession
    )
    
    return {
        "genome": genome_data.accession,
        "total_genes": len(genome_data.genes),
        "groups": summaries
    }
