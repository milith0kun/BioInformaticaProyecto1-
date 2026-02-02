"""
Analysis API Routes
Handles genomic analysis endpoints
"""
from fastapi import APIRouter, HTTPException, Request, BackgroundTasks
from typing import Optional
import os
import json

from app.models.schemas import (
    CodonAnalysisResponse, 
    GeneAnalysisResponse,
    ValidationResponse,
    CompleteAnalysisResponse,
    AnalysisStatus,
    StopCodonInfo,
    GeneComparison,
    SizeStatistics,
    GeneInfo,
    ValidationItem,
    ValidationStatus,
    PaginatedGenesResponse
)
from app.core.genome_parser import GenomeParser
from app.core.analyzers.codon_analyzer import CodonAnalyzer
from app.core.analyzers.gene_analyzer import GeneAnalyzer
from app.core.analyzers.validator import Validator
from app.core.analyzers.ai_validator import get_ai_validator

router = APIRouter()

# Global state for analysis results
_analysis_cache = {
    "codons": None,
    "genes": None,
    "validation": None,
    "genome_data": None,
    "status": "idle"
}


def get_genbank_path(request: Request) -> str:
    """Get path to primary GenBank file"""
    file_detector = request.app.state.file_detector
    
    if not file_detector.detected_files:
        # Try to scan ncbi_dataset folder
        project_root = file_detector.project_root
        ncbi_folder = os.path.join(project_root, "ncbi_dataset")
        if os.path.exists(ncbi_folder):
            file_detector.scan_directory(ncbi_folder)
    
    genbank_file = file_detector.get_genbank_file()
    
    if not genbank_file:
        raise HTTPException(
            status_code=404,
            detail="No se encontró archivo GenBank (.gbff). Asegúrese de que los datos estén extraídos."
        )
    
    return genbank_file.filepath


@router.post("/codons", response_model=CodonAnalysisResponse)
async def analyze_codons(request: Request):
    """
    Execute codon analysis on the genome
    """
    global _analysis_cache
    
    # Check cache
    if _analysis_cache["codons"] is not None:
        return _analysis_cache["codons"]
    
    genbank_path = get_genbank_path(request)
    cache_dir = request.app.state.cache_dir
    
    # Parse genome
    parser = GenomeParser(cache_dir)
    genome_data = parser.parse_genbank(genbank_path)
    _analysis_cache["genome_data"] = genome_data
    
    # Analyze codons
    analyzer = CodonAnalyzer()
    result = analyzer.analyze(genome_data.sequence, len(genome_data.genes))
    
    response = CodonAnalysisResponse(
        genome_length=result.genome_length,
        atg_count=result.atg_count,
        atg_density=result.atg_density,
        stop_codons={
            codon: StopCodonInfo(count=data["count"], percentage=data["percentage"])
            for codon, data in result.stop_codons.items()
        },
        gene_comparison=GeneComparison(
            annotated_genes=result.gene_comparison["annotated_genes"],
            atg_found=result.gene_comparison["atg_found"],
            difference=result.gene_comparison["difference"]
        )
    )
    
    _analysis_cache["codons"] = response
    return response


@router.post("/genes", response_model=GeneAnalysisResponse)
async def analyze_genes(request: Request):
    """
    Execute gene analysis on the genome
    """
    global _analysis_cache
    
    # Check cache
    if _analysis_cache["genes"] is not None:
        return _analysis_cache["genes"]
    
    genbank_path = get_genbank_path(request)
    cache_dir = request.app.state.cache_dir
    
    # Parse genome if not cached
    if _analysis_cache["genome_data"] is None:
        parser = GenomeParser(cache_dir)
        genome_data = parser.parse_genbank(genbank_path)
        _analysis_cache["genome_data"] = genome_data
    else:
        genome_data = _analysis_cache["genome_data"]
    
    # Analyze genes
    analyzer = GeneAnalyzer()
    result = analyzer.analyze(
        genome_data.genes,
        genome_data.length,
        genome_data.cds_count,
        genome_data.gc_content
    )
    
    response = GeneAnalysisResponse(
        total_genes=result.total_genes,
        total_cds=result.total_cds,
        genome_length=result.genome_length,
        gc_content=result.gc_content,
        gene_density=result.gene_density,
        size_statistics=SizeStatistics(
            mean=result.size_statistics.mean,
            median=result.size_statistics.median,
            min=result.size_statistics.min,
            max=result.size_statistics.max,
            std=result.size_statistics.std
        ),
        genes=[
            GeneInfo(
                locus_tag=g["locus_tag"],
                start=g["start"],
                end=g["end"],
                length=g["length"],
                strand=g["strand"],
                product=g["product"],
                gc_content=g["gc_content"]
            )
            for g in result.genes
        ]
    )
    
    _analysis_cache["genes"] = response
    return response


@router.get("/validate", response_model=ValidationResponse)
async def validate_results(request: Request):
    """
    Validate analysis results against reference values
    """
    global _analysis_cache
    
    # Ensure we have gene analysis results
    if _analysis_cache["genes"] is None:
        await analyze_genes(request)
    
    genes_result = _analysis_cache["genes"]
    
    # Validate
    validator = Validator()
    calculated_values = {
        "total_genes": genes_result.total_genes,
        "genome_length": genes_result.genome_length,
        "gc_content": genes_result.gc_content,
        "total_cds": genes_result.total_cds,
        "gene_density": genes_result.gene_density
    }
    
    result = validator.validate(calculated_values)
    
    response = ValidationResponse(
        items=[
            ValidationItem(
                metric=item.metric,
                calculated=item.calculated,
                reference=item.reference,
                deviation_percent=item.deviation_percent,
                status=ValidationStatus(item.status.value)
            )
            for item in result.items
        ],
        overall_status=ValidationStatus(result.overall_status.value)
    )
    
    _analysis_cache["validation"] = response
    return response


@router.get("/complete")
async def complete_analysis(request: Request):
    """
    Execute all analyses and return combined results
    """
    global _analysis_cache
    _analysis_cache["status"] = "running"
    
    try:
        # Run all analyses
        codons = await analyze_codons(request)
        genes = await analyze_genes(request)
        validation = await validate_results(request)
        
        _analysis_cache["status"] = "completed"
        
        return {
            "codons": codons,
            "genes": {
                "total_genes": genes.total_genes,
                "total_cds": genes.total_cds,
                "genome_length": genes.genome_length,
                "gc_content": genes.gc_content,
                "gene_density": genes.gene_density,
                "size_statistics": genes.size_statistics,
                "genes_count": len(genes.genes)  # Don't return all genes here
            },
            "validation": validation
        }
    except Exception as e:
        _analysis_cache["status"] = "error"
        raise HTTPException(status_code=500, detail=str(e))


@router.get("/status", response_model=AnalysisStatus)
async def get_analysis_status():
    """
    Get status of current analysis
    """
    status_map = {
        "idle": ("idle", 0, "Sin análisis en curso"),
        "running": ("running", 50, "Análisis en progreso..."),
        "completed": ("completed", 100, "Análisis completado"),
        "error": ("error", 0, "Error en el análisis")
    }
    
    status, progress, message = status_map.get(
        _analysis_cache["status"], 
        ("unknown", 0, "Estado desconocido")
    )
    
    return AnalysisStatus(status=status, progress=progress, message=message)


@router.get("/results/codons")
async def get_codon_results():
    """
    Get cached codon analysis results
    """
    if _analysis_cache["codons"] is None:
        raise HTTPException(status_code=404, detail="No hay resultados de análisis de codones. Ejecute el análisis primero.")
    return _analysis_cache["codons"]


@router.get("/results/genes", response_model=PaginatedGenesResponse)
async def get_gene_results(page: int = 1, page_size: int = 50, search: Optional[str] = None):
    """
    Get cached gene analysis results with pagination
    """
    if _analysis_cache["genes"] is None:
        raise HTTPException(status_code=404, detail="No hay resultados de análisis de genes. Ejecute el análisis primero.")
    
    genes_result = _analysis_cache["genes"]
    all_genes = genes_result.genes
    
    # Apply search filter
    if search:
        search = search.lower()
        all_genes = [
            g for g in all_genes
            if search in g.locus_tag.lower() or 
               (g.product and search in g.product.lower())
        ]
    
    # Calculate pagination
    total = len(all_genes)
    total_pages = (total + page_size - 1) // page_size
    start = (page - 1) * page_size
    end = start + page_size
    
    return PaginatedGenesResponse(
        genes=all_genes[start:end],
        total=total,
        page=page,
        page_size=page_size,
        total_pages=total_pages
    )


@router.get("/results/statistics")
async def get_statistics():
    """
    Get general statistics from all analyses
    """
    if _analysis_cache["genes"] is None or _analysis_cache["codons"] is None:
        raise HTTPException(status_code=404, detail="Ejecute el análisis completo primero.")
    
    genes = _analysis_cache["genes"]
    codons = _analysis_cache["codons"]
    
    return {
        "genome": {
            "length": genes.genome_length,
            "gc_content": genes.gc_content
        },
        "genes": {
            "total": genes.total_genes,
            "cds": genes.total_cds,
            "density": genes.gene_density,
            "avg_length": genes.size_statistics.mean
        },
        "codons": {
            "atg_count": codons.atg_count,
            "atg_density": codons.atg_density,
            "stop_codons": codons.stop_codons
        }
    }


@router.post("/clear-cache")
async def clear_cache():
    """
    Clear all cached analysis results
    """
    global _analysis_cache
    _analysis_cache = {
        "codons": None,
        "genes": None,
        "validation": None,
        "genome_data": None,
        "status": "idle"
    }
    return {"message": "Cache limpiado correctamente"}


@router.post("/ai-validation")
async def ai_validation(api_key: Optional[str] = None):
    """
    Validar resultados usando IA (Google Gemini)
    
    Args:
        api_key: API key de Google Gemini (opcional si está en env)
    """
    if _analysis_cache["genes"] is None or _analysis_cache["codons"] is None:
        raise HTTPException(
            status_code=404, 
            detail="Ejecute el análisis completo primero."
        )
    
    # Obtener validador de IA
    validator = get_ai_validator(api_key)
    
    if validator is None:
        raise HTTPException(
            status_code=400,
            detail="API key de Gemini no configurada. Proporcione api_key o configure GEMINI_API_KEY."
        )
    
    try:
        # Preparar datos para validación
        codon_data = {
            "atg_total": _analysis_cache["codons"].atg_count,
            "atg_density": _analysis_cache["codons"].atg_density,
            "taa_total": _analysis_cache["codons"].stop_codons.taa,
            "tag_total": _analysis_cache["codons"].stop_codons.tag,
            "tga_total": _analysis_cache["codons"].stop_codons.tga,
            "total_stop_codons": (
                _analysis_cache["codons"].stop_codons.taa +
                _analysis_cache["codons"].stop_codons.tag +
                _analysis_cache["codons"].stop_codons.tga
            )
        }
        
        gene_data = {
            "total_genes": _analysis_cache["genes"].total_genes,
            "total_cds": _analysis_cache["genes"].total_cds,
            "genome_length": _analysis_cache["genes"].genome_length,
            "gc_content": _analysis_cache["genes"].gc_content,
            "gene_density": _analysis_cache["genes"].gene_density
        }
        
        # Validar con IA
        codon_validation = validator.validate_codon_analysis(codon_data)
        gene_validation = validator.validate_gene_analysis(gene_data)
        
        # Validación comprehensiva
        full_results = {
            "codons": codon_data,
            "genes": gene_data
        }
        comprehensive = validator.comprehensive_validation(full_results)
        
        return {
            "codon_validation": {
                "is_valid": codon_validation.is_valid,
                "confidence": codon_validation.confidence,
                "interpretation": codon_validation.interpretation,
                "discrepancies": codon_validation.discrepancies,
                "recommendations": codon_validation.recommendations,
                "timestamp": codon_validation.timestamp
            },
            "gene_validation": {
                "is_valid": gene_validation.is_valid,
                "confidence": gene_validation.confidence,
                "interpretation": gene_validation.interpretation,
                "discrepancies": gene_validation.discrepancies,
                "recommendations": gene_validation.recommendations,
                "timestamp": gene_validation.timestamp
            },
            "comprehensive_validation": {
                "is_valid": comprehensive.is_valid,
                "confidence": comprehensive.confidence,
                "interpretation": comprehensive.interpretation,
                "discrepancies": comprehensive.discrepancies,
                "recommendations": comprehensive.recommendations,
                "timestamp": comprehensive.timestamp
            }
        }
        
    except Exception as e:
        raise HTTPException(
            status_code=500,
            detail=f"Error en validación IA: {str(e)}"
        )

