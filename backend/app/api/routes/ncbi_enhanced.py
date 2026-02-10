"""
Enhanced NCBI API Routes
Provides endpoints for:
- Protein listing and search
- Gene location/detail lookup
- Genome sequence viewer
- Complete codon usage (RSCU, CAI, Nc)
- Genome circular map data
- Literature references search
"""
from fastapi import APIRouter, HTTPException, Request, Query
from typing import Optional, List
import os

from app.core.ncbi_service import get_ncbi_service

router = APIRouter()


def _get_genbank_path(request: Request) -> str:
    """Get the active GenBank file path"""
    file_detector = request.app.state.file_detector

    if not file_detector.detected_files:
        project_root = file_detector.project_root
        ncbi_folder = os.path.join(project_root, "ncbi_dataset")
        if os.path.exists(ncbi_folder):
            file_detector.scan_directory(ncbi_folder)

    genbank_file = file_detector.get_genbank_file()
    if not genbank_file:
        raise HTTPException(
            status_code=404,
            detail="No se encontró archivo GenBank. Active un genoma primero."
        )
    return genbank_file.filepath


# ==================== PROTEINS ====================

@router.get("/proteins")
async def get_proteins(
    request: Request,
    page: int = 1,
    page_size: int = 50,
    search: Optional[str] = None
):
    """
    Get proteins extracted from the active genome's GenBank file
    """
    genbank_path = _get_genbank_path(request)
    service = get_ncbi_service()

    proteins = service.get_proteins_from_genbank(genbank_path, limit=5000)

    # Apply search
    if search:
        search_lower = search.lower()
        proteins = [
            p for p in proteins
            if search_lower in p.get("product", "").lower()
            or search_lower in p.get("locus_tag", "").lower()
            or search_lower in p.get("gene_name", "").lower()
            or search_lower in p.get("protein_id", "").lower()
        ]

    total = len(proteins)
    total_pages = (total + page_size - 1) // page_size
    start = (page - 1) * page_size
    end = start + page_size

    # Remove full_sequence from paginated results for performance
    page_proteins = []
    for p in proteins[start:end]:
        p_copy = {k: v for k, v in p.items() if k != "full_sequence"}
        page_proteins.append(p_copy)

    return {
        "proteins": page_proteins,
        "total": total,
        "page": page,
        "page_size": page_size,
        "total_pages": total_pages
    }


@router.get("/protein/{protein_id}")
async def get_protein_detail(protein_id: str, request: Request):
    """Get full detail of a specific protein including sequence"""
    genbank_path = _get_genbank_path(request)
    service = get_ncbi_service()

    proteins = service.get_proteins_from_genbank(genbank_path, limit=10000)
    for p in proteins:
        if p.get("protein_id") == protein_id or p.get("locus_tag") == protein_id:
            return p

    raise HTTPException(status_code=404, detail=f"Proteína no encontrada: {protein_id}")


# ==================== GENE LOCATION ====================

@router.get("/gene-at-position/{position}")
async def get_gene_at_position(position: int, request: Request):
    """
    Find which gene is at a specific genomic position
    """
    genbank_path = _get_genbank_path(request)
    service = get_ncbi_service()

    result = service.get_gene_at_position(genbank_path, position)
    return result


@router.get("/gene-detail/{locus_tag}")
async def get_gene_detail(locus_tag: str, request: Request):
    """
    Get detailed information about a specific gene
    """
    genbank_path = _get_genbank_path(request)
    service = get_ncbi_service()

    result = service.get_gene_detail(genbank_path, locus_tag)
    if result is None:
        raise HTTPException(status_code=404, detail=f"Gen no encontrado: {locus_tag}")

    return result


# ==================== SEQUENCE VIEWER ====================

@router.get("/sequence")
async def get_sequence_segment(
    request: Request,
    start: int = Query(0, ge=0),
    end: int = Query(1000, ge=1),
):
    """
    Get a segment of the genome sequence with annotations
    """
    genbank_path = _get_genbank_path(request)
    service = get_ncbi_service()

    # Limit segment size to prevent memory issues
    if end - start > 10000:
        end = start + 10000

    result = service.get_genome_sequence_segment(genbank_path, start, end)
    if "error" in result:
        raise HTTPException(status_code=500, detail=result["error"])

    return result


# ==================== CODON USAGE ====================

@router.get("/codon-usage")
async def get_complete_codon_usage(request: Request):
    """
    Get complete codon usage analysis (RSCU, Nc, GC3, amino acid frequencies)
    """
    genbank_path = _get_genbank_path(request)
    service = get_ncbi_service()

    result = service.calculate_complete_codon_usage(genbank_path)
    if "error" in result:
        raise HTTPException(status_code=500, detail=result["error"])

    return result


# ==================== GENOME MAP ====================

@router.get("/genome-map")
async def get_genome_map_data(request: Request):
    """
    Get data for circular genome map visualization
    Includes genes, GC content windows, and GC skew
    """
    genbank_path = _get_genbank_path(request)
    service = get_ncbi_service()

    result = service.get_genome_map_data(genbank_path)
    if "error" in result:
        raise HTTPException(status_code=500, detail=result["error"])

    return result


# ==================== LITERATURE ====================

@router.get("/literature/search")
async def search_literature(
    query: str = Query(..., description="Search query for PubMed"),
    max_results: int = Query(8, ge=1, le=20)
):
    """
    Search PubMed for relevant scientific literature
    """
    service = get_ncbi_service()
    results = service.search_literature(query, max_results)

    return {
        "query": query,
        "count": len(results),
        "references": results
    }


@router.get("/genbank-reference/{accession}")
async def get_genbank_reference(accession: str):
    """
    Get GenBank reference information for a genome accession
    """
    service = get_ncbi_service()
    result = service.get_genbank_reference(accession)

    return result


# ==================== NCBI SEARCH ====================

@router.get("/search-gene")
async def search_ncbi_gene(
    query: str = Query(..., description="Gene name or keyword to search"),
    organism: str = Query("", description="Optional organism filter"),
    max_results: int = Query(5, ge=1, le=20)
):
    """
    Search NCBI Gene database for gene information.
    Returns gene records with NCBI links.
    """
    service = get_ncbi_service()
    results = service.search_ncbi_gene(query, organism, max_results)

    return {
        "query": query,
        "organism": organism,
        "count": len(results),
        "genes": results
    }


@router.get("/search-nucleotide")
async def search_ncbi_nucleotide(
    query: str = Query(..., description="Accession or keyword to search"),
    max_results: int = Query(3, ge=1, le=10)
):
    """
    Search NCBI Nucleotide database.
    Returns sequences with GenBank links.
    """
    service = get_ncbi_service()
    results = service.search_ncbi_nucleotide(query, max_results)

    return {
        "query": query,
        "count": len(results),
        "nucleotides": results
    }


# ==================== CENTRAL DOGMA ====================

@router.get("/central-dogma/{locus_tag}")
async def get_central_dogma(locus_tag: str, request: Request):
    """
    Get central dogma data for a gene:
    DNA (5'→3' / 3'→5') → mRNA (transcription) → Protein (translation)
    """
    genbank_path = _get_genbank_path(request)
    service = get_ncbi_service()

    result = service.get_central_dogma_data(genbank_path, locus_tag)
    if result is None:
        raise HTTPException(status_code=404, detail=f"Gen CDS no encontrado: {locus_tag}")

    return result


# ==================== GC SLIDING WINDOW ====================

@router.get("/gc-window")
async def get_gc_sliding_window(
    request: Request,
    window_size: int = 5000,
    step: int = 1000
):
    """
    GC content sliding window analysis across the genome.
    Returns GC% and GC skew for visualization.
    """
    genbank_path = _get_genbank_path(request)
    service = get_ncbi_service()

    result = service.get_gc_sliding_window(genbank_path, window_size, step)
    if result is None:
        raise HTTPException(status_code=500, detail="Error calculando ventana deslizante GC")

    return result


# ==================== BLAST SEARCH ====================

@router.post("/blast/submit")
async def submit_blast(request: Request):
    """
    Submit a BLAST search to NCBI.
    Body should contain: { "sequence": "ATGCG...", "program": "blastn", "database": "nt" }
    """
    body = await request.json()
    sequence = body.get("sequence", "")
    program = body.get("program", "blastn")
    database = body.get("database", "nt")
    max_hits = body.get("max_hits", 10)

    if not sequence or len(sequence) < 10:
        raise HTTPException(status_code=400, detail="Secuencia demasiado corta (mínimo 10 bases)")

    service = get_ncbi_service()
    result = service.submit_blast_search(sequence, program, database, max_hits)

    if result.get("status") == "ERROR":
        raise HTTPException(status_code=500, detail=result.get("error", "Error BLAST"))

    return result


@router.get("/blast/results/{rid}")
async def get_blast_results(rid: str):
    """
    Check BLAST results for a given RID.
    """
    service = get_ncbi_service()
    result = service.get_blast_results(rid)
    return result


# ==================== tRNA / rRNA ANALYSIS ====================

@router.get("/rna-analysis")
async def get_rna_analysis(request: Request):
    """
    Analyze tRNA and rRNA genes from GenBank file.
    """
    genbank_path = _get_genbank_path(request)
    service = get_ncbi_service()

    result = service.get_trna_rrna_analysis(genbank_path)
    if result is None:
        raise HTTPException(status_code=500, detail="Error analizando tRNA/rRNA")

    return result


# ==================== GENOME SUMMARY ====================

@router.get("/genome-summary")
async def get_genome_summary(request: Request):
    """
    Quick genome summary with key stats for the active genome.
    """
    genbank_path = _get_genbank_path(request)
    try:
        from Bio import SeqIO
        from Bio.SeqUtils import gc_fraction
        records = list(SeqIO.parse(genbank_path, "genbank"))
        if not records:
            raise HTTPException(status_code=404, detail="No records found in GenBank file")
        record = max(records, key=lambda r: len(r.seq))

        cds_count = sum(1 for f in record.features if f.type == "CDS")
        gene_count = sum(1 for f in record.features if f.type == "gene")
        trna_count = sum(1 for f in record.features if f.type == "tRNA")
        rrna_count = sum(1 for f in record.features if f.type == "rRNA")

        genome_len = len(record.seq)
        gc = round(gc_fraction(record.seq) * 100, 2)

        return {
            "organism": record.annotations.get("organism", ""),
            "accession": record.id,
            "description": record.description,
            "genome_length": genome_len,
            "genome_length_mb": round(genome_len / 1e6, 3),
            "gc_content": gc,
            "total_features": len(record.features),
            "total_genes": gene_count,
            "total_cds": cds_count,
            "total_trna": trna_count,
            "total_rrna": rrna_count,
            "gene_density": round((gene_count / genome_len) * 1e6, 1) if genome_len > 0 else 0,
            "taxonomy": record.annotations.get("taxonomy", []),
            "references": [
                {
                    "title": ref.title,
                    "authors": ref.authors,
                    "journal": ref.journal
                }
                for ref in record.annotations.get("references", [])[:3]
            ]
        }
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


# ==================== FUNCTIONAL CATEGORIES (COG) ====================

@router.get("/functional-categories")
async def get_functional_categories(request: Request):
    """
    Classify genes into COG-like functional categories.
    """
    genbank_path = _get_genbank_path(request)
    service = get_ncbi_service()

    result = service.get_functional_categories(genbank_path)
    if result is None:
        raise HTTPException(status_code=500, detail="Error clasificando genes en categorías funcionales")

    return result


# ==================== CAI PER GENE ====================

@router.get("/cai-analysis")
async def get_cai_analysis(request: Request, top_n: int = 50):
    """
    Calculate Codon Adaptation Index (CAI) for each gene.
    Identifies highly and lowly expressed genes.
    """
    genbank_path = _get_genbank_path(request)
    service = get_ncbi_service()

    result = service.calculate_cai_per_gene(genbank_path, top_n)
    if result is None:
        raise HTTPException(status_code=500, detail="Error calculando CAI")

    return result


# ==================== PHYLOGENETIC TREE ====================

@router.post("/phylogenetic-tree")
async def get_phylogenetic_tree(request: Request):
    """
    Calculate phylogenetic distances and build UPGMA tree.
    Body should contain: { "genbank_paths": [...] } or will use all available genomes.
    """
    import os
    import glob

    # Try to get paths from request body
    try:
        body = await request.json()
        paths = body.get("genbank_paths", [])
    except:
        paths = []

    # If no paths provided, find all GenBank files
    if not paths:
        file_detector = request.app.state.file_detector
        project_root = file_detector.project_root

        # Search in multiple locations
        search_dirs = [
            os.path.join(project_root, "genomes"),      # Main genomes directory
            os.path.join(project_root, "ncbi_dataset"), # Legacy location
            os.path.join(project_root, "data"),         # Alternative location
            request.app.state.cache_dir                 # Cache directory
        ]

        # Include all common GenBank extensions including .gbff (NCBI format)
        extensions = ["*.gbff", "*.gbk", "*.gb", "*.genbank"]

        print(f"🔍 [PHYLO] Buscando genomas en proyecto: {project_root}")

        for search_dir in search_dirs:
            if os.path.exists(search_dir):
                for ext in extensions:
                    pattern = os.path.join(search_dir, "**", ext)
                    found = glob.glob(pattern, recursive=True)
                    paths.extend(found)
                    if found:
                        print(f"  ✓ Encontrados {len(found)} archivos {ext} en {os.path.basename(search_dir)}")
            else:
                print(f"  ✗ Directorio no existe: {os.path.basename(search_dir)}")

        # Remove duplicates
        paths = list(set(paths))

        print(f"📊 [PHYLO] Total de genomas encontrados: {len(paths)}")
        for p in paths[:10]:  # Mostrar solo los primeros 10
            print(f"  - {os.path.basename(os.path.dirname(p))}/{os.path.basename(p)}")
        if len(paths) > 10:
            print(f"  ... y {len(paths) - 10} más")

    if len(paths) == 0:
        raise HTTPException(
            status_code=400,
            detail=f"No se encontraron genomas. Descargue al menos 1 genoma desde la pestaña de archivos. Encontrados: {len(paths)}"
        )

    # Handle single genome case
    if len(paths) == 1:
        from Bio import SeqIO
        try:
            record = next(SeqIO.parse(paths[0], "genbank"))
            genome_name = record.description[:50]
            accession = record.id
            gc = sum(1 for base in str(record.seq).upper() if base in 'GC') / len(record.seq) * 100 if len(record.seq) > 0 else 0
            gene_count = sum(1 for f in record.features if f.type == "gene")

            return {
                "genomes": [{
                    "name": genome_name,
                    "accession": accession,
                    "gc": round(gc, 2),
                    "length": len(record.seq),
                    "gene_count": gene_count
                }],
                "distance_matrix": [[0.0]],
                "tree": {
                    "name": genome_name,
                    "height": 0,
                    "branch_length": 0,
                    "children": []
                },
                "labels": [genome_name],
                "method": "Single Genome (No clustering needed)",
                "single_genome": True
            }
        except Exception as e:
            raise HTTPException(status_code=500, detail=f"Error procesando genoma: {str(e)}")

    # Multiple genomes case
    service = get_ncbi_service()
    result = service.calculate_phylogenetic_distances(paths)

    if result is None:
        raise HTTPException(status_code=500, detail="Error calculando distancias filogenéticas")

    if "error" in result:
        raise HTTPException(status_code=400, detail=result["error"])

    return result
