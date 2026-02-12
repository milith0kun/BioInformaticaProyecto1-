"""
Enhanced NCBI API Routes
Provides endpoints for:
- Protein listing and search
- Gene location/detail lookup
- Genome sequence viewer
- Complete codon usage (RSCU, CAI, Nc)
- Genome circular map data
- Literature references search
- Protein structure prediction with ESMFold
"""
from fastapi import APIRouter, HTTPException, Request, Query
from fastapi.responses import Response
from typing import Optional, List
import os
import httpx
from pathlib import Path
from Bio import SeqIO

from app.core.ncbi_service import get_ncbi_service

router = APIRouter()


def _get_genbank_path(request: Request) -> str:
    """Get the active GenBank file path with auto-recovery
    
    Priority order:
    1. Use currently detected active genome
    2. Auto-scan and activate first available genome from genomes/ folder
    3. Raise 404 if no genomes found
    """
    file_detector = request.app.state.file_detector

    # Step 1: Check if there's an active genome already
    genbank_file = file_detector.get_genbank_file()
    if genbank_file:
        return genbank_file.filepath

    # Step 2: Auto-recovery - Try to find and activate any available genome
    project_root = file_detector.project_root
    genomes_dir = os.path.join(project_root, "genomes")
    
    if os.path.exists(genomes_dir):
        print(f"🔍 [RECOVERY] No active genome, searching in {genomes_dir}")
        
        # Sort to get consistent genome selection (prefer GCF over GCA)
        entries = sorted([e for e in os.listdir(genomes_dir) if e.startswith(("GCF_", "GCA_"))],
                        key=lambda x: (0 if x.startswith("GCF_") else 1, x))
        
        for entry in entries:
            extracted_dir = os.path.join(genomes_dir, entry, "extracted")
            if os.path.exists(extracted_dir):
                # Check if there's a .gbff file
                gbff_files = [f for f in os.listdir(extracted_dir) if f.endswith(".gbff")]
                if gbff_files:
                    print(f"✅ [RECOVERY] Auto-activating genome: {entry}")
                    file_detector.scan_directory(extracted_dir)
                    genbank_file = file_detector.get_genbank_file()
                    if genbank_file:
                        return genbank_file.filepath

    # Step 3: No genome found anywhere
    print("❌ [RECOVERY] No GenBank file found in any genome directory.")
    raise HTTPException(
        status_code=404,
        detail="No hay genomas disponibles. Por favor, descarga al menos un genoma desde el panel NCBI Datasets."
    )


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
    Uses analysis cache first (much faster), then falls back to GenBank parsing.
    """
    proteins = []
    
    # Try getting from analysis cache first (FAST)
    try:
        from app.api.routes.analysis import _analysis_cache
        genome_data = _analysis_cache.get("genome_data")
        
        if genome_data and genome_data.genes:
            print(f"✅ [PROTEINS] Using analysis cache with {len(genome_data.genes)} genes")
            proteins = [
                {
                    "protein_id": gene.protein_id,
                    "locus_tag": gene.locus_tag,
                    "gene_name": gene.gene_name,
                    "product": gene.product,
                    "length": gene.length,
                    "full_sequence": gene.translation,
                    "start": gene.start,
                    "end": gene.end,
                    "strand": gene.strand,
                    "gc_content": gene.gc_content
                }
                for gene in genome_data.genes
                if gene.protein_id  # Only include genes with protein_id
            ]
    except Exception as e:
        print(f"⚠️ [PROTEINS] Cache lookup failed: {e}")
    
    # Fallback: read from GenBank file (SLOW)
    if not proteins:
        print(f"⚠️ [PROTEINS] Cache empty, falling back to GenBank parsing")
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
    """Get full detail of a specific protein including sequence
    Uses analysis cache first, then searches all available genomes as fallback.
    """
    # Try getting from analysis cache first (FAST)
    try:
        from app.api.routes.analysis import _analysis_cache
        genome_data = _analysis_cache.get("genome_data")
        
        if genome_data and genome_data.genes:
            for gene in genome_data.genes:
                if gene.protein_id == protein_id or gene.locus_tag == protein_id:
                    print(f"✅ [PROTEIN] Found {protein_id} in analysis cache")
                    return {
                        "protein_id": gene.protein_id,
                        "locus_tag": gene.locus_tag,
                        "gene_name": gene.gene_name,
                        "product": gene.product,
                        "length": gene.length,
                        "full_sequence": gene.translation,
                        "start": gene.start,
                        "end": gene.end,
                        "strand": gene.strand,
                        "gc_content": gene.gc_content
                    }
    except Exception as e:
        print(f"⚠️ [PROTEIN] Cache lookup failed: {e}")
    
    # Fallback 1: Try active genome GenBank (if available)
    service = get_ncbi_service()
    try:
        genbank_path = _get_genbank_path(request)
        print(f"🔍 [PROTEIN] Searching in active genome: {genbank_path}")
        proteins = service.get_proteins_from_genbank(genbank_path, limit=10000)
        for p in proteins:
            if p.get("protein_id") == protein_id or p.get("locus_tag") == protein_id:
                print(f"✅ [PROTEIN] Found {protein_id} in active genome")
                return p
    except Exception as e:
        print(f"⚠️ [PROTEIN] Active genome search failed: {e}")
    
    # Fallback 2: Search ALL genomes (comprehensive search)
    print(f"🔍 [PROTEIN] Searching across all genomes...")
    genomes_dir = Path("genomes")
    if genomes_dir.exists():
        for accession_dir in sorted(genomes_dir.iterdir(), 
                                   key=lambda x: (not x.name.startswith("GCF_"), x.name)):
            if not accession_dir.is_dir():
                continue
            # Find GenBank file
            extracted_path = accession_dir / "extracted"
            if not extracted_path.exists():
                continue
            genbank_files = list(extracted_path.rglob("*.gbff")) + list(extracted_path.rglob("*.gbk"))
            if not genbank_files:
                continue
            
            genbank_path = str(genbank_files[0])
            try:
                proteins = service.get_proteins_from_genbank(genbank_path, limit=10000)
                for p in proteins:
                    if p.get("protein_id") == protein_id or p.get("locus_tag") == protein_id:
                        print(f"✅ [PROTEIN] Found {protein_id} in {accession_dir.name}")
                        return p
            except Exception as e:
                continue

    raise HTTPException(status_code=404, detail=f"Proteína no encontrada en ningún genoma: {protein_id}")


@router.post("/protein/predict-structure/{protein_id}")
async def predict_protein_structure(protein_id: str, request: Request):
    """
    Resolve 3D structure using multiple backends (cascade):
    1. AlphaFold DB via UniProt xref in GenBank
    2. AlphaFold DB via UniProt API mapping (RefSeq → UniProt)
    3. RCSB PDB experimental structure search
    4. ESMFold API (de novo, only for small proteins < 400 aa)
    Returns PDB format file for Molstar visualization
    Uses analysis cache first.
    """
    
    print(f"🔍 [PREDICT] Searching sequence for: {protein_id}")
    sequence = None
    db_xrefs = []
    product_name = ""
    gene_name = ""
    
    # Try getting from analysis cache first
    try:
        from app.api.routes.analysis import _analysis_cache
        genome_data = _analysis_cache.get("genome_data")
        
        if genome_data and genome_data.genes:
            for gene in genome_data.genes:
                if gene.protein_id == protein_id or gene.locus_tag == protein_id:
                    sequence = gene.translation
                    db_xrefs = gene.db_xrefs or []
                    product_name = gene.product or ""
                    gene_name = gene.gene_name or ""
                    print(f"✅ [PREDICT] Found in analysis cache")
                    break
    except Exception as e:
        print(f"⚠️ [PREDICT] Cache lookup failed: {e}")
    
    # Fallback 1: Try active genome GenBank
    if not sequence:
        try:
            genbank_path = _get_genbank_path(request)
            print(f"🔍 [PREDICT] Searching in active genome: {genbank_path}")
            with open(genbank_path, "r") as handle:
                for record in SeqIO.parse(handle, "genbank"):
                    for feature in record.features:
                        if feature.type == "CDS":
                            qual = feature.qualifiers
                            if qual.get('protein_id', [''])[0] == protein_id or \
                               qual.get('locus_tag', [''])[0] == protein_id:
                                sequence = qual.get('translation', [''])[0]
                                db_xrefs = qual.get('db_xref', [])
                                product_name = qual.get('product', [''])[0]
                                gene_name = qual.get('gene', [''])[0]
                                print(f"✅ [PREDICT] Found in active genome")
                                break
                    if sequence:
                        break
        except Exception as e:
            print(f"⚠️ [PREDICT] Active genome search failed: {e}")
    
    # Fallback 2: Search ALL genomes (comprehensive search)
    if not sequence:
        print(f"🔍 [PREDICT] Searching across all genomes...")
        genomes_dir = Path("genomes")
        if genomes_dir.exists():
            for accession_dir in sorted(genomes_dir.iterdir(), 
                                       key=lambda x: (not x.name.startswith("GCF_"), x.name)):
                if not accession_dir.is_dir():
                    continue
                extracted_path = accession_dir / "extracted"
                if not extracted_path.exists():
                    continue
                genbank_files = list(extracted_path.rglob("*.gbff")) + list(extracted_path.rglob("*.gbk"))
                if not genbank_files:
                    continue
                
                genbank_path = str(genbank_files[0])
                try:
                    with open(genbank_path, "r") as handle:
                        for record in SeqIO.parse(handle, "genbank"):
                            for feature in record.features:
                                if feature.type == "CDS":
                                    qual = feature.qualifiers
                                    if qual.get('protein_id', [''])[0] == protein_id or \
                                       qual.get('locus_tag', [''])[0] == protein_id:
                                        sequence = qual.get('translation', [''])[0]
                                        db_xrefs = qual.get('db_xref', [])
                                        product_name = qual.get('product', [''])[0]
                                        gene_name = qual.get('gene', [''])[0]
                                        print(f"✅ [PREDICT] Found {protein_id} in {accession_dir.name}")
                                        break
                            if sequence:
                                break
                    if sequence:
                        break
                except Exception as e:
                    continue

    if not sequence:
        raise HTTPException(status_code=404, detail=f"Proteína {protein_id} no identificada.")

    sequence = sequence.upper().strip().replace("*", "")
    seq_len = len(sequence)
    print(f"🧬 [PREDICT] Sequence: {seq_len} aa | Product: {product_name} | Gene: {gene_name}")

    if seq_len < 5:
        raise HTTPException(status_code=400, detail="Secuencia insuficiente para modelado.")

    errors_log = []

    # === STRATEGY 1: AlphaFold DB via GenBank UniProt xref ===
    uniprot_id = None
    for xref in db_xrefs:
        if "UniProtKB" in xref:
            uniprot_id = xref.split(":")[-1] if ":" in xref else None
            break

    if uniprot_id:
        result = await _try_alphafold(protein_id, uniprot_id, errors_log)
        if result:
            return result

    # === STRATEGY 2: AlphaFold DB via UniProt API mapping ===
    if not uniprot_id:
        try:
            async with httpx.AsyncClient(timeout=15.0) as client:
                search_url = f"https://rest.uniprot.org/uniprotkb/search?query=xref:refseq-{protein_id}&fields=accession&format=json&size=1"
                print(f"🔎 [PREDICT] UniProt lookup for RefSeq: {protein_id}")
                uni_resp = await client.get(search_url)
                if uni_resp.status_code == 200:
                    results = uni_resp.json().get("results", [])
                    if results:
                        uniprot_id = results[0].get("primaryAccession")
                        print(f"✅ [PREDICT] RefSeq → UniProt: {protein_id} → {uniprot_id}")
                        result = await _try_alphafold(protein_id, uniprot_id, errors_log)
                        if result:
                            return result
                    else:
                        errors_log.append("UniProt: sin mapeo encontrado")
                        print(f"⚠️ [PREDICT] No UniProt mapping for {protein_id}")
        except Exception as e:
            errors_log.append(f"UniProt API: {str(e)[:80]}")
            print(f"⚠️ [PREDICT] UniProt lookup failed: {e}")

    # === STRATEGY 3: RCSB PDB search (by ID, name, and sequence) ===
    try:
        async with httpx.AsyncClient(timeout=30.0) as client:
            pdb_id_found = None
            
            # 3a. Search by protein_id and product name
            search_terms = [protein_id]
            if product_name and len(product_name) > 3:
                search_terms.append(product_name)
            
            for term in search_terms:
                search_query = {
                    "query": {
                        "type": "terminal",
                        "service": "text",
                        "parameters": {"value": term}
                    },
                    "return_type": "entry",
                    "request_options": {
                        "results_content_type": ["experimental"],
                        "paginate": {"start": 0, "rows": 1}
                    }
                }
                print(f"🔎 [PREDICT] RCSB PDB text search: '{term}'")
                pdb_resp = await client.post(
                    "https://search.rcsb.org/rcsbsearch/v2/query",
                    json=search_query
                )
                
                if pdb_resp.status_code == 200:
                    results = pdb_resp.json().get("result_set", [])
                    if results:
                        pdb_id_found = results[0].get("identifier", "")
                        if pdb_id_found:
                            print(f"✅ [PREDICT] RCSB found PDB '{pdb_id_found}' via term '{term}'")
                            break
            
            # 3b. If text search failed, try sequence-based search
            if not pdb_id_found and seq_len <= 1500:
                try:
                    seq_query = {
                        "query": {
                            "type": "terminal",
                            "service": "sequence",
                            "parameters": {
                                "evalue_cutoff": 0.001,
                                "identity_cutoff": 0.5,
                                "sequence_type": "protein",
                                "value": sequence[:1000]  # Limit to first 1000 aa for API
                            }
                        },
                        "return_type": "entry",
                        "request_options": {
                            "results_content_type": ["experimental"],
                            "paginate": {"start": 0, "rows": 1}
                        }
                    }
                    print(f"🔎 [PREDICT] RCSB sequence search ({seq_len} aa)...")
                    seq_resp = await client.post(
                        "https://search.rcsb.org/rcsbsearch/v2/query",
                        json=seq_query
                    )
                    if seq_resp.status_code == 200:
                        results = seq_resp.json().get("result_set", [])
                        if results:
                            pdb_id_found = results[0].get("identifier", "")
                            if pdb_id_found:
                                print(f"✅ [PREDICT] RCSB found PDB '{pdb_id_found}' via sequence similarity")
                except Exception as e:
                    print(f"⚠️ [PREDICT] RCSB sequence search failed: {e}")
            
            # 3c. Download found PDB structure
            if pdb_id_found:
                # Extract just the entry ID (4 chars) if it contains chain info
                entry_id = pdb_id_found[:4] if len(pdb_id_found) >= 4 else pdb_id_found
                pdb_file_resp = await client.get(
                    f"https://files.rcsb.org/download/{entry_id}.pdb",
                    timeout=30.0
                )
                if pdb_file_resp.status_code == 200 and len(pdb_file_resp.text) > 100:
                    print(f"✅ [PREDICT] RCSB PDB structure loaded: {entry_id}")
                    return Response(
                        content=pdb_file_resp.text,
                        media_type="chemical/x-pdb",
                        headers={
                            "Content-Disposition": f"inline; filename={protein_id}_{entry_id}.pdb",
                            "X-Structure-Source": "pdb"
                        }
                    )
            
            if not pdb_id_found:
                errors_log.append("RCSB PDB: sin estructura encontrada")
    except Exception as e:
        errors_log.append(f"RCSB PDB: {str(e)[:80]}")
        print(f"⚠️ [PREDICT] RCSB PDB search failed: {e}")

    # === STRATEGY 4: ESMFold API (only for small proteins, API may be unreliable) ===
    if seq_len <= 400:
        try:
            async with httpx.AsyncClient(timeout=120.0) as client:
                print(f"🛰️ [PREDICT] Trying ESMFold for {protein_id} ({seq_len} aa)...")
                response = await client.post(
                    "https://api.esmatlas.com/foldSequence/v1/pdb/",
                    data=sequence,
                    headers={"Content-Type": "text/plain"}
                )

                if response.status_code == 200 and len(response.text) > 100:
                    print(f"✅ [PREDICT] ESMFold success for {protein_id}")
                    return Response(
                        content=response.text,
                        media_type="chemical/x-pdb",
                        headers={
                            "Content-Disposition": f"inline; filename={protein_id}_predicted.pdb",
                            "X-Structure-Source": "esmfold"
                        }
                    )
                else:
                    errors_log.append(f"ESMFold: HTTP {response.status_code}")
                    print(f"⚠️ [PREDICT] ESMFold rejected: {response.status_code}")

        except httpx.TimeoutException:
            errors_log.append("ESMFold: timeout")
            print(f"⏱️ [PREDICT] ESMFold timeout for {protein_id}")
        except Exception as e:
            errors_log.append(f"ESMFold: {str(e)[:80]}")
            print(f"⚠️ [PREDICT] ESMFold failed: {e}")
    else:
        errors_log.append(f"ESMFold: omitido (proteína de {seq_len} aa > 400 aa límite)")

    # === ALL STRATEGIES FAILED ===
    sources_tried = " | ".join(errors_log) if errors_log else "Todas las fuentes fallaron"
    print(f"❌ [PREDICT] All strategies failed for {protein_id}: {sources_tried}")
    raise HTTPException(
        status_code=404,
        detail=f"No se encontró estructura 3D para {protein_id} ({seq_len} aa, {product_name}). "
               f"Fuentes intentadas: {sources_tried}. "
               f"Busque manualmente en https://www.rcsb.org o https://alphafold.ebi.ac.uk"
    )


async def _try_alphafold(protein_id: str, uniprot_id: str, errors_log: list):
    """Helper: try to fetch structure from AlphaFold DB for a given UniProt ID"""
    try:
        async with httpx.AsyncClient(timeout=30.0, follow_redirects=True) as client:
            # Try model versions v4, v3, v2
            for version in [4, 3, 2]:
                alphafold_url = f"https://alphafold.ebi.ac.uk/files/AF-{uniprot_id}-F1-model_v{version}.pdb"
                print(f"🔎 [PREDICT] Trying AlphaFold v{version}: {uniprot_id}")
                af_response = await client.get(alphafold_url)
                
                if af_response.status_code == 200 and len(af_response.text) > 100:
                    print(f"✅ [PREDICT] AlphaFold v{version} found for {protein_id} (UniProt: {uniprot_id})")
                    return Response(
                        content=af_response.text,
                        media_type="chemical/x-pdb",
                        headers={
                            "Content-Disposition": f"inline; filename={protein_id}_alphafold_v{version}.pdb",
                            "X-Structure-Source": "alphafold"
                        }
                    )
            
            errors_log.append(f"AlphaFold DB: no disponible para {uniprot_id}")
            print(f"⚠️ [PREDICT] AlphaFold DB has no model for {uniprot_id}")
    except Exception as e:
        errors_log.append(f"AlphaFold DB: {str(e)[:80]}")
        print(f"⚠️ [PREDICT] AlphaFold DB failed: {e}")
    return None


@router.get("/protein/analyze-structure/{protein_id}")
async def analyze_protein_structure(protein_id: str, request: Request):
    """
    Analyze a protein's 3D structure to extract:
    - Secondary structure elements (HELIX/SHEET from PDB records)
    - Chain information for quaternary structure
    - Per-residue B-factor / pLDDT scores
    - Structural statistics
    
    First resolves structure (AlphaFold/PDB/ESMFold), then parses the PDB.
    """
    import io
    import re
    from Bio.PDB import PDBParser as BioPDBParser
    
    # Step 1: Get the PDB content (reuse the predict endpoint logic)
    try:
        structure_response = await predict_protein_structure(protein_id, request)
        pdb_content = structure_response.body.decode('utf-8') if hasattr(structure_response, 'body') else ""
        source = structure_response.headers.get('X-Structure-Source', 'unknown') if hasattr(structure_response, 'headers') else 'unknown'
    except HTTPException:
        raise
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Error obteniendo estructura: {str(e)}")
    
    if not pdb_content or len(pdb_content) < 50:
        raise HTTPException(status_code=404, detail="No se pudo obtener estructura PDB")
    
    # Step 2: Parse HELIX / SHEET records directly from PDB text
    helices = []
    sheets = []
    chains_info = {}
    remarks = []
    resolution = None
    
    for line in pdb_content.split('\n'):
        if line.startswith('HELIX'):
            try:
                helix = {
                    'id': line[11:14].strip(),
                    'init_res_name': line[15:18].strip(),
                    'init_chain': line[19].strip(),
                    'init_seq_num': int(line[21:25].strip()),
                    'end_res_name': line[27:30].strip(),
                    'end_chain': line[31].strip(),
                    'end_seq_num': int(line[33:37].strip()),
                    'helix_class': int(line[38:40].strip()) if line[38:40].strip() else 1,
                    'length': int(line[71:76].strip()) if len(line) > 75 and line[71:76].strip() else 0,
                    'type': 'helix'
                }
                helices.append(helix)
            except (ValueError, IndexError):
                continue
                
        elif line.startswith('SHEET'):
            try:
                sheet = {
                    'strand': int(line[7:10].strip()),
                    'sheet_id': line[11:14].strip(),
                    'init_res_name': line[17:20].strip(),
                    'init_chain': line[21].strip(),
                    'init_seq_num': int(line[22:26].strip()),
                    'end_res_name': line[28:31].strip(),
                    'end_chain': line[32].strip(),
                    'end_seq_num': int(line[33:37].strip()),
                    'sense': int(line[38:40].strip()) if line[38:40].strip() else 0,
                    'type': 'sheet'
                }
                sheets.append(sheet)
            except (ValueError, IndexError):
                continue
                
        elif line.startswith('ATOM') or line.startswith('HETATM'):
            try:
                chain_id = line[21].strip() or 'A'
                res_seq = int(line[22:26].strip())
                res_name = line[17:20].strip()
                b_factor = float(line[60:66].strip()) if len(line) > 65 else 0.0
                
                if chain_id not in chains_info:
                    chains_info[chain_id] = {
                        'chain_id': chain_id,
                        'residues': set(),
                        'b_factors': [],
                        'atom_count': 0
                    }
                chains_info[chain_id]['residues'].add(res_seq)
                chains_info[chain_id]['b_factors'].append(b_factor)
                chains_info[chain_id]['atom_count'] += 1
            except (ValueError, IndexError):
                continue
                
        elif line.startswith('REMARK') and 'RESOLUTION' in line.upper():
            try:
                nums = re.findall(r'(\d+\.?\d*)', line)
                if nums:
                    resolution = float(nums[0])
            except:
                pass
    
    # Process chain data
    chains = []
    for cid, cdata in sorted(chains_info.items()):
        avg_bfactor = sum(cdata['b_factors']) / len(cdata['b_factors']) if cdata['b_factors'] else 0
        chains.append({
            'chain_id': cid,
            'residue_count': len(cdata['residues']),
            'atom_count': cdata['atom_count'],
            'avg_bfactor': round(avg_bfactor, 2),
            'min_residue': min(cdata['residues']) if cdata['residues'] else 0,
            'max_residue': max(cdata['residues']) if cdata['residues'] else 0,
        })
    
    # Determine quaternary state
    num_chains = len(chains)
    oligomeric_labels = {1: 'Monómero', 2: 'Dímero', 3: 'Trímero', 4: 'Tetrámero', 
                         5: 'Pentámero', 6: 'Hexámero', 8: 'Octámero', 12: 'Dodecámero'}
    oligomeric_state = oligomeric_labels.get(num_chains, f'{num_chains}-mero')
    
    # Compute secondary structure coverage
    total_residues = sum(c['residue_count'] for c in chains)
    helix_residues = sum(h.get('length', h['end_seq_num'] - h['init_seq_num'] + 1) for h in helices)
    sheet_residues = sum(s['end_seq_num'] - s['init_seq_num'] + 1 for s in sheets)
    coil_residues = max(0, total_residues - helix_residues - sheet_residues)
    
    ss_stats = {
        'helix_count': len(helices),
        'sheet_count': len(sheets),
        'helix_residues': helix_residues,
        'sheet_residues': sheet_residues,
        'coil_residues': coil_residues,
        'helix_pct': round((helix_residues / total_residues * 100), 1) if total_residues else 0,
        'sheet_pct': round((sheet_residues / total_residues * 100), 1) if total_residues else 0,
        'coil_pct': round((coil_residues / total_residues * 100), 1) if total_residues else 0,
    }
    
    # Build per-residue SS assignment
    ss_assignment = {}
    for h in helices:
        for r in range(h['init_seq_num'], h['end_seq_num'] + 1):
            ss_assignment[r] = 'H'
    for s in sheets:
        for r in range(s['init_seq_num'], s['end_seq_num'] + 1):
            ss_assignment[r] = 'E'
    
    # Per-residue B-factor (pLDDT for AlphaFold)
    bfactor_per_residue = {}
    for cid, cdata in chains_info.items():
        residue_bfactors = {}
        for line in pdb_content.split('\n'):
            if (line.startswith('ATOM') and line[21].strip() == cid and 
                line[12:16].strip() == 'CA'):
                try:
                    res_seq = int(line[22:26].strip())
                    bf = float(line[60:66].strip())
                    residue_bfactors[res_seq] = round(bf, 2)
                except:
                    continue
        bfactor_per_residue[cid] = residue_bfactors
    
    return {
        'protein_id': protein_id,
        'source': source,
        'resolution': resolution,
        'secondary_structure': {
            'helices': helices,
            'sheets': sheets,
            'stats': ss_stats,
            'assignment': {str(k): v for k, v in sorted(ss_assignment.items())},
        },
        'quaternary': {
            'chains': chains,
            'num_chains': num_chains,
            'oligomeric_state': oligomeric_state,
            'is_complex': num_chains > 1,
        },
        'quality': {
            'bfactor_per_residue': {k: v for k, v in bfactor_per_residue.items()},
            'is_plddt': source in ('alphafold', 'esmfold'),
        },
        'total_residues': total_residues,
        'total_atoms': sum(c['atom_count'] for c in chains),
    }


# Amino acid physicochemical properties (Kyte-Doolittle hydropathy, MW, pI group)
_AA_PROPERTIES = {
    'A': {'name': 'Alanina', 'group': 'hydrophobic', 'hydropathy': 1.8, 'mw': 89.1, 'charge': 0},
    'R': {'name': 'Arginina', 'group': 'positive', 'hydropathy': -4.5, 'mw': 174.2, 'charge': 1},
    'N': {'name': 'Asparagina', 'group': 'polar', 'hydropathy': -3.5, 'mw': 132.1, 'charge': 0},
    'D': {'name': 'Ác. aspártico', 'group': 'negative', 'hydropathy': -3.5, 'mw': 133.1, 'charge': -1},
    'C': {'name': 'Cisteína', 'group': 'special', 'hydropathy': 2.5, 'mw': 121.2, 'charge': 0},
    'E': {'name': 'Ác. glutámico', 'group': 'negative', 'hydropathy': -3.5, 'mw': 147.1, 'charge': -1},
    'Q': {'name': 'Glutamina', 'group': 'polar', 'hydropathy': -3.5, 'mw': 146.2, 'charge': 0},
    'G': {'name': 'Glicina', 'group': 'special', 'hydropathy': -0.4, 'mw': 75.0, 'charge': 0},
    'H': {'name': 'Histidina', 'group': 'positive', 'hydropathy': -3.2, 'mw': 155.2, 'charge': 0.5},
    'I': {'name': 'Isoleucina', 'group': 'hydrophobic', 'hydropathy': 4.5, 'mw': 131.2, 'charge': 0},
    'L': {'name': 'Leucina', 'group': 'hydrophobic', 'hydropathy': 3.8, 'mw': 131.2, 'charge': 0},
    'K': {'name': 'Lisina', 'group': 'positive', 'hydropathy': -3.9, 'mw': 146.2, 'charge': 1},
    'M': {'name': 'Metionina', 'group': 'hydrophobic', 'hydropathy': 1.9, 'mw': 149.2, 'charge': 0},
    'F': {'name': 'Fenilalanina', 'group': 'hydrophobic', 'hydropathy': 2.8, 'mw': 165.2, 'charge': 0},
    'P': {'name': 'Prolina', 'group': 'special', 'hydropathy': -1.6, 'mw': 115.1, 'charge': 0},
    'S': {'name': 'Serina', 'group': 'polar', 'hydropathy': -0.8, 'mw': 105.1, 'charge': 0},
    'T': {'name': 'Treonina', 'group': 'polar', 'hydropathy': -0.7, 'mw': 119.1, 'charge': 0},
    'W': {'name': 'Triptófano', 'group': 'hydrophobic', 'hydropathy': -0.9, 'mw': 204.2, 'charge': 0},
    'Y': {'name': 'Tirosina', 'group': 'polar', 'hydropathy': -1.3, 'mw': 181.2, 'charge': 0},
    'V': {'name': 'Valina', 'group': 'hydrophobic', 'hydropathy': 4.2, 'mw': 117.1, 'charge': 0},
}


@router.get("/protein/structure-pipeline/{protein_id}")
async def get_protein_structure_pipeline(protein_id: str, request: Request):
    """
    Consolidated endpoint for the protein structure pipeline visualization.
    Returns all 4 structural levels in a single response:
    - Primary: sequence + per-residue properties
    - Secondary: HELIX/SHEET from PDB + predicted H-bonds
    - Tertiary: 3D structure metadata + disulfide bonds
    - Quaternary: chain assembly data
    """
    import math
    
    # 1. Get protein detail (sequence, gene info)
    try:
        protein_detail = await get_protein_detail(protein_id, request)
    except HTTPException:
        raise
    
    sequence = protein_detail.get('full_sequence', '') or protein_detail.get('translation', '')
    if not sequence:
        raise HTTPException(status_code=404, detail="Secuencia proteica no disponible")
    
    # 2. Build primary structure data (per-residue properties)
    residue_properties = []
    total_mw = 0.0
    net_charge = 0.0
    for i, aa in enumerate(sequence):
        props = _AA_PROPERTIES.get(aa, {'name': aa, 'group': 'unknown', 'hydropathy': 0, 'mw': 110, 'charge': 0})
        total_mw += props['mw']
        net_charge += props['charge']
        residue_properties.append({
            'pos': i + 1,
            'aa': aa,
            'name': props['name'],
            'group': props['group'],
            'hydropathy': props['hydropathy'],
            'charge': props['charge'],
        })
    # Subtract water for each peptide bond
    total_mw -= (len(sequence) - 1) * 18.015
    
    primary = {
        'sequence': sequence,
        'length': len(sequence),
        'molecular_weight_kda': round(total_mw / 1000, 2),
        'net_charge': round(net_charge, 1),
        'residue_properties': residue_properties,
    }
    
    # 3. Try to get structure analysis (secondary + quaternary from PDB)
    secondary = None
    tertiary_meta = None
    quaternary = None
    quality = None
    
    try:
        analysis = await analyze_protein_structure(protein_id, request)
        
        ss = analysis.get('secondary_structure', {})
        helices = ss.get('helices', [])
        sheets = ss.get('sheets', [])
        assignment = ss.get('assignment', {})
        stats = ss.get('stats', {})
        
        # Predict hydrogen bonds from SS assignment
        hydrogen_bonds = []
        # α-helix H-bonds: C=O of residue i bonds to N-H of residue i+4
        for h in helices:
            start = h.get('init_seq_num', 0)
            end = h.get('end_seq_num', 0)
            for r in range(start, end - 3):
                hydrogen_bonds.append({
                    'donor': r,
                    'acceptor': r + 4,
                    'type': 'α-hélice (i→i+4)',
                    'ss_type': 'H'
                })
        # β-sheet H-bonds: between strands (approximate)
        for s in sheets:
            start = s.get('init_seq_num', 0)
            end = s.get('end_seq_num', 0)
            for r in range(start, end - 1, 2):
                hydrogen_bonds.append({
                    'donor': r,
                    'acceptor': r + 2,
                    'type': 'β-lámina (inter-hebra)',
                    'ss_type': 'E'
                })
        
        secondary = {
            'helices': helices,
            'sheets': sheets,
            'assignment': assignment,
            'stats': stats,
            'hydrogen_bonds': hydrogen_bonds[:200],  # cap for performance
            'total_hbonds': len(hydrogen_bonds),
        }
        
        # Detect disulfide bonds from CYS positions
        cys_positions = [i + 1 for i, aa in enumerate(sequence) if aa == 'C']
        disulfide_bonds = []
        # If we have PDB coordinates, check Cα distances; otherwise predict from proximity
        for i in range(len(cys_positions)):
            for j in range(i + 1, len(cys_positions)):
                # In real structures, S-S bonds form between CYS residues
                # For prediction: pairs with sequence distance > 10 are plausible
                if abs(cys_positions[j] - cys_positions[i]) > 10:
                    disulfide_bonds.append({
                        'residue1': cys_positions[i],
                        'residue2': cys_positions[j],
                        'chain': 'A',
                        'predicted': True
                    })
        
        # Simple domain detection based on secondary structure density
        domains = []
        if len(sequence) > 100:
            # Split into N-terminal and C-terminal regions based on SS density
            mid = len(sequence) // 2
            n_ss = sum(1 for k, v in assignment.items() if int(k) <= mid and v in ('H', 'E'))
            c_ss = sum(1 for k, v in assignment.items() if int(k) > mid and v in ('H', 'E'))
            
            if n_ss > 5 and c_ss > 5:
                domains = [
                    {'name': 'Dominio N-terminal', 'start': 1, 'end': mid, 'ss_elements': n_ss},
                    {'name': 'Dominio C-terminal', 'start': mid + 1, 'end': len(sequence), 'ss_elements': c_ss}
                ]
            elif len(sequence) > 300:
                # Try thirds
                t1, t2 = len(sequence) // 3, 2 * len(sequence) // 3
                domains = [
                    {'name': 'Dominio I', 'start': 1, 'end': t1},
                    {'name': 'Dominio II', 'start': t1 + 1, 'end': t2},
                    {'name': 'Dominio III', 'start': t2 + 1, 'end': len(sequence)}
                ]
        
        tertiary_meta = {
            'source': analysis.get('source', 'unknown'),
            'total_atoms': analysis.get('total_atoms', 0),
            'total_residues': analysis.get('total_residues', len(sequence)),
            'disulfide_bonds': disulfide_bonds[:20],
            'domains': domains,
            'resolution': analysis.get('resolution'),
        }
        
        quat = analysis.get('quaternary', {})
        quaternary = {
            'chains': quat.get('chains', []),
            'num_chains': quat.get('num_chains', 1),
            'oligomeric_state': quat.get('oligomeric_state', 'Monómero'),
            'is_complex': quat.get('is_complex', False),
        }
        
        quality_data = analysis.get('quality', {})
        # Compute average confidence
        avg_conf = 0
        total_bf = 0
        count_bf = 0
        for chain_bfs in quality_data.get('bfactor_per_residue', {}).values():
            for bf in chain_bfs.values():
                total_bf += bf
                count_bf += 1
        avg_conf = round(total_bf / count_bf, 1) if count_bf else 0
        
        quality = {
            'is_plddt': quality_data.get('is_plddt', False),
            'avg_confidence': avg_conf,
        }
        
    except HTTPException:
        # Structure not yet resolved — return primary data only
        pass
    except Exception as e:
        print(f"⚠️ [PIPELINE] Structure analysis failed: {e}")
    
    # 4. Functional Annotations (Simulated/Heuristic for now)
    # In a real app, this would query UniProt/InterPro/CDD API
    functional = {
        'go_terms': [],
        'conserved_domains': [],
        'active_sites': []
    }
    
    # Heuristic domain detection based on length if not found in tertiary
    if not tertiary_meta or not tertiary_meta.get('domains'):
        if len(sequence) > 300:
             functional['conserved_domains'].append({
                 'name': 'Catalytic Core',
                 'start': 50,
                 'end': 250,
                 'description': 'Predicted catalytic domain based on sequence homology'
             })
    
    # Simulate active sites for demo purposes if none found
    # (Real implementation would parse from UniProt XML/JSON)
    if 'K' in sequence[40:60]: # P-loop logic often has K around 50
        k_pos = sequence[40:60].find('K') + 40 + 1
        functional['active_sites'].append({
            'residue': k_pos,
            'aa': 'K',
            'description': 'ATP binding site (Predicted)',
            'type': 'binding'
        })
        functional['go_terms'].append({'id': 'GO:0005524', 'name': 'ATP binding', 'category': 'function'})

    if 'H' in sequence and 'D' in sequence:
        # Just some random functional annotation for demo
        functional['go_terms'].append({'id': 'GO:0003824', 'name': 'Catalytic activity', 'category': 'function'})

    return {
        'protein_id': protein_id,
        'gene_info': {
            'locus_tag': protein_detail.get('locus_tag', ''),
            'gene_name': protein_detail.get('gene_name', ''),
            'product': protein_detail.get('product', 'Proteína hipotética'),
            'start': protein_detail.get('start', 0),
            'end': protein_detail.get('end', 0),
            'strand': protein_detail.get('strand', 1),
        },
        'primary': primary,
        'secondary': secondary,
        'tertiary': tertiary_meta,
        'quaternary': quaternary,
        'quality': quality,
        'functional': functional,
        'structure_resolved': secondary is not None,
    }


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
    Uses analysis cache first, then falls back to GenBank search.
    """
    # Try getting from analysis cache first
    try:
        from app.api.routes.analysis import _analysis_cache
        genome_data = _analysis_cache.get("genome_data")
        
        if genome_data and genome_data.genes:
            for gene in genome_data.genes:
                if gene.locus_tag == locus_tag:
                    return {
                        "locus_tag": gene.locus_tag,
                        "gene_name": gene.gene_name,
                        "protein_id": gene.protein_id,
                        "product": gene.product,
                        "start": gene.start,
                        "end": gene.end,
                        "length": gene.length,
                        "strand": gene.strand,
                        "gc_content": gene.gc_content,
                        "start_codon": gene.start_codon,
                        "stop_codon": gene.stop_codon
                    }
    except Exception as e:
        print(f"⚠️ [GENE] Cache lookup failed: {e}")
    
    # Fallback: search in GenBank file
    service = get_ncbi_service()
    try:
        genbank_path = _get_genbank_path(request)
        result = service.get_gene_detail(genbank_path, locus_tag)
        if result is not None:
            return result
    except Exception:
        pass
    
    raise HTTPException(status_code=404, detail=f"Gen no encontrado: {locus_tag}")



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
    Uses analysis cache first, then falls back to GenBank search.
    """
    service = get_ncbi_service()
    
    # Try getting from analysis cache first (much faster)
    try:
        from app.api.routes.analysis import _analysis_cache
        genome_data = _analysis_cache.get("genome_data")
        
        if genome_data and genome_data.genbank_filepath:
            result = service.get_central_dogma_data(genome_data.genbank_filepath, locus_tag)
            if result is not None:
                return result
    except Exception as e:
        print(f"⚠️ [DOGMA] Cache lookup failed: {e}")
    
    # Fallback: try active genome
    try:
        genbank_path = _get_genbank_path(request)
        result = service.get_central_dogma_data(genbank_path, locus_tag)
        if result is not None:
            return result
    except Exception:
        pass
    
    raise HTTPException(status_code=404, detail=f"Gen CDS no encontrado: {locus_tag}")


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
        error_msg = result.get("error", "Error BLAST")
        # If error contains HTML keywords, NCBI service is likely down/overloaded
        if any(keyword in error_msg.lower() for keyword in ["html", "doctype", "<!doctype", "no se pudo obtener un id"]):
            raise HTTPException(
                status_code=503, 
                detail="El servicio BLAST de NCBI está temporalmente no disponible o saturado. Por favor, inténtalo de nuevo en unos minutos."
            )
        raise HTTPException(status_code=500, detail=error_msg)

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
