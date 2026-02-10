"""
Interactive AI Chat Routes
Provides a conversational interface with a bioinformatics expert AI
that has access to NCBI APIs, genome data, and can cite scientific references.
Enhanced: searches GenBank, NCBI Gene, PubMed for real references.
"""
from fastapi import APIRouter, HTTPException, Request, Query
from typing import Optional, List, Dict, Any
from pydantic import BaseModel
from datetime import datetime
import os
import json
import anthropic

from app.core.ncbi_service import get_ncbi_service

router = APIRouter()


# In-memory chat history (per session)
_chat_sessions: Dict[str, List[Dict]] = {}


class ChatMessage(BaseModel):
    message: str
    session_id: str = "default"
    include_genome_context: bool = True
    genome_accession: Optional[str] = None


class ChatResponse(BaseModel):
    response: str
    references: List[Dict] = []
    ncbi_links: List[Dict] = []
    genome_context_used: bool = False
    timestamp: str = ""
    session_id: str = ""


def _get_genome_context(request: Request, genome_accession: Optional[str] = None) -> str:
    """Build comprehensive genome context from ALL available analysis data
    
    Args:
        request: FastAPI request object
        genome_accession: Optional genome accession to load specific genome data
    """
    # Import getter safely to handle module reloading/cache updates
    try:
        from app.api.routes.analysis import get_analysis_cache
        analysis_cache = get_analysis_cache()
    except ImportError:
        analysis_cache = {}

    context_parts = []
    
    # 1. Active Genome Identification
    active_genome = None
    
    # Priority 1: Use provided genome_accession to load from disk
    if genome_accession:
        try:
            # Try to load genome info from disk
            import os
            import json
            from pathlib import Path
            
            genomes_dir = Path("genomes")
            genome_dir = genomes_dir / genome_accession
            
            # Check if genome directory exists
            if genome_dir.exists():
                # Try to load genome_info.json if exists
                genome_info_path = genome_dir / "genome_info.json"
                if genome_info_path.exists():
                    with open(genome_info_path, 'r') as f:
                        genome_info_data = json.load(f)
                    # Create a simple dict/object from loaded data
                    active_genome = type('obj', (object,), genome_info_data)
                else:
                    # If no genome_info.json, at least we know the accession exists
                    # Try to extract info from GenBank files
                    extracted_dir = genome_dir / "extracted"
                    if extracted_dir.exists():
                        gbff_files = list(extracted_dir.glob("*.gbff"))
                        if gbff_files:
                            # Parse basic info from GenBank file
                            from Bio import SeqIO
                            try:
                                record = next(SeqIO.parse(str(gbff_files[0]), "genbank"))
                                active_genome = type('obj', (object,), {
                                    'accession': genome_accession,
                                    'organism_name': record.annotations.get('organism', 'Unknown'),
                                    'organism_common_name': record.annotations.get('source', 'Unknown'),
                                    'genome_size_mb': round(len(record.seq) / 1e6, 2),
                                    'gc_percent': round((record.seq.count('G') + record.seq.count('C')) / len(record.seq) * 100, 2),
                                    'gene_count': len([f for f in record.features if f.type == 'gene']),
                                    'strain': record.annotations.get('strain', 'Unknown'),
                                    'assembly_level': 'Complete' if 'complete' in record.description.lower() else 'Chromosome'
                                })
                            except Exception:
                                pass
        except Exception as e:
            print(f"Error loading genome from accession {genome_accession}: {e}")
    
    # Priority 2: Check downloader state (most reliable for metadata)
    if not active_genome and hasattr(request.app.state, 'ncbi_downloader') and request.app.state.ncbi_downloader.current_genome:
        active_genome = request.app.state.ncbi_downloader.current_genome
    
    # Priority 3: Fallback: Check file_detector if downloader state was lost (e.g. restart)
    if not active_genome and hasattr(request.app.state, 'file_detector') and request.app.state.file_detector.primary_file:
        # We have a file but no metadata loaded. Try to recover.
        try:
            potential_accession = request.app.state.file_detector.primary_file.filepath.split('/')[-3] # .../accession/extracted/file
            if potential_accession.startswith('GCF_') or potential_accession.startswith('GCA_'):
                 # Try to reload info silently
                 downloader = request.app.state.ncbi_downloader
                 if downloader:
                     active_genome = downloader.get_genome_info(potential_accession)
        except Exception:
            pass

    if active_genome:
        genome_info = (
            f"**GENOMA ACTIVO**:\n"
            f"- Organismo: {getattr(active_genome, 'organism_name', 'N/A')} ({getattr(active_genome, 'organism_common_name', 'N/A')})\n"
            f"- Accession: {getattr(active_genome, 'accession', 'N/A')}\n"
            f"- Tamaño: {getattr(active_genome, 'genome_size_mb', 'N/A')} Mb\n"
            f"- GC%: {getattr(active_genome, 'gc_percent', 'N/A')}%\n"
            f"- Genes Totales: {getattr(active_genome, 'gene_count', 'N/A')}\n"
            f"- Cepa: {getattr(active_genome, 'strain', 'N/A')}\n"
            f"- Nivel Ensamblaje: {getattr(active_genome, 'assembly_level', 'N/A')}\n"
        )
        context_parts.append(genome_info)
    else:
        context_parts.append("⚠️ **ADVERTENCIA**: No se ha detectado un genoma activo. El usuario debe seleccionar un genoma del panel lateral o ejecutar el análisis.")

    # 2. Gene Analysis Summary
    if analysis_cache.get("genes"):
        genes = analysis_cache["genes"]
        # Safe access to attributes using .get or getattr if it's an object
        total_genes = getattr(genes, "total_count", 0) if hasattr(genes, "total_count") else len(genes.get("genes", [])) if isinstance(genes, dict) else 0
        
        gene_context = (
            f"**ANÁLISIS DE GENES**:\n"
            f"- Total Genes Identificados: {total_genes}\n"
        )
        # Add more stats if available in the object structure
        if hasattr(genes, "density"):
             gene_context += f"- Densidad Génica: {genes.density:.2f} genes/Mb\n"
        
        context_parts.append(gene_context)

    # 3. Codon Usage Analysis
    if analysis_cache.get("codons"):
        codons = analysis_cache["codons"]
        # Extract key metrics if available
        codon_context = "**USO DE CODONES (Resumen)**:\n"
        if hasattr(codons, "global_stats"):
            stats = codons.global_stats
            codon_context += (
                f"- GC3s: {stats.get('gc3s', 'N/A')}\n"
                f"- Nc (Effective Number of Codons): {stats.get('nc', 'N/A')}\n"
                f"- CBI (Codon Bias Index): {stats.get('cbi', 'N/A')}\n"
            )
        context_parts.append(codon_context)

    # 4. Genomic Architecture (GC Window)
    if analysis_cache.get("structure"):
        structure = analysis_cache["structure"]
        # Handle dict access for structure (it's a dict in current implementation)
        if isinstance(structure, dict) and "gc_stats" in structure:
            gc_stats = structure["gc_stats"]
            struct_context = (
                f"**ARQUITECTURA GENÓMICA (Ventana Deslizante)**:\n"
                f"- GC Promedio: {gc_stats.get('mean', 0):.2f}%\n"
                f"- GC Mín/Máx: {gc_stats.get('min', 0):.2f}% - {gc_stats.get('max', 0):.2f}%\n"
                f"- Desviación Estándar GC: {gc_stats.get('std', 0):.4f}\n"
            )
            context_parts.append(struct_context)

    # 5. Comparative Phylogeny
    if analysis_cache.get("phylogeny"):
        phylo = analysis_cache["phylogeny"]
        if isinstance(phylo, dict) and "genomes" in phylo:
             phylo_context = (
                 f"**CONTEXTO FILOGENÉTICO**:\n"
                 f"- Genomas comparados: {len(phylo['genomes'])}\n"
                 f"- Método: {phylo.get('method', 'UPGMA')}\n"
                 f"- Grupos cercanos identificados en el árbol.\n"
             )
             context_parts.append(phylo_context)

    return "\n\n".join(context_parts) if context_parts else "No hay datos de análisis disponibles. El usuario debe activar un genoma y ejecutar las herramientas de análisis."


SYSTEM_PROMPT = """Eres el Dr. GenomicAI, un sistema experto en bioinformática y genómica microbiana diseñado para asistir en laboratorios de investigación avanzado.

ROL Y PERSONALIDAD:
- Actúas como un **Bioinformático Senior** especializado en bacteriología y biología molecular.
- Tu tono es **profesional, académico y riguroso**, pero accesible.
- Tienes acceso directo a la base de conocimientos de NCBI (GenBank, RefSeq, PubMed).
- Eres capaz de correlacionar datos genómicos con fenotipos biológicos (patogenicidad, metabolismo, resistencia).

BASE DE CONOCIMIENTO (Contexto Crítico):
- Tienes acceso al **GENOMA ACTIVO** del usuario (ver abajo). Úsalo para dar respuestas específicas.
- Entiendes métricas avanzadas: CAI (Codon Adaptation Index), RSCU, Nc (Effective Number of Codons), GC Skew (origen de replicación).
- Conoces los mecanismos de CRISPR-Cas, operones, regulones y factores sigma.
- Manejas taxonomía bacteriana y evolución molecular.

REGLAS DE INTERACCIÓN:
1. **Evidencia Científica**: Basa tus afirmaciones en principios biológicos sólidos. Si mencionas un gen, describe su función.
2. **Contextualización**: Si el usuario pregunta "¿Qué significa este GC?", responde usando el GC% real del genoma activo (ej. "El 50.8% de E. coli K-12 indica...").
3. **Citas NCBI**: Cuando menciones genes o proteínas, asume que existen en NCBI. Sugiere buscar "gen X en NCBI".
4. **Respuesta Estructurada**:
   - Usa **Markdown** rico (tablas para comparaciones, negritas para énfasis).
   - Divide explicaciones complejas en pasos lógicos.
   - Usa emojis científicos (🧬, 🦠, 🧪, 📊) para mejorar la legibilidad visual.
5. **Proactividad**: Sugiere qué análisis realizar a continuación (ej. "Dado este alto contenido GC, deberíamos revisar el uso de codones...").

LIMITACIONES:
- Si no tienes datos del genoma activo, pídele al usuario que cargue uno o ejecute el análisis.
- No Inventes referencias bibliográficas específicas (título/año) a menos que sean papers clásicos muy conocidos (ej. Watson & Crick 1953, Wright 1990).

IDIOMA:
- Responde siempre en **ESPAÑOL CIENTÍFICO** estándar.
"""


@router.post("/message")
async def send_chat_message(chat_msg: ChatMessage, request: Request):
    """
    Send a message to the AI bioinformatics expert.
    Searches PubMed, NCBI Gene, and GenBank for references.
    """
    api_key = os.getenv("CLAUDE_API_KEY")
    if not api_key:
        raise HTTPException(
            status_code=400,
            detail="API key de Claude no configurada. Set CLAUDE_API_KEY en .env"
        )

    session_id = chat_msg.session_id

    # Initialize session if needed
    if session_id not in _chat_sessions:
        _chat_sessions[session_id] = []

    # Build context
    genome_context = ""
    if chat_msg.include_genome_context:
        try:
            genome_context = _get_genome_context(request, chat_msg.genome_accession)
        except Exception as e:
            print(f"Error getting genome context: {e}")
            genome_context = "No hay datos de genoma disponibles actualmente."

    # Search NCBI for relevant references
    references = []
    ncbi_links = []
    ncbi_service = get_ncbi_service()

    # Auto-search literature for relevant queries
    search_keywords = _extract_search_keywords(chat_msg.message)
    if search_keywords:
        try:
            # PubMed search
            refs = ncbi_service.search_literature(search_keywords, max_results=3)
            references = refs
        except Exception:
            pass

        try:
            # NCBI Gene search
            organism = ""
            if hasattr(request.app.state, 'ncbi_downloader') and request.app.state.ncbi_downloader.current_genome:
                organism = request.app.state.ncbi_downloader.current_genome.organism_name
            gene_results = ncbi_service.search_ncbi_gene(search_keywords, organism=organism, max_results=3)
            for g in gene_results:
                ncbi_links.append({
                    "type": "gene",
                    "name": g["name"],
                    "description": g["description"],
                    "url": g["url"],
                    "organism": g["organism"]
                })
        except Exception:
            pass

        try:
            # NCBI Nucleotide search
            nuc_results = ncbi_service.search_ncbi_nucleotide(search_keywords, max_results=2)
            for n in nuc_results:
                ncbi_links.append({
                    "type": "nucleotide",
                    "name": n["accession"],
                    "description": n["title"][:100],
                    "url": n["url"],
                    "genbank_url": n["genbank_url"]
                })
        except Exception:
            pass

    # Build messages for Claude
    messages = []

    # Add conversation history (last 10 messages)
    for msg in _chat_sessions[session_id][-10:]:
        messages.append(msg)

    # Add user message with context
    user_content = chat_msg.message
    if genome_context:
        user_content = f"""[CONTEXTO DEL GENOMA ACTIVO - DATOS REALES DEL ANÁLISIS]
{genome_context}

[PREGUNTA DEL USUARIO]
{chat_msg.message}"""

    if references:
        refs_text = "\n[REFERENCIAS ENCONTRADAS EN PUBMED - PUEDES CITARLAS]\n"
        for ref in references:
            refs_text += f"- \"{ref.get('title', 'Sin título')}\" ({ref.get('year', 'N/A')}) PMID:{ref.get('pmid', '')} - {ref.get('journal', '')} - URL: {ref.get('url', '')}\n"
        user_content += refs_text

    if ncbi_links:
        links_text = "\n[ENLACES NCBI ENCONTRADOS - PUEDES REFERENCIAR]\n"
        for link in ncbi_links:
            links_text += f"- [{link['type'].upper()}] {link['name']}: {link['description']} - {link['url']}\n"
        user_content += links_text

    messages.append({"role": "user", "content": user_content})

    try:
        client = anthropic.Anthropic(api_key=api_key)
        response = client.messages.create(
            model="claude-3-5-haiku-20241022",
            max_tokens=4096,
            temperature=0.3,
            system=SYSTEM_PROMPT,
            messages=messages
        )

        assistant_response = response.content[0].text if response.content else "Sin respuesta"

        # Save to history
        _chat_sessions[session_id].append({"role": "user", "content": chat_msg.message})
        _chat_sessions[session_id].append({"role": "assistant", "content": assistant_response})

        # Keep history manageable
        if len(_chat_sessions[session_id]) > 40:
            _chat_sessions[session_id] = _chat_sessions[session_id][-20:]

        return {
            "response": assistant_response,
            "references": references,
            "ncbi_links": ncbi_links,
            "genome_context_used": bool(genome_context),
            "timestamp": datetime.now().isoformat(),
            "session_id": session_id
        }

    except anthropic.RateLimitError:
        raise HTTPException(status_code=429, detail="Rate limit excedido. Espere 1-2 minutos.")
    except anthropic.APIError as e:
        raise HTTPException(status_code=500, detail=f"Error de API: {str(e)[:200]}")
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Error: {str(e)[:200]}")


@router.get("/history/{session_id}")
async def get_chat_history(session_id: str = "default"):
    """Get chat history for a session"""
    history = _chat_sessions.get(session_id, [])
    return {
        "session_id": session_id,
        "messages": history,
        "count": len(history)
    }


@router.delete("/history/{session_id}")
async def clear_chat_history(session_id: str = "default"):
    """Clear chat history for a session"""
    if session_id in _chat_sessions:
        del _chat_sessions[session_id]
    return {"message": "Historial eliminado", "session_id": session_id}


@router.get("/suggestions")
async def get_chat_suggestions(request: Request):
    """Get contextual suggestions based on current analysis state"""
    try:
        from app.api.routes.analysis import get_analysis_cache
        analysis_cache = get_analysis_cache()

        suggestions = [
            {"text": "¿Qué puedes decirme sobre este genoma?", "icon": "🧬"},
            {"text": "Explica el dogma central: DNA → mRNA → Proteína", "icon": "🔬"},
            {"text": "¿Qué significa el contenido GC y por qué es importante?", "icon": "📊"},
        ]

        if analysis_cache.get("genes"):
            genes = analysis_cache["genes"]
            # Safe access to density
            density = 0
            if isinstance(genes, dict):
                 density = genes.get("gene_density")
            else:
                 density = getattr(genes, "gene_density", None)
            
            if density is None:
                density = 0
                 
            suggestions.extend([
                {"text": f"Analiza la densidad génica de {float(density):.2f} genes/Mb", "icon": "📈"},
                {"text": "¿Cuáles son los genes más largos y qué función tienen?", "icon": "🧪"},
                {"text": "Busca información sobre el gen dnaA en NCBI", "icon": "🔗"},
                {"text": "¿Cómo se compara este GC% con otros organismos similares?", "icon": "🔍"},
            ])

        if analysis_cache.get("codons"):
            suggestions.extend([
                {"text": "¿Por qué TAA es el codón de parada más frecuente en bacterias?", "icon": "🛑"},
                {"text": "Explica qué es el RSCU y cómo se interpreta", "icon": "📋"},
                {"text": "¿Qué indica el Nc (codones efectivos) sobre la expresión génica?", "icon": "📊"},
                {"text": "Explica las hebras 5'→3' y 3'→5' del DNA", "icon": "🧬"},
            ])

        if analysis_cache.get("compared_genomes") and len(analysis_cache["compared_genomes"]) > 1:
            n = len(analysis_cache["compared_genomes"])
            suggestions.extend([
                {"text": f"Compara los {n} genomas e identifica diferencias evolutivas", "icon": "⚖️"},
                {"text": "¿Qué implicaciones tienen las diferencias en GC% entre los genomas?", "icon": "🌳"},
                {"text": "Analiza las diferencias en densidad génica entre los genomas comparados", "icon": "📊"},
            ])

        return {"suggestions": suggestions}
    except Exception as e:
        print(f"Error generating suggestions: {e}")
        # Return basic suggestions on error instead of 500
        return {"suggestions": [
            {"text": "¿Qué puedes decirme sobre este genoma?", "icon": "🧬"},
            {"text": "Explica el dogma central de la biología molecular", "icon": "🔬"},
        ]}


def _extract_search_keywords(message: str) -> str:
    """Extract potential search keywords from a message for NCBI search"""
    # Keywords that suggest scientific topics worth searching
    scientific_terms = [
        'gen ', 'genes', 'proteína', 'protein', 'codon', 'codón',
        'mutación', 'mutation', 'resistencia', 'resistance',
        'virulencia', 'virulence', 'metabolismo', 'metabolism',
        'regulación', 'regulation', 'expresión', 'expression',
        'evolución', 'evolution', 'patogénesis', 'pathogenesis',
        'antibiótico', 'antibiotic', 'plasmido', 'plasmid',
        'crispr', 'recombinación', 'transcripción', 'replicación',
        'promotor', 'operón', 'ribosom', 'trna', 'rrna',
        'traducción', 'translation', 'transcription',
        'dogma', 'mrna', 'dna', 'arn', 'adn',
        'replicación', 'genbank', 'ncbi', 'pubmed',
    ]

    message_lower = message.lower()
    has_scientific = any(term in message_lower for term in scientific_terms)

    if has_scientific:
        # Build a search query
        words = message.split()
        important_words = [w for w in words if len(w) > 3 and w.lower() not in (
            'sobre', 'cómo', 'cuáles', 'cuales', 'porque', 'puedes',
            'explica', 'están', 'tiene', 'tiene', 'hacer', 'puede',
            'dime', 'qué', 'como', 'cuál', 'para', 'está', 'este',
            'esta', 'esos', 'esas', 'entre', 'cada', 'pero',
        )]
        if important_words:
            return " ".join(important_words[:5])

    return ""
