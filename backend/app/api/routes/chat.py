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


class ChatResponse(BaseModel):
    response: str
    references: List[Dict] = []
    ncbi_links: List[Dict] = []
    genome_context_used: bool = False
    timestamp: str = ""
    session_id: str = ""


def _get_genome_context(request: Request) -> str:
    """Build comprehensive genome context from ALL available analysis data"""
    from app.api.routes.analysis import _analysis_cache

    context_parts = []

    # Current genome info
    if hasattr(request.app.state, 'ncbi_downloader') and request.app.state.ncbi_downloader.current_genome:
        genome = request.app.state.ncbi_downloader.current_genome
        context_parts.append(f"""
GENOMA ACTIVO:
- Accession: {genome.accession}
- Organismo: {genome.organism_name}
- Cepa: {genome.strain}
- Tamaño: {genome.genome_size:,} bp
- Genes: {genome.gene_count}
- GC%: {genome.gc_percent}%
- Nivel ensamblaje: {genome.assembly_level}
- Categoría RefSeq: {genome.refseq_category}
- URL GenBank: https://www.ncbi.nlm.nih.gov/datasets/genome/{genome.accession}/
- URL Nucleotide: https://www.ncbi.nlm.nih.gov/nuccore/{genome.accession}
""")

    # Analysis results
    if _analysis_cache.get("genes"):
        genes = _analysis_cache["genes"]
        context_parts.append(f"""
ANÁLISIS DE GENES:
- Total genes: {genes.total_genes}
- Total CDS: {genes.total_cds}
- Longitud genoma: {genes.genome_length:,} bp
- GC content: {genes.gc_content}%
- Densidad génica: {genes.gene_density} genes/Mb
- Tamaño promedio gen: {genes.size_statistics.mean:.0f} bp
- Desviación estándar: {genes.size_statistics.std:.0f} bp
- Gen más largo: {genes.size_statistics.max:.0f} bp
- Gen más corto: {genes.size_statistics.min:.0f} bp
""")

    if _analysis_cache.get("codons"):
        codons = _analysis_cache["codons"]
        stop_info = ""
        try:
            stop_codons = codons.stop_codons if hasattr(codons, 'stop_codons') else {}
            for name, data in stop_codons.items():
                count = data.count if hasattr(data, 'count') else data.get('count', 'N/A')
                stop_info += f"  - {name}: {count}\n"
        except Exception:
            stop_info = "  - Datos no disponibles\n"

        context_parts.append(f"""
ANÁLISIS DE CODONES:
- Codones ATG: {codons.atg_count:,}
- Densidad ATG: {codons.atg_density}/kb
- Stop codons:
{stop_info}
- Genes anotados: {getattr(codons, 'annotated_genes', 'N/A')}
""")

    # Comparison data - ALL genomes
    if _analysis_cache.get("compared_genomes"):
        genomes = _analysis_cache["compared_genomes"]
        comparison_text = f"GENOMAS COMPARADOS ({len(genomes)} genomas):\n"
        for i, g in enumerate(genomes, 1):
            comparison_text += f"""
Genoma #{i}:
  - Accession: {g.get('accession', 'N/A')}
  - Organismo: {g.get('organism_name', 'N/A')}
  - Tamaño: {g.get('genome_length', 0):,} bp
  - GC content: {g.get('gc_content', 0)}%
  - Total genes: {g.get('total_genes', 0)}
  - Densidad génica: {g.get('gene_density', 0)} genes/Mb
  - URL: https://www.ncbi.nlm.nih.gov/datasets/genome/{g.get('accession', '')}/
"""
        context_parts.append(comparison_text)

    # Try to include codon usage summary
    try:
        genbank_path = None
        if hasattr(request.app.state, 'file_detector'):
            genbank_file = request.app.state.file_detector.get_genbank_file()
            if genbank_file:
                genbank_path = genbank_file.filepath

        if genbank_path:
            service = get_ncbi_service()
            codon_data = service.calculate_complete_codon_usage(genbank_path)
            if "error" not in codon_data:
                context_parts.append(f"""
USO DE CODONES COMPLETO:
- Total codones analizados: {codon_data.get('total_codons', 0):,}
- GC3 Content (3ra posición): {codon_data.get('gc3_content', 0)}%
- Nc (Codones Efectivos, Wright 1990): {codon_data.get('effective_number_of_codons', 0)}
- Top 5 codones más usados: {', '.join(f"{c['codon']}({c['count']})" for c in sorted(codon_data.get('codon_table', []), key=lambda x: -x['count'])[:5])}
- Top 5 RSCU más altos: {', '.join(f"{c['codon']}({c['rscu']})" for c in sorted(codon_data.get('codon_table', []), key=lambda x: -x['rscu'])[:5])}
""")
    except Exception:
        pass

    return "\n".join(context_parts) if context_parts else "No hay genoma activo ni análisis ejecutado."


SYSTEM_PROMPT = """Eres un biólogo molecular y bioinformático experto con PhD en genómica microbiana. Tu nombre es GenomicAI.

PERSONALIDAD:
- Eres entusiasta pero riguroso científicamente
- Explicas conceptos complejos de forma accesible
- SIEMPRE citas fuentes y referencias cuando es posible
- Usas terminología técnica pero la explicas
- Eres proactivo sugiriendo análisis adicionales

CAPACIDADES:
- Conoces en profundidad la genómica bacteriana
- Puedes interpretar datos de uso de codones, contenido GC, RSCU, Nc, densidad génica
- Conoces las bases de datos NCBI (GenBank, RefSeq, PubMed, Gene, Protein)
- Puedes comparar genomas y explicar significancia biológica
- Entiendes evolución molecular, regulación génica, metabolismo
- Puedes explicar el dogma central (DNA→mRNA→Proteína), transcripción, traducción
- Entiendes las hebras 5'→3' y 3'→5', codones, anticodones

INSTRUCCIONES CRÍTICAS:
1. Responde siempre en español
2. Usa formato Markdown para estructurar respuestas
3. **SIEMPRE** incluye referencias y enlaces cuando menciones datos científicos:
   - Para genes: [nombre_gen en NCBI](https://www.ncbi.nlm.nih.gov/gene/?term=GENE_NAME)
   - Para proteínas: [protein_id en NCBI](https://www.ncbi.nlm.nih.gov/protein/PROTEIN_ID)
   - Para genomas: [accession en GenBank](https://www.ncbi.nlm.nih.gov/nuccore/ACCESSION)
   - Para papers: incluye PMID y enlace a PubMed
4. Si el usuario pregunta sobre los datos del genoma activo, usa el contexto proporcionado
5. Sugiere análisis adicionales relevantes
6. Si no tienes certeza, indícalo honestamente
7. Para preguntas sobre genes específicos, proporciona contexto funcional
8. Usa emojis relevantes (🧬🔬📊🧪) con moderación
9. Cuando hables de procesos biológicos, especifica la dirección (5'→3', N-terminal→C-terminal)
10. NUNCA inventes PMIDs o referencias. Solo usa las que se proporcionan en el contexto.

FORMATO DE RESPUESTA:
- Usa **negrita** para términos importantes
- Usa listas con viñetas para comparaciones
- Incluye datos numéricos cuando sean relevantes
- Si hay referencias, cítalas con formato: [Autor et al., Año](PMID: XXXXX)
- Incluye enlaces a NCBI cuando mencionas genes, proteínas o genomas"""


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
            genome_context = _get_genome_context(request)
        except Exception:
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
    from app.api.routes.analysis import _analysis_cache

    suggestions = [
        {"text": "¿Qué puedes decirme sobre este genoma?", "icon": "🧬"},
        {"text": "Explica el dogma central: DNA → mRNA → Proteína", "icon": "🔬"},
        {"text": "¿Qué significa el contenido GC y por qué es importante?", "icon": "📊"},
    ]

    if _analysis_cache.get("genes"):
        genes = _analysis_cache["genes"]
        suggestions.extend([
            {"text": f"Analiza la densidad génica de {genes.gene_density} genes/Mb", "icon": "📈"},
            {"text": "¿Cuáles son los genes más largos y qué función tienen?", "icon": "🧪"},
            {"text": "Busca información sobre el gen dnaA en NCBI", "icon": "🔗"},
            {"text": "¿Cómo se compara este GC% con otros organismos similares?", "icon": "🔍"},
        ])

    if _analysis_cache.get("codons"):
        suggestions.extend([
            {"text": "¿Por qué TAA es el codón de parada más frecuente en bacterias?", "icon": "🛑"},
            {"text": "Explica qué es el RSCU y cómo se interpreta", "icon": "📋"},
            {"text": "¿Qué indica el Nc (codones efectivos) sobre la expresión génica?", "icon": "📊"},
            {"text": "Explica las hebras 5'→3' y 3'→5' del DNA", "icon": "🧬"},
        ])

    if _analysis_cache.get("compared_genomes") and len(_analysis_cache["compared_genomes"]) > 1:
        n = len(_analysis_cache["compared_genomes"])
        suggestions.extend([
            {"text": f"Compara los {n} genomas e identifica diferencias evolutivas", "icon": "⚖️"},
            {"text": "¿Qué implicaciones tienen las diferencias en GC% entre los genomas?", "icon": "🌳"},
            {"text": "Analiza las diferencias en densidad génica entre los genomas comparados", "icon": "📊"},
        ])

    return {"suggestions": suggestions}


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
