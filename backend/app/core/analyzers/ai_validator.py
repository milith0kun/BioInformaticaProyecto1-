"""
AI Validator Module
Valida resultados científicos usando Claude AI (Anthropic)

Implementa validación crítica de resultados con IA:
- Compara resultados con conocimiento científico establecido
- Identifica discrepancias o anomalías
- Proporciona interpretación biológica
- Genera recomendaciones

Autor: Sistema de Análisis Genómico
Fecha: Febrero 2026
"""

import os
import json
import re
from typing import Dict, List, Optional, Any
from dataclasses import dataclass
from datetime import datetime
import anthropic


@dataclass
class AIValidation:
    """Resultado de validación con IA"""
    is_valid: bool
    confidence: float
    interpretation: str
    discrepancies: List[str]
    recommendations: List[str]
    scientific_context: str
    key_findings: List[str]
    timestamp: str


class ClaudeValidator:
    """
    Validador científico usando Claude AI (Anthropic)
    """
    
    def __init__(self, api_key: str):
        """
        Inicializar validador
        
        Args:
            api_key: API key de Anthropic Claude
        """
        self.api_key = api_key
        self.client = anthropic.Anthropic(api_key=api_key)
        self.model = "claude-3-5-haiku-20241022"
        
    def validate_codon_analysis(self, results: Dict[str, Any]) -> AIValidation:
        """
        Validar análisis de codones usando IA
        """
        prompt = f"""Eres un bioinformático experto en genómica bacteriana. Analiza estos resultados de E. coli K-12 MG1655:

DATOS:
- ATG (inicio): {results.get('atg_total', 0):,}
- Densidad ATG: {results.get('atg_density', 0):.2f}/kb
- TAA: {results.get('taa_total', 0):,}
- TAG: {results.get('tag_total', 0):,}
- TGA: {results.get('tga_total', 0):,}

REFERENCIA (NC_000913.3):
- ~4,300 genes | 4.6 Mbp | GC: 50.8%

EVALÚA:
1. ¿Los ATG son consistentes? (Debe haber más ATG que genes por ATG internos)
2. ¿Proporción de stop codons normal? (E. coli típico: TAA ~64%, TGA ~30%, TAG ~6%)
3. ¿La densidad de codones es biológicamente razonable?

Responde en JSON:
{{"is_valid":true/false,"confidence":0-100,"interpretation":"resumen científico breve","scientific_context":"contexto biológico","key_findings":["hallazgo1","hallazgo2"],"discrepancies":["si hay"],"recommendations":["recomendaciones"]}}"""
        
        return self._call_and_parse(prompt)
    
    def validate_gene_analysis(self, results: Dict[str, Any]) -> AIValidation:
        """
        Validar análisis de genes usando IA
        """
        prompt = f"""Eres un bioinformático experto. Analiza estos resultados de E. coli K-12 MG1655:

DATOS:
- Genes: {results.get('total_genes', 0):,}
- CDS: {results.get('total_cds', 0):,}
- Genoma: {results.get('genome_length', 0):,} bp
- GC: {results.get('gc_content', 0):.2f}%
- Densidad: {results.get('gene_density', 0):.2f} genes/Mb

REFERENCIA (NC_000913.3):
- Genes: 4,651 | CDS: 4,318 | GC: 50.79% | Densidad: ~1,000 genes/Mb

EVALÚA:
1. Desviación de valores de referencia (aceptable <5%)
2. Consistencia genes/CDS (CDS debe ser ~93% de genes)
3. GC content típico para enterobacterias (45-55%)
4. Cobertura génica del genoma (~87% esperado)

Responde en JSON:
{{"is_valid":true/false,"confidence":0-100,"interpretation":"resumen científico breve","scientific_context":"contexto biológico","key_findings":["hallazgo1","hallazgo2"],"discrepancies":["si hay"],"recommendations":["recomendaciones"]}}"""
        
        return self._call_and_parse(prompt)
    
    def comprehensive_validation(self, full_results: Dict[str, Any]) -> AIValidation:
        """
        Validación comprehensiva de todos los análisis
        """
        codons = full_results.get('codons', {})
        genes = full_results.get('genes', {})
        
        prompt = f"""Eres un bioinformático senior. Realiza validación integral de E. coli K-12 MG1655:

CODONES:
- ATG: {codons.get('atg_total', 0):,} ({codons.get('atg_density', 0):.2f}/kb)
- TAA: {codons.get('taa_total', 0):,} | TAG: {codons.get('tag_total', 0):,} | TGA: {codons.get('tga_total', 0):,}

GENES:
- Total: {genes.get('total_genes', 0):,} | CDS: {genes.get('total_cds', 0):,}
- Genoma: {genes.get('genome_length', 0):,} bp | GC: {genes.get('gc_content', 0):.2f}%

VALIDACIÓN CIENTÍFICA:
1. ¿Los datos son internamente consistentes?
2. ¿Corresponden a E. coli K-12 MG1655 (RefSeq NC_000913.3)?
3. ¿Hay anomalías que requieran investigación?

CONTEXTO:
- E. coli tiene genoma compacto con ~87% codificante
- Ratio ATG/genes indica ATG internos (~15-20x más ATG que genes)
- Stop codons: TAA predomina en genes altamente expresados

Responde en JSON:
{{"is_valid":true/false,"confidence":0-100,"interpretation":"conclusión científica","scientific_context":"relevancia biológica","key_findings":["hallazgo clave 1","hallazgo clave 2","hallazgo clave 3"],"discrepancies":["discrepancia si existe"],"recommendations":["acción recomendada"]}}"""
        
        return self._call_and_parse(prompt)
    
    def _call_and_parse(self, prompt: str) -> AIValidation:
        """
        Llamar a Claude y parsear respuesta
        """
        try:
            message = self.client.messages.create(
                model=self.model,
                max_tokens=2048,
                temperature=0.1,
                messages=[{"role": "user", "content": prompt}]
            )
            
            if message.content and len(message.content) > 0:
                response = message.content[0].text
                return self._parse_ai_response(response)
            
            raise ValueError("Respuesta de Claude vacía")
            
        except anthropic.RateLimitError as e:
            return self._error_response(f"Rate limit excedido. Espere 1-2 minutos.")
        except anthropic.APIError as e:
            return self._error_response(f"Error de API: {str(e)[:100]}")
        except Exception as e:
            return self._error_response(str(e)[:100])
    
    def _parse_ai_response(self, response: str) -> AIValidation:
        """
        Parsear respuesta JSON de la IA con manejo robusto de errores
        """
        # Limpiar respuesta
        response = response.strip()
        
        # Remover bloques de código markdown
        if "```json" in response:
            response = response.split("```json")[-1]
        if "```" in response:
            response = response.split("```")[0]
        response = response.strip()
        
        # Encontrar el JSON válido usando regex
        json_match = re.search(r'\{[^{}]*(?:\{[^{}]*\}[^{}]*)*\}', response, re.DOTALL)
        if json_match:
            response = json_match.group()
        
        # Intentar parsear
        try:
            data = json.loads(response)
        except json.JSONDecodeError:
            # Si falla, intentar arreglar JSON común
            response = re.sub(r',\s*}', '}', response)  # Trailing commas
            response = re.sub(r',\s*]', ']', response)  # Trailing commas in arrays
            try:
                data = json.loads(response)
            except:
                return self._error_response("Error parseando respuesta de IA")
        
        return AIValidation(
            is_valid=data.get("is_valid", False),
            confidence=float(data.get("confidence", 0)),
            interpretation=data.get("interpretation", ""),
            scientific_context=data.get("scientific_context", ""),
            key_findings=data.get("key_findings", []),
            discrepancies=data.get("discrepancies", []),
            recommendations=data.get("recommendations", []),
            timestamp=datetime.now().isoformat()
        )
    
    def _error_response(self, message: str) -> AIValidation:
        """
        Generar respuesta de error
        """
        return AIValidation(
            is_valid=False,
            confidence=0.0,
            interpretation=f"Error: {message}",
            scientific_context="No disponible debido al error",
            key_findings=[],
            discrepancies=["Validación IA no completada"],
            recommendations=["Verificar resultados manualmente", "Reintentar validación"],
            timestamp=datetime.now().isoformat()
        )


# Singleton para reutilizar la instancia
_validator_instance: Optional[ClaudeValidator] = None


def get_ai_validator(api_key: Optional[str] = None) -> Optional[ClaudeValidator]:
    """
    Obtener instancia del validador de IA
    """
    global _validator_instance
    
    if _validator_instance is None:
        key = api_key or os.getenv("CLAUDE_API_KEY")
        if key:
            _validator_instance = ClaudeValidator(key)
    
    return _validator_instance
