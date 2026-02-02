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
        self.model = "claude-3-5-haiku-20241022"  # Claude 3.5 Haiku
        
    def validate_codon_analysis(self, results: Dict[str, Any]) -> AIValidation:
        """
        Validar análisis de codones usando IA
        
        Args:
            results: Resultados del análisis de codones
            
        Returns:
            Validación con interpretación de IA
        """
        # Construir prompt científico
        prompt = f"""
Eres un experto en bioinformática y genómica bacteriana. Analiza los siguientes resultados del análisis de codones de E. coli K-12 MG1655:

RESULTADOS DEL ANÁLISIS:
- Codones ATG (inicio): {results.get('atg_total', 0):,}
- Densidad ATG: {results.get('atg_density', 0):.2f} por kb
- Codones TAA (terminación): {results.get('taa_total', 0):,}
- Codones TAG (terminación): {results.get('tag_total', 0):,}
- Codones TGA (terminación): {results.get('tga_total', 0):,}
- Total codones terminación: {results.get('total_stop_codons', 0):,}

VALORES DE REFERENCIA ESPERADOS:
- Genes anotados: ~4,300
- Genoma: 4.6 Mbp
- Contenido GC: ~50.8%

TAREAS:
1. Valida si los números son científicamente consistentes
2. Compara la cantidad de ATG con genes esperados (debe haber más ATG que genes)
3. Analiza la proporción de codones de terminación (TAA > TAG > TGA esperado en E. coli)
4. Identifica cualquier anomalía o discrepancia
5. Proporciona interpretación biológica
6. Da un nivel de confianza (0-100%)

Responde SOLO en formato JSON:
{{
    "is_valid": true/false,
    "confidence": número entre 0-100,
    "interpretation": "interpretación detallada",
    "discrepancies": ["lista", "de", "discrepancias"],
    "recommendations": ["lista", "de", "recomendaciones"]
}}
"""
        
        try:
            response = self._call_claude(prompt)
            validation = self._parse_ai_response(response)
            return validation
            
        except Exception as e:
            return self._handle_ai_error(e)
    
    def validate_gene_analysis(self, results: Dict[str, Any]) -> AIValidation:
        """
        Validar análisis de genes usando IA
        
        Args:
            results: Resultados del análisis de genes
            
        Returns:
            Validación con interpretación de IA
        """
        prompt = f"""
Eres un experto en genómica de E. coli. Analiza estos resultados del análisis génico de E. coli K-12 MG1655:

RESULTADOS:
- Total de genes: {results.get('total_genes', 0):,}
- Total CDS: {results.get('total_cds', 0):,}
- Tamaño del genoma: {results.get('genome_length', 0):,} bp
- Contenido GC: {results.get('gc_content', 0):.2f}%
- Densidad génica: {results.get('gene_density', 0):.2f} genes/Mb

VALORES ESPERADOS (NC_000913.3):
- Genes: ~4,651
- CDS: ~4,318
- Genoma: 4,641,652 bp
- GC: ~50.8%
- Densidad: ~1,000 genes/Mb

TAREAS:
1. Compara con valores de referencia (desviación aceptable <5%)
2. Valida consistencia entre genes y CDS
3. Verifica que el contenido GC sea típico para E. coli
4. Evalúa si la densidad génica es normal
5. Identifica problemas o anomalías
6. Da nivel de confianza (0-100%)

Responde SOLO en formato JSON:
{{
    "is_valid": true/false,
    "confidence": número entre 0-100,
    "interpretation": "interpretación científica detallada",
    "discrepancies": ["lista de discrepancias encontradas"],
    "recommendations": ["lista de recomendaciones"]
}}
"""
        
        try:
            response = self._call_claude(prompt)
            validation = self._parse_ai_response(response)
            return validation
            
        except Exception as e:
            return self._handle_ai_error(e)
    
    def comprehensive_validation(self, full_results: Dict[str, Any]) -> AIValidation:
        """
        Validación comprehensiva de todos los análisis
        
        Args:
            full_results: Resultados completos del análisis
            
        Returns:
            Validación global con IA
        """
        prompt = f"""
Eres un bioinformático experto. Realiza una validación comprehensiva del análisis genómico completo de E. coli K-12 MG1655:

RESULTADOS COMPLETOS:
{json.dumps(full_results, indent=2)}

CONTEXTO CIENTÍFICO:
E. coli K-12 MG1655 es el organismo modelo más estudiado en biología molecular. 
RefSeq: NC_000913.3
Genome: 4,641,652 bp
Genes: ~4,651
CDS: ~4,318
GC Content: ~50.8%

TAREAS CRÍTICAS:
1. Evalúa la consistencia global de todos los resultados
2. Verifica que los datos sean científicamente válidos
3. Identifica cualquier resultado inesperado o anómalo
4. Proporciona interpretación biológica comprehensiva
5. Sugiere análisis adicionales si es necesario
6. Nivel de confianza global (0-100%)

Responde SOLO en formato JSON:
{{
    "is_valid": true/false,
    "confidence": número entre 0-100,
    "interpretation": "interpretación científica comprehensiva",
    "discrepancies": ["lista completa de discrepancias"],
    "recommendations": ["lista de recomendaciones científicas"]
}}
"""
        
        try:
            response = self._call_claude(prompt)
            validation = self._parse_ai_response(response)
            return validation
            
        except Exception as e:
            return self._handle_ai_error(e)
    
    def _call_claude(self, prompt: str) -> str:
        """
        Llamar a la API de Claude
        
        Args:
            prompt: Prompt para el modelo
            
        Returns:
            Respuesta del modelo
        """
        try:
            message = self.client.messages.create(
                model=self.model,
                max_tokens=4096,
                temperature=0.2,
                messages=[
                    {
                        "role": "user",
                        "content": prompt
                    }
                ]
            )
            
            # Extraer texto de la respuesta
            if message.content and len(message.content) > 0:
                return message.content[0].text
            
            raise ValueError("Respuesta de Claude vacía")
            
        except anthropic.RateLimitError as e:
            raise Exception(f"Rate limit excedido: {str(e)}")
        except anthropic.APIError as e:
            raise Exception(f"Error de API de Claude: {str(e)}")
        except Exception as e:
            raise Exception(f"Error al llamar a Claude: {str(e)}")
    
    def _handle_ai_error(self, error: Exception) -> AIValidation:
        """
        Manejar errores de la API de IA con mensajes específicos
        
        Args:
            error: Excepción capturada
            
        Returns:
            Validación con mensaje de error apropiado
        """
        error_msg = str(error)
        interpretation = f"Error en validación IA: {error_msg}"
        discrepancies = ["No se pudo conectar con IA"]
        recommendations = ["Verificar manualmente los resultados"]
        
        # Mensajes más específicos según el error
        if "rate limit" in error_msg.lower() or "429" in error_msg:
            interpretation = "⏱️ Límite de peticiones excedido. Intente nuevamente en unos minutos."
            discrepancies = ["API temporalmente no disponible por límite de tasa"]
            recommendations = [
                "Espere 1-2 minutos antes de reintentar",
                "Los resultados del análisis son válidos independientemente de la validación IA"
            ]
        elif "404" in error_msg:
            interpretation = "🔧 Modelo no disponible. Verifique la configuración."
            discrepancies = ["Modelo Claude no encontrado"]
            recommendations = ["Verificar el modelo claude-sonnet-4"]
        elif "authentication" in error_msg.lower() or "401" in error_msg or "403" in error_msg:
            interpretation = "🔑 API key inválida. Verifique su configuración en .env"
            discrepancies = ["Autenticación fallida con Claude"]
            recommendations = ["Verificar CLAUDE_API_KEY en archivo .env"]
            
        return AIValidation(
            is_valid=False,
            confidence=0.0,
            interpretation=interpretation,
            discrepancies=discrepancies,
            recommendations=recommendations,
            timestamp=datetime.now().isoformat()
        )
    
    def _parse_ai_response(self, response: str) -> AIValidation:
        """
        Parsear respuesta JSON de la IA
        
        Args:
            response: Respuesta de texto de la IA
            
        Returns:
            Objeto AIValidation
        """
        # Limpiar respuesta (remover markdown si existe)
        response = response.strip()
        if response.startswith("```json"):
            response = response[7:]
        if response.startswith("```"):
            response = response[3:]
        if response.endswith("```"):
            response = response[:-3]
        response = response.strip()
        
        # Parsear JSON
        data = json.loads(response)
        
        return AIValidation(
            is_valid=data.get("is_valid", False),
            confidence=float(data.get("confidence", 0)),
            interpretation=data.get("interpretation", ""),
            discrepancies=data.get("discrepancies", []),
            recommendations=data.get("recommendations", []),
            timestamp=datetime.now().isoformat()
        )


# Singleton para reutilizar la instancia
_validator_instance: Optional[ClaudeValidator] = None


def get_ai_validator(api_key: Optional[str] = None) -> Optional[ClaudeValidator]:
    """
    Obtener instancia del validador de IA
    
    Args:
        api_key: API key de Claude (opcional si está en env)
        
    Returns:
        Instancia del validador o None si no está configurado
    """
    global _validator_instance
    
    if _validator_instance is None:
        key = api_key or os.getenv("CLAUDE_API_KEY")
        if key:
            _validator_instance = ClaudeValidator(key)
    
    return _validator_instance
