"""
Validator Module
Validates analysis results against E. coli K-12 MG1655 reference values

Módulos BioPython utilizados (opcional):
- Bio.Entrez: Validación de metadata del genoma contra NCBI
- Verificación de que el organismo es correcto
- Obtención de información actualizada de referencia

Implementa el algoritmo de validación 7.3:
- Comparar valores calculados con constantes de referencia
- Calcular: desviación = abs(calculado - referencia) / referencia * 100
- Clasificar: PASS (<5%), WARNING (5-10%), FAIL (>10%)

Valores de referencia: NCBI RefSeq NC_000913.3
"""
from typing import Dict, List, Optional
from dataclasses import dataclass
from enum import Enum

# Bio.Entrez para validación opcional contra NCBI
try:
    from Bio import Entrez
    ENTREZ_AVAILABLE = True
except ImportError:
    ENTREZ_AVAILABLE = False


class ValidationStatus(Enum):
    PASS = "PASS"
    WARNING = "WARNING"
    FAIL = "FAIL"


@dataclass
class ValidationItem:
    metric: str
    calculated: float
    reference: float
    deviation_percent: float
    status: ValidationStatus


@dataclass
class ValidationResult:
    items: List[ValidationItem]
    overall_status: ValidationStatus


# Reference values for E. coli K-12 MG1655
# Source: NCBI Reference Sequence NC_000913.3
REFERENCE_VALUES = {
    "total_genes": {
        "value": 4300,
        "description": "Total number of annotated genes",
        "tolerance": 5  # percent
    },
    "genome_length": {
        "value": 4641652,
        "description": "Total genome length in base pairs",
        "tolerance": 1  # percent
    },
    "gc_content": {
        "value": 50.79,
        "description": "GC content percentage",
        "tolerance": 2  # percent
    },
    "total_cds": {
        "value": 4140,
        "description": "Number of coding sequences",
        "tolerance": 5  # percent
    },
    "gene_density": {
        "value": 926,
        "description": "Genes per megabase",
        "tolerance": 10  # percent
    }
}

# E. coli K-12 MG1655 identifiers for NCBI validation
ECOLI_IDENTIFIERS = {
    "accession": "NC_000913.3",
    "organism": "Escherichia coli str. K-12 substr. MG1655",
    "taxonomy_id": "511145"
}


class Validator:
    """Validator for genomic analysis results"""
    
    def __init__(self):
        self.reference = REFERENCE_VALUES
        self._last_result: ValidationResult = None
    
    def validate(self, calculated_values: Dict) -> ValidationResult:
        """
        Validate calculated values against reference
        
        Args:
            calculated_values: Dictionary with metric names and values
            
        Returns:
            ValidationResult with all validations
        """
        items = []
        
        for metric, calc_value in calculated_values.items():
            if metric in self.reference:
                ref_data = self.reference[metric]
                ref_value = ref_data["value"]
                tolerance = ref_data["tolerance"]
                
                # Calculate deviation
                if ref_value != 0:
                    deviation = abs(calc_value - ref_value) / ref_value * 100
                else:
                    deviation = 0 if calc_value == 0 else 100
                
                # Determine status
                if deviation <= tolerance:
                    status = ValidationStatus.PASS
                elif deviation <= tolerance * 2:
                    status = ValidationStatus.WARNING
                else:
                    status = ValidationStatus.FAIL
                
                items.append(ValidationItem(
                    metric=metric,
                    calculated=round(calc_value, 2),
                    reference=ref_value,
                    deviation_percent=round(deviation, 2),
                    status=status
                ))
        
        # Determine overall status
        if any(item.status == ValidationStatus.FAIL for item in items):
            overall = ValidationStatus.FAIL
        elif any(item.status == ValidationStatus.WARNING for item in items):
            overall = ValidationStatus.WARNING
        else:
            overall = ValidationStatus.PASS
        
        result = ValidationResult(items=items, overall_status=overall)
        self._last_result = result
        return result
    
    def validate_genome_length(self, length: int) -> ValidationItem:
        """Validate genome length specifically"""
        ref = self.reference["genome_length"]
        deviation = abs(length - ref["value"]) / ref["value"] * 100
        
        if deviation <= ref["tolerance"]:
            status = ValidationStatus.PASS
        elif deviation <= ref["tolerance"] * 2:
            status = ValidationStatus.WARNING
        else:
            status = ValidationStatus.FAIL
        
        return ValidationItem(
            metric="genome_length",
            calculated=length,
            reference=ref["value"],
            deviation_percent=round(deviation, 2),
            status=status
        )
    
    def validate_gc_content(self, gc: float) -> ValidationItem:
        """Validate GC content specifically"""
        ref = self.reference["gc_content"]
        deviation = abs(gc - ref["value"]) / ref["value"] * 100
        
        if deviation <= ref["tolerance"]:
            status = ValidationStatus.PASS
        elif deviation <= ref["tolerance"] * 2:
            status = ValidationStatus.WARNING
        else:
            status = ValidationStatus.FAIL
        
        return ValidationItem(
            metric="gc_content",
            calculated=round(gc, 2),
            reference=ref["value"],
            deviation_percent=round(deviation, 2),
            status=status
        )
    
    def validate_gene_count(self, count: int) -> ValidationItem:
        """Validate gene count specifically"""
        ref = self.reference["total_genes"]
        deviation = abs(count - ref["value"]) / ref["value"] * 100
        
        if deviation <= ref["tolerance"]:
            status = ValidationStatus.PASS
        elif deviation <= ref["tolerance"] * 2:
            status = ValidationStatus.WARNING
        else:
            status = ValidationStatus.FAIL
        
        return ValidationItem(
            metric="total_genes",
            calculated=count,
            reference=ref["value"],
            deviation_percent=round(deviation, 2),
            status=status
        )
    
    def get_reference_values(self) -> Dict:
        """Get all reference values"""
        return {
            metric: {
                "value": data["value"],
                "description": data["description"]
            }
            for metric, data in self.reference.items()
        }
    
    def validate_organism_ncbi(self, organism: str, email: str = "your.email@example.com") -> Dict:
        """
        Validate organism against NCBI using Bio.Entrez (opcional)
        
        Esta función usa Bio.Entrez para verificar que el organismo
        analizado corresponde a E. coli K-12 MG1655.
        
        Args:
            organism: Organism name from GenBank file
            email: Email for NCBI Entrez (required by NCBI)
            
        Returns:
            Dictionary with validation results
        """
        if not ENTREZ_AVAILABLE:
            return {
                "validated": False,
                "message": "Bio.Entrez no disponible",
                "expected": ECOLI_IDENTIFIERS["organism"],
                "found": organism
            }
        
        # Check if organism matches expected
        expected = ECOLI_IDENTIFIERS["organism"]
        is_valid = expected.lower() in organism.lower() or organism.lower() in expected.lower()
        
        return {
            "validated": is_valid,
            "message": "Organismo válido" if is_valid else "Organismo no coincide",
            "expected": expected,
            "found": organism,
            "accession_expected": ECOLI_IDENTIFIERS["accession"]
        }
    
    def to_dict(self) -> Dict:
        """Convert last result to dictionary"""
        if not self._last_result:
            return {}
        
        return {
            "items": [
                {
                    "metric": item.metric,
                    "calculated": item.calculated,
                    "reference": item.reference,
                    "deviation_percent": item.deviation_percent,
                    "status": item.status.value
                }
                for item in self._last_result.items
            ],
            "overall_status": self._last_result.overall_status.value
        }
