"""
Validator Module - Dynamic Validation System
Validates analysis results dynamically against calculated averages
"""
from typing import Dict, List, Optional
from dataclasses import dataclass
from enum import Enum
import numpy as np

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
    validation_type: str = "dynamic"
    reference_source: str = ""


DYNAMIC_TOLERANCES = {
    "total_genes": {"pass": 15, "warning": 30},
    "genome_length": {"pass": 15, "warning": 30},
    "gc_content": {"pass": 5, "warning": 10},
    "total_cds": {"pass": 15, "warning": 30},
    "gene_density": {"pass": 15, "warning": 30}
}

BACTERIAL_RANGES = {
    "genome_length": {"min": 500000, "max": 15000000, "typical_min": 1000000, "typical_max": 8000000},
    "gc_content": {"min": 20, "max": 80, "typical_min": 30, "typical_max": 70},
    "gene_density": {"min": 400, "max": 1500, "typical_min": 800, "typical_max": 1100},
    "total_genes": {"min": 400, "max": 12000, "typical_min": 1500, "typical_max": 7000},
    "total_cds": {"min": 400, "max": 12000, "typical_min": 1500, "typical_max": 7000}
}


class Validator:
    def __init__(self):
        self._last_result: ValidationResult = None
        self._genome_data_list: List[Dict] = []
        self._reference_values: Dict = {}
    
    def set_genome_data(self, genome_data_list: List[Dict]):
        self._genome_data_list = genome_data_list
        self._calculate_dynamic_reference()
    
    def _calculate_dynamic_reference(self):
        if not self._genome_data_list:
            return
        metrics = ["total_genes", "genome_length", "gc_content", "total_cds", "gene_density"]
        for metric in metrics:
            values = [g.get(metric, 0) for g in self._genome_data_list if metric in g]
            if values:
                self._reference_values[metric] = {
                    "value": float(np.mean(values)),
                    "std": float(np.std(values)) if len(values) > 1 else 0,
                    "min": float(min(values)),
                    "max": float(max(values)),
                    "count": len(values)
                }
    
    def validate(self, calculated_values: Dict) -> ValidationResult:
        items = []
        num_genomes = len(self._genome_data_list)
        
        if num_genomes > 1:
            validation_type = "multi"
            reference_source = f"Promedio de {num_genomes} genomas analizados"
            for metric, calc_value in calculated_values.items():
                if metric in self._reference_values:
                    ref_data = self._reference_values[metric]
                    ref_value = ref_data["value"]
                    tolerances = DYNAMIC_TOLERANCES.get(metric, {"pass": 15, "warning": 30})
                    if ref_value != 0:
                        deviation = abs(calc_value - ref_value) / ref_value * 100
                    else:
                        deviation = 0 if calc_value == 0 else 100
                    if deviation <= tolerances["pass"]:
                        status = ValidationStatus.PASS
                    elif deviation <= tolerances["warning"]:
                        status = ValidationStatus.WARNING
                    else:
                        status = ValidationStatus.FAIL
                    items.append(ValidationItem(
                        metric=metric,
                        calculated=round(calc_value, 2),
                        reference=round(ref_value, 2),
                        deviation_percent=round(deviation, 2),
                        status=status
                    ))
        else:
            validation_type = "single"
            reference_source = "Rangos tipicos bacterianos"
            for metric, calc_value in calculated_values.items():
                if metric in BACTERIAL_RANGES:
                    ranges = BACTERIAL_RANGES[metric]
                    typical_center = (ranges["typical_min"] + ranges["typical_max"]) / 2
                    typical_range = ranges["typical_max"] - ranges["typical_min"]
                    distance_from_center = abs(calc_value - typical_center)
                    deviation = (distance_from_center / (typical_range / 2)) * 100 if typical_range > 0 else 0
                    in_absolute = ranges["min"] <= calc_value <= ranges["max"]
                    in_typical = ranges["typical_min"] <= calc_value <= ranges["typical_max"]
                    if in_typical:
                        status = ValidationStatus.PASS
                        deviation = min(deviation, 25)
                    elif in_absolute:
                        status = ValidationStatus.WARNING
                    else:
                        status = ValidationStatus.FAIL
                    items.append(ValidationItem(
                        metric=metric,
                        calculated=round(calc_value, 2),
                        reference=round(typical_center, 2),
                        deviation_percent=round(min(deviation, 100), 2),
                        status=status
                    ))
        
        if any(item.status == ValidationStatus.FAIL for item in items):
            overall = ValidationStatus.FAIL
        elif any(item.status == ValidationStatus.WARNING for item in items):
            overall = ValidationStatus.WARNING
        else:
            overall = ValidationStatus.PASS
        
        result = ValidationResult(
            items=items,
            overall_status=overall,
            validation_type=validation_type,
            reference_source=reference_source
        )
        self._last_result = result
        return result
    
    def get_reference_values(self) -> Dict:
        if self._reference_values:
            return {
                metric: {"value": data["value"], "std": data.get("std", 0), "source": "dynamic_average"}
                for metric, data in self._reference_values.items()
            }
        else:
            return {
                metric: {"value": (data["typical_min"] + data["typical_max"]) / 2, "source": "bacterial_ranges"}
                for metric, data in BACTERIAL_RANGES.items()
            }
    
    def get_group_statistics(self) -> Dict:
        if not self._reference_values:
            return {}
        return {
            metric: {
                "mean": data["value"],
                "std_dev": data.get("std", 0),
                "min": data.get("min", 0),
                "max": data.get("max", 0),
                "count": data.get("count", 0),
                "coefficient_of_variation": (data.get("std", 0) / data["value"] * 100) if data["value"] != 0 else 0
            }
            for metric, data in self._reference_values.items()
        }
    
    def to_dict(self) -> Dict:
        if not self._last_result:
            return {}
        result = {
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
            "overall_status": self._last_result.overall_status.value,
            "validation_type": self._last_result.validation_type,
            "reference_source": self._last_result.reference_source
        }
        group_stats = self.get_group_statistics()
        if group_stats:
            result["group_statistics"] = group_stats
        return result
