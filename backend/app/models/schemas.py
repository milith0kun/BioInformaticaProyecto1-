"""
Pydantic schemas for API request/response validation
"""
from pydantic import BaseModel
from typing import List, Optional, Dict, Any
from enum import Enum


class ValidationStatus(str, Enum):
    PASS = "PASS"
    WARNING = "WARNING"
    FAIL = "FAIL"


class FileInfo(BaseModel):
    filename: str
    filepath: str
    extension: str
    size_bytes: int
    file_type: str
    is_primary: bool = False


class FileListResponse(BaseModel):
    files: List[FileInfo]
    total_count: int
    primary_file: Optional[str] = None


class StopCodonInfo(BaseModel):
    count: int
    percentage: float


class GeneComparison(BaseModel):
    annotated_genes: int
    atg_found: int
    difference: int


class CodonAnalysisResponse(BaseModel):
    genome_length: int
    atg_count: int
    atg_density: float
    stop_codons: Dict[str, StopCodonInfo]
    gene_comparison: GeneComparison


class SizeStatistics(BaseModel):
    mean: float
    median: float
    min: int
    max: int
    std: float


class GeneInfo(BaseModel):
    locus_tag: str
    start: int
    end: int
    length: int
    strand: int
    product: Optional[str] = None
    gc_content: float
    gene_name: Optional[str] = None
    protein_id: Optional[str] = None
    start_codon: Optional[str] = None
    stop_codon: Optional[str] = None
    has_introns: bool = False


class GeneAnalysisResponse(BaseModel):
    total_genes: int
    total_cds: int
    genome_length: int
    gc_content: float
    gene_density: float
    size_statistics: SizeStatistics
    genes: List[GeneInfo]


class ValidationItem(BaseModel):
    metric: str
    calculated: float
    reference: float
    deviation_percent: float
    status: ValidationStatus


class ValidationResponse(BaseModel):
    items: List[ValidationItem]
    overall_status: ValidationStatus


class CompleteAnalysisResponse(BaseModel):
    codons: CodonAnalysisResponse
    genes: GeneAnalysisResponse
    validation: ValidationResponse


class AnalysisStatus(BaseModel):
    status: str
    progress: int
    message: str


class ExportFormat(str, Enum):
    JSON = "json"
    CSV = "csv"
    EXCEL = "excel"
    PDF = "pdf"


class PaginatedGenesResponse(BaseModel):
    genes: List[GeneInfo]
    total: int
    page: int
    page_size: int
    total_pages: int
