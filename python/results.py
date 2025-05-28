from dataclasses import dataclass

__all__ = [
    "AlleleFrequency",
    "InfoSummary",
    "AlleleBalance",
    "ConcordanceRow",
    "HWEResult",
    "InbreedingCoefficient",
    "VariantClassification",
    "CrossSampleConcordanceRow",
    "AncestryAssignment",
    "DosageRow",
    "AncestryInference",
    "DistanceRow",
    "IndexEntry",
]


@dataclass
class AlleleFrequency:
    CHROM: str
    POS: int
    ID: str
    REF: str
    ALT: str
    Allele_Frequency: float


@dataclass
class InfoSummary:
    INFO_Field: str
    Mean: float
    Median: float
    Mode: float


@dataclass
class AlleleBalance:
    CHROM: str
    POS: int
    ID: str
    REF: str
    ALT: str
    Sample: str
    Allele_Balance: float


@dataclass
class ConcordanceRow:
    CHROM: str
    POS: int
    ID: str
    REF: str
    ALT: str
    SAMPLE1_GT: str
    SAMPLE2_GT: str
    Concordance: str


@dataclass
class HWEResult:
    CHROM: str
    POS: int
    ID: str
    REF: str
    ALT: str
    HWE_pvalue: float


@dataclass
class InbreedingCoefficient:
    Sample: str
    InbreedingCoefficient: float


@dataclass
class VariantClassification:
    CHROM: str
    POS: int
    ID: str
    REF: str
    ALT: str
    Classification: str


@dataclass
class CrossSampleConcordanceRow:
    CHROM: str
    POS: int
    ID: str
    REF: str
    ALT: str
    Num_Samples: int
    Unique_Normalized_Genotypes: int
    Concordance_Status: str


@dataclass
class AncestryAssignment:
    Sample: str
    Assigned_Population: str


@dataclass
class DosageRow:
    CHROM: str
    POS: int
    ID: str
    REF: str
    ALT: str
    Dosages: str


@dataclass
class AncestryInference:
    Sample: str
    Inferred_Population: str


@dataclass
class DistanceRow:
    CHROM: str
    POS: int
    PREV_POS: int
    DISTANCE: int


@dataclass
class IndexEntry:
    CHROM: str
    POS: int
    FILE_OFFSET: int
