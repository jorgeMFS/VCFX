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
    POS: str
    ID: str
    REF: str
    ALT: str
    Allele_Frequency: str


@dataclass
class InfoSummary:
    INFO_Field: str
    Mean: str
    Median: str
    Mode: str


@dataclass
class AlleleBalance:
    CHROM: str
    POS: str
    ID: str
    REF: str
    ALT: str
    Sample: str
    Allele_Balance: str


@dataclass
class ConcordanceRow:
    CHROM: str
    POS: str
    ID: str
    REF: str
    ALT: str
    SAMPLE1_GT: str
    SAMPLE2_GT: str
    Concordance: str


@dataclass
class HWEResult:
    CHROM: str
    POS: str
    ID: str
    REF: str
    ALT: str
    HWE_pvalue: str


@dataclass
class InbreedingCoefficient:
    Sample: str
    InbreedingCoefficient: str


@dataclass
class VariantClassification:
    CHROM: str
    POS: str
    ID: str
    REF: str
    ALT: str
    Classification: str


@dataclass
class CrossSampleConcordanceRow:
    CHROM: str
    POS: str
    ID: str
    REF: str
    ALT: str
    Num_Samples: str
    Unique_Normalized_Genotypes: str
    Concordance_Status: str


@dataclass
class AncestryAssignment:
    Sample: str
    Assigned_Population: str


@dataclass
class DosageRow:
    CHROM: str
    POS: str
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
    POS: str
    PREV_POS: str
    DISTANCE: str


@dataclass
class IndexEntry:
    CHROM: str
    POS: str
    FILE_OFFSET: str
