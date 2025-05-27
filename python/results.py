from dataclasses import dataclass

__all__ = ["AlleleFrequency", "InfoSummary"]


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
