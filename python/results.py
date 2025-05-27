from dataclasses import dataclass

__all__ = [
    "AlleleFrequency",
    "InfoSummary",
    "AlleleBalance",
    "ConcordanceRow",
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
