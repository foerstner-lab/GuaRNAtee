"""
Constants and configuration values used throughout GuaRNAtee.
Centralizes all magic numbers and string literals.
"""
from typing import List


class GFFColumns:
    """Standard GFF3 column names"""
    SEQID = "seqid"
    SOURCE = "source"
    TYPE = "type"
    START = "start"
    END = "end"
    SCORE = "score"
    STRAND = "strand"
    PHASE = "phase"
    ATTRIBUTES = "attributes"

    ALL = [SEQID, SOURCE, TYPE, START, END, SCORE, STRAND, PHASE, ATTRIBUTES]
    ESSENTIAL = [SEQID, START, END]


class LibraryModes:
    """RNA-seq library preparation types"""
    FULL_LENGTH = "FL"
    PAIRED_1 = "P1"
    PAIRED_2 = "P2"
    END_5_PRIME = "5E"
    END_3_PRIME = "3E"
    DIFF_5_PRIME = "d5E"

    ALL = [FULL_LENGTH, PAIRED_1, PAIRED_2, END_5_PRIME, END_3_PRIME, DIFF_5_PRIME]
    PAIRED_MODES = [PAIRED_1, PAIRED_2]
    END_MODES = [END_5_PRIME, END_3_PRIME, DIFF_5_PRIME]


class LibraryTypes:
    """Library type identifiers for output"""
    FULL_LENGTH = "full_length"
    PAIRED = "Paired"
    DUAL_LIB = "dual_lib"
    DIFFERENTIAL = "differential"
    TSS_PAIR = "pair"


class FileExtensions:
    """File format extensions"""
    GFF = ".gff"
    XLSX = ".xlsx"
    TSV = ".tsv"
    WIG = ".wig"
    FASTA = ".fasta"


class WigAnnotationColumns:
    """Column names for wiggle file annotations"""
    FILE_PATH = "file_path"
    CONDITION = "condition"
    REPLICATE = "replicate"
    STRAND = "strand"
    LIB_MODE = "lib_mode"
    FILE_PATH_5E = "file_path_5e"
    FILE_PATH_3E = "file_path_3e"
    FILE_PATH_D5E = "file_path_d5e"
    FILE_DESC = "file_desc"

    ALL = [FILE_PATH, CONDITION, REPLICATE, STRAND, LIB_MODE]


class StrandNotation:
    """Strand notation conversions"""
    PLUS = "+"
    MINUS = "-"
    FORWARD = "f"
    REVERSE = "r"
    FORWARD_LETTER = "F"
    REVERSE_LETTER = "R"


class AnnotationClasses:
    """RNA candidate classification types"""
    INTERGENIC = "intergenic"
    ORF_INT = "ORF_int"
    OTHERS = "others"

    ALL_GROUPS = [INTERGENIC, ORF_INT, OTHERS]


class ExportColumns:
    """Columns to drop in Excel export"""
    DROP_COLS = [
        "source", "type", "score", "phase", "id",
        "ss_id", "ss_diff_height", "ss_height", "ss_upstream", "ss_downstream",
        "ts_id", "ts_diff_height", "ts_height", "ts_upstream", "ts_downstream"
    ]


class StatisticsColumns:
    """Column names for statistics output"""
    SEQID = "seqid"
    STRAND = "strand"
    ORGANISM = "Organism"
    FILE_DESC = "file_desc"
    TSS_LIB_TYPE = "TSS_lib_type"
    TSS_LIB_WINDOWS_COUNT = "TSS_lib_windows_count"
    TSS_LIB_PEAKS_COUNT = "TSS_lib_peaks_count"
    TTS_LIB_WINDOWS_COUNT = "TTS_lib_windows_count"
    TTS_LIB_PEAKS_COUNT = "TTS_lib_peaks_count"
    PEAKS_CONNECTIONS_COUNT = "peaks_connections_count"


class AnnotationSource:
    """Annotation source identifier"""
    GUARNATEE = "GuaRNAtee"


class AnnotationType:
    """Annotation type identifier"""
    CANDIDATE = "candidate"


class ConfigKeys:
    """Configuration file keys"""
    MIN_HEIGHT = "min_height"
    MIN_LEN = "min_len"
    MAX_LEN = "max_len"
    READ_LENGTH = "read_length"
    MIN_STEP_FACTOR = "min_step_factor"
    MIN_ORF_UD_FRAG_RATIO = "min_orf_ud_frag_ratio"
    MIN_DISTANCE = "min_distance"
    MAX_TSS_LEN = "max_tss_len"
    MAX_TTS_LEN = "max_tts_len"
    MAX_OUTBOUND_TSS_TOLERANCE = "max_outbound_tss_tolerance"
    MAX_OUTBOUND_TTS_TOLERANCE = "max_outbound_tts_tolerance"
    MIN_MFE = "min_mfe"
    DETAILED_OUTPUT = "detailed_output"
    MERGE_SIMILARITY_RATIO = "merge_similarity_ratio"
