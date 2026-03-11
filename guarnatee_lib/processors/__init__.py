"""
Processing modules for different library types.
"""
from guarnatee_lib.processors.paired_processor import PairedLibraryProcessor
from guarnatee_lib.processors.differential_processor import DifferentialLibraryProcessor

__all__ = [
    'PairedLibraryProcessor',
    'DifferentialLibraryProcessor',
]
