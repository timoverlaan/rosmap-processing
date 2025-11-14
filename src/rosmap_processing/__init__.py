"""ROSMAP data processing pipeline.

A reproducible processing pipeline for ROSMAP Alzheimer's Disease 
single-cell RNA-seq datasets.
"""

__version__ = "0.1.3"
__author__ = "Timo Verlaan"
__email__ = "t.verlaan@tudelft.nl"

# Import key functions for convenience
from .utils.logging import setup_logging, get_logger

__all__ = [
    "__version__",
    "setup_logging",
    "get_logger",
]
