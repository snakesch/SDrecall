"""
SDrecall - Segmental Duplication aware variant calling

SDrecall is a bioinformatics tool designed to improve variant calling accuracy
in segmental duplication regions by using custom realignment strategies and
statistical filters to distinguish true variants from mapping errors.
"""

# Package version
__version__ = "0.1.0"

# Import key classes and functions for public API
from src.const import SDrecallPaths
from src.log import configure_logger, logger
from prepare_recall_regions import prepare_recall_regions
from realign_and_recall import SDrecall_per_sample

# Define what should be imported with "from SDrecall import *"
__all__ = [
    'SDrecallPaths',
    'configure_logger',
    'logger',
    'prepare_recall_regions',
    'SDrecall_per_sample'
]
