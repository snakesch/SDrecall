"""
Realignment and recall functionality for SDrecall.

This package contains modules for realigning reads in segmental duplication regions
and recalling variants from the realigned data.
"""

from .stat_realign_group_regions import stat_all_RG_region_size
from .realign_per_RG import imap_process_masked_bam
from .cal_edge_NM_values import calculate_NM_distribution_poisson
from .annotate_HP_tag_to_vars import annotate_vcf
from .prepare_masked_align_region import imap_prepare_masked_align_region_per_RG
from .slice_bam_by_cov import split_bam_by_cov

__all__ = [
    'stat_all_RG_region_size',
    'imap_process_masked_bam',
    'calculate_NM_distribution_poisson',
    'annotate_vcf',
    'imap_prepare_masked_align_region_per_RG',
    'split_bam_by_cov'
]
