"""
xbp1-junction-counter: reliable counting of the XBP1 IRE1alpha-dependent
26-nucleotide splice junction from RNA sequencing data.

Addresses three failure modes that are silent in standard analyses:
  1. coordinate degeneracy from the flanking 6-nt direct repeat
     (seven equivalent placements; canonical-only counting recovers ~33%)
  2. class-asymmetric estimator bias in exact-match counting
  3. a measurement floor set by library composition, not sequencing effort
"""
from .count import (count_bam, count_fastq, Counts, revcomp,
                    CHROM, INTRON_START, INTRON_END, INTRON_LEN,
                    N_PLACEMENTS, SHIFTS, ANCHOR, SPECIES,
                    SPLICED_SEAM, UNSPLICED_SEAM)
from .stats import wilson, required_reads, adequacy, Interval

__version__ = "1.0.0"
__all__ = ["count_bam", "count_fastq", "Counts", "revcomp", "wilson",
           "required_reads", "adequacy", "Interval", "CHROM", "INTRON_START",
           "INTRON_END", "INTRON_LEN", "N_PLACEMENTS", "SHIFTS", "ANCHOR",
           "SPECIES", "SPLICED_SEAM", "UNSPLICED_SEAM"]
