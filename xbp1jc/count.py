"""
Core junction counting for the XBP1 IRE1alpha-dependent 26-nt intron.

Two counting modes:

  BAM mode (recommended)  - parses CIGAR strings from a splice-aware alignment
                            and enumerates all seven equivalent placements of
                            the deletion created by the flanking direct repeat.
  FASTQ mode (no refs)    - exact-match k-mer counting against the two seams.
                            Requires no reference genome, but is upward-biased:
                            the unspliced seam is longer and extends further
                            into the intron, so an exact-match requirement
                            loses denominator faster than numerator.

The reported quantity is

    f(XBP1s) = N(spliced) / [ N(spliced) + N(unspliced) ]

reported as a percentage. Note the sign convention is INVERTED relative to
percentage-spliced-in (PSI): here the measured form is the one from which
sequence has been removed.

An "informative junction read" crosses the junction with at least ANCHOR
nucleotides aligned on each side. Reads lying wholly within an exon or wholly
within the intron report nothing about splicing status and are excluded from
both numerator and denominator.
"""
from __future__ import annotations

import gzip
import re
import subprocess
from collections import Counter
from dataclasses import dataclass, asdict

# --------------------------------------------------------------------------
# Locus definition, GRCh38
# --------------------------------------------------------------------------
CHROM = "chr22"
INTRON_START = 28796122       # 1-based, first intronic base
INTRON_END = 28796147         # 1-based, last intronic base
INTRON_LEN = INTRON_END - INTRON_START + 1          # 26
STRAND = "-"

# The intron is flanked by a 6-nt direct repeat (GCAGCA). The identical spliced
# product is therefore representable at seven distinct coordinates: the number
# of equivalent placements equals repeat length + 1.
REPEAT = "GCAGCA"
N_PLACEMENTS = len(REPEAT) + 1                       # 7
SHIFTS = tuple(range(N_PLACEMENTS))                  # 0 .. +6

# Minimum nucleotides that must be aligned either side of the junction.
ANCHOR = 8

# --------------------------------------------------------------------------
# Seam sequences (human). Used only in FASTQ mode.
# --------------------------------------------------------------------------
HUMAN_INTRON = "CTCAGACTACGTGCACCTCTGCAGCA"
SPLICED_SEAM = "TCCGCAGCAGGTGCAGGCC"                 # spans exon1|exon2
UNSPLICED_SEAM = "TCCGCAGCACTCAGACTACG"              # spans exon1|intron

# Mouse differs from human at one position only (position 10, T vs C), with
# identical flanks and the same 6-nt repeat, so the same seven placements apply.
MOUSE_INTRON = "CTCAGACTATGTGCACCTCTGCAGCA"
MOUSE_UNSPLICED_SEAM = "TCCGCAGCACTCAGACTATG"

SPECIES = {
    "human": {"intron": HUMAN_INTRON,
              "spliced": SPLICED_SEAM,
              "unspliced": UNSPLICED_SEAM},
    "mouse": {"intron": MOUSE_INTRON,
              "spliced": SPLICED_SEAM,
              "unspliced": MOUSE_UNSPLICED_SEAM},
}

_CIGAR = re.compile(r"(\d+)([MIDNSHP=X])")
_COMP = str.maketrans("ACGTN", "TGCAN")


def revcomp(s: str) -> str:
    return s.translate(_COMP)[::-1]


@dataclass
class Counts:
    """Junction counts for one sample."""
    spliced: int = 0
    unspliced: int = 0
    mode: str = ""
    shift_distribution: dict | None = None

    @property
    def informative(self) -> int:
        return self.spliced + self.unspliced

    @property
    def f_xbp1s(self) -> float | None:
        """Spliced fraction as a percentage, or None if no informative reads."""
        n = self.informative
        return None if n == 0 else 100.0 * self.spliced / n

    def to_dict(self) -> dict:
        d = asdict(self)
        d["informative"] = self.informative
        d["f_xbp1s_pct"] = self.f_xbp1s
        return d


# ==========================================================================
# BAM mode
# ==========================================================================
def count_bam(bam_path: str,
              samtools: str = "samtools",
              chrom: str = CHROM,
              intron_start: int = INTRON_START,
              intron_len: int = INTRON_LEN,
              anchor: int = ANCHOR,
              min_mapq: int = 255,
              policy: str = "enumerate",
              window: int = 5) -> Counts:
    """
    Count informative junction reads from a splice-aware alignment.

    policy:
      "enumerate" - count a deletion of the correct length at ANY of the seven
                    equivalent placements. This is the default and the only
                    policy that recovers all evidence.
      "window"    - count a deletion of the correct length whose start lies
                    within `window` nt of the annotation. Provided for
                    comparison; recovers ~84% in our corpus.
      "canonical" - count only an exact coordinate match. Provided for
                    comparison; recovers ~33% in our corpus.

    Requires samtools on PATH (or pass an explicit path).
    """
    if policy not in ("enumerate", "window", "canonical"):
        raise ValueError(f"unknown policy: {policy}")

    intron_end = intron_start + intron_len - 1
    region = f"{chrom}:{intron_start - 60}-{intron_end + 60}"
    proc = subprocess.run(
        [samtools, "view", "-q", str(min_mapq), "-F", "0x904", bam_path, region],
        capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError(f"samtools failed: {proc.stderr.strip()[:300]}")

    allowed = {intron_start + s for s in SHIFTS}
    shifts: Counter = Counter()
    spliced = unspliced = 0

    for rec in proc.stdout.splitlines():
        f = rec.split("\t")
        if len(f) < 6:
            continue
        ref_pos, cigar = int(f[3]), f[5]
        pos = ref_pos
        left_aligned = 0
        saw_junction = False

        ops = _CIGAR.findall(cigar)
        for i, (n_str, op) in enumerate(ops):
            n = int(n_str)
            if op == "N":
                if n == intron_len and left_aligned >= anchor:
                    # nucleotides aligned to the right of this gap
                    right_aligned = sum(int(m) for m, o in ops[i + 1:]
                                        if o in "M=X")
                    if right_aligned >= anchor:
                        ok = (pos in allowed if policy == "enumerate"
                              else abs(pos - intron_start) <= window
                              if policy == "window"
                              else pos == intron_start)
                        if ok:
                            spliced += 1
                            shifts[pos - intron_start] += 1
                            saw_junction = True
                pos += n
                left_aligned = 0
            elif op in "M=X":
                left_aligned += n
                pos += n
            elif op in "DP":
                pos += n
            # I, S, H consume query only

        if saw_junction:
            continue

        # Unspliced: contiguous alignment spanning an exon-intron boundary with
        # the required anchor either side. Checked at the 5' boundary.
        if _spans_boundary(ref_pos, cigar, intron_start, anchor):
            unspliced += 1

    return Counts(spliced=spliced, unspliced=unspliced,
                  mode=f"bam:{policy}", shift_distribution=dict(shifts))


def _spans_boundary(ref_pos: int, cigar: str, boundary: int,
                    anchor: int) -> bool:
    """
    True if a contiguous aligned block covers `boundary` with at least `anchor`
    aligned nucleotides on each side and no junction gap at the locus.
    """
    pos = ref_pos
    for n_str, op in _CIGAR.findall(cigar):
        n = int(n_str)
        if op in "M=X":
            block_start, block_end = pos, pos + n - 1
            if (block_start <= boundary - anchor
                    and block_end >= boundary + anchor - 1):
                return True
            pos += n
        elif op in "DNP":
            pos += n
    return False


# ==========================================================================
# FASTQ mode
# ==========================================================================
def count_fastq(paths: list[str], species: str = "human") -> Counts:
    """
    Exact-match k-mer counting, no reference required.

    Both the read and its reverse complement are tested, so the function is
    strand-agnostic. A read matching both seams is ambiguous and is discarded.

    NOTE: upward-biased relative to BAM mode. In our corpus this mode recovered
    95.1% of spliced but only 85.8% of unspliced reads, with the ratio above
    unity in 18 of 18 samples. Use BAM mode where alignment is feasible.
    """
    if species not in SPECIES:
        raise ValueError(f"unknown species: {species}; expected one of "
                         f"{sorted(SPECIES)}")
    s = SPECIES[species]
    sp_set = {s["spliced"], revcomp(s["spliced"])}
    un_set = {s["unspliced"], revcomp(s["unspliced"])}

    spliced = unspliced = 0
    for path in paths:
        opener = gzip.open if path.endswith(".gz") else open
        with opener(path, "rt") as fh:
            for i, line in enumerate(fh):
                if i % 4 != 1:
                    continue
                seq = line.strip()
                is_sp = any(sig in seq for sig in sp_set)
                is_un = any(sig in seq for sig in un_set)
                if is_sp and not is_un:
                    spliced += 1
                elif is_un and not is_sp:
                    unspliced += 1
    return Counts(spliced=spliced, unspliced=unspliced,
                  mode=f"fastq:kmer:{species}")
