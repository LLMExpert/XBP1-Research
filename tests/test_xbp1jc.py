"""
Tests for xbp1jc.

The headline test case is the chemical induction control from the published
corpus: thapsigargin 71.0% against DMSO vehicle 1.1%. A user who reproduces
those two numbers from the shipped synthetic reads has a working installation,
because the effect is far too large to arise by accident.
"""
import math

import pytest

from xbp1jc import (count_fastq, wilson, required_reads, adequacy, revcomp,
                    N_PLACEMENTS, INTRON_LEN, SHIFTS, SPECIES,
                    SPLICED_SEAM, UNSPLICED_SEAM)


# ---------------------------------------------------------------- locus facts
def test_intron_length_is_26():
    assert INTRON_LEN == 26


def test_seven_equivalent_placements():
    """Repeat length 6 implies 7 placements. This is the paper's central point."""
    assert N_PLACEMENTS == 7
    assert tuple(SHIFTS) == (0, 1, 2, 3, 4, 5, 6)


def test_mouse_intron_differs_from_human_at_one_position():
    h, m = SPECIES["human"]["intron"], SPECIES["mouse"]["intron"]
    assert len(h) == len(m) == 26
    diffs = [i for i, (a, b) in enumerate(zip(h, m)) if a != b]
    assert diffs == [9], f"expected a single difference at index 9, got {diffs}"


def test_both_species_share_the_spliced_seam():
    assert SPECIES["human"]["spliced"] == SPECIES["mouse"]["spliced"]


def test_unspliced_seam_is_longer_than_spliced():
    """
    This asymmetry is the mechanistic cause of the k-mer estimator bias:
    the longer seam is harder to match exactly, so denominator is lost faster
    than numerator and the ratio is inflated.
    """
    assert len(UNSPLICED_SEAM) > len(SPLICED_SEAM)


def test_revcomp_roundtrip():
    assert revcomp(revcomp(SPLICED_SEAM)) == SPLICED_SEAM
    assert revcomp("ACGT") == "ACGT"


# ---------------------------------------------------------------- Wilson
def test_wilson_returns_none_with_no_reads():
    """Zero informative reads means no estimate, not 0%."""
    assert wilson(0, 0) is None


def test_wilson_contains_point_estimate():
    iv = wilson(71, 29)
    assert iv.lo < iv.point < iv.hi
    assert iv.n == 100
    assert math.isclose(iv.point, 71.0)


def test_wilson_never_negative_at_zero_proportion():
    """The reason Wilson is used rather than the normal approximation."""
    iv = wilson(0, 300)
    assert iv.lo >= 0.0
    assert iv.point == 0.0
    assert iv.hi > 0.0


def test_wilson_narrows_with_depth():
    narrow = wilson(710, 290)
    wide = wilson(71, 29)
    assert (narrow.hi - narrow.lo) < (wide.hi - wide.lo)


def test_wilson_matches_published_thapsigargin_interval():
    """
    Thapsigargin arm from the corpus: 71.0% [66.9, 74.7] on n=517.
    367/517 = 70.98%.
    """
    iv = wilson(367, 150)
    assert iv.n == 517
    assert math.isclose(iv.point, 71.0, abs_tol=0.1)
    assert math.isclose(iv.lo, 66.9, abs_tol=0.4)
    assert math.isclose(iv.hi, 74.7, abs_tol=0.4)


# ---------------------------------------------------------------- power
def test_depth_standard_is_247_at_the_published_baseline():
    """
    The paper's depth standard: 247 informative junction reads per group to
    detect a doubling of the 8.3% pooled mock baseline at 80% power, alpha=0.05.
    """
    assert required_reads(8.3, fold_change=2.0) == 247


def test_higher_baseline_needs_fewer_reads():
    assert required_reads(30.0, 2.0) < required_reads(5.0, 2.0)


def test_larger_effect_needs_fewer_reads():
    assert required_reads(8.3, 4.0) < required_reads(8.3, 2.0)


def test_fold_change_of_one_is_rejected():
    with pytest.raises(ValueError):
        required_reads(8.3, fold_change=1.0)


def test_adequacy_verdicts():
    below = adequacy(10, 78, baseline_pct=8.3)      # n=88, the 36 hpi point
    assert below["observed"] == 88
    assert below["required"] == 247
    assert not below["adequate"]
    assert below["verdict"] == "bound"

    above = adequacy(48, 421, baseline_pct=8.3)     # n=469, the mock point
    assert above["observed"] == 469
    assert above["adequate"]
    assert above["verdict"] == "estimate"


# ---------------------------------------------------------------- FASTQ mode
def _write_fastq(path, reads):
    with open(path, "w") as fh:
        for i, seq in enumerate(reads):
            fh.write(f"@r{i}\n{seq}\n+\n{'I' * len(seq)}\n")


PAD = "GGTCTGCTGAGTCCGCAGCAGGTCTGCTGAGT"


def test_fastq_counts_spliced_and_unspliced(tmp_path):
    p = tmp_path / "x.fastq"
    _write_fastq(p, [PAD + SPLICED_SEAM + PAD] * 7
                 + [PAD + UNSPLICED_SEAM + PAD] * 3)
    c = count_fastq([str(p)])
    assert c.spliced == 7
    assert c.unspliced == 3
    assert c.informative == 10
    assert math.isclose(c.f_xbp1s, 70.0)


def test_fastq_is_strand_agnostic(tmp_path):
    p = tmp_path / "rc.fastq"
    _write_fastq(p, [revcomp(PAD + SPLICED_SEAM + PAD)] * 4)
    assert count_fastq([str(p)]).spliced == 4


def test_fastq_ignores_uninformative_reads(tmp_path):
    """Reads wholly within an exon carry no splicing information."""
    p = tmp_path / "u.fastq"
    _write_fastq(p, ["ACGT" * 20, "TTTT" * 20, PAD])
    c = count_fastq([str(p)])
    assert c.informative == 0
    assert c.f_xbp1s is None


def test_fastq_mouse_seam_not_matched_by_human_set(tmp_path):
    """Species matters: the mouse unspliced seam differs at one base."""
    p = tmp_path / "m.fastq"
    _write_fastq(p, [PAD + SPECIES["mouse"]["unspliced"] + PAD] * 5)
    assert count_fastq([str(p)], species="mouse").unspliced == 5
    assert count_fastq([str(p)], species="human").unspliced == 0


def test_unknown_species_rejected(tmp_path):
    p = tmp_path / "e.fastq"
    _write_fastq(p, [PAD])
    with pytest.raises(ValueError):
        count_fastq([str(p)], species="ferret")


# ------------------------------------------------- the headline test case
def test_thapsigargin_versus_vehicle(tmp_path):
    """
    Synthetic reconstruction of the published control pair.
    Thapsigargin 71.0% [66.9,74.7] n=517 versus DMSO 1.1% [0.4,3.3] n=266.
    The intervals must not overlap; this is the installation check.
    """
    tg = tmp_path / "thapsigargin.fastq"
    _write_fastq(tg, [PAD + SPLICED_SEAM + PAD] * 367
                 + [PAD + UNSPLICED_SEAM + PAD] * 150)
    dm = tmp_path / "dmso.fastq"
    _write_fastq(dm, [PAD + SPLICED_SEAM + PAD] * 3
                 + [PAD + UNSPLICED_SEAM + PAD] * 263)

    t = count_fastq([str(tg)])
    v = count_fastq([str(dm)])
    ti, vi = wilson(t.spliced, t.unspliced), wilson(v.spliced, v.unspliced)

    assert math.isclose(ti.point, 71.0, abs_tol=0.1)
    assert math.isclose(vi.point, 1.1, abs_tol=0.1)
    assert ti.lo > vi.hi, "control intervals must separate"
    assert adequacy(t.spliced, t.unspliced, 8.3)["adequate"]
