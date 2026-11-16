# xbp1-junction-counter

Reliable counting of the IRE1α-dependent 26-nucleotide splice junction in XBP1 mRNA from RNA sequencing data.

XBP1 splicing is the only direct RNA-level readout of unfolded-protein-response activation, which makes it attractive to quantify from RNA-seq. Three properties of the locus cause standard approaches to misreport it, and all three are silent — they produce plausible-looking numbers rather than errors.

| Problem | Consequence if ignored |
|---|---|
| The intron is flanked by a 6-nt direct repeat, so the identical spliced product can sit at **seven** different coordinates | Counting the canonical coordinate alone recovers **33%** of spliced reads; a ±5 nt window recovers **84%** |
| The unspliced seam is longer than the spliced seam, so exact-match counting loses denominator faster than numerator | k-mer counting recovered **95.1%** of spliced but **85.8%** of unspliced reads, inflating the ratio in **18/18** samples |
| Informative reads depend on library *composition*, not sequencing effort | A library can be rich enough for gene-level conclusions and incapable of supporting any junction-level estimate |

This package addresses each explicitly.

## Install

```bash
pip install git+https://github.com/REPLACE_ME/xbp1-junction-counter.git
```

Or from a clone, for development:

```bash
git clone https://github.com/REPLACE_ME/xbp1-junction-counter.git
cd xbp1-junction-counter
pip install -e ".[test]"
pytest
```

No dependencies beyond the Python standard library. BAM mode requires `samtools` on `PATH`; FASTQ mode requires nothing and no reference genome.

## Quick start

Check the locus definition:

```bash
xbp1jc --locus
```

Count from alignments, enumerating all seven placements (the default):

```bash
xbp1jc sample1.bam sample2.bam --baseline 8.3
```

```
sample         f(XBP1s)           95% CI         n  spliced/unspliced
thapsigargin.bam  71.0%        66.9-74.7       517  367/150
               required 247 reads for 2x vs 8.3% baseline: estimate
dmso.bam           1.1%          0.4-3.3       266  3/263
               required 247 reads for 2x vs 8.3% baseline: estimate
```

Count from raw reads with no reference:

```bash
xbp1jc --fastq --species human reads.fastq.gz
```

Machine-readable output:

```bash
xbp1jc sample.bam --baseline 8.3 --json
```

## What is being measured

```
f(XBP1s) = N(spliced) / [ N(spliced) + N(unspliced) ]
```

reported as a percentage, where both counts are **informative junction reads**: reads crossing the junction with at least 8 nt aligned on each side. Reads lying wholly within an exon or wholly within the intron say nothing about splicing status and are excluded from both numerator and denominator.

Two conventions worth noting:

- **The sign is inverted relative to PSI.** Here the measured form is the one from which sequence has been *removed*.
- **`n` always means informative junction reads**, never sequencing depth. The two are related by a factor that varied by more than an order of magnitude across the conditions we examined, which is exactly why the distinction matters.

Because numerator and denominator are the same read class at the same locus in the same library, library size and locus expression cancel. No size factor is applied or needed.

## The seven placements

The intron is flanked by the direct repeat `GCAGCA`. An aligner can therefore place the deletion at any of seven offsets and reconstruct an identical transcript; which one it reports depends on its tie-breaking rule, not on biology.

```
        ...TCCGCAGCA | CTCAGACTACGTGCACCTCTGCAGCA | GGTGCAGGCC...
           ^^^^^^                                   
           repeat                      repeat ^^^^^^
```

Number of equivalent placements = repeat length + 1 = **7**.

Three policies are available so the effect can be measured rather than assumed:

```bash
xbp1jc sample.bam --policy enumerate   # all 7 placements (default, recommended)
xbp1jc sample.bam --policy window      # +/-5 nt tolerance
xbp1jc sample.bam --policy canonical   # exact coordinate only
```

Inspect which placements were actually observed:

```bash
xbp1jc sample.bam --show-placements
```

Any intron flanked by a direct repeat has this property. The details here are XBP1-specific; the rule is not.

## BAM mode or FASTQ mode

**Use BAM mode where alignment is feasible.** It draws both counts from one alignment pass under one filter, so it carries no class asymmetry.

FASTQ mode exists for cases where no alignment is available. It is **upward-biased**, for a structural reason: the unspliced seam (`TCCGCAGCACTCAGACTACG`, 20 nt) is longer than the spliced seam (`TCCGCAGCAGGTGCAGGCC`, 19 nt) and extends further into the intron, so a read needs slightly more favourable positioning to match exactly. Denominator is lost faster than numerator. The bias is systematic, does not average out across replicates, and is largest in the low-coverage samples where it is hardest to notice.

Both modes test the read and its reverse complement, so both are strand-agnostic.

## The depth standard

Supply a baseline and the tool reports whether a sample's read count supports an estimate:

```bash
xbp1jc sample.bam --baseline 8.3 --fold 2.0
```

Against an 8.3% baseline, detecting a doubling at 80% power and α = 0.05 requires **247 informative junction reads per group**. Groups below that should be reported as bounds, not point estimates.

```python
from xbp1jc import required_reads
required_reads(8.3, fold_change=2.0)   # 247
required_reads(30.0, fold_change=2.0)  # fewer: higher baseline is easier
```

The threshold depends on the baseline, so recompute it per application rather than adopting 247 as a universal constant.

## Measurement floor

Where the target locus is a small part of the library, informative reads can be too few regardless of total sequencing. Two mechanisms, same consequence:

- **Viral displacement.** In SARS-CoV-2-infected cells, unique mapping to the human genome fell to 9.5–10.2% against 91–93% in mock. Yield fell to 11.9 informative junction reads per 100M read pairs against 357.3 for a related coronavirus in the same protocol.
- **Failed rRNA depletion.** In 15 mouse lung libraries with rRNA at 64.6–86.7%, median usable reads were 564,860 and the junction returned **zero** informative reads — while interferon-stimulated genes from the same libraries were recovered cleanly (Irf7 34.3×).

The second case is the instructive one: those libraries supported strong gene-level conclusions and no junction-level conclusion at all. **Assess the two separately.**

## Python API

```python
from xbp1jc import count_bam, count_fastq, wilson, required_reads, adequacy

c = count_bam("sample.bam", policy="enumerate")
print(c.spliced, c.unspliced, c.informative, c.f_xbp1s)
print(c.shift_distribution)        # which of the 7 placements were seen

iv = wilson(c.spliced, c.unspliced)
print(iv)                          # 71.0% (95% CI 66.9-74.7; n = 517)

print(adequacy(c.spliced, c.unspliced, baseline_pct=8.3))
# {'required': 247, 'observed': 517, 'adequate': True, 'verdict': 'estimate', ...}
```

`wilson()` returns `None` when there are no informative reads, rather than 0%. An arm with zero informative junction reads has no estimate, and reporting zero would assert something the data cannot support.

## Species

Human and mouse are included. The mouse intron differs from human at a single position (index 9, T vs C), has identical flanks and the same 6-nt repeat, so the same seven placements apply.

```bash
xbp1jc --fastq --species mouse reads.fastq.gz
```

## Reporting checklist

Four items make a value at this locus comparable with others. None requires extra experimental work.

1. Informative junction read count
2. An interval, preferably Wilson
3. The coordinate policy used
4. The estimator used

## Validation

```bash
pytest
```

22 tests, including reproduction of the published control pair — thapsigargin 71.0% [66.9, 74.7] against DMSO 1.1% [0.4, 3.3], intervals required to separate. The effect is far too large to arise by accident, so reproducing it confirms a working installation.

## Citation

If you use this software, please cite the accompanying paper and the archived release. See `CITATION.cff`.

## Licence

MIT. See `LICENSE`.
