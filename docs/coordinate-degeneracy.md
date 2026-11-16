# Why a 26-nt intron needs seven coordinates

## The repeat

The XBP1 intron in GRCh38 occupies `chr22:28,796,122-28,796,147` on the minus
strand. Immediately upstream of the intron, and again at its 5' end, sits the
same six-nucleotide sequence:

```
exon 1 ...TCCGCAGCA | CTCAGACTACGTGCACCTCTGCAGCA | GGTGCAGGCC... exon 2
              GCAGCA                      GCAGCA
              ^^^^^^ repeat               ^^^^^^ repeat
```

Because the sequence removed is identical whether you start the deletion at the
annotated position or up to six nucleotides later, there are seven distinct
coordinate representations of one biological product.

    number of equivalent placements = repeat length + 1 = 7

## Why aligners disagree

Nothing in the data distinguishes the seven. An aligner picks one by its
tie-breaking rule. That rule is deterministic for a given software version and
parameter set, so the resulting undercount is:

- reproducible within a study, and therefore easy to mistake for a real value
- potentially different between studies, and therefore a source of
  cross-study disagreement that looks biological

## Measured effect

Across 67 runs and 4,881 enumerated spliced reads:

| Policy | Recovered | Percentage |
|---|---|---|
| Canonical coordinate only | 1,606 | 32.9% |
| ±5 nt window | 4,104 | 84.1% |
| All seven placements | 4,881 | 100% |

A tolerance window helps but does not solve the problem: the placements span a
range wider than ±5 nt and their distribution is not centred on the annotation.

## Recommendation

Enumerate. The cost is negligible — seven integer comparisons per candidate
read — and it removes the dependency on aligner internals entirely.

```bash
xbp1jc sample.bam --policy enumerate --show-placements
```

Check `--show-placements` on your own data. A distribution concentrated away
from shift 0 means a canonical-only counter would have missed most of your
evidence.

## Generality

This is not an XBP1 quirk. Any intron flanked by a direct repeat admits
multiple equivalent placements, with the count given by the formula above.
XBP1 is simply a case where the repeat is long relative to the intron and the
locus is widely quantified.
