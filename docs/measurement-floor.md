# The measurement floor

## Statement

Junction-level measurement requires informative reads *at the junction*.
Where the target locus is a small fraction of the library, that count can be
too low to support any estimate, and **additional sequencing does not
necessarily fix it** because the limiting quantity is composition, not depth.

## Demonstration 1: viral displacement

In productive SARS-CoV-2 infection the virus progressively displaces the host
transcriptome. Same cell background, same protocol:

| Condition | Unique mapping to human genome | Informative reads per 100M pairs |
|---|---|---|
| Mock | 91-93% | — |
| HCoV-OC43 | 39-40% | 357.3 |
| SARS-CoV-2 | 9.5-10.2% | 11.9 |

Pooled across replicates the SARS-CoV-2 arm gave 4 spliced and 35 unspliced
reads from 327.7 million read pairs sequenced. The interval is too wide to be
informative.

The uncomfortable feature of this case: the process under study erodes the
material needed to study it. Interest in the measurement peaks exactly where
it becomes impossible.

## Demonstration 2: failed depletion

Fifteen mouse lung libraries, ribosomal RNA at 64.6-86.7%, median 564,860
usable reads per sample. Pooled across all fifteen: 1,088 reads mapped to the
Xbp1 locus, **0** crossed the junction informatively.

The locus was expressed — Xbp1 transcript abundance rose 2.3x on infection —
so this is absence of coverage, not absence of the gene. Positional analysis
excluded isoform-selective dropout: 97.6% of the 1,104-nt uncovered region is
sequence shared by both isoforms, so no isoform-specific mechanism could
produce the gap.

## Why case 2 is the important one

The same libraries yielded clean interferon-stimulated gene induction:

| Gene | Fold change vs mock |
|---|---|
| Irf7 | 34.3x |
| Isg15 | 24.6x |
| Cxcl9 | 12.9x |
| Ifit3b | 8.4x |

So these libraries support strong conclusions about abundant transcripts and
*no* conclusion about one specific junction.

**Gene-level adequacy is not junction-level adequacy.** Assess them
separately. A dataset that passed quality control, aligned well, and produced
a publishable differential-expression result can still be incapable of
answering a junction-level question.

## Making it operational

Fix the threshold before analysis:

```bash
xbp1jc *.bam --baseline 8.3 --fold 2.0
```

Anything flagged `bound` should be reported as a bound with its interval, not
as a point estimate. The alternative — deciding after the fact whether a
number looks reasonable — is how underpowered values enter the literature
indistinguishable from well-supported ones.

## Relationship to prediction-based methods

Methods that predict splicing from gene expression address the same constraint
from the other direction, and are complementary. Where junction reads are
sparse but present, prediction may supply an estimate counting cannot support,
and the depth standard identifies exactly that regime.

Where junction reads are *absent*, as in the fifteen libraries above, there is
nothing at that locus in that material against which a prediction could be
validated. The honest report is that the measurement was not made.
