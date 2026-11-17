# Data sources and tools

This repository holds the **code** (the `xbp1jc` junction counter). The **raw reads** are public in NCBI
GEO/SRA under the accessions below; the **run-level counts, per-study XBP1-locus BAMs, and software mirror**
are archived in the Zenodo data package (DOI: *to be assigned on deposit*).

## Datasets (NCBI GEO / SRA)

| GEO | SRA study | BioProject | study | runs used | role |
|-----|-----------|------------|-------|-----------|------|
| GSE103719 | SRP117282 | PRJNA406823 | Woo et al. 2018 | 8 | thapsigargin assay control |
| GSE147507 | SRP253951 | PRJNA615032 | Blanco-Melo et al. 2020 | 2 in reference set (+66 full series) | permissiveness / host genes |
| GSE157490 | SRP280103 | PRJNA661467 | Kim et al. 2021 | 28 | SARS-CoV-2 time course |
| GSE193169 | SRP353826 | PRJNA795179 | Nguyen et al. 2022 | 27 | MERS-CoV, SARS-CoV-2, HCoV-OC43, sgERN1 control |
| GSE254699 | PRJNA1071544 | PRJNA1071544 | Zhou et al. 2026 | 42 | SARS-CoV-2 mutant panel (14 arms) |
| GSE215433 | SRP403967 | — | Kehrer et al. 2023 | 18 | independent ORF6-deletion |

Per-run accession lists are in [`runs/`](runs/): `corpus_65_runs.tsv`, `panel_42_runs.tsv`,
`kehrer_18_runs.tsv`, `blanco_full_series_runs.tsv`. The private, pre-publication study GSE301842 is
embargoed and is **not** included.

Reads are retrieved from SRA by accession, e.g.:

```
prefetch SRR11517747 && fasterq-dump --split-files SRR11517747
```

## Tools and versions

| tool | version | role |
|------|---------|------|
| SRA Toolkit (prefetch, fasterq-dump) | 3.4.1 | read retrieval |
| STAR | 2.7.11b | alignment for the junction counts and locus BAMs |
| samtools | 1.21 / 1.24 | BAM handling |
| xbp1jc (this repo) | 1.0.0 | junction counting, policy `enumerate` (7 equivalent placements of the GCAGCA repeat), Wilson intervals |
| Salmon | 2.8.0 | gene-level TPM (Blanco-Melo host genes) |
| BLAST (blastn) | 2.17.0+ | ARSDA-style read-to-CDS matching for the host-gene TPM cross-check |
| minimap2 | 2.31-r1302 | viral-read counting (not the source of the published junction counts) |
| regtools | 1.0.0 | independent re-parse of the STAR BAM (benchmark) |
| rMATS | 4.4.0 | retained-intron benchmark |

Reference for the STAR benchmark: GRCh38 (no-alt), GENCODE v44.

## Reproducing a count

```
# align a run and count the XBP1 (chr22:28,796,122-28,796,147, minus strand) junction
python -m xbp1jc count --bam aligned.sortedByCoord.out.bam
```

The host-gene TPM in Tables 2-4 uses host-only normalisation (SumN over human genes); the ARSDA/BLAST
cross-check reproduces Prof. Xia's published values (e.g. S16 ATF4 ~1,348, DDIT3 ~493) under host-only
normalisation. The ARSDA per-run outputs are included in the Zenodo package for method-for-method checking.
