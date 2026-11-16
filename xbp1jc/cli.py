"""Command-line interface for xbp1jc."""
from __future__ import annotations

import argparse
import json
import sys

from .count import (count_bam, count_fastq, N_PLACEMENTS, INTRON_LEN,
                    CHROM, INTRON_START, INTRON_END, SPECIES)
from .stats import wilson, required_reads, adequacy


def _fmt_text(rows: list[dict], baseline: float | None,
              fold: float, show_shifts: bool) -> str:
    out = []
    w = max((len(r["sample"]) for r in rows), default=6)
    head = f"{'sample':<{w}}  {'f(XBP1s)':>9}  {'95% CI':>15}  {'n':>8}  spliced/unspliced"
    out.append(head)
    out.append("-" * len(head))
    for r in rows:
        if r["informative"] == 0:
            out.append(f"{r['sample']:<{w}}  {'n/a':>9}  {'no reads':>15}  "
                       f"{0:>8}  0/0")
            continue
        ci = f"{r['lo']:.1f}-{r['hi']:.1f}"
        out.append(f"{r['sample']:<{w}}  {r['f_xbp1s_pct']:>8.1f}%  {ci:>15}  "
                   f"{r['informative']:>8,}  {r['spliced']}/{r['unspliced']}")
        if r.get("adequacy"):
            a = r["adequacy"]
            flag = "" if a["adequate"] else "   <-- below threshold, report as bound"
            out.append(f"{'':<{w}}  required {a['required']:,} reads for "
                       f"{fold:g}x vs {baseline:g}% baseline: "
                       f"{a['verdict']}{flag}")
        if show_shifts and r.get("shift_distribution"):
            sd = r["shift_distribution"]
            parts = ", ".join(f"+{k}:{v}" for k, v in sorted(sd.items()))
            out.append(f"{'':<{w}}  placements: {parts}")
    return "\n".join(out)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(
        prog="xbp1jc",
        description="Count the XBP1 IRE1alpha-dependent 26-nt splice junction, "
                    "enumerating all equivalent coordinate placements.",
        epilog="f(XBP1s) = spliced / (spliced + unspliced), as a percentage. "
               "Note: sign convention is inverted relative to PSI.")
    ap.add_argument("inputs", nargs="*",
                    help="BAM files, or FASTQ files with --fastq")
    ap.add_argument("--fastq", action="store_true",
                    help="treat inputs as FASTQ and use k-mer mode "
                         "(no reference needed; upward-biased)")
    ap.add_argument("--species", default="human", choices=sorted(SPECIES),
                    help="seam set for FASTQ mode (default: human)")
    ap.add_argument("--policy", default="enumerate",
                    choices=["enumerate", "window", "canonical"],
                    help="coordinate policy in BAM mode (default: enumerate)")
    ap.add_argument("--anchor", type=int, default=8,
                    help="minimum aligned nt each side of the junction "
                         "(default: 8)")
    ap.add_argument("--min-mapq", type=int, default=255,
                    help="minimum MAPQ in BAM mode (default: 255, uniquely "
                         "mapped under STAR)")
    ap.add_argument("--samtools", default="samtools",
                    help="path to samtools (default: samtools on PATH)")
    ap.add_argument("--baseline", type=float, default=None,
                    help="baseline f(XBP1s) percentage against which to assess "
                         "depth adequacy, e.g. 8.3")
    ap.add_argument("--fold", type=float, default=2.0,
                    help="fold change to power for (default: 2.0)")
    ap.add_argument("--power", type=float, default=0.80)
    ap.add_argument("--alpha", type=float, default=0.05)
    ap.add_argument("--json", action="store_true", help="emit JSON")
    ap.add_argument("--show-placements", action="store_true",
                    help="report the distribution of observed placements")
    ap.add_argument("--locus", action="store_true",
                    help="print locus definition and exit")
    args = ap.parse_args(argv)

    if args.locus:
        print(f"GRCh38 {CHROM}:{INTRON_START:,}-{INTRON_END:,} minus strand")
        print(f"intron length      {INTRON_LEN} nt")
        print(f"equivalent placements  {N_PLACEMENTS} (shifts 0..+{N_PLACEMENTS-1})")
        print(f"flanking repeat    GCAGCA (6 nt direct repeat)")
        for name, s in SPECIES.items():
            print(f"\n{name}:")
            print(f"  spliced seam     {s['spliced']}")
            print(f"  unspliced seam   {s['unspliced']}")
            print(f"  intron           {s['intron']}")
        return 0

    if not args.inputs:
        ap.error("no input files given (use --locus to print locus definition)")

    rows = []
    for path in args.inputs:
        try:
            c = (count_fastq([path], species=args.species) if args.fastq
                 else count_bam(path, samtools=args.samtools,
                                anchor=args.anchor, min_mapq=args.min_mapq,
                                policy=args.policy))
        except Exception as e:
            print(f"error on {path}: {e}", file=sys.stderr)
            return 2
        row = c.to_dict()
        row["sample"] = path.rsplit("/", 1)[-1]
        iv = wilson(c.spliced, c.unspliced, alpha=args.alpha)
        if iv:
            row["lo"], row["hi"] = iv.lo, iv.hi
        if args.baseline is not None:
            row["adequacy"] = adequacy(c.spliced, c.unspliced, args.baseline,
                                       args.fold, args.power, args.alpha)
        rows.append(row)

    if args.json:
        print(json.dumps({"locus": f"{CHROM}:{INTRON_START}-{INTRON_END}",
                          "placements_enumerated": N_PLACEMENTS,
                          "results": rows}, indent=2))
    else:
        print(_fmt_text(rows, args.baseline, args.fold, args.show_placements))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
