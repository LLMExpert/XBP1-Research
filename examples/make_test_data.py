#!/usr/bin/env python3
"""
Generate the shipped test case: synthetic reads reproducing the published
thapsigargin / DMSO control pair.

    thapsigargin  71.0% [66.9, 74.7]  n = 517   (367 spliced / 150 unspliced)
    DMSO vehicle   1.1% [0.4, 3.3]    n = 266   (3 spliced / 263 unspliced)

The two intervals do not overlap. Reproducing them confirms a working
installation.

Usage:  python examples/make_test_data.py [outdir]
        xbp1jc --fastq --baseline 8.3 <outdir>/*.fastq
"""
import sys
from pathlib import Path

from xbp1jc import SPLICED_SEAM, UNSPLICED_SEAM

FLANK = "GGTCTGCTGAGTCCGCAGCAGGTCTGCTGAGT"


def write(path: Path, n_spliced: int, n_unspliced: int) -> None:
    with path.open("w") as fh:
        i = 0
        for seam, count in ((SPLICED_SEAM, n_spliced),
                            (UNSPLICED_SEAM, n_unspliced)):
            read = FLANK + seam + FLANK
            for _ in range(count):
                fh.write(f"@read{i}\n{read}\n+\n{'I' * len(read)}\n")
                i += 1
    print(f"{path}: {n_spliced} spliced + {n_unspliced} unspliced "
          f"= {n_spliced + n_unspliced} informative reads")


def main() -> None:
    out = Path(sys.argv[1] if len(sys.argv) > 1 else ".")
    out.mkdir(parents=True, exist_ok=True)
    write(out / "thapsigargin.fastq", 367, 150)
    write(out / "dmso.fastq", 3, 263)
    print(f"\nnow run:  xbp1jc --fastq --baseline 8.3 {out}/thapsigargin.fastq "
          f"{out}/dmso.fastq")


if __name__ == "__main__":
    main()
