#!/usr/bin/env python3
"""
pairs_to_reads.py -- how many mapped reads a simulated read pair yields.

A simulation that subsamples FASTQ files counts generated read *pairs*; what
reaches the ratio is mapped *reads* on the chromosomes being compared, after
the mapping-quality filter. The two differ by a factor that depends on the
genome, the aligner settings and the filter, so it has to be measured rather
than assumed, and a figure axis labelled in one unit cannot be compared with a
table column in the other.

This reads whichever record of the simulation you still have and reports the
factor, so the two can be put on the same scale.

Two input forms
---------------
A table, if your pipeline recorded the mapped count per replicate. The R
function in this project writes it as column `s`:

    pairs_to_reads.py --table sim_ivar_tab.csv --pairs read_count --mapped s

Or the idxstats files themselves, with the pair count taken from the file
name. Mapped reads are summed over the chromosomes named in the chromosome
file, matching what karyoratio counts:

    pairs_to_reads.py --idxstats 'sim/*.idxstats' \\
                      --chrom-file mmarfoi2.1_chroms.tsv \\
                      --pattern '_n(\\d+)_'

Requires: nothing beyond the standard library.
"""

from __future__ import annotations

import argparse
import csv
import glob
import re
import statistics as st
import sys
from pathlib import Path

__version__ = "1.0.0"


def read_chrom_file(path: Path) -> tuple[set[str], str | None]:
    """Sequence names karyoratio would count: the autosomes and the X."""
    keep: set[str] = set()
    x = None
    for raw in path.read_text().splitlines():
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split()
        if len(parts) < 2:
            continue
        name, role = parts[0], parts[1].lower()
        if role in ("autosome", "a", "auto"):
            keep.add(name)
        elif role == "x":
            keep.add(name)
            x = name
    return keep, x


def from_idxstats(pattern: str, chroms: Path, name_re: str) -> list[tuple[int, int]]:
    """(pairs, mapped) per file, pairs taken from the file name."""
    keep, _ = read_chrom_file(chroms)
    if not keep:
        sys.exit(f"[error] no autosomes or X in {chroms}")

    rx = re.compile(name_re)
    files = sorted(glob.glob(pattern))
    if not files:
        sys.exit(f"[error] no files matched {pattern!r}")

    out: list[tuple[int, int]] = []
    skipped = 0
    for f in files:
        m = rx.search(Path(f).name)
        if not m:
            skipped += 1
            continue
        pairs = int(m.group(1))
        mapped = 0
        for line in Path(f).read_text().splitlines():
            cols = line.split("\t")
            if len(cols) >= 3 and cols[0] in keep:
                mapped += int(cols[2])
        out.append((pairs, mapped))

    print(f"[info] {len(out)} files used, {skipped} skipped "
          f"(name did not match {name_re!r})", file=sys.stderr)
    return out


def from_table(path: Path, pairs_col: str, mapped_col: str) -> list[tuple[int, int]]:
    text = path.read_text()
    delim = "\t" if text.count("\t") > text.count(",") else ","
    rows = list(csv.DictReader(text.splitlines(), delimiter=delim))
    if not rows:
        sys.exit(f"[error] {path} has no rows")
    for col in (pairs_col, mapped_col):
        if col not in rows[0]:
            sys.exit(f"[error] column {col!r} not in {path}; columns are: "
                     + ", ".join(rows[0]))
    out = []
    for r in rows:
        try:
            out.append((int(float(r[pairs_col])), int(float(r[mapped_col]))))
        except (ValueError, TypeError):
            continue
    return out


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(
        description="Measure the generated-pair to mapped-read factor.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    p.add_argument("--table", type=Path,
                   help="CSV or TSV with one row per replicate")
    p.add_argument("--pairs", default="read_count",
                   help="column holding the generated pair count "
                        "(default: read_count)")
    p.add_argument("--mapped", default="s",
                   help="column holding the mapped read count (default: s)")
    p.add_argument("--idxstats",
                   help="glob for samtools idxstats files, quoted")
    p.add_argument("--chrom-file", type=Path,
                   help="karyoratio chromosome file, needed with --idxstats")
    p.add_argument("--pattern", default=r"(\d+)",
                   help="regex whose first group is the pair count in the "
                        "file name (default: the first run of digits)")
    p.add_argument("-o", "--output", type=Path,
                   help="write the per-level table here")
    p.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    args = p.parse_args(argv)

    if args.table:
        data = from_table(args.table, args.pairs, args.mapped)
    elif args.idxstats:
        if not args.chrom_file:
            p.error("--idxstats needs --chrom-file")
        data = from_idxstats(args.idxstats, args.chrom_file, args.pattern)
    else:
        p.error("give either --table or --idxstats")

    if not data:
        sys.exit("[error] nothing to summarise")

    levels: dict[int, list[int]] = {}
    for pairs, mapped in data:
        levels.setdefault(pairs, []).append(mapped)

    out = args.output.open("w") if args.output else sys.stdout
    try:
        print("pairs\tn\tmapped_median\tmapped_mean\tmapped_min\tmapped_max\t"
              "reads_per_pair", file=out)
        for pairs in sorted(levels):
            v = levels[pairs]
            factor = st.median(v) / pairs if pairs else float("nan")
            print(f"{pairs}\t{len(v)}\t{st.median(v):.0f}\t{st.mean(v):.0f}\t"
                  f"{min(v)}\t{max(v)}\t{factor:.4f}", file=out)
    finally:
        if args.output:
            out.close()

    factors = [st.median(v) / pairs for pairs, v in levels.items() if pairs]
    overall = st.median(factors)
    spread = (min(factors), max(factors))

    print(f"\n[result] {overall:.4f} mapped reads per generated pair "
          f"(range {spread[0]:.4f}-{spread[1]:.4f} across levels)",
          file=sys.stderr)
    print(f"         To convert a figure axis in pairs to mapped reads, "
          f"multiply by {overall:.3f};", file=sys.stderr)
    print(f"         to go the other way, divide.", file=sys.stderr)

    if spread[1] / max(spread[0], 1e-9) > 1.15:
        print("\n[warn] the factor varies by more than 15% across read counts, "
              "so a single\n"
              "       multiplier will not do. Convert level by level from the "
              "table above,\n"
              "       or report both units side by side.", file=sys.stderr)
    else:
        print(f"\n[note] the factor is stable across levels, so one multiplier "
              f"serves.", file=sys.stderr)

    return 0


if __name__ == "__main__":
    sys.exit(main())
