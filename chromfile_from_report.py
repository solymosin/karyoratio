#!/usr/bin/env python3
"""
chromfile_from_report.py -- write a karyoratio chromosome definition file.

Takes the roles of sequences from an NCBI sequence report, so nothing is
inferred from name patterns or length cutoffs, and matches them to whatever
naming the aligner index or BAM actually uses: RefSeq accessions, GenBank
accessions, UCSC-style names or bare chromosome numbers. Getting this wrong is
the usual way a prebuilt index causes trouble, because the mismatch surfaces
much later as an empty count rather than as an error.

Usage
-----
    # from an NCBI sequence report, matched to the names the index uses
    bowtie2-inspect -n my_index > names.txt
    chromfile_from_report.py --names names.txt \\
                             --report sequence_report.jsonl \\
                             --expect-autosomes 38 -o chroms.tsv

    # without a report, falling back to the names themselves
    chromfile_from_report.py --names names.txt --max-autosome 38 -o chroms.tsv

The names file is one sequence name per line; `bowtie2-inspect -n`,
`samtools idxstats | cut -f1` and `cut -f1 ref.fa.fai` all produce it.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

__version__ = "1.2.0"


def sam_name(name: str) -> str:
    """
    The name an aligner will write into the BAM.

    An index may store a full FASTA header -- accession followed by a
    description -- but SAM forbids whitespace in a reference name, so
    aligners emit only the part up to the first space. The chromosome file
    has to carry that short form, or it will match nothing.
    """
    parts = name.split()
    return parts[0] if parts else name


def build_lookup(names: list[str]) -> dict[str, str]:
    """
    Map every plausible spelling of a sequence to the name a BAM will carry.

    Index names vary more than one would hope. Some carry the description
    after the accession, some drop the version suffix, some differ only in
    case. Each variant is registered so a later lookup finds the short form.
    """
    lookup: dict[str, str] = {}

    def add(key: str, short: str) -> None:
        key = key.strip().lower()
        if key and key not in lookup:
            lookup[key] = short

    for name in names:
        short = sam_name(name)
        add(name, short)
        add(short, short)
        # Accessions with a version suffix: NC_051805.1 -> NC_051805
        if "." in short:
            add(short.rsplit(".", 1)[0], short)
    return lookup


def resolve(lookup: dict[str, str], *candidates: str | None) -> str | None:
    """First candidate that the index actually carries, as the index spells it."""
    for cand in candidates:
        if not cand:
            continue
        hit = lookup.get(str(cand).strip().lower())
        if hit:
            return hit
        if "." in str(cand):
            hit = lookup.get(str(cand).rsplit(".", 1)[0].strip().lower())
            if hit:
                return hit
    return None


def from_report(report: Path, lookup: dict[str, str]):
    """Roles from the sequence report, names as the index spells them."""
    autosomes: list[str] = []
    x = y = mt = None

    def pick(rec: dict) -> str | None:
        chrname = str(rec.get("chrName", ""))
        return resolve(
            lookup,
            rec.get("refseqAccession"), rec.get("genbankAccession"),
            rec.get("ucscStyleName"), f"chr{chrname}", chrname,
        )

    for line in report.read_text().splitlines():
        line = line.strip()
        if not line:
            continue
        rec = json.loads(line)
        if rec.get("role") != "assembled-molecule":
            continue
        used = pick(rec)
        if not used:
            continue
        loc = rec.get("assignedMoleculeLocationType", "")
        label = str(rec.get("chrName", "")).upper()
        if loc == "Mitochondrion" or label == "MT":
            mt = used
        elif loc == "Chromosome":
            if label == "X":
                x = used
            elif label == "Y":
                y = used
            else:
                autosomes.append(used)

    return autosomes, x, y, mt


def from_names(names: list[str], max_autosome: int):
    """
    Fall back to the names alone.

    Only numbered sequences up to `max_autosome` are taken as autosomes, so
    that unplaced scaffolds and the contigs of a co-mapped genome stay out.
    """
    autosomes: list[str] = []
    x = y = mt = None

    for n in names:
        short = sam_name(n)
        bare = short[3:] if short.lower().startswith("chr") else short
        upper = bare.upper()
        if upper == "X":
            x = short
        elif upper == "Y":
            y = short
        elif upper in ("M", "MT"):
            mt = short
        elif bare.isdigit() and 1 <= int(bare) <= max_autosome:
            autosomes.append(short)

    return autosomes, x, y, mt


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(
        description="Write a karyoratio chromosome definition file.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    p.add_argument("--names", required=True, type=Path,
                   help="one sequence name per line, as the index or BAM uses them")
    p.add_argument("--report", type=Path, default=None,
                   help="NCBI sequence_report.jsonl; without it the names "
                        "themselves are used")
    p.add_argument("--max-autosome", type=int, default=50,
                   help="highest chromosome number treated as an autosome in "
                        "the fallback (default: 50)")
    p.add_argument("--expect-autosomes", type=int, default=None,
                   help="warn if a different number of autosomes is found; "
                        "for the dog this is 38")
    p.add_argument("--comment", default=None,
                   help="header line written into the file")
    p.add_argument("-o", "--output", type=Path, default=None)
    p.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    args = p.parse_args(argv)

    names = [l.strip() for l in args.names.read_text().splitlines() if l.strip()]
    if not names:
        sys.exit(f"[error] {args.names} is empty")

    if args.report and args.report.exists():
        autosomes, x, y, mt = from_report(args.report, build_lookup(names))
        source = "NCBI sequence report"
        if x is None and not autosomes:
            print("[note] nothing in the sequence report matched the names, "
                  "so falling back\n       to the names themselves",
                  file=sys.stderr)
            autosomes, x, y, mt = from_names(names, args.max_autosome)
            source = "sequence names"
    else:
        if args.report:
            print(f"[note] {args.report} not found; falling back to names",
                  file=sys.stderr)
        autosomes, x, y, mt = from_names(names, args.max_autosome)
        source = "sequence names"

    if x is None or not autosomes:
        missing = "X chromosome" if x is None else "autosomes"
        print(f"[error] could not identify the {missing}.", file=sys.stderr)
        print(f"        {len(names)} names were given; the first few are:",
              file=sys.stderr)
        for n in names[:8]:
            print(f"          {n}", file=sys.stderr)
        print("        If these are not accessions, UCSC-style names or bare\n"
              "        chromosome numbers, write the chromosome file by hand:\n"
              "        two columns, <sequence name> and one of "
              "autosome/X/Y/MT.", file=sys.stderr)
        return 1

    out = args.output.open("w") if args.output else sys.stdout
    try:
        if args.comment:
            print(f"# {args.comment}", file=out)
        print(f"# roles from {source}; "
              f"{len(autosomes)} autosomes + X"
              + (", Y" if y else ", no Y in this assembly")
              + (", MT" if mt else ""), file=out)
        for acc in autosomes:
            print(f"{acc}\tautosome", file=out)
        print(f"{x}\tX", file=out)
        if y:
            print(f"{y}\tY", file=out)
        if mt:
            print(f"{mt}\tMT", file=out)
    finally:
        if args.output:
            out.close()

    print(f"  {len(autosomes)} autosomes, X = {x}"
          + (f", Y = {y}" if y else ", no Y")
          + (f", MT = {mt}" if mt else ""), file=sys.stderr)

    if args.expect_autosomes is not None and len(autosomes) != args.expect_autosomes:
        print(f"[warn] expected {args.expect_autosomes} autosomes but found "
              f"{len(autosomes)}; inspect the result before using it",
              file=sys.stderr)
        return 2

    return 0


if __name__ == "__main__":
    sys.exit(main())
