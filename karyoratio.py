#!/usr/bin/env python3
"""
karyoratio.py -- chromosome coverage ratios from an indexed BAM file.

By default computes Rx, the X-to-autosome read ratio of Mittnik et al.
(2016, PLoS ONE 11:e0163019), generalised to any reference genome with an
identifiable X chromosome, and estimates genetic sex from it.

The comparison itself is parameterised: any set of sequences can be pooled
as the target and compared against any other set as the reference. This
makes other karyotypic questions expressible -- a Y-to-autosome ratio, a
Z-to-autosome ratio in ZW systems, the relative copy number of a single
chromosome -- but note that the classification thresholds and the
correlation criterion were validated for Rx only. Outside the default
comparison the script computes the ratio and warns that no calibrated
decision rule exists for it.

Two differences from the original method:

  * read counts are taken from the BAM with an explicit mapping-quality
    filter (idxstats counts every mapped read regardless of MAPQ);
  * the correlation between chromosome length and read count, used as a
    validity criterion, is quantified with Spearman's rho rather than by
    linear regression.

Requires: pysam, numpy, scipy

Example
-------
    # default: Rx, sex assignment
    karyoratio.py --bam sample.bam --chrom-file chroms.tsv --min-mapq 30

    # ZW system: Z against the autosomes
    karyoratio.py --bam bird.bam --chrom-file bird.tsv \\
                  --target Z --reference autosome --labels ZW,n.a.,ZZ
"""

from __future__ import annotations

import argparse
import math
import sys
from collections import OrderedDict
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pysam
from scipy import stats

__version__ = "2.0.0"

# The comparison for which the thresholds below were established.
DEFAULT_TARGET = ("X",)
DEFAULT_REFERENCE = ("autosome",)

# Classification thresholds of Mittnik et al. (2016), applied to the bounds
# of the 95% confidence interval. Validated for Rx in diploid XY systems only.
DEFAULT_LOWER = 0.60
DEFAULT_UPPER = 0.80
DEFAULT_LABELS = ("XY", "n.a.", "XX")

CANONICAL_ROLES = {"autosome", "x", "y", "mt"}


# --------------------------------------------------------------------------- #
# Chromosome specification
# --------------------------------------------------------------------------- #

@dataclass
class ChromSpec:
    """Sequence name -> role. Roles are free-form; four are canonical."""
    roles: "OrderedDict[str, str]"

    def resolve(self, tokens: list[str]) -> list[str]:
        """
        Expand tokens into sequence names. A token is either a role label,
        which expands to every sequence carrying it, or a literal sequence
        name. Order is preserved and duplicates are dropped.
        """
        out: list[str] = []
        for token in tokens:
            by_role = [n for n, r in self.roles.items() if r == token]
            if by_role:
                out.extend(by_role)
            elif token in self.roles:
                out.append(token)
            else:
                known = sorted(set(self.roles.values()))
                raise KeyError(
                    f"{token!r} matches no role and no sequence name; "
                    f"roles present: {', '.join(known)}"
                )
        seen: set[str] = set()
        return [n for n in out if not (n in seen or seen.add(n))]

    def role_of(self, name: str) -> str:
        return self.roles.get(name, "?")


def read_chrom_file(path: Path) -> ChromSpec:
    """
    Read a two-column TSV of  <sequence name>\t<role>.

    Roles may be any label. Four are canonical and understood by the
    defaults: autosome, X, Y, MT. Anything after the second column, and any
    line starting with #, is ignored.

        NC_132780.1    autosome
        NC_132796.1    X
    """
    roles: "OrderedDict[str, str]" = OrderedDict()

    for lineno, raw in enumerate(path.read_text().splitlines(), start=1):
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split()
        if len(parts) < 2:
            raise ValueError(f"{path}:{lineno}: expected two columns, got {raw!r}")
        name, role = parts[0], parts[1]
        if role.lower() in CANONICAL_ROLES:
            role = {"autosome": "autosome", "x": "X", "y": "Y", "mt": "MT"}[role.lower()]
        if name in roles:
            raise ValueError(f"{path}:{lineno}: {name} listed more than once")
        roles[name] = role

    if not roles:
        raise ValueError(f"{path}: no sequences listed")
    return ChromSpec(roles=roles)


def write_chrom_template(
    bam_path: Path,
    x: str,
    y: str | None,
    mt: str | None,
    min_length: int,
    prefix: str | None,
    out,
) -> None:
    """
    Emit a chromosome definition file from the sequences in a BAM header.

    Everything at or above `min_length` (and matching `prefix`, if given) that
    is not the X, Y or mitochondrion is written as an autosome. Shorter
    sequences -- unplaced scaffolds, and the contigs of any co-mapped genome --
    are listed as comments so they can be inspected and promoted by hand.
    """
    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        seqs = list(zip(bam.references, bam.lengths))

    special = {name: role for name, role in
               ((x, "X"), (y, "Y"), (mt, "MT")) if name}
    missing = [n for n in special if n not in dict(seqs)]
    if missing:
        raise KeyError(f"{bam_path}: not in the BAM header: " + ", ".join(missing))

    autosomes, skipped = [], []
    for name, length in seqs:
        if name in special:
            continue
        if prefix and not name.startswith(prefix):
            skipped.append((name, length))
        elif length >= min_length:
            autosomes.append((name, length))
        else:
            skipped.append((name, length))

    print(f"# generated by karyoratio.py {__version__} from {bam_path.name}", file=out)
    print(f"# autosome criterion: length >= {min_length:,} bp"
          + (f", name starts with {prefix!r}" if prefix else ""), file=out)
    print("#", file=out)
    for name, length in autosomes:
        print(f"{name}\tautosome\t# {length:,} bp", file=out)
    for role, name in (("X", x), ("Y", y), ("MT", mt)):
        if name:
            print(f"{name}\t{role}", file=out)

    if skipped:
        print("#", file=out)
        print(f"# {len(skipped)} sequence(s) excluded; promote any that belong:", file=out)
        for name, length in skipped[:40]:
            print(f"# {name}\tautosome\t# {length:,} bp", file=out)
        if len(skipped) > 40:
            print(f"# ... and {len(skipped) - 40} more", file=out)


# --------------------------------------------------------------------------- #
# Counting
# --------------------------------------------------------------------------- #

def count_reads(
    bam_path: Path,
    names: list[str],
    min_mapq: int,
    require_proper_pair: bool,
    count_duplicates: bool,
) -> tuple[dict[str, int], dict[str, int]]:
    """
    Count filtered reads per reference sequence.

    Returns (counts, lengths); both keyed by reference name. Sequence lengths
    are taken from the BAM header.
    """
    counts: dict[str, int] = {}
    lengths: dict[str, int] = {}

    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        if not bam.has_index():
            raise RuntimeError(
                f"{bam_path}: no index found. Run `samtools index {bam_path}` first."
            )

        header_lengths = dict(zip(bam.references, bam.lengths))
        missing = [n for n in names if n not in header_lengths]
        if missing:
            raise KeyError(
                f"{bam_path}: these sequences are not in the BAM header: "
                + ", ".join(missing[:5])
                + (" ..." if len(missing) > 5 else "")
            )

        def keep(read: pysam.AlignedSegment) -> bool:
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                return False
            if read.is_qcfail:
                return False
            if read.is_duplicate and not count_duplicates:
                return False
            if read.mapping_quality < min_mapq:
                return False
            if require_proper_pair and read.is_paired and not read.is_proper_pair:
                return False
            return True

        for name in names:
            lengths[name] = header_lengths[name]
            counts[name] = bam.count(contig=name, read_callback=keep)

    return counts, lengths


# --------------------------------------------------------------------------- #
# Ratio
# --------------------------------------------------------------------------- #

@dataclass
class RatioResult:
    sample: str
    n_target: int
    n_reference: int
    n_mapped: int
    ratio: float
    CI_l: float
    CI_u: float
    rho: float
    rho_p: float
    assignment: str

    def row(self) -> str:
        return (
            f"{self.sample}\t{self.n_target}\t{self.n_reference}\t{self.n_mapped}\t"
            f"{self.ratio:.4f}\t{self.CI_l:.4f}\t{self.CI_u:.4f}\t"
            f"{self.rho:.4f}\t{self.rho_p:.3g}\t{self.assignment}"
        )

    @staticmethod
    def header() -> str:
        return ("sample\tn_target\tn_reference\tn_mapped\t"
                "ratio\tCI_l\tCI_u\trho\trho_p\tassignment")


def classify(
    ci_l: float,
    ci_u: float,
    n_mapped: int,
    rho: float,
    lower: float,
    upper: float,
    labels: tuple[str, str, str],
    min_reads: int = 0,
    min_rho: float = 0.0,
) -> str:
    """
    Assign a label from the confidence interval, subject to validity criteria.

    The interval must clear a threshold entirely: `labels[0]` if its upper
    bound is below `lower`, `labels[2]` if its lower bound is above `upper`,
    otherwise `labels[1]`.

    Two optional guards exist because thresholds separate cleanly only at
    adequate coverage. Below it the distributions overlap, so an interval can
    clear a threshold on a sample carrying too few reads to support any call.
    """
    if min_reads and n_mapped < min_reads:
        return labels[1]
    if min_rho and rho < min_rho:
        return labels[1]
    if ci_u < lower:
        return labels[0]
    if ci_l > upper:
        return labels[2]
    return labels[1]


def compute_ratio(
    counts: dict[str, int],
    lengths: dict[str, int],
    target: list[str],
    reference: list[str],
    sample: str,
    lower: float = DEFAULT_LOWER,
    upper: float = DEFAULT_UPPER,
    labels: tuple[str, str, str] = DEFAULT_LABELS,
    min_reads: int = 0,
    min_rho: float = 0.0,
) -> RatioResult:
    """
    Pooled target-to-reference coverage ratio.

    The target sequences are pooled: their reads and their lengths are summed,
    so the ratio refers to one unit rather than to a mean over several. The
    pooled target density is then compared with the density of each reference
    sequence separately,

        ratio_i = (m_T / s_T) / (m_i / s_i)

    and the reported value is the mean of ratio_i, with a confidence interval
    from the standard error of that mean.

    Note that the normalising totals used in the original formulation cancel
    in this quotient, so sequences present in neither set do not affect the
    result. What a sequence does affect is the mean, if it is in `reference`.
    """
    tgt = [n for n in target if counts.get(n, 0) > 0]
    ref = [n for n in reference if counts.get(n, 0) > 0]

    overlap = set(target) & set(reference)
    if overlap:
        raise ValueError(
            f"{sample}: sequence(s) in both target and reference: "
            + ", ".join(sorted(overlap))
        )
    if not tgt:
        raise ValueError(f"{sample}: no reads on any target sequence")
    if not ref:
        raise ValueError(f"{sample}: no reads on any reference sequence")

    m_t = float(sum(counts[n] for n in tgt))
    s_t = float(sum(lengths[n] for n in tgt))
    density_t = m_t / s_t

    m_ref = np.array([counts[n] for n in ref], dtype=float)
    s_ref = np.array([lengths[n] for n in ref], dtype=float)
    ratios = density_t / (m_ref / s_ref)

    value = float(ratios.mean())
    if ratios.size > 1:
        se = float(ratios.std(ddof=1) / math.sqrt(ratios.size))
    else:
        se = float("nan")
    ci_l, ci_u = value - 1.96 * se, value + 1.96 * se

    # Validity criterion over every sequence that entered the comparison.
    used = tgt + ref
    rho, rho_p = stats.spearmanr(
        [lengths[n] for n in used], [counts[n] for n in used]
    )

    n_mapped = int(m_t + m_ref.sum())

    return RatioResult(
        sample=sample,
        n_target=len(tgt),
        n_reference=len(ref),
        n_mapped=n_mapped,
        ratio=value,
        CI_l=ci_l,
        CI_u=ci_u,
        rho=float(rho),
        rho_p=float(rho_p),
        assignment=classify(ci_l, ci_u, n_mapped, float(rho),
                            lower=lower, upper=upper, labels=labels,
                            min_reads=min_reads, min_rho=min_rho),
    )


# --------------------------------------------------------------------------- #
# CLI
# --------------------------------------------------------------------------- #

def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Chromosome coverage ratios from a BAM file; by default "
                    "the X-to-autosome ratio (Rx) and a sex assignment.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    p.add_argument("--bam", nargs="+", required=True, type=Path,
                   help="one or more indexed BAM files")

    g = p.add_argument_group("chromosome specification")
    g.add_argument("--chrom-file", type=Path,
                   help="TSV of <sequence name> <role>; roles are free-form, "
                        "with autosome, X, Y and MT understood by the defaults")
    g.add_argument("--autosomes",
                   help="comma-separated autosome names (alternative to --chrom-file)")
    g.add_argument("--x", help="name of the X chromosome")
    g.add_argument("--y", default=None, help="name of the Y chromosome, if present")
    g.add_argument("--mt", default=None, help="name of the mitochondrial sequence")

    c = p.add_argument_group("comparison")
    c.add_argument("--target", default=",".join(DEFAULT_TARGET),
                   help="roles or sequence names forming the numerator; pooled "
                        f"if more than one (default: {','.join(DEFAULT_TARGET)})")
    c.add_argument("--reference", default=",".join(DEFAULT_REFERENCE),
                   help="roles or sequence names forming the denominator; each "
                        "compared separately (default: "
                        f"{','.join(DEFAULT_REFERENCE)})")
    c.add_argument("--lower", type=float, default=DEFAULT_LOWER,
                   help=f"CI upper bound below this gives the first label "
                        f"(default: {DEFAULT_LOWER})")
    c.add_argument("--upper", type=float, default=DEFAULT_UPPER,
                   help=f"CI lower bound above this gives the third label "
                        f"(default: {DEFAULT_UPPER})")
    c.add_argument("--labels", default=",".join(DEFAULT_LABELS),
                   help="three comma-separated labels: below, ambiguous, above "
                        f"(default: {','.join(DEFAULT_LABELS)})")

    f = p.add_argument_group("read filtering")
    f.add_argument("--min-mapq", type=int, default=30,
                   help="minimum mapping quality (default: 30)")
    f.add_argument("--proper-pairs-only", action="store_true",
                   help="count only properly paired reads; off by default, since "
                        "orphan reads carry real signal in low-coverage data")
    f.add_argument("--count-duplicates", action="store_true",
                   help="include reads flagged as duplicates (off by default)")

    v = p.add_argument_group("validity criteria")
    v.add_argument("--min-reads", type=int, default=0,
                   help="report the ambiguous label below this many mapped "
                        "reads, whatever the confidence interval says; "
                        "species-specific, so off by default")
    v.add_argument("--min-rho", type=float, default=0.0,
                   help="report the ambiguous label below this Spearman "
                        "correlation between sequence length and read count; "
                        "off by default, and calibrated for Rx only")

    t = p.add_argument_group("chromosome file generation")
    t.add_argument("--make-chrom-file", action="store_true",
                   help="instead of computing a ratio, write a chromosome "
                        "definition file built from the header of the first "
                        "BAM; requires --x")
    t.add_argument("--min-autosome-length", type=int, default=10_000_000,
                   help="sequences at or above this length become autosomes "
                        "(default: 10,000,000)")
    t.add_argument("--prefix", default=None,
                   help="only consider sequences whose name starts with this, "
                        "useful when the reference is a composite genome")

    p.add_argument("--quiet", action="store_true",
                   help="suppress the warning issued for non-default comparisons")
    p.add_argument("-o", "--output", type=Path, default=None,
                   help="write TSV here instead of stdout")
    p.add_argument("--per-chrom", type=Path, default=None,
                   help="also write per-chromosome counts to this TSV")
    p.add_argument("--version", action="version", version=f"%(prog)s {__version__}")

    args = p.parse_args(argv)

    if args.make_chrom_file:
        if not args.x:
            p.error("--make-chrom-file requires --x")
    elif args.chrom_file is None:
        if not args.autosomes or not args.x:
            p.error("give either --chrom-file, or both --autosomes and --x")

    args.labels = tuple(s.strip() for s in args.labels.split(","))
    if len(args.labels) != 3:
        p.error("--labels needs exactly three comma-separated values")
    return args


def warn_if_uncalibrated(args: argparse.Namespace) -> None:
    """
    The thresholds and the correlation criterion were established for Rx.
    Any departure from that comparison leaves them without support, so say so
    rather than returning a confident-looking label built on nothing.
    """
    target = tuple(s.strip() for s in args.target.split(","))
    reference = tuple(s.strip() for s in args.reference.split(","))

    changed = []
    if target != DEFAULT_TARGET or reference != DEFAULT_REFERENCE:
        changed.append(f"comparison ({','.join(target)} vs {','.join(reference)})")
    if (args.lower, args.upper) != (DEFAULT_LOWER, DEFAULT_UPPER):
        changed.append(f"thresholds ({args.lower}/{args.upper})")

    if not changed:
        return

    print(
        "[warn] this is not the Rx comparison the thresholds were validated for: "
        + "; ".join(changed) + ".\n"
        "       The ratio and its confidence interval are still computed as "
        "defined, but the\n"
        "       0.60/0.80 decision rule of Mittnik et al. does not carry over, "
        "and neither does\n"
        "       any --min-rho or --min-reads cutoff derived for Rx. Calibrate "
        "by simulation\n"
        "       before treating the assignment column as meaningful.",
        file=sys.stderr,
    )

    if args.min_rho:
        print(
            "[warn] --min-rho was calibrated for Rx, where the correlation "
            "behaves differently\n"
            "       between the sexes; its behaviour under this comparison is "
            "unknown.",
            file=sys.stderr,
        )


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)

    if args.make_chrom_file:
        out = args.output.open("w") if args.output else sys.stdout
        try:
            write_chrom_template(
                args.bam[0], args.x, args.y, args.mt,
                min_length=args.min_autosome_length,
                prefix=args.prefix,
                out=out,
            )
        except KeyError as exc:
            print(f"[error] {exc}", file=sys.stderr)
            return 1
        finally:
            if args.output:
                out.close()
        return 0

    if args.chrom_file:
        spec = read_chrom_file(args.chrom_file)
    else:
        roles: "OrderedDict[str, str]" = OrderedDict()
        for name in (c.strip() for c in args.autosomes.split(",")):
            if name:
                roles[name] = "autosome"
        roles[args.x] = "X"
        if args.y:
            roles[args.y] = "Y"
        if args.mt:
            roles[args.mt] = "MT"
        spec = ChromSpec(roles=roles)

    try:
        target = spec.resolve([s.strip() for s in args.target.split(",")])
        reference = spec.resolve([s.strip() for s in args.reference.split(",")])
    except KeyError as exc:
        print(f"[error] {exc}", file=sys.stderr)
        return 1

    if not args.quiet:
        warn_if_uncalibrated(args)

    names = target + [n for n in reference if n not in target]

    out = args.output.open("w") if args.output else sys.stdout
    per_chrom = args.per_chrom.open("w") if args.per_chrom else None
    if per_chrom:
        per_chrom.write("sample\tchrom\trole\tset\tlength\tmapped\n")

    try:
        print(RatioResult.header(), file=out)
        failed = 0

        for bam_path in args.bam:
            sample = bam_path.name
            for suffix in (".sorted.bam", ".bam"):
                if sample.endswith(suffix):
                    sample = sample[: -len(suffix)]
                    break
            try:
                counts, lengths = count_reads(
                    bam_path,
                    names,
                    min_mapq=args.min_mapq,
                    require_proper_pair=args.proper_pairs_only,
                    count_duplicates=args.count_duplicates,
                )
                result = compute_ratio(
                    counts, lengths, target, reference, sample,
                    lower=args.lower, upper=args.upper, labels=args.labels,
                    min_reads=args.min_reads, min_rho=args.min_rho,
                )
            except (RuntimeError, KeyError, ValueError) as exc:
                print(f"[warn] {exc}", file=sys.stderr)
                failed += 1
                continue

            print(result.row(), file=out)

            if per_chrom:
                for name in names:
                    which = "target" if name in target else "reference"
                    per_chrom.write(
                        f"{sample}\t{name}\t{spec.role_of(name)}\t{which}\t"
                        f"{lengths[name]}\t{counts[name]}\n"
                    )
    finally:
        if args.output:
            out.close()
        if per_chrom:
            per_chrom.close()

    return 1 if failed and failed == len(args.bam) else 0


if __name__ == "__main__":
    sys.exit(main())
