#!/usr/bin/env python3
"""
calibrate.py -- establish karyoratio thresholds for a new species by simulation.

The 0.60/0.80 decision rule of Mittnik et al. (2016) was derived from human
data. Before applying it to another species one needs to know two things:
whether the thresholds separate the sexes at all, and how many reads a sample
must carry for the separation to hold. This script answers both.

What it does
------------
1. Builds artificial female and male genomes from the reference, by repeating
   sequences to the copy number each sex carries.
2. Simulates reads from both with ART, optionally using an error profile
   derived from real data.
3. Maps each simulated set to the unmodified reference, once.
4. Subsamples the mapped reads in memory at a range of read counts, with
   replicates, and computes the ratio for each.
5. Reports, per read count and sex, how many replicates were assigned
   correctly, how many were left ambiguous, and how many were assigned to the
   wrong sex.

Subsampling after mapping rather than before is what makes this cheap: two
alignments are run instead of one per replicate. It also fixes the read count
to mapped reads, which is the quantity karyoratio reports, removing the
ambiguity between generated and mapped counts.

Requires: pysam, numpy, scipy, and on PATH: art_illumina, bowtie2 (or bwa),
samtools. art_profiler_illumina is needed only with --profile-from.

Example
-------
    calibrate.py --reference mMarFoi2.1.fa \\
                 --chrom-file mmarfoi2.1_chroms.tsv \\
                 --outdir calib_marten \\
                 --profile-from real_R1.fq.gz \\
                 --replicates 30 --threads 8
"""

from __future__ import annotations

import argparse
import importlib.util
import shutil
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pysam

__version__ = "1.0.0"

DEFAULT_COUNTS = (1000, 1500, 2000, 2500, 3000, 5000, 7500, 10000,
                  25000, 50000, 100000, 250000)


# --------------------------------------------------------------------------- #
# karyoratio is the single source of truth for the ratio and the decision rule
# --------------------------------------------------------------------------- #

def load_karyoratio(path: Path | None = None):
    candidate = path or (Path(__file__).resolve().parent / "karyoratio.py")
    if not candidate.exists():
        sys.exit(f"[error] karyoratio.py not found next to this script "
                 f"(looked at {candidate}); pass --karyoratio")
    spec = importlib.util.spec_from_file_location("karyoratio", candidate)
    module = importlib.util.module_from_spec(spec)
    sys.modules["karyoratio"] = module
    spec.loader.exec_module(module)
    return module


# --------------------------------------------------------------------------- #
# External tools
# --------------------------------------------------------------------------- #

def require(*tools: str) -> None:
    missing = [t for t in tools if shutil.which(t) is None]
    if missing:
        sys.exit("[error] not found on PATH: " + ", ".join(missing))


def run(cmd: list[str], log: Path | None = None) -> None:
    printable = " ".join(str(c) for c in cmd)
    print(f"[run] {printable[:160]}{'...' if len(printable) > 160 else ''}",
          file=sys.stderr)
    handle = log.open("w") if log else None
    try:
        subprocess.run([str(c) for c in cmd], check=True,
                       stdout=handle or subprocess.DEVNULL,
                       stderr=subprocess.STDOUT if handle else subprocess.DEVNULL)
    except subprocess.CalledProcessError as exc:
        sys.exit(f"[error] command failed ({exc.returncode}): {printable}"
                 + (f"\n        see {log}" if log else ""))
    finally:
        if handle:
            handle.close()


# --------------------------------------------------------------------------- #
# Stage 1: artificial genomes
# --------------------------------------------------------------------------- #

def build_genomes(
    reference: Path,
    autosomes: list[str],
    x: str,
    y: str | None,
    outdir: Path,
    skip_existing: bool,
) -> tuple[Path, Path]:
    """
    Write a female and a male genome at the copy numbers each sex carries.

    Female: every autosome twice, X twice.
    Male:   every autosome twice, X once, plus Y once if the reference has one.

    Repeating a sequence doubles its share of the simulated reads, which is the
    only property of ploidy the ratio depends on. The duplicated X of the
    female is a copy of one haplotype, so the simulation models copy number,
    not genetic variability -- it cannot speak to heterozygosity.
    """
    female = outdir / "genome_female.fa"
    male = outdir / "genome_male.fa"
    if skip_existing and female.exists() and male.exists():
        print("[skip] artificial genomes already present", file=sys.stderr)
        return female, male

    fasta = pysam.FastaFile(str(reference))
    present = set(fasta.references)
    for name in autosomes + [x] + ([y] if y else []):
        if name not in present:
            sys.exit(f"[error] {name} is not in {reference}")

    def emit(handle, name: str, copies: int) -> None:
        seq = fasta.fetch(name)
        for copy in range(copies):
            handle.write(f">{name}_copy{copy + 1}\n")
            for i in range(0, len(seq), 60):
                handle.write(seq[i:i + 60] + "\n")

    print(f"[stage 1] writing artificial genomes to {outdir}", file=sys.stderr)
    with female.open("w") as fh:
        for name in autosomes:
            emit(fh, name, 2)
        emit(fh, x, 2)
    with male.open("w") as fh:
        for name in autosomes:
            emit(fh, name, 2)
        emit(fh, x, 1)
        if y:
            emit(fh, y, 1)

    if not y:
        print("[note] the reference has no Y chromosome, so the male genome "
              "carries a single X\n"
              "       alongside the doubled autosomes. This is the correct "
              "construction for the\n"
              "       X-to-autosome ratio, which does not use the Y.",
              file=sys.stderr)
    return female, male


# --------------------------------------------------------------------------- #
# Stage 2: read simulation
# --------------------------------------------------------------------------- #

def build_profile(reads: Path, outdir: Path, skip_existing: bool) -> Path:
    """Derive an ART error profile from real reads."""
    prefix = outdir / "profile"
    first = Path(str(prefix) + "R1.txt")
    if skip_existing and first.exists():
        print("[skip] error profile already present", file=sys.stderr)
        return prefix

    require("art_profiler_illumina")
    folder = outdir / "profile_input"
    folder.mkdir(exist_ok=True)
    link = folder / reads.name
    if not link.exists():
        link.symlink_to(reads.resolve())

    suffix = "".join(reads.suffixes)
    print(f"[stage 2a] deriving error profile from {reads.name}", file=sys.stderr)
    run(["art_profiler_illumina", str(prefix), str(folder), suffix],
        log=outdir / "log_profiler.txt")
    return prefix


def simulate(
    genome: Path,
    label: str,
    outdir: Path,
    read_length: int,
    total_reads: int,
    fragment_mean: int,
    fragment_sd: int,
    profile: Path | None,
    seed: int,
    skip_existing: bool,
) -> tuple[Path, Path]:
    """Generate paired-end reads from one artificial genome."""
    prefix = outdir / f"sim_{label}_"
    r1, r2 = Path(f"{prefix}1.fq"), Path(f"{prefix}2.fq")
    if skip_existing and r1.exists() and r2.exists():
        print(f"[skip] simulated reads for {label} already present", file=sys.stderr)
        return r1, r2

    genome_bp = sum(len(pysam.FastaFile(str(genome)).fetch(r))
                    for r in pysam.FastaFile(str(genome)).references)
    # ART takes fold coverage; convert from the requested read total.
    fold = (total_reads * read_length) / genome_bp

    cmd = ["art_illumina", "-ss", "HS25", "-i", genome, "-p",
           "-l", read_length, "-f", f"{fold:.6f}",
           "-m", fragment_mean, "-s", fragment_sd,
           "-rs", seed, "-na", "-o", prefix]
    if profile:
        cmd = ["art_illumina", "-1", f"{profile}R1.txt", "-2", f"{profile}R2.txt",
               "-i", genome, "-p", "-l", read_length, "-f", f"{fold:.6f}",
               "-m", fragment_mean, "-s", fragment_sd,
               "-rs", seed, "-na", "-o", prefix]

    print(f"[stage 2b] simulating ~{total_reads:,} read pairs for {label} "
          f"({fold:.3f}x)", file=sys.stderr)
    run(cmd, log=outdir / f"log_art_{label}.txt")
    return r1, r2


# --------------------------------------------------------------------------- #
# Stage 3: mapping
# --------------------------------------------------------------------------- #

def build_index(reference: Path, aligner: str, outdir: Path,
                threads: int, skip_existing: bool) -> Path:
    index = outdir / "ref_index"
    marker = Path(str(index) + (".1.bt2" if aligner == "bowtie2" else ".bwt"))
    if skip_existing and marker.exists():
        print("[skip] aligner index already present", file=sys.stderr)
        return index

    print(f"[stage 3a] building {aligner} index", file=sys.stderr)
    if aligner == "bowtie2":
        run(["bowtie2-build", "--threads", threads, reference, index],
            log=outdir / "log_index.txt")
    else:
        run(["bwa", "index", "-p", index, reference], log=outdir / "log_index.txt")
    return index


def map_reads(r1: Path, r2: Path, index: Path, label: str, aligner: str,
              outdir: Path, threads: int, skip_existing: bool) -> Path:
    bam = outdir / f"sim_{label}.bam"
    if skip_existing and bam.exists() and Path(str(bam) + ".bai").exists():
        print(f"[skip] alignment for {label} already present", file=sys.stderr)
        return bam

    print(f"[stage 3b] mapping {label}", file=sys.stderr)
    log = outdir / f"log_map_{label}.txt"
    if aligner == "bowtie2":
        align = ["bowtie2", "--very-sensitive-local", "-p", str(threads),
                 "-x", str(index), "-1", str(r1), "-2", str(r2)]
    else:
        align = ["bwa", "mem", "-t", str(threads), str(index), str(r1), str(r2)]

    with log.open("w") as errlog:
        proc = subprocess.Popen(align, stdout=subprocess.PIPE, stderr=errlog)
        sort = subprocess.Popen(
            ["samtools", "sort", "-@", str(threads), "-o", str(bam)],
            stdin=proc.stdout, stderr=errlog)
        proc.stdout.close()
        sort.communicate()
        proc.wait()
    if sort.returncode or proc.returncode:
        sys.exit(f"[error] mapping failed for {label}; see {log}")

    run(["samtools", "index", bam])
    return bam


# --------------------------------------------------------------------------- #
# Stage 4: subsampling and ratio
# --------------------------------------------------------------------------- #

@dataclass
class Pool:
    """Every filtered read of one simulated sample, as a sequence index."""
    indices: np.ndarray
    names: list[str]
    lengths: dict[str, int]
    total_reads: int


def collect_pool(bam_path: Path, names: list[str], min_mapq: int) -> Pool:
    """
    Read the BAM once and record, for every passing read, which sequence it
    mapped to. Subsampling then happens in memory.
    """
    name_to_idx = {n: i for i, n in enumerate(names)}
    lengths: dict[str, int] = {}
    collected: list[np.ndarray] = []
    total = 0

    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        header_lengths = dict(zip(bam.references, bam.lengths))
        for name in names:
            if name not in header_lengths:
                sys.exit(f"[error] {name} not in {bam_path} header")
            lengths[name] = header_lengths[name]

        for name in names:
            idx = name_to_idx[name]
            n = bam.count(
                contig=name,
                read_callback=lambda r: not (
                    r.is_unmapped or r.is_secondary or r.is_supplementary
                    or r.is_qcfail or r.is_duplicate
                    or r.mapping_quality < min_mapq
                ),
            )
            collected.append(np.full(n, idx, dtype=np.int32))
            total += n

    return Pool(indices=np.concatenate(collected) if collected
                else np.empty(0, dtype=np.int32),
                names=names, lengths=lengths, total_reads=total)


def subsample_counts(pool: Pool, n: int, rng: np.random.Generator) -> dict[str, int]:
    """Draw n reads without replacement and tally them per sequence."""
    drawn = rng.choice(pool.indices, size=n, replace=False)
    tally = np.bincount(drawn, minlength=len(pool.names))
    return {name: int(tally[i]) for i, name in enumerate(pool.names)}


# --------------------------------------------------------------------------- #
# Driver
# --------------------------------------------------------------------------- #

def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Calibrate karyoratio thresholds for a species by simulation.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    p.add_argument("--reference", required=True, type=Path,
                   help="indexed reference FASTA (samtools faidx)")
    p.add_argument("--chrom-file", required=True, type=Path,
                   help="karyoratio chromosome definition file")
    p.add_argument("--outdir", required=True, type=Path)

    s = p.add_argument_group("simulation")
    s.add_argument("--profile-from", type=Path, default=None,
                   help="real FASTQ from which to derive an ART error profile; "
                        "without it ART's built-in HiSeq 2500 profile is used")
    s.add_argument("--read-length", type=int, default=151)
    s.add_argument("--total-reads", type=int, default=2_000_000,
                   help="read pairs to simulate per sex (default: 2,000,000)")
    s.add_argument("--fragment-mean", type=int, default=300)
    s.add_argument("--fragment-sd", type=int, default=50)
    s.add_argument("--seed", type=int, default=1)

    a = p.add_argument_group("alignment")
    a.add_argument("--aligner", choices=("bowtie2", "bwa"), default="bowtie2")
    a.add_argument("--threads", type=int, default=4)
    a.add_argument("--min-mapq", type=int, default=30)

    c = p.add_argument_group("subsampling")
    c.add_argument("--read-counts", default=",".join(str(c) for c in DEFAULT_COUNTS),
                   help="mapped-read counts to test, comma-separated")
    c.add_argument("--replicates", type=int, default=30)

    p.add_argument("--karyoratio", type=Path, default=None,
                   help="path to karyoratio.py (default: next to this script)")
    p.add_argument("--skip-existing", action="store_true",
                   help="reuse intermediate files from an earlier run")
    p.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    return p.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    kr = load_karyoratio(args.karyoratio)

    require("art_illumina", "samtools",
            "bowtie2" if args.aligner == "bowtie2" else "bwa")
    if args.aligner == "bowtie2":
        require("bowtie2-build")

    args.outdir.mkdir(parents=True, exist_ok=True)
    spec = kr.read_chrom_file(args.chrom_file)
    autosomes = spec.resolve(["autosome"])
    x = spec.resolve(["X"])[0]
    y_list = [n for n, r in spec.roles.items() if r == "Y"]
    y = y_list[0] if y_list else None

    counts = sorted(int(c) for c in args.read_counts.split(",") if c.strip())

    female_fa, male_fa = build_genomes(
        args.reference, autosomes, x, y, args.outdir, args.skip_existing)

    profile = (build_profile(args.profile_from, args.outdir, args.skip_existing)
               if args.profile_from else None)

    index = build_index(args.reference, args.aligner, args.outdir,
                        args.threads, args.skip_existing)

    pools: dict[str, Pool] = {}
    for label, genome in (("female", female_fa), ("male", male_fa)):
        r1, r2 = simulate(genome, label, args.outdir, args.read_length,
                          args.total_reads, args.fragment_mean, args.fragment_sd,
                          profile, args.seed, args.skip_existing)
        bam = map_reads(r1, r2, index, label, args.aligner, args.outdir,
                        args.threads, args.skip_existing)
        print(f"[stage 4a] collecting reads for {label}", file=sys.stderr)
        pools[label] = collect_pool(bam, autosomes + [x], args.min_mapq)
        print(f"           {pools[label].total_reads:,} reads pass the filters "
              f"({pools[label].total_reads / (2 * args.total_reads):.1%} of "
              f"simulated)", file=sys.stderr)

    usable = [n for n in counts
              if all(n <= p.total_reads for p in pools.values())]
    dropped = [n for n in counts if n not in usable]
    if dropped:
        print(f"[note] skipping read counts above the simulated pool: "
              f"{', '.join(str(d) for d in dropped)}", file=sys.stderr)
    if not usable:
        sys.exit("[error] no requested read count fits in the simulated pool; "
                 "raise --total-reads")

    rng = np.random.default_rng(args.seed)
    per_rep = args.outdir / "replicates.tsv"
    print(f"[stage 4b] {len(usable)} read counts x {args.replicates} replicates "
          f"x 2 sexes", file=sys.stderr)

    rows = []
    with per_rep.open("w") as fh:
        fh.write("sex\tread_count\treplicate\tratio\tCI_l\tCI_u\trho\trho_p\t"
                 "assignment\n")
        for sex, pool in pools.items():
            expected = "XX" if sex == "female" else "XY"
            for n in usable:
                for rep in range(1, args.replicates + 1):
                    counts_d = subsample_counts(pool, n, rng)
                    try:
                        res = kr.compute_ratio(
                            counts_d, pool.lengths, [x], autosomes,
                            f"{sex}_{n}_{rep}")
                    except ValueError:
                        continue
                    fh.write(f"{sex}\t{n}\t{rep}\t{res.ratio:.6f}\t"
                             f"{res.CI_l:.6f}\t{res.CI_u:.6f}\t{res.rho:.6f}\t"
                             f"{res.rho_p:.3g}\t{res.assignment}\n")
                    rows.append((sex, n, expected, res.assignment))

    summary = args.outdir / "summary.tsv"
    print(f"[stage 5] writing summary", file=sys.stderr)
    with summary.open("w") as fh:
        fh.write("sex\tread_count\tn\tcorrect\tambiguous\twrong\tpct_ambiguous\n")
        for sex in ("female", "male"):
            for n in usable:
                sel = [r for r in rows if r[0] == sex and r[1] == n]
                if not sel:
                    continue
                correct = sum(1 for r in sel if r[3] == r[2])
                amb = sum(1 for r in sel if r[3] == "n.a.")
                wrong = len(sel) - correct - amb
                fh.write(f"{sex}\t{n}\t{len(sel)}\t{correct}\t{amb}\t{wrong}\t"
                         f"{100 * amb / len(sel):.2f}\n")

    wrong_total = sum(1 for r in rows if r[3] != r[2] and r[3] != "n.a.")
    clean = [n for n in usable
             if all(r[3] == r[2] for r in rows if r[1] == n)]

    print(f"\n[result] per-replicate values: {per_rep}", file=sys.stderr)
    print(f"[result] summary:              {summary}", file=sys.stderr)
    print(f"[result] incorrect assignments: {wrong_total} of {len(rows)}",
          file=sys.stderr)
    if clean:
        print(f"[result] every replicate assigned correctly at and above "
              f"{min(clean):,} mapped reads", file=sys.stderr)
        print(f"         -> a reasonable --min-reads for this species",
              file=sys.stderr)
    else:
        print("[result] no read count classified every replicate; raise "
              "--total-reads or widen --read-counts", file=sys.stderr)
    if wrong_total:
        print("[warn] incorrect assignments occurred, so the 0.60/0.80 "
              "thresholds do not transfer\n"
              "       cleanly to this species. Inspect replicates.tsv before "
              "using them.", file=sys.stderr)

    return 0


if __name__ == "__main__":
    sys.exit(main())
