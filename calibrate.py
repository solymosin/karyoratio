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

__version__ = "2.0.0"

DEFAULT_COUNTS = (1000, 1500, 2000, 2500, 3000, 5000, 7500, 10000,
                  25000, 50000, 100000, 250000)

# Longest read each of ART's built-in error profiles will generate. Asking for
# more fails, with a message that does not make the cause obvious.
ART_SYSTEMS = {
    "GA1": 44, "GA2": 75, "HS10": 100, "HS20": 100, "HS25": 150,
    "HSXn": 150, "HSXt": 150, "MinS": 50, "MSv1": 250, "MSv3": 250,
    "NS50": 75,
}


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
        msg = [f"[error] command failed ({exc.returncode}): {printable}"]
        if log and log.exists():
            tail = log.read_text().strip().splitlines()[-12:]
            if tail:
                msg.append(f"        last lines of {log}:")
                msg.extend(f"          {line}" for line in tail)
        sys.exit("\n".join(msg))
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

    # Paired simulation needs a profile for each mate, and the profiler builds
    # them from whatever it finds in the folder. Given only R1 it writes only
    # profileR1.txt, and ART then fails on the missing R2. So the mate is
    # looked for alongside, under the usual naming conventions.
    mates = [reads]
    for a, b in (("_1.", "_2."), ("_R1", "_R2"), (".1.", ".2."),
                 ("_1_", "_2_"), ("R1_", "R2_")):
        if a in reads.name:
            candidate = reads.with_name(reads.name.replace(a, b, 1))
            if candidate.exists() and candidate != reads:
                mates.append(candidate)
                break

    for path in mates:
        link = folder / path.name
        if not link.exists():
            link.symlink_to(path.resolve())

    if len(mates) == 1:
        print(f"[note] no second mate found next to {reads.name}; the R1 "
              f"profile will be\n"
              f"       used for both mates. Name the pair _1/_2 or _R1/_R2 "
              f"to have each\n"
              f"       mate profiled separately.", file=sys.stderr)
    else:
        print(f"[note] profiling both mates: "
              f"{', '.join(p.name for p in mates)}", file=sys.stderr)

    # art_profiler_illumina builds the glob itself as "*." + suffix, so the
    # suffix must carry no leading dot. Everything after the first dot in the
    # name is taken, which handles sample.1.fastq.gz as well as reads.fq.
    name = reads.name
    dot = name.find(".")
    if dot < 0:
        sys.exit(f"[error] {reads} has no extension; art_profiler_illumina "
                 f"needs one to find the file")
    suffix = name[dot + 1:]

    print(f"[stage 2a] deriving error profile from {reads.name} "
          f"(extension {suffix!r})", file=sys.stderr)
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
    art_system: str = "HS25",
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

    cmd = ["art_illumina", "-ss", art_system, "-i", genome, "-p",
           "-l", read_length, "-f", f"{fold:.6f}",
           "-m", fragment_mean, "-s", fragment_sd,
           "-rs", seed, "-na", "-o", prefix]
    if profile:
        p1 = Path(f"{profile}R1.txt")
        p2 = Path(f"{profile}R2.txt")
        if not p1.exists():
            sys.exit(f"[error] no error profile at {p1}; the profiler step "
                     f"did not produce one")
        # With only one mate profiled, R1's profile stands in for R2.
        if not p2.exists():
            p2 = p1
        cmd = ["art_illumina", "-1", p1, "-2", p2,
               "-i", genome, "-p", "-l", read_length, "-f", f"{fold:.6f}",
               "-m", fragment_mean, "-s", fragment_sd,
               "-rs", seed, "-na", "-o", prefix]

    print(f"[stage 2b] simulating ~{total_reads // 2:,} read pairs "
          f"({total_reads:,} reads) for {label} "
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
    """
    The simulated sample as pair categories.

    Reads come in pairs, and the two mates of a pair almost always land on the
    same chromosome, so they are not two independent observations. Sampling
    reads individually would therefore understate the spread of the ratio --
    by close to a factor of two in practice. Subsampling here draws pairs, as
    FASTQ-level subsampling does.

    A category is the pair of chromosome indices the two mates occupy, which
    is usually (i, i); a mate that is unmapped or filtered out gives -1 and
    the pair contributes a single read.
    """
    categories: list[tuple[int, int]]
    counts: np.ndarray
    names: list[str]
    lengths: dict[str, int]
    total_pairs: int
    total_reads: int

    @property
    def reads_per_pair(self) -> float:
        return self.total_reads / self.total_pairs if self.total_pairs else 0.0

    def chrom_totals(self) -> dict[str, int]:
        """Reads per chromosome over the whole pool, undoing the pairing."""
        tally = np.zeros(len(self.names), dtype=np.int64)
        for (a, b), k in zip(self.categories, self.counts):
            tally[a] += k
            if b >= 0:
                tally[b] += k
        return {name: int(tally[i]) for i, name in enumerate(self.names)}


def _count_pairs_one(job: tuple[str, str, int, list[str]]) -> dict[tuple[int, int], int]:
    """
    Tally pair categories on one chromosome. Module level, so it can be forked.

    Only the first mate of each pair is examined; the chromosome of the other
    is read from its mate fields, which avoids having to find it in a
    coordinate-sorted file.
    """
    bam_path, contig, min_mapq, names = job
    index_of = {n: i for i, n in enumerate(names)}
    tally: dict[tuple[int, int], int] = {}

    with pysam.AlignmentFile(bam_path, "rb") as bam:
        ref_name = {}
        for read in bam.fetch(contig=contig):
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            if read.is_qcfail or read.is_duplicate:
                continue
            if read.mapping_quality < min_mapq:
                continue
            # Count each pair once. Unpaired reads are counted on their own.
            if read.is_paired and not read.is_read1:
                continue

            here = index_of[contig]
            there = -1
            if read.is_paired and read.next_reference_id >= 0:
                rid = read.next_reference_id
                if rid not in ref_name:
                    ref_name[rid] = bam.get_reference_name(rid)
                there = index_of.get(ref_name[rid], -1)

            key = (here, there)
            tally[key] = tally.get(key, 0) + 1

    return tally


def collect_pool(bam_path: Path, names: list[str], min_mapq: int,
                 jobs: int = 1) -> Pool:
    """
    Read the BAM once and record how many pairs occupy each combination of
    chromosomes. Subsampling then draws pairs from those categories.
    """
    lengths: dict[str, int] = {}
    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        header_lengths = dict(zip(bam.references, bam.lengths))
        for name in names:
            if name not in header_lengths:
                sys.exit(f"[error] {name} not in {bam_path} header")
            lengths[name] = header_lengths[name]

    work = [(str(bam_path), name, min_mapq, names) for name in names]
    merged: dict[tuple[int, int], int] = {}
    if jobs > 1 and len(work) > 1:
        import multiprocessing as mp
        with mp.Pool(min(jobs, len(work))) as pool:
            parts = pool.map(_count_pairs_one, work)
    else:
        parts = [_count_pairs_one(job) for job in work]
    for part in parts:
        for key, n in part.items():
            merged[key] = merged.get(key, 0) + n

    categories = sorted(merged)
    counts = np.array([merged[k] for k in categories], dtype=np.int64)
    total_pairs = int(counts.sum())
    total_reads = int(sum(n * (2 if b >= 0 else 1)
                          for (a, b), n in zip(categories, counts)))

    return Pool(categories=categories, counts=counts, names=names,
                lengths=lengths, total_pairs=total_pairs,
                total_reads=total_reads)


def subsample_counts(pool: Pool, n_reads: int,
                     rng: np.random.Generator) -> dict[str, int]:
    """
    Draw pairs without replacement until roughly n_reads reads are covered,
    and tally the reads per chromosome.

    The target is expressed in mapped reads, because that is what karyoratio
    reports and what a results table records; it is converted to a pair count
    using the ratio observed in this pool. Only the tally is needed, so the
    per-category counts are drawn straight from the multivariate
    hypergeometric distribution rather than pair by pair.
    """
    rpp = pool.reads_per_pair
    n_pairs = int(round(n_reads / rpp)) if rpp else 0
    n_pairs = max(1, min(n_pairs, pool.total_pairs))

    drawn = rng.multivariate_hypergeometric(pool.counts, n_pairs)

    tally = np.zeros(len(pool.names), dtype=np.int64)
    for (a, b), k in zip(pool.categories, drawn):
        if k:
            tally[a] += k
            if b >= 0:
                tally[b] += k
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
    p.add_argument("--reference", type=Path, default=None,
                   help="indexed reference FASTA (samtools faidx); "
                        "not needed with --from-bam")
    p.add_argument("--chrom-file", required=True, type=Path,
                   help="karyoratio chromosome definition file")
    p.add_argument("--outdir", required=True, type=Path)

    p.add_argument("--from-bam", action="append", default=[], metavar="SEX=PATH",
                   help="validate on real samples of known sex instead of "
                        "simulating, e.g. --from-bam male=dog_m.bam "
                        "--from-bam female=dog_f.bam. Stages 1-3 are skipped; "
                        "the subsampling and the summary are identical. "
                        "Repeatable.")

    s = p.add_argument_group("simulation")
    s.add_argument("--profile-from", type=Path, default=None,
                   help="real FASTQ from which to derive an ART error profile; "
                        "without it ART's built-in HiSeq 2500 profile is used")
    s.add_argument("--read-length", type=int, default=150)
    s.add_argument("--art-system", default="HS25", choices=sorted(ART_SYSTEMS),
                   help="ART built-in error profile, used when --profile-from "
                        "is not given (default: HS25, max 150 bp)")
    s.add_argument("--total-reads", type=int, default=2_000_000,
                   help="reads to simulate per sex, counting both mates, so "
                        "half that many pairs (default: 2,000,000)")
    s.add_argument("--fragment-mean", type=int, default=300)
    s.add_argument("--fragment-sd", type=int, default=50)
    s.add_argument("--seed", type=int, default=1)

    a = p.add_argument_group("alignment")
    a.add_argument("--aligner", choices=("bowtie2", "bwa"), default="bowtie2")
    a.add_argument("--index", type=Path, default=None,
                   help="basename of a prebuilt aligner index, skipping the "
                        "build. Its sequence names must match those in "
                        "--chrom-file; check with `bowtie2-inspect -n`")
    a.add_argument("--threads", type=int, default=4,
                   help="threads for bowtie2/bwa and samtools (default: 4)")
    a.add_argument("--jobs", type=int, default=0,
                   help="processes for counting sequences in parallel; "
                        "0 means use --threads")
    a.add_argument("--min-mapq", type=int, default=30)
    a.add_argument("--serial-simulation", action="store_true",
                   help="simulate the two sexes one after the other; by "
                        "default they run concurrently, since ART is "
                        "single-threaded")

    c = p.add_argument_group("subsampling")
    c.add_argument("--read-counts", default=",".join(str(c) for c in DEFAULT_COUNTS),
                   help="mapped-read counts to test, comma-separated")
    c.add_argument("--replicates", type=int, default=30)

    k = p.add_argument_group("comparison")
    k.add_argument("--target", default="X",
                   help="roles or sequence names forming the numerator, "
                        "pooled (default: X)")
    k.add_argument("--ref-chroms", default="autosome",
                   help="roles or sequence names forming the denominator; "
                        "this is karyoratio's --reference, renamed here "
                        "because --reference is the genome FASTA "
                        "(default: autosome)")
    k.add_argument("--lower", type=float, default=0.60,
                   help="CI upper bound below this gives the first label")
    k.add_argument("--upper", type=float, default=0.80,
                   help="CI lower bound above this gives the third label")
    k.add_argument("--labels", default="XY,n.a.,XX",
                   help="three comma-separated labels: below, ambiguous, above")
    k.add_argument("--male-label", default="XY",
                   help="which label a male sample should receive")
    k.add_argument("--female-label", default="XX",
                   help="which label a female sample should receive")

    p.add_argument("--karyoratio", type=Path, default=None,
                   help="path to karyoratio.py (default: next to this script)")
    p.add_argument("--skip-existing", action="store_true",
                   help="reuse intermediate files from an earlier run")
    p.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    args = p.parse_args(argv)

    args.labels = tuple(s.strip() for s in args.labels.split(","))
    if len(args.labels) != 3:
        p.error("--labels needs exactly three comma-separated values")
    for role, lab in (("--male-label", args.male_label),
                      ("--female-label", args.female_label)):
        if lab not in args.labels:
            p.error(f"{role} {lab!r} is not among --labels "
                    f"({', '.join(args.labels)})")
    if args.male_label == args.female_label:
        p.error("--male-label and --female-label must differ")

    # ART refuses a read length its built-in profile cannot supply, and the
    # message it gives does not make the cause obvious.
    if not args.from_bam and not args.profile_from:
        cap = ART_SYSTEMS[args.art_system]
        if args.read_length > cap:
            suitable = sorted(s for s, m in ART_SYSTEMS.items()
                              if m >= args.read_length)
            p.error(
                f"--read-length {args.read_length} is more than ART's "
                f"{args.art_system} profile can generate ({cap} bp). "
                + (f"Either lower it to {cap}, pick a profile that reaches it "
                   f"({', '.join(suitable)}), " if suitable
                   else "Either lower it, ")
                + "or give --profile-from with real reads of the length you "
                  "want, which removes the limit."
            )
    return args


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    kr = load_karyoratio(args.karyoratio)

    if not args.from_bam:
        if args.reference is None:
            sys.exit("[error] --reference is required unless --from-bam is used")
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

    try:
        target = spec.resolve([s.strip() for s in args.target.split(",")])
        reference = spec.resolve([s.strip() for s in args.ref_chroms.split(",")])
    except KeyError as exc:
        sys.exit(f"[error] {exc}")
    overlap = set(target) & set(reference)
    if overlap:
        sys.exit("[error] sequence(s) in both target and reference: "
                 + ", ".join(sorted(overlap)))
    needed = target + [n for n in reference if n not in target]

    default_cmp = (args.target.strip() == "X"
                   and args.ref_chroms.strip() == "autosome"
                   and (args.lower, args.upper) == (0.60, 0.80))
    if not default_cmp:
        print("[warn] this is not the Rx comparison the 0.60/0.80 thresholds "
              "were validated for.\n"
              "       The point of calibrating is to establish thresholds for "
              "it, so the\n"
              "       summary below is only as good as the --lower/--upper "
              "you supplied.\n"
              "       Read the separation reported at the end, then re-run "
              "with thresholds\n"
              "       taken from it.", file=sys.stderr)

    counts = sorted(int(c) for c in args.read_counts.split(",") if c.strip())
    pools: dict[str, Pool] = {}
    # Each real sample is kept separate: pooling individuals would blend
    # coverage profiles that differ, and the per-sample spread is part of
    # what the validation is meant to reveal.
    sex_of: dict[str, str] = {}

    if args.from_bam:
        print("[mode] validating on real samples of known sex", file=sys.stderr)
        for item in args.from_bam:
            if "=" not in item:
                sys.exit(f"[error] --from-bam expects SEX=PATH, got {item!r}")
            label, path = item.split("=", 1)
            label = label.strip().lower()
            if label not in ("male", "female"):
                sys.exit(f"[error] --from-bam sex must be male or female, "
                         f"got {label!r}")
            bam = Path(path)
            if not bam.exists():
                sys.exit(f"[error] {bam} does not exist")
            key = bam.name
            for suffix in (".sorted.bam", ".bam"):
                if key.endswith(suffix):
                    key = key[: -len(suffix)]
                    break
            if key in pools:
                sys.exit(f"[error] two samples share the name {key!r}")
            print(f"[stage 4a] collecting reads for {key} ({label})",
                  file=sys.stderr)
            pools[key] = collect_pool(bam, needed, args.min_mapq,
                                      jobs=args.jobs or args.threads)
            sex_of[key] = label
            print(f"           {pools[key].total_reads:,} reads in "
                  f"{pools[key].total_pairs:,} pairs "
                  f"({pools[key].reads_per_pair:.3f} reads per pair)",
                  file=sys.stderr)
        if len(set(sex_of.values())) < 2:
            print("[note] only one sex supplied; the summary will cover it "
                  "alone, which\n"
                  "       shows the ambiguous rate but cannot demonstrate "
                  "separation.", file=sys.stderr)
    else:
        sex_of = {"female": "female", "male": "male"}
        female_fa, male_fa = build_genomes(
            args.reference, autosomes, x, y, args.outdir, args.skip_existing)

        profile = (build_profile(args.profile_from, args.outdir,
                                 args.skip_existing)
                   if args.profile_from else None)

        if args.index:
            index = args.index
            marker = Path(str(index) + (".1.bt2" if args.aligner == "bowtie2"
                                        else ".bwt"))
            if not marker.exists():
                sys.exit(f"[error] no {args.aligner} index at {index} "
                         f"(expected {marker})")
            print(f"[stage 3a] using the prebuilt index at {index}",
                  file=sys.stderr)
        else:
            index = build_index(args.reference, args.aligner, args.outdir,
                                args.threads, args.skip_existing)

        sim_args = dict(outdir=args.outdir, read_length=args.read_length,
                        total_reads=args.total_reads,
                        fragment_mean=args.fragment_mean,
                        fragment_sd=args.fragment_sd, profile=profile,
                        skip_existing=args.skip_existing)

        # ART is single-threaded, so the two sexes are simulated side by side.
        # Mapping is not, and gets every thread in turn.
        if args.serial_simulation:
            reads = {label: simulate(genome, label, seed=args.seed + i,
                                     art_system=args.art_system, **sim_args)
                     for i, (label, genome) in
                     enumerate((("female", female_fa), ("male", male_fa)))}
        else:
            from concurrent.futures import ThreadPoolExecutor
            print("[stage 2b] simulating both sexes concurrently", file=sys.stderr)
            with ThreadPoolExecutor(max_workers=2) as pool:
                futures = {
                    label: pool.submit(simulate, genome, label,
                                       seed=args.seed + i,
                                       art_system=args.art_system, **sim_args)
                    for i, (label, genome) in
                    enumerate((("female", female_fa), ("male", male_fa)))
                }
                reads = {label: fut.result() for label, fut in futures.items()}

        for label in ("female", "male"):
            r1, r2 = reads[label]
            bam = map_reads(r1, r2, index, label, args.aligner, args.outdir,
                            args.threads, args.skip_existing)
            print(f"[stage 4a] collecting reads for {label}", file=sys.stderr)
            pools[label] = collect_pool(bam, needed, args.min_mapq,
                                        jobs=args.jobs or args.threads)
            p = pools[label]
            # --total-reads counts both mates, so ART was asked for half as
            # many pairs as that number.
            gen_pairs = args.total_reads / 2
            print(f"           {p.total_reads:,} reads in {p.total_pairs:,} "
                  f"pairs; {p.total_pairs / gen_pairs:.1%} of generated pairs "
                  f"survive, giving {p.total_reads / gen_pairs:.3f} mapped "
                  f"reads per generated pair", file=sys.stderr)

    usable = [n for n in counts
              if all(n <= p.total_reads for p in pools.values())]
    dropped = [n for n in counts if n not in usable]
    if dropped:
        print(f"[note] skipping read counts above the simulated pool: "
              f"{', '.join(str(d) for d in dropped)}", file=sys.stderr)
    if not usable:
        sys.exit("[error] no requested read count fits in the available pool; "
                 + ("use samples with more reads, or lower --read-counts"
                    if args.from_bam else "raise --total-reads"))

    rng = np.random.default_rng(args.seed)
    per_rep = args.outdir / "replicates.tsv"
    print(f"[stage 4b] {len(usable)} read counts x {args.replicates} replicates "
          f"x 2 sexes", file=sys.stderr)

    rows = []
    with per_rep.open("w") as fh:
        fh.write("sample\tsex\tread_count\treplicate\tratio\tCI_l\tCI_u\t"
                 "rho\trho_p\tassignment\n")
        for sample_key, pool in pools.items():
            sex = sex_of[sample_key]
            expected = (args.female_label if sex == "female"
                        else args.male_label)
            for n in usable:
                for rep in range(1, args.replicates + 1):
                    counts_d = subsample_counts(pool, n, rng)
                    try:
                        res = kr.compute_ratio(
                            counts_d, pool.lengths, target, reference,
                            f"{sex}_{n}_{rep}",
                            lower=args.lower, upper=args.upper,
                            labels=args.labels)
                    except ValueError:
                        continue
                    fh.write(f"{sample_key}\t{sex}\t{n}\t{rep}\t"
                             f"{res.ratio:.6f}\t{res.CI_l:.6f}\t"
                             f"{res.CI_u:.6f}\t{res.rho:.6f}\t"
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
                amb = sum(1 for r in sel if r[3] == args.labels[1])
                wrong = len(sel) - correct - amb
                fh.write(f"{sex}\t{n}\t{len(sel)}\t{correct}\t{amb}\t{wrong}\t"
                         f"{100 * amb / len(sel):.2f}\n")

    if args.from_bam and len(pools) > 2:
        print("\n[result] full-pool ratio per sample", file=sys.stderr)
        for key, pool in pools.items():
            try:
                full = kr.compute_ratio(
                    pool.chrom_totals(),
                    pool.lengths, target, reference, key,
                    lower=args.lower, upper=args.upper, labels=args.labels)
                print(f"         {key:<28} {sex_of[key]:<7} "
                      f"ratio={full.ratio:.4f}  rho={full.rho:.3f}  "
                      f"{full.assignment}", file=sys.stderr)
            except ValueError as exc:
                print(f"         {key:<28} {exc}", file=sys.stderr)
        print("         Spread between samples of the same sex is the thing "
              "to look at: a\n"
              "         threshold has to clear the least favourable of them, "
              "not the mean.",
              file=sys.stderr)

    # Where the two sexes actually sit, which is what thresholds should be
    # read off. Needed whenever the comparison is not the calibrated one.
    try:
        import csv
        vals: dict[tuple[str, int], list[float]] = {}
        with per_rep.open() as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                vals.setdefault((row["sex"], int(row["read_count"])),
                                []).append(float(row["ratio"]))
        print("\n[result] observed separation (ratio range per sex)",
              file=sys.stderr)
        print(f"         {'reads':>9}  {'female':>21}  {'male':>21}",
              file=sys.stderr)
        for n in usable:
            f = vals.get(("female", n), [])
            m = vals.get(("male", n), [])
            fs = f"{min(f):.3f} - {max(f):.3f}" if f else "-"
            ms = f"{min(m):.3f} - {max(m):.3f}" if m else "-"
            print(f"         {n:>9,}  {fs:>21}  {ms:>21}", file=sys.stderr)
        if vals:
            top = max(usable)
            f, m = vals.get(("female", top), []), vals.get(("male", top), [])
            if f and m:
                lo, hi = (m, f) if max(m) < min(f) else (f, m)
                if max(lo) < min(hi):
                    print(f"         At {top:,} reads the groups do not "
                          f"overlap: {max(lo):.3f} < {min(hi):.3f}.",
                          file=sys.stderr)
                    print(f"         Thresholds for this comparison belong "
                          f"between them.", file=sys.stderr)
                else:
                    print(f"         At {top:,} reads the groups overlap; no "
                          f"threshold separates them.", file=sys.stderr)
    except Exception as exc:
        print(f"[note] could not summarise the separation: {exc}",
              file=sys.stderr)

    wrong_total = sum(1 for r in rows
                      if r[3] != r[2] and r[3] != args.labels[1])

    print(f"\n[result] per-replicate values: {per_rep}", file=sys.stderr)
    print(f"[result] summary:              {summary}", file=sys.stderr)
    print(f"[result] incorrect assignments: {wrong_total} of {len(rows)}",
          file=sys.stderr)
    # The first depth with no unclassified replicate is not a property of the
    # method: run more replicates and rare events surface, pushing it upward.
    # A rate is stable, so the depth at which the rate drops below a target is
    # reported instead.
    print("[result] depth at which the unclassified rate falls below:",
          file=sys.stderr)
    print(f"         {'sex':<8} {'<10%':>10} {'<5%':>10} {'<1%':>10}",
          file=sys.stderr)
    for sex in ("female", "male"):
        sel = [(n, sum(1 for r in rows
                       if r[0] == sex and r[1] == n and r[3] == args.labels[1])
                   / max(1, sum(1 for r in rows if r[0] == sex and r[1] == n))
                   * 100)
               for n in usable]
        sel = [(n, p) for n, p in sel if any(r[0] == sex and r[1] == n
                                             for r in rows)]
        if not sel:
            continue
        cells = []
        for target in (10, 5, 1):
            ans = "-"
            for n, pct in reversed(sel):
                if pct >= target:
                    break
                ans = f"{n:,}"
            cells.append(ans)
        print(f"         {sex:<8} {cells[0]:>10} {cells[1]:>10} {cells[2]:>10}",
              file=sys.stderr)
    if args.replicates < 50:
        print(f"         With {args.replicates} replicates a rate below "
              f"roughly {100/args.replicates:.0f}% cannot be\n"
              f"         resolved; raise --replicates before quoting the "
              f"<1% column.", file=sys.stderr)
    if wrong_total:
        bad = sorted({r[1] for r in rows
                      if r[3] != r[2] and r[3] != args.labels[1]})
        print(f"[warn] {wrong_total} incorrect assignment(s), at read counts "
              f"{', '.join(f'{b:,}' for b in bad)}.", file=sys.stderr)
        print(f"       The highest of these is {max(bad):,}. A threshold set "
              f"above that is unaffected; below it the thresholds do not "
              f"hold here.", file=sys.stderr)

    return 0


if __name__ == "__main__":
    sys.exit(main())
