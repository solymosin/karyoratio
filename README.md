# karyoratio

Chromosome coverage ratios computed directly from an indexed BAM file.

By default it computes $R_x$, the X-to-autosome read ratio of Mittnik et al.
(2016, *PLoS ONE* 11:e0163019), and estimates genetic sex from it. The
comparison itself is parameterised, so other karyotypic questions can be
expressed with the same machinery -- but the calibrated decision rule applies
to $R_x$ alone, and the script says so when you depart from it.

## Installation

```bash
pip install pysam numpy scipy
chmod +x karyoratio.py
```

## Usage

```bash
# default: Rx and a sex assignment
./karyoratio.py --bam sample.bam --chrom-file chroms.tsv --min-mapq 30

# several samples, with per-chromosome counts kept
./karyoratio.py --bam *.bam --chrom-file chroms.tsv \
                -o rx.tsv --per-chrom rx_per_chrom.tsv
```

### Chromosome definition file

Two whitespace-separated columns, `<sequence name>` and `<role>`. Roles are
free-form labels; four are canonical and understood by the defaults:
`autosome`, `X`, `Y`, `MT`. Lines starting with `#` are ignored, as is
anything after the second column.

```
NC_132780.1	autosome
NC_132781.1	autosome
NC_132796.1	X
```

The sequence names must match the BAM header exactly.

Two example files are provided: `vulvul3_chroms.tsv` (red fox, complete) and
`mmarfoi2.1_chroms.tsv` (stone marten, sex chromosomes and mitochondrion only).

### Generating the file from a BAM

Rather than transcribing accessions by hand, build the file from the header of
a BAM mapped to the reference in question:

```bash
./karyoratio.py --bam marten_sample.bam --make-chrom-file \
                --x OZ199613.1 --y OZ199623.1 --mt OZ199624.1 \
                --prefix OZ --min-autosome-length 10000000 \
                -o mmarfoi2.1_chroms.tsv
```

Every sequence at or above `--min-autosome-length` becomes an autosome;
shorter ones are listed at the end as comments so they can be inspected and
promoted by hand if any belong. `--prefix` restricts consideration to names
with a given prefix, which keeps the contigs of a co-mapped genome out when a
composite reference was used.

**Check the result before using it.** The criterion is a length cutoff, which
works when the assembly is chromosome-length and the unplaced scaffolds are
clearly smaller, but it is a heuristic, not an authority on the karyotype.

## Output

One row per BAM:

| column | meaning |
|---|---|
| `sample` | BAM file name with the extension stripped |
| `n_target` | target sequences carrying reads (pooled into one unit) |
| `n_reference` | reference sequences carrying reads |
| `n_mapped` | reads passing the filters across both sets |
| `ratio` | mean of the pooled target density over each reference density |
| `CI_l`, `CI_u` | $ratio \pm 1.96 \times SE$ |
| `rho`, `rho_p` | Spearman correlation of sequence length against read count |
| `assignment` | one of the three `--labels`, default `XY` / `n.a.` / `XX` |

Under the defaults `ratio` is $R_x$ and `assignment` applies the thresholds of
Mittnik et al. to the CI bounds: `XY` if the upper bound is below 0.60, `XX`
if the lower bound is above 0.80, `n.a.` otherwise.

`rho` is a validity criterion, not part of the estimate. Read counts should
scale with sequence length; a low value means the mapping is too sparse or too
uneven for the ratio to be meaningful. Male samples give systematically lower
values than females, because one X chromosome departs from the expected
length-count relationship.

## How the ratio is computed

The target sequences are **pooled**: their reads and their lengths are summed,
so the ratio refers to one unit rather than to a mean over several. The pooled
target density is then compared with the density of each reference sequence
separately,

```
ratio_i = (m_T / s_T) / (m_i / s_i)
```

and the reported value is the mean of `ratio_i` across the reference set.

The normalising totals of the original formulation cancel in this quotient, so
**a sequence present in neither set has no effect on the result**. What matters
is membership of the reference set: a sequence there contributes one term to
the mean. This is why the mitochondrion must be kept out. Its copy number per
cell is orders of magnitude above that of the nuclear genome, so as a reference
term it drags the mean down sharply -- not because it distorts a normalisation,
but because it is a denominator that does not belong.

The Y chromosome is excluded from the defaults for a different reason: its
coverage depends on the very thing the method is trying to determine.

## Other comparisons

`--target` and `--reference` take roles or literal sequence names. The three
`--labels` and the two thresholds are likewise free, so a different karyotypic
question can be expressed directly:

```bash
# ZW system: the female is the heterogametic sex, so the labels swap
./karyoratio.py --bam bird.bam --chrom-file bird_chroms.tsv \
                --target Z --reference autosome \
                --labels ZW,n.a.,ZZ

# relative copy number of one chromosome against the rest
./karyoratio.py --bam sample.bam --chrom-file chroms.tsv \
                --target chr21 --reference autosome \
                --lower 1.25 --upper 1.40 \
                --labels disomic,n.a.,trisomic
```

> **The thresholds do not generalise.** The 0.60/0.80 rule was established for
> $R_x$ in diploid XY systems, and nothing about it transfers to another
> comparison. The same holds for any `--min-rho` or `--min-reads` cutoff
> derived for $R_x$: the correlation criterion in particular behaves
> differently between the sexes, and its behaviour under a different comparison
> is simply unknown. The script computes the ratio and its interval as defined
> whatever you ask of it, and prints a warning whenever the comparison or the
> thresholds depart from the defaults. Treat the `assignment` column as
> meaningful only after calibrating.

## The thresholds alone are not enough

The 0.60/0.80 rule separates the sexes cleanly once coverage is adequate. Below
that, the male and female distributions overlap, and a confidence interval can
clear a threshold on a sample that carries far too few reads to support any
assignment. Applying the rule on its own will then return a confident-looking
call that the data do not justify.

Two guards are available, both inactive by default:

```bash
./karyoratio.py --bam sample.bam --chrom-file chroms.tsv \
                --min-reads 5800 --min-rho 0.90
```

`--min-reads` returns the ambiguous label below a given number of mapped reads,
whatever the interval says. `--min-rho` does the same below a given
length-count correlation.

There is no universal default for either, because the usable minimum depends on
the chromosome count and size distribution of the species. Establish it by
calibration and then set it explicitly.

Note that `--min-rho` behaves asymmetrically between the sexes: male samples
give systematically lower correlations than females at the same coverage. A
threshold strict enough to exclude poor female samples may discard good male
ones.

## Working from a BAM rather than from idxstats

`samtools idxstats` counts every mapped read regardless of mapping quality.
Counting from the BAM allows an explicit MAPQ filter, which matters when the
reference is a composite genome (for example host plus human, used to absorb
contaminating reads competitively). Reads flagged as secondary, supplementary,
QC-fail or duplicate are excluded by default.

Orphan reads -- paired reads whose mate did not map properly -- **are** counted
by default. In low-coverage data they carry real signal, and excluding them can
silently drop the only reads covering a region. Use `--proper-pairs-only` to
exclude them.

## Calibrating for a new species

The thresholds of Mittnik et al. were established on human data. Before relying
on them for another species, two things need establishing: whether they separate
the sexes at all, and how many reads a sample must carry for the separation to
hold. `calibrate.py` answers both by simulation.

```bash
./calibrate.py --reference mMarFoi2.1.fa \
               --chrom-file mmarfoi2.1_chroms.tsv \
               --outdir calib_marten \
               --profile-from real_R1.fq.gz \
               --replicates 30 --threads 8
```

It builds artificial female and male genomes from the reference at the copy
number each sex carries, simulates reads from both with ART, maps each set to
the unmodified reference, then subsamples the mapped reads in memory at a range
of read counts and computes the ratio for each replicate.

Subsampling after mapping rather than before is what makes this affordable: two
alignments instead of one per replicate. It also fixes the read count to
*mapped* reads, which is what `karyoratio` reports, removing the ambiguity
between generated and mapped counts that otherwise creeps into these
calibrations.

Two files are written. `replicates.tsv` holds every individual result, for
plotting the ratio against read count by sex. `summary.tsv` gives, per read
count and sex, how many replicates were assigned correctly, how many were left
ambiguous, and how many were assigned to the wrong sex. The script prints the
lowest read count at which every replicate was classified correctly, which is
the figure to pass to `--min-reads`.

Requires `art_illumina`, `bowtie2` (or `bwa`) and `samtools` on PATH;
`art_profiler_illumina` only when `--profile-from` is used. Without a profile,
ART's built-in HiSeq 2500 model is used, which is adequate for a first pass but
worth replacing with a profile from the real data before drawing conclusions.

### Interpreting the result

Two outcomes matter. If any replicate is assigned to the wrong sex, the
thresholds do not transfer and need adjusting; the script says so explicitly.
If none is, the thresholds hold, and the remaining question is the read count
below which replicates start falling into the ambiguous interval.

In our own calibration for *Martes foina* and *Vulpes vulpes* the thresholds
never produced an incorrect assignment, but below a species-specific read count
a proportion of samples fell into the 0.60--0.80 interval and could not be
classified.

Note that a male genome is constructed with a single X alongside the doubled
autosomes, and a Y only if the reference has one. The X-to-autosome ratio does
not use the Y, so its absence -- as in the red fox and dog references -- costs
nothing here.

## Related tools

Two existing tools address the same question for human data:

* [SexChrLab/SexInference](https://github.com/SexChrLab/SexInference) computes
  X/chr19, Y/chr19 and Y/X read ratios with fixed thresholds. The chromosome
  names are hard-coded, a Y chromosome is required, and normalisation is
  against a single reference chromosome.
* [wdecoster's `infer_gender.py`](https://gist.github.com/wdecoster/a4182f1e200bc6eb8f2f1643678aae53)
  takes the Y/X ratio outside the pseudoautosomal regions, with human
  coordinates hard-coded.

Both are human-specific and depend on a Y chromosome. `karyoratio` differs in
being species-agnostic, in requiring no Y chromosome -- which many non-model
reference genomes lack, including the red fox and the dog -- in normalising
against the full autosome set rather than a single chromosome, and in reporting
a confidence interval and a validity criterion rather than a point estimate.

## Citation

If you use this script, please cite Mittnik et al. (2016) for the method.

```
Mittnik A, Wang C-C, Svoboda J, Krause J (2016) A molecular approach to the
sexing of the triple burial at the Upper Paleolithic site of Dolní Věstonice.
PLoS ONE 11(10): e0163019.
```
