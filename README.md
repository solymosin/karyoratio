# karyoratio

Chromosome coverage ratios computed directly from an indexed BAM file.

By default it computes $R_x$, the X-to-autosome read ratio of Mittnik et al.
(2016, *PLoS ONE* 11:e0163019), and estimates genetic sex from it. The
comparison itself is parameterised, so other karyotypic questions can be
expressed with the same machinery -- but the calibrated decision rule applies
to $R_x$ alone, and the script says so when you depart from it.

A companion script calibrates the thresholds for a new species by simulation,
and can put real samples of known sex through the identical procedure.

The method, the calibration and a worked application to stone marten and red
fox faecal metagenomes are described in the accompanying preprint:
[Sexing carnivores from faecal metagenomes: calibrating $R_x$ for stone marten
and red fox](https://www.biorxiv.org/content/10.64898/2026.11.04.764345v1)
(bioRxiv, 2026).

## Contents

| file | what it does |
|---|---|
| `karyoratio.py` | the estimator: ratio, interval, assignment |
| `calibrate.py` | calibration by simulation, or validation on real samples |
| `chromfile_from_report.py` | builds a chromosome file from an NCBI sequence report |
| `pairs_to_reads.py` | measures the generated-pair to mapped-read factor |
| `run_calibration.sh` | calibration wrapper for a local genome and index |
| `run_marten_fox_calibration.sh` | the calibration run reported in the preprint |
| `*_chroms.tsv` | chromosome definitions for the stone marten and red fox |

## Installation

```bash
pip install pysam numpy scipy
chmod +x karyoratio.py calibrate.py chromfile_from_report.py pairs_to_reads.py
```

Calibration additionally needs `art_illumina`, `bowtie2` (or `bwa`) and
`samtools` on PATH, and `art_profiler_illumina` if an error profile is derived
from real reads.

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

The sequence names must match the BAM header exactly. Note that SAM forbids
whitespace in a reference name, so an aligner writes only the part up to the
first space even when the index stores a full FASTA header; the chromosome
file has to carry that short form.

Two example files are provided: `vulvul3_chroms.tsv` (red fox) and
`mmarfoi2.1_chroms.tsv` (stone marten).

### Building the file

From a BAM, by length:

```bash
./karyoratio.py --bam sample.bam --make-chrom-file \
                --x OZ199613.1 --y OZ199623.1 --mt OZ199624.1 \
                --prefix OZ --min-autosome-length 10000000 \
                -o chroms.tsv
```

Or from an NCBI sequence report, which states the role of every sequence so
that nothing is inferred from names or lengths:

```bash
bowtie2-inspect -n my_index > names.txt
./chromfile_from_report.py --names names.txt \
    --report sequence_report.jsonl --expect-autosomes 38 -o chroms.tsv
```

The second form matches roles to whatever naming the index uses: RefSeq
accessions, GenBank accessions, UCSC-style names or bare chromosome numbers.
Either way, check the result; both are heuristics about an assembly, not
authorities on a karyotype.

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

Under the defaults `ratio` is $R_x$ and `assignment` follows the decision rule
of Mittnik et al., which amounts to a pair of one-sided tests: `XY` is returned
when $R_x \geq 0.60$ is rejected (`CI_u` < 0.60) and `XX` when
$R_x \leq 0.80$ is rejected (`CI_l` > 0.80), each at a nominal one-sided level
of 2.5 per cent; otherwise `n.a.`. `CI_l` and `CI_u` bound the nominal 95 per
cent confidence interval, $R_x \pm 1.96 \times SE$, with SE the standard error
of the mean over the autosomes, as in the original method; the interval itself
makes no probability statement about where a sample's ratio lies.

The nominal level is not attained. Every per-autosome term shares the same X
count, so the sampling noise of that count, which dominates at low depth, never
enters the SE. In simulation the spread of $R_x$ between replicates was roughly
2.5 to 8 times the reported SE at depths up to 100,000 read pairs. The error
rates of the rule are therefore those established by calibration (below), not
the nominal ones.

`rho` is a validity criterion, not part of the estimate. Read counts should
scale with sequence length; a low value means the mapping is too sparse or too
uneven for the ratio to be meaningful. Male samples give systematically lower
values than females, because one X chromosome departs from the expected
length-count relationship.

`example_output.tsv` shows the tool run on seven faecal metagenomes.

## How the ratio is computed

The target sequences are **pooled**: their reads and their lengths are summed,
so the ratio refers to one unit rather than to a mean over several. The pooled
target density is then compared with the density of each reference sequence
separately,

```
ratio_i = (m_T / s_T) / (m_i / s_i)
```

and the reported value is the mean of `ratio_i` across the reference set.

The X chromosome carries a slightly smaller share of uniquely mappable sequence
than the autosomes, so a MAPQ filter lowers the ratio below the value the
karyotype prescribes. In simulated female genomes of the stone marten and the
red fox the ratio is 1.00 without a threshold and 0.96 at MAPQ 30, and it falls
further as the threshold rises. Calibrate at the threshold you intend to use,
and read samples against the calibrated centre rather than against 1.00.

The normalising totals of the original formulation cancel in this quotient, so
**a sequence present in neither set has no effect on the result**. What matters
is membership of the reference set: a sequence there contributes one term to
the mean. This is why the mitochondrion must be kept out. Its copy number per
cell is orders of magnitude above that of the nuclear genome, so as a reference
term it drags the mean down sharply -- not because it distorts a normalisation,
but because it is a denominator that does not belong.

The Y chromosome is excluded from the defaults for a different reason: its
coverage depends on the very thing the method is trying to determine.

## Depth, and when not to trust an assignment

The 0.60/0.80 rule separates the sexes cleanly once coverage is adequate. Below
that, the male and female distributions overlap, and the interval can
clear a threshold on a sample that carries far too few reads to support any
assignment. Applying the rule on its own will then return a confident-looking
call that the data do not justify.

Two guards are available, both inactive by default:

```bash
./karyoratio.py --bam sample.bam --chrom-file chroms.tsv --min-reads 7500
```

**`--min-reads` is the one to set.** Every incorrect assignment we have
observed in simulation came from a sample with too few reads rather than from
an unusual distribution of them.

`--min-rho` returns the ambiguous label below a given length-count correlation.
It is a diagnostic rather than a primary filter: it catches a sample that has
ample reads distributed oddly across the chromosomes -- uneven coverage, or a
mis-specified chromosome file -- but it is an expensive way to guard against
low depth.

There is no universal default for either, because the usable minimum depends on
the species. Establish it by calibration and set it explicitly.

## Calibrating for a new species

```bash
./calibrate.py --reference genome.fa --chrom-file chroms.tsv \
               --index bowtie2_index --outdir calib \
               --profile-from real_R1.fq.gz \
               --replicates 200 --threads 16
```

It builds artificial female and male genomes from the reference at the copy
number each sex carries, simulates reads from both with ART, maps each set to
the unmodified reference, then subsamples the mapped reads and computes the
ratio for each replicate.

Mapping once per sex and subsampling afterwards is what makes this affordable:
two alignments instead of one per replicate. **Subsampling draws pairs, not
individual reads**, because the two mates of a pair almost always land on the
same chromosome and are therefore not two independent observations; sampling
reads singly understates the spread by close to a factor of two. The target is
expressed in mapped reads, which is what `karyoratio` reports, and converted
internally using the pair-to-read ratio observed in the pool.

The index should be built from the host genome alone. Competitive mapping
against a host-plus-human reference belongs to the analysis of real samples,
where it absorbs contaminating reads; it has no place in a calibration, where
every simulated read is host-derived.

### Reading the result

`replicates.tsv` holds every individual result; `summary.tsv` gives, per depth
and sex, how many replicates were assigned correctly, how many were left
unclassified, and how many were assigned to the wrong sex.

The script reports the depth at which the unclassified **rate** falls below 10,
5 and 1 per cent. It does not report the first depth with no unclassified
replicate, because that is not a property of the method: run more replicates
and rare events surface, pushing it upward without anything having changed.
Resolving a rate near 1 per cent needs on the order of 200 replicates.

An incorrect assignment is a different matter from an unclassified one, and the
depths at which any occurred are listed separately.

### Validating on real samples

Simulation establishes where the thresholds break down under ideal conditions.
Real samples of known sex carry everything it leaves out: reference bias,
repeat structure, contamination, and above all variation between individuals.
The same machinery handles both:

```bash
./calibrate.py --from-bam male=m1.bam --from-bam female=f1.bam \
               --from-bam male=m2.bam --from-bam female=f2.bam \
               --chrom-file chroms.tsv --outdir downsample
```

Stages 1--3 are skipped and each sample is kept separate, since the spread
between individuals is part of what the validation is meant to reveal.

### What the calibration found

Two species, 200 replicates per level, depth in mapped reads:

| | <10% | <5% | <1% |
|---|---|---|---|
| stone marten, simulated | 4,000--5,000 | 7,500 | 15,000 |
| red fox, simulated | 4,000--5,000 | 7,500 | 10,000--15,000 |

### Converting between generated pairs and mapped reads

A simulation that subsamples FASTQ files counts generated read pairs; what
reaches the ratio is mapped reads. The factor is stable but has to be measured:

```bash
./pairs_to_reads.py --table sim_results.csv --pairs read_count --mapped s
```

For the stone marten it was 1.944--1.946 and for the red fox 1.921--1.923
mapped reads per generated pair. A figure axis in one unit cannot be compared
with a table column in the other.

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

## Other comparisons

`--target` and `--reference` take roles or literal sequence names. The three
`--labels` and the two thresholds are likewise free, so a different karyotypic
question can be expressed directly:

```bash
# ZW system: the female is the heterogametic sex, so the labels swap
./karyoratio.py --bam bird.bam --chrom-file bird_chroms.tsv \
                --target Z --reference autosome --labels ZW,n.a.,ZZ
```

> **The thresholds do not generalise.** The 0.60/0.80 rule was established for
> $R_x$ in diploid XY systems, and nothing about it transfers to another
> comparison. The script computes the ratio and its interval as defined
> whatever you ask of it, and warns whenever the comparison or the thresholds
> depart from the defaults. `calibrate.py` accepts the same `--target`,
> `--ref-chroms`, `--lower`, `--upper` and `--labels`, and reports where the
> two groups actually sit so thresholds can be read off rather than guessed.

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
reference genomes lack, including the red fox -- in normalising against the
full autosome set rather than a single chromosome, in reporting an
interval and a validity criterion rather than a point estimate, and in shipping
the means to calibrate itself for a species it has not seen.

## Citation

If you use this software, please cite Mittnik et al. (2016) for the method and
the preprint below for the implementation and its calibration.

```
Mittnik A, Wang C-C, Svoboda J, Krause J (2016) A molecular approach to the
sexing of the triple burial at the Upper Paleolithic site of Dolní Věstonice.
PLoS ONE 11(10): e0163019.
```

```
Vujić V, Mayer C, Helyes G, Maróti-Agóts Á, Maróti G, Csabai I, Solymosi N (2026)
Sexing carnivores from faecal metagenomes: calibrating Rx for stone marten and
red fox.
bioRxiv 2026.11.04.764345. https://www.biorxiv.org/content/10.64898/2026.11.04.764345v1
```
