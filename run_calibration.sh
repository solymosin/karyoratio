#!/usr/bin/env bash
#
# run_calibration.sh -- calibrate karyoratio thresholds for a species whose
# reference and aligner index you already have locally.
#
# The one decision worth getting right is which index to map the simulated
# reads to. If the real samples were mapped competitively against a composite
# reference -- host plus human, to absorb contaminating reads -- then the
# calibration has to use that same index, or it measures a pipeline you never
# ran. The simulated reads still come from the host genome alone; only the
# target of the alignment is composite.
#
# Usage:
#   REF=... INDEX=... CHROMS=... OUTDIR=... ./run_calibration.sh
#
# Required:
#   REF      host reference FASTA, indexed with samtools faidx. Host only --
#            the artificial male and female genomes are built from it, and a
#            composite FASTA would put human sequence into them.
#   INDEX    bowtie2 index basename. Composite, if that is what the real
#            samples were mapped against.
#   CHROMS   karyoratio chromosome file listing the host chromosomes only.
#   OUTDIR   where to write the results.
#
# Optional:
#   PROFILE_FROM   real FASTQ for an ART error profile (recommended)
#   READ_COUNTS    mapped-read counts to test
#   REPLICATES     default 30
#   READ_LENGTH    default 150; above that ART needs PROFILE_FROM
#   TOTAL_READS    read pairs simulated per sex, default 2000000
#   THREADS        default 8
#
# Examples:
#
#   # stone marten, against the marten+human index the samples were mapped to
#   REF=mMarFoi2.1.fa \
#   INDEX=idx/marten_human \
#   CHROMS=mmarfoi2.1_chroms.tsv \
#   OUTDIR=calib_marten \
#   PROFILE_FROM=sample1_R1.fq.gz READ_LENGTH=151 THREADS=16 \
#     ./run_calibration.sh
#
#   # red fox
#   REF=VulVul3.fa \
#   INDEX=idx/fox_human \
#   CHROMS=vulvul3_chroms.tsv \
#   OUTDIR=calib_fox \
#   PROFILE_FROM=sample2_R1.fq.gz READ_LENGTH=151 THREADS=16 \
#     ./run_calibration.sh
#
set -euo pipefail

KARYODIR="${KARYODIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"

REF="${REF:-}"
INDEX="${INDEX:-}"
CHROMS="${CHROMS:-}"
OUTDIR="${OUTDIR:-}"

THREADS="${THREADS:-8}"
REPLICATES="${REPLICATES:-30}"
TOTAL_READS="${TOTAL_READS:-2000000}"
READ_LENGTH="${READ_LENGTH:-150}"
READ_COUNTS="${READ_COUNTS:-250,500,750,1000,1500,2000,2500,3000,5000,7500,10000,25000,50000,100000,250000}"
PROFILE_FROM="${PROFILE_FROM:-}"
MIN_MAPQ="${MIN_MAPQ:-30}"

usage () {
  sed -n '3,50p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'
  exit 1
}

for v in REF INDEX CHROMS OUTDIR; do
  [[ -n "${!v}" ]] || { echo "missing: $v" >&2; echo >&2; usage; }
done

echo "== checking inputs =="
MISSING=()
for tool in art_illumina bowtie2 samtools python3; do
  command -v "$tool" >/dev/null 2>&1 || MISSING+=("$tool")
done
[[ -n "$PROFILE_FROM" ]] && {
  command -v art_profiler_illumina >/dev/null 2>&1 || MISSING+=("art_profiler_illumina")
}
if [[ ${#MISSING[@]} -gt 0 ]]; then
  echo "missing from PATH: ${MISSING[*]}" >&2
  echo "  conda install -c bioconda samtools art bowtie2" >&2
  exit 1
fi

for f in karyoratio.py calibrate.py; do
  [[ -f "$KARYODIR/$f" ]] || { echo "not found: $KARYODIR/$f" >&2; exit 1; }
done
[[ -f "$REF" ]] || { echo "not found: $REF" >&2; exit 1; }
[[ -f "$REF.fai" ]] || { echo "not indexed: run samtools faidx $REF" >&2; exit 1; }
[[ -f "$CHROMS" ]] || { echo "not found: $CHROMS" >&2; exit 1; }
[[ -f "${INDEX}.1.bt2" || -f "${INDEX}.1.bt2l" ]] \
  || { echo "no bowtie2 index at $INDEX" >&2; exit 1; }

# Every chromosome named in the file has to exist in the index, or the counts
# will silently come back empty.
echo "== checking that the chromosome file matches the index =="
bowtie2-inspect -n "$INDEX" | awk '{print $1}' | sort -u > "/tmp/.kr_idx_$$"
MISS=$(awk '$1 !~ /^#/ && NF>=2 {print $1}' "$CHROMS" \
       | sort -u | comm -23 - "/tmp/.kr_idx_$$" | head -5)
rm -f "/tmp/.kr_idx_$$"
if [[ -n "$MISS" ]]; then
  echo "these sequences from $CHROMS are not in the index:" >&2
  echo "$MISS" >&2
  echo "The chromosome file and the index must use the same names." >&2
  exit 1
fi
echo "   all chromosomes present"

# The artificial genomes are built from REF, so it must not already contain
# the second genome of a composite reference.
HOST_N=$(grep -c '^>' "$REF" || true)
IDX_N=$(bowtie2-inspect -n "$INDEX" | wc -l)
echo "   $REF has $HOST_N sequences; the index has $IDX_N"
if [[ "$IDX_N" -gt "$HOST_N" ]]; then
  echo "   (the index is larger, so it looks composite -- which is right if"
  echo "    that is what the real samples were mapped against)"
fi

mkdir -p "$OUTDIR"

echo
echo "== running the calibration =="
echo "   replicates : $REPLICATES"
echo "   read counts: $READ_COUNTS"
echo "   simulated  : $TOTAL_READS pairs per sex, $READ_LENGTH bp"
if [[ -n "$PROFILE_FROM" ]]; then
  echo "   error profile from $PROFILE_FROM"
else
  echo "   ART's built-in HiSeq 2500 model (max 150 bp)"
fi
echo

PROFILE_ARG=()
[[ -n "$PROFILE_FROM" ]] && PROFILE_ARG=(--profile-from "$PROFILE_FROM")

python3 "$KARYODIR/calibrate.py" \
    --reference "$REF" \
    --chrom-file "$CHROMS" \
    --index "$INDEX" \
    --outdir "$OUTDIR" \
    --karyoratio "$KARYODIR/karyoratio.py" \
    --read-counts "$READ_COUNTS" \
    --replicates "$REPLICATES" \
    --total-reads "$TOTAL_READS" \
    --read-length "$READ_LENGTH" \
    --min-mapq "$MIN_MAPQ" \
    --threads "$THREADS" \
    --skip-existing \
    "${PROFILE_ARG[@]}"

echo
echo "== summary =="
column -t -s$'\t' "$OUTDIR/summary.tsv" 2>/dev/null || cat "$OUTDIR/summary.tsv"

cat <<EOF

Written to $OUTDIR:

  replicates.tsv   every individual result
  summary.tsv      correct / ambiguous / wrong per read count and sex

What this does and does not measure:

  Every replicate is drawn from one simulated genome, so the spread reported
  here is the spread of read sampling alone. Variation between individuals --
  how far a given animal sits from the reference -- is not in it, and on real
  samples of another species it was not small. Treat the minimum read count
  from this as a floor rather than as a sufficient condition.
EOF
