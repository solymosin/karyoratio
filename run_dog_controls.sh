#!/usr/bin/env bash
#
# run_dog_controls.sh -- validate karyoratio on real dog samples of known sex.
#
# Downloads two low-coverage (~1x) whole genome sequencing runs from the
# Darwin's Ark project, maps them to the dog reference, and checks both the
# sex assignment at full depth and the behaviour under downsampling.
#
#   SRR30353844  male    (BioSample SAMN43320406, sex = male)
#   SRR30353845  female  (verify the BioSample annotation before relying on it)
#
# Real controls complement the simulation rather than replacing it. The
# simulation establishes where the thresholds break down under ideal
# conditions; real samples carry everything the simulation leaves out --
# reference bias, repeat structure, contamination, library artefacts. An
# estimator that works on both is on firmer ground than one validated on
# either alone.
#
# Run run_dog_calibration.sh first: this script reuses its reference,
# chromosome file and bowtie2 index.
#
# Usage:
#   ./run_dog_controls.sh [workdir]
#
set -euo pipefail

WORKDIR="${1:-dog_calibration}"
KARYODIR="${KARYODIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
THREADS="${THREADS:-8}"
MIN_MAPQ="${MIN_MAPQ:-30}"
REPLICATES="${REPLICATES:-30}"
READ_COUNTS="${READ_COUNTS:-1000,1500,2000,2500,3000,5000,7500,10000,25000,50000,100000,250000}"

MALE_RUN="${MALE_RUN:-SRR30353844}"
FEMALE_RUN="${FEMALE_RUN:-SRR30353845}"

echo "== checking dependencies =="
MISSING=()
for tool in fasterq-dump samtools bowtie2 python3; do
  command -v "$tool" >/dev/null 2>&1 || MISSING+=("$tool")
done
if [[ ${#MISSING[@]} -gt 0 ]]; then
  echo "missing from PATH: ${MISSING[*]}" >&2
  echo "  conda install -c bioconda sra-tools samtools bowtie2" >&2
  exit 1
fi
echo "all present"

cd "$WORKDIR"

REF="reference.fa"
CHROMS="canfam_chroms.tsv"
INDEX="calib/ref_index"

for f in "$REF" "$CHROMS"; do
  [[ -f "$f" ]] || { echo "not found: $WORKDIR/$f -- run run_dog_calibration.sh first" >&2; exit 1; }
done
[[ -f "${INDEX}.1.bt2" ]] || {
  echo "no bowtie2 index at $WORKDIR/$INDEX -- run run_dog_calibration.sh first" >&2
  exit 1
}

mkdir -p controls
cd controls

# --------------------------------------------------------------------------- #
# 1. Confirm the annotated sex, rather than trusting the accession
# --------------------------------------------------------------------------- #

if command -v esearch >/dev/null 2>&1; then
  echo
  echo "== BioSample sex annotation =="
  for run in "$MALE_RUN" "$FEMALE_RUN"; do
    printf '  %s: ' "$run"
    esearch -db sra -query "$run" 2>/dev/null \
      | elink -target biosample 2>/dev/null \
      | efetch -format xml 2>/dev/null \
      | xtract -pattern BioSample \
          -block Attribute -if Attribute@attribute_name -equals sex \
            -element Attribute 2>/dev/null \
      || echo "could not retrieve"
    echo
  done
  echo "  If either disagrees with the label used below, stop and fix it."
else
  echo
  echo "[note] Entrez Direct not installed, so the BioSample sex annotation is"
  echo "       not being checked. Verify it by hand before trusting the result."
fi

# --------------------------------------------------------------------------- #
# 2. Download and map
# --------------------------------------------------------------------------- #

map_run () {
  local run="$1" label="$2"
  local bam="dog_${label}.bam"

  if [[ -f "$bam.bai" ]]; then
    echo "== $label ($run): alignment already present, skipping =="
    return
  fi

  if [[ ! -f "${run}_1.fastq" ]]; then
    echo "== $label ($run): downloading, roughly 2.6 Gbases =="
    fasterq-dump --split-files --threads "$THREADS" --progress "$run"
  fi

  echo "== $label ($run): mapping =="
  bowtie2 --very-sensitive-local -p "$THREADS" \
      -x "../$INDEX" -1 "${run}_1.fastq" -2 "${run}_2.fastq" \
      2> "log_bowtie2_${label}.txt" \
    | samtools sort -@ "$THREADS" -o "$bam" -
  samtools index "$bam"

  echo "   alignment rate: $(grep 'overall alignment rate' "log_bowtie2_${label}.txt" || echo 'see log')"

  # The FASTQs are large and no longer needed once the BAM exists.
  rm -f "${run}_1.fastq" "${run}_2.fastq"
}

map_run "$MALE_RUN" male
map_run "$FEMALE_RUN" female

# --------------------------------------------------------------------------- #
# 3. Assignment at full depth
# --------------------------------------------------------------------------- #

echo
echo "== assignment at full depth =="
python3 "$KARYODIR/karyoratio.py" \
    --bam dog_male.bam dog_female.bam \
    --chrom-file "../$CHROMS" \
    --min-mapq "$MIN_MAPQ" \
    --per-chrom per_chrom.tsv \
    -o full_depth.tsv
column -t -s$'\t' full_depth.tsv 2>/dev/null || cat full_depth.tsv

echo
echo "   dog_male should be XY and dog_female XX. Anything else means the"
echo "   estimator does not work on this species as configured, and nothing"
echo "   below is worth reading until that is resolved."

# --------------------------------------------------------------------------- #
# 4. Behaviour under downsampling, on real data
# --------------------------------------------------------------------------- #

echo
echo "== downsampling the real samples =="
echo "   $REPLICATES replicates at each of: $READ_COUNTS"
echo

python3 "$KARYODIR/calibrate.py" \
    --from-bam male=dog_male.bam \
    --from-bam female=dog_female.bam \
    --chrom-file "../$CHROMS" \
    --karyoratio "$KARYODIR/karyoratio.py" \
    --outdir downsample \
    --read-counts "$READ_COUNTS" \
    --replicates "$REPLICATES" \
    --min-mapq "$MIN_MAPQ"

echo
echo "== summary =="
column -t -s$'\t' downsample/summary.tsv 2>/dev/null || cat downsample/summary.tsv

cat <<EOF

Files written under $(pwd):

  full_depth.tsv           the assignment at full coverage
  per_chrom.tsv            per-chromosome read counts, useful for spotting a
                           mis-specified chromosome file
  downsample/replicates.tsv  every individual result under downsampling
  downsample/summary.tsv     correct / ambiguous / wrong per read count

Comparing this with the simulation:

  The two summary tables answer the same question from different directions.
  Where the real samples need more reads than the simulated ones to reach a
  clean call, the gap is what the simulation failed to model -- mapping bias,
  repeats, contamination. Take the larger of the two minima as --min-reads.

  A real sample assigned to the wrong sex at any depth is more serious than an
  ambiguous one, and worth tracing through per_chrom.tsv before going further.
EOF
