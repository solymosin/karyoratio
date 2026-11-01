#!/usr/bin/env bash
#
# run_dog_calibration.sh -- calibrate karyoratio thresholds for the dog.
#
# Uses the prebuilt bowtie2 index of ROS_Cfam_1.0 published in the bowtie2
# index zone, which removes the only expensive step: building an index for a
# 2.4 Gb genome takes one to three hours, downloading one takes minutes.
#
# The reference FASTA is reconstructed from the index with bowtie2-inspect
# rather than downloaded separately. This guarantees that the sequence names
# used to simulate reads are exactly those the aligner reports, which is the
# usual way a prebuilt index goes wrong.
#
# The dog karyotype is 2n = 78: 38 autosomes plus X and Y. ROS_Cfam_1.0 was
# assembled from a male and does carry a Y, unlike the older CanFam3.1 and
# unlike the red fox reference. That makes no difference to the calibration:
# the X-to-autosome ratio does not use the Y, which is excluded from the
# comparison by default, though it is included in the simulated male genome.
#
# Usage:
#   ./run_dog_calibration.sh [workdir]
#
# SMOKE=1        fast end-to-end check before committing to the full run
# BUILD_INDEX=1  download the genome from NCBI and build the index locally
#
set -euo pipefail

WORKDIR="${1:-dog_calibration}"
ACCESSION="GCF_014441545.1"
ASSEMBLY="ROS_Cfam_1.0"
INDEX_URL="https://genome-idx.s3.amazonaws.com/bt/ROS_Cfam_1.0.zip"

KARYODIR="${KARYODIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"

THREADS="${THREADS:-8}"
REPLICATES="${REPLICATES:-30}"
TOTAL_READS="${TOTAL_READS:-2000000}"
READ_LENGTH="${READ_LENGTH:-150}"
READ_COUNTS="${READ_COUNTS:-1000,1500,2000,2500,3000,5000,7500,10000,25000,50000,100000,250000}"
PROFILE_FROM="${PROFILE_FROM:-}"
BUILD_INDEX="${BUILD_INDEX:-0}"

if [[ "${SMOKE:-0}" == "1" ]]; then
  echo "### SMOKE TEST: small read counts, 3 replicates ###"
  TOTAL_READS=200000
  REPLICATES=3
  READ_COUNTS="1000,5000,20000"
fi

# --------------------------------------------------------------------------- #
# 0. Dependencies
# --------------------------------------------------------------------------- #

echo "== checking dependencies =="
MISSING=()
for tool in samtools art_illumina bowtie2 bowtie2-inspect python3 curl unzip; do
  command -v "$tool" >/dev/null 2>&1 || MISSING+=("$tool")
done
if [[ "$BUILD_INDEX" == "1" ]]; then
  command -v bowtie2-build >/dev/null 2>&1 || MISSING+=("bowtie2-build")
  command -v datasets >/dev/null 2>&1 || MISSING+=("datasets")
fi
if [[ -n "$PROFILE_FROM" ]]; then
  command -v art_profiler_illumina >/dev/null 2>&1 || MISSING+=("art_profiler_illumina")
fi
if [[ ${#MISSING[@]} -gt 0 ]]; then
  echo "missing from PATH: ${MISSING[*]}" >&2
  echo "  conda install -c bioconda samtools art bowtie2" >&2
  echo "  conda install -c conda-forge ncbi-datasets-cli   # only with BUILD_INDEX=1" >&2
  exit 1
fi
for f in karyoratio.py calibrate.py chromfile_from_report.py; do
  [[ -f "$KARYODIR/$f" ]] || { echo "not found: $KARYODIR/$f" >&2; exit 1; }
done
python3 -c "import pysam, numpy, scipy" 2>/dev/null \
  || { echo "python needs pysam, numpy and scipy" >&2; exit 1; }
echo "all present"

mkdir -p "$WORKDIR"
cd "$WORKDIR"
mkdir -p calib

REF="reference.fa"
INDEX="calib/ref_index"
CHROMS="canfam_chroms.tsv"

# --------------------------------------------------------------------------- #
# 1. Aligner index and reference
# --------------------------------------------------------------------------- #

if [[ "$BUILD_INDEX" == "1" ]]; then

  if [[ ! -f "$REF.fai" ]]; then
    echo "== downloading $ASSEMBLY ($ACCESSION) from NCBI =="
    datasets download genome accession "$ACCESSION" \
        --include genome,seq-report --filename dog.zip
    unzip -o -q dog.zip -d ncbi && rm -f dog.zip
    FNA=$(find ncbi -name "*_genomic.fna" | head -1)
    [[ -n "$FNA" ]] || { echo "no genomic FASTA in the download" >&2; exit 1; }
    mv "$FNA" "$REF"
    samtools faidx "$REF"
  fi

  if [[ ! -f "${INDEX}.1.bt2" ]]; then
    echo "== building the bowtie2 index; expect one to three hours =="
    bowtie2-build --threads "$THREADS" "$REF" "$INDEX"
  fi

else

  # The index is considered usable only when the whole set is present; a
  # partial set produces a confusing complaint about an "older version".
  index_complete () {
    local base="$1" ext
    for ext in 1 2 3 4 rev.1 rev.2; do
      [[ -f "${base}.${ext}.bt2" || -f "${base}.${ext}.bt2l" ]] || return 1
    done
    return 0
  }

  if ! index_complete "$INDEX"; then
    if [[ ! -d bt2_index ]]; then
      echo "== downloading the prebuilt bowtie2 index for $ASSEMBLY =="
      echo "   $INDEX_URL"
      curl -L --fail --progress-bar -o bt2_index.zip "$INDEX_URL"
      unzip -o -q bt2_index.zip -d bt2_index && rm -f bt2_index.zip
    else
      echo "== linking the already downloaded index =="
    fi

    # Match only the forward .1 file: '*.1.bt2*' would also match
    # '*.rev.1.bt2', which would truncate the basename to '....rev'.
    FIRST=$(find bt2_index -name "*.1.bt2" -o -name "*.1.bt2l" \
            | grep -v '\.rev\.' | head -1)
    [[ -n "$FIRST" ]] || { echo "no bowtie2 index inside the archive" >&2; exit 1; }
    BASE="${FIRST%.1.bt2}"; BASE="${BASE%.1.bt2l}"
    echo "   index basename: $BASE"

    rm -f calib/ref_index.*
    for part in "$BASE".*.bt2 "$BASE".*.bt2l; do
      [[ -e "$part" ]] || continue
      ln -sf "$(cd "$(dirname "$part")" && pwd)/$(basename "$part")" \
             "calib/ref_index.${part#"$BASE".}"
    done

    echo "   linked: $(ls calib/ref_index.* 2>/dev/null | wc -l) files"
    if ! index_complete "$INDEX"; then
      echo "the index set is incomplete after linking:" >&2
      ls -l calib/ref_index.* >&2 || true
      echo "expected .1 .2 .3 .4 .rev.1 .rev.2; inspect bt2_index/" >&2
      exit 1
    fi
  fi

  if [[ ! -f "$REF.fai" ]]; then
    echo "== reconstructing the reference FASTA from the index =="
    echo "   (this guarantees the sequence names match what the aligner reports)"
    bowtie2-inspect "$INDEX" > "$REF"
    samtools faidx "$REF"
  fi

fi

echo "   sequences in the index: $(bowtie2-inspect -n "$INDEX" | wc -l)"

# --------------------------------------------------------------------------- #
# 2. Chromosome definition file
# --------------------------------------------------------------------------- #
# Roles come from the NCBI sequence report, matched to whatever naming the
# index uses: accessions, UCSC-style names or bare chromosome numbers.

if [[ -f "$CHROMS" ]]; then
  echo "== chromosome file already present, skipping =="
else
  echo "== building $CHROMS =="

  REPORT=$(find ncbi -name "sequence_report.jsonl" 2>/dev/null | head -1 || true)
  if [[ -z "$REPORT" ]] && command -v datasets >/dev/null 2>&1; then
    echo "   fetching the sequence report from NCBI"
    datasets download genome accession "$ACCESSION" \
        --include seq-report --filename seqrep.zip
    unzip -o -q seqrep.zip -d ncbi && rm -f seqrep.zip
    REPORT=$(find ncbi -name "sequence_report.jsonl" | head -1)
  fi
  [[ -n "$REPORT" ]] || echo "   no sequence report; falling back to name patterns"

  bowtie2-inspect -n "$INDEX" > index_names.txt

  set +e
  python3 "$KARYODIR/chromfile_from_report.py" \
      --names index_names.txt \
      ${REPORT:+--report "$REPORT"} \
      --max-autosome 38 \
      --expect-autosomes 38 \
      --comment "Canis lupus familiaris, $ASSEMBLY ($ACCESSION)" \
      -o "$CHROMS"
  RC=$?
  set -e
  # Exit 2 is the "unexpected autosome count" warning, which is survivable;
  # anything else, including a missing script, is not.
  if [[ $RC -ne 0 && $RC -ne 2 ]] || [[ ! -s "$CHROMS" ]]; then
    rm -f "$CHROMS"
    echo "could not build $CHROMS (exit $RC); inspect index_names.txt" >&2
    exit 1
  fi

  N_AUT=$(grep -c $'\tautosome$' "$CHROMS")
  echo "   sanity check: the dog karyotype is 2n = 78, so 38 autosomes expected"
  echo "   autosomes found: $N_AUT"
  if [[ "$N_AUT" != "38" ]]; then
    echo "   [warn] that is not 38; inspect $CHROMS before continuing"
  fi
fi

# --------------------------------------------------------------------------- #
# 3. Calibration
# --------------------------------------------------------------------------- #

echo
echo "== running the calibration =="
echo "   replicates : $REPLICATES"
echo "   read counts: $READ_COUNTS"
echo "   simulated  : $TOTAL_READS pairs per sex"
echo

PROFILE_ARG=()
if [[ -n "$PROFILE_FROM" ]]; then
  PROFILE_ARG=(--profile-from "$PROFILE_FROM")
  echo "   error profile derived from $PROFILE_FROM"
else
  echo "   no --profile-from given, so ART's built-in HiSeq 2500 model is used."
  echo "   Set PROFILE_FROM=/path/to/real_R1.fq.gz to match your own data."
fi
echo

python3 "$KARYODIR/calibrate.py" \
    --reference "$REF" \
    --chrom-file "$CHROMS" \
    --index "$INDEX" \
    --outdir calib \
    --karyoratio "$KARYODIR/karyoratio.py" \
    --read-counts "$READ_COUNTS" \
    --replicates "$REPLICATES" \
    --total-reads "$TOTAL_READS" \
    --read-length "$READ_LENGTH" \
    --threads "$THREADS" \
    --skip-existing \
    "${PROFILE_ARG[@]}"

# --------------------------------------------------------------------------- #
# 4. What came out
# --------------------------------------------------------------------------- #

echo
echo "== summary =="
column -t -s$'\t' calib/summary.tsv 2>/dev/null || cat calib/summary.tsv

cat <<EOF

Files written under $(pwd)/calib:

  replicates.tsv   every individual result; plot ratio against read_count,
                   coloured by sex, with the 0.60 and 0.80 thresholds marked
  summary.tsv      correct / ambiguous / wrong counts per read count and sex

Reading the result:

  * any value in the 'wrong' column means the 0.60/0.80 thresholds do not
    transfer to the dog and need adjusting before use
  * for the depth to pass to karyoratio as --min-reads, read the rate table
    printed above rather than the first depth with no unclassified replicate;
    that depth rises with the number of replicates and is not a property of
    the method

Next, validate on real samples of known sex:

  ./run_dog_controls.sh $WORKDIR
EOF
