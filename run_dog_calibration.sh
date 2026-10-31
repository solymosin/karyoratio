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
# The dog karyotype is 2n = 78: 38 autosomes plus X. The assembly carries no
# Y chromosome, which is why the X-to-autosome ratio is the right estimator
# here -- a Y-based one is not available.
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
READ_LENGTH="${READ_LENGTH:-151}"
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
for f in karyoratio.py calibrate.py; do
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

  if [[ ! -f "${INDEX}.1.bt2" ]]; then
    echo "== downloading the prebuilt bowtie2 index for $ASSEMBLY =="
    echo "   $INDEX_URL"
    curl -L --fail --progress-bar -o bt2_index.zip "$INDEX_URL"
    unzip -o -q bt2_index.zip -d bt2_index && rm -f bt2_index.zip

    FIRST=$(find bt2_index -name "*.1.bt2*" | head -1)
    [[ -n "$FIRST" ]] || { echo "no bowtie2 index inside the archive" >&2; exit 1; }
    BASE="${FIRST%.1.bt2*}"
    echo "   index basename: $BASE"
    for part in "$BASE".*; do
      ln -sf "$(cd "$(dirname "$part")" && pwd)/$(basename "$part")" \
             "calib/ref_index.${part#"$BASE".}"
    done
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

  python3 - "${REPORT:-}" index_names.txt "$CHROMS" <<'PY'
import json, sys
report, names_file, out = sys.argv[1], sys.argv[2], sys.argv[3]
names = [l.strip() for l in open(names_file) if l.strip()]
nameset = set(names)

autosomes, x, y, mt = [], None, None, None

def pick(rec):
    """Which of this record's identifiers does the index actually use?"""
    chrname = str(rec.get("chrName", ""))
    for cand in (rec.get("refseqAccession"), rec.get("genbankAccession"),
                 rec.get("ucscStyleName"), f"chr{chrname}", chrname):
        if cand and cand in nameset:
            return cand
    return None

if report:
    for line in open(report):
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
else:
    # Without a sequence report, fall back to the names themselves, keeping
    # only the 38 autosomes the dog karyotype specifies.
    for n in names:
        bare = n[3:] if n.lower().startswith("chr") else n
        if bare.upper() == "X":
            x = n
        elif bare.upper() == "Y":
            y = n
        elif bare.upper() in ("M", "MT"):
            mt = n
        elif bare.isdigit() and 1 <= int(bare) <= 38:
            autosomes.append(n)

if x is None:
    sys.exit("could not identify the X chromosome; inspect index_names.txt "
             "and write the chromosome file by hand")
if not autosomes:
    sys.exit("could not identify any autosome; inspect index_names.txt")

with open(out, "w") as fh:
    fh.write("# Canis lupus familiaris, ROS_Cfam_1.0 (GCF_014441545.1)\n")
    fh.write(f"# {len(autosomes)} autosomes + X"
             + (", Y" if y else ", no Y in this assembly")
             + (", MT\n" if mt else "\n"))
    for acc in autosomes:
        fh.write(f"{acc}\tautosome\n")
    fh.write(f"{x}\tX\n")
    if y:
        fh.write(f"{y}\tY\n")
    if mt:
        fh.write(f"{mt}\tMT\n")

print(f"  {len(autosomes)} autosomes, X = {x}"
      + (f", Y = {y}" if y else ", no Y")
      + (f", MT = {mt}" if mt else ""))
PY

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
  * the lowest read count whose 'ambiguous' column is zero for both sexes is
    the figure to pass to karyoratio as --min-reads

Next, validate on real samples of known sex:

  ./run_dog_controls.sh $WORKDIR
EOF
