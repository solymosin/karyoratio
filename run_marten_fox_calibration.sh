#!/usr/bin/env bash
#
# run_marten_fox_calibration.sh -- calibrate karyoratio for the stone marten
# and the red fox from genomes already present on the server.
#
#   ../genomes/GCA_964304585.1_mMarFoi2.1_genomic.fna.gz
#   ../genomes/GCF_048418805.1_VulVul3_genomic.fna.gz
#
# Both are plain gzip, which samtools faidx and pysam cannot index, so each is
# decompressed once into the working directory.
#
# The index is built from the host genome alone. Competitive mapping against a
# host-plus-human reference belongs to the analysis of real samples, where it
# absorbs contaminating human reads; it has no role here, because the
# simulated reads are all host-derived and the human contigs would simply sit
# empty. What the calibration asks is whether the thresholds hold for this
# karyotype, and that question is about the host genome.
#
# Building an index for a 2.4 Gb genome takes one to three hours per species.
# MARTEN_INDEX and FOX_INDEX can point at an existing host-only index to skip
# it.
#
# Usage:
#   ./run_marten_fox_calibration.sh [workdir]
#
# SMOKE=1        quick end-to-end check first
# SPECIES=marten only one species
#
set -euo pipefail

WORKDIR="${1:-mf_calibration}"
KARYODIR="${KARYODIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
GENOMEDIR="${GENOMEDIR:-../genomes}"

MARTEN_GZ="${MARTEN_GZ:-$GENOMEDIR/GCA_964304585.1_mMarFoi2.1_genomic.fna.gz}"
FOX_GZ="${FOX_GZ:-$GENOMEDIR/GCF_048418805.1_VulVul3_genomic.fna.gz}"

MARTEN_CHROMS="${MARTEN_CHROMS:-$KARYODIR/mmarfoi2.1_chroms.tsv}"
FOX_CHROMS="${FOX_CHROMS:-$KARYODIR/vulvul3_chroms.tsv}"

MARTEN_INDEX="${MARTEN_INDEX:-}"
FOX_INDEX="${FOX_INDEX:-}"

# A real FASTQ per species gives ART an error profile matched to the data and
# lifts the 150 bp ceiling of the built-in profiles.
MARTEN_PROFILE="${MARTEN_PROFILE:-}"
FOX_PROFILE="${FOX_PROFILE:-}"

THREADS="${THREADS:-8}"
REPLICATES="${REPLICATES:-30}"
TOTAL_READS="${TOTAL_READS:-2000000}"
READ_LENGTH="${READ_LENGTH:-150}"
READ_COUNTS="${READ_COUNTS:-250,500,750,1000,1500,2000,2500,3000,5000,7500,10000,25000,50000,100000,250000}"
MIN_MAPQ="${MIN_MAPQ:-30}"
SPECIES="${SPECIES:-marten fox}"

if [[ "${SMOKE:-0}" == "1" ]]; then
  echo "### SMOKE TEST ###"
  TOTAL_READS=200000
  REPLICATES=3
  READ_COUNTS="1000,5000,20000"
fi

echo "== checking dependencies =="
MISSING=()
for tool in samtools art_illumina bowtie2 python3; do
  command -v "$tool" >/dev/null 2>&1 || MISSING+=("$tool")
done
for f in karyoratio.py calibrate.py; do
  [[ -f "$KARYODIR/$f" ]] || { echo "not found: $KARYODIR/$f" >&2; exit 1; }
done
if [[ ${#MISSING[@]} -gt 0 ]]; then
  echo "missing from PATH: ${MISSING[*]}" >&2
  exit 1
fi
python3 -c "import pysam, numpy, scipy" 2>/dev/null \
  || { echo "python needs pysam, numpy and scipy" >&2; exit 1; }
echo "all present"

# Resolve every input path before changing directory, so relative paths are
# read as the user meant them -- relative to where the script was invoked, not
# to the working directory it creates.
abspath () {
  case "$1" in
    "" ) echo "" ;;
    /* ) echo "$1" ;;
    *  ) echo "$PWD/$1" ;;
  esac
}
MARTEN_GZ=$(abspath "$MARTEN_GZ")
FOX_GZ=$(abspath "$FOX_GZ")
MARTEN_CHROMS=$(abspath "$MARTEN_CHROMS")
FOX_CHROMS=$(abspath "$FOX_CHROMS")
MARTEN_INDEX=$(abspath "$MARTEN_INDEX")
FOX_INDEX=$(abspath "$FOX_INDEX")
MARTEN_PROFILE=$(abspath "$MARTEN_PROFILE")
FOX_PROFILE=$(abspath "$FOX_PROFILE")

mkdir -p "$WORKDIR"
cd "$WORKDIR"

# --------------------------------------------------------------------------- #

run_species () {
  local name="$1" gz="$2" chroms="$3" index="$4" profile="$5"
  local ref="${name}.fa"

  echo
  echo "########## $name ##########"

  [[ -f "$gz" ]] || { echo "genome not found: $gz" >&2; return 1; }
  [[ -f "$chroms" ]] || { echo "chromosome file not found: $chroms" >&2; return 1; }

  # 1. Decompress and index. Plain gzip is not seekable, so a copy is needed.
  if [[ ! -f "$ref.fai" ]]; then
    echo "== decompressing $(basename "$gz") =="
    gzip -dc "$gz" > "$ref"
    echo "== samtools faidx =="
    samtools faidx "$ref"
  else
    echo "== reference already prepared =="
  fi
  echo "   $(grep -c . "$ref.fai") sequences, $(awk '{s+=$2} END {printf "%.2f Gb", s/1e9}' "$ref.fai")"

  # 2. The chromosome file has to name sequences that exist.
  local missing
  missing=$(awk '$1 !~ /^#/ && NF>=2 {print $1}' "$chroms" | sort -u \
            | comm -23 - <(cut -f1 "$ref.fai" | sort -u) | head -5)
  if [[ -n "$missing" ]]; then
    echo "these sequences from $chroms are not in $ref:" >&2
    echo "$missing" >&2
    return 1
  fi
  echo "   chromosome file matches the reference"
  echo "   $(awk '$2=="autosome"' "$chroms" | wc -l) autosomes, X = $(awk '$2=="X" {print $1}' "$chroms")"

  # 3. Index: use the one given, or build a host-only one.
  local idx
  if [[ -n "$index" ]]; then
    if [[ ! -f "${index}.1.bt2" && ! -f "${index}.1.bt2l" ]]; then
      echo "no bowtie2 index at $index" >&2
      return 1
    fi
    idx="$index"
    echo "== using the index at $idx =="
    local n_idx n_ref
    n_idx=$(bowtie2-inspect -n "$idx" | wc -l)
    n_ref=$(grep -c . "$ref.fai")
    echo "   index has $n_idx sequences, the reference $n_ref"
    if [[ "$n_idx" -gt "$n_ref" ]]; then
      echo "   [warn] the index holds more sequences than the host reference."
      echo "          A composite index belongs to the analysis of real"
      echo "          samples, not here; the calibration should run against"
      echo "          the host genome alone."
    fi
  else
    idx="${name}_index"
    if [[ ! -f "${idx}.1.bt2" && ! -f "${idx}.1.bt2l" ]]; then
      echo "== building the bowtie2 index; expect one to three hours =="
      command -v bowtie2-build >/dev/null 2>&1 \
        || { echo "bowtie2-build not on PATH" >&2; return 1; }
      bowtie2-build --threads "$THREADS" "$ref" "$idx"
    else
      echo "== index already built =="
    fi
  fi

  # 4. Calibrate.
  local profile_arg=()
  if [[ -n "$profile" ]]; then
    [[ -f "$profile" ]] || { echo "profile FASTQ not found: $profile" >&2; return 1; }
    profile_arg=(--profile-from "$profile")
    echo "   error profile from $(basename "$profile")"
  else
    echo "   ART's built-in HiSeq 2500 model (max 150 bp); set ${name^^}_PROFILE"
    echo "   to a real FASTQ to match your own data and allow 151 bp"
  fi

  echo
  echo "== calibrating $name =="
  python3 "$KARYODIR/calibrate.py" \
      --reference "$ref" \
      --chrom-file "$chroms" \
      --index "$idx" \
      --outdir "calib_${name}" \
      --karyoratio "$KARYODIR/karyoratio.py" \
      --read-counts "$READ_COUNTS" \
      --replicates "$REPLICATES" \
      --total-reads "$TOTAL_READS" \
      --read-length "$READ_LENGTH" \
      --min-mapq "$MIN_MAPQ" \
      --threads "$THREADS" \
      --skip-existing \
      "${profile_arg[@]}"

  echo
  echo "== $name summary =="
  column -t -s$'\t' "calib_${name}/summary.tsv" 2>/dev/null \
    || cat "calib_${name}/summary.tsv"
}

for sp in $SPECIES; do
  case "$sp" in
    marten) run_species marten "$MARTEN_GZ" "$MARTEN_CHROMS" \
                        "$MARTEN_INDEX" "$MARTEN_PROFILE" ;;
    fox)    run_species fox    "$FOX_GZ"    "$FOX_CHROMS" \
                        "$FOX_INDEX"    "$FOX_PROFILE" ;;
    *) echo "unknown species: $sp" >&2; exit 1 ;;
  esac
done

# --------------------------------------------------------------------------- #

echo
echo "########## across species ##########"
cat <<'EOF'
The depth at which the unclassified rate falls below a given level. A rate is
reported rather than the first depth with no unclassified replicate, because
that depth is not a property of the method: run more replicates and rare
events surface, pushing it upward without anything having changed.

EOF
printf '%-10s %-8s %10s %10s %10s\n' species sex "<10%" "<5%" "<1%"
for sp in $SPECIES; do
  f="calib_${sp}/summary.tsv"
  [[ -f "$f" ]] || continue
  for sex in female male; do
    awk -F'\t' -v sp="$sp" -v sx="$sex" '
      $1==sx {n++; cnt[n]=$2; rate[n]=$7}
      END {
        split("10 5 1", want, " ")
        out = sprintf("%-10s %-8s", sp, sx)
        for (w=1; w<=3; w++) {
          ans="-"
          # the lowest depth from which the rate stays under the target
          for (i=n; i>=1; i--) { if (rate[i]+0 >= want[w]+0) break; ans=cnt[i] }
          out = out sprintf(" %10s", ans)
        }
        print out
      }' "$f"
  done
done

cat <<'EOF'

Incorrect assignments, as distinct from unclassified ones, are listed above
for each species; they appeared only at the lowest depths tested.

How far to trust these numbers:

  The precision of the ratio depends on how many autosomes stand in the
  denominator, since it is their mean that is reported. A species with fewer,
  larger chromosomes gives a noisier estimate at the same depth, which is one
  reason a threshold calibrated on one species does not carry to another.

  Every replicate here is drawn from a single simulated genome, so the spread
  is the spread of read sampling alone. Variation between individuals is not
  in it, and on thirty real dogs it was not small -- the female ratio ranged
  from 0.86 to 1.01. A minimum read count from simulation is therefore a
  floor, not a sufficient condition, and is not comparable with one obtained
  from real samples of known sex.
EOF
