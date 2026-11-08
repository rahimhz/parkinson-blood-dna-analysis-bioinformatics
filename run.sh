#!/usr/bin/env bash
# Download, validate, test, and optionally analyze processed GSE165081 M-values.
set -Eeuo pipefail

ACCESSION="GSE165081"
BASE_URL="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE165nnn/${ACCESSION}/suppl/"
MATRIX_FILE="GSE165081_M.csv.gz"
SHEET_FILE="GSE165081_Sample_Sheet.PD-CN.csv.gz"
ROOT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=scripts/run_options.sh
source "${ROOT_DIR}/scripts/run_options.sh"

FORCE_DOWNLOAD=0
REQUESTED_DESIGN=auto
DESIGN=
SAMPLE_COLUMN=
GROUP_COLUMN=
PAIR_COLUMN=
PAIR_OPTION_SET=0
CURRENT_DOWNLOAD_TEMP=
TEMP_DIR=
RUN_DIR=
LOG_FILE=

usage() {
  cat <<'USAGE'
Usage: ./run.sh [options]
  -h            Show help and exit.
  -f            Force re-download the two supplementary GEO files.
  -d DESIGN     Design: auto (default), independent, or adjusted.
  -s COLUMN     Actual sample-ID header if detection is ambiguous.
  -g COLUMN     Actual group header if detection is ambiguous.
  -p COLUMN     Verified one-to-one pair-ID header; selects paired analysis.

Auto mode validates inputs, writes a validation summary and archive, then
exits successfully without statistical testing. -p requires verified pairs.
-d independent requests an exploratory unadjusted Welch test.
-d adjusted requests OLS with PD status, age, and sex covariates.
USAGE
}

while getopts ':hfd:s:g:p:' option; do
  case "$option" in
    h) usage; exit 0 ;;
    f) FORCE_DOWNLOAD=1 ;;
    d) REQUESTED_DESIGN="$OPTARG" ;;
    s) SAMPLE_COLUMN="$OPTARG" ;;
    g) GROUP_COLUMN="$OPTARG" ;;
    p) PAIR_COLUMN="$OPTARG"; PAIR_OPTION_SET=1 ;;
    :) echo "ERROR: option -$OPTARG requires a value." >&2; usage >&2; exit 2 ;;
    \?) echo "ERROR: unknown option -$OPTARG." >&2; usage >&2; exit 2 ;;
  esac
done
shift "$((OPTIND - 1))"
if (( $# > 0 )); then
  echo "ERROR: unexpected positional arguments: $*" >&2
  usage >&2
  exit 2
fi
if (( PAIR_OPTION_SET == 1 )) && [[ -z "$PAIR_COLUMN" ]]; then
  echo 'ERROR: -p requires a non-empty pair-column name.' >&2
  exit 2
fi
DESIGN="$(resolve_design "$REQUESTED_DESIGN" "$PAIR_COLUMN")" || exit $?

cd "$ROOT_DIR"
mkdir -p data/raw data/processed results/runs logs
RUN_ID="$(date -u '+%Y%m%dT%H%M%SZ')-$$"
RUN_DIR="$(mktemp -d "results/runs/${RUN_ID}.XXXXXX")"
RUN_ID="$(basename "$RUN_DIR")"
LOG_FILE="logs/run_${RUN_ID}.log"
CHECKSUM_FILE="logs/SHA256SUMS_${RUN_ID}.txt"
ARCHIVE_FILE="${RUN_DIR}/results_${RUN_ID}.tar.gz"
exec > >(tee -a "$LOG_FILE") 2>&1

cleanup() {
  local status=$?
  if [[ -n "$CURRENT_DOWNLOAD_TEMP" && -e "$CURRENT_DOWNLOAD_TEMP" ]]; then
    rm -f -- "$CURRENT_DOWNLOAD_TEMP"
  fi
  if [[ -n "$TEMP_DIR" && -d "$TEMP_DIR" ]]; then
    rm -rf -- "$TEMP_DIR"
  fi
  return "$status"
}
on_error() {
  local status="$1"
  trap - ERR
  log "ERROR: command failed (exit ${status}); see ${LOG_FILE}."
  exit "$status"
}
on_signal() {
  local status="$1"
  log "ERROR: interrupted (exit ${status})."
  exit "$status"
}
trap cleanup EXIT
trap 'on_error "$?"' ERR
trap 'on_signal 130' INT
trap 'on_signal 143' TERM

log() { printf '[%s] %s\n' "$(date -u '+%Y-%m-%dT%H:%M:%SZ')" "$*"; }
fail() { log "ERROR: $*"; exit 1; }

check_dependencies() {
  local command_name
  for command_name in wget gzip sha256sum awk tr head grep tar tee sort cp mktemp find mkdir chmod mv rm basename date dirname; do
    command -v "$command_name" >/dev/null 2>&1 || fail "Required command not found: ${command_name}."
  done
  command -v python >/dev/null 2>&1 || fail 'Python not found; activate environment.yml first.'
  python -c 'import numpy, pandas, scipy, statsmodels, matplotlib' ||
    fail 'Python packages missing; activate environment.yml first.'
}

download_file() {
  local filename="$1"
  local destination="data/raw/${filename}"
  local temporary="${destination}.part"
  if [[ -s "$destination" && "$FORCE_DOWNLOAD" -eq 0 ]]; then
    log "Reusing existing download: ${destination}"
  else
    rm -f -- "$temporary"
    CURRENT_DOWNLOAD_TEMP="$temporary"
    log "Downloading ${filename} from official NCBI GEO FTP."
    wget --https-only --tries=3 --timeout=30 --output-document="$temporary" "${BASE_URL}${filename}"
    [[ -s "$temporary" ]] || fail "Downloaded file is empty: ${filename}"
    gzip -t "$temporary" || fail "gzip integrity check failed: ${filename}"
    mv -- "$temporary" "$destination"
    CURRENT_DOWNLOAD_TEMP=
  fi
  gzip -t "$destination" || fail "gzip integrity check failed: ${destination}"
  chmod 0444 "$destination"
}

verify_downloads() {
  log 'Calculating SHA-256 checksums.'
  (cd data/raw && sha256sum "$MATRIX_FILE" "$SHEET_FILE") > "$ROOT_DIR/$CHECKSUM_FILE"
  log 'Verifying local checksum manifest.'
  (cd data/raw && sha256sum --check "$ROOT_DIR/$CHECKSUM_FILE")
}

clean_sample_sheet() {
  local source="data/raw/${SHEET_FILE}"
  local destination="data/processed/sample_sheet_clean.csv"
  local temporary="${destination}.tmp"
  log 'Decompressing the sample sheet without rewriting CSV records.'
  # Preserve quoted fields and embedded newlines; pandas handles CRLF and blank rows.
  gzip -dc "$source" > "$temporary"
  [[ -s "$temporary" ]] || fail 'Cleaned sample sheet is empty.'
  mv -- "$temporary" "$destination"
  grep -q '[^[:space:]]' "$destination" || fail 'Sample sheet has no non-empty content.'
  log 'Sample-sheet preview (header and up to four rows):'
  head -n 5 "$destination"
  log 'Non-empty physical line count, including header:'
  awk 'NF { count++ } END { print count+0 }' "$destination"
}

run_tests() {
  log 'Running offline tests on synthetic data and CLI options.'
  python -m unittest discover -s tests -v
}

run_analysis() {
  local arguments=(--matrix "data/raw/${MATRIX_FILE}" --samples data/processed/sample_sheet_clean.csv --output-dir "$RUN_DIR" --design "$DESIGN")
  [[ -z "$SAMPLE_COLUMN" ]] || arguments+=(--sample-column "$SAMPLE_COLUMN")
  [[ -z "$GROUP_COLUMN" ]] || arguments+=(--group-column "$GROUP_COLUMN")
  [[ -z "$PAIR_COLUMN" ]] || arguments+=(--pair-column "$PAIR_COLUMN")
  log "Running analysis with design=${DESIGN}."
  python scripts/analyze.py "${arguments[@]}"
}

archive_run() {
  local members expected
  TEMP_DIR="$(mktemp -d)"
  mkdir -p "${TEMP_DIR}/run"
  while IFS= read -r -d '' file; do
    [[ "$file" == "$ARCHIVE_FILE" ]] && continue
    relative="${file#${RUN_DIR}/}"
    mkdir -p "${TEMP_DIR}/run/$(dirname "$relative")"
    cp -- "$file" "${TEMP_DIR}/run/${relative}"
  done < <(find "$RUN_DIR" -type f -print0)
  [[ -n "$(find "${TEMP_DIR}/run" -type f -print -quit)" ]] || fail 'No run output files to archive.'
  # Write archive beside, not inside, source run directory to avoid self-inclusion.
  ARCHIVE_FILE="results/runs/${RUN_ID}.tar.gz"
  tar -czf "$ARCHIVE_FILE" -C "$TEMP_DIR" run
  members="$(tar -tzf "$ARCHIVE_FILE" | awk '!/\/$/ { sub(/^run\//, ""); print }' | sort)"
  expected="$(find "${TEMP_DIR}/run" -type f -printf '%P\n' | sort)"
  [[ "$members" == "$expected" ]] || fail 'Archive members do not match current run regular files.'
  log 'Verified archive members:'
  printf '%s\n' "$members"
  rm -rf -- "$TEMP_DIR"
  TEMP_DIR=
}

main() {
  check_dependencies
  log "Starting ${ACCESSION} processed M-value workflow; design=${DESIGN}."
  download_file "$MATRIX_FILE"
  download_file "$SHEET_FILE"
  verify_downloads
  clean_sample_sheet
  run_tests
  run_analysis
  archive_run
  log 'Pipeline completed successfully.'
  log "Log: ${LOG_FILE}"
  log "Checksum manifest: ${CHECKSUM_FILE}"
  log "Archive: ${ARCHIVE_FILE}"
}
main "$@"
