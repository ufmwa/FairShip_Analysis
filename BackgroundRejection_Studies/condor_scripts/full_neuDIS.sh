#!/usr/bin/env bash
# Läuft neuDIS-Kanäle über viele job_* in einem Eingabeordner.
# Optional: nur Jobs aus einer Jobliste-Datei abarbeiten.
# Optional: nur einen Kanal rechnen (partialreco|fullreco|leptonrho).
#
# Usage:
#   bash full_neuDIS.sh <INPDIR> <OUTBASE> <SCRIPTDIR> [JOBLIST] [CHANNEL] [CASE]
#
# JOBLIST:
#   Datei mit einem Job pro Zeile, z.B.:
#     job_1770
#     9950
#   Leere Zeilen und Zeilen beginnend mit # werden ignoriert.
#
# CHANNEL:
#   partialreco | fullreco | leptonrho | all
#   Default: all
#
# CASE:
#   all | vesselCase | heliumCase | caveCase | ubtCase
#   Default: all

set -u
set -o pipefail

if [ $# -lt 3 ] || [ $# -gt 6 ]; then
  echo "Usage: $0 <INPDIR> <OUTBASE> <SCRIPTDIR> [JOBLIST] [CHANNEL] [CASE]" >&2
  exit 1
fi

INPDIR="$1"
OUTBASE="$2"
SCRIPTDIR="$3"
JOBLIST="${4:-}"
CHANNEL_SEL="${5:-all}"
CASE_SEL="${6:-all}"

# --- Kanal-Auswahl prüfen ---
case "$CHANNEL_SEL" in
  all|partialreco|fullreco|leptonrho) ;;
  *)
    echo "Invalid CHANNEL '$CHANNEL_SEL'. Use: all|partialreco|fullreco|leptonrho" >&2
    exit 1
    ;;
esac

# --- Case-Auswahl prüfen ---
case "$CASE_SEL" in
  all|vesselCase|heliumCase|caveCase|ubtCase) ;;
  *)
    echo "Invalid CASE '$CASE_SEL'. Use: all|vesselCase|heliumCase|caveCase|ubtCase" >&2
    exit 1
    ;;
esac

# --- Hilfsfunktion: 1 Job x 1 Kanal ---
run_one() {
  local JOB="$1"
  local CHANNEL="$2"
  local FLAG

  case "$CHANNEL" in
    partialreco) FLAG="--partialreco" ;;
    fullreco)    FLAG="--fullreco" ;;
    leptonrho)   FLAG="--leptonrho" ;;
    *) echo "Unknown channel: $CHANNEL" >&2; return 2 ;;
  esac

  local CASEARGS=()
  if [ "$CASE_SEL" != "all" ]; then
    CASEARGS=(--case "$CASE_SEL")
  fi

  echo ">>> [$JOB][$CHANNEL] start $(date)"
  rm -f selectionparameters_*.root selection_summary_*.csv

  if ! python "$SCRIPTDIR/BackgroundRejection_Studies/run_neuDIS.py" \
        -p "$INPDIR" -i "$JOB" "$FLAG" "${CASEARGS[@]}" ; then
    echo "!!! [$JOB][$CHANNEL] FAILED" >&2
    return 3
  fi

  local OUTDIR="$OUTBASE/neuDIS/$CHANNEL/$JOB"
  mkdir -p "$OUTDIR"
  cp -f selectionparameters_*.root selection_summary_*.csv "$OUTDIR"/ 2>/dev/null || true
  rm -f selectionparameters_*.root selection_summary_*.csv

  echo "<<< [$JOB][$CHANNEL] done  $(date)"
}

# --- Logging ---
LOGDIR="$OUTBASE/logs"
mkdir -p "$LOGDIR"

# --- Jobliste bauen ---
declare -a JOBS=()

if [ -n "$JOBLIST" ]; then
  if [ ! -f "$JOBLIST" ]; then
    echo "JOBLIST file not found: $JOBLIST" >&2
    exit 1
  fi

  while IFS= read -r JOB || [[ -n "$JOB" ]]; do
    JOB="${JOB//$'\r'/}"          # CR entfernen
    JOB="${JOB%%#*}"              # Kommentare abschneiden
    JOB="$(echo "$JOB" | xargs)"  # whitespace trim
    [[ -z "$JOB" ]] && continue

    # optional: akzeptiere auch "9950" und mache "job_9950" draus
    [[ "$JOB" =~ ^[0-9]+$ ]] && JOB="job_$JOB"

    if [[ -d "$INPDIR/$JOB" ]]; then
      JOBS+=("$JOB")
    else
      echo "Skipping (not found in INPDIR): $JOB" >&2
    fi
  done < "$JOBLIST"

  if [ ${#JOBS[@]} -eq 0 ]; then
    echo "No valid jobs found in JOBLIST: $JOBLIST" >&2
    exit 1
  fi
  echo "Using JOBLIST '$JOBLIST' with ${#JOBS[@]} jobs."
else
  shopt -s nullglob
  JOBPATHS=("$INPDIR"/job_*)
  if [ ${#JOBPATHS[@]} -eq 0 ]; then
    echo "Keine job_* Ordner in $INPDIR gefunden." >&2
    exit 1
  fi
  for JP in "${JOBPATHS[@]}"; do
    [ -d "$JP" ] || continue
    JOBS+=("$(basename "$JP")")
  done
  echo "Discovered ${#JOBS[@]} jobs in $INPDIR."
fi

# --- Welche Kanäle laufen? ---
declare -a CHANNELS=()
if [ "$CHANNEL_SEL" = "all" ]; then
  CHANNELS=(partialreco fullreco leptonrho)
else
  CHANNELS=("$CHANNEL_SEL")
fi

echo "Channels to run: ${CHANNELS[*]}"

# --- Jobs abarbeiten ---
for JOB in "${JOBS[@]}"; do
  JP="$INPDIR/$JOB"
  if [ ! -d "$JP" ]; then
    echo "WARN: job dir missing, skip: $JP" >&2
    continue
  fi

  {
    echo "===== $JOB ====="
    for CH in "${CHANNELS[@]}"; do
      run_one "$JOB" "$CH"
    done
    echo "===== $JOB DONE ====="
  } > "$LOGDIR/${JOB}.log" 2>&1
done

echo "Alle Jobs fertig. Reports/Logs unter: $LOGDIR"
echo "Ergebnisse unter: $OUTBASE/neuDIS/{partialreco,fullreco,leptonrho}/job_*/"
