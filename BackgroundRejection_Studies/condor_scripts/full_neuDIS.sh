#!/usr/bin/env bash
# Läuft alle neuDIS-Kanäle über alle job_* in einem Eingabeordner.
# Usage:
#   bash neuDIS_alljobs_local.sh <INPDIR> <OUTBASE> <SCRIPTDIR>
# Example:
#   bash neuDIS_alljobs_local.sh \
# /work/nouachem/FairShip_Analysis/BackgroundRejection_Studies/condor_scripts/full_neuDIS.sh
    # /storage/9/rquishpe/ship/NeutrinoDIS_2024helium_noCavern/test
    # /ceph/nouachem/Analysis_Sebastian/test
    # /work/nouachem/FairShip_Analysis

set -u
set -o pipefail

if [ $# -ne 3 ]; then
  echo "Usage: $0 <INPDIR> <OUTBASE> <SCRIPTDIR>" >&2
  exit 1
fi

INPDIR="$1"     # enthält viele job_* Ordner -> geht an run_neuDIS.py -p
OUTBASE="$2"    # lokaler Output-Basisordner für Ergebnisse (statt EOSDIR)
SCRIPTDIR="$3"  # Repo-Wurzel, die 'BackgroundRejection_Studies/' enthält

# --- Hilfsfunktion: 1 Job x 1 Kanal ---
run_one() {
  local JOB="$1"         # z.B. job_000123 -> geht an -i
  local CHANNEL="$2"     # partialreco | fullreco | leptonrho
  local FLAG
  case "$CHANNEL" in
    partialreco) FLAG="--partialreco" ;;
    fullreco)    FLAG="--fullreco" ;;
    leptonrho)   FLAG="--leptonrho" ;;
    *) echo "Unknown channel: $CHANNEL" >&2; return 2 ;;
  esac

  echo ">>> [$JOB][$CHANNEL] start $(date)"
  # Aufräumen, falls vom letzten Lauf was rumliegt:
  rm -f selectionparameters_*.root selection_summary_*.csv

  # Python-Analyse
  if ! python "$SCRIPTDIR/BackgroundRejection_Studies/run_neuDIS.py" \
        -p "$INPDIR" -i "$JOB" "$FLAG" ; then
    echo "!!! [$JOB][$CHANNEL] FAILED" >&2
    return 3
  fi

  # Ergebnisse einsortieren (Merge-kompatible Struktur!)
  local OUTDIR="$OUTBASE/neuDIS/$CHANNEL/$JOB"
  mkdir -p "$OUTDIR"
  # Pro Job werden mehrere Kategorien geschrieben: all / heliumCase / vesselCase
  cp selectionparameters_*.root selection_summary_*.csv "$OUTDIR"/ || true
  rm -f selectionparameters_*.root selection_summary_*.csv

  echo "<<< [$JOB][$CHANNEL] done  $(date)"
}

# --- Logging ---
LOGDIR="$OUTBASE/logs"
mkdir -p "$LOGDIR"

# --- Alle job_* durchgehen ---
shopt -s nullglob
JOBPATHS=("$INPDIR"/job_*)
if [ ${#JOBPATHS[@]} -eq 0 ]; then
  echo "Keine job_* Ordner in $INPDIR gefunden." >&2
  exit 1
fi

for JP in "${JOBPATHS[@]}"; do
  [ -d "$JP" ] || continue
  JOB="$(basename "$JP")"
  {
    echo "===== $JOB ====="
    run_one "$JOB" partialreco
    run_one "$JOB" fullreco
    run_one "$JOB" leptonrho
    echo "===== $JOB DONE ====="
  } > "$LOGDIR/${JOB}.log" 2>&1
done

echo "Alle Jobs fertig. Reports/Logs unter: $LOGDIR"
echo "Ergebnisse unter: $OUTBASE/neuDIS/{partialreco,fullreco,leptonrho}/job_*/"
