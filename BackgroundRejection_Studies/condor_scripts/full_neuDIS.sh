#!/usr/bin/env bash
# Läuft neuDIS über alle job_* in <INPDIR> – nur partialreco & fullreco.

set -euo pipefail

if [ $# -ne 3 ]; then
  echo "Usage: $0 <INPDIR> <OUTBASE> <SCRIPTDIR>" >&2
  echo "  INPDIR   = INPUT mit job_* (Sim-Ausgabe von run_simScript/ship_reco)"
  echo "  OUTBASE  = Zielbasis für Analyse-Outputs (hier entsteht neuDIS/.../job_*)"
  echo "  SCRIPTDIR= Repo-Root mit BackgroundRejection_Studies/"
  exit 1
fi

INPDIR="$1"
OUTBASE="$2"
SCRIPTDIR="$3"

echo "[INFO] INPDIR   = $INPDIR"
echo "[INFO] OUTBASE  = $OUTBASE"
echo "[INFO] SCRIPTDIR= $SCRIPTDIR"

# --- job_* überall bis Tiefe 2 finden (deckt .../<run>/job_* und .../<grp>/<run>/job_* ab)
mapfile -t JOBPATHS < <(find "$INPDIR" -mindepth 1 -maxdepth 2 -type d -name 'job_*' | sort)
if [ ${#JOBPATHS[@]} -eq 0 ]; then
  echo "[ERROR] Keine job_* Ordner unter $INPDIR (bis 2 Ebenen) gefunden." >&2
  exit 1
fi
echo "[INFO] Found ${#JOBPATHS[@]} job folders"

LOGDIR="$OUTBASE/logs"
mkdir -p "$LOGDIR"

run_one() {
  local JOB="$1"         # job_XXXXXX
  local JOBPARENT="$2"   # Verzeichnis, das den job_* enthält
  local CHANNEL="$3"     # partialreco | fullreco

  local FLAG
  case "$CHANNEL" in
    partialreco) FLAG="--partialreco" ;;
    fullreco)    FLAG="--fullreco" ;;
    *) echo "Unknown channel: $CHANNEL" >&2; return 2 ;;
  esac

  # Sanity: existiert eine Geometrie im Input?
  if ! compgen -G "$JOBPARENT/$JOB/geofile_full*.root" >/dev/null && \
     ! compgen -G "$JOBPARENT/$JOB/*_rec.root" >/dev/null; then
    echo "[WARN] $JOB: keine Geometrie/REC im Input ($JOBPARENT/$JOB) – skip."
    return 0
  fi

  echo ">>> [$JOB][$CHANNEL] start $(date)"
  rm -f selectionparameters_*.root selection_summary_*.csv

  # WICHTIG: $FLAG NICHT quoten (kann mehrere Tokens enthalten)
  if ! python "$SCRIPTDIR/BackgroundRejection_Studies/run_neuDIS.py" \
        -p "$JOBPARENT" -i "$JOB" $FLAG ; then
    echo "!!! [$JOB][$CHANNEL] FAILED" >&2
    return 3
  fi

  local OUTDIR="$OUTBASE/neuDIS/$CHANNEL/$JOB"
  mkdir -p "$OUTDIR"
  cp selectionparameters_*.root selection_summary_*.csv "$OUTDIR"/ || true
  rm -f selectionparameters_*.root selection_summary_*.csv
  echo "<<< [$JOB][$CHANNEL] done $(date)"
}

for JP in "${JOBPATHS[@]}"; do
  [ -d "$JP" ] || continue
  JOB="$(basename "$JP")"
  JOBPARENT="$(dirname "$JP")"
  {
    echo "===== $JOB ====="
    run_one "$JOB" "$JOBPARENT" partialreco
    run_one "$JOB" "$JOBPARENT" fullreco
    echo "===== $JOB DONE ====="
  } > "$LOGDIR/${JOB}.log" 2>&1
done

echo "Alle Jobs fertig. Logs: $LOGDIR"
echo "Ergebnisse: $OUTBASE/neuDIS/{partialreco,fullreco}/job_*/"


