
PATH_TO_FILES=$1
KEY=$2

#######################################################################################
source /cvmfs/ship.cern.ch/24.10/setUp.sh 
source /afs/cern.ch/user/a/anupamar/HTCondor/configfiles/config_ECN3_2024.sh #alienv load FairShip/latest-master-release > config_<version>.sh
echo 'config sourced'
#######################################################################################

OUT="combined_selectionparameters_${KEY}.root"
FLAVOURS=( all heliumCase vesselCase )

TMPDIR=$(mktemp -d)
trap 'rm -rf "$TMPDIR"' EXIT

[[ -f $OUT ]] && rm -f "$OUT"

###############################################################################
# 3) loop over flavours
###############################################################################
for FLAV in "${FLAVOURS[@]}"; do
  PAT="selectionparameters_${FLAV}.root"

  echo "=== [$FLAV] looking for $PAT ==="
  echo "find \"$PATH_TO_FILES\" -type f -path \"*/job_*/*\" -name \"$PAT\" -print0"
  
  mapfile -d '' FILES < <(find "$PATH_TO_FILES" -type f -path "*/job_*/*" -name "$PAT" -print0)

  echo ">>> found ${#FILES[@]} file(s)"
  if (( ${#FILES[@]} > 0 )); then
      printf '    %s\n' "${FILES[@]:0:3}"
      [[ ${#FILES[@]} -gt 3 ]] && echo "    …"
  fi

  (( ${#FILES[@]} )) || { echo "↪ skipping $FLAV (none found)" ; echo ; continue ; }


  #echo ">>> hadd command:"
  #echo "    hadd -j 8 -f \"$TMPDIR/$FLAV.root\"  <${#FILES[@]} input files>"
  #hadd -j 8 -f "${TMPDIR}/${FLAV}.root" "${FILES[@]}"
  # -------------------------------------------------------------------
  # build the list file instead of expanding a huge "$@" on the command line
  listfile="$TMPDIR/${FLAV}_files.txt"
  printf '%s\n' "${FILES[@]}" > "$listfile"

  echo ">>> hadd command:"
  echo "    hadd -j 8 -f \"$TMPDIR/$FLAV.root\" @$listfile"
  hadd -j 8 -f "${TMPDIR}/${FLAV}.root" @"$listfile"
  # -------------------------------------------------------------------

  # --- inside the FLAV loop ---------------------------------------------
  destdir="${FLAV}"          # e.g. all, heliumCase, vesselCase
  srcfile="${TMPDIR}/${FLAV}.root"

  # create the sub-directory once – does nothing if it is already there
  rootmkdir "${OUT}:${destdir}"
  echo ">>> rootcp → /${FLAV} in $OUT"
  # now copy *everything* from the source file into that directory
  rootcp -r "${srcfile}:/"   "${OUT}:${destdir}"
  echo
done

###############################################################################
# 2) done
###############################################################################
echo "✅  Output written to  $OUT"
