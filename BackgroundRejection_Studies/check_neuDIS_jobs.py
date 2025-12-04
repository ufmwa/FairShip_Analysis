#!/usr/bin/env python3
import os
import re
import sys

# Default base directory; kann als erstes Argument überschrieben werden
BASE_DIR = (
    "/storage/9/rquishpe/ship/"
    "NeutrinoDIS_2024helium_noCavern_neuDIS/7974071/neuDIS/fullreco"
)

if len(sys.argv) > 1:
    BASE_DIR = sys.argv[1]

job_dir_re = re.compile(r"^job_(\d+)$")
csv_re = re.compile(
    r"^selection_summary_(all|heliumCase|vesselCase)_job_(\d+)\.csv$"
)

# Per-directory Diagnosen
empty_jobs = []        # job_N directory has no files at all
misplaced_jobs = []    # directory contains CSVs for other job IDs
incomplete_jobs = []   # not all three expected CSVs for its own job ID

# Globale Info
jobs_expected = set()          # job IDs from directory names (job_*)
jobs_with_any_csv = set()      # job IDs that appear in at least one CSV filename
csv_map = {}                   # (job_id, variant) -> [file paths]


def check_job_dir(path, dirname):
    m = job_dir_re.match(dirname)
    if not m:
        return
    job_id_dir = m.group(1)  # ID from directory name
    jobs_expected.add(job_id_dir)

    entries = [e for e in os.scandir(path) if e.is_file()]
    csv_files = [e for e in entries if e.name.endswith(".csv")]
    root_files = [e for e in entries if e.name.endswith(".root")]

    # Fall 1: komplett leer
    if not csv_files and not root_files:
        empty_jobs.append(job_id_dir)
        return

    foreign_indices = set()
    own_csv_names = set()

    for entry in csv_files:
        fname = entry.name
        m2 = csv_re.match(fname)
        if not m2:
            # andere CSVs ignorieren wir
            continue
        variant, job_id_file = m2.groups()

        # global Tracking
        jobs_with_any_csv.add(job_id_file)
        key = (job_id_file, variant)
        csv_map.setdefault(key, []).append(entry.path)

        # per-directory Diagnose
        if job_id_file == job_id_dir:
            own_csv_names.add(fname)
        else:
            foreign_indices.add(job_id_file)

    # Erwartete CSV-Namen für dieses job_N im eigenen Ordner:
    expected_csv = {
        f"selection_summary_all_job_{job_id_dir}.csv",
        f"selection_summary_heliumCase_job_{job_id_dir}.csv",
        f"selection_summary_vesselCase_job_{job_id_dir}.csv",
    }

    # Fremde Job-IDs in diesem Ordner?
    if foreign_indices:
        misplaced_jobs.append(job_id_dir)

    # Fehlen eigene CSVs im eigenen Ordner?
    if not expected_csv.issubset(own_csv_names):
        missing = expected_csv - own_csv_names
        if missing:
            incomplete_jobs.append(job_id_dir)


def fmt_list(ids):
    """
    Hilfsfunktion: nimmt eine Iterable von string-Job-IDs
    und gibt einen String im Python-Listenformat zurück, z.B. '[1, 2, 10]'.
    """
    ids_sorted = sorted(set(ids), key=int)
    if not ids_sorted:
        return "[]"
    return "[" + ", ".join(ids_sorted) + "]"


def main():
    if not os.path.isdir(BASE_DIR):
        sys.stderr.write(f"Base dir does not exist: {BASE_DIR}\n")
        sys.exit(1)

    for entry in sorted(os.scandir(BASE_DIR), key=lambda e: e.name):
        if not entry.is_dir():
            continue
        check_job_dir(entry.path, entry.name)

    # ===== Coverage- & Varianten-Analyse =====
    jobs_expected_sorted = sorted(jobs_expected, key=int)
    jobs_with_any_csv_sorted = sorted(jobs_with_any_csv, key=int)

    # Erwartete Varianten (Fälle) pro Job
    expected_variants = {"all", "heliumCase", "vesselCase"}

    # Für jede erwartete Job-ID: welche Varianten haben wir irgendwo gefunden?
    job_variants = {job_id: set() for job_id in jobs_expected}
    for (job_id, variant), paths in csv_map.items():
        if job_id in job_variants and paths:
            job_variants[job_id].add(variant)

    # Jobs ohne CSVs (gar nichts gefunden)
    jobs_with_no_csv = sorted(
        [job_id for job_id, variants in job_variants.items()
         if len(variants) == 0],
        key=int,
    )

    # Jobs mit teilweise vorhandenen CSVs (mind. 1, aber nicht alle 3)
    jobs_with_partial_csv = sorted(
        [job_id for job_id, variants in job_variants.items()
         if 0 < len(variants) < len(expected_variants)],
        key=int,
    )

    # Jobs mit allen 3 CSVs
    jobs_with_complete_csv = sorted(
        [job_id for job_id, variants in job_variants.items()
         if variants == expected_variants],
        key=int,
    )

    # Zur Kompatibilität: jobs_missing = Jobs ohne CSVs
    jobs_missing = jobs_with_no_csv

    # Doppelte CSVs: gleiche (job_id, variant) mit mehr als einem Pfad
    duplicate_keys = {
        key: paths for key, paths in csv_map.items() if len(paths) > 1
    }
    jobs_with_duplicates = sorted(
        {job_id for (job_id, variant) in duplicate_keys.keys()},
        key=int,
    )

    # ===== Menschliche Zusammenfassung =====
    print("Base directory:", BASE_DIR)
    print()

    print("==== SUMMARY (human readable) ====")
    print(f"Total job_* directories (expected jobs): {len(jobs_expected_sorted)}")
    print(f"Jobs with any CSV (>=1 of 3):            {len(jobs_with_any_csv_sorted)}")
    print(f"Jobs with all 3 CSVs:                   {len(jobs_with_complete_csv)}")
    print(f"Jobs with partial CSVs (some missing):  {len(jobs_with_partial_csv)}")
    print(f"Jobs with no CSVs at all:               {len(jobs_with_no_csv)}")
    print(f"Jobs missing (alias of no CSVs):        {len(jobs_missing)}")
    print(f"Jobs with duplicated CSVs:              {len(jobs_with_duplicates)}")
    print(f"Empty job dirs (no files at all):       {len(set(empty_jobs))}")
    print(f"Misplaced job dirs (foreign IDs):       {len(set(misplaced_jobs))}")
    print(f"Incomplete job dirs (in own folder):    {len(set(incomplete_jobs))}")
    print()

    # ===== Python-kompatible Listen-Ausgabe =====
    print("==== PYTHON LISTS ====")
    print(f"jobs_expected = {fmt_list(jobs_expected)}")
    print(f"jobs_with_any_csv = {fmt_list(jobs_with_any_csv)}")
    print(f"jobs_with_complete_csv = {fmt_list(jobs_with_complete_csv)}")
    print(f"jobs_with_partial_csv = {fmt_list(jobs_with_partial_csv)}")
    print(f"jobs_with_no_csv = {fmt_list(jobs_with_no_csv)}")
    print(f"jobs_missing = {fmt_list(jobs_missing)}")
    print(f"jobs_with_duplicates = {fmt_list(jobs_with_duplicates)}")
    print(f"empty_jobs = {fmt_list(empty_jobs)}")
    print(f"misplaced_jobs = {fmt_list(misplaced_jobs)}")
    print(f"incomplete_jobs = {fmt_list(incomplete_jobs)}")


if __name__ == "__main__":
    main()


