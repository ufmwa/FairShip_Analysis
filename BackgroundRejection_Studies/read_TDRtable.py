import sys
import pandas as pd
from tabulate import tabulate  # pip install tabulate

def norm(s: str) -> str:
    """Normalize spacing to make matching robust."""
    return " ".join(s.split())

def unique_order(seq):
    """Keep first occurrence, preserve order."""
    seen = set()
    out = []
    for x in seq:
        if x not in seen:
            seen.add(x)
            out.append(x)
    return out

def extract_unique_rows(input_file, rows_of_interest):
    # normalize target tags (and keep order without duplicates)
    targets = [norm(t) for t in rows_of_interest]
    targets = unique_order(targets)

    found = {}  # tag -> row (first one wins)
    with open(input_file, "r") as f:
        for line in f:
            parts = [p.strip() for p in line.split("│")]
            if len(parts) < 4:
                continue
            tag_raw = parts[1]
            tag = norm(tag_raw)
            if tag in targets and tag not in found:
                found[tag] = {
                    "Tag": tag_raw,  # keep original formatting for display
                    "nEvents_generated": parts[2],
                    "nEvents15y": parts[3],
                }

    # build dataframe in the same order as targets, skipping any missing
    rows = [found[t] for t in targets if t in found]
    return pd.DataFrame(rows)

def print_block(title, input_file, rows):
    print(title)
    df = extract_unique_rows(input_file, rows)
    if df.empty:
        print("(no matching rows found)")
    else:
        print(tabulate(df, headers="keys", tablefmt="github", showindex=False))
    print()  # blank line between blocks

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print(f"Usage: python {sys.argv[0]} <results.txt>")
        sys.exit(1)

    input_file = sys.argv[1]

    print("=================================================================")
    print_block("WITH PID", input_file, [
        "reconstructed",
        "preselection+PID",
        "preselection+AdvSBT@90MeV+PID",
        "preselection+AdvSBT@45MeV+PID",
        "preselection+[AdvSBT + GNNSBT ]@45MeV+ PID",
    ])

    print("-----------------------------------------------------------------")
    print_block("WITH PID + UBT", input_file, [
        "reconstructed",
        "preselection+UBT+PID",
        "preselection+UBT+AdvSBT@90MeV+PID",
        "preselection+UBT+AdvSBT@45MeV+PID",
        "preselection+UBT+[AdvSBT + GNNSBT ]@45MeV+ PID",
    ])

    print("=================================================================")
    print_block("WITHOUT PID", input_file, [
        "reconstructed",
        "preselection",
        "preselection+AdvSBT@90MeV",
        "preselection+AdvSBT@45MeV",
        "preselection+[AdvSBT + GNNSBT ]@45MeV",
    ])

    print("-----------------------------------------------------------------")
    print_block("WITHOUT PID + UBT", input_file, [
        "reconstructed",
        "preselection+UBT",
        "preselection+UBT+AdvSBT@90MeV",
        "preselection+UBT+AdvSBT@45MeV",
        "preselection+UBT+[AdvSBT + GNNSBT ]@45MeV",
    ])
    print("=================================================================")
