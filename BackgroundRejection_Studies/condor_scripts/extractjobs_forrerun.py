import re

logfile = "error_messages.txt"
outfile = "joblists_neuDIS_ECN3_2024_runs5_rerun_fullreco.whitespace"

with open(logfile, "r") as f:
    content = f.read()

# Extract all job_XXXX
jobs = re.findall(r'job_\d+', content)

# Make unique while preserving order
seen = set()
unique_jobs = []
for job in jobs:
    if job not in seen:
        seen.add(job)
        unique_jobs.append(job)

# Group into chunks of 5
chunks = [unique_jobs[i:i+5] for i in range(0, len(unique_jobs), 5)]

# Write to file
with open(outfile, "w") as out:
    for group in chunks:
        out.write(" ".join(group) + "\n")

print(f"✅ Extracted {len(unique_jobs)} unique job IDs → saved to {outfile}")
