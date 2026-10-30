"""Rebuild panelB_data_regrouped.csv from panelB_data.csv and sra_run_metadata.csv.

Only Treatment-selection runs change:
  - PRJEB55794 (monoclonal-antibody treatment run) -> Clinical surveillance
  - PRJNA1035781, PRJNA692078 (drug-selection passages) -> Drug selection (cell culture)
All other rows are copied unchanged, with the BioProject added from the run metadata.
"""
import csv
import sys

IN_DATA = "panelB_data.csv"
IN_META = "sra_run_metadata.csv"
OUT = sys.argv[1] if len(sys.argv) > 1 else "panelB_data_regrouped.csv"

bioproject = {}
with open(IN_META, newline="") as f:
    for r in csv.DictReader(f):
        bioproject[r["Run"]] = r["bioproject"]

with open(IN_DATA, newline="") as f:
    rows = list(csv.DictReader(f))

fields = list(rows[0].keys()) + ["bioproject"]
with open(OUT, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=fields)
    w.writeheader()
    for r in rows:
        bp = bioproject.get(r["Run"], "")
        ctx = r["context"]
        if ctx == "Treatment-selection":
            ctx = ("Clinical surveillance" if bp == "PRJEB55794"
                   else "Drug selection (cell culture)")
        out = dict(r, context=ctx, bioproject=bp)
        w.writerow(out)
print(f"wrote {OUT}: {len(rows)} rows")
