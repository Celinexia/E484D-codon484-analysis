"""Count GenBank genomes carrying each amino-acid change at spike codon 484.

Queries the LAPIS open endpoint behind CoV-Spectrum (sars_cov-2_nextstrain_open)
and saves each raw response, so the counts can be re-checked later.

Usage: python3 query_lapis_codon484.py [OUT_DIR]
Requires curl (the network here presents a custom TLS CA that curl trusts).
"""
import csv
import datetime
import json
import os
import subprocess
import sys

BASE = "https://lapis.cov-spectrum.org/open/v2/sample/aggregated"
OUT = sys.argv[1] if len(sys.argv) > 1 else "lapis_raw"
os.makedirs(OUT, exist_ok=True)

# Manuscript Methods, LAPIS query on 12 Oct 2026 (for comparison)
MS = {"total": 9438404, "E484A": 4306700, "E484K": 370685, "E484Q": 8237,
      "E484D": 229, "E484G": 170, "E484V": 0}

QUERIES = [("total", "")] + [(f"E484{aa}", f"?aminoAcidMutations=S:E484{aa}")
                             for aa in ("A", "K", "Q", "D", "G", "V")]

def fetch(query_string):
    url = BASE + query_string
    raw = subprocess.run(["curl", "-sS", "--fail", "--max-time", "120", url],
                         capture_output=True, check=True, text=True).stdout
    return url, json.loads(raw)

rows = []
accessed = datetime.date.today().isoformat()
for name, qs in QUERIES:
    url, resp = fetch(qs)
    with open(os.path.join(OUT, f"{name}.json"), "w") as f:
        json.dump({"url": url, "accessed": accessed, "response": resp}, f, indent=2)
    count = resp["data"][0]["count"]
    rows.append({"name": name, "count": count, "ms_12oct": MS[name],
                 "dataVersion": resp["info"]["dataVersion"],
                 "lapisVersion": resp["info"]["lapisVersion"],
                 "requestId": resp["info"]["requestId"], "url": url})

total = next(r["count"] for r in rows if r["name"] == "total")
for r in rows:
    if r["name"] != "total":
        r["freq"] = f"{r['count'] / total:.6%}"
with open(os.path.join(OUT, "summary.csv"), "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys()) + ["freq"])
    w.writeheader(); w.writerows(rows)

print(f"accessed {accessed}, dataVersion {rows[0]['dataVersion']}")
for r in rows:
    print(f"{r['name']:8s} now={r['count']:>10,}  12Oct={r['ms_12oct']:>10,}  "
          f"{'' if r['name']=='total' else r['freq']}")
